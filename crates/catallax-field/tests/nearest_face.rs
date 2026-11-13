//! F8 ownership on the public target parts, checked against the B-rep built
//! from the same recipe.
//!
//! Owners at random B-rep surface points must be a name of the B-rep face.
//! Every B-rep edge must lie on an owner border within the display tolerance:
//! points half the tolerance away on each side are owned by the two adjoining
//! faces, so the border crosses between them. At the edge point itself every
//! name the query keeps must be within twice the width and include a name of
//! either face: a refusal naming both, a corner refusal naming the two planes
//! that the rounded point lies on exactly, or a proved side when binary64
//! rounding moved the point measurably off the edge. Answers and logical work
//! must be identical with 1, 2 and 4 workers.
use std::time::{Duration, Instant};

use catallax_field::{
    Direction, Plan,
    distance::{NearestOptions, Ownership, Prepared},
};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Origin, Recipe, Role, StepId};
use catallax_shapes::{CurveKind, Point3, RPoint3};
use catallax_topology::{FaceId, Model};

/// The verification-request allowance, per query.
const ALLOWANCE: u64 = 100_000_000;
/// Half the default 0.01 mm display tolerance.
const STEP: f64 = 0.005;
/// Surface samples keep this far from every B-rep edge polyline, so mesh chords
/// and carrier projection cannot carry them across a trim.
const CLEARANCE: f64 = 0.05;

const P1: &str = include_str!("../../../parts/P1-bracket.json");
const P2: &str = include_str!("../../../parts/P2-collar.json");
const P3: &str = include_str!("../../../parts/P3-filleted-boss.json");
const P4: &str = include_str!("../../../parts/P4-shelled-tray.json");

fn point(p: [f64; 3]) -> Point3 {
    let [x, y, z] = p.map(|c| Rational::from_f64_exact(c).unwrap());
    RPoint3::new(x, y, z).to_real()
}

#[derive(Clone, Debug)]
enum Expect {
    /// Any name of one B-rep face.
    Owner(Vec<Origin>),
    /// The names of the two B-rep faces that meet at an edge.
    Border([Vec<Origin>; 2]),
}

#[derive(Clone, Debug)]
struct Sample {
    point: [f64; 3],
    expect: Expect,
}

struct Brep {
    model: Model,
    faces: Vec<FaceId>,
    mesh: Mesh,
}

fn brep(recipe: &Recipe) -> Brep {
    let (built, work) = with_operation_budget(ALLOWANCE, || {
        catallax_features::evaluate(recipe, &mut Decider::default()).unwrap()
    });
    assert!(!work.exhausted);
    assert!(built.model.check(built.result).is_empty());
    let mesh = mesh_body(
        &built.model,
        built.result,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(mesh.marks.is_empty());
    Brep {
        faces: built.model.faces_of_body(built.result),
        model: built.model,
        mesh: mesh.value,
    }
}

// Deterministic sample order, independent of the platform's hashing.
fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn unit(state: &mut u64) -> f64 {
    (splitmix(state) >> 11) as f64 / (1u64 << 53) as f64
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(b, a);
    let t = (dot(sub(p, a), d) / dot(d, d).max(f64::MIN_POSITIVE)).clamp(0.0, 1.0);
    let q: [f64; 3] = std::array::from_fn(|i| a[i] + t * d[i]);
    dot(sub(p, q), sub(p, q)).sqrt()
}

/// Area-weighted random points of the B-rep display mesh, at least `CLEARANCE`
/// from every edge polyline, projected onto their B-rep face's exact carrier.
fn surface_samples(brep: &Brep, count: usize) -> Vec<Sample> {
    let positions = &brep.mesh.positions;
    let segments: Vec<[[f64; 3]; 2]> = brep
        .mesh
        .edges
        .iter()
        .flat_map(|e| {
            e.vertices
                .windows(2)
                .map(|w| [positions[w[0]], positions[w[1]]])
        })
        .collect();
    let mut cumulative = Vec::with_capacity(brep.mesh.triangles.len());
    let mut total = 0.0;
    for t in &brep.mesh.triangles {
        let [a, b, c] = t.vertices.map(|v| positions[v]);
        let (u, v) = (sub(b, a), sub(c, a));
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        total += dot(n, n).sqrt() / 2.0;
        cumulative.push(total);
    }
    let mut dec = Decider::default();
    let mut state = 0xf8;
    let mut samples = Vec::with_capacity(count);
    let mut attempts = 0;
    while samples.len() < count {
        attempts += 1;
        assert!(attempts < 100 * count, "too few clear surface points");
        let target = unit(&mut state) * total;
        let index = cumulative
            .partition_point(|&c| c <= target)
            .min(cumulative.len() - 1);
        let triangle = &brep.mesh.triangles[index];
        let [a, b, c] = triangle.vertices.map(|v| positions[v]);
        let s = unit(&mut state).sqrt();
        let r = unit(&mut state);
        let q: [f64; 3] =
            std::array::from_fn(|i| (1.0 - s) * a[i] + s * (1.0 - r) * b[i] + s * r * c[i]);
        if segments
            .iter()
            .any(|[a, b]| segment_distance(q, *a, *b) < CLEARANCE)
        {
            continue;
        }
        let surface = brep
            .model
            .geom()
            .surface(brep.model.face(triangle.face).surface);
        let projected = catallax_shapes::distance::project(surface, &point(q), &mut dec).unwrap();
        let p = [&projected.x, &projected.y, &projected.z].map(|x| {
            dec.enclose_budgeted_at(x, 128)
                .unwrap()
                .mid()
                .to_f64_nearest()
        });
        samples.push(Sample {
            point: p,
            expect: Expect::Owner(brep.model.face_origins(triangle.face).to_vec()),
        });
    }
    samples
}

#[derive(Default, Debug)]
struct EdgeCount {
    checked: usize,
    short_arcs: usize,
    seams: usize,
    shared: usize,
}

/// For each B-rep edge between two faces with distinct names, `per_edge`
/// points on the edge: interior polyline vertices of a curved edge, evenly
/// spaced points of a straight one, or the end vertices of a short arc. Each gives the edge point itself and one
/// point `STEP` from it toward the centroid of an adjoining triangle of each face.
fn edge_samples(brep: &Brep, per_edge: usize) -> (Vec<Sample>, EdgeCount) {
    let positions = &brep.mesh.positions;
    let mut samples = Vec::new();
    let mut count = EdgeCount::default();
    for edge in &brep.mesh.edges {
        let mut faces = Vec::new();
        for coedge in brep.model.coedges_of_edge(edge.edge) {
            // The model also keeps the recipe's intermediate bodies.
            if let Some((_, face)) = brep.model.owner_of_coedge(coedge)
                && brep.faces.contains(&face)
                && !faces.contains(&face)
            {
                faces.push(face);
            }
        }
        if faces.len() == 1 {
            // A periodic seam joins one face to itself: no owner border.
            count.seams += 1;
            continue;
        }
        assert_eq!(faces.len(), 2, "a manifold edge joins two faces");
        let names = [0, 1].map(|i| brep.model.face_origins(faces[i]).to_vec());
        if names[0].iter().any(|n| names[1].contains(n)) {
            // Faces that share a name have no border between names.
            count.shared += 1;
            continue;
        }
        count.checked += 1;
        let n = edge.vertices.len();
        let straight = brep
            .model
            .geom()
            .curve(brep.model.edge(edge.edge).curve)
            .kind()
            == CurveKind::Line;
        for j in 0..per_edge {
            // The polyline segment holding the sample, and its position on it.
            // A curved edge short enough to need no interior vertex is sampled
            // at its exact end vertices instead of a chord point off the curve.
            let (segment, e) = if !straight && n == 2 {
                count.short_arcs += usize::from(j == 0);
                ((0, 1), positions[edge.vertices[j % 2]])
            } else if straight {
                let [a, b] = [edge.vertices[0], edge.vertices[n - 1]].map(|v| positions[v]);
                let t = (j + 1) as f64 / (per_edge + 1) as f64;
                let k = ((t * (n - 1) as f64) as usize).min(n - 2);
                (
                    (k, k + 1),
                    std::array::from_fn(|i| a[i] + t * (b[i] - a[i])),
                )
            } else {
                let k = 1 + j * (n - 2) / per_edge;
                ((k, k), positions[edge.vertices[k]])
            };
            let ends = [edge.vertices[segment.0], edge.vertices[segment.1]];
            samples.push(Sample {
                point: e,
                expect: Expect::Border(names.clone()),
            });
            for (face, name) in faces.iter().zip(&names) {
                let triangle = brep
                    .mesh
                    .triangles
                    .iter()
                    .find(|t| t.face == *face && ends.iter().all(|v| t.vertices.contains(v)))
                    .expect("an adjoining triangle on each side");
                let c: [f64; 3] = std::array::from_fn(|i| {
                    triangle
                        .vertices
                        .iter()
                        .map(|&w| positions[w][i])
                        .sum::<f64>()
                        / 3.0
                });
                let d = sub(c, e);
                let length = dot(d, d).sqrt();
                samples.push(Sample {
                    point: std::array::from_fn(|i| e[i] + STEP * d[i] / length),
                    expect: Expect::Owner(name.clone()),
                });
            }
        }
    }
    (samples, count)
}

struct Outcome {
    answer: Ownership,
    work: u64,
    elapsed: Duration,
}

fn query(prepared: &Prepared<'_>, p: [f64; 3]) -> Outcome {
    let started = Instant::now();
    let (answer, work) = with_operation_budget(ALLOWANCE, || {
        prepared.nearest_face(
            &point(p),
            &NearestOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted, "{p:?} exhausted the allowance");
    Outcome {
        answer: answer.unwrap(),
        work: work.used,
        elapsed: started.elapsed(),
    }
}

// Everything a worker could change: owner, bounds, witness, margin, cells, work.
fn fingerprint(outcome: &Outcome) -> String {
    match &outcome.answer {
        Ownership::Proved(o) => format!(
            "proved {:?} face={} [{:?}, {:?}] witness={:?} direction={:?} margin={:?} next={:?} cells={} work={}",
            o.origin,
            o.face.index(),
            o.lower,
            o.upper,
            o.witness.as_rational(),
            o.direction,
            o.margin,
            o.limiting.as_ref().map(|r| (&r.origin, &r.lower, &r.upper)),
            o.cells,
            outcome.work
        ),
        Ownership::Refused(a) => format!(
            "refused {:?} spread={:?} cells={} work={}",
            a.contenders, a.spread, a.cells, outcome.work
        ),
    }
}

fn run(prepared: &Prepared<'_>, samples: &[Sample], workers: usize) -> Vec<Outcome> {
    let mut slots: Vec<Option<Outcome>> = (0..samples.len()).map(|_| None).collect();
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|worker| {
                scope.spawn(move || {
                    samples
                        .iter()
                        .enumerate()
                        .skip(worker)
                        .step_by(workers)
                        .map(|(i, s)| (i, query(prepared, s.point)))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        for handle in handles {
            for (i, outcome) in handle.join().unwrap() {
                slots[i] = Some(outcome);
            }
        }
    });
    slots.into_iter().map(Option::unwrap).collect()
}

#[derive(Default, Debug)]
struct Tally {
    owners: usize,
    owner_matches: usize,
    // Surface points where two names of one merged B-rep face tie at their seam.
    merged_refusals: usize,
    borders: usize,
    border_matches: usize,
    // Edge points refused naming both sides; refused naming one side, at a
    // corner whose rounded coordinates lie exactly on two planes but provably
    // off a curved third face; or proved for one side, when binary64 rounding
    // moved the point measurably off the edge.
    border_refusals: usize,
    border_corners: usize,
    border_sides: usize,
    // Edge-point refusals that also keep a witness-less buried candidate.
    border_buried: usize,
    misses: Vec<String>,
    max_work: u64,
    total: Duration,
    slowest: Duration,
}

fn judge(samples: &[Sample], outcomes: &[Outcome], width: &Rational) -> Tally {
    let close = width.mul_pow2(1);
    let mut tally = Tally::default();
    for (sample, outcome) in samples.iter().zip(outcomes) {
        tally.max_work = tally.max_work.max(outcome.work);
        tally.total += outcome.elapsed;
        tally.slowest = tally.slowest.max(outcome.elapsed);
        let matched = match (&sample.expect, &outcome.answer) {
            (Expect::Owner(names), answer) => {
                tally.owners += 1;
                match answer {
                    Ownership::Proved(o) if names.contains(&o.origin) => {
                        tally.owner_matches += 1;
                        true
                    }
                    Ownership::Refused(a)
                        if names.len() > 1
                            && a.contenders.iter().all(|c| names.contains(&c.origin))
                            && a.spread.as_ref().is_some_and(|s| s <= &close) =>
                    {
                        tally.owner_matches += 1;
                        tally.merged_refusals += 1;
                        true
                    }
                    _ => false,
                }
            }
            (Expect::Border(sides), answer) => {
                tally.borders += 1;
                // Contenders with an admitted witness. A refusal may also keep a
                // buried candidate whose source face passes through the point;
                // it has no witness and cannot be excluded there.
                let witnessed = |a: &catallax_field::distance::Ambiguity| {
                    a.contenders
                        .iter()
                        .filter(|c| c.upper.is_some())
                        .map(|c| c.origin.clone())
                        .collect::<Vec<_>>()
                };
                match answer {
                    Ownership::Refused(a)
                        if a.contenders
                            .iter()
                            .all(|c| c.upper.as_ref().is_none_or(|u| u <= &close))
                            && sides
                                .iter()
                                .any(|s| witnessed(a).iter().any(|n| s.contains(n))) =>
                    {
                        tally.border_matches += 1;
                        let kept = witnessed(a);
                        if sides.iter().all(|s| kept.iter().any(|n| s.contains(n))) {
                            tally.border_refusals += 1;
                        } else {
                            tally.border_corners += 1;
                        }
                        if kept.len() < a.contenders.len() {
                            tally.border_buried += 1;
                        }
                        true
                    }
                    Ownership::Proved(o)
                        if sides.iter().any(|s| s.contains(&o.origin)) && o.upper <= close =>
                    {
                        tally.border_matches += 1;
                        tally.border_sides += 1;
                        true
                    }
                    _ => false,
                }
            }
        };
        if !matched {
            tally.misses.push(format!(
                "{:?}: expected {:?}, got {}",
                sample.point,
                sample.expect,
                fingerprint(outcome)
            ));
        }
    }
    tally
}

fn check_part(name: &str, text: &str, surface: usize, per_edge: usize, workers: &[usize]) {
    let recipe = Recipe::parse(text).unwrap();
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let prepared = plan.prepare_distances(&mut dec).unwrap();
    prepared.prepare_nearest(&mut dec).unwrap();
    let brep = brep(&recipe);
    let mut samples = surface_samples(&brep, surface);
    let (edges, count) = edge_samples(&brep, per_edge);
    samples.extend(edges);
    let width = NearestOptions::default().width;
    let mut reference: Option<Vec<String>> = None;
    for &workers in workers {
        let started = Instant::now();
        let outcomes = run(&prepared, &samples, workers);
        let wall = started.elapsed();
        let tally = judge(&samples, &outcomes, &width);
        // Informational shared-machine timings, not a speed qualification.
        eprintln!(
            "{name} ownership workers={workers}: owners {}/{} ({} merged-name ties), \
             edge points {}/{} ({} refused naming both, {} corner refusals naming one, \
             {} proved one side; {} keep a buried candidate), edges {:?}, \
             max work {}, query sum {:?}, slowest {:?}, wall {wall:?}",
            tally.owner_matches,
            tally.owners,
            tally.merged_refusals,
            tally.border_matches,
            tally.borders,
            tally.border_refusals,
            tally.border_corners,
            tally.border_sides,
            tally.border_buried,
            count,
            tally.max_work,
            tally.total,
            tally.slowest,
        );
        let prints: Vec<String> = outcomes.iter().map(fingerprint).collect();
        // Raw evidence for the record; not part of the check.
        if let Some(dir) = std::env::var_os("CATALLAX_OWNERSHIP_EVIDENCE") {
            let lines: Vec<String> = samples
                .iter()
                .zip(&prints)
                .map(|(s, p)| format!("{:?} {:?} => {p}", s.point, s.expect))
                .collect();
            std::fs::write(
                std::path::Path::new(&dir)
                    .join(format!("{name}-{surface}-{per_edge}-w{workers}.txt")),
                lines.join("\n"),
            )
            .unwrap();
        }
        assert!(tally.misses.is_empty(), "{name}: {:#?}", tally.misses);
        assert_eq!(tally.owners, surface + 2 * count.checked * per_edge);
        assert_eq!(tally.borders, count.checked * per_edge);
        match &reference {
            None => reference = Some(prints),
            Some(first) => assert_eq!(
                first, &prints,
                "{name}: answers differ with {workers} workers"
            ),
        }
    }
}

/// A sample of the P3 pass mark sized for every check run.
#[test]
fn p3_owners_match_brep_faces_and_edges_lie_on_owner_borders() {
    check_part("P3", P3, 24, 1, &[1, 2, 4]);
}

/// The full pass mark: 1,000 B-rep surface points and four points per edge on
/// P3, then the other analytic target parts.
#[test]
#[ignore = "full F8 pass mark, thousands of proved queries; run on demand"]
fn full_ownership_pass_mark() {
    check_part("P3", P3, 1000, 4, &[1, 2, 4]);
    check_part("P2", P2, 200, 2, &[1, 2, 4]);
    check_part("P4", P4, 200, 2, &[1, 2, 4]);
    check_part("P1", P1, 200, 2, &[1, 2, 4]);
}

/// The other analytic target parts that build quickly, at every check run.
#[test]
fn p2_and_p4_owners_match_brep_faces_and_edges_lie_on_owner_borders() {
    check_part("P2", P2, 12, 1, &[1, 2, 4]);
    check_part("P4", P4, 12, 1, &[1, 2, 4]);
}

fn p3_prepared(plan: &Plan) -> Prepared<'_> {
    plan.prepare_distances(&mut Decider::default()).unwrap()
}

#[test]
fn ownership_reports_distance_direction_margin_and_refusals() {
    let plan = Plan::compile(&Recipe::parse(P3).unwrap()).unwrap();
    let prepared = p3_prepared(&plan);
    let bottom = Origin {
        step: StepId::new("block").unwrap(),
        role: Role::Start,
    };

    // Below the block: three millimetres from the bottom face, straight up.
    let Ownership::Proved(owner) = query(&prepared, [5.0, 6.0, -3.0]).answer else {
        panic!("a point below the bottom face has one owner");
    };
    assert_eq!(owner.origin, bottom);
    assert_eq!(owner.lower, Rational::from_i64(3));
    assert_eq!(owner.upper, Rational::from_i64(3));
    let Direction::Unique(direction) = owner.direction else {
        panic!("an off-surface point has a direction");
    };
    for (axis, expected) in direction.iter().zip([0.0, 0.0, 1.0]) {
        assert!(axis.contains(expected), "{direction:?}");
    }
    // The limiting name's lower bound exceeds the owner's upper bound; the
    // margin is half that proven gap.
    let next = owner.limiting.unwrap();
    assert!(next.lower > owner.upper);
    assert_eq!(
        owner.margin.unwrap(),
        (&next.lower - &owner.upper).mul_pow2(-1)
    );

    // On the bottom face itself the witness is the point; no direction.
    let Ownership::Proved(owner) = query(&prepared, [5.0, 6.0, 0.0]).answer else {
        panic!("a bottom-face point has one owner");
    };
    assert_eq!(owner.origin, bottom);
    assert_eq!(owner.upper, Rational::zero());
    assert_eq!(owner.direction, Direction::Undefined);

    // On the edge between the bottom and the y_min side the two names are
    // equidistant: a typed refusal naming both, never a guess.
    let Ownership::Refused(ambiguity) = query(&prepared, [5.0, 0.0, 0.0]).answer else {
        panic!("an edge point has no single owner");
    };
    let names: Vec<_> = ambiguity
        .contenders
        .iter()
        .map(|c| c.origin.clone())
        .collect();
    assert_eq!(names.len(), 2, "{names:?}");
    assert!(names.contains(&bottom));
    assert!(
        names
            .iter()
            .any(|n| matches!(&n.role, Role::Side { segment, .. } if segment.as_str() == "y_min"))
    );
    assert_eq!(ambiguity.spread, Some(Rational::zero()));
    assert_eq!(ambiguity.refusal().code, "field.ownership-ambiguous");

    // Outside a convex edge the nearest point is the edge, so the refusal
    // covers the whole wedge, not only the edge itself.
    let Ownership::Refused(ambiguity) = query(&prepared, [5.0, -2.0, -2.0]).answer else {
        panic!("a point whose nearest point is an edge has no single owner");
    };
    assert_eq!(ambiguity.contenders.len(), 2);

    // Known limit: inside the part, close to the buried boss bottom disk, that
    // registry candidate has no witness and its source disk is nearer than any
    // surviving face, so the query refuses rather than guess.
    let Ownership::Refused(ambiguity) = query(&prepared, [20.0, 15.0, 9.5]).answer else {
        panic!("a buried source face near the point is not excluded");
    };
    eprintln!("buried-disk limit: {ambiguity:?}");
    assert!(ambiguity.contenders.iter().any(|c| c.upper.is_none()
        && c.origin.step.as_str() == "boss"
        && c.origin.role == Role::Start));
    assert_eq!(ambiguity.spread, None);

    // Invalid controls are a refusal, not a default.
    let controls = NearestOptions {
        max_cells: 1,
        ..NearestOptions::default()
    };
    let refusal = prepared
        .nearest_face(&point([1.0, 1.0, 1.0]), &controls, &mut Decider::default())
        .unwrap_err();
    assert_eq!(refusal.code, "field.ownership-options");
}

#[test]
fn ownership_exhaustion_is_a_refusal() {
    let plan = Plan::compile(&Recipe::parse(P3).unwrap()).unwrap();
    let prepared = p3_prepared(&plan);
    let (answer, work) = with_operation_budget(1_000, || {
        prepared.nearest_face(
            &point([5.0, 6.0, -3.0]),
            &NearestOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted);
    assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
}
