//! Extruded slot top boundary: a straight strip runs past its rail ends and is
//! cut by the two cylindrical end walls, with exact volumes, terminal curves,
//! radius-independent names, an immutable input after refusal and a closed
//! display at the application's default limits.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::OpenOptions,
    io::Write,
    time::{Duration, Instant},
};

use catallax_boolean::{Containment, classify_point};
use catallax_features::{BuildError, Built, evaluate, evaluate_incremental};
use catallax_mesh::{Mesh, MeshOptions, PreviewOptions, mesh_body, preview_mesh_body};
use catallax_numbers::{
    BudgetPolicy, Decider, Outcome, Rational, Real, StrictSign, with_operation_budget,
};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::{Curve, Point3, RPoint3, Surface};
use catallax_topology::{Loop, Sense};
use serde_json::json;

fn recipe(segments: &[&str], radius: &str, version: u8) -> Recipe {
    let mut steps = vec![
        json!({"id":"profile","op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
            "loops":[[
                {"kind":"line","id":"bottom","to":["20","-5"]},
                {"kind":"arc","id":"right","centre":["20","0"],"to":["20","5"],"turn":"ccw"},
                {"kind":"line","id":"top","to":["-20","5"]},
                {"kind":"arc","id":"left","centre":["-20","0"],"to":["-20","-5"],"turn":"ccw"}]]}),
        json!({"id":"body","op":"extrude","profile":"profile","distance":"10"}),
    ];
    let mut result = "body";
    if !segments.is_empty() {
        let edges: Vec<_> = segments
            .iter()
            .map(|segment| {
                json!({"between":[{"step":"body","role":"end"},
                    {"step":"body","role":{"side":{"segment":segment}}}]})
            })
            .collect();
        steps.push(json!({"id":"round","op":"fillet","body":"body","edges":edges,"radius":radius}));
        result = "round";
    }
    Recipe::parse(
        &json!({"format":format!("catallax-recipe/{version}"),"name":"slot tangent ends",
            "steps":steps,"result":result})
        .to_string(),
    )
    .unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 100_000_000,
        ..Default::default()
    })
}
fn zero(value: Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(&value), Outcome::ProvenCoincident(_))
}
fn positive(value: &Real, dec: &mut Decider) -> bool {
    matches!(dec.sign(value), Outcome::Proven(StrictSign::Positive))
}
fn int(v: i64) -> Real {
    Real::from_i64(v)
}
fn export(built: &Built) -> Vec<u8> {
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    bytes
}
fn assert_volume(built: &Built, want: &str) {
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &Default::default(),
        &mut decider(),
    )
    .unwrap()
    .value;
    let want = Rational::parse(want).unwrap();
    let tolerance = Rational::parse("1e-5").unwrap();
    assert!(&volume.hi - &volume.lo <= tolerance, "{volume:?}");
    assert!(
        &want - &tolerance <= volume.lo && volume.hi <= &want + &tolerance,
        "{} not within 1e-5 of {}",
        volume.lo.to_f64_nearest(),
        want.to_f64_nearest()
    );
}
fn classify(built: &mut Built, points: &[([&str; 3], Containment)]) {
    for (p, want) in points {
        let [x, y, z] = p.map(|v| Rational::parse(v).unwrap());
        let point = RPoint3::new(x, y, z);
        let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
        assert!(
            matches!(actual, Outcome::Proven(actual) if actual == *want),
            "at {p:?}: {actual:?}"
        );
    }
}
/// The sorted origins of every face, the blend's included.
fn origins(built: &Built) -> Vec<Vec<catallax_recipe::Origin>> {
    let mut all: Vec<_> = built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .map(|f| built.model.face_origins(f).to_vec())
        .collect();
    all.sort();
    all
}

#[test]
fn a_straight_strip_ends_on_both_cylindrical_end_walls() {
    for (radius, volume) in [
        (1, "4776.244436115815"),
        (2, "4747.891984181232"),
        (4, "4630.739296215008"),
    ] {
        let mut built = evaluate(&recipe(&["top"], &radius.to_string(), 3), &mut decider())
            .unwrap_or_else(|e| panic!("r={radius}: {e}"));
        let model = &built.model;
        assert_eq!(model.check(built.result), vec![]);
        let lumps = &model.body(built.result).lumps;
        assert_eq!(lumps.len(), 1);
        assert_eq!(model.lump(lumps[0]).shells.len(), 1);
        let faces = model.faces_of_body(built.result);
        let mut kinds = BTreeMap::new();
        let mut blends = Vec::new();
        for face in &faces {
            let surface = model.geom().surface(model.face(*face).surface);
            *kinds.entry(format!("{:?}", surface.kind())).or_insert(0) += 1;
            if model
                .face_origins(*face)
                .iter()
                .any(|o| matches!(o.role, Role::Blend { .. }))
            {
                blends.push(*face);
            }
        }
        assert_eq!(
            kinds,
            BTreeMap::from([("Cylinder".into(), 3), ("Plane".into(), 4)])
        );
        let [blend] = blends[..] else {
            panic!("one blend face, not {blends:?}");
        };
        let mut dec = decider();
        let r = int(radius);
        let Surface::Cylinder(c) = model.geom().surface(model.face(blend).surface) else {
            panic!("the strip is a cylinder");
        };
        assert!(zero(&c.radius_sq - &r.square(), &mut dec));
        assert!(zero(c.axis.y.clone(), &mut dec) && zero(c.axis.z.clone(), &mut dec));
        assert!(zero(&c.origin.y - &(int(5) - r.clone()), &mut dec));
        assert!(zero(&c.origin.z - &(int(10) - r.clone()), &mut dec));

        // Every vertex lies on every face it bounds; every edge has two
        // opposite uses.
        let mut uses: BTreeMap<_, Vec<Sense>> = BTreeMap::new();
        for face in &faces {
            let surface = model.face(*face).surface;
            for lp in &model.face(*face).loops {
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    continue;
                };
                for c in cs {
                    let c = model.coedge(*c);
                    uses.entry(c.edge).or_default().push(c.sense);
                    let e = model.edge(c.edge);
                    for v in [e.start, e.end].into_iter().flatten() {
                        let p = model.geom().position(model.vertex(v).point);
                        assert!(matches!(
                            model.geom().side_of(surface, &p, &mut dec),
                            Outcome::ProvenCoincident(_)
                        ));
                    }
                }
            }
        }
        for senses in uses.values() {
            assert_eq!(senses.len(), 2);
            assert_ne!(senses[0], senses[1]);
        }

        // Both terminal curves are blend/end-wall intersections. They run
        // from the wall contact (+-20, 5, 10-r) to the top contact
        // (+-(20 + sqrt(10r - r^2)), 5-r, 10), beyond the straight segment.
        let mut terminals = 0;
        let mut ends = BTreeSet::new();
        let Loop::Coedges(cs) = model.loop_(model.face(blend).loops[0]) else {
            panic!("the strip has a coedge boundary");
        };
        for c in cs {
            let e = model.edge(model.coedge(*c).edge);
            if matches!(model.geom().curve(e.curve), Curve::Intersection(_)) {
                terminals += 1;
            }
            for v in [e.start, e.end].into_iter().flatten() {
                let p: Point3 = model.geom().position(model.vertex(v).point);
                let side = if positive(&p.x, &mut dec) { 1 } else { -1 };
                let offset = &(&p.x * &int(side)) - &int(20);
                let wall = zero(&p.y - &int(5), &mut dec)
                    && zero(&p.z - &(int(10) - r.clone()), &mut dec)
                    && zero(offset.clone(), &mut dec);
                let top = zero(&p.y - &(int(5) - r.clone()), &mut dec)
                    && zero(&p.z - &int(10), &mut dec)
                    && positive(&offset, &mut dec)
                    && zero(
                        &offset.square() - &(&int(10 * radius) - &r.square()),
                        &mut dec,
                    );
                assert!(wall ^ top, "r={radius}: {p:?}");
                ends.insert((side, top));
            }
        }
        assert_eq!(terminals, 2);
        assert_eq!(ends.len(), 4);
        assert_volume(&built, volume);
        // The strip extends past x=+-20, while the opposite edge stays sharp.
        classify(
            &mut built,
            &[
                (["0", "49/10", "99/10"], Containment::Outside),
                (["22", "9/2", "99/10"], Containment::Outside),
                (["0", "-49/10", "99/10"], Containment::Inside),
                (["249/10", "0", "99/10"], Containment::Inside),
                (["0", "0", "1"], Containment::Inside),
            ],
        );
    }
}

/// `catallax_app::DEFAULT_BUILD_WORK`: construction and display share it.
const APPLICATION_WORK: u64 = 20_000_000;
/// Speed target for one display mesh. Clock time is reported and never fails a check.
const MESH_TARGET: Duration = Duration::from_secs(1);

/// SolidDraft's default display: the default chord and a 15 degree normal limit.
fn display_options() -> MeshOptions {
    MeshOptions {
        angular_deviation: Some(Rational::parse("0.2617993877991494").unwrap()),
        ..MeshOptions::default()
    }
}

/// A passing test's output is hidden, so CI also gets the warning in its step summary.
fn report_warning(warning: &str) {
    eprintln!("warning: {warning}");
    let summary = std::env::var_os("GITHUB_STEP_SUMMARY")
        .and_then(|path| OpenOptions::new().append(true).open(path).ok());
    if let Some(mut summary) = summary {
        let _ = writeln!(summary, "- warning: {warning}");
    }
}

/// Builds and meshes one recipe inside a single application allowance, as
/// `catallax_app::build` does.
fn display(recipe: &Recipe, options: &MeshOptions, case: &str) -> (Built, Mesh) {
    let mut dec = Decider::default();
    let ((built, mesh, elapsed), work) = with_operation_budget(APPLICATION_WORK, || {
        let built = evaluate(recipe, &mut dec).unwrap_or_else(|e| panic!("{case}: {e}"));
        let start = Instant::now();
        let mesh = mesh_body(&built.model, built.result, options, &mut dec);
        (built, mesh, start.elapsed())
    });
    assert!(!work.exhausted, "{case} exhausted the work cap: {work:?}");
    if elapsed >= MESH_TARGET {
        report_warning(&format!(
            "{case} mesh took {elapsed:?}, target {MESH_TARGET:?}"
        ));
    }
    let mesh = mesh.unwrap_or_else(|e| panic!("{case}: {e:?}")).value;
    (built, mesh)
}

#[test]
fn the_rounded_slot_meshes_closed_at_the_default_display_limits() {
    let options = display_options();
    let (_, plain) = display(&recipe(&[], "0", 3), &options, "plain slot");
    for radius in ["1", "4"] {
        let case = format!("r={radius}");
        let (built, mesh) = display(&recipe(&["top"], radius, 3), &options, &case);
        let model = &built.model;
        assert!(mesh.unavailable_controls.is_empty(), "{case}");
        let faces = model.faces_of_body(built.result);
        let meshed: BTreeSet<_> = mesh.faces.iter().map(|f| f.face).collect();
        assert_eq!(
            meshed,
            faces.iter().copied().collect::<BTreeSet<_>>(),
            "{case}"
        );
        assert_eq!(mesh.edges.len(), model.edges_of_body(built.result).len());

        // Closed and oriented: every segment has two opposite uses, including
        // every shared edge polyline, and the genus-0 surface has V - E + F = 2.
        let mut incidence = BTreeMap::<_, (usize, i32)>::new();
        for triangle in &mesh.triangles {
            for i in 0..3 {
                let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
                let uses = incidence.entry((a.min(b), a.max(b))).or_default();
                uses.0 += 1;
                uses.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(incidence.values().all(|uses| *uses == (2, 0)), "{case}");
        assert_eq!(mesh.triangles.len(), 2 * mesh.positions.len() - 4, "{case}");
        for edge in &mesh.edges {
            assert!(edge.max_chord_error <= options.chord_tolerance, "{case}");
            for pair in edge.vertices.windows(2) {
                let key = (pair[0].min(pair[1]), pair[0].max(pair[1]));
                assert_eq!(incidence.get(&key), Some(&(2, 0)), "{case}");
            }
        }
        let polylines: BTreeMap<_, _> = mesh.edges.iter().map(|e| (e.edge, e)).collect();
        for face in &mesh.faces {
            assert!(face.max_chord_error <= options.chord_tolerance, "{case}");
            assert_eq!(face.max_angular_deviation, options.angular_deviation);
            // Both end walls and the strip are ruled strips between their
            // boundary samples: a disc with S boundary segments and no
            // interior vertex has S - 2 triangles.
            let definition = model.face(face.face);
            if !matches!(
                model.geom().surface(definition.surface),
                Surface::Cylinder(_)
            ) {
                continue;
            }
            let [boundary] = definition.loops[..] else {
                panic!("{case}: a cylindrical face has one loop");
            };
            let Loop::Coedges(cs) = model.loop_(boundary) else {
                panic!("{case}: a cylindrical face has a coedge boundary");
            };
            let segments: usize = cs
                .iter()
                .map(|c| polylines[&model.coedge(*c).edge].vertices.len() - 1)
                .sum();
            assert_eq!(face.triangles.len(), segments - 2, "{case}");
        }
        // At these radii the strip cuts more triangles from the top cap than it
        // and its two wall cuts add. When an end wall's cap could not join two
        // edges, each wall alone took about 60,000 triangles at the default
        // chord and never met the normal limit.
        assert!(
            mesh.triangles.len() <= plain.triangles.len(),
            "{case}: {} triangles, unrounded {}",
            mesh.triangles.len(),
            plain.triangles.len()
        );

        let (preview, work) = with_operation_budget(APPLICATION_WORK, || {
            preview_mesh_body(
                model,
                built.result,
                &PreviewOptions::from(&options),
                &mut Decider::default(),
            )
        });
        assert!(!work.exhausted, "{case} preview: {work:?}");
        let preview = preview.unwrap_or_else(|e| panic!("{case} preview: {e:?}"));
        assert_eq!(preview.faces.len(), faces.len());
        assert!(preview.faces.iter().all(|f| !f.triangles.is_empty()));
    }
}

#[test]
fn strip_names_survive_radius_edits_and_the_rounded_edge_is_gone() {
    let mut previous = None;
    for radius in ["1", "4", "1"] {
        let built = evaluate(&recipe(&["top"], radius, 3), &mut decider()).unwrap();
        let names = origins(&built);
        if let Some(previous) = &previous {
            assert_eq!(previous, &names);
        }
        previous = Some(names);
    }
    // The rounded edge no longer exists, although both its faces survive.
    let mut again = recipe(&["top"], "1", 3);
    let mut second = again.steps[2].clone();
    second.id = StepId::new("again").unwrap();
    let catallax_recipe::Operation::Fillet { body, .. } = &mut second.op else {
        unreachable!();
    };
    *body = StepId::new("round").unwrap();
    again.steps.push(second);
    again.result = StepId::new("again").unwrap();
    match evaluate(&again, &mut decider()) {
        Err(BuildError::Geometry { step, refusal }) => {
            assert_eq!(step.as_str(), "again");
            assert_eq!(refusal.code, "fillet.lost-edge");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn oversized_slot_rounds_refuse_without_touching_the_retained_input() {
    let base = evaluate_incremental(&recipe(&[], "0", 3), None, &mut decider()).unwrap();
    let bytes = export(base.built());
    let names = origins(base.built());
    for (segments, code) in [
        (
            &["bottom", "right", "top", "left"][..],
            "fillet.unsupported-contact",
        ),
        (&["top"][..], "fillet.unsupported-overflow"),
    ] {
        match evaluate_incremental(
            &recipe(segments, "6", 3),
            Some(base.clone()),
            &mut decider(),
        ) {
            Err(BuildError::Geometry { step, refusal }) => {
                assert_eq!(step.as_str(), "round");
                assert_eq!(refusal.code, code, "{segments:?}");
            }
            other => panic!("{segments:?}: {other:?}"),
        }
        let retained = base.built();
        assert_eq!(retained.model.check(retained.result), vec![]);
        assert_eq!(export(retained), bytes);
        assert_eq!(origins(retained), names);
    }
    assert_volume(base.built(), "4785.398163397448");
}

const LOOP: [&str; 4] = ["bottom", "right", "top", "left"];

fn kinds(built: &Built) -> BTreeMap<String, usize> {
    let mut kinds = BTreeMap::new();
    for face in built.model.faces_of_body(built.result) {
        let surface = built.model.geom().surface(built.model.face(face).surface);
        *kinds.entry(format!("{:?}", surface.kind())).or_insert(0) += 1;
    }
    kinds
}
fn is_blend(built: &Built, face: catallax_topology::FaceId) -> bool {
    built
        .model
        .face_origins(face)
        .iter()
        .any(|o| matches!(o.role, Role::Blend { .. }))
}
fn blend_names(built: &Built) -> BTreeSet<catallax_recipe::Origin> {
    let mut names = BTreeSet::new();
    for face in built.model.faces_of_body(built.result) {
        if is_blend(built, face) {
            names.extend(built.model.face_origins(face).iter().cloned());
        }
    }
    names
}
/// Origin sets of the faces that are not blends, as a set.
fn input_origins(built: &Built) -> BTreeSet<Vec<catallax_recipe::Origin>> {
    built
        .model
        .faces_of_body(built.result)
        .into_iter()
        .filter(|f| !is_blend(built, *f))
        .map(|f| built.model.face_origins(f).to_vec())
        .collect()
}
/// Every edge between two blend pieces is a shared exact meridian circle,
/// used once by a cylinder strip and once by a toroidal piece, with equal
/// unit normals at both ends: a tangent join with no corner face.
fn assert_tangent_meridians(built: &Built) -> usize {
    let model = &built.model;
    let mut owners: BTreeMap<_, Vec<_>> = BTreeMap::new();
    for face in model.faces_of_body(built.result) {
        for lp in &model.face(face).loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                continue;
            };
            for c in cs {
                owners.entry(model.coedge(*c).edge).or_default().push(face);
            }
        }
    }
    let mut dec = decider();
    let mut meridians = 0;
    for (edge, faces) in owners {
        let [a, b] = faces[..] else {
            panic!("edge {edge:?} has {} uses", faces.len());
        };
        if !(is_blend(built, a) && is_blend(built, b)) {
            continue;
        }
        meridians += 1;
        let e = model.edge(edge);
        assert!(matches!(model.geom().curve(e.curve), Curve::Circle(_)));
        let mut surfaces: Vec<_> = [a, b]
            .map(|f| model.geom().surface(model.face(f).surface).kind())
            .into_iter()
            .map(|k| format!("{k:?}"))
            .collect();
        surfaces.sort();
        assert_eq!(surfaces[0], "Cylinder");
        assert!(surfaces[1].starts_with("Torus"), "{surfaces:?}");
        for v in [e.start, e.end].into_iter().flatten() {
            let p = model.geom().position(model.vertex(v).point);
            let [na, nb] = [a, b].map(|f| {
                let face = model.face(f);
                let n = model
                    .geom()
                    .surface(face.surface)
                    .normal_at(&p, &mut dec)
                    .unwrap();
                if face.sense == Sense::Same {
                    n
                } else {
                    n.scaled(&int(-1))
                }
            });
            assert!(zero(na.cross(&nb).norm_squared(), &mut dec));
            assert!(positive(&na.dot(&nb), &mut dec));
        }
    }
    meridians
}

#[test]
fn the_top_boundary_rounds_as_one_tangent_loop_at_every_radius() {
    let base = evaluate(&recipe(&[], "0", 3), &mut decider()).unwrap();
    let inputs = input_origins(&base);
    assert_eq!(inputs.len(), 6);
    let mut explicit = None;
    for (radius, volume) in [
        (1, "4761.789286491508"),
        (2, "4692.167397995440"),
        (4, "4422.113039563441"),
    ] {
        let r = radius.to_string();
        for (segments, version) in [(&LOOP[..], 3), (&LOOP[..], 4), (&["right"][..], 4)] {
            let mut built = evaluate(&recipe(segments, &r, version), &mut decider())
                .unwrap_or_else(|e| panic!("r={radius} {segments:?} v{version}: {e}"));
            assert_eq!(built.model.check(built.result), vec![]);
            assert_eq!(
                kinds(&built),
                BTreeMap::from([
                    ("Cylinder".into(), 4),
                    ("Plane".into(), 4),
                    // Radius 4 leaves the ring torus for its outer sheet.
                    ((if radius == 4 { "TorusSheet" } else { "Torus" }).into(), 2)
                ])
            );
            assert_eq!(assert_tangent_meridians(&built), 4);
            assert_eq!(input_origins(&built), inputs);
            assert_volume(&built, volume);
            let ([x, ..], [tx, ..]) = (
                [format!("{}/10", 10 * (25 - radius) - 1)],
                [format!("{}/10", 10 * (25 - radius) + 5)],
            );
            let side = format!("{}", 10 - radius);
            classify(
                &mut built,
                &[
                    (["0", "49/10", "99/10"], Containment::Outside),
                    (["0", "-49/10", "99/10"], Containment::Outside),
                    (["249/10", "0", "99/10"], Containment::Outside),
                    (["0", "0", "1"], Containment::Inside),
                    (["0", "49/10", "1/10"], Containment::Inside),
                    // At z=10-r the side contact keeps the full end radius 5;
                    // just below the top the capsule radius is about 5-r.
                    (["249/10", "0", side.as_str()], Containment::Inside),
                    (["251/10", "0", side.as_str()], Containment::Outside),
                    ([x.as_str(), "0", "9999/1000"], Containment::Inside),
                    ([tx.as_str(), "0", "9999/1000"], Containment::Outside),
                ],
            );
            let names = blend_names(&built);
            if segments.len() == 1 {
                // The continued loop keeps the one selected arc's name.
                assert_eq!(names.len(), 1);
                let arc = names.iter().next().unwrap();
                assert!(
                    explicit
                        .as_ref()
                        .is_some_and(|e: &BTreeSet<_>| e.contains(arc))
                );
            } else {
                assert_eq!(names.len(), 4);
                if let Some(previous) = &explicit {
                    assert_eq!(previous, &names, "radius edits keep the four names");
                }
                explicit = Some(names);
            }
        }
    }
}

#[test]
fn permuted_and_repeated_references_build_the_same_loop() {
    let reference = evaluate(&recipe(&LOOP, "4", 4), &mut decider()).unwrap();
    let permuted = evaluate(
        &recipe(&["left", "top", "right", "bottom", "top"], "4", 4),
        &mut decider(),
    )
    .unwrap();
    assert_eq!(origins(&permuted), origins(&reference));
    assert_eq!(export(&permuted), export(&reference));
}

#[test]
fn two_straight_strips_with_meeting_carriers_round_disjoint_sectors() {
    for (radius, volume) in [(1, "4767.090708834183"), (4, "4476.080429032569")] {
        let built = evaluate(
            &recipe(&["bottom", "top"], &radius.to_string(), 3),
            &mut decider(),
        )
        .unwrap_or_else(|e| panic!("r={radius}: {e}"));
        assert_eq!(built.model.check(built.result), vec![]);
        assert_eq!(
            kinds(&built),
            BTreeMap::from([("Cylinder".into(), 4), ("Plane".into(), 4)])
        );
        assert_eq!(blend_names(&built).len(), 2);
        assert_volume(&built, volume);
    }
}
