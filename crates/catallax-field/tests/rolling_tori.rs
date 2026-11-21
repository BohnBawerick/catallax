//! 4.2a step 6: torus pairs. The corpus rows `torus-slice` and
//! `torus-intersection` round the edge where a plane normal to a torus's
//! axis cuts it, and a cut below the torus's centre plane does the same with
//! the offset plane below the core circle.
//!
//! The field derives each edge's ball path itself: the plane's offset meets
//! the torus offset by the ball radius (the torus of tube radius `a - r`) in
//! a circle about the axis, a ring piece. Where the offset plane touches the
//! offset tube's top, as in `torus-intersection`, both edges share one
//! centre circle and their fans meet in the cylinder over the core circle.
//! Its membership must match each witness: the corpus row's and points
//! derived here, each checked against the B-rep first. Every display vertex
//! of the B-rep's blend faces must lie within the field's guaranteed distance
//! range of a piece's face at 1 nm, with the same records at 1, 2 and 4
//! workers. The corpus row's boundary witnesses lie on a piece's region
//! boundary with a distance range from 0.
use catallax_features::evaluate;
use catallax_field::rolling_region::piece::Piece;
use catallax_field::rolling_region::{Batch, Range, Relation, distance_batch};
use catallax_field::{Location, Plan, Workspace};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, OperationAllowance, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Point3, RPoint3};
use serde_json::{Value, json};

fn ratio(n: i64, d: i64) -> Rational {
    Rational::from_ratio(n, d).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn parse(value: &Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}
fn number(value: &Value) -> Rational {
    Rational::parse(value.as_str().unwrap()).unwrap()
}
fn corpus_row(id: &str) -> Value {
    include_str!("../../../parts/fillets/cases.jsonl")
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == id)
        .unwrap()
}

/// The corpus row's recipe with the block's bottom moved to `cut` and the
/// fillet's radius set to `radius`.
fn cut_torus(cut: &str, radius: &str) -> Recipe {
    let mut value = corpus_row("torus-slice")["recipe"].clone();
    let steps = value["steps"].as_array_mut().unwrap();
    steps[1]["pose"]["translation"] = json!(["-20", "-20", cut]);
    steps[3]["radius"] = json!(radius);
    value["name"] = json!(format!("torus cut at {cut}, radius {radius}"));
    parse(&value)
}

type Points = [([i64; 2], [i64; 2], [i64; 2])];
fn points(list: &Points) -> Vec<[Rational; 3]> {
    list.iter()
        .map(|(x, y, z)| [ratio(x[0], x[1]), ratio(y[0], y[1]), ratio(z[0], z[1])])
        .collect()
}

/// One torus case: its recipe, the ring pieces the field compiles, points
/// inside and outside, and points on the boundary.
struct Case {
    name: String,
    recipe: Recipe,
    pieces: usize,
    inside: Vec<[Rational; 3]>,
    outside: Vec<[Rational; 3]>,
    on: Vec<[Rational; 3]>,
}
impl Case {
    fn new(name: &str, recipe: Recipe, inside: &Points, outside: &Points) -> Self {
        Self {
            name: name.to_string(),
            recipe,
            pieces: 2,
            inside: points(inside),
            outside: points(outside),
            on: Vec::new(),
        }
    }
    /// A corpus row with its witnesses added to the derived points.
    fn corpus(id: &str, inside: &Points, outside: &Points) -> Self {
        let row = corpus_row(id);
        let mut case = Self::new(id, parse(&row["recipe"]), inside, outside);
        for witness in row["witnesses"].as_array().unwrap() {
            let p: [Rational; 3] = std::array::from_fn(|i| number(&witness["point"][i]));
            match witness["state"].as_str().unwrap() {
                "TopAbs_IN" => case.inside.push(p),
                "TopAbs_OUT" => case.outside.push(p),
                "TopAbs_ON" => case.on.push(p),
                other => panic!("{id}: witness state {other}"),
            }
        }
        case
    }
}

fn cases() -> Vec<Case> {
    vec![
        // Cut at z = 6/5, r = 1/5: centre circles of radius 6 +- 4 sqrt(2) / 5
        // at z = 7/5. Derived: inside the outer ball and in another meridian
        // (3/5, 4/5), cut outside it there and opposite, inside the inner
        // ball.
        Case::corpus(
            "torus-slice",
            &[
                ([71, 10], [0, 1], [7, 5]),
                ([0, 1], [71, 10], [7, 5]),
                ([213, 50], [142, 25], [7, 5]),
                ([49, 10], [0, 1], [7, 5]),
            ],
            &[
                ([111, 25], [148, 25], [5, 4]),
                ([-37, 5], [0, 1], [5, 4]),
                ([0, 1], [-23, 5], [5, 4]),
            ],
        ),
        // Cut at z = 1, r = 1/2: one centre circle of radius 6 at z = 3/2.
        // Derived: inside the tube near its bottom and to either side, in
        // another meridian; cut outside it on both sides of the core
        // circle's cylinder.
        Case::corpus(
            "torus-intersection",
            &[
                ([6, 1], [0, 1], [6, 5]),
                ([31, 5], [0, 1], [9, 5]),
                ([18, 5], [24, 5], [3, 2]),
            ],
            &[
                ([27, 4], [0, 1], [3, 2]),
                ([0, 1], [-21, 4], [3, 2]),
                ([-15, 2], [0, 1], [5, 4]),
            ],
        ),
        // Cut at z = -1, below the core circle's plane, r = 1/5: the offset
        // plane at z = -4/5 meets the offset tube in circles of radius
        // 6 +- sqrt(65) / 5, and the edges are at 6 +- sqrt(3). Derived:
        // inside each ball, cut outside them, and uncut material.
        Case::new(
            "torus cut below its centre",
            cut_torus("-1", "1/5"),
            &[
                ([38, 5], [0, 1], [-4, 5]),
                ([22, 5], [0, 1], [-4, 5]),
                ([6, 1], [0, 1], [0, 1]),
            ],
            &[
                ([771, 100], [0, 1], [-99, 100]),
                ([0, 1], [-429, 100], [-99, 100]),
            ],
        ),
    ]
}

fn brep_location(built: &mut catallax_features::Built, p: &[Rational; 3]) -> Location {
    use catallax_boolean::{Containment, classify_point};
    use catallax_numbers::Outcome;
    let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
    match classify_point(&mut built.model, built.result, &point, &mut decider()) {
        Outcome::Proven(Containment::Inside) => Location::Inside,
        Outcome::Proven(Containment::Outside) => Location::Outside,
        Outcome::ProvenCoincident(_) => Location::Boundary,
        other => panic!("the B-rep does not settle {p:?}: {other:?}"),
    }
}

fn location(plan: &Plan, p: &[Rational; 3]) -> Location {
    let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
    plan.proven_checked(&point, &mut Workspace::default(), &mut decider())
        .unwrap()
        .location
}

fn real_point(p: &[Rational; 3]) -> Point3 {
    Point3::new(
        p[0].clone().into(),
        p[1].clone().into(),
        p[2].clone().into(),
    )
}

#[test]
fn the_field_rounds_each_torus_edge_and_agrees_on_membership() {
    let mut failures = Vec::new();
    for case in cases() {
        let mut dec = decider();
        let plan = match Plan::compile_checked(&case.recipe, &mut dec) {
            Ok(plan) => plan,
            Err(e) => {
                failures.push(format!("{}: {e:?}", case.name));
                continue;
            }
        };
        let pieces = plan.rolling_pieces().count();
        if pieces != case.pieces {
            failures.push(format!("{}: {pieces} rolling pieces", case.name));
        }
        let mut built = evaluate(&case.recipe, &mut decider()).unwrap();
        for (points, want) in [
            (&case.inside, Location::Inside),
            (&case.outside, Location::Outside),
        ] {
            for p in points {
                assert_eq!(
                    brep_location(&mut built, p),
                    want,
                    "{}: the B-rep at {p:?}",
                    case.name
                );
                let found = location(&plan, p);
                if found != want {
                    failures.push(format!("{} at {p:?}: {found:?}", case.name));
                }
            }
        }
        // A boundary witness lies on the region boundary of one piece at
        // least, at distance 0 from its face.
        let width = ratio(1, 1_000_000_000_000);
        for p in &case.on {
            assert_eq!(
                brep_location(&mut built, p),
                Location::Boundary,
                "{}: the B-rep at {p:?}",
                case.name
            );
            let point = real_point(p);
            let mut on = false;
            for piece in plan.rolling_pieces() {
                if piece.classify(&point, &mut decider()).unwrap() == Relation::Boundary {
                    let range = piece.distance(&point, &width, &mut decider()).unwrap();
                    on = on || range.lower == Rational::zero();
                }
            }
            if !on {
                failures.push(format!("{} at {p:?}: on no piece's face", case.name));
            }
        }
        eprintln!(
            "{}: built with {} work, {} witnesses",
            case.name,
            dec.stats().work,
            case.inside.len() + case.outside.len() + case.on.len()
        );
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// The display vertices of the B-rep's blend faces made by `step`. Every
/// such face here is an exact torus, so the part's tolerance at each vertex
/// is 1 nm.
fn blend_vertices(recipe: &Recipe, step: &str, chord: &str) -> Vec<Point3> {
    let built = evaluate(recipe, &mut decider()).unwrap();
    let mesh = mesh_body(
        &built.model,
        built.result,
        &MeshOptions {
            chord_tolerance: Rational::parse(chord).unwrap(),
            ..MeshOptions::default()
        },
        &mut decider(),
    )
    .unwrap()
    .value;
    let mut indices = std::collections::BTreeSet::new();
    for face in &mesh.faces {
        if !built
            .model
            .face_origins(face.face)
            .iter()
            .any(|o| o.step.as_str() == step && matches!(o.role, Role::Blend { .. }))
        {
            continue;
        }
        for triangle in &mesh.triangles[face.triangles.clone()] {
            indices.extend(triangle.vertices);
        }
    }
    indices
        .into_iter()
        .map(|v| {
            let [x, y, z] = mesh.positions[v].map(|c| Rational::from_f64_exact(c).unwrap());
            Point3::new(x.into(), y.into(), z.into())
        })
        .collect()
}

/// Upper bounds on the distance from each point to the nearest of the
/// pieces' faces, from guaranteed ranges at most `width` wide, and the most
/// work one query took. Every 128th point is also queried at 1, 2 and 4
/// workers, whose records must agree.
fn nearest_witnesses(
    name: &str,
    pieces: &[&Piece],
    points: &[Point3],
    width: &Rational,
) -> (Vec<Option<Rational>>, u64) {
    const CAP: u64 = 20_000;
    let allowance = |n: usize| OperationAllowance::reserve(n as u64 * CAP).unwrap();
    let mut nearest = vec![None::<Rational>; points.len()];
    let mut most = 0;
    for &piece in pieces {
        let index: Vec<usize> = (0..points.len())
            .filter(|&i| nearest[i].as_ref().is_none_or(|u| u > width))
            .collect();
        if index.is_empty() {
            break;
        }
        let mine: Vec<Point3> = index.iter().map(|&i| points[i].clone()).collect();
        let sample: Vec<Point3> = mine.iter().step_by(128).cloned().collect();
        let all = distance_batch(
            piece,
            &mine,
            width,
            None,
            allowance(mine.len()),
            &mut decider(),
        )
        .unwrap();
        let mut batches: Vec<Batch<Range>> = Vec::new();
        for workers in [1, 2, 4] {
            batches.push(
                distance_batch(
                    piece,
                    &sample,
                    width,
                    Some(workers),
                    allowance(sample.len()),
                    &mut decider(),
                )
                .unwrap(),
            );
        }
        assert_eq!(batches[0].records, batches[1].records, "{name}");
        assert_eq!(batches[0].records, batches[2].records, "{name}");
        for (&i, record) in index.iter().zip(&all.records) {
            most = most.max(record.work);
            let Ok(range) = &record.result else {
                assert!(record.exhausted, "{name}: {record:?}");
                continue;
            };
            if let Some(upper) = &range.upper
                && nearest[i].as_ref().is_none_or(|best| upper < best)
            {
                nearest[i] = Some(upper.clone());
            }
        }
    }
    (nearest, most)
}

/// Every B-rep display vertex of the rounds' blend faces lies within the
/// field's guaranteed range of a ring piece's face, at 1 nm.
#[test]
fn torus_blend_vertices_lie_on_the_field_face() {
    let tolerance = ratio(1, 1_000_000);
    for case in cases() {
        let start = std::time::Instant::now();
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        let points = blend_vertices(&case.recipe, "round", "1/50");
        assert!(!points.is_empty(), "{}: no blend vertices", case.name);
        let pieces: Vec<&Piece> = plan.rolling_pieces().collect();
        let (nearest, most) = nearest_witnesses(&case.name, &pieces, &points, &tolerance);
        let mut worst = Rational::zero();
        for found in &nearest {
            let found = found
                .clone()
                .unwrap_or_else(|| panic!("{}: no witness", case.name));
            assert!(found <= tolerance, "{}: {found:?}", case.name);
            worst = worst.max(found);
        }
        eprintln!(
            "{}: {} blend vertices in {:?}, worst field distance {}, at most {most} work a query",
            case.name,
            nearest.len(),
            start.elapsed(),
            worst.to_f64_nearest()
        );
    }
}
