//! 4.2a step 2: a plane with a cylinder in any pose, a cone or a sphere, and
//! straight runs between supports parallel to one line.
//!
//! The field compiles these fillets from its own path pieces. Its membership
//! must match the witnesses of the fillet corpus (checked against OpenCASCADE)
//! and of the step's test files (checked against the B-rep). Where the B-rep
//! builds, every display vertex of a B-rep blend face must lie within the
//! field's guaranteed distance range of the field's own face, at the 1 nm
//! own-part tolerance, with the same records at 1, 2 and 4 workers.
use catallax_features::evaluate;
use catallax_field::rolling_region::piece::Piece;
use catallax_field::rolling_region::{Batch, Range, distance_batch};
use catallax_field::{Location, Plan, Workspace};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{BudgetPolicy, Decider, OperationAllowance, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Point3, RPoint3};
use serde_json::{Value, json};

/// A half turn about y of `2 atan(7/24)`: its z axis is `(336, 0, 527) / 625`.
const TILT: [&str; 4] = ["24/25", "0", "7/25", "0"];
const LEVEL: [&str; 4] = ["1", "0", "0", "0"];

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000_000,
        ..Default::default()
    })
}
fn parse(value: Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}

/// One recipe with points whose membership the corpus or a B-rep test fixed.
struct Case {
    name: String,
    recipe: Recipe,
    inside: Vec<[Rational; 3]>,
    outside: Vec<[Rational; 3]>,
    /// The B-rep builds it and no earlier case has the same recipe, so its
    /// display vertices are checked.
    brep: bool,
    /// A corpus row, meshed at the default display chord of 0.01 mm; the
    /// test files' recipes are meshed at 0.05 mm to keep the test short.
    corpus: bool,
    /// The elementary route of step 1 keeps it: no rolling piece.
    elementary: bool,
}
fn case(name: &str, recipe: Recipe, inside: &[[&str; 3]], outside: &[[&str; 3]]) -> Case {
    Case {
        name: name.into(),
        recipe,
        inside: inside.iter().map(|p| p.map(q)).collect(),
        outside: outside.iter().map(|p| p.map(q)).collect(),
        brep: true,
        corpus: false,
        elementary: false,
    }
}

fn corpus(id: &str) -> Case {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../parts/fillets/cases.jsonl"
    ))
    .unwrap();
    let row: Value = text
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).unwrap())
        .find(|row| row["id"] == id)
        .unwrap();
    let mut inside = Vec::new();
    let mut outside = Vec::new();
    for witness in row["witnesses"].as_array().into_iter().flatten() {
        let point = std::array::from_fn(|i| q(witness["point"][i].as_str().unwrap()));
        match witness["state"].as_str().unwrap() {
            "TopAbs_IN" => inside.push(point),
            "TopAbs_OUT" => outside.push(point),
            other => panic!("{other}"),
        }
    }
    Case {
        name: id.into(),
        recipe: parse(row["recipe"].clone()),
        inside,
        outside,
        brep: row["expected"]["status"] == "built",
        corpus: true,
        elementary: false,
    }
}

// The oblique file: a cylinder cut, bossed or bored by a tilted slab.
fn slab(at: &str, distance: &str) -> [Value; 2] {
    [
        json!({"id":"cut-section","op":"sketch",
            "plane":{"rotation":TILT,"translation":["0","0",at]},
            "loops":[[{"kind":"line","id":"s","to":["30","-30"]},
                {"kind":"line","id":"e","to":["30","30"]},
                {"kind":"line","id":"n","to":["-30","30"]},
                {"kind":"line","id":"w","to":["-30","-30"]}]]}),
        json!({"id":"cut","op":"extrude","profile":"cut-section","distance":distance}),
    ]
}
fn oblique(kind: &str, radius: &str) -> Recipe {
    let wall = |height: &str, z: &str| {
        json!({"id":"wall","op":"cylinder","radius":"5","height":height,
            "pose":{"rotation":LEVEL,"translation":["0","0",z]}})
    };
    let (mut steps, boolean, role) = match kind {
        "rim" => (vec![wall("16", "0")], ("subtract", "wall", "cut"), "start"),
        "boss" => (vec![wall("22", "-10")], ("union", "wall", "cut"), "end"),
        _ => (vec![wall("100", "-60")], ("subtract", "cut", "wall"), "end"),
    };
    steps.extend(if kind == "rim" {
        slab("8", "30")
    } else {
        slab("-25", "527/25")
    });
    steps.push(json!({"id":"joined","op":"boolean","kind":boolean.0,
        "target":boolean.1,"tool":boolean.2}));
    steps.push(
        json!({"id":"round","op":"fillet","body":"joined","radius":radius,
        "edges":[{"between":[{"step":"wall","role":{"side":{"segment":"wall"}}},
            {"step":"cut","role":role}]}]}),
    );
    parse(
        json!({"format":"catallax-recipe/3","name":"oblique round","steps":steps,
        "result":"round"}),
    )
}

// The oblique-cone file: the corpus frustum cut by the tilted plane, mirrored
// in x = 0 after the round when asked.
fn cut_cone(mirror: bool) -> Recipe {
    let mut steps = vec![
        json!({"id":"wall","op":"cone","base_radius":"10","top_radius":"4","height":"24",
            "pose":{"rotation":LEVEL,"translation":["0","0","0"]}}),
    ];
    steps.extend(slab("8", "30"));
    steps
        .push(json!({"id":"joined","op":"boolean","kind":"subtract","target":"wall","tool":"cut"}));
    steps.push(
        json!({"id":"round","op":"fillet","body":"joined","radius":"1/2",
        "edges":[{"between":[{"step":"wall","role":{"side":{"segment":"wall"}}},
            {"step":"cut","role":"start"}]}]}),
    );
    let mut result = "round";
    if mirror {
        steps.push(json!({"id":"reflected","op":"mirror","body":"round",
            "plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}));
        result = "reflected";
    }
    parse(
        json!({"format":"catallax-recipe/4","name":"oblique cone round","steps":steps,
        "result":result}),
    )
}

// The generator file: cylinders of radius 5 and height 12 along z.
fn cylinder(id: &str, x: &str) -> Value {
    json!({"id":id,"op":"cylinder","radius":"5","height":"12",
        "pose":{"rotation":LEVEL,"translation":[x,"0","0"]}})
}
fn generator(steps: Vec<Value>, radius: &str, edge: Value) -> Recipe {
    let mut steps = steps;
    let body = steps.last().unwrap()["id"].clone();
    steps.push(json!({"id":"round","op":"fillet","body":body,"radius":radius,"edges":[edge]}));
    parse(
        json!({"format":"catallax-recipe/3","name":"generator round","steps":steps,
        "result":"round"}),
    )
}
fn walls() -> Value {
    json!({"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
        {"step":"b","role":{"side":{"segment":"wall"}}}]})
}
fn pair(kind: &str) -> Vec<Value> {
    vec![
        cylinder("a", "-3"),
        cylinder("b", "3"),
        json!({"id":"joined","op":"boolean","kind":kind,"target":"a","tool":"b"}),
    ]
}
fn flat() -> Value {
    json!({"between":[{"step":"a","role":{"side":{"segment":"wall"}}},
        {"step":"block","role":{"side":{"segment":"n"}}}]})
}
/// The cylinder at the origin intersected with a block whose north side is a
/// flat at y = 3, spanning x from -12 to `east`.
fn flatted(east: &str) -> Vec<Value> {
    vec![
        cylinder("a", "0"),
        json!({"id":"block-section","op":"sketch",
            "plane":{"rotation":LEVEL,"translation":["0","0","0"]},
            "loops":[[{"kind":"line","id":"s","to":[east,"-10"]},
                {"kind":"line","id":"e","to":[east,"3"]},
                {"kind":"line","id":"n","to":["-12","3"]},
                {"kind":"line","id":"w","to":["-12","-10"]}]]}),
        json!({"id":"block","op":"extrude","profile":"block-section","distance":"12"}),
        json!({"id":"joined","op":"boolean","kind":"intersect","target":"a","tool":"block"}),
    ]
}

fn cases() -> Vec<Case> {
    let mut cases: Vec<Case> = [
        "cone-top",
        "cone-base",
        "sphere-dome",
        "oblique-cylinder",
        "narrow-cone",
        "oblique-cone",
    ]
    .into_iter()
    .map(corpus)
    .collect();
    // Radius 1 is the corpus row `oblique-cylinder`, whose vertices are
    // already checked.
    for radius in ["1", "2"] {
        let mut rim = case(
            &format!("oblique rim r={radius}"),
            oblique("rim", radius),
            &[["0", "0", "4"], ["0", "4", "6"], ["-3", "-3", "7"]],
            &[["0", "0", "9"]],
        );
        rim.brep = radius != "1";
        cases.push(rim);
    }
    let mut blend = case(
        "oblique rim r=1, at the blend",
        oblique("rim", "1"),
        &[["24/5", "0", "47/10"], ["-43/10", "0", "991/100"]],
        &[["99/20", "0", "239/50"], ["-24/5", "0", "541/50"]],
    );
    blend.brep = false;
    cases.push(blend);
    cases.push(case(
        "oblique boss",
        oblique("boss", "1"),
        &[["101/20", "0", "-79/25"], ["-26/5", "0", "71/20"]],
        &[["26/5", "0", "-77/25"], ["-57/10", "0", "223/50"]],
    ));
    cases.push(case(
        "oblique bore",
        oblique("bore", "1"),
        &[["57/10", "0", "-223/50"], ["-26/5", "0", "77/25"]],
        &[["26/5", "0", "-71/20"], ["-101/20", "0", "79/25"]],
    ));
    let cone_inside = [
        ["0", "0", "4"],
        ["0", "5", "5"],
        ["-5", "-3", "6"],
        ["47/5", "0", "37/20"],
        ["-673/100", "0", "1197/100"],
    ];
    let cone_outside = [
        ["0", "0", "9"],
        ["10", "0", "1"],
        ["0", "0", "-1/10"],
        ["4749/500", "0", "1919/1000"],
        ["-6821/1000", "0", "12201/1000"],
    ];
    cases.push(case(
        "cut cone",
        cut_cone(false),
        &cone_inside,
        &cone_outside,
    ));
    let mirrored = |points: &[[&'static str; 3]]| -> Vec<[String; 3]> {
        points
            .iter()
            .map(|[x, y, z]| {
                let x = (-q(x)).to_string();
                [x, y.to_string(), z.to_string()]
            })
            .collect()
    };
    let (inside, outside) = (mirrored(&cone_inside), mirrored(&cone_outside));
    let refs = |v: &[[String; 3]]| -> Vec<[Rational; 3]> {
        v.iter().map(|p| p.each_ref().map(|s| q(s))).collect()
    };
    cases.push(Case {
        name: "mirrored cut cone".into(),
        recipe: cut_cone(true),
        inside: refs(&inside),
        outside: refs(&outside),
        brep: true,
        corpus: false,
        elementary: false,
    });
    cases.push(case(
        "parallel cylinders, intersect",
        generator(pair("intersect"), "5/4", walls()),
        &[["0", "17/5", "6"], ["0", "-17/5", "6"]],
        &[["0", "18/5", "6"], ["0", "-18/5", "6"]],
    ));
    cases.push(case(
        "parallel cylinders, union",
        generator(pair("union"), "11/8", walls()),
        &[["0", "41/10", "6"], ["0", "-41/10", "6"]],
        &[["0", "22/5", "6"], ["0", "-22/5", "6"]],
    ));
    cases.push(case(
        "parallel cylinders, crescent",
        generator(pair("subtract"), "1", walls()),
        &[["-1", "4", "6"]],
        &[],
    ));
    cases.push(case(
        "a flat's two generators",
        generator(flatted("12"), "7/4", flat()),
        &[["18/5", "13/5", "6"], ["-18/5", "13/5", "6"]],
        &[["39/10", "29/10", "6"], ["-39/10", "29/10", "6"]],
    ));
    let mut rational = case(
        "one generator, rational centre",
        generator(flatted("0"), "7/4", flat()),
        &[["-18/5", "13/5", "6"], ["-3", "1", "6"], ["-1", "1", "6"]],
        &[["-39/10", "29/10", "6"], ["1", "1", "6"]],
    );
    rational.elementary = true;
    cases.push(rational);
    cases.push(case(
        "one generator, irrational centre",
        generator(flatted("0"), "3/4", flat()),
        &[["-3", "1", "6"], ["-1", "1", "6"]],
        &[["1", "1", "6"]],
    ));
    cases
}

fn location(plan: &Plan, p: &[Rational; 3]) -> Location {
    let point = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
    plan.proven_checked(&point, &mut Workspace::default(), &mut decider())
        .unwrap()
        .location
}

#[test]
fn the_field_builds_every_step_two_fillet_and_agrees_on_membership() {
    for case in cases() {
        let plan = Plan::compile_checked(&case.recipe, &mut decider())
            .unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
        assert_eq!(
            plan.rolling_pieces().next().is_none(),
            case.elementary,
            "{}: rolling pieces",
            case.name
        );
        for p in &case.inside {
            assert_eq!(
                location(&plan, p),
                Location::Inside,
                "{} at {p:?}",
                case.name
            );
        }
        for p in &case.outside {
            assert_eq!(
                location(&plan, p),
                Location::Outside,
                "{} at {p:?}",
                case.name
            );
        }
    }
}

/// The cases whose B-rep blend vertices are checked, one test per group so
/// they run side by side.
const VERTEX_GROUPS: &[&[&str]] = &[
    &["cone-top", "cone-base", "narrow-cone"],
    &["sphere-dome"],
    &["oblique-cylinder"],
    &["oblique rim r=2"],
    &["oblique boss"],
    &["oblique bore"],
    &["cut cone"],
    &["mirrored cut cone"],
    &[
        "parallel cylinders, intersect",
        "parallel cylinders, union",
        "parallel cylinders, crescent",
        "a flat's two generators",
        "one generator, irrational centre",
    ],
];
macro_rules! vertex_tests {
    ($($test:ident = $group:literal;)+) => {
        $(
            #[test]
            fn $test() {
                blend_vertices_lie_on_the_field_face(VERTEX_GROUPS[$group]);
            }
        )+
    };
}
vertex_tests! {
    cone_blend_vertices_lie_on_the_field_face = 0;
    sphere_blend_vertices_lie_on_the_field_face = 1;
    oblique_cylinder_blend_vertices_lie_on_the_field_face = 2;
    oblique_rim_blend_vertices_lie_on_the_field_face = 3;
    oblique_boss_blend_vertices_lie_on_the_field_face = 4;
    oblique_bore_blend_vertices_lie_on_the_field_face = 5;
    cut_cone_blend_vertices_lie_on_the_field_face = 6;
    mirrored_cut_cone_blend_vertices_lie_on_the_field_face = 7;
    generator_blend_vertices_lie_on_the_field_face = 8;
}

#[test]
fn every_built_case_has_its_vertices_checked() {
    let listed: std::collections::BTreeSet<String> = VERTEX_GROUPS
        .iter()
        .flat_map(|group| group.iter().map(|name| name.to_string()))
        .collect();
    let wanted: std::collections::BTreeSet<String> = cases()
        .into_iter()
        .filter(|c| c.brep && !c.elementary)
        .map(|c| c.name)
        .collect();
    assert_eq!(listed, wanted);
}

/// Every B-rep display vertex of a blend face lies within the field's
/// guaranteed range of a field face, at the 1 nm own-part tolerance.
fn blend_vertices_lie_on_the_field_face(names: &[&str]) {
    // The vertex lies within the field's bound once a face witness within the
    // tolerance is proven, so the range need be no narrower than that.
    let tolerance = q("1/1000000");
    let width = tolerance.clone();
    for case in cases()
        .into_iter()
        .filter(|c| names.contains(&c.name.as_str()))
    {
        assert!(case.brep && !case.elementary, "{}", case.name);
        let start = std::time::Instant::now();
        let plan = Plan::compile_checked(&case.recipe, &mut decider()).unwrap();
        let compiled = start.elapsed();
        let built = evaluate(&case.recipe, &mut decider())
            .unwrap_or_else(|e| panic!("{}: {e:?}", case.name));
        let brep = start.elapsed();
        let mesh = mesh_body(
            &built.model,
            built.result,
            &MeshOptions {
                chord_tolerance: if case.corpus {
                    MeshOptions::default().chord_tolerance
                } else {
                    q("1/20")
                },
                ..MeshOptions::default()
            },
            &mut decider(),
        )
        .unwrap()
        .value;
        let mut indices = std::collections::BTreeSet::new();
        for face in &mesh.faces {
            if !built.model.face_origins(face.face).iter().any(|o| {
                let mut role = &o.role;
                while let Role::Copy { of, .. } = role {
                    role = &of.role;
                }
                matches!(role, Role::Blend { .. })
            }) {
                continue;
            }
            for triangle in &mesh.triangles[face.triangles.clone()] {
                indices.extend(triangle.vertices);
            }
        }
        let points: Vec<Point3> = indices
            .into_iter()
            .map(|v| {
                let [x, y, z] = mesh.positions[v].map(|c| Rational::from_f64_exact(c).unwrap());
                Point3::new(x.into(), y.into(), z.into())
            })
            .collect();
        assert!(!points.is_empty(), "{}: no blend vertices", case.name);
        let pieces: Vec<_> = plan.rolling_pieces().collect();
        let mut nearest = vec![None::<Rational>; points.len()];
        // Every vertex once on all workers; every 128th also at 1, 2 and 4
        // workers, whose records and work must agree.
        // A range at a far piece needs deep refinement, so each query has a
        // work cap; an exhausted query only yields no witness, never a wrong
        // one. The matching piece settled every vertex within 7,781 units on
        // circles and 574 on ellipses when this was written.
        const CAP: u64 = 20_000;
        let allowance = |n: usize| OperationAllowance::reserve(n as u64 * CAP).unwrap();
        let mut used = 0;
        // A piece whose box lies farther than the tolerance from a vertex
        // cannot hold its witness; a mirror leaves the original piece in the
        // plan beside its copy.
        let near = |piece: &Piece, p: &Point3| {
            let Some(bounds) = piece.bounds() else {
                return true;
            };
            let [x, y, z] = [&p.x, &p.y, &p.z].map(|c| c.rational_ref().unwrap().clone());
            let gap = [x, y, z]
                .iter()
                .zip(bounds)
                .map(|(c, [lo, hi])| {
                    let g = if c < lo {
                        lo - c
                    } else if c > hi {
                        c - hi
                    } else {
                        Rational::zero()
                    };
                    &g * &g
                })
                .fold(Rational::zero(), |sum, g| sum + g);
            gap <= &tolerance * &tolerance
        };
        // Newest piece first: a mirror or pattern copy is compiled after its
        // original, and a vertex leaves the search once it has a witness.
        for piece in pieces.iter().rev() {
            let index: Vec<usize> = (0..points.len())
                .filter(|&i| {
                    nearest[i].as_ref().is_none_or(|u| u > &tolerance) && near(piece, &points[i])
                })
                .collect();
            if index.is_empty() {
                continue;
            }
            let mine: Vec<Point3> = index.iter().map(|&i| points[i].clone()).collect();
            let sample: Vec<Point3> = mine.iter().step_by(128).cloned().collect();
            let all = distance_batch(
                *piece,
                &mine,
                &width,
                None,
                allowance(mine.len()),
                &mut decider(),
            )
            .unwrap();
            let mut batches: Vec<Batch<Range>> = Vec::new();
            for workers in [1, 2, 4] {
                batches.push(
                    distance_batch(
                        *piece,
                        &sample,
                        &width,
                        Some(workers),
                        allowance(sample.len()),
                        &mut decider(),
                    )
                    .unwrap(),
                );
            }
            assert_eq!(batches[0].records, batches[1].records, "{}", case.name);
            assert_eq!(batches[0].records, batches[2].records, "{}", case.name);
            for (&i, record) in index.iter().zip(&all.records) {
                let Ok(range) = &record.result else {
                    assert!(record.exhausted, "{}: {record:?}", case.name);
                    continue;
                };
                if range.upper.as_ref().is_some_and(|u| u <= &tolerance) {
                    used = used.max(record.work);
                }
                if let Some(upper) = &range.upper
                    && nearest[i].as_ref().is_none_or(|best| upper < best)
                {
                    nearest[i] = Some(upper.clone());
                }
            }
        }
        let queried = start.elapsed();
        eprintln!(
            "{}: field compile {compiled:?}, B-rep and mesh {:?}, queries {:?}, \
             most work for a witness {used}",
            case.name,
            brep - compiled,
            queried - brep
        );
        let worst = nearest
            .iter()
            .map(|u| {
                u.clone()
                    .unwrap_or_else(|| panic!("{}: no witness", case.name))
            })
            .max()
            .unwrap();
        eprintln!(
            "{}: {} blend vertices, worst field distance {}",
            case.name,
            points.len(),
            worst.to_f64_nearest()
        );
        assert!(worst <= tolerance, "{}: {worst:?}", case.name);
    }
}
