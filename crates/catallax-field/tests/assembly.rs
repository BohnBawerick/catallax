//! Analytic controls and public complete recipes; no saved-answer corpus.
use catallax_field::{
    Plan,
    assembly::{self, Mass, Options, PlacedPart, Relation},
    mass,
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::{Pose, RVector3, Rotation};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn placed(plan: &Plan, xyz: [&str; 3], rotation: [&str; 4]) -> PlacedPart {
    let [w, rx, ry, rz] = rotation.map(q);
    let [x, y, z] = xyz.map(q);
    let pose = Pose::from_rational(
        Rotation::from_quaternion(w, rx, ry, rz).unwrap(),
        RVector3::new(x, y, z),
    );
    PlacedPart::new(plan, &pose, &mut Decider::default()).unwrap()
}
fn block() -> Plan {
    Plan::compile(&Recipe::parse(r#"{"format":"catallax-recipe/3","name":"box","steps":[{"id":"box","op":"block","size":["1","2","3"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"box"}"#).unwrap()).unwrap()
}
fn options() -> Options {
    Options {
        gap_width: q("1/1000000"),
        max_pairs: 64,
    }
}
fn compare(a: &PlacedPart, b: &PlacedPart) -> Relation {
    assembly::relation(a, b, &options(), &mut Decider::default())
}
fn gap(result: Relation, lo: &str, hi: &str) {
    let Relation::Clear { gap } = result else {
        panic!("expected clear: {result:?}")
    };
    assert!(gap.lo <= q(lo) && gap.hi >= q(hi));
    assert!(gap.lo > q("0") && gap.width() <= options().gap_width);
}
#[test]
fn boxes_clash_touch_clear_and_exact_diagonal_gap() {
    let plan = block();
    let a = placed(&plan, ["0", "0", "0"], ["1", "0", "0", "0"]);
    let b = placed(&plan, ["4", "6", "0"], ["1", "0", "0", "0"]);
    gap(compare(&a, &b), "5", "5");
    gap(compare(&b, &a), "5", "5");
    assert!(matches!(compare(&a, &a), Relation::Clash { .. }));
    let touch = placed(&plan, ["1", "0", "0"], ["1", "0", "0", "0"]);
    assert!(matches!(compare(&a,&touch),Relation::Unresolved(r) if r.code == "assembly.contact"));
    let tiny = placed(&plan, ["1.000000000001", "0", "0"], ["1", "0", "0", "0"]);
    gap(compare(&a, &tiny), "0.000000000001", "0.000000000001");
}
#[test]
fn common_rigid_transform_preserves_gap_and_moves_centroid() {
    let a = placed(&block(), ["10", "20", "30"], ["1", "0", "0", "1"]);
    let b = placed(&block(), ["4", "24", "30"], ["1", "0", "0", "1"]);
    gap(compare(&a, &b), "5", "5");
    let m = assembly::uniform_mass(
        &a,
        &q("2"),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(m.grams.lo, q("12"));
    assert!(m.grams.is_exact());
    let c = m.centre_of_gravity.unwrap();
    for (c, v) in c.iter().zip(["9", "20.5", "31.5"]) {
        assert!(c.is_exact());
        assert_eq!(c.lo, q(v));
    }
}
#[test]
fn curved_control_and_irrational_gap_are_enclosed() {
    let sphere = Plan::compile(&Recipe::parse(r#"{"format":"catallax-recipe/3","name":"sphere","steps":[{"id":"s","op":"sphere","radius":"1","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"s"}"#).unwrap()).unwrap();
    let a = placed(&sphere, ["0", "0", "0"], ["1", "0", "0", "0"]);
    let b = placed(&sphere, ["3", "0", "0"], ["1", "0", "0", "0"]);
    gap(compare(&a, &b), "1", "1");
    let c = placed(&sphere, ["1", "0", "0"], ["1", "0", "0", "0"]);
    assert!(matches!(compare(&a, &c), Relation::Clash { .. }));
    let a = placed(&block(), ["0", "0", "0"], ["1", "0", "0", "0"]);
    let b = placed(&block(), ["2", "3", "0"], ["1", "0", "0", "0"]);
    let Relation::Clear { gap } = compare(&a, &b) else {
        panic!("expected irrational clear gap")
    };
    assert!(&gap.lo * &gap.lo <= q("2"));
    assert!(&gap.hi * &gap.hi >= q("2"));
    assert!(gap.lo > q("0") && gap.width() <= options().gap_width);
}
#[test]
fn invalid_options_and_material_refuse() {
    let a = placed(&block(), ["0", "0", "0"], ["1", "0", "0", "0"]);
    for density in ["0", "-1"] {
        assert_eq!(
            assembly::uniform_mass(
                &a,
                &q(density),
                &mass::Accuracy::default(),
                &mass::Options::default(),
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
            "assembly.density"
        );
    }
    let bad = Options {
        gap_width: q("0"),
        max_pairs: 1,
    };
    assert!(
        matches!(assembly::relation(&a,&a,&bad,&mut Decider::default()),Relation::Unresolved(r) if r.code == "assembly.options")
    );
}
#[test]
fn overlapping_envelopes_do_not_prove_clash_or_clearance() {
    let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"hollow","steps":[{"id":"outer","op":"block","size":["10","10","10"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"inner","op":"block","size":["8","8","8"],"pose":{"rotation":["1","0","0","0"],"translation":["1","1","1"]}},{"id":"shell","op":"boolean","kind":"subtract","target":"outer","tool":"inner"}],"result":"shell"}"#).unwrap();
    let shell = placed(
        &Plan::compile(&recipe).unwrap(),
        ["0", "0", "0"],
        ["1", "0", "0", "0"],
    );
    let inside = placed(&block(), ["3", "3", "3"], ["1", "0", "0", "0"]);
    let bounded = Options {
        gap_width: q("1/1000000"),
        max_pairs: 1,
    };
    assert!(
        matches!(assembly::relation(&shell,&inside,&bounded,&mut Decider::default()),Relation::Unresolved(r) if r.code == "assembly.pair-budget")
    );
}

const PARTS: [(&str, &str); 6] = [
    ("P1", include_str!("../../../parts/P1-bracket.json")),
    ("P2", include_str!("../../../parts/P2-collar.json")),
    ("P3", include_str!("../../../parts/P3-filleted-boss.json")),
    ("P4", include_str!("../../../parts/P4-shelled-tray.json")),
    ("P5", include_str!("../../../parts/P5-flanged-housing.json")),
    ("P6", include_str!("../../../parts/P6-lattice-bracket.json")),
];
fn target_record(text: &str) -> String {
    let plan = match Plan::compile(&Recipe::parse(text).unwrap()) {
        Ok(p) => p,
        Err(r) => return format!("compile: {}", r.code),
    };
    let part = match PlacedPart::new(&plan, &Pose::identity(), &mut Decider::default()) {
        Ok(p) => p,
        Err(r) => return format!("placement: {}", r.code),
    };
    let mass = assembly::uniform_mass(
        &part,
        &q("0.00785"),
        &mass::Accuracy::default(),
        &mass::Options { max_cells: 128 },
        &mut Decider::default(),
    );
    let mass = match mass {
        Ok(Mass {
            grams,
            centre_of_gravity,
        }) => format!("mass={grams:?}; cog={centre_of_gravity:?}"),
        Err(r) => format!("mass unresolved: {}", r.code),
    };
    let other = placed(&plan, ["100", "0", "0"], ["1", "0", "0", "0"]);
    let relation = assembly::relation(
        &part,
        &other,
        &Options {
            gap_width: q("1/1000000"),
            max_pairs: 16,
        },
        &mut Decider::default(),
    );
    format!("{mass}; relation={relation:?}")
}
#[test]
fn complete_public_targets_are_identical_at_one_two_four_workers() {
    let baseline: Vec<_> = PARTS
        .iter()
        .map(|(name, text)| {
            let r = target_record(text);
            eprintln!("{name}: {r}");
            r
        })
        .collect();
    for (record, expected) in baseline.iter().zip([
        "mass unresolved: field.print-coverage; relation=Unresolved",
        "mass unresolved: field.print-coverage; relation=Unresolved",
        "mass unresolved: field.print-coverage; relation=Unresolved",
        "mass unresolved: field.mass-budget; relation=Clear",
        "compile: fillet.unsupported-overflow",
        "mass unresolved: field.mass-budget; relation=Clear",
    ]) {
        assert!(record.starts_with(expected), "{record}");
    }
    for workers in [1usize, 2, 4] {
        let result = std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers)
                .map(|worker| {
                    let start = worker * PARTS.len() / workers;
                    let end = (worker + 1) * PARTS.len() / workers;
                    let chunk = &PARTS[start..end];
                    assert!(!chunk.is_empty());
                    scope.spawn(move || {
                        chunk
                            .iter()
                            .map(|(_, text)| target_record(text))
                            .collect::<Vec<_>>()
                    })
                })
                .collect();
            assert_eq!(handles.len(), workers, "spawned worker count");
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        assert_eq!(result, baseline, "{workers} workers");
    }
}

#[test]
fn analytical_answers_and_original_proven_query_survive_parallel_calls() {
    use catallax_field::Workspace;
    use catallax_shapes::RPoint3;
    let plan = block();
    let query = || {
        let mut dec = Decider::default();
        let result = plan
            .proven_checked(
                &RPoint3::new(q("0.5"), q("1"), q("1.5")),
                &mut Workspace::default(),
                &mut dec,
            )
            .unwrap();
        (result, dec.stats().work)
    };
    let before = query();
    let record = |x: &str| {
        let a = placed(&plan, ["0", "0", "0"], ["1", "0", "0", "0"]);
        let b = placed(&plan, [x, "0", "0"], ["1", "0", "0", "0"]);
        let mut dec = Decider::default();
        let r = assembly::relation(&a, &b, &options(), &mut dec);
        (format!("{r:?}"), dec.stats().work)
    };
    let inputs = ["0", "0.5", "1", "1.000000000001", "4", "-4", "10", "100"];
    let expected: Vec<_> = inputs.iter().map(|x| record(x)).collect();
    for workers in [1usize, 2, 4] {
        let got = std::thread::scope(|s| {
            let jobs: Vec<_> = inputs
                .chunks(inputs.len().div_ceil(workers))
                .map(|chunk| {
                    let record = &record;
                    s.spawn(move || chunk.iter().map(|x| record(x)).collect::<Vec<_>>())
                })
                .collect();
            jobs.into_iter()
                .flat_map(|h| h.join().unwrap())
                .collect::<Vec<_>>()
        });
        assert_eq!(got, expected);
    }
    assert_eq!(before, query());
}

#[test]
fn union_cover_proves_global_gap_inside_overlapping_outer_boxes() {
    let recipe=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"L","steps":[{"id":"x","op":"block","size":["8","2","2"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"y","op":"block","size":["2","8","2"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"l","op":"boolean","kind":"union","target":"x","tool":"y"}],"result":"l"}"#).unwrap();
    let a = placed(
        &Plan::compile(&recipe).unwrap(),
        ["0", "0", "0"],
        ["1", "0", "0", "0"],
    );
    let b = placed(&block(), ["4", "4", "0"], ["1", "0", "0", "0"]);
    let options = Options {
        gap_width: q("1/1000000"),
        max_pairs: 512,
    };
    gap(
        assembly::relation(&a, &b, &options, &mut Decider::default()),
        "2",
        "2",
    );
    gap(
        assembly::relation(&b, &a, &options, &mut Decider::default()),
        "2",
        "2",
    );
}

#[test]
fn empty_operand_has_no_finite_gap_or_centre_of_gravity() {
    let recipe=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"empty","steps":[{"id":"box","op":"block","size":["1","2","3"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"empty","op":"boolean","kind":"subtract","target":"box","tool":"box"}],"result":"empty"}"#).unwrap();
    let empty = placed(
        &Plan::compile(&recipe).unwrap(),
        ["0", "0", "0"],
        ["1", "0", "0", "0"],
    );
    let a = placed(&block(), ["0", "0", "0"], ["1", "0", "0", "0"]);
    assert!(matches!(compare(&empty,&a),Relation::Unresolved(r) if r.code == "assembly.empty"));
    let m = assembly::uniform_mass(
        &empty,
        &q("2"),
        &mass::Accuracy::default(),
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(m.grams.is_exact() && m.grams.lo.is_zero());
    assert!(m.centre_of_gravity.is_none());
}
