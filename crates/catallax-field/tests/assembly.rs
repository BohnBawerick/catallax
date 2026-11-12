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
        // Helper threads' adopted work reaches the request's operation count,
        // not the writer's decider statistics.
        let (r, work) = catallax_numbers::with_operation_budget(u64::MAX / 2, || {
            assembly::relation(&a, &b, &options(), &mut Decider::default())
        });
        (format!("{r:?}"), work.used)
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

fn target(text: &str) -> Plan {
    Plan::compile(&Recipe::parse(text).unwrap()).unwrap()
}
fn sized(size: [&str; 3]) -> Plan {
    let [x, y, z] = size;
    target(&format!(
        r#"{{"format":"catallax-recipe/3","name":"box","steps":[{{"id":"box","op":"block","size":["{x}","{y}","{z}"],"pose":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}}}],"result":"box"}}"#
    ))
}
const SAME: [&str; 4] = ["1", "0", "0", "0"];
/// Not axis permutations: cos = 3/5 about x and about z.
const ABOUT_X: [&str; 4] = ["2", "1", "0", "0"];
const ABOUT_Z: [&str; 4] = ["2", "0", "0", "1"];

#[derive(Debug)]
enum Expect {
    /// Clear, with the enclosure containing this exact gap.
    Clear(&'static str),
    Clash,
    Unresolved(&'static str),
}
struct Case {
    label: &'static str,
    a: PlacedPart,
    b: PlacedPart,
    pairs: usize,
    width: &'static str,
    expect: Expect,
}
impl Case {
    fn run(&self) -> Relation {
        let options = Options {
            gap_width: q(self.width),
            max_pairs: self.pairs,
        };
        assembly::relation(&self.a, &self.b, &options, &mut Decider::default())
    }
    fn check(&self, result: &Relation) {
        let label = self.label;
        match (&self.expect, result) {
            (Expect::Clear(value), Relation::Clear { gap }) => {
                assert!(gap.lo > q("0"), "{label}: {gap:?}");
                assert!(gap.contains(&q(value)), "{label}: {gap:?}");
                assert!(gap.width() <= q(self.width), "{label}: {gap:?}");
            }
            (Expect::Clash, Relation::Clash { witness }) => {
                let point = catallax_shapes::RPoint3::new(
                    witness[0].clone(),
                    witness[1].clone(),
                    witness[2].clone(),
                );
                for part in [&self.a, &self.b] {
                    let at = part
                        .plan()
                        .proven_checked(
                            &point,
                            &mut catallax_field::Workspace::default(),
                            &mut Decider::default(),
                        )
                        .unwrap();
                    assert_eq!(at.location, catallax_field::Location::Inside, "{label}");
                }
            }
            (Expect::Unresolved(code), Relation::Unresolved(r)) => {
                assert_eq!(&r.code, code, "{label}: {r:?}")
            }
            (expect, result) => panic!("{label}: expected {expect:?}, got {result:?}"),
        }
    }
}
fn case(
    label: &'static str,
    (a, b): (PlacedPart, PlacedPart),
    pairs: usize,
    width: &'static str,
    expect: Expect,
) -> Case {
    Case {
        label,
        a,
        b,
        pairs,
        width,
        expect,
    }
}
fn copies(plan: &Plan, at: &str, rotation: [&str; 4]) -> (PlacedPart, PlacedPart) {
    (
        placed(plan, ["0", "0", "0"], SAME),
        placed(plan, [at, "0", "0"], rotation),
    )
}
// Second copies along x. The x extents are P1 60 mm, P4 40 mm and P6 12 mm,
// each bounded by planar faces normal to x.
fn cases() -> Vec<Case> {
    use Expect::*;
    let p1 = target(PARTS[0].1);
    let p2 = target(PARTS[1].1);
    let p3 = target(PARTS[2].1);
    let p4 = target(PARTS[3].1);
    let p6 = target(PARTS[5].1);
    let small = sized(["1", "2", "3"]);
    let fine = "1/1000000";
    vec![
        case(
            "P1 touching",
            copies(&p1, "60", SAME),
            512,
            fine,
            Unresolved("assembly.contact"),
        ),
        case(
            "P1 clear by 1 um",
            copies(&p1, "60.001", SAME),
            512,
            fine,
            Clear("1/1000"),
        ),
        case(
            "P1 interpenetrating",
            copies(&p1, "30", SAME),
            4096,
            fine,
            Clash,
        ),
        case("P1 coincident", copies(&p1, "0", SAME), 512, fine, Clash),
        case(
            "P1 rotated about x",
            copies(&p1, "60.001", ABOUT_X),
            512,
            fine,
            Clear("1/1000"),
        ),
        case(
            "P4 touching",
            copies(&p4, "40", SAME),
            512,
            fine,
            Unresolved("assembly.contact"),
        ),
        case(
            "P4 clear by 1 um",
            copies(&p4, "40.001", SAME),
            512,
            fine,
            Clear("1/1000"),
        ),
        case(
            "P4 interpenetrating",
            copies(&p4, "20", SAME),
            4096,
            fine,
            Clash,
        ),
        case(
            "P4 rotated about x",
            copies(&p4, "40.001", ABOUT_X),
            512,
            fine,
            Clear("1/1000"),
        ),
        case(
            "block buried in the P4 floor",
            (
                placed(&p4, ["0", "0", "0"], SAME),
                placed(&sized(["1", "2", "1"]), ["0", "0", "0.5"], SAME),
            ),
            64,
            fine,
            Clash,
        ),
        // Nested in the tray's cavity, 2 mm above its floor. The floor plane is
        // off the cover's dyadic grid, so the lower bound converges slowly; the
        // field reports the reached enclosure and refuses.
        case(
            "block nested in the P4 cavity",
            (
                placed(&p4, ["0", "0", "0"], SAME),
                placed(&small, ["0", "0", "4"], SAME),
            ),
            512,
            fine,
            Unresolved("assembly.pair-budget"),
        ),
        case(
            "P6 touching",
            copies(&p6, "12", SAME),
            512,
            fine,
            Unresolved("assembly.contact"),
        ),
        case(
            "P6 clear by 1 um",
            copies(&p6, "12.001", SAME),
            64,
            fine,
            Clear("1/1000"),
        ),
        case(
            "P6 interpenetrating",
            copies(&p6, "6", SAME),
            4096,
            fine,
            Clash,
        ),
        case(
            "P6 rotated about x",
            copies(&p6, "12.001", ABOUT_X),
            64,
            fine,
            Clear("1/1000"),
        ),
        case("P2 coincident", copies(&p2, "0", SAME), 512, fine, Clash),
        // Loose recipe bounds around the revolved collar: no certified upper
        // bound within the limit, so no clear answer.
        case(
            "P2 100 mm apart",
            copies(&p2, "100", SAME),
            512,
            fine,
            Unresolved("assembly.pair-budget"),
        ),
        case("P3 coincident", copies(&p3, "0", SAME), 512, fine, Clash),
        case(
            "P3 interpenetrating",
            copies(&p3, "10", SAME),
            4096,
            fine,
            Clash,
        ),
        case(
            "P3 100 mm apart",
            copies(&p3, "100", SAME),
            512,
            fine,
            Unresolved("assembly.pair-budget"),
        ),
        // A turned block's vertical edge at x = 3/5 faces an upright block.
        case(
            "turned block clear by 1 um",
            (
                placed(&small, ["0", "0", "0"], ABOUT_Z),
                placed(&small, ["0.601", "0", "0"], SAME),
            ),
            512,
            "1/20",
            Clear("1/1000"),
        ),
        case(
            "turned block interpenetrating",
            (
                placed(&small, ["0", "0", "0"], ABOUT_Z),
                placed(&small, ["0.5", "0", "0"], SAME),
            ),
            512,
            fine,
            Clash,
        ),
        // A touching turned edge has no rational contact proposal on the cover's grid.
        case(
            "turned block touching",
            (
                placed(&small, ["0", "0", "0"], ABOUT_Z),
                placed(&small, ["0.6", "0", "0"], SAME),
            ),
            512,
            fine,
            Unresolved("assembly.pair-budget"),
        ),
    ]
}

#[test]
fn placed_copies_answer_clash_clear_or_a_typed_reason() {
    for case in cases() {
        let result = case.run();
        eprintln!("{}: {result:?}", case.label);
        case.check(&result);
    }
}

#[test]
fn every_thread_count_gives_the_same_relation_and_work() {
    use std::num::NonZeroUsize;
    let picked = [
        "P1 clear by 1 um",
        "P4 interpenetrating",
        "P6 touching",
        "P2 100 mm apart",
        "turned block clear by 1 um",
    ];
    let cases: Vec<_> = cases()
        .into_iter()
        .filter(|c| picked.contains(&c.label))
        .collect();
    assert_eq!(cases.len(), picked.len());
    let record = |case: &Case| {
        let (result, work) = catallax_numbers::with_operation_budget(u64::MAX / 2, || case.run());
        assert!(!work.exhausted);
        (format!("{result:?}"), work.used)
    };
    for wait in [true, false] {
        catallax_numbers::set_run_ahead_wait(wait);
        for case in &cases {
            let mut seen = None;
            for threads in [Some(1), Some(2), Some(4), None] {
                catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
                let now = record(case);
                case.check(&case.run());
                match &seen {
                    Some(first) => assert_eq!(first, &now, "{} at {threads:?} threads", case.label),
                    None => seen = Some(now),
                }
            }
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
}

#[test]
fn every_thread_count_spends_a_short_allowance_the_same_way() {
    use std::num::NonZeroUsize;
    let picked = ["P1 clear by 1 um", "P4 interpenetrating"];
    let cases: Vec<_> = cases()
        .into_iter()
        .filter(|c| picked.contains(&c.label))
        .collect();
    assert_eq!(cases.len(), picked.len());
    let record = |case: &Case, limit: u64| {
        let (result, work) = catallax_numbers::with_operation_budget(limit, || case.run());
        (format!("{result:?}"), work.used, work.exhausted)
    };
    for case in &cases {
        catallax_numbers::set_threads(NonZeroUsize::new(1));
        let (_, full, exhausted) = record(case, u64::MAX / 2);
        assert!(!exhausted);
        // An allowance short of the whole search stops it at the same item,
        // with the same refusal and the same work, on every thread count.
        for limit in [full / 4, full / 2, full - 1] {
            for wait in [true, false] {
                catallax_numbers::set_run_ahead_wait(wait);
                let mut seen = None;
                for threads in [Some(1), Some(2), Some(4), None] {
                    catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
                    let now = record(case, limit);
                    assert!(now.2, "{} limit {limit}: {now:?}", case.label);
                    match &seen {
                        Some(first) => assert_eq!(
                            first, &now,
                            "{} limit {limit} wait {wait} at {threads:?} threads",
                            case.label
                        ),
                        None => seen = Some(now),
                    }
                }
            }
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
}
