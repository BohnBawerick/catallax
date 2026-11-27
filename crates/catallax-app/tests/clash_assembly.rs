//! Roadmap 3.3, step 4: a test assembly of the public target parts P1 to P6,
//! with pins, a bolt, a block and a ball, placed distant, near, turned,
//! nested, touching and clashing in one model. Every pair gets clash, clear by
//! a stated gap, or unresolved with a typed reason, the same at 1, 2 and 4
//! threads. The field kernel's relation is the yardstick (Checkpoint 2
//! decision Q16): where both kernels answer, they agree.
use catallax_app::queries::{BodyRelation, body_relation, field_only_envelope};
use catallax_features::{Built, evaluate};
use catallax_field::{
    Plan,
    assembly::{self as field, PlacedPart, Relation},
};
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{Decider, Rational, Refusal, set_threads, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{Isometry, Pose, RVector3, Rotation};
use catallax_topology::{BodyId, transformed_body};
use std::num::NonZeroUsize;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

/// The parts, each a recipe of its own; the instances below place them.
fn parts() -> Vec<(&'static str, String)> {
    let primitive = |id: &str, op: &str, size: &str| {
        format!(
            r#"{{"format":"catallax-recipe/3","name":"{id}","steps":[{{"id":"{id}","op":"{op}",{size},"pose":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}}}],"result":"{id}"}}"#
        )
    };
    vec![
        ("p1", include_str!("../../../parts/P1-bracket.json").into()),
        ("p2", include_str!("../../../parts/P2-collar.json").into()),
        (
            "p3",
            include_str!("../../../parts/P3-filleted-boss.json").into(),
        ),
        (
            "p4",
            include_str!("../../../parts/P4-shelled-tray.json").into(),
        ),
        (
            "p5",
            include_str!("../../../parts/P5-flanged-housing.json").into(),
        ),
        (
            "p6",
            include_str!("../../../parts/P6-lattice-bracket.json").into(),
        ),
        (
            "pin7",
            primitive("pin", "cylinder", r#""radius":"7","height":"34""#),
        ),
        (
            "pin35",
            primitive("pin", "cylinder", r#""radius":"7/2","height":"20""#),
        ),
        (
            "bolt",
            primitive("bolt", "cylinder", r#""radius":"5/2","height":"9""#),
        ),
        (
            "block",
            primitive("block", "block", r#""size":["8","8","4"]"#),
        ),
        ("ball", primitive("ball", "sphere", r#""radius":"1/10""#)),
    ]
}

// Recipe fields that name an earlier step.
const READS: [&str; 7] = ["id", "body", "target", "tool", "profile", "region", "step"];

fn prefixed(value: &mut serde_json::Value, key: Option<&str>, ids: &[String], prefix: &str) {
    match value {
        serde_json::Value::Object(map) => {
            for (k, v) in map.iter_mut() {
                prefixed(v, Some(k), ids, prefix);
            }
        }
        serde_json::Value::Array(items) => {
            for v in items {
                prefixed(v, key, ids, prefix);
            }
        }
        serde_json::Value::String(s)
            if key.is_some_and(|k| READS.contains(&k)) && ids.contains(s) =>
        {
            *s = format!("{prefix}-{s}");
        }
        _ => (),
    }
}

/// One recipe holding every part, each step id prefixed with its part's name,
/// so all parts live in one model. Format 4, which P5 needs, lets fillets roll
/// over; P3, the only other part with fillets, builds the same faces and exact
/// volume under it as under format 3.
fn combined(parts: &[(&str, String)]) -> (Recipe, Vec<StepId>) {
    let mut steps = Vec::new();
    let mut results = Vec::new();
    for (prefix, text) in parts {
        let v: serde_json::Value = serde_json::from_str(text).unwrap();
        let ids: Vec<String> = v["steps"]
            .as_array()
            .unwrap()
            .iter()
            .map(|s| s["id"].as_str().unwrap().to_owned())
            .collect();
        for step in v["steps"].as_array().unwrap() {
            let mut step = step.clone();
            prefixed(&mut step, None, &ids, prefix);
            steps.push(step);
        }
        results.push(StepId::new(format!("{prefix}-{}", v["result"].as_str().unwrap())).unwrap());
    }
    let recipe = serde_json::json!({
        "format": "catallax-recipe/4",
        "name": "clash assembly",
        "steps": steps,
        "result": results[0].to_string(),
    });
    (Recipe::parse(&recipe.to_string()).unwrap(), results)
}

fn pose(quaternion: [&str; 4], to: [&str; 3]) -> Pose {
    let [w, x, y, z] = quaternion.map(q);
    Pose::from_rational(
        Rotation::from_quaternion(w, x, y, z).unwrap(),
        RVector3::new(q(to[0]), q(to[1]), q(to[2])),
    )
}
const IDENTITY: [&str; 4] = ["1", "0", "0", "0"];
/// A quarter turn about x, which takes y to z.
const UPRIGHT: [&str; 4] = ["1", "1", "0", "0"];
/// About z, with cosine 3/5 and sine 4/5.
const TURNED: [&str; 4] = ["2", "0", "0", "1"];
/// About the axis (1, 1, 1), by 2 acos(2 / sqrt 7).
const TUMBLED: [&str; 4] = ["2", "1", "1", "1"];

/// Where each part goes, in millimetres.
fn placements() -> Vec<(&'static str, &'static str, Pose)> {
    vec![
        ("tray P4", "p4", pose(IDENTITY, ["0", "0", "0"])),
        (
            "collar P2 standing in the tray",
            "p2",
            pose(UPRIGHT, ["0", "0", "3"]),
        ),
        (
            "lattice bracket P6 beside the tray",
            "p6",
            pose(IDENTITY, ["21", "-4", "0"]),
        ),
        (
            "housing P5 upright beyond the tray",
            "p5",
            pose(UPRIGHT, ["0", "46", "0"]),
        ),
        (
            "pin of radius 7 in the housing bore",
            "pin7",
            pose(IDENTITY, ["0", "46", "-2"]),
        ),
        (
            "boss P3 turned beside the tray",
            "p3",
            pose(TURNED, ["-30", "-60", "0"]),
        ),
        (
            "boss P3 turned into the tray corner",
            "p3",
            pose(TURNED, ["-10", "-50", "0"]),
        ),
        ("bracket P1", "p1", pose(IDENTITY, ["-30", "-130", "0"])),
        (
            "pin of radius 7/2 in the bracket bore",
            "pin35",
            pose(IDENTITY, ["0", "-105", "-2"]),
        ),
        (
            "bolt of radius 5/2 through a bracket hole",
            "bolt",
            pose(IDENTITY, ["-20", "-118", "-2"]),
        ),
        (
            "block resting on the bracket plate",
            "block",
            pose(IDENTITY, ["10", "-110", "5"]),
        ),
        (
            "lattice bracket P6 tumbled beside the housing",
            "p6",
            pose(TUMBLED, ["30", "46", "10"]),
        ),
        // At cavity height in the first P6's bore, off its axis toward +x, -y.
        (
            "ball in the lattice bracket's bore",
            "ball",
            pose(IDENTITY, ["109/4", "3/4", "3/2"]),
        ),
    ]
}

struct Instance {
    label: &'static str,
    body: BodyId,
    /// The B-rep carries this part as a field-only envelope.
    envelope: Option<Refusal>,
    field: Result<PlacedPart, Refusal>,
}

struct Assembly {
    built: Built,
    instances: Vec<Instance>,
}

fn assembly() -> Assembly {
    let parts = parts();
    let (recipe, results) = combined(&parts);
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap_or_else(|e| panic!("{e}"));
    let plans: Vec<Result<Plan, Refusal>> = parts
        .iter()
        .map(|(_, text)| Plan::compile(&Recipe::parse(text).unwrap()))
        .collect();
    let instances = placements()
        .into_iter()
        .map(|(label, name, pose)| {
            let part = parts.iter().position(|(n, _)| *n == name).unwrap();
            let source = built.bodies[&results[part]];
            let body = transformed_body(
                &mut built.model,
                source,
                &Isometry::from(&pose),
                Clone::clone,
            );
            Instance {
                label,
                body,
                envelope: field_only_envelope(&built, &results[part]),
                field: plans[part]
                    .clone()
                    .and_then(|plan| PlacedPart::new(&plan, &pose, &mut Decider::default())),
            }
        })
        .collect();
    Assembly { built, instances }
}

/// What a B-rep answer looks like, exactly, for comparison across runs.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Answered {
    Clash(Enclosure, Vec<&'static str>),
    Clear(Enclosure, Vec<&'static str>),
    Unresolved(&'static str),
}

/// The B-rep answer and work for one pair, under the application allowance.
/// A part the B-rep carries as a field-only envelope refuses first, as it
/// does through `ApplicationBody`.
fn brep(assembly: &Assembly, i: usize, j: usize) -> (Answered, u64) {
    let (x, y) = (&assembly.instances[i], &assembly.instances[j]);
    if let Some(refusal) = x.envelope.as_ref().or(y.envelope.as_ref()) {
        return (Answered::Unresolved(refusal.code), 0);
    }
    let (answer, work) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
        body_relation(
            &assembly.built.model,
            x.body,
            y.body,
            &Accuracy::default(),
            &mut Decider::default(),
        )
    });
    let answered = match answer {
        Ok(a) => {
            let marks = a.marks.iter().map(|m| m.code).collect();
            match a.value {
                BodyRelation::Clash { overlap } => Answered::Clash(overlap, marks),
                BodyRelation::Clear { gap } => Answered::Clear(gap, marks),
            }
        }
        Err(r) => Answered::Unresolved(r.refusal.code),
    };
    (answered, work.used)
}

/// What one pair answers in the B-rep.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Expect<S = &'static str> {
    /// Clear by this exact gap.
    Exact(S),
    /// Clear within the requested accuracy, unmarked: the gap's lower end to
    /// ten significant digits.
    Near(S),
    /// Clear by a proven gap wider than the requested accuracy, marked: its
    /// lower end to six significant digits.
    Wide(S),
    /// A proven positive overlap; marked when a ball and the boxes bound it.
    Clash { marked: bool },
    /// A typed refusal.
    Unresolved(&'static str),
}

use Expect::{Clash, Exact, Near, Unresolved, Wide};

fn expect(answer: &Answered) -> Expect<String> {
    let marked = |marks: &Vec<&str>| marks.contains(&"measure.accuracy-not-reached");
    match answer {
        Answered::Clash(overlap, marks) => {
            assert!(overlap.lo > q("0"), "{overlap:?}");
            Clash {
                marked: marked(marks),
            }
        }
        Answered::Clear(gap, marks) => {
            assert!(gap.lo > q("0"), "{gap:?}");
            if gap.lo == gap.hi {
                Exact(gap.lo.to_string())
            } else if marked(marks) {
                Wide(gap.lo.to_decimal_bounds(6).0)
            } else {
                assert!(
                    gap.width() <= &Accuracy::default().relative * &gap.lo,
                    "{gap:?}"
                );
                Near(gap.lo.to_decimal_bounds(10).0)
            }
        }
        Answered::Unresolved(code) => Unresolved(code),
    }
}

/// Every pair (i, j) with i < j, in order.
fn pairs(n: usize) -> Vec<(usize, usize)> {
    (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .collect()
}

/// A P6 instance: its lattice is field-only, so the B-rep holds an envelope.
const ENVELOPE: Expect = Unresolved("app.relation-envelope");
/// A P5 instance: its fillet is fitted, a tolerant participant until 3.0B.
const TOLERANT: Expect = Unresolved("app.relation-construction-bound");
const TOLERANT_ANSWER: Answered = Answered::Unresolved("app.relation-construction-bound");

/// What each pair answers, as the B-rep stated it at 1, 2 and 4 threads.
const EXPECTED: [(usize, usize, Expect); 78] = [
    (0, 1, Exact("1")),  // tray P4 / collar P2 standing in the tray
    (0, 2, ENVELOPE),    // tray P4 / lattice bracket P6 beside the tray
    (0, 3, TOLERANT),    // tray P4 / housing P5 upright beyond the tray
    (0, 4, Exact("24")), // tray P4 / pin of radius 7 in the housing bore
    (0, 5, Unresolved("numbers.operation-budget")), // tray P4 / boss P3 turned beside the tray
    (0, 6, Unresolved("numbers.operation-budget")), // tray P4 / boss P3 turned into the tray corner
    (0, 7, Exact("75")), // tray P4 / bracket P1
    (0, 8, Exact("173/2")), // tray P4 / pin of radius 7/2 in the bracket bore
    (0, 9, Near("100.5319182")), // tray P4 / bolt of radius 5/2 through a bracket hole
    (0, 10, Exact("87")), // tray P4 / block resting on the bracket plate
    (0, 11, ENVELOPE),   // tray P4 / lattice bracket P6 tumbled beside the housing
    (0, 12, Exact("143/20")), // tray P4 / ball in the lattice bracket's bore
    (1, 2, ENVELOPE),    // collar P2 standing in the tray / lattice bracket P6 beside the tray
    (1, 3, TOLERANT),    // collar P2 standing in the tray / housing P5 upright beyond the tray
    (1, 4, Exact("27")), // collar P2 standing in the tray / pin of radius 7 in the housing bore
    (1, 5, Unresolved("numbers.operation-budget")), // collar P2 standing in the tray / boss P3 turned beside the tray
    (1, 6, Unresolved("numbers.operation-budget")), // collar P2 standing in the tray / boss P3 turned into the tray corner
    (1, 7, Exact("78")),                            // collar P2 standing in the tray / bracket P1
    (1, 8, Exact("179/2")), // collar P2 standing in the tray / pin of radius 7/2 in the bracket bore
    (1, 9, Near("105.1829144")), // collar P2 standing in the tray / bolt of radius 5/2 through a bracket hole
    (1, 10, Near("90.4890238")), // collar P2 standing in the tray / block resting on the bracket plate
    (1, 11, ENVELOPE), // collar P2 standing in the tray / lattice bracket P6 tumbled beside the housing
    (1, 12, Near("15.36374277")), // collar P2 standing in the tray / ball in the lattice bracket's bore
    (2, 3, ENVELOPE), // lattice bracket P6 beside the tray / housing P5 upright beyond the tray
    (2, 4, ENVELOPE), // lattice bracket P6 beside the tray / pin of radius 7 in the housing bore
    (2, 5, ENVELOPE), // lattice bracket P6 beside the tray / boss P3 turned beside the tray
    (2, 6, ENVELOPE), // lattice bracket P6 beside the tray / boss P3 turned into the tray corner
    (2, 7, ENVELOPE), // lattice bracket P6 beside the tray / bracket P1
    (2, 8, ENVELOPE), // lattice bracket P6 beside the tray / pin of radius 7/2 in the bracket bore
    (2, 9, ENVELOPE), // lattice bracket P6 beside the tray / bolt of radius 5/2 through a bracket hole
    (2, 10, ENVELOPE), // lattice bracket P6 beside the tray / block resting on the bracket plate
    (2, 11, ENVELOPE), // lattice bracket P6 beside the tray / lattice bracket P6 tumbled beside the housing
    (2, 12, ENVELOPE), // lattice bracket P6 beside the tray / ball in the lattice bracket's bore
    (3, 4, TOLERANT),  // housing P5 upright beyond the tray / pin of radius 7 in the housing bore
    (3, 5, TOLERANT),  // housing P5 upright beyond the tray / boss P3 turned beside the tray
    (3, 6, TOLERANT),  // housing P5 upright beyond the tray / boss P3 turned into the tray corner
    (3, 7, TOLERANT),  // housing P5 upright beyond the tray / bracket P1
    (3, 8, TOLERANT),  // housing P5 upright beyond the tray / pin of radius 7/2 in the bracket bore
    (3, 9, TOLERANT), // housing P5 upright beyond the tray / bolt of radius 5/2 through a bracket hole
    (3, 10, TOLERANT), // housing P5 upright beyond the tray / block resting on the bracket plate
    (3, 11, ENVELOPE), // housing P5 upright beyond the tray / lattice bracket P6 tumbled beside the housing
    (3, 12, TOLERANT), // housing P5 upright beyond the tray / ball in the lattice bracket's bore
    (4, 5, Wide("48.9999")), // pin of radius 7 in the housing bore / boss P3 turned beside the tray
    (4, 6, Wide("38.9999")), // pin of radius 7 in the housing bore / boss P3 turned into the tray corner
    (4, 7, Exact("129")),    // pin of radius 7 in the housing bore / bracket P1
    (4, 8, Exact("281/2")), // pin of radius 7 in the housing bore / pin of radius 7/2 in the bracket bore
    (4, 9, Near("155.7150114")), // pin of radius 7 in the housing bore / bolt of radius 5/2 through a bracket hole
    (4, 10, Near("141.3374531")), // pin of radius 7 in the housing bore / block resting on the bracket plate
    (4, 11, ENVELOPE), // pin of radius 7 in the housing bore / lattice bracket P6 tumbled beside the housing
    (4, 12, Near("45.72163382")), // pin of radius 7 in the housing bore / ball in the lattice bracket's bore
    (5, 6, Clash { marked: true }), // boss P3 turned beside the tray / boss P3 turned into the tray corner
    (5, 7, Exact("30")),            // boss P3 turned beside the tray / bracket P1
    (5, 8, Wide("41.5752")), // boss P3 turned beside the tray / pin of radius 7/2 in the bracket bore
    (5, 9, Wide("55.536")), // boss P3 turned beside the tray / bolt of radius 5/2 through a bracket hole
    (5, 10, Wide("44.9444")), // boss P3 turned beside the tray / block resting on the bracket plate
    (5, 11, ENVELOPE), // boss P3 turned beside the tray / lattice bracket P6 tumbled beside the housing
    (5, 12, Wide("35.5964")), // boss P3 turned beside the tray / ball in the lattice bracket's bore
    (6, 7, Exact("40")), // boss P3 turned into the tray corner / bracket P1
    (6, 8, Wide("51.9085")), // boss P3 turned into the tray corner / pin of radius 7/2 in the bracket bore
    (6, 9, Wide("65.5305")), // boss P3 turned into the tray corner / bolt of radius 5/2 through a bracket hole
    (6, 10, Wide("51.9999")), // boss P3 turned into the tray corner / block resting on the bracket plate
    (6, 11, ENVELOPE), // boss P3 turned into the tray corner / lattice bracket P6 tumbled beside the housing
    (6, 12, Wide("15.1031")), // boss P3 turned into the tray corner / ball in the lattice bracket's bore
    (7, 8, Exact("1/2")),     // bracket P1 / pin of radius 7/2 in the bracket bore
    (7, 9, Clash { marked: false }), // bracket P1 / bolt of radius 5/2 through a bracket hole
    (7, 10, Unresolved("assembly.contact")), // bracket P1 / block resting on the bracket plate
    (7, 11, ENVELOPE),        // bracket P1 / lattice bracket P6 tumbled beside the housing
    (7, 12, Exact("1813/20")), // bracket P1 / ball in the lattice bracket's bore
    (8, 9, Near("17.85372088")), // pin of radius 7/2 in the bracket bore / bolt of radius 5/2 through a bracket hole
    (8, 10, Exact("13/2")), // pin of radius 7/2 in the bracket bore / block resting on the bracket plate
    (8, 11, ENVELOPE), // pin of radius 7/2 in the bracket bore / lattice bracket P6 tumbled beside the housing
    (8, 12, Near("105.6045099")), // pin of radius 7/2 in the bracket bore / ball in the lattice bracket's bore
    (9, 10, Near("28.54834939")), // bolt of radius 5/2 through a bracket hole / block resting on the bracket plate
    (9, 11, ENVELOPE), // bolt of radius 5/2 through a bracket hole / lattice bracket P6 tumbled beside the housing
    (9, 12, Near("125.2050272")), // bolt of radius 5/2 through a bracket hole / ball in the lattice bracket's bore
    (10, 11, ENVELOPE), // block resting on the bracket plate / lattice bracket P6 tumbled beside the housing
    (10, 12, Near("103.1248758")), // block resting on the bracket plate / ball in the lattice bracket's bore
    (11, 12, ENVELOPE), // lattice bracket P6 tumbled beside the housing / ball in the lattice bracket's bore
];

#[test]
fn every_pair_gets_one_of_three_answers_at_every_thread_count() {
    // Each count builds its own assembly. The geometry store builds a curve's
    // equations on first use and keeps them, also through the clones a query
    // works on, so a second query on one model can spend less work.
    let runs: Vec<(Assembly, Vec<(Answered, u64)>)> = [1, 2, 4]
        .map(|n| {
            set_threads(NonZeroUsize::new(n));
            let assembly = assembly();
            let records = pairs(assembly.instances.len())
                .into_iter()
                .map(|(i, j)| brep(&assembly, i, j))
                .collect();
            (assembly, records)
        })
        .into();
    set_threads(None);
    let (assembly, records) = &runs[0];
    let pairs = pairs(assembly.instances.len());
    let label = |i: usize, j: usize| {
        format!(
            "{} / {}",
            assembly.instances[i].label, assembly.instances[j].label
        )
    };
    let observed: Vec<String> = pairs
        .iter()
        .zip(records)
        .map(|(&(i, j), (answer, _))| format!("({i}, {j}, {:?})", expect(answer)))
        .collect();
    let table: String = observed
        .iter()
        .zip(pairs.iter().zip(records))
        .map(|(line, (&(i, j), (_, work)))| {
            format!("    {line}, // {}, {work} work\n", label(i, j))
        })
        .collect();
    println!("{table}");
    for (k, &(i, j)) in pairs.iter().enumerate() {
        let at: Vec<_> = runs.iter().map(|r| &r.1[k]).collect();
        assert!(at.iter().all(|r| *r == at[0]), "{}: {at:?}", label(i, j));
    }
    let expected: Vec<String> = EXPECTED
        .iter()
        .map(|(i, j, e)| format!("({i}, {j}, {e:?})"))
        .collect();
    assert_eq!(observed, expected, "observed:\n{table}");
}

/// Pairs where the field confirms the B-rep's clash or gap, of 78.
const FIELD_AGREES: usize = 30;
/// Pairs the B-rep leaves unresolved and the field answers, of 78.
const FIELD_ONLY: usize = 7;

#[test]
fn the_field_agrees_wherever_both_kernels_answer() {
    let assembly = assembly();
    let (mut agreed, mut field_only) = (Vec::new(), Vec::new());
    for (i, j) in pairs(assembly.instances.len()) {
        let (Ok(a), Ok(b)) = (&assembly.instances[i].field, &assembly.instances[j].field) else {
            continue;
        };
        let (answer, _) = brep(&assembly, i, j);
        // A tolerant participant (P5) is out of this milestone until 3.0B,
        // and its pairs have no B-rep answer to check; the field's P5 route
        // also spends about 16 s on each of them without an answer.
        if answer == TOLERANT_ANSWER {
            continue;
        }
        // A sixteenth of the B-rep's gap tells the faces at the minimum
        // from any other pair of faces, without a fine search.
        let gap_width = match &answer {
            Answered::Clear(gap, _) => gap.lo.mul_pow2(-4),
            _ => q("1/2"),
        };
        let (relation, _) = with_operation_budget(catallax_app::DEFAULT_BUILD_WORK, || {
            field::relation(
                a,
                b,
                &field::Options {
                    gap_width,
                    max_pairs: 1024,
                },
                &mut Decider::default(),
            )
        });
        let label = format!(
            "{} / {}",
            assembly.instances[i].label, assembly.instances[j].label
        );
        match (&answer, &relation) {
            (Answered::Clear(x, _), Relation::Clear { gap: y }) => {
                assert!(
                    x.lo <= y.hi && y.lo <= x.hi,
                    "{label}: B-rep {x:?}, field {y:?}"
                );
                agreed.push(label);
            }
            (Answered::Clash(..), Relation::Clash { .. }) => agreed.push(label),
            (Answered::Clear(..), Relation::Clash { .. })
            | (Answered::Clash(..), Relation::Clear { .. }) => {
                panic!("{label}: B-rep {answer:?}, field {relation:?}")
            }
            (Answered::Unresolved(_), Relation::Clear { .. } | Relation::Clash { .. }) => {
                field_only.push(label);
            }
            _ => (),
        }
    }
    println!(
        "agreed on {}: {agreed:#?}\nfield only on {}: {field_only:#?}",
        agreed.len(),
        field_only.len()
    );
    assert!(
        agreed.len() >= FIELD_AGREES && field_only.len() >= FIELD_ONLY,
        "agreed on {}: {agreed:#?}\nfield only on {}: {field_only:#?}",
        agreed.len(),
        field_only.len()
    );
}

#[test]
fn a_lattice_part_relates_through_the_field_not_its_envelope() {
    // P6's lattice is field-only. Its B-rep holds the lattice region as a
    // solid envelope, which at cavity height in the bore is solid where the
    // part has struts and open cells. A ball of radius 1/10 sits there, off
    // the bore's axis by sqrt(2)/4 toward +x, -y.
    let parts = [
        (
            "p6",
            include_str!("../../../parts/P6-lattice-bracket.json").to_owned(),
        ),
        ("ball", parts().pop().unwrap().1),
    ];
    let (recipe, results) = combined(&parts);
    let mut built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let at = pose(IDENTITY, ["25/4", "19/4", "3/2"]);
    let lattice = built.bodies[&results[0]];
    let ball = transformed_body(
        &mut built.model,
        built.bodies[&results[1]],
        &Isometry::from(&at),
        Clone::clone,
    );
    // The envelope's gap is the ball's distance to the bore wall,
    // 2/5 - sqrt(2)/4 = 0.0464466094...
    let envelope = match body_relation(
        &built.model,
        lattice,
        ball,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    {
        BodyRelation::Clear { gap } => gap,
        other => panic!("{other:?}"),
    };
    assert!(envelope.lo >= q("0.0464466094") && envelope.hi <= q("0.0464466095"));
    // The field holds the struts: no material lies that near the ball.
    let place = |text: &str, pose: &Pose| {
        let plan = Plan::compile(&Recipe::parse(text).unwrap()).unwrap();
        PlacedPart::new(&plan, pose, &mut Decider::default()).unwrap()
    };
    let relation = field::relation(
        &place(&parts[0].1, &Pose::identity()),
        &place(&parts[1].1, &at),
        &field::Options {
            gap_width: q("1/8"),
            max_pairs: 20_000,
        },
        &mut Decider::default(),
    );
    let Relation::Clear { gap } = relation else {
        panic!("{relation:?}");
    };
    assert!(gap.lo > envelope.hi, "field {gap:?}, envelope {envelope:?}");
    // So a relation never uses the envelope in the part's place.
    assert_eq!(
        field_only_envelope(&built, &results[0]).map(|r| r.code),
        Some("app.relation-envelope")
    );
    assert!(field_only_envelope(&built, &results[1]).is_none());
}

#[test]
fn application_queries_refuse_a_lattice_envelope() {
    let build = |text: &str| {
        catallax_app::build(
            &Recipe::parse(text).unwrap(),
            &catallax_mesh::MeshOptions::default(),
            &mut Decider::default(),
        )
        .unwrap()
    };
    let lattice = build(include_str!("../../../parts/P6-lattice-bracket.json"));
    let bore = lattice.built.bodies[&StepId::new("bore").unwrap()];
    let accuracy = Accuracy::default();
    let refusals = [
        lattice
            .body_relation(bore, &accuracy, &mut Decider::default())
            .map(|_| ())
            .unwrap_err(),
        lattice
            .overlap_volume(bore, &accuracy, &mut Decider::default())
            .map(|_| ())
            .unwrap_err(),
        lattice
            .contact_faces(bore, &mut Decider::default())
            .map(|_| ())
            .unwrap_err(),
    ];
    for refusal in refusals {
        assert_eq!(refusal.refusal.code, "app.relation-envelope");
        assert!(refusal.refusal.detail.contains("infill"), "{refusal}");
        assert_eq!(refusal.operands.len(), 2);
    }
    // A part without field-only steps answers as before: the collar and its
    // ring before the chamfers overlap.
    let collar = build(include_str!("../../../parts/P2-collar.json"));
    let ring = collar.built.bodies[&StepId::new("ring").unwrap()];
    let answer = collar.body_relation(ring, &accuracy, &mut Decider::default());
    assert!(
        answer
            .as_ref()
            .err()
            .is_none_or(|r| r.refusal.code != "app.relation-envelope"),
        "{answer:?}"
    );
}
