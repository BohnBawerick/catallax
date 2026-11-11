//! M10's spline-with-spline stage on R13's controls: the 2026 field method
//! with an exact proof layer, under the unchanged 20M operation allowance.

use catallax_numbers::{Decider, Rational, Real, Refusal, with_operation_budget};
use catallax_shapes::{
    BSplineSurface, PairCertificate, PatchPairIntersection, Point3, RationalPatch,
};
use std::num::NonZeroUsize;

type Q = Rational;
fn q(n: i64) -> Q {
    Q::from_i64(n)
}
fn r(text: &str) -> Q {
    Q::parse(text).unwrap()
}

#[path = "../../catallax-features/examples/r13/wavy_fixture.rs"]
mod wavy_fixture;

const ALLOWANCE: u64 = 20_000_000;

#[derive(Clone, Copy, Debug)]
enum Pose {
    Base,
    Moved,
    Turned,
    Scaled,
}
const POSES: [Pose; 4] = [Pose::Base, Pose::Moved, Pose::Turned, Pose::Scaled];

impl Pose {
    fn apply(self, [x, y, z]: [Q; 3]) -> [Q; 3] {
        match self {
            Pose::Base => [x, y, z],
            Pose::Moved => [&x + &r("17/3"), &y - &q(11), &z + &r("2/7")],
            Pose::Turned => [
                &(&r("3/5") * &x) + &(&r("4/5") * &z),
                y,
                &(&r("-4/5") * &x) + &(&r("3/5") * &z),
            ],
            Pose::Scaled => {
                let s = r("7/3");
                [&s * &x, &s * &y, &s * &z]
            }
        }
    }
    fn name(self) -> &'static str {
        match self {
            Pose::Base => "base",
            Pose::Moved => "moved",
            Pose::Turned => "turned",
            Pose::Scaled => "scaled",
        }
    }
}

/// A bicubic patch on knots 0,0,0,0,1,2,3,3,3,3 with controls `x = i`,
/// `y = j` and the given heights, in the given pose.
fn patch(height: impl Fn(usize, usize) -> Q, pose: Pose, dec: &mut Decider) -> RationalPatch {
    weighted(height, pose, Real::one(), dec)
}

/// The same patch with every weight `weight`: the same surface, carried by
/// homogeneous coordinates whose weight is not one.
fn weighted(
    height: impl Fn(usize, usize) -> Q,
    pose: Pose,
    weight: Real,
    dec: &mut Decider,
) -> RationalPatch {
    let knots: Vec<Q> = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3].map(q).to_vec();
    let mut controls = Vec::with_capacity(36);
    for i in 0..6 {
        for j in 0..6 {
            let [x, y, z] = pose.apply([q(i as i64), q(j as i64), height(i, j)]);
            controls.push(Point3::new(x.into(), y.into(), z.into()));
        }
    }
    let carrier = BSplineSurface::new(
        [3, 3],
        [knots.clone(), knots],
        [6, 6],
        controls,
        vec![weight; 36],
        dec,
    )
    .unwrap();
    RationalPatch::new(carrier, dec).unwrap()
}

/// `H = N2 + N3`: its controls are 0, 0, 1, 1, 0, 0.
fn hill(i: usize) -> Q {
    q(i64::from(i == 2 || i == 3))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Family {
    Wavy,
    Loop,
    Contact,
    SmallLoop,
    Disjoint,
    Open,
    TinyLoop,
}

impl Family {
    fn name(self) -> &'static str {
        match self {
            Family::Wavy => "m25-wavy",
            Family::Loop => "loop",
            Family::Contact => "contact",
            Family::SmallLoop => "small-loop",
            Family::Disjoint => "disjoint",
            Family::Open => "open-control",
            Family::TinyLoop => "tiny-loop",
        }
    }
    fn level(self) -> Q {
        let peak = r("225/256");
        match self {
            Family::Loop => r("3/10"),
            Family::Contact => peak,
            Family::SmallLoop => &peak - &r("1/1000000"),
            Family::Disjoint => &peak + &r("1/1000000"),
            Family::TinyLoop => &peak - &r("1/10000000000"),
            Family::Wavy | Family::Open => q(0),
        }
    }
    fn patches(self, pose: Pose, dec: &mut Decider) -> [RationalPatch; 2] {
        match self {
            Family::Wavy => {
                let nets = wavy_fixture::pair();
                [
                    patch(|i, j| nets[0][i][j][2].clone(), pose, dec),
                    patch(|i, j| nets[1][i][j][2].clone(), pose, dec),
                ]
            }
            Family::Open => [
                patch(|i, j| &hill(i) * &hill(j), pose, dec),
                patch(
                    |i, j| &(&hill(i) * &hill(j)) + &(&q(10) * &(&q(i as i64) - &r("5/2"))),
                    pose,
                    dec,
                ),
            ],
            _ => {
                let c = self.level();
                [
                    patch(|i, j| &hill(i) * &hill(j), pose, dec),
                    patch(|i, j| (&(&hill(i) * &hill(j)) + &c).mul_pow2(-1), pose, dec),
                ]
            }
        }
    }
}

/// R13's six families; with the four poses, its 24 pairs.
const R13: [Family; 6] = [
    Family::Wavy,
    Family::Loop,
    Family::Contact,
    Family::SmallLoop,
    Family::Disjoint,
    Family::Open,
];

/// Each branch: whether it is closed, and its arcs written out.
type Branches = Vec<(bool, Vec<String>)>;

/// What a run proved, in a form that compares across thread counts.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Summary {
    answer: Result<Branches, (String, String)>,
    work: u64,
}

fn summary(result: &Result<PatchPairIntersection, Refusal>, work: u64) -> Summary {
    let answer = match result {
        Ok(found) => Ok(found
            .branches()
            .iter()
            .map(|branch| {
                (
                    branch.is_closed(),
                    branch
                        .arcs()
                        .iter()
                        .map(|arc| {
                            format!(
                                "{} {:?} {} {:?} {:?} {:?} {:?}",
                                arc.first_span,
                                arc.first_cell,
                                arc.second_span,
                                arc.second_cell,
                                arc.certificate,
                                arc.ends[0].enclosure(),
                                arc.ends[1].enclosure()
                            )
                        })
                        .collect(),
                )
            })
            .collect()),
        Err(refusal) => Err((refusal.code.to_string(), refusal.detail.clone())),
    };
    Summary { answer, work }
}

fn run(family: Family, pose: Pose) -> (Result<PatchPairIntersection, Refusal>, u64) {
    let mut dec = Decider::default();
    let [a, b] = family.patches(pose, &mut dec);
    let (result, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        a.patch_meeting(&b, &mut dec)
    });
    assert!(
        !work.exhausted || result.is_err(),
        "an exhausted allowance cannot answer"
    );
    (result, work.used)
}

fn shape(result: &Result<PatchPairIntersection, Refusal>) -> String {
    match result {
        Ok(found) => {
            let closed = found.branches().iter().filter(|b| b.is_closed()).count();
            let open = found.branches().len() - closed;
            let s = found.search();
            format!(
                "{open} open, {closed} closed; {} arcs, {} cells, {} cover, {} excluded, field {}/{}, cone {}, {} face systems, depth {}",
                s.arcs,
                s.cells,
                s.cover,
                s.excluded,
                s.field_first,
                s.field_second,
                s.normal_cone,
                s.face_systems,
                s.depth
            )
        }
        Err(refusal) => format!("{}: {}", refusal.code, refusal.detail),
    }
}

fn check(family: Family, pose: Pose, result: &Result<PatchPairIntersection, Refusal>) {
    let branches = |expected_open: usize, expected_closed: usize| {
        let found = result
            .as_ref()
            .unwrap_or_else(|e| panic!("{}: {} {}", family.name(), e.code, e.detail));
        let closed = found.branches().iter().filter(|b| b.is_closed()).count();
        assert_eq!(
            (found.branches().len() - closed, closed),
            (expected_open, expected_closed),
            "{}",
            family.name()
        );
    };
    match family {
        Family::Wavy => branches(2, 0),
        Family::Loop | Family::SmallLoop | Family::TinyLoop => branches(0, 1),
        Family::Disjoint => branches(0, 0),
        Family::Open => branches(1, 0),
        Family::Contact => {
            let refusal = result
                .as_ref()
                .expect_err("an exact tangent contact refuses");
            assert_eq!(
                refusal.code, "intersect.spline-pair-tangent",
                "{}",
                refusal.detail
            );
            // The refusal names the hilltop in the caller's frame.
            let [x, y, z] = pose.apply([r("5/2"), r("5/2"), r("225/256")]);
            let point = format!("contact at ({x}, {y}, {z})");
            assert!(refusal.detail.contains(&point), "{}", refusal.detail);
        }
    }
}

/// The 24 pairs and the tiny loop in every pose, each inside the allowance.
/// `--nocapture` prints the work table.
#[test]
fn r13_pairs_and_the_tiny_loop_inside_the_allowance() {
    for family in R13.into_iter().chain([Family::TinyLoop]) {
        for pose in POSES {
            let (result, work) = run(family, pose);
            println!(
                "{}/{}: {work} units; {}",
                family.name(),
                pose.name(),
                shape(&result)
            );
            check(family, pose, &result);
            assert!(work <= ALLOWANCE);
        }
    }
}

#[test]
fn every_thread_count_gives_the_same_answer_and_work() {
    catallax_numbers::set_run_ahead_wait(true);
    let mut cases: Vec<(Family, Pose)> = R13.iter().map(|f| (*f, Pose::Base)).collect();
    cases.push((Family::Loop, Pose::Turned));
    cases.push((Family::TinyLoop, Pose::Base));
    for (family, pose) in cases {
        let mut seen: Option<Summary> = None;
        for threads in [Some(1), Some(2), Some(4), None] {
            catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
            let (result, work) = run(family, pose);
            let now = summary(&result, work);
            if let Some(first) = &seen {
                assert_eq!(
                    first,
                    &now,
                    "{}/{} at {threads:?} threads",
                    family.name(),
                    pose.name()
                );
            } else {
                seen = Some(now);
            }
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
}

#[test]
fn loop_arcs_carry_exact_certificates_and_join_exactly() {
    let (result, _) = run(Family::Loop, Pose::Base);
    let found = result.unwrap();
    let [branch] = found.branches() else {
        panic!("one loop");
    };
    assert!(branch.is_closed());
    let arcs = branch.arcs();
    for (k, arc) in arcs.iter().enumerate() {
        assert!(matches!(
            arc.certificate,
            PairCertificate::FieldFirst(_)
                | PairCertificate::FieldSecond(_)
                | PairCertificate::NormalCone(_)
        ));
        let next = &arcs[(k + 1) % arcs.len()];
        // Consecutive arcs share their join point's enclosure overlap.
        let a = arc.ends[1].enclosure();
        let b = next.ends[0].enclosure();
        for axis in 0..4 {
            assert!(a[axis][0] <= b[axis][1] && b[axis][0] <= a[axis][1]);
        }
    }
    // On this family both patches share the chart, so u = s and v = t.
    for arc in arcs {
        let e = arc.ends[0].enclosure();
        assert!(e[0][0] <= e[2][1] && e[2][0] <= e[0][1]);
        assert!(e[1][0] <= e[3][1] && e[3][0] <= e[1][1]);
    }
}

#[test]
fn a_constant_weight_keeps_the_loop() {
    let mut dec = Decider::default();
    let [plain, b] = Family::Loop.patches(Pose::Base, &mut dec);
    let heavy = weighted(
        |i, j| &hill(i) * &hill(j),
        Pose::Base,
        Real::from_i64(2),
        &mut dec,
    );
    let mut shapes = Vec::new();
    for a in [plain, heavy] {
        let (result, work) = with_operation_budget(ALLOWANCE, || {
            let mut dec = Decider::default();
            a.patch_meeting(&b, &mut dec)
        });
        assert!(!work.exhausted);
        check(Family::Loop, Pose::Base, &result);
        let found = result.unwrap();
        let arcs: Vec<_> = found.branches()[0]
            .arcs()
            .iter()
            .map(|arc| (arc.first_cell.clone(), arc.second_span))
            .collect();
        shapes.push(arcs);
    }
    assert_eq!(shapes[0], shapes[1]);
}

#[test]
fn an_exhausted_allowance_refuses() {
    let mut dec = Decider::default();
    let [a, b] = Family::Loop.patches(Pose::Base, &mut dec);
    let (result, work) = with_operation_budget(1, || {
        let mut dec = Decider::default();
        a.patch_meeting(&b, &mut dec)
    });
    assert!(work.exhausted);
    assert!(result.is_err());
}
