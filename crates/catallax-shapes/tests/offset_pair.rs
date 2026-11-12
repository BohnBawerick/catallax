//! Exact offset inputs to M10's pair search: a plane control and the three
//! edges of the M11 offset-distance experiment, under the unchanged 20M
//! operation allowance.

use catallax_numbers::{Decider, MPoly, Rational, Real, Refusal, with_operation_budget};
use catallax_shapes::{BSplineSurface, OffsetPatch, PatchPairIntersection, Point3, RationalPatch};
use std::num::NonZeroUsize;

type Q = Rational;
const ALLOWANCE: u64 = 20_000_000;

fn q(n: i64) -> Q {
    Q::from_i64(n)
}
fn r(text: &str) -> Q {
    Q::parse(text).unwrap()
}

fn binomial(n: usize, k: usize) -> Q {
    (0..k).fold(q(1), |out, i| {
        &(&out * &q((n - i) as i64)) / &q((i + 1) as i64)
    })
}

/// The Bezier patch on `[0,1]^2` of the given degrees whose coordinates are
/// the polynomials `p` in `(u, v)`, with unit weights.
fn patch(p: &[MPoly; 3], degree: [usize; 2], dec: &mut Decider) -> RationalPatch {
    let mut controls = Vec::new();
    for a in 0..=degree[0] {
        for b in 0..=degree[1] {
            let mut point: [Q; 3] = std::array::from_fn(|_| q(0));
            for (value, poly) in point.iter_mut().zip(p) {
                for (m, c) in poly.terms() {
                    let [i, j] = [m.exponents()[0] as usize, m.exponents()[1] as usize];
                    assert!(i <= degree[0] && j <= degree[1]);
                    if i > a || j > b {
                        continue;
                    }
                    let f = &(&binomial(a, i) / &binomial(degree[0], i))
                        * &(&binomial(b, j) / &binomial(degree[1], j));
                    *value = &*value + &(&f * &c);
                }
            }
            let [x, y, z] = point;
            controls.push(Point3::new(x.into(), y.into(), z.into()));
        }
    }
    let knots = |d: usize| -> Vec<Q> { (0..2 * d + 2).map(|i| q(i64::from(i > d))).collect() };
    let count = controls.len();
    let carrier = BSplineSurface::new(
        degree,
        [knots(degree[0]), knots(degree[1])],
        [degree[0] + 1, degree[1] + 1],
        controls,
        vec![Real::one(); count],
        dec,
    )
    .unwrap();
    RationalPatch::new(carrier, dec).unwrap()
}

/// M11's three support pairs, with their degrees.
fn supports(index: usize) -> (&'static str, [MPoly; 3], [MPoly; 3], [usize; 2]) {
    let k = |n: i64| MPoly::constant(2, q(n));
    let (u, v) = (MPoly::var(2, 0), MPoly::var(2, 1));
    let one = k(1);
    let a = [
        &(&k(-1) + &u.scale(&q(2))) + &(&v * &(&one - &u)),
        &k(-1) + &(&u * &v),
        v.scale(&q(2)),
    ];
    let b = [
        &one - &(&u * &v),
        &(&k(-1) + &u.scale(&q(2))) + &v,
        v.scale(&q(2)),
    ];
    let c = [
        &(&one - &u.scale(&q(2))) - &v,
        &one + &(&v * &(&one - &u.scale(&q(2)))),
        v.scale(&q(2)),
    ];
    match index {
        0 => ("m27-loft-edge", a, b, [1, 1]),
        1 => ("loft-next-edge", b, c, [1, 1]),
        _ => {
            let bump = (&(&(&u * &(&one - &u)) * &v) * &(&one - &v)).scale(&r("1/5"));
            let mut first = a;
            let mut second = b;
            first[2] = &first[2] + &bump;
            second[2] = &second[2] - &bump;
            ("biquadratic-edge", first, second, [2, 2])
        }
    }
}

/// The inward offsets `P - N/(10 |N|)` of one M11 edge's supports.
fn offsets(index: usize, dec: &mut Decider) -> (&'static str, [OffsetPatch; 2]) {
    let (name, a, b, degree) = supports(index);
    let radius = -r("1/10");
    (
        name,
        [
            patch(&a, degree, dec).offset(radius.clone()),
            patch(&b, degree, dec).offset(radius),
        ],
    )
}

/// Each branch: whether it is closed, and its arcs written out.
type Branches = Vec<(bool, Vec<String>)>;

/// What a run proved, in a form that compares across thread counts.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Record {
    answer: Result<Branches, (String, String)>,
    used: u64,
    exhausted: bool,
}

fn record(result: &Result<PatchPairIntersection, Refusal>, used: u64, exhausted: bool) -> Record {
    let answer = match result {
        Ok(found) => Ok(found
            .branches()
            .iter()
            .map(|branch| {
                let arcs = branch
                    .arcs()
                    .iter()
                    .map(|arc| {
                        format!(
                            "{:?} {:?} {:?} {:?} {:?}",
                            arc.first_cell,
                            arc.second_cell,
                            arc.certificate,
                            arc.ends[0].enclosure(),
                            arc.ends[1].enclosure()
                        )
                    })
                    .collect();
                (branch.is_closed(), arcs)
            })
            .collect()),
        Err(refusal) => Err((refusal.code.to_string(), refusal.detail.clone())),
    };
    Record {
        answer,
        used,
        exhausted,
    }
}

fn shape(result: &Result<PatchPairIntersection, Refusal>) -> String {
    match result {
        Ok(found) => {
            let closed = found.branches().iter().filter(|b| b.is_closed()).count();
            let s = found.search();
            let mut text = format!(
                "{} open, {closed} closed; {} arcs, {} cells, {} excluded, field {}/{}, cone {}, {} face systems, depth {}",
                found.branches().len() - closed,
                s.arcs,
                s.cells,
                s.excluded,
                s.field_first,
                s.field_second,
                s.normal_cone,
                s.face_systems,
                s.depth
            );
            for branch in found.branches() {
                let [start, end] = branch.ends();
                let at = |e: &[[Q; 2]; 4]| e.each_ref().map(|i| i[0].to_f64_nearest());
                text.push_str(&format!(
                    "; branch {:?} -> {:?}",
                    at(start.enclosure()),
                    at(end.enclosure())
                ));
            }
            text
        }
        Err(refusal) => format!("{}: {}", refusal.code, refusal.detail),
    }
}

fn run(index: usize, allowance: u64) -> (Record, String) {
    let mut dec = Decider::default();
    let (name, [a, b]) = offsets(index, &mut dec);
    let (result, work) = with_operation_budget(allowance, || {
        let mut dec = Decider::default();
        a.offset_meeting(&b, &mut dec)
    });
    let text = format!("{name}: {} units; {}", work.used, shape(&result));
    (record(&result, work.used, work.exhausted), text)
}

/// Offsets of two perpendicular planes are planes: `z = -1/10` and
/// `y = 1/10` meet along `v = 1/10`, `t = 2/5`, `s = (2u + 1)/4`, from the
/// first patch's border `u = 0` to `u = 1`.
#[test]
fn plane_offsets_meet_on_their_exact_line() {
    let mut dec = Decider::default();
    let k = |n: &str| MPoly::constant(2, r(n));
    let (u, v) = (MPoly::var(2, 0), MPoly::var(2, 1));
    let a = patch(&[u.clone(), v.clone(), k("0")], [1, 1], &mut dec);
    let b = patch(
        &[&u.scale(&q(2)) - &k("1/2"), k("0"), &v - &k("1/2")],
        [1, 1],
        &mut dec,
    );
    let (result, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        a.offset(-r("1/10"))
            .offset_meeting(&b.offset(-r("1/10")), &mut dec)
    });
    println!("planes: {} units; {}", work.used, shape(&result));
    let found = result.unwrap();
    let [branch] = found.branches() else {
        panic!("one branch");
    };
    assert!(!branch.is_closed());
    for arc in branch.arcs() {
        for end in &arc.ends {
            let e = end.enclosure();
            assert!(e[1][0] <= r("1/10") && r("1/10") <= e[1][1]);
            assert!(e[3][0] <= r("2/5") && r("2/5") <= e[3][1]);
            // s = (2u + 1)/4 on the line.
            let s = |u: &Q| &(&(u * &q(2)) + &q(1)) / &q(4);
            assert!(s(&e[0][0]) <= e[2][1] && e[2][0] <= s(&e[0][1]));
        }
    }
    let ends = branch.ends().map(|end| end.enclosure()[0].clone());
    assert!(ends.iter().any(|u| u == &[q(0), q(0)]));
    assert!(ends.iter().any(|u| u == &[q(1), q(1)]));
    // The reported point is the offset point `(u, 1/10, -1/10)`.
    for end in branch.ends() {
        let point = found.point(end);
        let u = end.enclosure()[0][0].to_f64_nearest();
        for (value, expected) in [(&point.x, u), (&point.y, 0.1), (&point.z, -0.1)] {
            assert!((value.to_f64() - expected).abs() < 1e-9);
        }
    }
    assert!(work.used <= ALLOWANCE);
}

/// The three M11 edges, each one request under the unchanged allowance, with
/// the same records at 1, 2 and 4 workers and at the default.
#[test]
#[ignore = "M11 opt-in: offset spines through the exact offset input"]
fn m11_offset_spines() {
    catallax_numbers::set_run_ahead_wait(true);
    let mut seen: Option<Vec<Record>> = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let mut records = Vec::new();
        for index in 0..3 {
            let (record, text) = run(index, ALLOWANCE);
            println!("M11 workers={threads:?} {text}");
            assert!(record.used <= ALLOWANCE);
            records.push(record);
        }
        match &seen {
            Some(first) => assert_eq!(first, &records, "{threads:?} workers"),
            None => seen = Some(records),
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
    // Published: each spine is one open branch of three certified arcs.
    for (record, used) in seen.unwrap().iter().zip([2_332_939, 2_666_401, 6_793_048]) {
        let branches = record.answer.as_ref().expect("a proved spine");
        assert_eq!(branches.len(), 1);
        assert!(!branches[0].0);
        assert_eq!(branches[0].1.len(), 3);
        assert_eq!((record.used, record.exhausted), (used, false));
    }
}

#[test]
fn an_exhausted_allowance_refuses() {
    let (record, _) = run(0, 1);
    assert!(record.exhausted);
    let refusal = record.answer.expect_err("one unit cannot answer");
    assert_eq!(refusal.0, "numbers.operation-budget");
}
