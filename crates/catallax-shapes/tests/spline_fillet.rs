//! M11: rolling-ball fillet faces between two rational spline faces, on the
//! proven spines of the M11 edges, under the unchanged 20M operation
//! allowance.

use catallax_numbers::{
    Ball, Decider, MPoly, Rational, Real, Refusal, Sign, with_operation_budget,
};
use catallax_shapes::{
    BSplineSurface, Point3, RailSpan, RationalPatch, SPLINE_FILLET_DEFINITION,
    SPLINE_FILLET_SUPPORT_SPANS, SPLINE_FILLET_TARGET_NM, SplineFillet, SplineFilletDefinition,
};
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

/// M11's three edges: their names, support pairs, and the range of the first
/// support's `v` each face covers, inside the part of the spine where both
/// contacts lie on their faces.
fn edge(index: usize, dec: &mut Decider) -> (&'static str, [RationalPatch; 2], [Q; 2]) {
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
    let (name, first, second, degree, start) = match index {
        0 => ("m27-loft-edge", a, b, [1, 1], r("1/32")),
        1 => ("loft-next-edge", b, c, [1, 1], r("3/64")),
        _ => {
            let bump = (&(&(&u * &(&one - &u)) * &v) * &(&one - &v)).scale(&r("1/5"));
            let mut first = a;
            let mut second = b;
            first[2] = &first[2] + &bump;
            second[2] = &second[2] - &bump;
            ("biquadratic-edge", first, second, [2, 2], r("1/32"))
        }
    };
    (
        name,
        [patch(&first, degree, dec), patch(&second, degree, dec)],
        [start, q(1)],
    )
}

/// A ball of radius 1/10 inside the material of the edge, along the first
/// support's `v`.
fn definition(supports: [RationalPatch; 2], range: [Q; 2]) -> SplineFilletDefinition {
    SplineFilletDefinition::new(
        supports,
        r("1/10"),
        [Sign::Negative, Sign::Negative],
        1,
        range,
    )
    .unwrap()
}

/// A proved face: its spans, its bound and its two rails.
type Face = (usize, Q, [Vec<RailSpan>; 2]);

/// What a fit proved, in a form that compares across thread counts.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Record {
    answer: Result<Face, (String, String)>,
    used: u64,
    exhausted: bool,
}

fn run(index: usize, allowance: u64) -> (Record, Option<SplineFillet>, String) {
    let mut dec = Decider::default();
    let (name, supports, range) = edge(index, &mut dec);
    let definition = definition(supports, range);
    let (result, work) = with_operation_budget(allowance, || {
        let mut dec = Decider::default();
        SplineFillet::fit(definition, &mut dec)
    });
    let answer = match &result {
        Ok(fillet) => Ok((
            fillet.spans(),
            fillet.bound().upper.clone(),
            [fillet.rail(0).to_vec(), fillet.rail(1).to_vec()],
        )),
        Err(refusal) => Err((refusal.code.to_string(), refusal.detail.clone())),
    };
    let text = match &result {
        Ok(fillet) => format!(
            "{name}: {} units; {} spans, bound {} nm",
            work.used,
            fillet.spans(),
            nanometres(&fillet.bound().upper)
        ),
        Err(e) => format!("{name}: {} units; {}: {}", work.used, e.code, e.detail),
    };
    let record = Record {
        answer,
        used: work.used,
        exhausted: work.exhausted,
    };
    (record, result.ok(), text)
}

/// A bound in millimetres as whole nanometres, rounded up.
fn nanometres(upper: &Q) -> i64 {
    let scaled = upper * &q(1_000_000);
    scaled.ceil().to_string().parse().unwrap()
}

/// Runs `edges` at 1, 2 and 4 workers and the default, checks every record
/// and carrier agrees, and returns the records with one fit per edge.
fn matrix(edges: &[usize]) -> Vec<(Record, Option<SplineFillet>)> {
    catallax_numbers::set_run_ahead_wait(true);
    let mut seen: Option<Vec<(Record, Option<SplineFillet>)>> = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let mut records = Vec::new();
        for &index in edges {
            let (record, fillet, text) = run(index, ALLOWANCE);
            println!("M11 fillet workers={threads:?} {text}");
            assert!(record.used <= ALLOWANCE);
            records.push((record, fillet));
        }
        match &seen {
            Some(first) => {
                for ((a, fa), (b, fb)) in first.iter().zip(&records) {
                    assert_eq!(a, b, "{threads:?} workers");
                    if let (Some(fa), Some(fb)) = (fa, fb) {
                        assert!(
                            fa.carrier().same_definition(fb.carrier()),
                            "{threads:?} workers"
                        );
                    }
                }
            }
            None => seen = Some(records),
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
    seen.unwrap()
}

// The M27 loft edge's supports in binary64, for an independent check: the
// point and its two partials.
fn first(u: f64, v: f64) -> [[f64; 3]; 3] {
    [
        [-1.0 + 2.0 * u + v * (1.0 - u), -1.0 + u * v, 2.0 * v],
        [2.0 - v, v, 0.0],
        [1.0 - u, u, 2.0],
    ]
}
fn second(u: f64, v: f64) -> [[f64; 3]; 3] {
    [
        [1.0 - u * v, -1.0 + 2.0 * u + v, 2.0 * v],
        [-v, 2.0, 0.0],
        [-u, 1.0, 2.0],
    ]
}

/// The ball's centre `P - N/(10 |N|)` over a support point.
fn centre(jet: [[f64; 3]; 3]) -> [f64; 3] {
    let [p, a, b] = jet;
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let length = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    std::array::from_fn(|i| p[i] - 0.1 * n[i] / length)
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

fn det(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// The ball's centre on the M27 edge at `v = tau` of the first support, as
/// `(u_A, u_B, v_B)` by binary64 Newton from `guess`.
fn spine_point(tau: f64, guess: [f64; 3]) -> [f64; 3] {
    let residual = |x: [f64; 3]| -> [f64; 3] {
        let (a, b) = (centre(first(x[0], tau)), centre(second(x[1], x[2])));
        std::array::from_fn(|i| a[i] - b[i])
    };
    let mut x = guess;
    for _ in 0..50 {
        let f = residual(x);
        let mut jacobian = [[0.0; 3]; 3];
        for column in 0..3 {
            let step = 1e-7;
            let (mut up, mut down) = (x, x);
            up[column] += step;
            down[column] -= step;
            let (fu, fd) = (residual(up), residual(down));
            for (row, values) in jacobian.iter_mut().enumerate() {
                values[column] = (fu[row] - fd[row]) / (2.0 * step);
            }
        }
        let whole = det(jacobian);
        let delta: [f64; 3] = std::array::from_fn(|column| {
            let mut replaced = jacobian;
            for (row, values) in replaced.iter_mut().enumerate() {
                values[column] = f[row];
            }
            det(replaced) / whole
        });
        x = std::array::from_fn(|i| x[i] - delta[i]);
        if delta.iter().all(|d| d.abs() < 1e-15) {
            break;
        }
    }
    x
}

#[test]
fn definitions_are_checked() {
    let mut dec = Decider::default();
    let (_, supports, range) = edge(0, &mut dec);
    let refused = |radius: Q, sides: [Sign; 2], axis: usize, range: [Q; 2]| {
        SplineFilletDefinition::new(supports.clone(), radius, sides, axis, range)
            .expect_err("an unfittable definition")
            .code
    };
    let sides = [Sign::Negative, Sign::Negative];
    assert_eq!(
        refused(q(0), sides, 1, range.clone()),
        SPLINE_FILLET_DEFINITION
    );
    assert_eq!(
        refused(r("1/10"), [Sign::Zero, Sign::Negative], 1, range.clone()),
        SPLINE_FILLET_DEFINITION
    );
    assert_eq!(
        refused(r("1/10"), sides, 2, range.clone()),
        SPLINE_FILLET_DEFINITION
    );
    assert_eq!(
        refused(r("1/10"), sides, 1, [r("1/2"), r("1/4")]),
        SPLINE_FILLET_DEFINITION
    );
    assert_eq!(
        refused(r("1/10"), sides, 1, [q(0), q(2)]),
        SPLINE_FILLET_DEFINITION
    );
}

/// A support of two knot spans is refused before any work on the spine: the
/// M27 edge's first support with the knot `u = 1/2` inserted.
#[test]
fn a_support_of_several_spans_refuses() {
    let mut dec = Decider::default();
    let (_, [first, second], range) = edge(0, &mut dec);
    let corner = |x: &str, y: &str, z: &str| Point3::new(r(x).into(), r(y).into(), r(z).into());
    let split = BSplineSurface::new(
        [1, 1],
        [
            vec![q(0), q(0), r("1/2"), q(1), q(1)],
            vec![q(0), q(0), q(1), q(1)],
        ],
        [3, 2],
        vec![
            corner("-1", "-1", "0"),
            corner("0", "-1", "2"),
            corner("0", "-1", "0"),
            corner("1/2", "-1/2", "2"),
            corner("1", "-1", "0"),
            corner("1", "0", "2"),
        ],
        vec![Real::one(); 6],
        &mut dec,
    )
    .unwrap();
    let split = RationalPatch::new(split, &mut dec).unwrap();
    assert_eq!(first.domain(), split.domain());
    let (result, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        SplineFillet::fit(definition([split, second], range), &mut dec)
    });
    let refusal: Refusal = result.expect_err("a support of two spans");
    assert_eq!(refusal.code, SPLINE_FILLET_SUPPORT_SPANS);
    assert!(work.used < 100_000, "{} units", work.used);
}

#[test]
fn an_exhausted_allowance_refuses() {
    let (record, _, _) = run(0, 1_000_000);
    assert!(record.exhausted);
    let refusal = record
        .answer
        .expect_err("a million units cannot fit the face");
    assert_eq!(refusal.0, "numbers.operation-budget");
}

/// The M27 loft edge's fillet face, at 1, 2 and 4 workers and the default.
/// Each carrier point checked lies on the ball round an independently
/// computed centre, and each rail on its contact, within the proved bound.
#[test]
fn m27_loft_edge_fillet_face() {
    let mut records = matrix(&[0]);
    let (record, fillet) = records.remove(0);
    let fillet = fillet.expect("a proved face");
    let (spans, upper, _) = record.answer.clone().expect("a proved face");
    assert!(upper <= Q::from_ratio(SPLINE_FILLET_TARGET_NM, 1_000_000).unwrap());
    // Published: 13 spans, a bound below 78 nm, and the work.
    assert_eq!(
        (spans, nanometres(&upper), record.used, record.exhausted),
        (13, 78, 5_252_596, false)
    );
    let mut dec = Decider::default();
    assert_eq!(fillet.remeasure(&mut dec).unwrap(), upper);
    let bound = upper.to_f64_nearest() + 1e-9;
    let at = |value: f64| Ball::from_rational(&Q::from_f64_exact(value).unwrap(), 128);
    let mut x = [0.95, 0.02, 0.03];
    for k in 0..=32 {
        let tau = 1.0 / 32.0 + (1.0 - 1.0 / 32.0) * f64::from(k) / 32.0;
        x = spine_point(tau, x);
        let c = centre(first(x[0], tau));
        let contacts = [first(x[0], tau)[0], second(x[1], x[2])[0]];
        for j in 0..=8u32 {
            let theta = f64::from(j) / 8.0;
            let point = fillet
                .carrier()
                .evaluate(&at(tau), &at(theta), &dec)
                .unwrap()
                .map(|b| b.mid().to_f64_nearest());
            let off = (distance(point, c) - 0.1).abs();
            assert!(off <= bound, "tau {tau} theta {theta}: {off} off the ball");
            if j % 8 == 0 {
                let rail = distance(point, contacts[(j / 8) as usize]);
                assert!(
                    rail <= bound,
                    "tau {tau} rail {}: {rail} off the contact",
                    j / 8
                );
            }
        }
    }
}

/// M11's three edges at 1, 2 and 4 workers and the default.
#[test]
#[ignore = "M11 opt-in: fillet faces on the three M11 edges"]
fn m11_fillet_faces() {
    let records = matrix(&[0, 1, 2]);
    // Published: spans, the bound in whole nanometres, and the work.
    for ((record, _), published) in records.iter().zip([
        (13, 78, 5_252_596),
        (16, 5, 6_970_870),
        (22, 80, 16_432_916),
    ]) {
        let (spans, upper, _) = record.answer.as_ref().expect("a proved face");
        assert_eq!((*spans, nanometres(upper), record.used), published);
        assert!(!record.exhausted);
    }
}
