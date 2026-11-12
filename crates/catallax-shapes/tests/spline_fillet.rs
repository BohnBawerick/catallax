//! M11: rolling-ball fillet faces between two rational spline faces, on the
//! proven spines of the M11 edges, under the unchanged 20M operation
//! allowance.

use catallax_numbers::{
    Ball, BudgetPolicy, Decider, MPoly, Rational, Real, Refusal, Sign, with_operation_budget,
};
use catallax_shapes::{
    BSplineSurface, FilletCut, Point3, Pose, RailBorder, RailSpan, RationalPatch, Rotation,
    SPLINE_FILLET_DEFINITION, SPLINE_FILLET_SUPPORT_SPANS, SPLINE_FILLET_TARGET_NM, SplineFillet,
    SplineFilletDefinition, SplineSupport, Vector3,
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
    fitted(name, definition(supports, range), allowance)
}

/// Fits one definition under an allowance and records the answer.
fn fitted(
    name: &str,
    definition: SplineFilletDefinition,
    allowance: u64,
) -> (Record, Option<SplineFillet>, String) {
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

/// Runs `edges` through `run` at 1, 2 and 4 workers and the default,
/// checks every record and carrier agrees, and returns the records with one
/// fit per edge.
fn matrix(
    edges: &[usize],
    run: impl Fn(usize) -> (Record, Option<SplineFillet>, String),
) -> Vec<(Record, Option<SplineFillet>)> {
    catallax_numbers::set_run_ahead_wait(true);
    let mut seen: Option<Vec<(Record, Option<SplineFillet>)>> = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        let mut records = Vec::new();
        for &index in edges {
            let (record, fillet, text) = run(index);
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
                        assert_eq!(fa.borders(), fb.borders(), "{threads:?} workers");
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
    let mut records = matrix(&[0], |index| run(index, ALLOWANCE));
    let (record, fillet) = records.remove(0);
    let fillet = fillet.expect("a proved face");
    let (spans, upper, _) = record.answer.clone().expect("a proved face");
    assert!(upper <= Q::from_ratio(SPLINE_FILLET_TARGET_NM, 1_000_000).unwrap());
    // Published: 13 spans, a bound below 78 nm, and the work.
    assert_eq!(
        (spans, nanometres(&upper), record.used, record.exhausted),
        (13, 78, 5_252_648, false)
    );
    let mut dec = Decider::default();
    assert_eq!(fillet.remeasure(&mut dec).unwrap(), upper);
    rails_on_their_faces(&fillet);
    on_the_ball(&fillet, 1.0 / 32.0, 1.0);
}

/// Checks carrier points over `tau` in `[from, to]` of an M27 face against
/// an independent binary64 centre: each on the ball, and each rail on its
/// contact, within the proved bound.
fn on_the_ball(fillet: &SplineFillet, from: f64, to: f64) {
    let dec = Decider::default();
    let bound = fillet.bound().upper.to_f64_nearest() + 1e-9;
    let at = |value: f64| Ball::from_rational(&Q::from_f64_exact(value).unwrap(), 128);
    let mut x = spine_point(from, [0.95, 0.02, 0.03]);
    for k in 0..=32 {
        let tau = from + (to - from) * f64::from(k) / 32.0;
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

/// Checks each returned rail cubic at exact points along every span: both
/// support parameters stay in the faces' domain `[0, 1]`.
fn rails_on_their_faces(fillet: &SplineFillet) {
    for side in 0..2 {
        for span in fillet.rail(side) {
            for parameter in &span.parameters {
                for k in 0..=64 {
                    let t = r(&format!("{k}/64"));
                    let s = &q(1) - &t;
                    let value = (0..4).fold(q(0), |sum, i| {
                        let basis = &(&binomial(3, i) * &t.pow(i as i32).unwrap())
                            * &s.pow(3 - i as i32).unwrap();
                        &sum + &(&basis * &parameter[i])
                    });
                    assert!(
                        value >= q(0) && value <= q(1),
                        "rail {side} leaves [0, 1] near tau {}",
                        span.tau[0]
                    );
                }
            }
        }
    }
}

/// M11's three edges at 1, 2 and 4 workers and the default.
#[test]
#[ignore = "M11 opt-in: fillet faces on the three M11 edges"]
fn m11_fillet_faces() {
    let records = matrix(&[0, 1, 2], |index| run(index, ALLOWANCE));
    // Published: spans, the bound in whole nanometres, and the work.
    for ((record, _), published) in records.iter().zip([
        (13, 78, 5_252_648),
        (16, 5, 6_970_934),
        (22, 80, 16_433_004),
    ]) {
        let (spans, upper, _) = record.answer.as_ref().expect("a proved face");
        assert_eq!((*spans, nanometres(upper), record.used), published);
        assert!(!record.exhausted);
    }
    for (_, fillet) in &records {
        rails_on_their_faces(fillet.as_ref().expect("a proved face"));
    }
}

/// A loft edge between the caps `z = 0` and `z = 2`, on supports extended
/// past their faces in `v` by a quarter each way: the range runs an eighth
/// past both caps, with a border where each rail meets `v = 0` and `v = 1`
/// of its support. The fillet on a loft body uses the same margins.
fn capped(index: usize) -> (&'static str, SplineFilletDefinition) {
    let mut dec = Decider::default();
    let (name, supports, _) = edge(index, &mut dec);
    let domain = [[q(0), q(1)], [r("-1/4"), r("5/4")]];
    let supports = supports.map(|patch| patch.extended(domain.clone(), &mut dec).unwrap());
    let borders = [(0, q(0)), (0, q(1)), (1, q(0)), (1, q(1))]
        .map(|(support, value)| RailBorder {
            support,
            parameter: 1,
            value,
        })
        .to_vec();
    let definition = definition(supports, [r("-1/8"), r("9/8")])
        .with_borders(borders)
        .unwrap();
    (name, definition)
}

fn run_capped(index: usize) -> (Record, Option<SplineFillet>, String) {
    let (name, definition) = capped(index);
    fitted(name, definition, ALLOWANCE)
}

/// Each cap's cut of a capped loft face, and the work of each.
fn cap_cuts(fillet: &SplineFillet) -> [(FilletCut, u64); 2] {
    [0, 2].map(|height| {
        let (cut, work) = with_operation_budget(ALLOWANCE, || {
            let mut dec = Decider::default();
            fillet.plane_cut([q(0), q(0), q(1)], q(height), &mut dec)
        });
        (cut.unwrap().expect("a cut"), work.used)
    })
}

/// The rail parameters at a node of `tau`.
fn rail_at(fillet: &SplineFillet, side: usize, tau: &Q) -> [Q; 2] {
    let span = fillet
        .rail(side)
        .iter()
        .find(|span| &span.tau[0] == tau)
        .expect("a span from the node");
    [span.parameters[0][0].clone(), span.parameters[1][0].clone()]
}

/// The M27 loft edge's face between the loft's caps, on extended supports,
/// at 1, 2 and 4 workers and the default. Each rail meets each cap at a
/// node, at a rational point, and each cap cuts the face once.
#[test]
fn m27_loft_edge_face_between_the_caps() {
    let mut records = matrix(&[0], run_capped);
    let (record, fillet) = records.remove(0);
    let fillet = fillet.expect("a proved face");
    let (spans, upper, _) = record.answer.clone().expect("a proved face");
    let [(bottom, bottom_work), (top, top_work)] = cap_cuts(&fillet);
    println!(
        "M27 between the caps: borders {:?}, cuts {:?} and {:?}, cut work {bottom_work} and {top_work}",
        fillet
            .borders()
            .iter()
            .map(Q::to_f64_nearest)
            .collect::<Vec<_>>(),
        bottom.spans,
        top.spans,
    );
    // Published: 18 spans, a bound of 27 nm rounded up, and the work.
    assert_eq!(
        (spans, nanometres(&upper), record.used, record.exhausted),
        (18, 27, 8_199_428, false)
    );
    // The first rail's borders are `tau = 0` and `1`; the second's lie at
    // binary64 nodes where its `v` is exactly 0 and 1.
    let borders = fillet.borders();
    assert_eq!(&borders[..2], &[q(0), q(1)]);
    assert!((borders[2].to_f64_nearest() - 0.020075084665161).abs() < 1e-12);
    assert!((borders[3].to_f64_nearest() - 1.024668970491692).abs() < 1e-12);
    assert_eq!(rail_at(&fillet, 1, &borders[2])[1], q(0));
    assert_eq!(rail_at(&fillet, 1, &borders[3])[1], q(1));
    assert_eq!(rail_at(&fillet, 0, &borders[0])[1], q(0));
    assert_eq!(bottom.ends, [q(0), borders[2].clone()]);
    assert_eq!(top.ends, [q(1), borders[3].clone()]);
    assert_eq!(
        (bottom.before, top.before),
        (Sign::Negative, Sign::Negative)
    );
    assert_eq!((bottom.spans, top.spans), ([0, 2], [14, 16]));
    assert_eq!((bottom_work, top_work), (4_808, 4_808));
    on_the_ball(&fillet, -1.0 / 8.0, 9.0 / 8.0);
}

/// Planes that miss the face, or cut it along `tau`, are told apart.
#[test]
fn a_plane_along_the_edge_is_refused() {
    let (record, fillet, _) = run(0, ALLOWANCE);
    assert!(record.answer.is_ok());
    let fillet = fillet.unwrap();
    let mut dec = Decider::default();
    let missed = fillet
        .plane_cut([q(0), q(0), q(1)], q(10), &mut dec)
        .unwrap();
    assert!(missed.is_none());
    let along = fillet
        .plane_cut([q(1), q(0), q(0)], r("19/20"), &mut dec)
        .expect_err("a plane along the edge");
    assert_eq!(along.code, catallax_shapes::SPLINE_FILLET_CUT);
}

#[test]
fn borders_and_extensions_are_checked() {
    let mut dec = Decider::default();
    let (_, supports, range) = edge(0, &mut dec);
    let border = |support: usize, parameter: usize, value: Q| RailBorder {
        support,
        parameter,
        value,
    };
    let with = |borders: Vec<RailBorder>| {
        definition(supports.clone(), range.clone())
            .with_borders(borders)
            .map(|_| ())
            .map_err(|e| e.code)
    };
    assert_eq!(
        with(vec![border(1, 1, q(2))]),
        Err(SPLINE_FILLET_DEFINITION)
    );
    assert_eq!(
        with(vec![border(0, 1, q(0))]),
        Err(SPLINE_FILLET_DEFINITION)
    );
    assert_eq!(
        with(vec![border(2, 0, q(0))]),
        Err(SPLINE_FILLET_DEFINITION)
    );
    assert_eq!(with(vec![border(1, 0, r("1/2"))]), Ok(()));
    // The second rail's `u` stays near 0.05, so it never meets `u = 1/2`.
    let definition = definition(supports.clone(), range.clone())
        .with_borders(vec![border(1, 0, r("1/2"))])
        .unwrap();
    let (result, _) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        SplineFillet::fit(definition, &mut dec)
    });
    let refusal = result.expect_err("no crossing");
    assert_eq!(refusal.code, SPLINE_FILLET_DEFINITION);
    assert!(refusal.detail.contains("at 0 places"), "{}", refusal.detail);
    // An extension keeps the patch's own domain inside the new one.
    let refusal = supports[0]
        .extended([[q(0), r("1/2")], [q(0), q(1)]], &mut dec)
        .expect_err("a smaller domain");
    assert_eq!(refusal.code, "shapes.rational-patch-domain");
}

/// Both loft edges between the caps at 1, 2 and 4 workers and the default.
#[test]
#[ignore = "M11 opt-in: fillet faces between the loft's caps"]
fn m11_faces_between_the_caps() {
    let records = matrix(&[0, 1], run_capped);
    // Published: spans, the bound in whole nanometres, and the work.
    for ((record, fillet), published) in records
        .iter()
        .zip([(18, 27, 8_199_428), (18, 28, 9_685_010)])
    {
        let (spans, upper, _) = record.answer.as_ref().expect("a proved face");
        assert_eq!((*spans, nanometres(upper), record.used), published);
        assert!(!record.exhausted);
        let fillet = fillet.as_ref().unwrap();
        let [(bottom, _), (top, _)] = cap_cuts(fillet);
        let borders = fillet.borders();
        println!(
            "M11 between the caps: borders {:?}",
            borders.iter().map(Q::to_f64_nearest).collect::<Vec<_>>()
        );
        assert_eq!(bottom.ends, [q(0), borders[2].clone()]);
        assert_eq!(top.ends, [q(1), borders[3].clone()]);
    }
}

/// The M27 face as a surface of a body: chart preimages, exact rail points,
/// normals and sides, line hits, and a rigid move that keeps them.
#[test]
fn m27_face_as_a_surface() {
    let (_, fillet, _) = run(0, ALLOWANCE);
    let fillet = fillet.expect("a proved face");
    let mut dec = Decider::default();
    let at = |value: &Q| Ball::from_rational(value, 128);
    let inside = |ball: &Ball, value: &Q| &ball.lo() <= value && value <= &ball.hi();
    let (tau, theta) = (r("1/2"), r("1/3"));
    let point = fillet
        .carrier()
        .evaluate(&at(&tau), &at(&theta), &dec)
        .unwrap();
    let surface = SplineSupport::SplineFillet(fillet.clone());
    let found = surface
        .parameters_of(&point, &dec)
        .unwrap()
        .expect("a preimage");
    assert!(
        inside(&found[0], &tau) && inside(&found[1], &theta),
        "{found:?}"
    );
    // A rail point at a node is exactly on the face, and its normal is the
    // support's, turned toward the ball inside the material.
    let span = &fillet.rail(0)[3];
    let [u, v] = [&span.parameters[0][0], &span.parameters[1][0]];
    let support = &fillet.supports()[0];
    let rail = support
        .point_at(&[u.clone().into(), v.clone().into()], &mut dec)
        .unwrap();
    assert_eq!(surface.side_at(&rail, &mut dec).unwrap(), Sign::Zero);
    let normal = surface.normal_at(&rail, &mut dec).unwrap();
    let own = support.normal_at(&rail, &mut dec).unwrap();
    assert!(normal.dot(&own).to_f64() < 0.0);
    let far = Point3::new(q(0).into(), q(0).into(), q(10).into());
    assert_eq!(surface.side_at(&far, &mut dec).unwrap(), Sign::Positive);
    // A line through the face along its binary64 normal there meets it once.
    let mid = point.each_ref().map(|b| b.mid());
    let step = r("1/1024");
    let along = |du: &Q, dv: &Q| {
        let near = fillet
            .carrier()
            .evaluate(&at(&(&tau + du)), &at(&(&theta + dv)), &dec)
            .unwrap();
        std::array::from_fn::<f64, 3, _>(|i| (near[i].mid() - mid[i].clone()).to_f64_nearest())
    };
    let [a, b] = [along(&step, &q(0)), along(&q(0), &step)];
    let n = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    let rational = |x: f64| Q::from_f64_exact(x).unwrap();
    let origin = Point3::new(
        mid[0].clone().into(),
        mid[1].clone().into(),
        mid[2].clone().into(),
    );
    let direction = Vector3::new(
        rational(n[0]).into(),
        rational(n[1]).into(),
        rational(n[2]).into(),
    );
    let hits = surface
        .line_hits(&origin, &direction, &mut dec)
        .unwrap()
        .expect("a line solver");
    assert_eq!(hits.len(), 1);
    let gap = [
        (&hits[0].x, &origin.x),
        (&hits[0].y, &origin.y),
        (&hits[0].z, &origin.z),
    ]
    .iter()
    .map(|(a, b)| {
        let d = (*a - *b).to_f64();
        d * d
    })
    .sum::<f64>();
    assert!(gap.sqrt() < 1e-12, "{}", gap.sqrt());
    let away = Vector3::new(q(1).into(), q(0).into(), q(0).into());
    assert!(
        surface
            .line_hits(&far, &away, &mut dec)
            .unwrap()
            .unwrap()
            .is_empty()
    );
    // Moved by a rational rotation and translation, the face keeps its proofs
    // and answers in the new frame.
    let rotation = Rotation::from_quaternion(q(1), q(0), q(0), q(1)).unwrap();
    let pose = Pose::from_rational(
        rotation,
        catallax_shapes::RVector3 {
            x: r("1/3"),
            y: q(2),
            z: r("-1/2"),
        },
    );
    let moved = fillet.transformed(&pose);
    assert!(moved.shares_definition_nodes(&fillet));
    let placed = moved
        .carrier()
        .evaluate(&at(&tau), &at(&theta), &dec)
        .unwrap();
    let found = SplineSupport::SplineFillet(moved.clone())
        .parameters_of(&placed, &dec)
        .unwrap()
        .expect("a preimage");
    assert!(
        inside(&found[0], &tau) && inside(&found[1], &theta),
        "{found:?}"
    );
    let moved_rail = pose.apply_point(&rail);
    assert!(moved.on_boundary(&moved_rail, &mut dec).unwrap());
    assert!(
        moved
            .plane_cut([q(0), q(0), q(1)], q(10), &mut dec)
            .unwrap()
            .is_none()
    );
}

/// A moved face's cut is reported in the caller's frame. The capped M27 face
/// is translated, and also turned a third of a turn about `(1, 1, 1)`, which
/// takes `z` to `x`. Its bottom cap's cut returns the caller's plane, and both
/// of the cut's ends, the placed rails' exact points, lie on it.
#[test]
fn a_moved_face_is_cut_in_the_callers_frame() {
    let (_, fillet, _) = run_capped(0);
    let fillet = fillet.expect("a proved face");
    let degree = fillet.carrier().degrees()[0];
    let shift = [q(3), q(-2), q(10)];
    for (rotation, normal) in [
        (Rotation::identity(), [q(0), q(0), q(1)]),
        (
            Rotation::from_quaternion(q(1), q(1), q(1), q(1)).unwrap(),
            [q(1), q(0), q(0)],
        ),
    ] {
        let pose = Pose::from_rational(
            rotation,
            catallax_shapes::RVector3 {
                x: shift[0].clone(),
                y: shift[1].clone(),
                z: shift[2].clone(),
            },
        );
        let moved = fillet.transformed(&pose);
        let offset = normal
            .iter()
            .zip(&shift)
            .fold(q(0), |sum, (n, t)| sum + n * t);
        let mut dec = Decider::default();
        let cut = moved
            .plane_cut(normal.clone(), offset.clone(), &mut dec)
            .unwrap()
            .expect("a cut");
        assert_eq!((&cut.normal, &cut.offset), (&normal, &offset));
        let nodes: Vec<Q> = moved.nodes().cloned().collect();
        for (side, tau) in cut.ends.iter().enumerate() {
            let node = nodes.iter().position(|n| n == tau).expect("a node");
            let point = &moved.carrier().controls()[node * degree * 3 + 2 * side];
            let value = [&point.x, &point.y, &point.z]
                .into_iter()
                .zip(&cut.normal)
                .fold(q(0), |sum, (x, n)| sum + &x.as_rational().unwrap() * n);
            assert_eq!(value, cut.offset, "rail {side}");
        }
    }
}

/// An extension reserves its exact arithmetic before it does it. On the
/// biquadratic patch, an endpoint of about 12,700 bits is squared, which a
/// per-decision allowance of 10,000 units does not cover, so the extension
/// stops at that reservation. With a passed deadline it stops at once.
#[test]
fn an_extension_reserves_its_arithmetic_first() {
    let mut dec = Decider::default();
    let (_, supports, _) = edge(2, &mut dec);
    let far = -q(3).pow(8000).unwrap();
    let mut small = Decider::new(BudgetPolicy {
        max_work: 10_000,
        ..Default::default()
    });
    let refusal = supports[0]
        .extended([[q(0), q(1)], [far, q(2)]], &mut small)
        .expect_err("a reservation refuses");
    assert_eq!(
        (refusal.code, refusal.detail.as_str()),
        (
            "numbers.budget-exhausted",
            "exact geometric proof arithmetic stopped"
        )
    );
    let mut late = Decider::default();
    late.set_deadline(Some(std::time::Instant::now()));
    let refusal = supports[0]
        .extended([[q(0), q(1)], [r("-1/4"), r("5/4")]], &mut late)
        .expect_err("a passed deadline");
    assert_eq!(refusal.code, "numbers.deadline");
}
