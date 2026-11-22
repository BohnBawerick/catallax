//! 4.2a step 7.3: the field's own proof of the rolling ball's spine on the
//! three M11 edges.
//!
//! The supports are those of `catallax-shapes/tests/spline_fillet.rs`: the M27
//! loft edge, the next loft edge and the biquadratic edge, a ball of radius
//! 1/10 inside the material of each, `tau` the first support's `v`. The third
//! check solves the spine again in binary64 with other unknowns (the three
//! contact parameters, the centre written with a square root, a central
//! difference Jacobian), from code that shares nothing with the field's
//! proof, and finds every solution inside the proved tube.

use std::num::NonZeroUsize;
use std::sync::Mutex;
use std::time::Instant;

use catallax_field::spline_spine::sweep::{Sweep, cap_continuation};
use catallax_field::spline_spine::{Spine, SpineDefinition, TUBE_CAP};
use catallax_numbers::{Decider, MPoly, Rational, Real, Sign, with_operation_budget};
use catallax_shapes::{BSplineSurface, Point3, RPoint3, RVector3, RationalPatch};

type Q = Rational;

const ALLOWANCE: u64 = 20_000_000;
static THREAD_COUNT: Mutex<()> = Mutex::new(());

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

/// M11's three edges: name and support pair.
fn edge(index: usize, dec: &mut Decider) -> (&'static str, [RationalPatch; 2]) {
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
    let (name, first, second, degree) = match index {
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
    };
    (
        name,
        [patch(&first, degree, dec), patch(&second, degree, dec)],
    )
}

/// A ball of radius 1/10 inside the material, `tau` the first support's `v`,
/// seeded half way up, at a part tolerance of 1 um.
fn definition(index: usize, dec: &mut Decider) -> (&'static str, SpineDefinition) {
    let (name, supports) = edge(index, dec);
    (
        name,
        SpineDefinition {
            supports,
            radius: r("1/10"),
            sides: [Sign::Negative, Sign::Negative],
            axis: 1,
            seed: r("1/2"),
            tolerance: r("1/1000"),
        },
    )
}

/// What a proof found, in a form that compares across worker counts.
#[derive(Clone, Debug, PartialEq)]
struct Record {
    answer: Result<String, (String, String)>,
    used: u64,
}

fn run(definition: &SpineDefinition) -> (Record, Option<Spine>) {
    let (result, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        Spine::prove(definition, &mut dec)
    });
    let answer = match &result {
        Ok(spine) => Ok(format!(
            "{:?} {:?} {}",
            spine.spans(),
            spine.ends(),
            spine.eps()
        )),
        Err(e) => Err((e.code.to_string(), e.detail.clone())),
    };
    (
        Record {
            answer,
            used: work.used,
        },
        result.ok(),
    )
}

fn proved(index: usize) -> (&'static str, Spine, u64) {
    let mut dec = Decider::default();
    let (name, definition) = definition(index, &mut dec);
    let start = Instant::now();
    let (record, spine) = run(&definition);
    let spine = spine.unwrap_or_else(|| panic!("{name}: {:?}", record.answer));
    println!(
        "{name}: {} units, {:.3} s; {} spans, tube {:.3e} mm; ends tau {:?} and {:?}",
        record.used,
        start.elapsed().as_secs_f64(),
        spine.spans().len(),
        spine.eps(),
        spine.ends()[0].point[4],
        spine.ends()[1].point[4],
    );
    (name, spine, record.used)
}

#[test]
fn each_m11_spine_is_proved_inside_its_tube() {
    // The B-rep's M10 spines end on the second support's `v = 0` and on the
    // first support's `v = 1` (docs/M11-OFFSET-PROOF.md, "The spine by the M10
    // engine"); the M27 edge's ends lie near
    // `(u, v, s, t) = (0.94997, 0.02008, 0.05470, 0)` and
    // `(0.94181, 1, 0.03853, 0.97520)`, and the other two start near
    // `v = 0.0407` and `0.0199`.
    let starts = [0.0201, 0.0407, 0.0199];
    for index in 0..3 {
        let (name, spine, used) = proved(index);
        assert!(used <= ALLOWANCE, "{name}: {used}");
        assert!(spine.eps() <= TUBE_CAP, "{name}: {}", spine.eps());
        assert!(spine.eps() <= spine.target(), "{name}");
        let [first, last] = spine.ends();
        assert_eq!((first.parameter, first.border), (3, 0.0), "{name}");
        assert_eq!((last.parameter, last.border), (1, 1.0), "{name}");
        assert!(
            (first.point[4].mid() - starts[index]).abs() < 1e-4,
            "{name}: {:?}",
            first.point[4]
        );
        if index == 0 {
            let near = |x: catallax_numbers::Interval, value: f64| (x.mid() - value).abs() < 1e-5;
            let [_, _, _, u, v, s, t, _, _] = first.point;
            assert!(near(u, 0.94997) && near(v, 0.02008) && near(s, 0.05470));
            assert_eq!(t, catallax_numbers::Interval::point(0.0).unwrap());
            let [_, _, _, u, v, s, t, _, _] = last.point;
            assert!(near(u, 0.94181) && near(s, 0.03853) && near(t, 0.97520));
            assert_eq!(v, catallax_numbers::Interval::point(1.0).unwrap());
        }
        // The spans cover the spine from end to end.
        let spans = spine.spans();
        assert!(spans[0].tau[0] <= first.point[4].lo());
        assert_eq!(spans.last().unwrap().tau[1], 1.0);
    }
}

#[test]
fn the_spine_is_the_same_at_one_two_and_four_workers() {
    let _guard = THREAD_COUNT.lock().unwrap();
    let mut dec = Decider::default();
    let mut seen: Option<Vec<Record>> = None;
    for threads in [Some(1), Some(2), Some(4), None] {
        catallax_numbers::set_threads(threads.and_then(NonZeroUsize::new));
        // Give every metered run a fresh support so it pays for the same
        // Bézier extraction instead of reading a cache warmed by the first.
        let records: Vec<Record> = (0..3).map(|i| run(&definition(i, &mut dec).1).0).collect();
        match &seen {
            Some(first) => {
                for (index, (expected, found)) in first.iter().zip(&records).enumerate() {
                    assert_eq!(
                        expected.used, found.used,
                        "edge {index}, {threads:?} workers"
                    );
                    assert!(
                        expected.answer == found.answer,
                        "edge {index}, {threads:?} workers"
                    );
                }
            }
            None => seen = Some(records),
        }
    }
    catallax_numbers::set_threads(None);
}

/// A support's point and partials in binary64, written out by hand from the
/// polynomials above.
type Surface = fn(f64, f64) -> [[f64; 3]; 3];

fn support_a(u: f64, v: f64) -> [[f64; 3]; 3] {
    [
        [-1.0 + 2.0 * u + v * (1.0 - u), -1.0 + u * v, 2.0 * v],
        [2.0 - v, v, 0.0],
        [1.0 - u, u, 2.0],
    ]
}

fn support_b(u: f64, v: f64) -> [[f64; 3]; 3] {
    [
        [1.0 - u * v, -1.0 + 2.0 * u + v, 2.0 * v],
        [-v, 2.0, 0.0],
        [-u, 1.0, 2.0],
    ]
}

fn support_c(u: f64, v: f64) -> [[f64; 3]; 3] {
    [
        [1.0 - 2.0 * u - v, 1.0 + v * (1.0 - 2.0 * u), 2.0 * v],
        [-2.0, -2.0 * v, 0.0],
        [-1.0, 1.0 - 2.0 * u, 2.0],
    ]
}

/// `(1/5) u (1 - u) v (1 - v)` and its partials.
fn bump(u: f64, v: f64) -> [f64; 3] {
    [
        0.2 * u * (1.0 - u) * v * (1.0 - v),
        0.2 * (1.0 - 2.0 * u) * v * (1.0 - v),
        0.2 * u * (1.0 - u) * (1.0 - 2.0 * v),
    ]
}

fn bumped_a(u: f64, v: f64) -> [[f64; 3]; 3] {
    let mut out = support_a(u, v);
    let [h, hu, hv] = bump(u, v);
    out[0][2] += h;
    out[1][2] += hu;
    out[2][2] += hv;
    out
}

fn bumped_b(u: f64, v: f64) -> [[f64; 3]; 3] {
    let mut out = support_b(u, v);
    let [h, hu, hv] = bump(u, v);
    out[0][2] -= h;
    out[1][2] -= hu;
    out[2][2] -= hv;
    out
}

/// The ball's centre `P - N / (10 |N|)` over a support point.
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

fn det(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// The ball's centre at `v = tau` of the first support, as
/// `(u_A, u_B, v_B)` by binary64 Newton from `guess`, and the centre.
fn spine_point(first: Surface, second: Surface, tau: f64, guess: [f64; 3]) -> ([f64; 3], [f64; 3]) {
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
    (x, centre(first(x[0], tau)))
}

#[test]
fn an_independent_binary64_solve_lies_in_every_tube() {
    let pairs: [(Surface, Surface); 3] = [
        (support_a, support_b),
        (support_b, support_c),
        (bumped_a, bumped_b),
    ];
    for (index, (first, second)) in pairs.into_iter().enumerate() {
        let (name, spine, _) = proved(index);
        let start = spine.ends()[0].point;
        let mut guess = [start[3].mid(), start[5].mid(), start[6].mid()];
        let from = start[4].hi();
        let samples = 1000;
        let mut worst: f64 = 0.0;
        for k in 0..=samples {
            let tau = from + (1.0 - from) * k as f64 / samples as f64;
            let (x, c) = spine_point(first, second, tau, guess);
            guess = x;
            let tube = spine
                .at(tau)
                .unwrap_or_else(|| panic!("{name}: no span at {tau}"));
            // The binary64 solve is good to a few units in the last place;
            // 1e-12 is far above that and far below the tube.
            let slack = 1e-12;
            for (value, bounds) in c.iter().zip(&tube[0..3]).chain([
                (&x[0], &tube[3]),
                (&x[1], &tube[5]),
                (&x[2], &tube[6]),
            ]) {
                assert!(
                    bounds.lo() - slack <= *value && *value <= bounds.hi() + slack,
                    "{name}: tau {tau}: {value} outside {bounds:?}"
                );
            }
            for (value, bounds) in c.iter().zip(&tube[0..3]) {
                worst = worst.max((value - bounds.mid()).abs());
            }
        }
        println!(
            "{name}: {} binary64 centres inside the tube, farthest {worst:.3e} mm from its middle",
            samples + 1
        );
    }
}

#[test]
fn an_unreachable_tolerance_refuses_with_the_width_reached() {
    let mut dec = Decider::default();
    let (_, mut definition) = definition(0, &mut dec);
    definition.tolerance = r("1/1000000000000000000");
    let (record, _) = run(&definition);
    let (code, detail) = record.answer.unwrap_err();
    assert_eq!(code, "field.spline-spine-eps", "{detail}");
    let (width, _) = detail
        .strip_prefix("the tube reached ")
        .and_then(|reached| reached.split_once(" nm "))
        .unwrap_or_else(|| panic!("{detail}"));
    let width = width.parse::<f64>().unwrap();
    let limit_nm = definition.tolerance.to_f64_nearest() * 0.5e6;
    assert!(width.is_finite() && width > limit_nm, "{detail}");
    println!("{detail}");
}

#[test]
fn support_conversion_reserves_basis_arithmetic_before_computing_it() {
    let mut dec = Decider::default();
    let (_, definition) = definition(0, &mut dec);
    for support in &definition.supports {
        assert_eq!(support.bezier_spans(&dec).unwrap().count(), 1);
    }
    let (result, work) = with_operation_budget(7, || Spine::prove(&definition, &mut dec));
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
    assert_eq!(work.uncharged, 0, "{work:?}");
}

#[test]
fn part_tolerance_can_require_tighter_spans() {
    let mut dec = Decider::default();
    let (_, mut definition) = definition(0, &mut dec);
    let default = run(&definition).1.unwrap();
    definition.tolerance = r("1/250000000");
    let (record, tighter) = run(&definition);
    let tighter = tighter.unwrap_or_else(|| panic!("{:?}", record.answer));
    assert!(tighter.eps() <= 2e-9 && tighter.eps() < default.eps());
    assert!(tighter.spans().len() > default.spans().len());
}

#[test]
fn definitions_are_checked() {
    let mut dec = Decider::default();
    let (_, base) = definition(0, &mut dec);
    let refused = |change: &dyn Fn(&mut SpineDefinition)| {
        let mut definition = base.clone();
        change(&mut definition);
        run(&definition).0.answer.unwrap_err()
    };
    for (code, detail) in [
        refused(&|d| d.radius = q(0)),
        refused(&|d| d.sides[1] = Sign::Zero),
        refused(&|d| d.axis = 2),
        refused(&|d| d.seed = r("3/2")),
        refused(&|d| d.tolerance = q(0)),
    ] {
        assert_eq!(code, "field.spline-spine-definition", "{detail}");
    }
    // A ball outside the material of both supports: their offsets do not meet
    // over the squares.
    let (code, detail) = refused(&|d| d.sides = [Sign::Positive, Sign::Positive]);
    assert_eq!(code, "field.spline-spine", "{detail}");
    println!("outside both: {detail}");
}

#[test]
fn a_sweep_distance_accounts_for_the_short_arc_and_its_ends() {
    let mut dec = Decider::default();
    let z = MPoly::constant(2, q(0));
    let u = MPoly::var(2, 0);
    let v = MPoly::var(2, 1);
    let a = patch(&[z.clone(), u.clone(), v.scale(&q(2))], [1, 1], &mut dec);
    let b = patch(
        &[u, z, &v + &MPoly::constant(2, r("1/2"))],
        [1, 1],
        &mut dec,
    );
    let definition = SpineDefinition {
        supports: [a, b],
        radius: r("1/10"),
        sides: [Sign::Positive, Sign::Negative],
        axis: 1,
        seed: r("1/2"),
        tolerance: r("1/1000"),
    };
    let spine = Spine::prove(&definition, &mut dec).unwrap();
    let sweep = Sweep::new(&spine, &mut dec).unwrap();
    let width = r("1/1000000");
    // A 3-4-5 direction in the short arc at z=1, the same direction below
    // the end, and a point opposite the arc whose nearest point is a rail.
    for ([x, y, z], expected) in [
        (["1/25", "1/50", "1"], q(0)),
        (["1/25", "1/50", "1/4"], r("1/4")),
        (["0", "1/5", "1"], r("1/10")),
    ] {
        let p = RPoint3::new(r(x), r(y), r(z));
        let result = sweep.distance(&p, &width, &mut dec).unwrap();
        assert!(result.converged, "{p:?}: {result:?}");
        assert!(
            result.lower <= expected && result.upper.as_ref().unwrap() >= &expected,
            "{p:?}: expected {expected}, got {result:?}"
        );
    }
    let past = RPoint3::new(r("1/25"), r("1/50"), r("301/200"));
    let exact = r("1/200");
    let bound = r("1229/250000");
    let coarse = sweep
        .distance(&past, &(&r("1/100") - &bound), &mut dec)
        .unwrap();
    assert!(coarse.lower <= bound && coarse.upper.unwrap() >= exact);
    let result = sweep.distance_at_most(&past, &bound, &mut dec).unwrap();
    assert!(bound < result.lower && result.lower <= exact, "{result:?}");
    assert!(result.upper.unwrap() >= exact);
}

#[test]
fn independently_sampled_arcs_on_the_three_m11_edges_lie_in_the_sweep_range() {
    let mut dec = Decider::default();
    let pairs: [(Surface, Surface); 3] = [
        (support_a, support_b),
        (support_b, support_c),
        (bumped_a, bumped_b),
    ];
    let width = r("1/1000000");
    for (index, (first, second)) in pairs.into_iter().enumerate() {
        let (name, spine, _) = proved(index);
        let sweep = Sweep::new(&spine, &mut dec).unwrap_or_else(|e| panic!("{name}: {e:?}"));
        let start = spine.ends()[0].point;
        let [from, to] = [start[4].hi(), spine.ends()[1].point[4].lo()];
        let mut guess = [start[3].mid(), start[5].mid(), start[6].mid()];
        let mut durations = Vec::new();
        for i in 0..=100 {
            let tau = from + (to - from) * f64::from(i) / 100.0;
            let (x, centre) = spine_point(first, second, tau, guess);
            guess = x;
            let a = first(x[0], tau)[0];
            let b = second(x[1], x[2])[0];
            for angle in [0.0, 0.25, 0.5, 0.75, 1.0] {
                let radial: [f64; 3] = std::array::from_fn(|j| {
                    (a[j] - centre[j]) * (1.0 - angle) + (b[j] - centre[j]) * angle
                });
                let l = radial.iter().map(|x| x * x).sum::<f64>().sqrt();
                let p: [Q; 3] = std::array::from_fn(|j| {
                    Q::from_f64_exact(centre[j] + radial[j] * 0.1 / l).unwrap()
                });
                let [x, y, z] = p;
                let p = RPoint3::new(x, y, z);
                let before = Instant::now();
                let range = sweep.distance(&p, &width, &mut dec).unwrap();
                durations.push(before.elapsed());
                assert!(
                    range.lower <= r("1/1000000000000"),
                    "{name} {tau}: {range:?}"
                );
                assert!(
                    range.upper.as_ref().is_some_and(|u| u <= &width),
                    "{name} {tau}: {range:?}"
                );
            }
        }
        durations.sort_unstable();
        println!(
            "{name}: {} distance queries, median {:?}",
            durations.len(),
            durations[durations.len() / 2]
        );
    }
}

#[test]
fn the_m27_cap_keeps_the_blend_past_the_first_rail_end() {
    let mut dec = Decider::default();
    let (_, definition) = definition(0, &mut dec);
    let [original, continued] = cap_continuation(&definition, &mut dec).unwrap();
    let caps = [
        (
            RPoint3::new(q(0), q(0), q(0)),
            RVector3::new(q(0), q(0), q(1)),
        ),
        (
            RPoint3::new(q(0), q(0), q(2)),
            RVector3::new(q(0), q(0), q(-1)),
        ),
    ];
    let sweep = Sweep::between_caps(&continued, &caps, &mut dec).unwrap();
    let finite = Sweep::new(&original, &mut dec).unwrap();
    // tau = 1.01 is beyond the first rail's top border. The second rail
    // is still below the cap here. Its point must survive the cap cut.
    let tau = 1.01;
    let (parameters, _) = spine_point(support_a, support_b, tau, [0.94, 0.04, 0.98]);
    let rail = support_b(parameters[1], parameters[2])[0];
    assert!(rail[2] < 2.0);
    let [x, y, z] = rail.map(|x| Q::from_f64_exact(x).unwrap());
    let p = RPoint3::new(x, y, z);
    let width = r("1/1000000");
    let trimmed = sweep.distance(&p, &width, &mut dec).unwrap();
    let stopped = finite.distance(&p, &width, &mut dec).unwrap();
    assert!(
        trimmed.upper.as_ref().is_some_and(|d| d <= &width),
        "{trimmed:?}"
    );
    assert!(stopped.lower > r("1/1000"), "{stopped:?}");
}

#[test]
fn capped_loft_edges_time_independent_points_near_both_cuts() {
    let mut dec = Decider::default();
    let width = r("1/1000000");
    let caps = [
        (RPoint3::origin(), RVector3::new(q(0), q(0), q(1))),
        (
            RPoint3::new(q(0), q(0), q(2)),
            RVector3::new(q(0), q(0), q(-1)),
        ),
    ];
    let pairs: [(Surface, Surface); 2] = [(support_a, support_b), (support_b, support_c)];
    for (index, (first, second)) in pairs.into_iter().enumerate() {
        let (name, definition) = definition(index, &mut dec);
        let [_, continued] = cap_continuation(&definition, &mut dec).unwrap();
        let sweep = Sweep::between_caps(&continued, &caps, &mut dec).unwrap();

        for (cap, direction) in [(0.0_f64, 1.0), (2.0, -1.0)] {
            for clipped in [false, true] {
                let mut times = Vec::new();
                for i in 1..10 {
                    let angle = f64::from(i) / 10.0;
                    let height = if clipped { cap } else { cap + direction * 0.2 };
                    let section = |tau| {
                        let ([u, s, t], centre) =
                            spine_point(first, second, tau, [0.94, 0.05, tau]);
                        let (a, b) = (first(u, tau)[0], second(s, t)[0]);
                        let radial: [f64; 3] = std::array::from_fn(|j| {
                            (a[j] - centre[j]) * (1.0 - angle) + (b[j] - centre[j]) * angle
                        });
                        let length = radial.iter().map(|x| x * x).sum::<f64>().sqrt();
                        let point: [f64; 3] =
                            std::array::from_fn(|j| centre[j] + radial[j] * 0.1 / length);
                        (point, [a[2], b[2]], centre[2])
                    };
                    let (mut lo, mut hi) = (height / 2.0 - 0.15, height / 2.0 + 0.15);
                    assert!(section(lo).0[2] < height && section(hi).0[2] > height);
                    for _ in 0..44 {
                        let tau = lo.midpoint(hi);
                        if section(tau).0[2] < height {
                            lo = tau;
                        } else {
                            hi = tau;
                        }
                    }
                    let (point, rails, centre_z) = section(lo.midpoint(hi));
                    assert!((point[2] - height).abs() < 1e-12);
                    if clipped {
                        if (rails[0] - cap) * (rails[1] - cap) >= 0.0 {
                            continue;
                        }
                    } else {
                        assert!(centre_z > 0.1 && centre_z < 1.9);
                    }
                    let [x, y, z] =
                        [point[0], point[1], height].map(|x| Q::from_f64_exact(x).unwrap());
                    let p = RPoint3::new(x, y, z);
                    let before = Instant::now();
                    let range = sweep.distance(&p, &width, &mut dec).unwrap();
                    times.push(before.elapsed());
                    assert!(
                        range.lower <= r("1/1000000000000"),
                        "{name} {p:?}: {range:?}"
                    );
                    assert!(
                        range.upper.as_ref().is_some_and(|u| u <= &width),
                        "{name} {p:?}: {range:?}"
                    );
                }
                assert!(times.len() >= 3, "{name}: cap {cap}, clipped {clipped}");
                times.sort_unstable();
                println!(
                    "{name}: cap {cap}, clipped {clipped}, {} points, median {:?}",
                    times.len(),
                    times[times.len() / 2]
                );
            }
        }
    }
}

fn near_tangent(slope: Q, dec: &mut Decider) -> SpineDefinition {
    let u = MPoly::var(2, 0);
    let v = MPoly::var(2, 1);
    let v2 = v.scale(&q(2));
    let a = patch(
        &[u.clone(), (&u * &u).scale(&r("1/16")), v2.clone()],
        [2, 1],
        dec,
    );
    // The skew along the edge gives the two contact rails different border
    // stations without changing the second face's plane or the shared edge.
    let b = patch(
        &[u.scale(&q(-1)), u.scale(&slope), &v2 + &u.scale(&r("1/8"))],
        [1, 1],
        dec,
    );
    SpineDefinition {
        supports: [a, b],
        radius: r("1/10"),
        sides: [Sign::Negative, Sign::Positive],
        axis: 1,
        seed: r("1/2"),
        tolerance: r("1/1000"),
    }
}

#[test]
fn a_plane_offset_can_enter_the_spline_spine_system_at_zero_offset() {
    let guard = THREAD_COUNT.lock().unwrap();
    let mut dec = Decider::default();
    let mut record = None;
    for workers in [1, 2, 4] {
        catallax_numbers::set_threads(NonZeroUsize::new(workers));
        let mixed = near_tangent(r("15/8"), &mut dec);
        let (spine, work) =
            with_operation_budget(ALLOWANCE, || Spine::prove_spline_plane(&mixed, 1, &mut dec));
        let spine = spine.unwrap();
        assert!(!work.exhausted);
        let signature = format!("{:?} {:?}", spine.spans(), spine.ends());
        if let Some((_, seen_work, seen)) = &record {
            assert_eq!((&work.used, &signature), (seen_work, seen));
        }
        record = Some((spine, work.used, signature));
    }
    catallax_numbers::set_threads(None);
    drop(guard);
    let (spine, work, _) = record.unwrap();
    let [x, y, z, ..] = spine.at(0.5).unwrap();
    // The plane is 15x + 8y = 0. Offset its exact unit normal by 1/10.
    let p = |x| catallax_numbers::Interval::point(x).unwrap();
    assert!((x * p(15.0) + y * p(8.0)).contains(1.7));
    assert!(z.contains(1.0));
    assert!(spine.eps() <= spine.target());
    let sweep = Sweep::new(&spine, &mut dec).unwrap();
    let ([cx, cy], [a, b]) = mixed_contact([15.0, 8.0]);
    let radial = [a[0] - cx + b[0] - cx, a[1] - cy + b[1] - cy];
    let scale = 0.1 / radial[0].hypot(radial[1]);
    let exact = |x| Q::from_f64_exact(x).unwrap();
    let width = r("1/1000000");
    let mut times = Vec::new();
    for i in 1..100 {
        let p = RPoint3::new(
            exact(cx + scale * radial[0]),
            exact(cy + scale * radial[1]),
            r("1/50") * q(i),
        );
        let before = Instant::now();
        let range = sweep.distance(&p, &width, &mut dec).unwrap();
        times.push(before.elapsed());
        assert!(range.lower <= r("1/1000000000000"));
        assert!(range.upper.is_some_and(|u| u <= width));
    }
    times.sort_unstable();
    println!(
        "spline-plane tube {:.6e} mm, {work} work, median {:?} at {} independent points",
        spine.eps(),
        times[times.len() / 2],
        times.len()
    );
}

fn mixed_contact(plane: [f64; 2]) -> ([f64; 2], [[f64; 2]; 2]) {
    let length = plane[0].hypot(plane[1]);
    let centre = |u: f64| {
        let slope = u / 8.0;
        let length = (1.0 + slope * slope).sqrt();
        [u - 0.1 * slope / length, u * u / 16.0 + 0.1 / length]
    };
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    while hi > lo.next_up() {
        let u = lo.midpoint(hi);
        let [x, y] = centre(u);
        if plane[0] * x + plane[1] * y < 0.1 * length {
            lo = u;
        } else {
            hi = u;
        }
    }
    let [cx, cy] = centre(lo);
    let feet = [
        [lo, lo * lo / 16.0],
        [cx - 0.1 * plane[0] / length, cy - 0.1 * plane[1] / length],
    ];
    ([cx, cy], feet)
}

#[test]
fn near_tangent_sweep_has_a_deterministic_certified_cutoff() {
    let mut dec = Decider::default();
    let mut last = None;
    for power in 0..=24 {
        let slope = Q::one().mul_pow2(-power);
        let definition = near_tangent(slope.clone(), &mut dec);
        let (result, work) = with_operation_budget(ALLOWANCE, || {
            let spine = Spine::prove(&definition, &mut dec)?;
            let sweep = Sweep::new(&spine, &mut dec)?;
            Ok::<_, catallax_numbers::Refusal>((spine.eps(), sweep.tail_bound()))
        });
        let degrees = slope.to_f64_nearest().atan().to_degrees();
        if power == 14 {
            let (_, tail) = result.as_ref().expect("the 2^-14 sweep must prove");
            let bound = tail.as_ref().expect("expected a certified rail bound");
            let (_, [a, b]) = mixed_contact([slope.to_f64_nearest(), 1.0]);
            let chord = (a[0] - b[0]).hypot(a[1] - b[1]);
            assert!(
                chord > 0.0 && &Q::from_f64_exact(chord).unwrap() <= bound,
                "rail chord {chord:e} mm exceeds the bound {bound} mm"
            );
        }
        match result {
            Ok((width, tail)) => {
                last = Some(power);
                println!("{degrees:e} degrees: {} work, {width:e} mm", work.used);
                if let Some(bound) = &tail {
                    assert!(bound <= &r("1/10000"));
                    println!("  labelled near-tangent rail bound {bound} mm");
                }
            }
            Err(refusal) => {
                println!(
                    "{degrees:e} degrees: {} work, {} {}",
                    work.used, refusal.code, refusal.detail
                );
            }
        }
    }
    assert!(
        last.is_some_and(|power| power >= 14),
        "the certified cutoff regressed"
    );
    let tangent = near_tangent(q(0), &mut dec);
    let (result, _) = with_operation_budget(ALLOWANCE, || {
        let spine = Spine::prove(&tangent, &mut dec)?;
        Sweep::new(&spine, &mut dec)?;
        Ok::<_, catallax_numbers::Refusal>(())
    });
    let refusal = result.unwrap_err();
    assert!(refusal.code.starts_with("field.spline-"), "{refusal:?}");
}
