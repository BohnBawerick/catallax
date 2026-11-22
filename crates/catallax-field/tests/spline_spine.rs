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
use std::time::Instant;

use catallax_field::spline_spine::{Spine, SpineDefinition, TUBE_CAP};
use catallax_numbers::{Decider, MPoly, Rational, Real, Sign, with_operation_budget};
use catallax_shapes::{BSplineSurface, Point3, RationalPatch};

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
    assert!(detail.contains("the tube reached"), "{detail}");
    println!("{detail}");
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
