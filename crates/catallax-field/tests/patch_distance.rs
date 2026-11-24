//! 4.2a step 7.1: certified point-to-patch distance on the M11 support pairs.
//!
//! The supports are those of `catallax-shapes/tests/spline_fillet.rs`: the M27
//! loft edge's bilinear pair and the biquadratic edge's pair. Each query's
//! enclosure is checked against an independent binary64 search (a dense grid,
//! then zoomed grids round the best sample, as SMLib's `FindRoughGuess` does),
//! which shares no code with `PatchDistance`.

use std::time::{Duration, Instant};

use catallax_field::patch_distance::{Enclosure, PatchDistance, WIDTH};
use catallax_numbers::{Decider, MPoly, Rational, Real};
use catallax_shapes::{BSplineSurface, Point3, RPoint3, RationalPatch};

type Q = Rational;

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

/// The supports of the M27 loft edge (`biquadratic = false`) or of the
/// biquadratic edge, as in `spline_fillet.rs`.
fn supports(biquadratic: bool, dec: &mut Decider) -> [RationalPatch; 2] {
    let k = |n: i64| MPoly::constant(2, q(n));
    let (u, v) = (MPoly::var(2, 0), MPoly::var(2, 1));
    let one = k(1);
    let mut a = [
        &(&k(-1) + &u.scale(&q(2))) + &(&v * &(&one - &u)),
        &k(-1) + &(&u * &v),
        v.scale(&q(2)),
    ];
    let mut b = [
        &one - &(&u * &v),
        &(&k(-1) + &u.scale(&q(2))) + &v,
        v.scale(&q(2)),
    ];
    let degree = if biquadratic {
        let bump = (&(&(&u * &(&one - &u)) * &v) * &(&one - &v)).scale(&r("1/5"));
        a[2] = &a[2] + &bump;
        b[2] = &b[2] - &bump;
        [2, 2]
    } else {
        [1, 1]
    };
    [patch(&a, degree, dec), patch(&b, degree, dec)]
}

/// The same supports in binary64: `which` 0 is the first, 1 the second.
fn point(biquadratic: bool, which: usize, u: f64, v: f64) -> [f64; 3] {
    let bump = if biquadratic {
        0.2 * u * (1.0 - u) * v * (1.0 - v)
    } else {
        0.0
    };
    if which == 0 {
        [-1.0 + 2.0 * u + v * (1.0 - u), -1.0 + u * v, 2.0 * v + bump]
    } else {
        [1.0 - u * v, -1.0 + 2.0 * u + v, 2.0 * v - bump]
    }
}

fn squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| (a[i] - b[i]) * (a[i] - b[i])).sum()
}

/// The third check: the smallest squared distance a binary64 search finds.
/// A 64 by 64 grid, then twenty 9 by 9 grids, each a quarter as wide, round
/// the best sample so far.
fn sampled(biquadratic: bool, which: usize, p: [f64; 3]) -> f64 {
    let at = |u: f64, v: f64| squared(point(biquadratic, which, u, v), p);
    let mut best = (f64::INFINITY, 0.0, 0.0);
    for i in 0..=64 {
        for j in 0..=64 {
            let (u, v) = (f64::from(i) / 64.0, f64::from(j) / 64.0);
            let value = at(u, v);
            if value < best.0 {
                best = (value, u, v);
            }
        }
    }
    let mut half = 1.0 / 64.0;
    for _ in 0..20 {
        let (_, cu, cv) = best;
        for i in -4..=4 {
            for j in -4..=4 {
                let u = (cu + f64::from(i) * half / 4.0).clamp(0.0, 1.0);
                let v = (cv + f64::from(j) * half / 4.0).clamp(0.0, 1.0);
                let value = at(u, v);
                if value < best.0 {
                    best = (value, u, v);
                }
            }
        }
        half /= 4.0;
    }
    best.0
}

/// `count` query points with dyadic coordinates: the first half within 0.3 of
/// the shared edge `(1, v - 1, 2 v)`, the rest in the box round both supports.
/// A fixed linear congruential sequence (Knuth's MMIX constants) draws them.
fn points(count: usize) -> Vec<RPoint3> {
    let mut state: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) as i64 % 4096
    };
    let unit = |n: i64| Q::from_ratio(n, 4096).unwrap();
    (0..count)
        .map(|index| {
            let (a, b, c) = (unit(next()), unit(next()), unit(next()));
            if index < count / 2 {
                let v = unit(next());
                let near = |t: Q| &(&t - &r("1/2")) * &r("3/5");
                RPoint3::new(
                    &q(1) + &near(a),
                    &(&v - &q(1)) + &near(b),
                    &(&v * &q(2)) + &near(c),
                )
            } else {
                let spread = |t: Q, lo: i64, size: i64| &q(lo) + &(&t * &q(size));
                RPoint3::new(spread(a, -2, 4), spread(b, -2, 4), spread(c, -1, 4))
            }
        })
        .collect()
}

/// One query's answer and its work, in a form that compares across runs.
#[derive(Clone, Debug, PartialEq)]
struct Answer {
    enclosure: Enclosure,
    work: u64,
}

fn query(field: &PatchDistance, p: &RPoint3) -> (Answer, Duration) {
    let mut dec = Decider::default();
    let start = Instant::now();
    let enclosure = field.distance(p, &mut dec).expect("a patch distance");
    let time = start.elapsed();
    let work = dec.stats().work;
    (Answer { enclosure, work }, time)
}

/// Every query on `workers` threads, each taking every `workers`-th point.
fn run(fields: &[PatchDistance], points: &[RPoint3], workers: usize) -> Vec<Answer> {
    let mut out: Vec<Option<Answer>> = vec![None; points.len() * fields.len()];
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|worker| {
                scope.spawn(move || {
                    let mut mine = Vec::new();
                    for (index, p) in points.iter().enumerate().skip(worker).step_by(workers) {
                        for (which, field) in fields.iter().enumerate() {
                            mine.push((index * fields.len() + which, query(field, p).0));
                        }
                    }
                    mine
                })
            })
            .collect();
        for handle in handles {
            for (slot, answer) in handle.join().unwrap() {
                out[slot] = Some(answer);
            }
        }
    });
    out.into_iter().map(Option::unwrap).collect()
}

/// Checks every enclosure against the binary64 search, and returns the
/// sorted query times and work.
fn check(biquadratic: bool, count: usize) -> (Vec<Duration>, Vec<u64>, usize) {
    let mut dec = Decider::default();
    let fields = supports(biquadratic, &mut dec).map(|s| PatchDistance::new(&s, &mut dec).unwrap());
    let points = points(count);
    let mut times = Vec::new();
    let mut work = Vec::new();
    let mut cells = 0;
    let slack = 1e-12;
    for p in &points {
        let fast = [&p.x, &p.y, &p.z].map(|x| x.to_f64_nearest());
        for (which, field) in fields.iter().enumerate() {
            let (answer, time) = query(field, p);
            let enclosure = answer.enclosure;
            assert!(
                enclosure.converged,
                "{p:?} on support {which}: {enclosure:?}"
            );
            assert!(enclosure.distance.hi() - enclosure.distance.lo() <= WIDTH);
            let found = sampled(biquadratic, which, fast);
            // The search's value is the squared distance to a point of the
            // patch, so at least the true minimum; it reaches the minimum to
            // within the slack when it finds the nearest point's basin.
            assert!(
                found >= enclosure.squared.lo() - slack,
                "{p:?} on support {which}: sampled {found} below {enclosure:?}"
            );
            assert!(
                found <= enclosure.squared.hi() + slack,
                "{p:?} on support {which}: sampled {found} above {enclosure:?}"
            );
            times.push(time);
            work.push(answer.work);
            cells = cells.max(enclosure.cells);
        }
    }
    times.sort();
    work.sort();
    (times, work, cells)
}

fn report(name: &str, (times, work, cells): (Vec<Duration>, Vec<u64>, usize)) {
    let total: Duration = times.iter().sum();
    let at = |q: usize| times[(times.len() - 1) * q / 100];
    println!(
        "{name}: {} queries, {:.1} s in all, median {:?}, 99th percentile {:?}, \
         slowest {:?}; work median {}, largest {}; at most {cells} cells",
        times.len(),
        total.as_secs_f64(),
        at(50),
        at(99),
        times[times.len() - 1],
        work[work.len() / 2],
        work[work.len() - 1],
    );
}

#[test]
fn patch_distances_enclose_the_sampled_minimum() {
    report("m27 pair", check(false, 200));
    report("biquadratic pair", check(true, 200));
}

#[test]
fn patch_distances_are_the_same_at_every_worker_count() {
    let mut dec = Decider::default();
    let fields = supports(true, &mut dec).map(|s| PatchDistance::new(&s, &mut dec).unwrap());
    let points = points(64);
    let one = run(&fields, &points, 1);
    assert_eq!(one, run(&fields, &points, 2));
    assert_eq!(one, run(&fields, &points, 4));
}

/// A point on a support is at distance zero, within the width.
#[test]
fn a_point_on_a_support_is_at_distance_zero() {
    let mut dec = Decider::default();
    let [first, _] = supports(false, &mut dec);
    let field = PatchDistance::new(&first, &mut dec).unwrap();
    // a(1/4, 1/2) = (-1 + 1/2 + 3/8, -1 + 1/8, 1).
    let p = RPoint3::new(r("-1/8"), r("-7/8"), q(1));
    let (answer, _) = query(&field, &p);
    assert!(answer.enclosure.converged);
    assert_eq!(
        answer.enclosure.distance.lo().abs().to_bits(),
        0.0_f64.to_bits()
    );
    assert!(answer.enclosure.distance.hi() <= WIDTH);
}

/// The step 7.1 spike of the decisions: 10^4 points on the M27 pair, timed.
/// Run in release: `cargo test --release -p catallax-field --test
/// patch_distance -- --ignored --nocapture`.
#[test]
#[ignore = "timing spike, run in release"]
fn spike_ten_thousand_points() {
    report("m27 pair", check(false, 10_000));
    report("biquadratic pair", check(true, 10_000));
}
