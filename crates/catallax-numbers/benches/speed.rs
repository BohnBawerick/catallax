#![allow(missing_docs)]
#[path = "../tests/common/mod.rs"]
mod common;
use catallax_numbers::*;
use common::*;
use std::collections::BTreeMap;
use std::hint::black_box;
use std::time::{Duration, Instant};

fn decisions(label: &str, values: Vec<Real>) {
    let mut d = Decider::default();
    let mut paths: BTreeMap<&str, (u64, Duration)> = BTreeMap::new();
    let mut filter = 0;
    let mut max_bits = 0;
    let mut work = 0;
    for value in &values {
        // Prime the public approximation cache outside the measured decision.
        black_box(value.approx());
        let start = Instant::now();
        let _ = black_box(d.sign(black_box(value)));
        let elapsed = start.elapsed();
        let s = d.take_stats();
        filter += s.filter;
        max_bits = max_bits.max(s.max_bits);
        work += s.work;
        let path = if s.filter > 0 {
            "filter"
        } else if s.refined > 0 {
            "refined"
        } else if s.exact > 0 {
            "exact"
        } else if s.coincident_exact > 0 {
            "coincident_exact"
        } else if s.coincident_separation > 0 {
            "coincident_separation"
        } else if s.not_exact > 0 {
            "not_exact"
        } else {
            "refused"
        };
        let entry = paths.entry(path).or_default();
        entry.0 += 1;
        entry.1 += elapsed;
    }
    println!(
        "decisions family={label} count={} filter_share={} max_bits={max_bits} work={work}",
        values.len(),
        filter as f64 / values.len() as f64
    );
    for (path, (count, time)) in paths {
        println!(
            "path family={label} path={path} count={count} ns_per_decision={}",
            time.as_nanos() / u128::from(count)
        );
    }
}
fn orient2(p: [Real; 2], q: [Real; 2], r: [Real; 2]) -> Real {
    (&q[0] - &p[0]) * (&r[1] - &p[1]) - (&q[1] - &p[1]) * (&r[0] - &p[0])
}
fn orient3(a: [Real; 3], b: [Real; 3], c: [Real; 3], d: [Real; 3]) -> Real {
    let b: Vec<_> = (0..3).map(|i| &b[i] - &a[i]).collect();
    let c: Vec<_> = (0..3).map(|i| &c[i] - &a[i]).collect();
    let d: Vec<_> = (0..3).map(|i| &d[i] - &a[i]).collect();
    &b[0] * &(&c[1] * &d[2] - &c[2] * &d[1]) - &b[1] * &(&c[0] * &d[2] - &c[2] * &d[0])
        + &b[2] * &(&c[0] * &d[1] - &c[1] * &d[0])
}
fn main() {
    if std::env::var_os("CATALLAX_NUMBERS_RUN_BENCH").as_deref() != Some(std::ffi::OsStr::new("1"))
    {
        println!(
            "Disabled against the interface skeleton. Set CATALLAX_NUMBERS_RUN_BENCH=1 to measure piece 1."
        );
        return;
    }
    let mut rng = Rng(0xbeac4);
    let mut two = Vec::new();
    let mut three = Vec::new();
    for _ in 0..10_000 {
        let mut coord = || Real::from_i64((rng.next() % 2_000_001) as i64 - 1_000_000);
        two.push(orient2(
            std::array::from_fn(|_| coord()),
            std::array::from_fn(|_| coord()),
            std::array::from_fn(|_| coord()),
        ));
        three.push(orient3(
            std::array::from_fn(|_| coord()),
            std::array::from_fn(|_| coord()),
            std::array::from_fn(|_| coord()),
            std::array::from_fn(|_| coord()),
        ));
    }
    decisions("general_position_2d", two);
    decisions("general_position_3d", three);
    let mut near = Vec::new();
    for i in 0..256 {
        for j in 0..256 {
            let point =
                |v: i64| Real::from_rational(r("1/2") + Rational::from_i64(v).mul_pow2(-53));
            near.push(orient2(
                [point(i), point(j)],
                [real("12"), real("12")],
                [real("24"), real("24")],
            ));
        }
    }
    decisions("near_orientation", near);
    for n in ["1000000", "100000000000000000000"] {
        let n = real(n);
        decisions(
            "concave_sqrt",
            vec![(&n + &real("1")).sqrt() + (&n - &real("1")).sqrt() - real("2") * n.sqrt()],
        );
    }
    decisions("near_sqrt", vec![near_tie()]);
    decisions("hidden_zeros", hidden_zeros());
    let mut products = Vec::new();
    for _ in 0..128 {
        let x = eq(&rng.rational().abs());
        let y = eq(&rng.rational().abs());
        products.push(
            Real::from_rational(x.clone()).sqrt() * Real::from_rational(y.clone()).sqrt()
                - Real::from_rational(x * y).sqrt(),
        );
    }
    decisions("sqrt_product_zeros", products);
    for degree in [10usize, 20, 21, 40, 60, 80, 100] {
        for (family, p) in [
            (
                "wilkinson",
                UPoly::from_roots(
                    &(1..=degree)
                        .map(|i| Rational::from_i64(i as i64))
                        .collect::<Vec<_>>(),
                ),
            ),
            ("mignotte", mignotte(degree)),
            ("chebyshev", chebyshev(degree)),
        ] {
            let mut d = Decider::default();
            let start = Instant::now();
            let outcome = d.real_roots(black_box(&p));
            let time = start.elapsed();
            let status = match &outcome {
                Outcome::Proven(v) => format!("roots:{}", v.len()),
                Outcome::Refused(e) => e.code.into(),
                other => format!("{other:?}"),
            };
            let _ = black_box(outcome);
            println!(
                "roots family={family} degree={degree} ns={} max_bits={} work={} outcome={status}",
                time.as_nanos(),
                d.stats().max_bits,
                d.stats().work
            );
        }
    }
    for positive in [false, true] {
        let mut d = Decider::default();
        let system = sphere_system();
        let candidate = sphere_box(positive);
        let start = Instant::now();
        let outcome = d.certify_root(black_box(&system), black_box(&candidate));
        println!(
            "certify spheres positive={positive} ns={} max_bits={} outcome={outcome:?}",
            start.elapsed().as_nanos(),
            d.stats().max_bits
        );
    }
    for (digits, bits) in [(1000, 3322), (10000, 33220)] {
        let start = Instant::now();
        let pi = Ball::pi(black_box(bits));
        let elapsed = start.elapsed();
        println!(
            "pi digits={digits} bits={bits} ns={} radius={}",
            elapsed.as_nanos(),
            pi.rad()
        );
        black_box(pi);
    }
    let mut facts = FactBook::new();
    for i in 0u32..100_000 {
        facts.record(Fact {
            kind: FactKind::Coaxial,
            a: i,
            b: i + 1,
            at: None,
            origin: FactOrigin::Construction("benchmark chain"),
        });
    }
    let start = Instant::now();
    for i in 0u32..1000 {
        black_box(facts.holds(FactKind::Coaxial, black_box(i), black_box(i + 1000)));
    }
    println!(
        "facts count=100000 queries=1000 chain=1000 ns_per_query={}",
        start.elapsed().as_nanos() / 1000
    );
}
