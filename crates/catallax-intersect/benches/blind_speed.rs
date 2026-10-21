//! Non-gating piece 4 timings, enabled explicitly after implementation.
#[path = "../tests/common/mod.rs"]
mod common;
use common::curves::*;
use common::*;
use std::{hint::black_box, time::Instant};

fn report(label: &str, start: Instant, d: &Decider) {
    let stats = d.stats();
    let total = stats.decisions.max(1) as f64;
    println!(
        "{label}: ns={} decisions={} fact_share={} filter_share={} max_precision_bits={} work={}",
        start.elapsed().as_nanos(),
        stats.decisions,
        stats.coincident_fact as f64 / total,
        stats.filter as f64 / total,
        stats.max_bits,
        stats.work
    );
}
fn main() {
    // cargo test --all-targets must remain green while the implementation is todo!().
    if std::env::var_os("CATALLAX_INTERSECT_BENCH").as_deref() != Some(std::ffi::OsStr::new("1")) {
        println!("Set CATALLAX_INTERSECT_BENCH=1 to run the blind piece 4 timing suite.");
        return;
    }
    let mut rng = Rng(0x04be11);
    for case in cases::pairs() {
        for pose in [Pose::identity(), rng.pose()] {
            let mut g = GeomStore::new();
            let a = g.add_surface(case.a.transformed(&pose));
            let b = g.add_surface(case.b.transformed(&pose));
            let mut d = Decider::default();
            let start = Instant::now();
            black_box(intersect_surfaces(&mut g, a, b, &mut d)).ok();
            report(&format!("surface {} pose={pose:?}", case.name), start, &d);
        }
    }
    for case in hits() {
        let mut g = GeomStore::new();
        let c = case.curve.add(&mut g);
        let s = g.add_surface(case.surface);
        let mut d = Decider::default();
        let start = Instant::now();
        black_box(intersect_curve_surface(&mut g, c, s, &mut d)).ok();
        report(
            &format!(
                "curve-surface {} curve={:?} surface={:?}",
                case.name,
                g.curve(c).kind(),
                g.surface(s).kind()
            ),
            start,
            &d,
        );
    }
    for carrier in [
        Carrier::Analytic(line28()),
        Carrier::Analytic(k5()),
        Carrier::Analytic(ellipse4()),
        viviani(true),
        loop25(),
    ] {
        let mut g = GeomStore::new();
        let c = carrier.add(&mut g);
        let point = g.curve(c).parameter_origin(&g);
        let p = match g.curve(c) {
            Curve::Intersection(branch) => branch.witness,
            _ => explicit(
                &mut g,
                rational_point(&point).expect("these analytic parameter origins are rational"),
            ),
        };
        if let Curve::Intersection(branch) = g.curve(c) {
            let surfaces = branch.surfaces;
            for s in surfaces {
                let mut d = Decider::default();
                let start = Instant::now();
                let _ = black_box(intersect_curve_surface(&mut g, c, s, &mut d));
                report("curve-surface branch own support (fact)", start, &d);
            }
        }
        let mut d = Decider::default();
        let start = Instant::now();
        let _ = black_box(on_curve(&g, c, p, &mut d));
        report(&format!("on_curve {:?}", g.curve(c).kind()), start, &d);
        let mut d = Decider::default();
        let start = Instant::now();
        let _ = black_box(same_curve(&g, c, c, &mut d));
        report(&format!("same_curve {:?}", g.curve(c).kind()), start, &d);
        let (from, to) = if matches!(g.curve(c), Curve::Line(_)) {
            let q = explicit(&mut g, pi([5, 0, 0]));
            let mut d = Decider::default();
            let start = Instant::now();
            let _ = black_box(order_on_curve(&g, c, p, q, &mut d));
            report("order_on_curve Line", start, &d);
            (Some(p), Some(q))
        } else {
            let mut d = Decider::default();
            let start = Instant::now();
            let _ = black_box(cyclic_order(&g, c, [p, p, p], &mut d));
            report(&format!("cyclic_order {:?}", g.curve(c).kind()), start, &d);
            (None, None)
        };
        let mut d = Decider::default();
        let start = Instant::now();
        let result = trace(&g, c, from, to, &mut d);
        report(&format!("trace {:?}", g.curve(c).kind()), start, &d);
        if let Ok(trace) = result {
            for bits in [64, 1000] {
                for order in [4, 16] {
                    let start = Instant::now();
                    let mut count = 0;
                    for piece in trace.value.pieces() {
                        let (a, b) = piece.range();
                        let t = Ball::from_rational(&((&a + &b) / n(2)), bits);
                        black_box(piece.taylor(&t, order, &Decider::default()).unwrap());
                        count += 1;
                    }
                    println!(
                        "taylor {:?}: order={order} bits={bits} pieces={count} ns={}",
                        g.curve(c).kind(),
                        start.elapsed().as_nanos()
                    );
                }
            }
        }
    }
}
