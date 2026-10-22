#![allow(missing_docs)]
#[path = "../tests/common/mod.rs"]
mod common;
use catallax_numbers::*;
use catallax_shapes::*;
use common::cases::*;
use common::*;
use std::hint::black_box;
use std::time::Instant;

fn measure<T>(label: &str, n: u64, mut f: impl FnMut() -> T) {
    let start = Instant::now();
    for _ in 0..n {
        black_box(f());
    }
    println!(
        "operation={label} samples={n} ns_per_call={}",
        start.elapsed().as_nanos() / u128::from(n)
    );
}
fn decisions<T: std::fmt::Debug>(
    label: &str,
    n: u64,
    mut f: impl FnMut(&mut Decider) -> Outcome<T>,
) {
    let mut dec = Decider::default();
    let mut outcomes = [0u64; 4];
    let start = Instant::now();
    for _ in 0..n {
        match black_box(f(&mut dec)) {
            Outcome::Proven(_) => outcomes[0] += 1,
            Outcome::ProvenCoincident(_) => outcomes[1] += 1,
            Outcome::NotExact(_) => outcomes[2] += 1,
            Outcome::Refused(e) => {
                outcomes[3] += 1;
                println!(
                    "refusal operation={label} code={} detail={}",
                    e.code, e.detail
                );
            }
        }
    }
    let elapsed = start.elapsed();
    let s = dec.stats();
    let divisor = s.decisions.max(1) as f64;
    println!(
        "operation={label} samples={n} ns_per_call={} proven={} coincident={} not_exact={} refused={} decisions={} filter_share={} fact_share={} max_bits={} work={}",
        elapsed.as_nanos() / u128::from(n),
        outcomes[0],
        outcomes[1],
        outcomes[2],
        outcomes[3],
        s.decisions,
        s.filter as f64 / divisor,
        s.coincident_fact as f64 / divisor,
        s.max_bits,
        s.work
    );
}
fn main() {
    if std::env::var_os("CATALLAX_SHAPES_RUN_BENCH").as_deref() != Some(std::ffi::OsStr::new("1")) {
        println!(
            "Disabled against interface bodies. Set CATALLAX_SHAPES_RUN_BENCH=1 to measure piece 2."
        );
        return;
    }
    let mut rng = Rng(0x2beac4);
    let pose = rng.pose();
    for surface in surfaces() {
        let kind = surface.kind();
        let mut store = GeomStore::new();
        let s = store.add_surface(surface.clone());
        let ids: Vec<_> = (0..1000)
            .map(|_| explicit(&mut store, rp(std::array::from_fn(|_| rng.rational()))))
            .collect();
        let mut i = 0;
        decisions(&format!("side.general.{kind:?}"), 1000, |d| {
            let p = ids[i % ids.len()];
            i += 1;
            store.side(s, p, d)
        });
        decisions(&format!("normals.{kind:?}"), 1000, |d| {
            store.normals_agree(s, s, ids[0], d)
        });
        decisions(&format!("same_point.{kind:?}"), 1000, |d| {
            store.same_point(ids[0], ids[1], d)
        });
        measure(&format!("implicit.{kind:?}"), 1000, || surface.implicit());
        measure(&format!("transformed.{kind:?}"), 1000, || {
            surface.transformed(black_box(&pose))
        });
        for bits in [64, 1000] {
            let u = ball("1", bits);
            let v = ball("1/2", bits);
            let point = surface
                .evaluate(&u, &v, &catallax_numbers::Decider::default())
                .unwrap();
            measure(&format!("evaluate.{kind:?}.{bits}"), 32, || {
                surface
                    .evaluate(
                        black_box(&u),
                        black_box(&v),
                        &catallax_numbers::Decider::default(),
                    )
                    .unwrap()
            });
            measure(&format!("derivatives.{kind:?}.{bits}"), 32, || {
                surface
                    .derivatives(
                        black_box(&u),
                        black_box(&v),
                        &catallax_numbers::Decider::default(),
                    )
                    .unwrap()
            });
            let mut none = 0;
            measure(&format!("parameters_of.{kind:?}.{bits}"), 32, || {
                let out = surface
                    .parameters_of(black_box(&point), &catallax_numbers::Decider::default())
                    .unwrap();
                none += u64::from(out.is_none());
                out
            });
            println!("parameters_of kind={kind:?} bits={bits} none={none}");
        }
    }
    for (surface, point, expected) in side_cases().into_iter().skip(11).take(10) {
        let mut store = GeomStore::new();
        let kind = surface.kind();
        let s = store.add_surface(surface);
        let p = explicit(&mut store, point);
        decisions(&format!("side.near_1e-30.{kind:?}.{expected}"), 1000, |d| {
            store.side(s, p, d)
        });
    }
    let mut store = GeomStore::new();
    let s = store.add_surface(cylinder("1"));
    let x = real("1/2").sqrt();
    for (label, offset) in [("on", r("0")), ("front", r("1e-30")), ("back", -r("1e-30"))] {
        let p = Point3::new(&x + &Real::from_rational(offset), x.clone(), real("0"));
        decisions(&format!("side.square_root.{label}"), 32, |d| {
            store.side_of(s, &p, d)
        });
    }
    let root = sphere_root(true);
    let spheres: Vec<_> = [p("0", "0", "0"), p("6", "0", "0"), p("0", "6", "0")]
        .into_iter()
        .map(|center| {
            store.add_surface(Surface::Sphere(
                Sphere::new(
                    (center).into(),
                    (r("5")).into(),
                    (v("0", "0", "1")).into(),
                    (v("1", "0", "0")).into(),
                    &mut catallax_numbers::Decider::default(),
                )
                .unwrap(),
            ))
        })
        .collect();
    let root_id = store.add_point(PointDef::Meet(Meet {
        of: MeetOf::Surfaces(spheres.clone()),
        position: Position::Isolated(root),
    }));
    let point = store.position(root_id);
    for s in &spheres {
        decisions("side.certified_root", 32, |d| store.side_of(*s, &point, d));
    }
    decisions("side.fact", 10000, |d| store.side(spheres[0], root_id, d));
    let explicit = explicit(&mut store, p("8", "9", "10"));
    let mut counter = 0;
    decisions("side.mixed", 10000, |d| {
        counter += 1;
        store.side(
            spheres[0],
            if counter % 2 == 0 { root_id } else { explicit },
            d,
        )
    });
    for curve in curves() {
        let kind = curve.kind();
        let c = store.add_curve(curve.clone());
        decisions(&format!("point_on_curve.{kind:?}"), 1000, |d| {
            store.point_on_curve(c, explicit, d)
        });
        decisions(&format!("same_point.{kind:?}"), 1000, |d| {
            store.same_point(root_id, explicit, d)
        });
        measure(&format!("curve.transformed.{kind:?}"), 1000, || {
            curve.transformed(&pose)
        });
        for bits in [64, 1000] {
            let t = ball("1", bits);
            measure(&format!("curve.evaluate.{kind:?}.{bits}"), 32, || {
                curve
                    .evaluate(&t, &catallax_numbers::Decider::default())
                    .unwrap()
            });
        }
    }
    let c = store.add_curve(Curve::Intersection(IntersectionCurve {
        surfaces: [spheres[0], spheres[1]],
        witness: root_id,
    }));
    decisions("point_on_curve.Intersection", 1000, |d| {
        store.point_on_curve(c, explicit, d)
    });
    decisions("same_point.Intersection", 1000, |d| {
        store.same_point(root_id, explicit, d)
    });
    for bits in [64, 1000] {
        let width = r("1").mul_pow2(-bits);
        let mut d = Decider::default();
        let mut refused = 0;
        measure(&format!("enclose_point.2^-{bits}"), 16, || {
            let out = store.enclose_point(root_id, &width, &mut d);
            if let Err(e) = &out {
                refused += 1;
                println!("enclose_refusal code={} detail={}", e.code, e.detail);
            }
            out
        });
        println!(
            "enclose bits={bits} refused={refused} stats={:?}",
            d.stats()
        );
    }
    measure("Rotation.from_quaternion", 1000, || {
        Rotation::from_quaternion(r("1"), r("2"), r("3"), r("4"))
    });
}
