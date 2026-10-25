//! Non-gating measurement of requested point-distance work and latency on
//! complete cylinder bands and caps: direct analytic reduction against the
//! general complete search on the same face.

use catallax_measure::distance::PreparedFace;
use catallax_numbers::{Ball, Decider, Rational};
use catallax_shapes::{
    Point3, Pose, RPoint3,
    distance::{self, Domain, Membership, Options},
};
use catallax_topology::Model;
use std::time::Instant;

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn point(x: i64, y: i64, z: i64) -> Point3 {
    RPoint3::new(r(x), r(y), r(z)).into()
}

fn main() {
    if std::env::var_os("CATALLAX_SPHERICAL_DISTANCE_BENCH").is_some() {
        spherical_queries();
        return;
    }
    let mut model = Model::new();
    let body = catallax_topology::cylinder(&mut model, r(2), r(3), &Pose::identity()).unwrap();
    let options = Options::default();
    let mut dec = Decider::default();
    let bounds: [Ball; 3] = [(-2, 2), (-2, 2), (0, 3)]
        .map(|(a, b)| Ball::from_rational(&r(a), 128).hull(&Ball::from_rational(&r(b), 128)));
    for (i, face) in model.faces_of_body(body).into_iter().enumerate() {
        let surface = model.geom().surface(model.face(face).surface).clone();
        // Prepared before the closure below borrows the model mutably.
        let prepared = PreparedFace::new(&model, face, &mut dec).unwrap();
        for (name, at) in [
            ("side", point(5, 0, 1)),
            ("rim", point(5, 0, 7)),
            ("axis", point(0, 0, 7)),
            ("cap", point(1, 1, -4)),
        ] {
            let start = Instant::now();
            let direct = prepared
                .point_to_face(&mut model, &at, &options, &mut dec)
                .unwrap();
            let direct_ns = start.elapsed().as_nanos();
            // The general search, with the exact trim of this face, as a reference.
            let domain = Domain::new(&surface, bounds.clone(), &mut dec).unwrap();
            let start = Instant::now();
            let mut inside = |p: &Point3, _: &[Ball; 3], d: &mut Decider| {
                Ok(Membership {
                    witness: catallax_boolean::face_contains(&mut model, face, p, d)?,
                    excluded: false,
                })
            };
            let general = distance::point_distance(&domain, &at, &options, &mut inside, &mut dec);
            let general_ns = start.elapsed().as_nanos();
            println!(
                "face={i} query={name} direct_cells={} direct_ns={direct_ns} converged={} general_cells={:?} general_ns={general_ns}",
                direct.cells,
                direct.converged,
                general.map(|g| g.cells)
            );
        }
    }
}

// This mode uses PreparedFace unchanged on both builds, so the baseline measures
// its actual boundary-excluding fallback rather than a weaker substitute adapter.
fn spherical_queries() {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Circle, Curve, GeomKey, RVector3, Sphere, Surface};
    use catallax_topology::{Coedge, Edge, Face, Loop, Sense};
    let selected = std::env::var("CATALLAX_DISTANCE_CASE").ok();
    for (name, rims) in [
        ("cap", vec![(3, Sense::Same)]),
        ("band", vec![(-3, Sense::Same), (3, Sense::Reversed)]),
    ] {
        let mut model = Model::new();
        let mut build_dec = Decider::default();
        let vector = |x, y, z| RVector3::new(r(x), r(y), r(z)).to_real();
        let surface = model.geom_mut().add_surface(Surface::Sphere(
            Sphere::new(
                point(0, 0, 0),
                r(5).into(),
                vector(0, 0, 1),
                vector(1, 0, 0),
                &mut build_dec,
            )
            .unwrap(),
        ));
        let loops = rims
            .into_iter()
            .map(|(z, sense)| {
                let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                    center: point(0, 0, z),
                    normal: vector(0, 0, 1),
                    radius_sq: r(25 - z * z).into(),
                    ref_dir: vector(1, 0, 0),
                }));
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Curve(curve),
                    b: GeomKey::Surface(surface),
                    at: None,
                    origin: FactOrigin::Construction("public-spherical-distance-bench"),
                });
                let edge = model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: None,
                    end: None,
                });
                let c = model.add_coedge(Coedge { edge, sense });
                model.add_loop(Loop::Coedges(vec![c]))
            })
            .collect();
        let face = model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops,
        });
        let options = Options {
            width: Rational::from_ratio(1, 1000).unwrap(),
            max_cells: 4096,
        };
        let start = Instant::now();
        let prepared = PreparedFace::new(&model, face, &mut build_dec).unwrap();
        let prepare_ns = start.elapsed().as_nanos();
        for (query, at) in [
            ("radial", point(6, 0, 8)),
            ("rim", point(7, 0, -8)),
            ("axis", point(0, 0, -5)),
        ] {
            let case = format!("{name}-{query}");
            if selected.as_ref().is_some_and(|s| s != &case) {
                continue;
            }
            // A fresh decider for each query: no prior query/root-cache advantage.
            let mut dec = Decider::default();
            let start = Instant::now();
            let answer = prepared.point_to_face(&mut model, &at, &options, &mut dec);
            let query_ns = start.elapsed().as_nanos();
            match answer {
                Ok(b) => println!(
                    "case={case} prepare_ns={prepare_ns} prepare_work={} query_ns={query_ns} query_work={} cells={} converged={} lower={} upper={} width={} max_cells={} max_work={} max_bits={}",
                    build_dec.stats().work,
                    dec.stats().work,
                    b.cells,
                    b.converged,
                    b.lower,
                    b.upper.map_or("none".into(), |q| q.to_string()),
                    options.width,
                    options.max_cells,
                    dec.policy().max_work,
                    dec.policy().max_precision_bits
                ),
                Err(e) => println!(
                    "case={case} prepare_ns={prepare_ns} prepare_work={} query_ns={query_ns} query_work={} refusal={} width={} max_cells={} max_work={} max_bits={}",
                    build_dec.stats().work,
                    dec.stats().work,
                    e.code,
                    options.width,
                    options.max_cells,
                    dec.policy().max_work,
                    dec.policy().max_precision_bits
                ),
            }
        }
    }
}
