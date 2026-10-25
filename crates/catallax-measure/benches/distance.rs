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
