//! A tapered bore exercises the reversed cone face and both periodic seam uses.
#[allow(dead_code)]
#[path = "../../catallax-topology/tests/common/periodic.rs"]
mod fixture;

use catallax_numbers::{Ball, Decider, Real};
use catallax_shapes::{Cone, Curve, Point3, Surface, Vector3};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;
use fixture::q;

fn tapered_bore() -> n::NumericalDefinition {
    let mut d = fixture::bore();
    d.name = "authored-conical-bore".into();
    d.vertices[9].position = [q("1.5"), q("0"), q("1")];
    for edge in &mut d.edges {
        if edge.key == "rim-1" {
            let Curve::Circle(circle) = &mut edge.curve else {
                unreachable!()
            };
            circle.radius_sq = q("2.25").into();
        }
        if edge.vertices == [8, 9] {
            let Curve::Line(line) = &mut edge.curve else {
                unreachable!()
            };
            line.dir = Vector3::new(q("0.5").into(), q("0").into(), q("1").into());
        }
    }
    for face in &mut d.faces {
        if matches!(face.surface, Surface::Cylinder(_)) {
            face.surface = Surface::Cone(Cone {
                apex: Point3::new(q("0").into(), q("0").into(), q("-2").into()),
                axis: Vector3::new(q("0").into(), q("0").into(), q("1").into()),
                slope: q("0.5").into(),
                ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
            });
            for uv in face.loops.iter_mut().flatten().flat_map(|t| &mut t.uv) {
                uv[1] = &uv[1] + &q("2");
            }
        } else {
            for trim in face.loops.iter_mut().flatten() {
                if d.edges[trim.edge].key != "rim-1" {
                    continue;
                }
                let p = d.vertices[9]
                    .position
                    .each_ref()
                    .map(|q| Ball::from_rational(q, 128));
                let uv = face
                    .surface
                    .parameters_of(&p, &Decider::default())
                    .unwrap()
                    .unwrap()
                    .map(|v| v.mid());
                trim.uv = [uv.clone(), uv];
            }
        }
    }
    d
}

#[test]
fn inner_cone_seams_export_from_every_cyclic_loop_start() {
    for start in 0..4 {
        let mut d = tapered_bore();
        d.faces.last_mut().unwrap().loops[0].rotate_left(start);
        let body = n::build(d, n::NumericalOptions::default(), &mut Decider::default()).unwrap();
        let file = output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(
            std::str::from_utf8(file.bytes())
                .unwrap()
                .contains("CONICAL_SURFACE(")
        );
        assert!(
            file.report()
                .checks
                .iter()
                .any(|c| c == "file.cone-step-v-generator-distance")
        );
        file.check_written(file.bytes(), &mut Decider::default())
            .unwrap();
        if let Ok(directory) = std::env::var("CATALLAX_CONE_EVIDENCE") {
            std::fs::create_dir_all(&directory).unwrap();
            std::fs::write(
                std::path::Path::new(&directory).join(format!("cone-bore-loop-{start}.step")),
                file.bytes(),
            )
            .unwrap();
        }
    }
}
