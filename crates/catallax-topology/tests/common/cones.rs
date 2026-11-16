//! Authored conical bores with positive and signed negative axial charts.
#[allow(dead_code)]
#[path = "periodic.rs"]
mod periodic;
use catallax_numbers::{Ball, Decider, Real};
use catallax_shapes::{Cone, Curve, Point3, Surface, Vector3};
use catallax_topology::numerical as n;
use periodic::q;

pub fn positive_bore() -> n::NumericalDefinition {
    tapered_bore(false)
}
pub fn negative_height_bore() -> n::NumericalDefinition {
    tapered_bore(true)
}

fn tapered_bore(negative: bool) -> n::NumericalDefinition {
    let top = q(if negative { "0.5" } else { "1.5" });
    let apex = q(if negative { "2" } else { "-2" });
    let mut d = periodic::bore();
    d.name = "authored-conical-bore".into();
    d.vertices[9].position = [top.clone(), q("0"), q("1")];
    for edge in &mut d.edges {
        if edge.key == "rim-1" {
            let Curve::Circle(circle) = &mut edge.curve else {
                unreachable!()
            };
            circle.radius_sq = (&top * &top).into();
        }
        if edge.vertices == [8, 9] {
            let Curve::Line(line) = &mut edge.curve else {
                unreachable!()
            };
            line.dir = Vector3::new((&top - &q("1")).into(), q("0").into(), q("1").into());
        }
    }
    for face in &mut d.faces {
        if matches!(face.surface, Surface::Cylinder(_)) {
            face.surface = Surface::Cone(Cone {
                apex: Point3::new(q("0").into(), q("0").into(), apex.clone().into()),
                axis: Vector3::new(q("0").into(), q("0").into(), q("1").into()),
                slope: q("0.5").into(),
                ref_dir: Vector3::new(
                    q(if negative { "-1" } else { "1" }).into(),
                    Real::zero(),
                    Real::zero(),
                ),
            });
            for uv in face.loops.iter_mut().flatten().flat_map(|t| &mut t.uv) {
                uv[1] = &uv[1] - &apex;
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
