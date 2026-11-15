//! Finite cylindrical pcurves with source vertices slightly off their carriers.
#[path = "../../catallax-topology/tests/common/analytic.rs"]
mod fixture;

use catallax_numbers::{Ball, Decider, Real};
use catallax_shapes::{Curve, Point3, RVector3, Surface, Vector3};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;
use fixture::q;

fn rotated(mut d: n::NumericalDefinition) -> n::NumericalDefinition {
    let c = q("20/29");
    let s = q("21/29");
    let rotate = |x: &Real, y: &Real| {
        [
            &(x * &Real::from(c.clone())) - &(y * &Real::from(s.clone())),
            &(x * &Real::from(s.clone())) + &(y * &Real::from(c.clone())),
        ]
    };
    let point = |p: &Point3| {
        let [x, y] = rotate(&p.x, &p.y);
        Point3::new(x, y, p.z.clone())
    };
    let vector = |v: &Vector3| {
        let [x, y] = rotate(&v.x, &v.y);
        Vector3::new(x, y, v.z.clone())
    };
    for vertex in &mut d.vertices {
        let [x, y, _] = vertex.position.clone();
        vertex.position[0] = &(&x * &c) - &(&y * &s);
        vertex.position[1] = &(&x * &s) + &(&y * &c);
    }
    for edge in &mut d.edges {
        match &mut edge.curve {
            Curve::Line(line) => {
                line.origin = point(&line.origin);
                line.dir = vector(&line.dir);
            }
            Curve::Circle(circle) => {
                circle.center = point(&circle.center);
                circle.normal = vector(&circle.normal);
                circle.ref_dir = vector(&circle.ref_dir);
            }
            _ => unreachable!(),
        }
    }
    for face in &mut d.faces {
        match &mut face.surface {
            Surface::Plane(plane) => {
                plane.origin = point(&plane.origin);
                plane.normal = vector(&plane.normal);
                plane.ref_dir = vector(&plane.ref_dir);
            }
            Surface::Cylinder(cylinder) => {
                cylinder.origin = point(&cylinder.origin);
                cylinder.axis = vector(&cylinder.axis);
                cylinder.ref_dir = vector(&cylinder.ref_dir);
            }
            _ => unreachable!(),
        }
    }
    d
}

fn chart(phase: usize) -> n::NumericalDefinition {
    let mut d = fixture::quarter();
    let surface = d.faces.last_mut().unwrap();
    let Surface::Cylinder(cylinder) = &mut surface.surface else {
        unreachable!()
    };
    let reference = [[1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0]][phase];
    cylinder.ref_dir = RVector3::new(
        reference[0].into(),
        reference[1].into(),
        reference[2].into(),
    )
    .to_real();
    let offset = &(Ball::pi(128).mid() / q("2"))
        * &catallax_numbers::Rational::from_i64(phase.try_into().unwrap());
    for uv in surface.loops.iter_mut().flatten().flat_map(|t| &mut t.uv) {
        uv[0] = &uv[0] - &offset;
    }
    d
}

#[test]
fn finite_cylinder_exchange_with_toleranced_source_vertex_incidence() {
    let displacements = ["-1e-10", "-1e-12", "-1e-14", "0", "1e-14", "1e-12", "1e-10"];
    for oblique in [false, true] {
        for phase in 0..4 {
            for vertex in [1, 2, 4, 5] {
                for (noise, displacement) in displacements.iter().enumerate() {
                    let mut definition = chart(phase);
                    if oblique {
                        definition = rotated(definition);
                    }
                    definition.vertices[vertex].position[2] =
                        &definition.vertices[vertex].position[2] + &q(displacement);
                    let body = n::build(
                        definition,
                        n::NumericalOptions::default(),
                        &mut Decider::default(),
                    )
                    .unwrap();
                    let file = output::build(
                        &body,
                        &NumericalStepOptions::default(),
                        &mut Decider::default(),
                    )
                    .unwrap();
                    file.check_written(file.bytes(), &mut Decider::default())
                        .unwrap();
                    assert_eq!(file.report().counts, [5, 9, 6]);
                    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
                        let name = format!(
                            "incidence-oblique-{oblique}-phase-{phase}-vertex-{vertex}-noise-{noise}.step"
                        );
                        std::fs::write(std::path::Path::new(&path).join(name), file.bytes())
                            .unwrap();
                    }
                }
            }
        }
    }
}

#[test]
fn source_vertex_displacement_above_tolerance_is_refused() {
    for displacement in ["-0.02", "0.02"] {
        let mut definition = chart(1);
        definition.vertices[1].position[2] = q(displacement);
        assert!(
            n::build(
                definition,
                n::NumericalOptions::default(),
                &mut Decider::default(),
            )
            .is_err()
        );
    }
}

#[test]
fn tiny_tilted_carriers_keep_existing_flat_pcurves() {
    for oblique in [false, true] {
        for phase in 0..4 {
            for (noise, tilt) in ["-1e-14", "-1e-16", "1e-16", "1e-14"].iter().enumerate() {
                let mut definition = chart(phase);
                for edge in &mut definition.edges {
                    if let Curve::Circle(circle) = &mut edge.curve {
                        circle.normal = RVector3::new(q("0"), q(tilt), q("1")).to_real();
                    }
                }
                let face = definition.faces.last_mut().unwrap();
                let Surface::Cylinder(cylinder) = &mut face.surface else {
                    unreachable!()
                };
                let offset = q("37.3578654636477");
                cylinder.axis = RVector3::new(q("0"), q("0"), q("-1")).to_real();
                cylinder.origin.z = Real::from(-&offset);
                for uv in face
                    .loops
                    .iter_mut()
                    .flatten()
                    .flat_map(|trim| &mut trim.uv)
                {
                    uv[0] = -&uv[0];
                    uv[1] = -(&uv[1] + &offset);
                }
                if oblique {
                    definition = rotated(definition);
                }
                let body = n::build(
                    definition,
                    n::NumericalOptions::default(),
                    &mut Decider::default(),
                )
                .unwrap();
                let file = output::build(
                    &body,
                    &NumericalStepOptions::default(),
                    &mut Decider::default(),
                )
                .unwrap();
                file.check_written(file.bytes(), &mut Decider::default())
                    .unwrap();
                assert_eq!(file.report().counts, [5, 9, 6]);
                if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
                    let name = format!("tilted-oblique-{oblique}-phase-{phase}-noise-{noise}.step");
                    std::fs::write(std::path::Path::new(&path).join(name), file.bytes()).unwrap();
                }
            }
        }
    }
}

#[test]
fn circle_tilt_above_source_tolerance_is_refused() {
    for tilt in ["-0.02", "0.02"] {
        let mut definition = chart(1);
        for edge in &mut definition.edges {
            if let Curve::Circle(circle) = &mut edge.curve {
                circle.normal = RVector3::new(q("0"), q(tilt), q("1")).to_real();
            }
        }
        assert!(
            n::build(
                definition,
                n::NumericalOptions::default(),
                &mut Decider::default(),
            )
            .is_err()
        );
    }
}
