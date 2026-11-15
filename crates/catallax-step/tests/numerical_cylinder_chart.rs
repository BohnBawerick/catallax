//! Equivalent full-cylinder charts must preserve the same finite material.
#[path = "../../catallax-topology/tests/common/periodic.rs"]
mod fixture;

use catallax_numbers::{Ball, Decider, Rational};
use catallax_shapes::{RVector3, Surface};
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::{Sense, numerical as n};

fn flip(sense: Sense) -> Sense {
    match sense {
        Sense::Same => Sense::Reversed,
        Sense::Reversed => Sense::Same,
    }
}

fn outer_cylinder() -> n::NumericalDefinition {
    let source = fixture::bore();
    let edge_indices = ["rim-0", "rim-1", "cylinder-u-seam"].map(|name| {
        source
            .edges
            .iter()
            .position(|edge| edge.key == name)
            .unwrap()
    });
    let vertices = [
        source.edges[edge_indices[0]].vertices[0],
        source.edges[edge_indices[1]].vertices[0],
    ];
    let edges = edge_indices.map(|index| {
        let mut edge = source.edges[index].clone();
        edge.vertices = edge
            .vertices
            .map(|v| vertices.iter().position(|i| *i == v).unwrap());
        edge
    });
    let mut faces = Vec::new();
    for index in [0, 1] {
        let mut cap = source.faces[index].clone();
        let mut boundary = cap.loops[1].clone();
        for trim in &mut boundary {
            trim.edge = edge_indices.iter().position(|e| *e == trim.edge).unwrap();
            trim.sense = flip(trim.sense);
            trim.uv.swap(0, 1);
        }
        cap.loops = vec![boundary];
        faces.push(cap);
    }
    let mut wall = source.faces.last().unwrap().clone();
    wall.sense = flip(wall.sense);
    for boundary in &mut wall.loops {
        boundary.reverse();
        for trim in boundary {
            trim.edge = edge_indices.iter().position(|e| *e == trim.edge).unwrap();
            trim.sense = flip(trim.sense);
            trim.uv.swap(0, 1);
        }
    }
    faces.push(wall);
    n::NumericalDefinition {
        name: "authored-periodic-outer-cylinder".into(),
        vertices: vertices.map(|i| source.vertices[i].clone()).to_vec(),
        edges: edges.to_vec(),
        faces,
    }
}

#[test]
fn full_cylinder_exchange_keeps_equivalent_surface_charts() {
    for outer in [false, true] {
        for axis in [-1_i64, 1] {
            for quarter_turns in -3_i64..=3 {
                for reverse_seam in [false, true] {
                    for start in 0..4 {
                        let mut definition = if outer {
                            outer_cylinder()
                        } else {
                            fixture::bore()
                        };
                        let phase = &(Ball::pi(128).mid() / fixture::q("2"))
                            * &Rational::from_i64(quarter_turns);
                        for face in &mut definition.faces {
                            if let Surface::Cylinder(cylinder) = &mut face.surface {
                                cylinder.axis =
                                    RVector3::new(0.into(), 0.into(), axis.into()).to_real();
                                let reference = match quarter_turns.rem_euclid(4) {
                                    0 => [1, 0, 0],
                                    1 => [0, axis, 0],
                                    2 => [-1, 0, 0],
                                    _ => [0, -axis, 0],
                                };
                                cylinder.ref_dir = RVector3::new(
                                    reference[0].into(),
                                    reference[1].into(),
                                    reference[2].into(),
                                )
                                .to_real();
                                for uv in face
                                    .loops
                                    .iter_mut()
                                    .flatten()
                                    .flat_map(|trim| &mut trim.uv)
                                {
                                    uv[0] = &(&uv[0] * &Rational::from_i64(axis)) - &phase;
                                    uv[1] = &uv[1] * &Rational::from_i64(axis);
                                }
                            }
                        }
                        if reverse_seam {
                            let edge = definition
                                .edges
                                .iter_mut()
                                .find(|e| e.key == "cylinder-u-seam")
                                .unwrap();
                            let catallax_shapes::Curve::Line(line) = &mut edge.curve else {
                                unreachable!()
                            };
                            line.dir = RVector3::new(0.into(), 0.into(), (-1).into()).to_real();
                            edge.parameters[1] = fixture::q("-1");
                        }
                        definition.faces.last_mut().unwrap().loops[0].rotate_left(start);
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
                        assert_eq!(
                            file.report().counts,
                            if outer { [3, 3, 2] } else { [7, 15, 10] }
                        );
                        if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
                            let name = format!(
                                "cylinder-chart-outer-{outer}-axis-{axis}-phase-{quarter_turns}-reverse-{reverse_seam}-start-{start}.step"
                            );
                            std::fs::write(std::path::Path::new(&path).join(name), file.bytes())
                                .unwrap();
                        }
                    }
                }
            }
        }
    }
}
