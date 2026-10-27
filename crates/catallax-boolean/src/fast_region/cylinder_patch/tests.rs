use super::*;
use catallax_numbers::{Fact, FactOrigin, Real, with_operation_budget};
use catallax_shapes::{Circle, IntersectionCurve, Line, PointDef, Surface};
use catallax_topology::{Face, Vertex};

fn point([x, y, z]: [i64; 3], placed: bool) -> Point3 {
    let [x, y, z] = if placed {
        [z + 11, x - 7, y + 13]
    } else {
        [x, y, z]
    };
    Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}

fn patch(
    model: &mut Model,
    major: bool,
    normals: [i64; 2],
    reversed: bool,
    placed: bool,
    mismatched_caps: bool,
    incidence: bool,
) -> FaceId {
    let mut dec = Decider::default();
    let origin = point([0, 0, 0], placed);
    let axis = point([0, 0, 3], placed).minus(&origin);
    let ref_dir = point([2, 0, 0], placed).minus(&origin);
    let surface = model.geom_mut().add_surface(Surface::Cylinder(Cylinder {
        origin,
        axis,
        radius_sq: Real::from_i64(25),
        ref_dir: ref_dir.clone(),
    }));
    let y = if major { -5 } else { 5 };
    let positions = [[5, 0, 0], [0, y, 0], [0, y, 4], [5, 0, 4]].map(|p| point(p, placed));
    let vertices = positions.each_ref().map(|position| {
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(position.clone()));
        model.add_vertex(Vertex { point })
    });
    let mut coedges = Vec::new();
    for index in 0..4 {
        let next = (index + 1) % 4;
        let (curve, sense) = if index % 2 == 0 {
            let z = if index == 0 { 0 } else { 4 };
            let center = point([0, 0, z], placed);
            let normal = point([0, 0, z + 3 * normals[index / 2]], placed).minus(&center);
            let curve = Curve::Circle(
                Circle::new(
                    center,
                    normal,
                    Real::from_i64(25),
                    ref_dir.clone(),
                    &mut dec,
                )
                .unwrap(),
            );
            let forward = (index == 0) == (normals[index / 2] > 0);
            let forward = if index == 2 && mismatched_caps {
                !forward
            } else {
                forward
            };
            (
                curve,
                if forward {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
            )
        } else {
            (
                Curve::Line(Line {
                    origin: positions[index].clone(),
                    dir: positions[next].minus(&positions[index]),
                }),
                Sense::Same,
            )
        };
        let curve = model.geom_mut().add_curve(curve);
        if incidence {
            for vertex in [vertices[index], vertices[next]] {
                let point = model.vertex(vertex).point;
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Point(point),
                    b: GeomKey::Curve(curve),
                    at: None,
                    origin: FactOrigin::Construction("boolean.test.cylinder-patch"),
                });
            }
        }
        let edge = model.add_edge(Edge {
            curve,
            sense,
            start: Some(vertices[index]),
            end: Some(vertices[next]),
        });
        coedges.push(model.add_coedge(Coedge {
            edge,
            sense: if reversed {
                Sense::Reversed
            } else {
                Sense::Same
            },
        }));
    }
    if reversed {
        coedges.reverse();
    }
    let loop_ = model.add_loop(Loop::Coedges(coedges));
    model.add_face(Face {
        surface,
        sense: if reversed {
            Sense::Reversed
        } else {
            Sense::Same
        },
        loops: vec![loop_],
    })
}

fn intersection_generators(model: &mut Model, face: FaceId, placed: bool) -> FaceId {
    let mut face = model.face(face).clone();
    let Surface::Cylinder(cylinder) = model.geom().surface(face.surface).clone() else {
        unreachable!()
    };
    let Loop::Coedges(original) = model.loop_(face.loops[0]).clone() else {
        unreachable!()
    };
    let offset = point([3, 4, 0], placed).minus(&point([0, 0, 0], placed));
    let mut coedges = Vec::new();
    for id in original {
        let mut coedge = model.coedge(id).clone();
        let mut edge = model.edge(coedge.edge).clone();
        if let Curve::Line(line) = model.geom().curve(edge.curve).clone() {
            let start = model.vertex(edge.start.unwrap()).point;
            let end = model.vertex(edge.end.unwrap()).point;
            let other = model.geom_mut().add_surface(Surface::Cylinder(Cylinder {
                origin: line.origin.plus(&offset.scaled(&Real::from_i64(-1))),
                ..cylinder.clone()
            }));
            edge.curve = model
                .geom_mut()
                .add_curve(Curve::Intersection(IntersectionCurve {
                    surfaces: [face.surface, other],
                    witness: start,
                }));
            let mut dec = Decider::default();
            let [natural, _] =
                catallax_intersect::tangent_at(model.geom(), edge.curve, start, &mut dec).unwrap();
            edge.sense = if sign(&natural.dot(&line.dir), &mut dec, CODE).unwrap() == Sign::Positive
            {
                Sense::Same
            } else {
                Sense::Reversed
            };
            for point in [start, end] {
                crate::geometry::record_incident(
                    model,
                    GeomKey::Point(point),
                    GeomKey::Curve(edge.curve),
                );
            }
            coedge.edge = model.add_edge(edge);
        }
        coedges.push(model.add_coedge(coedge));
    }
    face.loops = vec![model.add_loop(Loop::Coedges(coedges))];
    model.add_face(face)
}

#[test]
fn product_trim_preserves_major_minor_arcs_caps_and_generators() {
    for major in [false, true] {
        for normals in [[1, 1], [1, -1], [-1, 1], [-1, -1]] {
            for reversed in [false, true] {
                for (placed, intersections) in
                    [(false, false), (true, false), (false, true), (true, true)]
                {
                    let mut model = Model::new();
                    let face = patch(&mut model, major, normals, reversed, placed, false, true);
                    let face = if intersections {
                        intersection_generators(&mut model, face, placed)
                    } else {
                        face
                    };
                    let between = if major {
                        Location::Inside
                    } else {
                        Location::Outside
                    };
                    let end_y = if major { -5 } else { 5 };
                    for (position, expected) in [
                        ([3, 4, 2], Location::Inside),
                        ([-3, 4, 2], between),
                        ([-3, -4, 2], between),
                        ([3, -4, 2], Location::Outside),
                        ([5, 0, 2], Location::Boundary),
                        ([0, end_y, 2], Location::Boundary),
                        ([3, 4, 0], Location::Boundary),
                        ([3, 4, 4], Location::Boundary),
                        ([3, 4, -1], Location::Outside),
                        ([3, 4, 5], Location::Outside),
                        ([5, 0, -1], Location::Outside),
                        ([3, -4, 0], Location::Outside),
                    ] {
                        let query = model
                            .geom_mut()
                            .add_point(PointDef::Explicit(point(position, placed)));
                        let Surface::Cylinder(cylinder) =
                            model.geom().surface(model.face(face).surface)
                        else {
                            unreachable!()
                        };
                        assert_eq!(
                            locate(&model, face, cylinder, query, &mut Decider::default()).unwrap(),
                            Some(expected),
                            "{position:?}, major={major}, normals={normals:?}, reversed={reversed}, placed={placed}, intersections={intersections}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn inconsistent_cap_travel_and_face_orientation_keep_the_fallback() {
    for bad_caps in [false, true] {
        let mut model = Model::new();
        let face = patch(&mut model, false, [1, 1], false, false, bad_caps, true);
        let face = if bad_caps {
            face
        } else {
            let mut wrong = model.face(face).clone();
            wrong.sense = Sense::Reversed;
            model.add_face(wrong)
        };
        let query = model
            .geom_mut()
            .add_point(PointDef::Explicit(point([3, 4, 2], false)));
        let Surface::Cylinder(cylinder) = model.geom().surface(model.face(face).surface) else {
            unreachable!()
        };
        assert_eq!(
            locate(&model, face, cylinder, query, &mut Decider::default()).unwrap(),
            None
        );
    }
}

#[test]
fn missing_incidence_uses_the_ordinary_exact_membership_query() {
    let mut model = Model::new();
    let face = patch(&mut model, false, [1, 1], false, false, false, false);
    let query = model
        .geom_mut()
        .add_point(PointDef::Explicit(point([3, 4, 2], false)));
    let Surface::Cylinder(cylinder) = model.geom().surface(model.face(face).surface) else {
        unreachable!()
    };
    assert_eq!(
        locate(&model, face, cylinder, query, &mut Decider::default()).unwrap(),
        None
    );
    assert_eq!(
        crate::region::locate(&mut model, face, query, &mut Decider::default()).unwrap(),
        Location::Inside
    );
}

#[test]
fn native_patch_obeys_the_callers_remaining_work() {
    let mut model = Model::new();
    let face = patch(&mut model, false, [1, 1], false, false, false, true);
    let query = model
        .geom_mut()
        .add_point(PointDef::Explicit(point([3, 4, 2], false)));
    let (answer, work) = with_operation_budget(0, || {
        crate::region::locate(&mut model, face, query, &mut Decider::default())
    });
    assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
}
