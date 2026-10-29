//! Closed-form winding for planar line/circle loops and complete circular boundaries.
//! Other loops use certified traces, with the same orientation convention.

use crate::geometry::{product, refusal, sign};
use crate::region::Location;
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Curve, PointId, Surface};
use catallax_topology::{FaceId, Loop, LoopId, Model, Sense};

mod cylinder_patch;

fn ring(model: &Model, l: LoopId) -> Option<(catallax_topology::EdgeId, Sense)> {
    let Loop::Coedges(cs) = model.loop_(l) else {
        return None;
    };
    if cs.len() != 1 {
        return None;
    }
    let c = model.coedge(cs[0]);
    let e = model.edge(c.edge);
    if e.start.is_some() || e.end.is_some() {
        return None;
    }
    Some((c.edge, product(e.sense, c.sense)))
}

fn ring_side(
    model: &Model,
    l: LoopId,
    p: PointId,
    dec: &mut Decider,
) -> Result<Option<(Sign, Sense)>, Refusal> {
    let Some((edge, travel)) = ring(model, l) else {
        return Ok(None);
    };
    let CurveForm::Circle { center, normal, .. } =
        curve_form(model.geom(), model.edge(edge).curve, dec)?
    else {
        return Ok(None);
    };
    let value = model.geom().position(p).minus(&center).dot(&normal);
    Ok(Some((
        sign(&value, dec, "boolean.classify-ambiguous")?,
        travel,
    )))
}

fn chord_winding(
    (ax, ay): (&Real, &Real),
    (bx, by): (&Real, &Real),
    (qx, qy): (&Real, &Real),
    dec: &mut Decider,
) -> Result<Option<(i64, Sign)>, Refusal> {
    let dx = qx - ax;
    let dy = qy - ay;
    let cross = &(bx - ax) * &dy - &(by - ay) * &dx;
    let cross = sign(&cross, dec, "boolean.classify-ambiguous")?;
    let ex = qx - bx;
    if cross == Sign::Zero {
        let between = &dx * &ex + &dy * &(qy - by);
        if sign(&between, dec, "boolean.classify-ambiguous")? != Sign::Positive {
            return Ok(None);
        }
    }
    let left = sign(&dx, dec, "boolean.classify-ambiguous")?;
    let right = sign(&ex, dec, "boolean.classify-ambiguous")?;
    let value = if left != Sign::Negative && right == Sign::Negative && cross == Sign::Negative {
        -1
    } else if right != Sign::Negative && left == Sign::Negative && cross == Sign::Positive {
        1
    } else {
        0
    };
    Ok(Some((value, cross)))
}

fn native_planar_loops(model: &Model, loops: impl IntoIterator<Item = LoopId>) -> bool {
    loops.into_iter().all(|l| {
        let Loop::Coedges(coedges) = model.loop_(l) else {
            return false;
        };
        coedges.iter().all(|&coedge| {
            matches!(
                model
                    .geom()
                    .curve(model.edge(model.coedge(coedge).edge).curve),
                Curve::Line(_) | Curve::Circle(_)
            )
        })
    })
}

fn plane_winding(
    model: &Model,
    surface: &catallax_shapes::Plane,
    l: LoopId,
    p: PointId,
    dec: &mut Decider,
    native_arcs: bool,
) -> Result<Option<i64>, Refusal> {
    let query = model.geom().position(p);
    if let Some((edge, travel)) = ring(model, l) {
        let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = curve_form(model.geom(), model.edge(edge).curve, dec)?
        else {
            return Ok(None);
        };
        let side = sign(
            &(query.minus(&center).norm_squared() - radius_sq),
            dec,
            "boolean.classify-ambiguous",
        )?;
        if side == Sign::Zero {
            return Ok(None);
        }
        if side == Sign::Positive {
            return Ok(Some(0));
        }
        let agree = sign(
            &surface.normal.dot(&normal),
            dec,
            "boolean.classify-ambiguous",
        )? == Sign::Positive;
        return Ok(Some(if agree == (travel == Sense::Same) {
            1
        } else {
            -1
        }));
    }
    let Loop::Coedges(cs) = model.loop_(l) else {
        return Ok(None);
    };
    let u = surface.ref_dir.clone();
    let v = surface.normal.cross(&surface.ref_dir);
    let origin = surface.origin.clone();
    let d = query.minus(&origin);
    let qx = d.dot(&u);
    let qy = d.dot(&v);
    let project = |vertex| {
        let d = model
            .geom()
            .position(model.vertex(vertex).point)
            .minus(&origin);
        (d.dot(&u), d.dot(&v))
    };
    let mut previous = None;
    let mut result = 0;
    for &coedge in cs {
        let coedge = model.coedge(coedge);
        let edge = model.edge(coedge.edge);
        let form = curve_form(model.geom(), edge.curve, dec)?;
        match &form {
            CurveForm::Line { .. } => {}
            CurveForm::Circle { .. } if native_arcs => {}
            _ => return Ok(None),
        }
        let (Some(start), Some(end)) = (edge.start, edge.end) else {
            return Ok(None);
        };
        let (start, end) = if coedge.sense == Sense::Same {
            (start, end)
        } else {
            (end, start)
        };
        let (ax, ay) = match previous.take() {
            Some((vertex, coordinates)) if vertex == start => coordinates,
            _ => project(start),
        };
        let (bx, by) = project(end);
        let Some((chord, cross)) = chord_winding((&ax, &ay), (&bx, &by), (&qx, &qy), dec)? else {
            return Ok(None);
        };
        result += chord;
        if let CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } = form
        {
            let side = sign(
                &(query.minus(&center).norm_squared() - radius_sq),
                dec,
                "boolean.classify-ambiguous",
            )?;
            if side != Sign::Positive {
                let agree = sign(
                    &surface.normal.dot(&normal),
                    dec,
                    "boolean.classify-ambiguous",
                )?;
                if agree == Sign::Zero {
                    return Ok(None);
                }
                let ccw =
                    (agree == Sign::Positive) == (product(edge.sense, coedge.sense) == Sense::Same);
                let closed = if cross != Sign::Zero {
                    false
                } else if start == end {
                    true
                } else {
                    let dx = &bx - &ax;
                    let dy = &by - &ay;
                    sign(&(&dx * &dx + &dy * &dy), dec, "boolean.classify-ambiguous")? == Sign::Zero
                };
                // Arc minus chord bounds the disk segment on the arc side.
                // Its winding needs only signs, not a new intersection witness.
                if closed || cross == if ccw { Sign::Negative } else { Sign::Positive } {
                    if side == Sign::Zero {
                        return Ok(None);
                    }
                    result += if ccw { 1 } else { -1 };
                }
            }
        }
        previous = Some((end, (bx, by)));
    }
    Ok(Some(result))
}

fn signature(
    model: &Model,
    surface: &Surface,
    loops: &[LoopId],
    point: PointId,
    dec: &mut Decider,
    native_arcs: bool,
) -> Result<Option<i64>, Refusal> {
    match surface {
        Surface::Plane(plane) => {
            let mut value = 0;
            for &l in loops {
                let Some(w) = plane_winding(model, plane, l, point, dec, native_arcs)? else {
                    return Ok(None);
                };
                value += w;
            }
            Ok(Some(value))
        }
        Surface::Sphere(_) | Surface::Cylinder(_) | Surface::Cone(_) => {
            if loops.len() != 1 {
                return Ok(None);
            }
            if matches!(model.loop_(loops[0]), Loop::Vertex(_)) {
                return Ok(Some(0));
            }
            let Some((side, _)) = ring_side(model, loops[0], point, dec)? else {
                return Ok(None);
            };
            match side {
                Sign::Positive => Ok(Some(1)),
                Sign::Negative => Ok(Some(0)),
                Sign::Zero => Err(refusal(
                    "boolean.classify-ambiguous",
                    "cell probe is on a circular boundary",
                )),
            }
        }
        _ => Ok(None),
    }
}

pub(crate) fn signatures(
    model: &Model,
    surface: &Surface,
    loop_groups: &[Vec<LoopId>],
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Vec<i64>>, Refusal> {
    // A later unsupported loop/group would discard the arc decisions. Keep
    // the legacy predicate order instead of proving expensive unused residuals.
    let native_arcs = matches!(surface, Surface::Plane(_))
        && native_planar_loops(model, loop_groups.iter().flatten().copied());
    let mut values = Vec::with_capacity(loop_groups.len());
    for loops in loop_groups {
        let Some(value) = signature(model, surface, loops, point, dec, native_arcs)? else {
            return Ok(None);
        };
        values.push(value);
    }
    Ok(Some(values))
}

pub(crate) fn locate(
    model: &Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let f = model.face(face);
    let surface = model.geom().surface(f.surface);
    if let Surface::Cylinder(cylinder) = surface
        && let Some(location) = cylinder_patch::locate(model, face, cylinder, point, dec)?
    {
        return Ok(Some(location));
    }
    match surface {
        Surface::Plane(plane) => {
            let native_arcs = native_planar_loops(model, f.loops.iter().copied());
            let mut value = 0;
            for &l in &f.loops {
                let Some(w) = plane_winding(model, plane, l, point, dec, native_arcs)? else {
                    return Ok(None);
                };
                value += w;
            }
            if f.sense == Sense::Reversed {
                value = -value;
            }
            match value {
                0 => Ok(Some(Location::Outside)),
                1 => Ok(Some(Location::Inside)),
                _ => Err(refusal(
                    "boolean.classify-ambiguous",
                    "polygon winding is not a unit face region",
                )),
            }
        }
        Surface::Sphere(_) | Surface::Cylinder(_) | Surface::Cone(_) => {
            if let Surface::Cone(cone) = surface
                && sign(
                    &model
                        .geom()
                        .position(point)
                        .minus(&cone.apex)
                        .dot(&cone.axis),
                    dec,
                    "boolean.classify-ambiguous",
                )? == Sign::Negative
            {
                return Ok(Some(Location::Outside));
            }
            let mut boundary = false;
            for &l in &f.loops {
                if let Loop::Vertex(vertex) = model.loop_(l) {
                    if crate::geometry::coincident(
                        model
                            .geom()
                            .same_point(model.vertex(*vertex).point, point, dec),
                        "boolean.classify-ambiguous",
                    )? {
                        return Ok(Some(Location::Boundary));
                    }
                    continue;
                }
                let Some((side, travel)) = ring_side(model, l, point, dec)? else {
                    return Ok(None);
                };
                let front = travel == f.sense;
                match side {
                    Sign::Zero => boundary = true,
                    Sign::Positive if !front => return Ok(Some(Location::Outside)),
                    Sign::Negative if front => return Ok(Some(Location::Outside)),
                    _ => {}
                }
            }
            Ok(Some(if boundary {
                Location::Boundary
            } else {
                Location::Inside
            }))
        }
        _ => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Circle, Curve, GeomKey, Line, Plane, Point3, PointDef, Vector3};
    use catallax_topology::{Coedge, Edge, Face, Vertex};

    fn point(x: i64, y: i64) -> Point3 {
        Point3::new(Real::from_i64(x), Real::from_i64(y), Real::zero())
    }

    fn vector(x: i64, y: i64, z: i64) -> Vector3 {
        Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
    }

    fn arc_face(
        model: &mut Model,
        radius: i64,
        end: [i64; 2],
        through_origin: bool,
        reversed: bool,
        dec: &mut Decider,
    ) -> FaceId {
        let surface = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(point(0, 0), vector(0, 0, 3), vector(2, 0, 0), dec).unwrap(),
        ));
        let mut positions = vec![point(radius, 0)];
        if end != [radius, 0] {
            positions.push(point(end[0], end[1]));
        }
        if through_origin {
            positions.push(point(0, 0));
        }
        let vertices: Vec<_> = positions
            .iter()
            .map(|p| {
                let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
                model.add_vertex(Vertex { point })
            })
            .collect();
        let mut coedges = Vec::new();
        for index in 0..positions.len() {
            let next = (index + 1) % positions.len();
            let curve = if index == 0 {
                Curve::Circle(
                    Circle::new(
                        point(0, 0),
                        vector(0, 0, 1),
                        Real::from_i64(radius * radius),
                        vector(1, 0, 0),
                        dec,
                    )
                    .unwrap(),
                )
            } else {
                Curve::Line(Line {
                    origin: positions[index].clone(),
                    dir: positions[next].minus(&positions[index]),
                })
            };
            let curve = model.geom_mut().add_curve(curve);
            for (a, b) in [
                (
                    GeomKey::Point(model.vertex(vertices[index]).point),
                    GeomKey::Curve(curve),
                ),
                (
                    GeomKey::Point(model.vertex(vertices[next]).point),
                    GeomKey::Curve(curve),
                ),
                (GeomKey::Curve(curve), GeomKey::Surface(surface)),
            ] {
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a,
                    b,
                    at: None,
                    origin: FactOrigin::Construction("boolean.test.circular-arc"),
                });
            }
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
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

    fn locations(model: &mut Model, face: FaceId, expected: &[([i64; 2], Location)]) {
        for &(xy, expected) in expected {
            let p = model
                .geom_mut()
                .add_point(PointDef::Explicit(point(xy[0], xy[1])));
            assert_eq!(
                crate::region::locate(model, face, p, &mut Decider::default()).unwrap(),
                expected,
                "point {xy:?}",
            );
        }
    }

    #[test]
    fn minor_and_major_arcs_keep_the_selected_segment_and_boundary() {
        for reversed in [false, true] {
            let mut model = Model::new();
            let mut dec = Decider::default();
            let minor = arc_face(&mut model, 5, [0, 5], false, reversed, &mut dec);
            locations(
                &mut model,
                minor,
                &[
                    ([3, 3], Location::Inside),
                    ([2, 2], Location::Outside),
                    ([1, 4], Location::Boundary),
                    ([3, 4], Location::Boundary),
                    ([4, -3], Location::Outside),
                ],
            );
            let major = arc_face(&mut model, 5, [0, -5], false, reversed, &mut dec);
            locations(
                &mut model,
                major,
                &[
                    ([0, 0], Location::Inside),
                    ([3, -3], Location::Outside),
                    ([3, 4], Location::Boundary),
                    ([4, -3], Location::Outside),
                ],
            );
        }
    }

    #[test]
    fn an_artificial_arc_chord_is_not_a_face_boundary() {
        let mut model = Model::new();
        let face = arc_face(&mut model, 5, [0, 5], true, false, &mut Decider::default());
        locations(
            &mut model,
            face,
            &[
                ([2, 3], Location::Inside),
                ([1, 1], Location::Inside),
                ([-1, 1], Location::Outside),
                ([0, 3], Location::Boundary),
                ([4, 3], Location::Boundary),
            ],
        );
    }

    #[test]
    fn seamed_full_arcs_keep_holes_and_both_boundary_orientations() {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let outer = arc_face(&mut model, 10, [10, 0], false, false, &mut dec);
        let inner = arc_face(&mut model, 5, [5, 0], false, true, &mut dec);
        let mut loops = model.face(outer).loops.clone();
        loops.extend_from_slice(&model.face(inner).loops);
        let face = model.add_face(Face {
            surface: model.face(outer).surface,
            sense: Sense::Same,
            loops,
        });
        locations(
            &mut model,
            face,
            &[
                ([0, 0], Location::Outside),
                ([7, 0], Location::Inside),
                ([12, 0], Location::Outside),
                ([5, 0], Location::Boundary),
                ([0, -10], Location::Boundary),
            ],
        );
    }

    #[test]
    fn mixed_conic_boundary_uses_its_certificate_within_the_remaining_work() {
        use catallax_intersect::{CurveSurface, intersect_curve_surface};
        use catallax_numbers::{Outcome, Rational, StrictSign, with_operation_budget};
        use catallax_shapes::{Cylinder, Ellipse};

        let r = |n, d| Real::from(Rational::from_i64(n) / Rational::from_i64(d));
        let origin = Point3::new(r(1, 7), r(2, 11), r(3, 13));
        let u = Vector3::new(r(1, 3), r(2, 3), r(2, 3));
        let v = Vector3::new(r(2, 3), r(1, 3), r(-2, 3));
        let normal = u.cross(&v);
        let at = |x, y| {
            origin
                .plus(&u.scaled(&Real::from_i64(x)))
                .plus(&v.scaled(&Real::from_i64(y)))
        };
        let mut model = Model::new();
        let mut dec = Decider::default();
        let surface = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(origin.clone(), normal.clone(), u.clone(), &mut dec).unwrap(),
        ));
        let circle = model.geom_mut().add_curve(Curve::Circle(
            Circle::new(
                origin.clone(),
                normal.clone(),
                Real::from_i64(25),
                u.clone(),
                &mut dec,
            )
            .unwrap(),
        ));
        let a = model.geom_mut().add_point(PointDef::Explicit(at(5, 0)));
        let b = model.geom_mut().add_point(PointDef::Explicit(at(0, 5)));
        let va = model.add_vertex(Vertex { point: a });
        let vb = model.add_vertex(Vertex { point: b });
        let line = model.geom_mut().add_curve(Curve::Line(Line {
            origin: at(0, 5),
            dir: at(5, 0).minus(&at(0, 5)),
        }));
        let coedges = [(circle, va, vb), (line, vb, va)]
            .into_iter()
            .map(|(curve, start, end)| {
                let edge = model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: Some(start),
                    end: Some(end),
                });
                model.add_coedge(Coedge {
                    edge,
                    sense: Sense::Same,
                })
            })
            .collect();
        let outer = model.add_loop(Loop::Coedges(coedges));
        let ellipse = model.geom_mut().add_curve(Curve::Ellipse(
            Ellipse::new(
                at(3, 3),
                normal.clone(),
                u.clone(),
                r(1, 16),
                r(1, 64),
                &mut dec,
            )
            .unwrap(),
        ));
        let hole_edge = model.add_edge(Edge {
            curve: ellipse,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        let hole_coedge = model.add_coedge(Coedge {
            edge: hole_edge,
            sense: Sense::Reversed,
        });
        let hole = model.add_loop(Loop::Coedges(vec![hole_coedge]));
        for (a, b) in [
            (GeomKey::Point(a), GeomKey::Curve(circle)),
            (GeomKey::Point(b), GeomKey::Curve(circle)),
            (GeomKey::Point(a), GeomKey::Curve(line)),
            (GeomKey::Point(b), GeomKey::Curve(line)),
            (GeomKey::Curve(circle), GeomKey::Surface(surface)),
            (GeomKey::Curve(line), GeomKey::Surface(surface)),
            (GeomKey::Curve(ellipse), GeomKey::Surface(surface)),
        ] {
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a,
                b,
                at: None,
                origin: FactOrigin::Construction("boolean.test.mixed-conic"),
            });
        }
        let face = model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![outer, hole],
        });
        // In this orthonormal frame the section has x²+y²=25 and xy=1.
        // Its algebraic coordinates need not be re-eliminated to prove the
        // already-certified circle incidence of a boundary query.
        let cutter = model.geom_mut().add_surface(Surface::Cylinder(
            Cylinder::from_radius_sq(
                origin.clone(),
                u.plus(&v).plus(&normal),
                Real::from_i64(16),
                u.minus(&v),
                &mut dec,
            )
            .unwrap(),
        ));
        let hits = intersect_curve_surface(model.geom_mut(), circle, cutter, &mut dec).unwrap();
        assert!(hits.marks.is_empty());
        let CurveSurface::Meets(hits) = hits.value else {
            panic!("the circle and elliptic section must meet");
        };
        let query = hits
            .into_iter()
            .find(|hit| {
                let delta = model.geom().position(hit.point).minus(&origin);
                matches!(
                    dec.sign(&delta.dot(&u)),
                    Outcome::Proven(StrictSign::Positive)
                ) && matches!(
                    dec.sign(&delta.dot(&v)),
                    Outcome::Proven(StrictSign::Positive)
                )
            })
            .unwrap()
            .point;
        let (answer, work) = with_operation_budget(10_000, || {
            crate::region::locate(&mut model, face, query, &mut Decider::default())
        });
        assert_eq!(answer.unwrap(), Location::Boundary);
        assert!(!work.exhausted);
    }
}
