//! Complete planar paths in arrangements containing retained finite sections.

use crate::curve_pair::{Meetings, meetings};
use crate::geometry::{natural_ends, on_arc};
use catallax_intersect::{is_bezier_extrusion_section, on_curve, same_curve, tangent_at};
use catallax_numbers::{Decider, Outcome, Real, Refusal, Sign, StrictSign, ZeroProof};
use catallax_shapes::{Curve, CurveId, Line, Point3, PointId, Surface, SurfaceId, Vector3};
use catallax_topology::{EdgeId, Model};

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(p) if p.proof != ZeroProof::SeparationBound => Ok(Sign::Zero),
        Outcome::ProvenCoincident(_) => Err(Refusal::new(
            "numbers.geometric-gap-unproven",
            "a finite planar probe needs exact incidence",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

pub(super) fn adjacent(
    model: &Model,
    surface: SurfaceId,
    source: (PointId, &Point3, CurveId),
    q: &Point3,
    left: &Vector3,
    barriers: &mut [EdgeId],
    dec: &mut Decider,
) -> Result<Option<bool>, Refusal> {
    let Surface::Plane(plane) = model.geom().surface(surface) else {
        return Ok(None);
    };
    let mut polynomial = matches!(
        model.geom().curve(source.2),
        Curve::PatchSection(_) | Curve::BSpline(_)
    ) || is_bezier_extrusion_section(model.geom(), source.2, dec)?;
    if !polynomial {
        for edge in barriers.iter() {
            if matches!(
                model.geom().curve(model.edge(*edge).curve),
                Curve::PatchSection(_) | Curve::BSpline(_)
            ) || is_bezier_extrusion_section(model.geom(), model.edge(*edge).curve, dec)?
            {
                polynomial = true;
                break;
            }
        }
    }
    if !polynomial {
        return Ok(None);
    }
    let before_separation = dec.stats().coincident_separation;
    let (source_id, p, source_curve) = source;
    if sign(
        &model.geom().position(source_id).minus(p).norm_squared(),
        dec,
    )? != Sign::Zero
    {
        return Err(Refusal::new(
            "boolean.bezier-probe-source",
            "the source point differs from its retained arc probe",
        ));
    }
    match on_curve(model.geom(), source_curve, source_id, dec) {
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {}
        Outcome::Refused(error) => return Err(error),
        Outcome::NotExact(error) => return Err(Refusal::from(error)),
        _ => {
            return Err(Refusal::new(
                "boolean.bezier-probe-source",
                "the source point is not proved on its curve",
            ));
        }
    }
    let tangent = tangent_at(model.geom(), source_curve, source_id, dec)?[1].clone();
    if sign(&left.norm_squared(), dec)? != Sign::Positive
        || sign(&left.dot(&tangent), dec)? != Sign::Zero
        || sign(&left.dot(&plane.normal), dec)? != Sign::Zero
    {
        return Ok(Some(false));
    }
    if sign(&p.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
        || sign(&q.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
    {
        return Ok(Some(false));
    }
    let delta = q.minus(p);
    let squared = delta.norm_squared();
    if sign(&squared, dec)? != Sign::Positive || sign(&left.dot(&delta), dec)? != Sign::Positive {
        return Ok(Some(false));
    }
    // A positive left derivative chooses the source arc's local left cell.
    // Complete finite barrier meetings then certify the whole remaining path.
    let mut fork = model.clone();
    let path = fork
        .geom_mut()
        .add_curve(Curve::Line(Line::new(p.clone(), delta.clone(), dec)?));
    let mut blocked_index = None;
    for (index, &edge_id) in barriers.iter().enumerate() {
        let edge = fork.edge(edge_id).clone();
        let (from, to) = natural_ends(&fork, edge_id);
        let station = |point: &Point3| point.minus(p).dot(&delta) / squared.clone();
        let blocked = match meetings(&mut fork, path, edge.curve, dec)? {
            Meetings::Coincident => {
                let (Some(from), Some(to)) = (from, to) else {
                    return Ok(Some(false));
                };
                let from = station(&fork.geom().position(from));
                let to = station(&fork.geom().position(to));
                let before =
                    sign(&from, dec)? == Sign::Negative && sign(&to, dec)? == Sign::Negative;
                let after = sign(&(&from - &Real::one()), dec)? == Sign::Positive
                    && sign(&(&to - &Real::one()), dec)? == Sign::Positive;
                !before && !after
            }
            Meetings::Points(points) => {
                let mut blocked = false;
                for point in points {
                    if !on_arc(
                        &fork,
                        edge.curve,
                        from,
                        to,
                        point,
                        dec,
                        "boolean.imprint-ambiguous",
                    )? {
                        continue;
                    }
                    let parameter = station(&fork.geom().position(point));
                    match sign(&parameter, dec)? {
                        Sign::Negative => {}
                        Sign::Positive => {
                            if sign(&(&parameter - &Real::one()), dec)? != Sign::Positive {
                                blocked = true;
                                break;
                            }
                        }
                        Sign::Zero => {
                            if edge.curve != source_curve
                                && !matches!(same_curve(fork.geom(),edge.curve,source_curve,dec),Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound)
                            {
                                blocked = true;
                                break;
                            }
                        }
                    }
                }
                blocked
            }
        };
        if blocked {
            blocked_index = Some(index);
            break;
        }
    }
    if let Some(index) = blocked_index {
        barriers.swap(0, index);
        return Ok(Some(false));
    }
    if dec.stats().coincident_separation != before_separation {
        return Err(Refusal::new(
            "numbers.geometric-gap-unproven",
            "a finite planar probe used an unresolved separation allowance",
        ));
    }
    Ok(Some(true))
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_intersect::{Component, SurfacePair, intersect_surfaces};
    use catallax_numbers::Rational;
    use catallax_shapes::{BezierExtrusion, Plane, PointDef, RPoint3, RVector3};
    use catallax_topology::{Edge, Sense, Vertex};

    fn p(x: i64, y: i64, z: i64) -> Point3 {
        RPoint3::new(x.into(), y.into(), z.into()).to_real()
    }
    fn v(x: i64, y: i64, z: i64) -> Vector3 {
        RVector3::new(x.into(), y.into(), z.into()).to_real()
    }
    fn edge(model: &mut Model, curve: CurveId, a: Point3, b: Point3) -> EdgeId {
        let vertices = [a, b].map(|point| {
            let point = model.geom_mut().add_point(PointDef::Explicit(point));
            model.add_vertex(Vertex { point })
        });
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(vertices[0]),
            end: Some(vertices[1]),
        })
    }
    fn line(model: &mut Model, a: Point3, b: Point3) -> EdgeId {
        let curve = model.geom_mut().add_curve(Curve::Line(
            Line::new(a.clone(), b.minus(&a), &mut Decider::default()).unwrap(),
        ));
        edge(model, curve, a, b)
    }
    fn fixture() -> (Model, SurfaceId, PointId, CurveId, Vec<EdgeId>) {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let wall = BezierExtrusion::new(
            vec![(0, 0), (0, 4), (4, 4), (4, 0)]
                .into_iter()
                .map(|(x, y)| RPoint3::new(x.into(), y.into(), 0.into()))
                .collect(),
            RVector3::new(0.into(), 0.into(), 2.into()),
            &mut dec,
        )
        .unwrap();
        let wall = model.geom_mut().add_surface(Surface::BSpline(wall.into()));
        let plane = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(p(0, 0, 1), v(0, 0, 1), v(1, 0, 0), &mut dec).unwrap(),
        ));
        let SurfacePair::Meets(parts) = intersect_surfaces(model.geom_mut(), plane, wall, &mut dec)
            .unwrap()
            .value
        else {
            panic!("section expected")
        };
        let Component::Curve { curve, .. } = parts[0] else {
            panic!("curve expected")
        };
        let curved = edge(&mut model, curve, p(0, 0, 1), p(4, 0, 1));
        let base = line(&mut model, p(0, 0, 1), p(4, 0, 1));
        let source = model.geom_mut().add_point(PointDef::Explicit(p(2, 3, 1)));
        (model, plane, source, curve, vec![curved, base])
    }
    #[test]
    fn complete_paths_stop_at_the_baseline_and_added_arrangement_barriers() {
        let (mut model, plane, source, curve, mut barriers) = fixture();
        let position = model.geom().position(source);
        let left = v(0, -1, 0);
        for (target, clear) in [
            (p(2, 1, 1), true),
            (p(2, 0, 1), false),
            (p(2, -1, 1), false),
            (p(2, 4, 1), false),
        ] {
            assert_eq!(
                adjacent(
                    &model,
                    plane,
                    (source, &position, curve),
                    &target,
                    &left,
                    &mut barriers,
                    &mut Decider::default()
                )
                .unwrap(),
                Some(clear)
            );
        }
        barriers.push(line(&mut model, p(1, 2, 1), p(3, 2, 1)));
        assert_eq!(
            adjacent(
                &model,
                plane,
                (source, &position, curve),
                &p(2, 1, 1),
                &left,
                &mut barriers,
                &mut Decider::default()
            )
            .unwrap(),
            Some(false)
        );
        let above = Point3::new(
            2.into(),
            Rational::from_ratio(5, 2).unwrap().into(),
            1.into(),
        );
        assert_eq!(
            adjacent(
                &model,
                plane,
                (source, &position, curve),
                &above,
                &left,
                &mut barriers,
                &mut Decider::default()
            )
            .unwrap(),
            Some(true)
        );
    }
    #[test]
    fn collinear_barriers_are_clipped_to_their_physical_segments() {
        let (mut model, plane, source, curve, mut barriers) = fixture();
        let position = model.geom().position(source);
        barriers.push(line(&mut model, p(2, 10, 1), p(2, 11, 1)));
        assert_eq!(
            adjacent(
                &model,
                plane,
                (source, &position, curve),
                &p(2, 1, 1),
                &v(0, -1, 0),
                &mut barriers,
                &mut Decider::default()
            )
            .unwrap(),
            Some(true)
        );
        barriers.push(line(&mut model, p(2, 1, 1), p(2, 2, 1)));
        assert_eq!(
            adjacent(
                &model,
                plane,
                (source, &position, curve),
                &p(2, 1, 1),
                &v(0, -1, 0),
                &mut barriers,
                &mut Decider::default()
            )
            .unwrap(),
            Some(false)
        );
    }
    #[test]
    fn a_curve_meeting_outside_its_physical_trim_does_not_block_the_path() {
        let (mut model, plane, _, curve, mut barriers) = fixture();
        let end = Point3::new(
            Rational::from_ratio(5, 8).unwrap().into(),
            Rational::from_ratio(9, 4).unwrap().into(),
            1.into(),
        );
        barriers[0] = edge(&mut model, curve, p(0, 0, 1), end);
        let source_curve = model.edge(barriers[1]).curve;
        let position = p(2, 0, 1);
        let source = model
            .geom_mut()
            .add_point(PointDef::Explicit(position.clone()));
        assert_eq!(
            adjacent(
                &model,
                plane,
                (source, &position, source_curve),
                &p(2, 4, 1),
                &v(0, 1, 0),
                &mut barriers,
                &mut Decider::default()
            )
            .unwrap(),
            Some(true)
        );
    }

    #[test]
    fn native_planar_boundary_certifies_the_path_without_an_implicit_curve_equation() {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let plane = model.geom_mut().add_surface(Surface::Plane(
            Plane::new(p(0, 0, 0), v(0, 0, 1), v(1, 0, 0), &mut dec).unwrap(),
        ));
        let curve = model.geom_mut().add_curve(Curve::BSpline(
            catallax_shapes::BSplineCurve::new(
                2,
                vec![0, 0, 0, 1, 1, 1]
                    .into_iter()
                    .map(Rational::from_i64)
                    .collect(),
                vec![p(0, 0, 0), p(2, 4, 0), p(4, 0, 0)],
                vec![Real::one(); 3],
                &mut dec,
            )
            .unwrap(),
        ));
        let roof = edge(&mut model, curve, p(0, 0, 0), p(4, 0, 0));
        let floor = line(&mut model, p(0, 0, 0), p(4, 0, 0));
        let source = model.geom_mut().add_point(PointDef::Explicit(p(2, 2, 0)));
        let position = model.geom().position(source);
        let mut barriers = vec![roof, floor];
        for (target, clear) in [
            (p(2, 1, 0), true),
            (p(2, 0, 0), false),
            (p(2, -1, 0), false),
            (p(2, 3, 0), false),
        ] {
            assert_eq!(
                adjacent(
                    &model,
                    plane,
                    (source, &position, curve),
                    &target,
                    &v(0, -1, 0),
                    &mut barriers,
                    &mut dec
                )
                .unwrap(),
                Some(clear)
            );
        }
    }
}
