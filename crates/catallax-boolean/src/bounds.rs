//! Conservative finite face bounds used only to prove irrelevant support pairs.
//! A strictly nonzero support enclosure excludes an intersection; a box overlap
//! never claims a meeting or changes the analytic intersection method.

use crate::{
    geometry::sign,
    region::{Location, locate},
};
use catallax_numbers::{Ball, Decider, FactKind, Rational, Real, Refusal, Sign};
use catallax_shapes::{Circle, Curve, GeomKey, Plane, Point3, PointDef, Surface};
use catallax_topology::{FaceId, Loop, Model};
use std::collections::BTreeMap;

fn enclose(x: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    if let Some(q) = x.rational_ref() {
        Ok(Ball::from_rational(q, 128))
    } else {
        dec.enclose_budgeted_at(x, 128)
    }
}

fn support_box(
    surface: &Surface,
    bounds: &[Ball; 3],
    dec: &mut Decider,
) -> Result<Option<Ball>, Refusal> {
    let mut local = surface.clone();
    let origin = match &mut local {
        Surface::Plane(s) => &mut s.origin,
        Surface::Cylinder(s) => &mut s.origin,
        Surface::Cone(s) => &mut s.apex,
        Surface::Sphere(s) => &mut s.center,
        Surface::Torus(s) => &mut s.center,
        Surface::EllipticCylinder(s) => &mut s.origin,
        _ => return Ok(None),
    };
    let offset = [
        enclose(&origin.x, dec)?,
        enclose(&origin.y, dec)?,
        enclose(&origin.z, dec)?,
    ];
    *origin = Point3::origin();
    // Translate before expanding the support polynomial. Evaluating its global
    // expansion loses local gaps through interval cancellation at large offsets.
    let bounds: [Ball; 3] = std::array::from_fn(|i| &bounds[i] - &offset[i]);
    local
        .implicit()
        .map(|polynomial| polynomial.eval_ball(&bounds, dec))
        .transpose()
}

fn sphere_box(center: &Point3, radius: &Real, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    let radius = enclose(radius, dec)?;
    let mut result = Vec::with_capacity(3);
    for x in [&center.x, &center.y, &center.z] {
        let c = enclose(x, dec)?;
        result.push((&c - &radius).hull(&(&c + &radius)));
    }
    Ok(result.try_into().expect("three coordinates"))
}
fn zero(x: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(x, dec, "boolean.bound-ambiguous")? == Sign::Zero)
}

fn circle_box(circle: &Circle, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    let normal_squared = circle.normal.norm_squared();
    let coordinates = [&circle.center.x, &circle.center.y, &circle.center.z];
    let normal = [&circle.normal.x, &circle.normal.y, &circle.normal.z];
    let mut bounds = Vec::with_capacity(3);
    for axis in 0..3 {
        let factor = &Real::one() - &(&normal[axis].square() / &normal_squared);
        let width = if factor.rational_ref().is_some_and(Rational::is_zero) {
            Real::zero()
        } else {
            (&circle.radius_sq * &factor).sqrt()
        };
        let width = enclose(&width, dec)?;
        let center = enclose(coordinates[axis], dec)?;
        bounds.push((&center - &width).hull(&(&center + &width)));
    }
    Ok(bounds.try_into().expect("three coordinates"))
}

fn planar_box(
    model: &mut Model,
    face: FaceId,
    plane: &Plane,
    dec: &mut Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    let mut bounds: Option<[Ball; 3]> = None;
    for lp in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return Ok(None);
        };
        for coedge in coedges {
            let edge = model.edge(model.coedge(*coedge).edge);
            let next = match model.geom().curve(edge.curve) {
                Curve::Circle(circle) => circle_box(circle, dec)?,
                Curve::Line(_) => {
                    let (Some(start), Some(end)) = (edge.start, edge.end) else {
                        return Ok(None);
                    };
                    let a = model.geom().position(model.vertex(start).point);
                    let b = model.geom().position(model.vertex(end).point);
                    let mut result = Vec::with_capacity(3);
                    for (a, b) in [(&a.x, &b.x), (&a.y, &b.y), (&a.z, &b.z)] {
                        result.push(enclose(a, dec)?.hull(&enclose(b, dec)?));
                    }
                    result.try_into().expect("three coordinates")
                }
                _ => return Ok(None),
            };
            if let Some(bounds) = &mut bounds {
                for i in 0..3 {
                    bounds[i] = bounds[i].hull(&next[i]);
                }
            } else {
                bounds = Some(next);
            }
        }
    }
    let Some(bounds) = bounds else {
        return Ok(None);
    };
    // All boundaries lie in this box. Prove the unbounded planar component is
    // outside using a point farther from the origin than the entire box. Every
    // exterior point is in that component, so every retained cell is bounded.
    let mut radius = Rational::zero();
    for (axis, coordinate) in [&plane.origin.x, &plane.origin.y, &plane.origin.z]
        .into_iter()
        .enumerate()
    {
        let origin = enclose(coordinate, dec)?;
        radius = radius
            + (&bounds[axis].lo() - &origin.hi())
                .abs()
                .max((&bounds[axis].hi() - &origin.lo()).abs());
    }
    let norm = plane.ref_dir.norm_squared();
    let length = if norm.rational_ref() == Some(&Rational::one()) {
        Real::one()
    } else {
        norm.sqrt()
    };
    let sample = plane.origin.plus(
        &plane
            .ref_dir
            .scaled(&(&Real::from(radius + Rational::one()) / &length)),
    );
    let point = model.geom_mut().add_point(PointDef::Explicit(sample));
    if locate(model, face, point, dec)? != Location::Outside {
        return Ok(None);
    }
    Ok(Some(bounds))
}

fn face_box(
    model: &mut Model,
    face: FaceId,
    planar: bool,
    dec: &mut Decider,
) -> Result<Option<[Ball; 3]>, Refusal> {
    let surface = model.geom().surface(model.face(face).surface).clone();
    let (origin, axis, reference) = match &surface {
        Surface::Plane(plane) if planar => return planar_box(model, face, plane, dec),
        Surface::Sphere(s) => return sphere_box(&s.center, &s.radius, dec).map(Some),
        Surface::Torus(s) => return sphere_box(&s.center, &(&s.major + &s.minor), dec).map(Some),
        Surface::Cylinder(s) => (&s.origin, &s.axis, &s.ref_dir),
        Surface::Cone(s) => (&s.apex, &s.axis, &s.ref_dir),
        _ => return Ok(None),
    };
    let mut ends = Vec::with_capacity(2);
    for l in &model.face(face).loops {
        match model.loop_(*l) {
            Loop::Vertex(v) => {
                let PointDef::Explicit(p) = model.geom().point(model.vertex(*v).point) else {
                    return Ok(None);
                };
                let same_apex = !surface.has_rational_definition()
                    && p.x.same_definition(&origin.x)
                    && p.y.same_definition(&origin.y)
                    && p.z.same_definition(&origin.z);
                if !matches!(&surface, Surface::Cone(_))
                    || (!same_apex && !zero(&p.minus(origin).norm_squared(), dec)?)
                {
                    return Ok(None);
                }
                ends.push((p.clone(), Real::zero()));
            }
            Loop::Coedges(cs) if cs.len() == 1 => {
                let edge = model.edge(model.coedge(cs[0]).edge);
                if edge.start.is_some() || edge.end.is_some() {
                    return Ok(None);
                }
                let Curve::Circle(c) = model.geom().curve(edge.curve) else {
                    return Ok(None);
                };
                // Incidence on a circular carrier and coaxiality already prove
                // this entire ring is its section, including the retained nappe.
                // Keep the cheap rational path and avoid eliminating algebraic
                // construction identities that the producer has recorded.
                let constructed_section = !surface.has_rational_definition()
                    && [FactKind::Incident, FactKind::Coaxial]
                        .into_iter()
                        .all(|kind| {
                            model.geom().facts().holds(
                                kind,
                                GeomKey::Curve(edge.curve),
                                GeomKey::Surface(model.face(face).surface),
                            )
                        });
                if constructed_section {
                    dec.note_fact_decision();
                } else {
                    if !zero(&c.center.minus(origin).cross(axis).norm_squared(), dec)?
                        || !zero(&c.normal.cross(axis).norm_squared(), dec)?
                    {
                        return Ok(None);
                    }
                    let expected = match &surface {
                        Surface::Cylinder(s) => s.radius_sq.clone(),
                        Surface::Cone(s) => {
                            let axial = c.center.minus(&s.apex).dot(&s.axis);
                            if sign(&axial, dec, "boolean.bound-ambiguous")? == Sign::Negative {
                                return Ok(None);
                            }
                            &s.slope.square() * &(&axial.square() / &s.axis.norm_squared())
                        }
                        _ => unreachable!(),
                    };
                    if !zero(&(&expected - &c.radius_sq), dec)? {
                        return Ok(None);
                    }
                }
                let radius = if let Some(q) = c.radius_sq.rational_ref() {
                    if *q == Rational::one() {
                        Real::one()
                    } else {
                        c.radius_sq.sqrt()
                    }
                } else {
                    c.radius_sq.sqrt()
                };
                ends.push((c.center.clone(), radius));
            }
            _ => return Ok(None),
        }
    }
    if ends.len() != 2 {
        return Ok(None);
    }
    // Two complete coaxial circles cut a cylinder/cone nappe into axial bands.
    // Proving the middle band is the trimmed face excludes the unbounded cells.
    let half = Real::from(Rational::from_ratio(1, 2).expect("fixed denominator"));
    let center = ends[0].0.plus(&ends[1].0.minus(&ends[0].0).scaled(&half));
    let radius = &(&ends[0].1 + &ends[1].1) * &half;
    let norm = reference.norm_squared();
    let length = if norm.rational_ref() == Some(&Rational::one()) {
        Real::one()
    } else {
        norm.sqrt()
    };
    let sample = center.plus(&reference.scaled(&(&radius / &length)));
    let point = model.geom_mut().add_point(PointDef::Explicit(sample));
    if locate(model, face, point, dec)? != Location::Inside {
        return Ok(None);
    }
    let a = sphere_box(&ends[0].0, &ends[0].1, dec)?;
    let b = sphere_box(&ends[1].0, &ends[1].1, dec)?;
    Ok(Some(std::array::from_fn(|i| a[i].hull(&b[i]))))
}

fn band_radius_upper(model: &Model, face: FaceId, dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut radius = Rational::zero();
    for lp in &model.face(face).loops {
        if let Loop::Coedges(coedges) = model.loop_(*lp) {
            for coedge in coedges {
                let edge = model.edge(model.coedge(*coedge).edge);
                let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
                    unreachable!("face_box certified a complete circular band");
                };
                radius = radius.max(enclose(&circle.radius_sq.sqrt(), dec)?.hi());
            }
        }
    }
    Ok(radius)
}

// The caller certified this cone face as the band between complete rings.
// Bound it by a tube and compare parallel axes, independent of world boxes.
fn cone_band_misses_surface(
    model: &mut Model,
    face: FaceId,
    other: FaceId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if !matches!(
        model.geom().surface(model.face(face).surface),
        Surface::Cone(_)
    ) {
        return Ok(false);
    }
    if matches!(
        model.geom().surface(model.face(other).surface),
        Surface::Cone(_)
    ) && face_box(model, other, false, dec)?.is_none()
    {
        return Ok(false);
    }
    let Surface::Cone(cone) = model.geom().surface(model.face(face).surface) else {
        unreachable!()
    };
    let (origin, axis, radius, cylinder_inner) =
        match model.geom().surface(model.face(other).surface) {
            Surface::Cylinder(cylinder) => {
                let radius = enclose(&cylinder.radius_sq.sqrt(), dec)?;
                (
                    &cylinder.origin,
                    &cylinder.axis,
                    radius.hi(),
                    Some(radius.lo()),
                )
            }
            Surface::Cone(other_cone) => (
                &other_cone.apex,
                &other_cone.axis,
                band_radius_upper(model, other, dec)?,
                None,
            ),
            _ => return Ok(false),
        };
    if !zero(&cone.axis.cross(axis).norm_squared(), dec)? {
        return Ok(false);
    }
    let distance_sq = cone.apex.minus(origin).cross(axis).norm_squared() / axis.norm_squared();
    let cone_radius = band_radius_upper(model, face, dec)?;
    let radii = &cone_radius + &radius;
    if enclose(&distance_sq, dec)?.lo() > &radii * &radii {
        return Ok(true);
    }
    // Only a cylinder has a constant-radius interior that excludes this tube.
    if let Some(inner) = cylinder_inner {
        return Ok(enclose(&distance_sq.sqrt(), dec)?.hi() + cone_radius < inner);
    }
    Ok(false)
}

pub(crate) struct FaceBounds {
    algebraic: bool,
    boxes: BTreeMap<FaceId, Option<[Ball; 3]>>,
}

impl FaceBounds {
    pub(crate) fn new(algebraic: bool) -> Self {
        Self {
            algebraic,
            boxes: BTreeMap::new(),
        }
    }

    pub(crate) fn separated(
        &mut self,
        model: &mut Model,
        a: FaceId,
        b: FaceId,
        dec: &mut Decider,
    ) -> bool {
        for face in [a, b] {
            if let std::collections::btree_map::Entry::Vacant(entry) = self.boxes.entry(face) {
                // Bounds are optional pruning proofs. Unsupported or unsettled
                // trims retain the ordinary intersection path.
                entry.insert(face_box(model, face, true, dec).ok().flatten());
            }
        }
        let (Some(a), Some(b)) = (&self.boxes[&a], &self.boxes[&b]) else {
            return false;
        };
        (0..3).any(|i| a[i].hi() < b[i].lo() || b[i].hi() < a[i].lo())
    }

    pub(crate) fn excludes(
        &mut self,
        model: &mut Model,
        a: FaceId,
        b: FaceId,
        dec: &mut Decider,
    ) -> bool {
        if self.separated(model, a, b, dec) {
            return true;
        }
        for (bounded, other) in [(a, b), (b, a)] {
            if let Some(bounds) = &self.boxes[&bounded]
                && let Ok(Some(value)) =
                    support_box(model.geom().surface(model.face(other).surface), bounds, dec)
                && !value.contains_zero()
            {
                return true;
            }
        }
        false
    }

    pub(crate) fn disjoint(
        &mut self,
        model: &mut Model,
        a: FaceId,
        b: FaceId,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for (bounded, other) in [(a, b), (b, a)] {
            // Retain the native rational fallback and its accounting. Broader
            // finite-plane bounds and memoization belong to the algebraic
            // construction path that performs support pruning up front.
            let bounds = if self.algebraic {
                if let std::collections::btree_map::Entry::Vacant(entry) = self.boxes.entry(bounded)
                {
                    entry.insert(face_box(model, bounded, true, dec)?);
                }
                self.boxes
                    .get(&bounded)
                    .expect("cached face bounds")
                    .clone()
            } else {
                face_box(model, bounded, false, dec)?
            };
            if let Some(bounds) = bounds {
                let Some(polynomial) = model.geom().surface(model.face(other).surface).implicit()
                else {
                    continue;
                };
                let value = polynomial.eval_ball(&bounds, dec)?;
                if !value.contains_zero() {
                    return Ok(true);
                }
                if cone_band_misses_surface(model, bounded, other, dec)? {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Pose, Vector3};
    use catallax_topology::{Face, Sense, cylinder};

    #[test]
    fn cone_tube_separation_preserves_tangency_and_crossings() {
        use catallax_shapes::{RVector3, Rotation};
        let r = Rational::from_i64;
        for (quaternion, offset) in [([1, 0, 0, 0], 0), ([1, 2, 3, 4], 1_000_000)] {
            let rotation = Rotation::from_quaternion(
                r(quaternion[0]),
                r(quaternion[1]),
                r(quaternion[2]),
                r(quaternion[3]),
            )
            .unwrap();
            let pose = Pose::from_rational(rotation, RVector3::new(r(offset), r(11), r(-13)));
            let mut model = Model::new();
            let cone = catallax_topology::cone(&mut model, r(4), r(2), r(2), &pose).unwrap();
            let cone_face = model
                .faces_of_body(cone)
                .into_iter()
                .find(|f| {
                    matches!(
                        model.geom().surface(model.face(*f).surface),
                        Surface::Cone(_)
                    )
                })
                .unwrap();
            let mut dec = Decider::default();
            assert!(
                face_box(&mut model, cone_face, false, &mut dec)
                    .unwrap()
                    .is_some()
            );
            for (distance, radius, excluded) in [
                (8, 3, true),
                (7, 3, false),
                (6, 3, false),
                (7, 12, true),
                (8, 12, false),
                (9, 12, false),
            ] {
                let mut placed = pose.clone();
                placed.translation = placed.translation.plus(&pose.apply_vector(&Vector3::new(
                    r(distance).into(),
                    Real::zero(),
                    Real::zero(),
                )));
                let cylinder = cylinder(&mut model, r(radius), r(2), &placed).unwrap();
                let wall = model
                    .faces_of_body(cylinder)
                    .into_iter()
                    .find(|f| {
                        matches!(
                            model.geom().surface(model.face(*f).surface),
                            Surface::Cylinder(_)
                        )
                    })
                    .unwrap();
                assert_eq!(
                    cone_band_misses_surface(&mut model, cone_face, wall, &mut dec).unwrap(),
                    excluded,
                    "distance {distance}, radius {radius}, quaternion {quaternion:?}"
                );
                if radius == 3 {
                    let other =
                        catallax_topology::cone(&mut model, r(3), r(1), r(2), &placed).unwrap();
                    let wall = model
                        .faces_of_body(other)
                        .into_iter()
                        .find(|f| {
                            matches!(
                                model.geom().surface(model.face(*f).surface),
                                Surface::Cone(_)
                            )
                        })
                        .unwrap();
                    assert_eq!(
                        cone_band_misses_surface(&mut model, cone_face, wall, &mut dec).unwrap(),
                        excluded,
                        "cone distance {distance}, quaternion {quaternion:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn support_exclusion_keeps_local_gaps_under_common_translation() {
        let r = Rational::from_i64;
        for offset in [0, 1_000_000] {
            let mut model = Model::new();
            let mut pose = Pose::identity();
            pose.translation.x = Real::from_i64(offset);
            let body = cylinder(&mut model, r(5), r(2), &pose).unwrap();
            let cap = model.faces_of_body(body)[0];
            let mut bounds = FaceBounds::new(false);
            let mut dec = Decider::default();
            for (distance, excluded) in [(15, true), (10, false), (9, false)] {
                pose.translation.x = Real::from_i64(offset + distance);
                let other = cylinder(&mut model, r(5), r(2), &pose).unwrap();
                let wall = model
                    .faces_of_body(other)
                    .into_iter()
                    .find(|face| {
                        matches!(
                            model.geom().surface(model.face(*face).surface),
                            Surface::Cylinder(_)
                        )
                    })
                    .unwrap();
                let support = model.add_face(Face {
                    surface: model.face(wall).surface,
                    sense: Sense::Same,
                    loops: vec![],
                });
                assert_eq!(
                    bounds.excludes(&mut model, cap, support, &mut dec),
                    excluded,
                    "offset {offset}, distance {distance}"
                );
            }
        }
    }

    #[test]
    fn coincident_caps_prune_only_strictly_separated_finite_trims() {
        let r = Rational::from_i64;
        let mut model = Model::new();
        let first = cylinder(&mut model, r(2), r(3), &Pose::identity()).unwrap();
        let cap = model.faces_of_body(first)[0];
        let mut bounds = FaceBounds::new(false);
        let mut dec = Decider::default();
        for (x, separated) in [(5, true), (4, false), (3, false), (0, false)] {
            let mut pose = Pose::identity();
            pose.translation.x = Real::from_i64(x);
            let other = cylinder(&mut model, r(2), r(3), &pose).unwrap();
            let other_cap = model.faces_of_body(other)[0];
            assert_eq!(
                bounds.separated(&mut model, cap, other_cap, &mut dec),
                separated,
                "cap offset {x}"
            );
        }
        let unbounded = model.add_face(Face {
            surface: model.face(cap).surface,
            sense: model.face(cap).sense,
            loops: vec![],
        });
        assert!(!bounds.separated(&mut model, cap, unbounded, &mut dec));
    }

    #[test]
    fn planar_bounds_exclude_only_supports_outside_the_finite_disk() {
        let r = Rational::from_i64;
        let mut model = Model::new();
        let body = cylinder(&mut model, r(2), r(3), &Pose::identity()).unwrap();
        let cap = model
            .faces_of_body(body)
            .into_iter()
            .find(|f| {
                matches!(
                    model.geom().surface(model.face(*f).surface),
                    Surface::Plane(_)
                )
            })
            .unwrap();
        let mut dec = Decider::default();
        let mut bounds = FaceBounds::new(true);
        for (x, excluded) in [(3, true), (2, false), (1, false)] {
            let surface = model.geom_mut().add_surface(Surface::Plane(Plane {
                origin: Point3::new(Real::from_i64(x), Real::zero(), Real::zero()),
                normal: Vector3::new(Real::one(), Real::zero(), Real::zero()),
                ref_dir: Vector3::new(Real::zero(), Real::one(), Real::zero()),
            }));
            let other = model.add_face(Face {
                surface,
                sense: Sense::Same,
                loops: vec![],
            });
            assert_eq!(
                bounds.disjoint(&mut model, cap, other, &mut dec).unwrap(),
                excluded
            );
        }
    }
}
