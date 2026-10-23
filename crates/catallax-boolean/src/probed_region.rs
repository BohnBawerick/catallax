//! A surface probe follows a section through the query to its first
//! transverse boundary hit. Incidence, arc membership and order all use
//! the earlier pieces' public predicates.

use crate::geometry::{coincident, natural_ends, on_arc, oriented, product, refusal, sign};
use crate::region::Location;
use catallax_intersect::{
    Component, CurveForm, CurveSurface, Cyclic, SurfacePair, curve_form, curve_is_closed,
    cyclic_order, intersect_curve_surface, intersect_surfaces, on_curve, order_on_curve,
    tangent_at,
};
use catallax_numbers::{
    Decider, Fact, FactKind, FactOrigin, Outcome, Refusal, Sign, StrictOrder, ZeroProof,
};
use catallax_shapes::{
    Curve, CurveId, GeomKey, MeetOf, Plane, Point3, PointDef, PointId, RVector3, Surface,
    SurfaceId, Vector3,
};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};
use std::collections::{BTreeMap, BTreeSet};

const PROBE_CONSTRUCTION: &str = "boolean.surface-probe";

fn plane_through(
    model: &mut Model,
    incident: GeomKey,
    origin: Point3,
    normal: Vector3,
    dec: &mut Decider,
) -> Result<Option<SurfaceId>, Refusal> {
    if sign(&normal.norm_squared(), dec, "boolean.classify-ambiguous")? == Sign::Zero {
        return Ok(None);
    }
    let reference = crate::geometry::perpendicular(&normal, dec, "boolean.classify-ambiguous")?;
    let surface = model
        .geom_mut()
        .add_surface(Surface::Plane(Plane::new(origin, normal, reference, dec)?));
    // Every caller constructs the plane through this point or curve exactly.
    model.geom_mut().record_fact(Fact {
        kind: FactKind::Incident,
        a: incident,
        b: GeomKey::Surface(surface),
        at: None,
        origin: FactOrigin::Construction(PROBE_CONSTRUCTION),
    });
    Ok(Some(surface))
}
fn auxiliary_geometry(model: &Model) -> (BTreeSet<SurfaceId>, BTreeSet<CurveId>) {
    let facts = model.geom().facts();
    let mut planes = BTreeSet::new();
    let mut inputs = BTreeSet::new();
    for fact in facts.iter() {
        if fact.origin == FactOrigin::Construction(PROBE_CONSTRUCTION) {
            if let GeomKey::Surface(surface) = fact.b {
                planes.insert(surface);
            }
            if let GeomKey::Curve(curve) = fact.a {
                inputs.insert(curve);
            }
        }
    }
    // Curves used to construct probe planes remain inputs. Only sections
    // constructed from those planes are auxiliary probe geometry.
    let curves = facts
        .iter()
        .filter_map(|fact| match (fact.kind, fact.a, fact.b) {
            (FactKind::Incident, GeomKey::Curve(curve), GeomKey::Surface(surface))
                if planes.contains(&surface) && !inputs.contains(&curve) =>
            {
                Some(curve)
            }
            _ => None,
        })
        .collect();
    (planes, curves)
}
fn boundary(model: &Model, face: FaceId) -> Vec<(EdgeId, Sense)> {
    let mut uses: BTreeMap<EdgeId, i32> = BTreeMap::new();
    for &l in &model.face(face).loops {
        if let Loop::Coedges(cs) = model.loop_(l) {
            for &c in cs {
                let c = model.coedge(c);
                *uses.entry(c.edge).or_default() += if c.sense == Sense::Same { 1 } else { -1 };
            }
        }
    }
    uses.into_iter()
        .filter_map(|(edge, count)| {
            if count > 0 {
                Some((edge, Sense::Same))
            } else if count < 0 {
                Some((edge, Sense::Reversed))
            } else {
                None
            }
        })
        .collect()
}
fn first_hit(
    model: &mut Model,
    face: FaceId,
    point: PointId,
    probe_surface: SurfaceId,
    curve: CurveId,
    edges: &[(EdgeId, Sense)],
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let line = matches!(
        curve_form(model.geom(), curve, dec)?,
        CurveForm::Line { .. }
    );
    let open = !curve_is_closed(model.geom(), curve, dec)?;
    let f = model.face(face).clone();
    let mut first: Option<(PointId, Sign)> = None;
    for &(edge_id, use_sense) in edges {
        let edge = model.edge(edge_id).clone();
        let (from, to) = natural_ends(model, edge_id);
        let hits = if line {
            // A generator already lies on the face. Intersect that line with
            // the boundary carrier instead of solving the branch against the
            // entire probe plane, which unnecessarily produces a quartic.
            match crate::curve_pair::meetings(model, curve, edge.curve, dec).map_err(
                |mut error| {
                    if error.code == "boolean.imprint-ambiguous" {
                        error.code = "boolean.classify-ambiguous";
                    }
                    error
                },
            )? {
                crate::curve_pair::Meetings::Coincident => return Ok(None),
                crate::curve_pair::Meetings::Points(points) => points,
            }
        } else {
            let answer = intersect_curve_surface(model.geom_mut(), edge.curve, probe_surface, dec)?;
            if !answer.marks.is_empty() {
                return Err(refusal(
                    "boolean.classify-ambiguous",
                    "surface probe required unresolved curve incidence",
                ));
            }
            match answer.value {
                CurveSurface::Disjoint => continue,
                CurveSurface::Contained => return Ok(None),
                CurveSurface::Meets(hits) => {
                    hits.into_iter().map(|candidate| candidate.point).collect()
                }
            }
        };
        for hit in hits {
            if !coincident(
                on_curve(model.geom(), curve, hit, dec),
                "boolean.classify-ambiguous",
            )? || !on_arc(
                model,
                edge.curve,
                from,
                to,
                hit,
                dec,
                "boolean.classify-ambiguous",
            )? {
                continue;
            }
            if coincident(
                model.geom().same_point(point, hit, dec),
                "boolean.classify-ambiguous",
            )? {
                return Ok(Some(Location::Boundary));
            }
            let mut vertex = false;
            for endpoint in [from, to].into_iter().flatten() {
                if coincident(
                    model.geom().same_point(endpoint, hit, dec),
                    "boolean.classify-ambiguous",
                )? {
                    vertex = true;
                    break;
                }
            }
            if vertex {
                return Ok(None);
            }
            if open {
                match order_on_curve(model.geom(), curve, point, hit, dec) {
                    Outcome::Proven(StrictOrder::Less) => {}
                    Outcome::Proven(_) => continue,
                    Outcome::ProvenCoincident(_) => return Ok(Some(Location::Boundary)),
                    Outcome::NotExact(n) => {
                        return Err(refusal("boolean.classify-ambiguous", n.detail));
                    }
                    Outcome::Refused(e) => return Err(e),
                }
            }
            let probe_tangent = tangent_at(model.geom(), curve, hit, dec)?[0].clone();
            let boundary_tangent = oriented(
                &tangent_at(model.geom(), edge.curve, hit, dec)?[1],
                product(edge.sense, use_sense),
            );
            let normal = oriented(&model.geom().normal(f.surface, hit, dec)?, f.sense);
            let crossing = sign(
                &probe_tangent.dot(&normal.cross(&boundary_tangent)),
                dec,
                "boolean.classify-ambiguous",
            )?;
            if crossing == Sign::Zero {
                continue;
            }
            let earlier = if let Some((previous, _)) = first {
                if coincident(
                    model.geom().same_point(previous, hit, dec),
                    "boolean.classify-ambiguous",
                )? {
                    return Ok(None);
                }
                if open {
                    match order_on_curve(model.geom(), curve, hit, previous, dec) {
                        Outcome::Proven(StrictOrder::Less) => true,
                        Outcome::Proven(_) => false,
                        Outcome::ProvenCoincident(_) => return Ok(None),
                        Outcome::NotExact(n) => {
                            return Err(refusal("boolean.classify-ambiguous", n.detail));
                        }
                        Outcome::Refused(e) => return Err(e),
                    }
                } else {
                    match cyclic_order(model.geom(), curve, [point, hit, previous], dec) {
                        Outcome::Proven(Cyclic::Forward) => true,
                        Outcome::Proven(Cyclic::Backward) => false,
                        Outcome::ProvenCoincident(_) => return Ok(None),
                        Outcome::NotExact(n) => {
                            return Err(refusal("boolean.classify-ambiguous", n.detail));
                        }
                        Outcome::Refused(e) => return Err(e),
                    }
                }
            } else {
                true
            };
            if earlier {
                first = Some((hit, crossing));
            }
        }
    }
    Ok(first.map(|(_, crossing)| {
        if crossing == Sign::Negative {
            Location::Inside
        } else {
            Location::Outside
        }
    }))
}
fn try_surface(
    model: &mut Model,
    face: FaceId,
    point: PointId,
    surface: SurfaceId,
    edges: &[(EdgeId, Sense)],
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    if !matches!(model.geom().surface(surface), Surface::Plane(_)) {
        return Ok(None);
    }
    let face_surface = model.face(face).surface;
    let pair = intersect_surfaces(model.geom_mut(), face_surface, surface, dec)?;
    if !pair.marks.is_empty() {
        return Err(refusal(
            "boolean.classify-ambiguous",
            "surface probe section required an unresolved decision",
        ));
    }
    let SurfacePair::Meets(parts) = pair.value else {
        return Ok(None);
    };
    if let [Component::Curve { curve, .. }] = parts.as_slice()
        && model.geom().facts().holds(
            FactKind::Incident,
            GeomKey::Point(point),
            GeomKey::Surface(face_surface),
        )
        && model.geom().facts().holds(
            FactKind::Incident,
            GeomKey::Point(point),
            GeomKey::Surface(surface),
        )
    {
        // The unmarked pair exhausts the common zero set. With one component,
        // both surface incidences identify its curve, without eliminating a
        // second algebraic coordinate representation of the same point.
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Curve(*curve),
            at: None,
            origin: FactOrigin::Proof(ZeroProof::Algebraic),
        });
    }
    for part in parts {
        if let Component::Curve { curve, .. } = part
            && coincident(
                on_curve(model.geom(), curve, point, dec),
                "boolean.classify-ambiguous",
            )?
            && let Some(location) = first_hit(model, face, point, surface, curve, edges, dec)?
        {
            return Ok(Some(location));
        }
    }
    Ok(None)
}
fn known(model: &Model, point: PointId) -> (Vec<SurfaceId>, Vec<CurveId>) {
    let (auxiliary_planes, auxiliary_curves) = auxiliary_geometry(model);
    let mut surfaces = Vec::new();
    let mut curves = Vec::new();
    if let PointDef::Meet(m) = model.geom().point(point) {
        match &m.of {
            MeetOf::Surfaces(s) => surfaces.extend(s),
            MeetOf::CurveSurface(c, s) => {
                curves.push(*c);
                surfaces.push(*s);
            }
            MeetOf::Curves(a, b) => curves.extend([a, b]),
        }
    }
    for fact in model.geom().facts().facts_about(GeomKey::Point(point)) {
        if fact.kind == FactKind::Incident && fact.a == GeomKey::Point(point) {
            match fact.b {
                GeomKey::Curve(c) => curves.push(c),
                GeomKey::Surface(s) => surfaces.push(s),
                _ => {}
            }
        }
    }
    // Keep proof facts, but do not feed derived probe geometry back into the
    // next query's seed set. Its field can be more expensive than the carrier's.
    curves.retain(|curve| !auxiliary_curves.contains(curve));
    let mut seen = BTreeSet::new();
    let mut index = 0;
    while index < curves.len() {
        let curve = curves[index];
        index += 1;
        if !seen.insert(curve) {
            continue;
        }
        if let Curve::Intersection(c) = model.geom().curve(curve) {
            surfaces.extend(c.surfaces);
        }
        for fact in model.geom().facts().facts_about(GeomKey::Curve(curve)) {
            if fact.kind == FactKind::Incident
                && fact.a == GeomKey::Curve(curve)
                && let GeomKey::Surface(s) = fact.b
            {
                surfaces.push(s);
            }
            if fact.kind == FactKind::Coincident {
                let other = if fact.a == GeomKey::Curve(curve) {
                    fact.b
                } else {
                    fact.a
                };
                if let GeomKey::Curve(c) = other
                    && !seen.contains(&c)
                    && !auxiliary_curves.contains(&c)
                {
                    curves.push(c);
                }
            }
        }
    }
    surfaces.retain(|surface| !auxiliary_planes.contains(surface));
    surfaces.sort();
    surfaces.dedup();
    curves.sort();
    curves.dedup();
    (surfaces, curves)
}

pub(crate) fn locate(
    model: &mut Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<Option<Location>, Refusal> {
    let edges = boundary(model, face);
    if edges.is_empty() {
        return Ok(Some(Location::Inside));
    }
    let (surfaces, curves) = known(model, point);
    for surface in surfaces {
        if !coincident(
            model
                .geom()
                .same_surface(surface, model.face(face).surface, dec),
            "boolean.classify-ambiguous",
        )? && let Some(location) = try_surface(model, face, point, surface, &edges, dec)?
        {
            return Ok(Some(location));
        }
    }
    let position = model.geom().position(point);
    let exact_probe = position.as_rational().is_some()
        || !model
            .geom()
            .surface(model.face(face).surface)
            .has_rational_definition();
    let attempts = crate::geometry::DIRECTION_ATTEMPTS
        .min((edges.len() as u32).saturating_mul(4).saturating_add(4));
    if exact_probe {
        let p = position;
        let axial = match model.geom().surface(model.face(face).surface) {
            Surface::Cylinder(c) => Some((c.origin.clone(), c.axis.clone())),
            Surface::EllipticCylinder(c) => Some((c.origin.clone(), c.axis.clone())),
            _ => None,
        };
        if let Some((origin, axis)) = axial {
            let normal = axis.cross(&p.minus(&origin));
            if let Some(surface) = plane_through(model, GeomKey::Point(point), origin, normal, dec)?
                && let Some(location) = try_surface(model, face, point, surface, &edges, dec)?
            {
                return Ok(Some(location));
            }
        }
        let normal = model
            .geom()
            .surface(model.face(face).surface)
            .normal_at(&p, dec)?;
        let u = crate::geometry::perpendicular(&normal, dec, "boolean.classify-ambiguous")?;
        let v = normal.cross(&u);
        for index in 0..attempts {
            let (a, b) = crate::geometry::probe_parameters(index);
            let n = u.plus(&v.scaled(&a.into())).plus(&normal.scaled(&b.into()));
            if let Some(surface) = plane_through(model, GeomKey::Point(point), p.clone(), n, dec)?
                && let Some(location) = try_surface(model, face, point, surface, &edges, dec)?
            {
                return Ok(Some(location));
            }
        }
    } else {
        for curve in curves {
            match model.geom().curve(curve).clone() {
                Curve::Line(line) => {
                    for index in 0..attempts {
                        let (a, b) = crate::geometry::probe_parameters(index);
                        let direction = RVector3::new(crate::geometry::r(1), a, b).to_real();
                        let normal = line.dir.cross(&direction);
                        if let Some(surface) = plane_through(
                            model,
                            GeomKey::Curve(curve),
                            line.origin.clone(),
                            normal,
                            dec,
                        )? && let Some(location) =
                            try_surface(model, face, point, surface, &edges, dec)?
                        {
                            return Ok(Some(location));
                        }
                    }
                }
                Curve::Circle(c) => {
                    if let Some(surface) =
                        plane_through(model, GeomKey::Curve(curve), c.center, c.normal, dec)?
                        && let Some(location) =
                            try_surface(model, face, point, surface, &edges, dec)?
                    {
                        return Ok(Some(location));
                    }
                }
                Curve::Ellipse(c) => {
                    if let Some(surface) =
                        plane_through(model, GeomKey::Curve(curve), c.center, c.normal, dec)?
                        && let Some(location) =
                            try_surface(model, face, point, surface, &edges, dec)?
                    {
                        return Ok(Some(location));
                    }
                }
                _ => {}
            }
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::unit;
    use catallax_numbers::{Rational, StrictSign};
    use catallax_shapes::{Pose, RPoint3};

    #[test]
    fn multi_component_probe_does_not_infer_a_branch() {
        for x in [-1, 1] {
            let mut model = Model::new();
            let body = catallax_topology::cylinder(
                &mut model,
                Rational::one(),
                Rational::from_i64(2),
                &Pose::identity(),
            )
            .unwrap();
            let face = model
                .faces_of_body(body)
                .into_iter()
                .find(|face| {
                    matches!(
                        model.geom().surface(model.face(*face).surface),
                        Surface::Cylinder(_)
                    )
                })
                .unwrap();
            let carrier = model.face(face).surface;
            let position = RPoint3::new(Rational::from_i64(x), Rational::zero(), Rational::one());
            let point = model
                .geom_mut()
                .add_point(PointDef::Explicit((position.clone()).into()));
            let mut dec = Decider::default();
            let Outcome::ProvenCoincident(proof) = model.geom().side(carrier, point, &mut dec)
            else {
                panic!("the rational query lies on the cylinder");
            };
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Point(point),
                b: GeomKey::Surface(carrier),
                at: None,
                origin: FactOrigin::Proof(proof.proof),
            });
            let plane = plane_through(
                &mut model,
                GeomKey::Point(point),
                RPoint3::new(Rational::zero(), Rational::zero(), Rational::zero()).into(),
                unit(1),
                &mut dec,
            )
            .unwrap()
            .unwrap();
            let edges = boundary(&model, face);
            assert!(matches!(
                try_surface(&mut model, face, point, plane, &edges, &mut dec).unwrap(),
                Some(Location::Inside)
            ));
            let generators: BTreeSet<_> = model
                .geom()
                .facts()
                .iter()
                .filter_map(|fact| match fact.a {
                    GeomKey::Curve(curve) if fact.kind == FactKind::Incident => Some(curve),
                    _ => None,
                })
                .filter(|curve| {
                    model.geom().facts().holds(
                        FactKind::Incident,
                        GeomKey::Curve(*curve),
                        GeomKey::Surface(carrier),
                    ) && model.geom().facts().holds(
                        FactKind::Incident,
                        GeomKey::Curve(*curve),
                        GeomKey::Surface(plane),
                    )
                })
                .collect();
            let mut memberships = BTreeSet::new();
            for curve in generators {
                let CurveForm::Line { origin, dir } =
                    curve_form(model.geom(), curve, &mut dec).unwrap()
                else {
                    panic!("a meridian cuts a cylinder in two generators");
                };
                let distance = position.to_real().minus(&origin).cross(&dir).norm_squared();
                let on = match dec.sign(&distance) {
                    Outcome::ProvenCoincident(_) => true,
                    Outcome::Proven(StrictSign::Positive) => false,
                    other => panic!("line distance must be decided exactly: {other:?}"),
                };
                let got = on_curve(model.geom(), curve, point, &mut dec);
                assert!(
                    matches!(
                        (on, got),
                        (true, Outcome::ProvenCoincident(_)) | (false, Outcome::Proven(_))
                    ),
                    "a parent-surface fact must not choose the wrong component"
                );
                memberships.insert(on);
            }
            assert_eq!(memberships, BTreeSet::from([false, true]));
        }
    }
}
