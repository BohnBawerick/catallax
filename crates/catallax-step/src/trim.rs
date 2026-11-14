//! Exact trim queries used by export-only periodic charts. The scratch store
//! shares the model's immutable geometry; no model topology or geometry changes.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use catallax_intersect::{Contact, CurveForm, CurveHit, CurveSurface, Cyclic};
use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictOrder, StrictSign};
use catallax_shapes::{
    BSplineCurve, Curve, CurveId, GeomStore, Plane, Point3, PointDef, PointId, Surface, SurfaceId,
    Vector3,
};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};

use crate::number;

/// Blind query planes tried after every chord-aimed plane has failed.
const BLIND_ATTEMPTS: usize = 32;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Location {
    Outside,
    Boundary,
    Inside,
}

#[derive(Clone)]
pub(crate) struct Boundary {
    pub edge: EdgeId,
    pub curve: CurveId,
    pub from: Option<PointId>,
    pub to: Option<PointId>,
    pub closed: bool,
    pub forward: bool,
}

pub(crate) struct Geometry {
    pub store: GeomStore,
    planes: Vec<([Real; 4], SurfaceId)>,
    framed: bool,
    hits: BTreeMap<(CurveId, SurfaceId), Rc<CurveSurface>>,
    forms: BTreeMap<CurveId, Rc<CurveForm>>,
    closed: BTreeMap<CurveId, bool>,
    boundaries: BTreeMap<FaceId, Rc<Vec<Boundary>>>,
    parameters: BTreeMap<(CurveId, PointId), Real>,
}

pub(crate) fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    }
}

pub(crate) fn try_sort<T>(
    values: &mut [T],
    mut compare: impl FnMut(&T, &T) -> Result<Ordering, Refusal>,
) -> Result<(), Refusal> {
    fn sift<T>(
        values: &mut [T],
        mut root: usize,
        end: usize,
        compare: &mut impl FnMut(&T, &T) -> Result<Ordering, Refusal>,
    ) -> Result<(), Refusal> {
        while root < end / 2 {
            let mut child = root * 2 + 1;
            if child + 1 < end && compare(&values[child], &values[child + 1])? == Ordering::Less {
                child += 1;
            }
            if compare(&values[root], &values[child])? != Ordering::Less {
                break;
            }
            values.swap(root, child);
            root = child;
        }
        Ok(())
    }
    let len = values.len();
    for root in (0..len / 2).rev() {
        sift(values, root, len, &mut compare)?;
    }
    for end in (1..len).rev() {
        values.swap(0, end);
        sift(values, 0, end, &mut compare)?;
    }
    Ok(())
}

pub(crate) fn fillet_face(model: &Model, face: FaceId) -> bool {
    model.face_origins(face).iter().any(|origin| {
        let mut role = &origin.role;
        while let catallax_recipe::Role::Copy { of, .. } = role {
            role = &of.role;
        }
        matches!(role, catallax_recipe::Role::Blend { .. })
    })
}

impl Geometry {
    pub fn new(model: &Model, framed: bool) -> Self {
        Self {
            store: model.geom().clone(),
            planes: Vec::new(),
            framed,
            hits: BTreeMap::new(),
            forms: BTreeMap::new(),
            closed: BTreeMap::new(),
            boundaries: BTreeMap::new(),
            parameters: BTreeMap::new(),
        }
    }

    pub fn point(&mut self, point: Point3) -> PointId {
        self.store.add_point(PointDef::Explicit(point))
    }

    pub fn plane(
        &mut self,
        origin: Point3,
        normal: Vector3,
        reference: Vector3,
        dec: &mut Decider,
    ) -> Result<SurfaceId, Refusal> {
        dec.check_deadline()?;
        let terms = [
            normal.x.clone(),
            normal.y.clone(),
            normal.z.clone(),
            -origin.minus(&Point3::origin()).dot(&normal),
        ];
        let mut pivot = None;
        for component in &terms[..3] {
            if sign(component, dec)? != Sign::Zero {
                pivot = Some(component);
                break;
            }
        }
        let pivot = pivot.ok_or_else(|| number::budget("chart plane normal is not separated"))?;
        let key = terms.each_ref().map(|value| value / pivot);
        // Only construction equality is used for this optional cache. Do not
        // re-eliminate unrelated algebraic planes merely to discover an alias.
        if let Some((_, id)) = self
            .planes
            .iter()
            .find(|(old, _)| old.iter().zip(&key).all(|(a, b)| a.same_definition(b)))
        {
            return Ok(*id);
        }
        let id = self.store.add_surface(Surface::Plane(Plane {
            origin,
            normal,
            ref_dir: reference,
        }));
        self.planes.push((key, id));
        Ok(id)
    }

    pub fn intersections(
        &mut self,
        curve: CurveId,
        plane: SurfaceId,
        dec: &mut Decider,
    ) -> Result<Rc<CurveSurface>, Refusal> {
        if let Some(result) = self.hits.get(&(curve, plane)) {
            return Ok(result.clone());
        }
        let intersect = if self.framed {
            catallax_intersect::intersect_curve_surface_framed
        } else {
            catallax_intersect::intersect_curve_surface
        };
        let result = match intersect(&mut self.store, curve, plane, dec) {
            Ok(result) => result.value,
            Err(error) if error.code == "intersect.spline-plane" => {
                let (Curve::BSpline(spline), Surface::Plane(cut)) = (
                    self.store.curve(curve).clone(),
                    self.store.surface(plane).clone(),
                ) else {
                    return Err(error);
                };
                self.spline_plane(curve, &spline, &cut, dec)?
            }
            Err(error) => return Err(error),
        };
        let result = Rc::new(result);
        self.hits.insert((curve, plane), result.clone());
        Ok(result)
    }

    // Algebraic finite spline trims need the same exact polynomial roots as
    // rational trims. Keep every knot interval; a partially contained spline
    // has no single CurveSurface result and cannot certify a chart cut.
    fn spline_plane(
        &mut self,
        curve: CurveId,
        spline: &BSplineCurve,
        cut: &Plane,
        dec: &mut Decider,
    ) -> Result<CurveSurface, Refusal> {
        let mut hits = Vec::<CurveHit>::new();
        let mut contained = 0;
        let mut spans = 0;
        let before = dec.stats().coincident_separation;
        let normal = [&cut.normal.x, &cut.normal.y, &cut.normal.z];
        let offset = cut.origin.minus(&Point3::origin()).dot(&cut.normal);
        for (span, [lo, hi]) in spline.spans() {
            spans += 1;
            let homogeneous = spline.query_span_polynomials(span, dec)?;
            let equation = (0..3).fold(catallax_numbers::RealPoly::zero(1), |sum, i| {
                sum + homogeneous[i].scale(normal[i])
            }) - homogeneous[3].scale(&offset);
            if equation.is_zero(dec)? {
                contained += 1;
                continue;
            }
            let roots = equation.real_roots(dec)?;
            for parameter in roots {
                let lower = sign(&(&parameter - &Real::from(lo.clone())), dec)?;
                let upper = sign(&(&Real::from(hi.clone()) - &parameter), dec)?;
                if lower == Sign::Negative || upper == Sign::Negative {
                    continue;
                }
                // A knot needs both one-sided tangents. Conservatively reject
                // it as a chart crossing; another meridian can avoid it.
                let tangent = lower == Sign::Zero
                    || upper == Sign::Zero
                    || sign(
                        &equation
                            .partial(0)
                            .eval_real(std::slice::from_ref(&parameter)),
                        dec,
                    )? == Sign::Zero;
                let point = self.point(spline.point_at(&parameter, dec)?);
                let mut duplicate = false;
                for hit in &hits {
                    duplicate |= sign(&(&parameter - &self.parameters[&(curve, hit.point)]), dec)?
                        == Sign::Zero;
                }
                if !duplicate {
                    self.parameters.insert((curve, point), parameter);
                    hits.push(CurveHit {
                        point,
                        contact: if tangent {
                            Contact::Tangent
                        } else {
                            Contact::Crossing
                        },
                    });
                }
            }
        }
        if dec.stats().coincident_separation != before {
            return Err(number::budget("spline chart cut needs exact finite roots"));
        }
        if contained == spans {
            Ok(CurveSurface::Contained)
        } else if contained > 0 {
            Err(number::budget(
                "a spline chart cut contains only some knot spans",
            ))
        } else if hits.is_empty() {
            Ok(CurveSurface::Disjoint)
        } else {
            Ok(CurveSurface::Meets(hits))
        }
    }

    pub fn form(&mut self, curve: CurveId, dec: &mut Decider) -> Result<Rc<CurveForm>, Refusal> {
        if let Some(form) = self.forms.get(&curve) {
            return Ok(form.clone());
        }
        let form = Rc::new(catallax_intersect::curve_form(&self.store, curve, dec)?);
        self.forms.insert(curve, form.clone());
        Ok(form)
    }

    pub fn closed(&mut self, curve: CurveId, dec: &mut Decider) -> Result<bool, Refusal> {
        if let Some(closed) = self.closed.get(&curve) {
            return Ok(*closed);
        }
        let closed = catallax_intersect::curve_is_closed(&self.store, curve, dec)?;
        self.closed.insert(curve, closed);
        Ok(closed)
    }

    pub fn same_point(&self, a: PointId, b: PointId, dec: &mut Decider) -> Result<bool, Refusal> {
        if a == b {
            return Ok(true);
        }
        match self.store.same_point(a, b, dec) {
            Outcome::Proven(_) => Ok(false),
            Outcome::ProvenCoincident(_) => Ok(true),
            Outcome::NotExact(n) => Err(Refusal::from(n)),
            Outcome::Refused(r) => Err(r),
        }
    }

    pub fn boundaries(
        &mut self,
        model: &Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Rc<Vec<Boundary>>, Refusal> {
        if let Some(boundaries) = self.boundaries.get(&face) {
            return Ok(boundaries.clone());
        }
        let mut boundaries = Vec::new();
        for loop_id in &model.face(face).loops {
            let Loop::Coedges(coedges) = model.loop_(*loop_id) else {
                continue;
            };
            for coedge_id in coedges {
                let coedge = model.coedge(*coedge_id);
                let edge = model.edge(coedge.edge);
                let from = edge.start.map(|v| model.vertex(v).point);
                let to = edge.end.map(|v| model.vertex(v).point);
                let (from, to) = if edge.sense == Sense::Same {
                    (from, to)
                } else {
                    (to, from)
                };
                let closed = self.closed(edge.curve, dec)?;
                let forward = ((edge.sense == Sense::Same) == (coedge.sense == Sense::Same))
                    == (model.face(face).sense == Sense::Same);
                boundaries.push(Boundary {
                    edge: coedge.edge,
                    curve: edge.curve,
                    from,
                    to,
                    closed,
                    forward,
                });
            }
        }
        let boundaries = Rc::new(boundaries);
        self.boundaries.insert(face, boundaries.clone());
        Ok(boundaries)
    }

    pub fn on_arc(
        &self,
        edge: &Boundary,
        point: PointId,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        if let Some(parameter) = self.parameters.get(&(edge.curve, point))
            && !edge.closed
            && let Curve::BSpline(spline) = self.store.curve(edge.curve)
            && let (Some(from), Some(to)) = (edge.from, edge.to)
        {
            let bound = |point, dec: &mut Decider| -> Result<Real, Refusal> {
                let position = self.store.position(point);
                let [lo, hi] = spline.domain();
                if sign(&position.minus(spline.control(0)).norm_squared(), dec)? == Sign::Zero {
                    Ok(Real::from(lo.clone()))
                } else if sign(
                    &position
                        .minus(spline.control(spline.control_count() - 1))
                        .norm_squared(),
                    dec,
                )? == Sign::Zero
                {
                    Ok(Real::from(hi.clone()))
                } else {
                    Ok(spline
                        .unique_parameter_of(&position, dec)?
                        .ok_or_else(|| number::budget("spline trim endpoint has no parameter"))?
                        .value()
                        .clone())
                }
            };
            let from = sign(&(parameter - &bound(from, dec)?), dec)?;
            let to = sign(&(&bound(to, dec)? - parameter), dec)?;
            return Ok(if from == Sign::Negative || to == Sign::Negative {
                Location::Outside
            } else if from == Sign::Zero || to == Sign::Zero {
                Location::Boundary
            } else {
                Location::Inside
            });
        }
        let (Some(from), Some(to)) = (edge.from, edge.to) else {
            assert!(edge.from.is_none() && edge.to.is_none() && edge.closed);
            return Ok(Location::Inside);
        };
        if self.same_point(from, point, dec)? || self.same_point(to, point, dec)? {
            return Ok(Location::Boundary);
        }
        if edge.closed {
            if self.same_point(from, to, dec)? {
                return Ok(Location::Inside);
            }
            match catallax_intersect::cyclic_order(&self.store, edge.curve, [from, point, to], dec)
            {
                Outcome::Proven(Cyclic::Forward) => Ok(Location::Inside),
                Outcome::Proven(Cyclic::Backward) => Ok(Location::Outside),
                Outcome::ProvenCoincident(_) => Err(number::budget("trim arc order is not strict")),
                Outcome::NotExact(n) => Err(Refusal::from(n)),
                Outcome::Refused(r) => Err(r),
            }
        } else {
            for (a, b) in [(from, point), (point, to)] {
                match catallax_intersect::order_on_curve(&self.store, edge.curve, a, b, dec) {
                    Outcome::Proven(StrictOrder::Less) => {}
                    Outcome::Proven(StrictOrder::Greater) => return Ok(Location::Outside),
                    Outcome::ProvenCoincident(_) => return Ok(Location::Boundary),
                    Outcome::NotExact(n) => return Err(Refusal::from(n)),
                    Outcome::Refused(r) => return Err(r),
                }
            }
            Ok(Location::Inside)
        }
    }

    pub fn tangent(
        &self,
        curve: CurveId,
        point: PointId,
        dec: &mut Decider,
    ) -> Result<Vector3, Refusal> {
        if let Some(parameter) = self.parameters.get(&(curve, point))
            && let Curve::BSpline(spline) = self.store.curve(curve)
        {
            return spline.first_at(parameter, dec);
        }
        Ok(catallax_intersect::tangent_at(&self.store, curve, point, dec)?[1].clone())
    }

    /// Locates `point` on a spherical face, or `None` without a certificate.
    // A sphere face bounded only by whole sphere band rails, decided as
    // `catallax_topology::band_cap` places the point: a rail walked forward
    // keeps the face in the cap its band's turn axis points to.
    fn band_caps(
        &mut self,
        model: &Model,
        boundaries: &[Boundary],
        query: &Point3,
        center: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<Location>, Refusal> {
        let mut inside = true;
        for boundary in boundaries {
            if !boundary.closed {
                return Ok(None);
            }
            let Some(above) =
                catallax_topology::band_cap(model, boundary.edge, center, query, dec)?
            else {
                return Ok(None);
            };
            inside &= above == boundary.forward;
        }
        Ok((!boundaries.is_empty()).then_some(if inside {
            Location::Inside
        } else {
            Location::Outside
        }))
    }

    pub fn sphere_contains(
        &mut self,
        model: &Model,
        face: FaceId,
        point: PointId,
        dec: &mut Decider,
    ) -> Result<Option<Location>, Refusal> {
        let Surface::Sphere(sphere) = self.store.surface(model.face(face).surface) else {
            unreachable!("sphere membership on another carrier");
        };
        let center = sphere.center.clone();
        let axis = sphere.axis.clone();
        let reference = sphere.ref_dir.clone();
        let radius = sphere.radius.clone();
        let boundaries = self.boundaries(model, face, dec)?;
        if boundaries.is_empty() {
            return Ok(Some(Location::Inside));
        }
        let query = self.store.position(point);
        let loop_count = model
            .face(face)
            .loops
            .iter()
            .filter(|id| matches!(model.loop_(**id), Loop::Coedges(_)))
            .count();
        // A full circular loop separates the sphere by its exact cutting plane.
        // Its natural circle direction puts the positive halfspace on the left.
        if loop_count == boundaries.len() {
            let mut circular = true;
            let mut boundary = false;
            for edge in boundaries.iter() {
                let whole = match (edge.from, edge.to) {
                    (None, None) => true,
                    (Some(a), Some(b)) => self.same_point(a, b, dec)?,
                    _ => false,
                };
                if !whole {
                    circular = false;
                    break;
                }
                let form = self.form(edge.curve, dec)?;
                let CurveForm::Circle { center, normal, .. } = form.as_ref() else {
                    circular = false;
                    break;
                };
                match sign(&query.minus(center).dot(normal), dec)? {
                    Sign::Zero => boundary = true,
                    Sign::Positive if !edge.forward => return Ok(Some(Location::Outside)),
                    Sign::Negative if edge.forward => return Ok(Some(Location::Outside)),
                    _ => {}
                }
            }
            if circular {
                return Ok(Some(if boundary {
                    Location::Boundary
                } else {
                    Location::Inside
                }));
            }
        }
        if let Some(location) = self.band_caps(model, &boundaries, &query, &center, dec)? {
            return Ok(Some(location));
        }
        for edge in boundaries.iter() {
            match catallax_intersect::on_curve(&self.store, edge.curve, point, dec) {
                Outcome::ProvenCoincident(_) => {
                    if self.on_arc(edge, point, dec)? != Location::Outside {
                        return Ok(Some(Location::Boundary));
                    }
                }
                Outcome::Proven(_) => {}
                Outcome::NotExact(n) => return Err(Refusal::from(n)),
                Outcome::Refused(r) => return Err(r),
            }
        }
        let query_vector = query.minus(&center);
        let across = axis.cross(&reference);
        // Admission is optional: an uncertified point keeps the legacy export.
        // Bound the blind fallback so that failure stays cheap.
        let attempts = boundaries.len().saturating_add(BLIND_ATTEMPTS);
        // Aim rational planes toward known trim chords first. Blind carrier
        // rays can all miss a small spherical patch. Enclose only to choose
        // the candidate: arc membership and crossing orientation stay exact.
        // Keeping algebraic vertices out of the plane coefficients also avoids
        // raising the intersection eliminant's coefficient field.
        let epsilon = dec.enclose(&radius, &Rational::one())?.hi().abs()
            * Rational::from_ratio(1, 100_000_000).unwrap();
        let mut approximate = BTreeMap::<PointId, Vector3>::new();
        for attempt in 0..attempts {
            dec.check_deadline()?;
            let (ray, reverse) = if attempt < boundaries.len() {
                let (Some(from), Some(to)) = (boundaries[attempt].from, boundaries[attempt].to)
                else {
                    continue;
                };
                if from == to {
                    continue;
                }
                for point in [from, to] {
                    if let std::collections::btree_map::Entry::Vacant(entry) =
                        approximate.entry(point)
                    {
                        let p = self.store.position(point).minus(&center);
                        let coordinate =
                            |value: &Real, dec: &mut Decider| -> Result<Real, Refusal> {
                                let enclosure = dec.enclose(value, &epsilon)?;
                                Ok(number::rational(&enclosure.mid(), 8).value.into())
                            };
                        entry.insert(Vector3::new(
                            coordinate(&p.x, dec)?,
                            coordinate(&p.y, dec)?,
                            coordinate(&p.z, dec)?,
                        ));
                    }
                }
                (approximate[&from].plus(&approximate[&to]), false)
            } else {
                let serial = attempt - boundaries.len();
                let index = serial / 2;
                let magnitude = index.div_ceil(2) as i64;
                let a = Real::from_i64(if index % 2 == 0 {
                    -magnitude
                } else {
                    magnitude
                });
                (
                    reference
                        .plus(&across.scaled(&a))
                        .plus(&axis.scaled(&(&a * &a))),
                    serial % 2 != 0,
                )
            };
            let mut normal = query_vector.cross(&ray);
            if reverse {
                normal = normal.scaled(&Real::from_i64(-1));
            }
            if sign(&normal.norm_squared(), dec)? == Sign::Zero {
                continue;
            }
            let plane = self.plane(center.clone(), normal.clone(), query_vector.clone(), dec)?;
            let mut nearest: Option<(PointId, usize, Vector3)> = None;
            let mut regular = true;
            for (index, boundary) in boundaries.iter().enumerate() {
                let hits = self.intersections(boundary.curve, plane, dec)?;
                match hits.as_ref() {
                    CurveSurface::Contained => {
                        regular = false;
                        break;
                    }
                    CurveSurface::Disjoint => {}
                    CurveSurface::Meets(hits) => {
                        for hit in hits {
                            let location = self.on_arc(boundary, hit.point, dec)?;
                            if location == Location::Outside {
                                continue;
                            }
                            let p = self.store.position(hit.point).minus(&center);
                            if sign(&normal.dot(&query_vector.cross(&p)), dec)? != Sign::Positive {
                                continue;
                            }
                            if location == Location::Boundary {
                                regular = false;
                                break;
                            }
                            if hit.contact != Contact::Crossing {
                                continue;
                            }
                            if let Some((_, _, previous)) = &nearest {
                                match sign(&normal.dot(&p.cross(previous)), dec)? {
                                    Sign::Negative => continue,
                                    Sign::Zero => {
                                        regular = false;
                                        break;
                                    }
                                    Sign::Positive => {}
                                }
                            }
                            nearest = Some((hit.point, index, p));
                        }
                    }
                }
                if !regular {
                    break;
                }
            }
            if !regular {
                continue;
            }
            if let Some((target, index, position)) = nearest {
                let tangent = catallax_intersect::tangent_at(
                    &self.store,
                    boundaries[index].curve,
                    target,
                    dec,
                )?;
                let side = sign(&position.cross(&tangent[1]).dot(&query_vector), dec)?;
                if side == Sign::Zero {
                    continue;
                }
                return Ok(Some(
                    if (side == Sign::Positive) == boundaries[index].forward {
                        Location::Inside
                    } else {
                        Location::Outside
                    },
                ));
            }
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::BudgetPolicy;
    use catallax_shapes::RPoint3;

    #[test]
    fn algebraic_finite_spline_cuts_keep_endpoints_and_contacts() {
        let mut dec = Decider::default();
        let mut geometry = Geometry::new(&Model::new(), true);
        let root = Real::from_i64(2).sqrt();
        let spline = BSplineCurve::new(
            2,
            vec![0.into(), 0.into(), 0.into(), 1.into(), 1.into(), 1.into()],
            vec![
                Point3::new(-root.clone(), root.clone(), Real::zero()),
                Point3::new(Real::zero(), -root.clone(), Real::zero()),
                Point3::new(root.clone(), root, Real::zero()),
            ],
            vec![Real::one(); 3],
            &mut dec,
        )
        .unwrap();
        let curve = geometry.store.add_curve(Curve::BSpline(spline));
        for (normal, origin, expected) in [
            (
                Vector3::new(1.into(), 0.into(), 0.into()),
                Point3::origin(),
                Contact::Crossing,
            ),
            (
                Vector3::new(0.into(), 1.into(), 0.into()),
                Point3::origin(),
                Contact::Tangent,
            ),
        ] {
            let plane = geometry
                .plane(
                    origin,
                    normal,
                    Vector3::new(0.into(), 0.into(), 1.into()),
                    &mut dec,
                )
                .unwrap();
            let hits = geometry.intersections(curve, plane, &mut dec).unwrap();
            let CurveSurface::Meets(hits) = hits.as_ref() else {
                panic!("finite cut");
            };
            assert_eq!(hits.len(), 1);
            assert_eq!(hits[0].contact, expected);
        }
        let spline = match geometry.store.curve(curve) {
            Curve::BSpline(s) => s.clone(),
            _ => unreachable!(),
        };
        let plane = geometry
            .plane(
                spline.control(0).clone(),
                Vector3::new(1.into(), 0.into(), 0.into()),
                Vector3::new(0.into(), 1.into(), 0.into()),
                &mut dec,
            )
            .unwrap();
        let hits = geometry.intersections(curve, plane, &mut dec).unwrap();
        let CurveSurface::Meets(hits) = hits.as_ref() else {
            panic!("endpoint cut");
        };
        assert_eq!(hits.len(), 1);
        assert!(
            sign(
                &geometry
                    .store
                    .position(hits[0].point)
                    .minus(spline.control(0))
                    .norm_squared(),
                &mut dec
            )
            .unwrap()
                == Sign::Zero
        );
    }

    // The sphere-cylinder fixture mirrored in y = 0 is the same solid with
    // left-handed band frames. Its sphere face still keeps the equator and
    // still loses both poles.
    #[test]
    fn a_mirrored_sphere_band_keeps_its_caps() {
        std::thread::Builder::new()
            .stack_size(64 * 1024 * 1024)
            .spawn(|| {
                let mut value = include_str!("../../../parts/fillets/cases.jsonl")
                    .lines()
                    .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
                    .find(|case| case["id"] == "sphere-cylinder")
                    .unwrap()["recipe"]
                    .clone();
                value["steps"]
                    .as_array_mut()
                    .unwrap()
                    .push(serde_json::json!({
                        "id": "flipped", "op": "mirror", "body": "round",
                        "plane": {"point": ["0", "0", "0"], "normal": ["0", "1", "0"]},
                        "keep_original": false
                    }));
                value["result"] = "flipped".into();
                let recipe = catallax_recipe::Recipe::parse(&value.to_string()).unwrap();
                let mut dec = Decider::new(BudgetPolicy {
                    max_work: 200_000_000,
                    ..Default::default()
                });
                let built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
                let model = &built.model;
                let face = model
                    .faces_of_body(built.result)
                    .into_iter()
                    .find(|&face| {
                        matches!(
                            model.geom().surface(model.face(face).surface),
                            Surface::Sphere(_)
                        )
                    })
                    .unwrap();
                let mut geometry = Geometry::new(model, true);
                for ([x, y, z], expected) in [
                    ([5, 0, 0], Location::Inside),
                    ([-5, 0, 0], Location::Inside),
                    ([0, -5, 0], Location::Inside),
                    ([0, 0, 5], Location::Outside),
                    ([0, 0, -5], Location::Outside),
                ] {
                    let point = geometry.point(RPoint3::new(x.into(), y.into(), z.into()).into());
                    let found = geometry
                        .sphere_contains(model, face, point, &mut dec)
                        .unwrap();
                    assert!(found == Some(expected), "at {:?}", [x, y, z]);
                }
            })
            .unwrap()
            .join()
            .unwrap();
    }
}

/// A starved decision must reach the caller with its code and bound intact,
/// whatever code this wrapper reports. See docs/NOT-EXACT-INVENTORY.md.
#[cfg(test)]
mod undecided_evidence {
    use super::*;
    use catallax_numbers::{BudgetPolicy, Decider, NotExact, Outcome, Real};

    fn starved() -> Decider {
        Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 0,
        })
    }

    fn undecided(value: &Real) -> NotExact {
        match starved().sign(value) {
            Outcome::NotExact(decision) => decision,
            other => panic!("an empty allowance must stop undecided: {other:?}"),
        }
    }

    #[test]
    fn step_sign_keeps_the_decision() {
        let value = Real::from(2).sqrt();
        let decision = undecided(&value);
        let refusal = sign(&value, &mut starved()).unwrap_err();
        assert_eq!(refusal.code, decision.code);
        assert_eq!(refusal.undecided(), Some(&decision));
    }
}
