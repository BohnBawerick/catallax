//! Exact trim queries used by export-only periodic charts. The scratch store
//! shares the model's immutable geometry; no model topology or geometry changes.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::rc::Rc;

use catallax_intersect::{Contact, CurveForm, CurveSurface, Cyclic};
use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign, StrictOrder, StrictSign};
use catallax_shapes::{
    CurveId, GeomStore, Plane, Point3, PointDef, PointId, Surface, SurfaceId, Vector3,
};
use catallax_topology::{EdgeId, FaceId, Loop, Model, Sense};

use crate::number;

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
    hits: BTreeMap<(CurveId, SurfaceId), Rc<CurveSurface>>,
    forms: BTreeMap<CurveId, Rc<CurveForm>>,
    closed: BTreeMap<CurveId, bool>,
    boundaries: BTreeMap<FaceId, Rc<Vec<Boundary>>>,
}

pub(crate) fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
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

impl Geometry {
    pub fn new(model: &Model) -> Self {
        Self {
            store: model.geom().clone(),
            planes: Vec::new(),
            hits: BTreeMap::new(),
            forms: BTreeMap::new(),
            closed: BTreeMap::new(),
            boundaries: BTreeMap::new(),
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
        let result = Rc::new(
            catallax_intersect::intersect_curve_surface(&mut self.store, curve, plane, dec)?.value,
        );
        self.hits.insert((curve, plane), result.clone());
        Ok(result)
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
            Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
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
                Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
                Outcome::Refused(r) => Err(r),
            }
        } else {
            for (a, b) in [(from, point), (point, to)] {
                match catallax_intersect::order_on_curve(&self.store, edge.curve, a, b, dec) {
                    Outcome::Proven(StrictOrder::Less) => {}
                    Outcome::Proven(StrictOrder::Greater) => return Ok(Location::Outside),
                    Outcome::ProvenCoincident(_) => return Ok(Location::Boundary),
                    Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
                    Outcome::Refused(r) => return Err(r),
                }
            }
            Ok(Location::Inside)
        }
    }

    pub fn sphere_contains(
        &mut self,
        model: &Model,
        face: FaceId,
        point: PointId,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let Surface::Sphere(sphere) = self.store.surface(model.face(face).surface) else {
            unreachable!("sphere membership on another carrier");
        };
        let center = sphere.center.clone();
        let axis = sphere.axis.clone();
        let reference = sphere.ref_dir.clone();
        let radius = sphere.radius.clone();
        let boundaries = self.boundaries(model, face, dec)?;
        if boundaries.is_empty() {
            return Ok(Location::Inside);
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
                    Sign::Positive if !edge.forward => return Ok(Location::Outside),
                    Sign::Negative if edge.forward => return Ok(Location::Outside),
                    _ => {}
                }
            }
            if circular {
                return Ok(if boundary {
                    Location::Boundary
                } else {
                    Location::Inside
                });
            }
        }
        for edge in boundaries.iter() {
            match catallax_intersect::on_curve(&self.store, edge.curve, point, dec) {
                Outcome::ProvenCoincident(_) => {
                    if self.on_arc(edge, point, dec)? != Location::Outside {
                        return Ok(Location::Boundary);
                    }
                }
                Outcome::Proven(_) => {}
                Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
                Outcome::Refused(r) => return Err(r),
            }
        }
        let query_vector = query.minus(&center);
        let across = axis.cross(&reference);
        let attempts = boundaries.len().saturating_add(1).saturating_mul(128);
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
                return Ok(if (side == Sign::Positive) == boundaries[index].forward {
                    Location::Inside
                } else {
                    Location::Outside
                });
            }
        }
        Err(number::budget(
            "no transverse spherical trim membership certificate",
        ))
    }
}
