//! Analytic extrema candidates with exact trim membership. Traced branch extrema
//! use a budgeted second-derivative enclosure, never display geometry.

use crate::{
    Accuracy, Enclosure,
    integrate::{Stop, Work},
};
use catallax_intersect::{
    CurveForm, Cyclic, TracePiece, curve_form, cyclic_order, order_on_curve, trace,
};
use catallax_numbers::{
    Answer, Answered, BudgetPolicy, Decider, MarkLog, Outcome, Rational, Real, Refusal,
};
use catallax_shapes::{Point3, PointDef, Surface, Vector3};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};
use std::collections::{BTreeMap, BTreeSet};

/// An inclusive axis-aligned box containing the exact geometry, in model coordinates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundingBox {
    /// Inclusive lower coordinate bounds.
    pub min: [Rational; 3],
    /// Inclusive upper coordinate bounds.
    pub max: [Rational; 3],
}

#[derive(Clone)]
struct Extrema {
    min: [Enclosure; 3],
    max: [Enclosure; 3],
}
impl Extrema {
    fn merge(&mut self, other: &Self) {
        for i in 0..3 {
            self.min[i].lo = self.min[i].lo.clone().min(other.min[i].lo.clone());
            self.min[i].hi = self.min[i].hi.clone().min(other.min[i].hi.clone());
            self.max[i].lo = self.max[i].lo.clone().max(other.max[i].lo.clone());
            self.max[i].hi = self.max[i].hi.clone().max(other.max[i].hi.clone());
        }
    }
    fn reached(&self, accuracy: &Accuracy) -> bool {
        self.min
            .iter()
            .chain(&self.max)
            .all(|v| accuracy.reached(v))
    }
    fn finish(self, accuracy: &Accuracy) -> Answer<BoundingBox> {
        let mut marks = MarkLog::new();
        accuracy.mark(&mut marks, &self.min.each_ref(), "bounding minimum");
        accuracy.mark(&mut marks, &self.max.each_ref(), "bounding maximum");
        Answer::with_marks(
            BoundingBox {
                min: self.min.map(|v| v.lo),
                max: self.max.map(|v| v.hi),
            },
            marks,
        )
    }
}
fn merge(total: &mut Option<Extrema>, values: Extrema) {
    if let Some(total) = total {
        total.merge(&values);
    } else {
        *total = Some(values);
    }
}
fn budget() -> Refusal {
    Refusal::new(
        "measure.budget-exhausted",
        "no extrema enclosure within the measurement budget",
    )
}
fn stop(error: Stop) -> Refusal {
    match error {
        Stop::Refusal(r) => r,
        _ => budget(),
    }
}

struct Budget {
    policy: BudgetPolicy,
    bits: u32,
    remaining: u64,
}
impl Budget {
    fn new(dec: &Decider) -> Result<Self, Refusal> {
        let policy = *dec.policy();
        let bits = policy.max_precision_bits.min(256);
        if bits == 0 || policy.max_work == 0 {
            return Err(budget());
        }
        Ok(Self {
            policy,
            bits,
            remaining: policy.max_work,
        })
    }
    fn work<T>(
        &mut self,
        dec: &Decider,
        f: impl FnOnce(&mut Work) -> Result<T, Stop>,
    ) -> Result<T, Refusal> {
        let mut work = Work::new(&self.policy, self.bits, self.remaining, dec).map_err(stop)?;
        let result = f(&mut work);
        self.remaining = work.remaining();
        result.map_err(stop)
    }
    fn point(&mut self, point: &Point3, dec: &Decider) -> Result<Extrema, Refusal> {
        self.work(dec, |work| {
            let values = [&point.x, &point.y, &point.z]
                .map(|x| work.real(x).map(|b| Enclosure::from_ball(&b)));
            let [x, y, z] = values;
            let min = [x?, y?, z?];
            Ok(Extrema {
                max: min.clone(),
                min,
            })
        })
    }
}

struct Query<'a> {
    model: &'a Model,
    scratch: Option<Model>,
    budget: Budget,
    edges: BTreeMap<EdgeId, Extrema>,
}
impl<'a> Query<'a> {
    fn new(model: &'a Model, dec: &Decider) -> Result<Self, Refusal> {
        Ok(Self {
            model,
            scratch: None,
            budget: Budget::new(dec)?,
            edges: BTreeMap::new(),
        })
    }
    fn face(
        &mut self,
        id: FaceId,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Extrema, Refusal> {
        let surface = self.model.geom().surface(self.model.face(id).surface);
        if self.model.face(id).loops.is_empty() {
            let (center, radii) = match surface {
                Surface::Sphere(s) => (
                    &s.center,
                    std::array::from_fn::<_, 3, _>(|_| s.radius.clone()),
                ),
                Surface::Torus(t) => {
                    let axis = [&t.axis.x, &t.axis.y, &t.axis.z];
                    let squared = t.axis.norm_squared();
                    (
                        &t.center,
                        std::array::from_fn(|i| {
                            &(&t.major * &(Real::from(1) - &axis[i].square() / &squared).sqrt())
                                + &t.minor
                        }),
                    )
                }
                _ => {
                    return Err(Refusal::new(
                        "measure.unbounded-face",
                        "noncompact carrier has no finite trim",
                    ));
                }
            };
            return self.budget.work(dec, |work| {
                let coordinates = [&center.x, &center.y, &center.z];
                let mut min = std::array::from_fn(|_| Enclosure::exact(Rational::zero()));
                let mut max = min.clone();
                for i in 0..3 {
                    min[i] = Enclosure::from_ball(&work.real(&(coordinates[i] - &radii[i]))?);
                    max[i] = Enclosure::from_ball(&work.real(&(coordinates[i] + &radii[i]))?);
                }
                Ok(Extrema { min, max })
            });
        }
        let mut total = None;
        let mut edges = BTreeSet::new();
        for lp in &self.model.face(id).loops {
            match self.model.loop_(*lp) {
                Loop::Vertex(v) => merge(
                    &mut total,
                    self.budget.point(
                        &self.model.geom().position(self.model.vertex(*v).point),
                        dec,
                    )?,
                ),
                Loop::Coedges(cs) => {
                    for c in cs {
                        edges.insert(self.model.coedge(*c).edge);
                    }
                }
            }
        }
        for edge in edges {
            if !self.edges.contains_key(&edge) {
                let result = self.edge(edge, accuracy, dec)?;
                self.edges.insert(edge, result);
            }
            if let Some(total) = &mut total {
                total.merge(&self.edges[&edge]);
            } else {
                total = Some(self.edges[&edge].clone());
            }
        }
        for candidate in interior_candidates(surface, dec)? {
            let scratch = self.scratch.get_or_insert_with(|| self.model.clone());
            if catallax_boolean::face_contains(scratch, id, &candidate, dec)? {
                merge(&mut total, self.budget.point(&candidate, dec)?);
            }
        }
        total.ok_or_else(|| {
            Refusal::new(
                "measure.unbounded-face",
                "face has no finite extrema candidates",
            )
        })
    }
    fn edge(
        &mut self,
        id: EdgeId,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Extrema, Refusal> {
        let edge = self.model.edge(id);
        let (mut from, mut to) = (
            edge.start.map(|v| self.model.vertex(v).point),
            edge.end.map(|v| self.model.vertex(v).point),
        );
        if edge.sense == Sense::Reversed {
            std::mem::swap(&mut from, &mut to);
        }
        let form = curve_form(self.model.geom(), edge.curve, dec)?;
        if matches!(form, CurveForm::Branch) {
            let traced = trace(self.model.geom(), edge.curve, from, to, dec)?;
            return traced_extrema(traced.value.pieces(), accuracy, &mut self.budget, dec);
        }
        let closed = matches!(form, CurveForm::Circle { .. } | CurveForm::Ellipse { .. });
        if !closed && (from.is_none() || to.is_none()) {
            return Err(Refusal::new(
                "measure.unbounded-face",
                "open edge lacks finite endpoints",
            ));
        }
        let mut total = None;
        for endpoint in [from, to].into_iter().flatten() {
            merge(
                &mut total,
                self.budget
                    .point(&self.model.geom().position(endpoint), dec)?,
            );
        }
        for candidate in curve_candidates(form, dec)? {
            let inside = if closed && from == to {
                true
            } else {
                let (Some(from), Some(to)) = (from, to) else {
                    return Err(Refusal::new(
                        "measure.unbounded-face",
                        "open edge lacks finite endpoints",
                    ));
                };
                let scratch = self.scratch.get_or_insert_with(|| self.model.clone());
                let point = scratch
                    .geom_mut()
                    .add_point(PointDef::Explicit(candidate.clone()));
                if closed {
                    match cyclic_order(scratch.geom(), edge.curve, [from, point, to], dec) {
                        Outcome::Proven(Cyclic::Forward) | Outcome::ProvenCoincident(_) => true,
                        Outcome::Proven(Cyclic::Backward) => false,
                        Outcome::Refused(r) => return Err(r),
                        _ => return Err(unsettled()),
                    }
                } else {
                    let first = order(scratch.geom(), edge.curve, from, point, dec)?;
                    let last = order(scratch.geom(), edge.curve, point, to, dec)?;
                    first == 0 || last == 0 || first == last
                }
            };
            if inside {
                merge(&mut total, self.budget.point(&candidate, dec)?);
            }
        }
        total.ok_or_else(|| {
            Refusal::new(
                "measure.unbounded-face",
                "edge has no finite extrema candidates",
            )
        })
    }
}
fn unsettled() -> Refusal {
    Refusal::new(
        "measure.bounds-not-separated",
        "an extrema membership predicate is not exact",
    )
}
fn order(
    geom: &catallax_shapes::GeomStore,
    curve: catallax_shapes::CurveId,
    a: catallax_shapes::PointId,
    b: catallax_shapes::PointId,
    dec: &mut Decider,
) -> Result<i8, Refusal> {
    match order_on_curve(geom, curve, a, b, dec) {
        Outcome::Proven(catallax_numbers::StrictOrder::Less) => Ok(-1),
        Outcome::Proven(catallax_numbers::StrictOrder::Greater) => Ok(1),
        Outcome::ProvenCoincident(_) => Ok(0),
        Outcome::Refused(r) => Err(r),
        _ => Err(unsettled()),
    }
}
fn nonzero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(_) => Ok(true),
        Outcome::ProvenCoincident(_) => Ok(false),
        Outcome::Refused(r) => Err(r),
        _ => Err(unsettled()),
    }
}
fn axis(i: usize) -> Vector3 {
    Vector3::new(
        Real::from((i == 0) as i64),
        Real::from((i == 1) as i64),
        Real::from((i == 2) as i64),
    )
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(Real::from(1) / v.norm_squared().sqrt()))
}
fn pair(out: &mut Vec<Point3>, center: &Point3, offset: &Vector3) {
    out.push(center.plus(offset));
    out.push(center.plus(&offset.scaled(&Real::from(-1))));
}
fn interior_candidates(surface: &Surface, dec: &mut Decider) -> Result<Vec<Point3>, Refusal> {
    let mut points = Vec::new();
    match surface {
        Surface::Sphere(s) => {
            for i in 0..3 {
                pair(&mut points, &s.center, &axis(i).scaled(&s.radius));
            }
        }
        Surface::Torus(t) => {
            let normal = unit(&t.axis);
            let n = [&normal.x, &normal.y, &normal.z];
            for (i, component) in n.into_iter().enumerate() {
                let perpendicular = axis(i).minus(&normal.scaled(component));
                let radial = if nonzero(&perpendicular.norm_squared(), dec)? {
                    unit(&perpendicular)
                } else {
                    unit(&t.ref_dir)
                };
                let major = radial.scaled(&t.major);
                let minor = axis(i).scaled(&t.minor);
                pair(&mut points, &t.center.plus(&major), &minor);
                pair(
                    &mut points,
                    &t.center.plus(&major.scaled(&Real::from(-1))),
                    &minor,
                );
            }
        }
        Surface::Plane(_)
        | Surface::Cylinder(_)
        | Surface::Cone(_)
        | Surface::EllipticCylinder(_) => {}
        _ => {
            return Err(Refusal::new(
                "measure.unsupported-surface",
                "surface has no extrema definition",
            ));
        }
    }
    Ok(points)
}
fn curve_candidates(form: CurveForm, dec: &mut Decider) -> Result<Vec<Point3>, Refusal> {
    let mut points = Vec::new();
    match form {
        CurveForm::Line { .. } => {}
        CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ..
        } => {
            let n = [&normal.x, &normal.y, &normal.z];
            let norm2 = normal.norm_squared();
            for (i, component) in n.into_iter().enumerate() {
                let direction = axis(i).minus(&normal.scaled(&(component / &norm2)));
                let squared = direction.norm_squared();
                if nonzero(&squared, dec)? {
                    pair(
                        &mut points,
                        &center,
                        &direction.scaled(&(&radius_sq / &squared).sqrt()),
                    );
                }
            }
        }
        CurveForm::Ellipse {
            center,
            normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            let first = major_dir.scaled(&(&major_sq / &major_dir.norm_squared()).sqrt());
            let second = normal.cross(&major_dir);
            let second = second.scaled(&(&minor_sq / &second.norm_squared()).sqrt());
            for (a, b) in [&first.x, &first.y, &first.z]
                .into_iter()
                .zip([&second.x, &second.y, &second.z])
            {
                let squared = &a.square() + &b.square();
                if nonzero(&squared, dec)? {
                    pair(
                        &mut points,
                        &center,
                        &first
                            .scaled(a)
                            .plus(&second.scaled(b))
                            .scaled(&(Real::from(1) / squared.sqrt())),
                    );
                }
            }
        }
        CurveForm::Parabola {
            vertex,
            normal,
            axis_dir,
            focal,
        } => {
            let first = unit(&axis_dir);
            let second = unit(&normal.cross(&axis_dir));
            for (a, b) in [&first.x, &first.y, &first.z]
                .into_iter()
                .zip([&second.x, &second.y, &second.z])
            {
                if nonzero(a, dec)? {
                    let t = &(&Real::from(-2) * &focal) * &(b / a);
                    let x = &t.square() / &(&Real::from(4) * &focal);
                    points.push(vertex.plus(&first.scaled(&x).plus(&second.scaled(&t))));
                }
            }
        }
        CurveForm::Hyperbola {
            center,
            normal,
            axis_dir,
            transverse_sq,
            conjugate_sq,
            branch,
        } => {
            let first = unit(&axis_dir);
            let second = unit(&normal.cross(&axis_dir));
            let a = &transverse_sq.sqrt()
                * &Real::from(if branch == catallax_numbers::StrictSign::Positive {
                    1
                } else {
                    -1
                });
            let b = conjugate_sq.sqrt();
            for (x, y) in [&first.x, &first.y, &first.z]
                .into_iter()
                .zip([&second.x, &second.y, &second.z])
            {
                let denominator = &a * x;
                if !nonzero(&denominator, dec)? {
                    continue;
                }
                let k = -(&(&b * y) / &denominator);
                let separation = Real::from(1) - k.square();
                match dec.sign(&separation) {
                    Outcome::Proven(catallax_numbers::StrictSign::Positive) => {
                        let scale = Real::from(1) / separation.sqrt();
                        points.push(
                            center.plus(
                                &first
                                    .scaled(&(&a * &scale))
                                    .plus(&second.scaled(&(&(&b * &k) * &scale))),
                            ),
                        );
                    }
                    Outcome::Proven(_) | Outcome::ProvenCoincident(_) => {}
                    Outcome::Refused(r) => return Err(r),
                    _ => return Err(unsettled()),
                }
            }
        }
        _ => {
            return Err(Refusal::new(
                "measure.unsupported-curve",
                "curve has no finite extrema definition",
            ));
        }
    }
    Ok(points)
}

struct Panel {
    piece: usize,
    lo: Rational,
    hi: Rational,
    extrema: Extrema,
}
fn panel(
    pieces: &[TracePiece],
    piece: usize,
    lo: Rational,
    hi: Rational,
    budget: &mut Budget,
    dec: &Decider,
) -> Result<Panel, Refusal> {
    let (a, b, range) = budget.work(dec, |work| {
        Ok((
            work.rational(&lo)?,
            work.rational(&hi)?,
            work.interval(&lo, &hi)?,
        ))
    })?;
    let a = pieces[piece].evaluate(&a, dec)?;
    let b = pieces[piece].evaluate(&b, dec)?;
    let second = pieces[piece].taylor(&range, 2, dec)?.pop().unwrap();
    let width = &hi - &lo;
    let factor = (&width * &width).mul_pow2(-2);
    let error: [Rational; 3] = second
        .each_ref()
        .map(|v| &v.lo().abs().max(v.hi().abs()) * &factor);
    let min = std::array::from_fn(|i| Enclosure {
        lo: a[i].lo().min(b[i].lo()) - error[i].clone(),
        hi: a[i].hi().min(b[i].hi()),
    });
    let max = std::array::from_fn(|i| Enclosure {
        lo: a[i].lo().max(b[i].lo()),
        hi: a[i].hi().max(b[i].hi()) + error[i].clone(),
    });
    Ok(Panel {
        piece,
        lo,
        hi,
        extrema: Extrema { min, max },
    })
}
fn traced_extrema(
    pieces: &[TracePiece],
    accuracy: &Accuracy,
    budget: &mut Budget,
    dec: &Decider,
) -> Result<Extrema, Refusal> {
    let mut panels = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        let (lo, hi) = piece.range();
        panels.push(panel(pieces, i, lo, hi, budget, dec)?);
    }
    loop {
        dec.check_deadline()?;
        let mut total = panels[0].extrema.clone();
        for p in &panels[1..] {
            total.merge(&p.extrema);
        }
        if total.reached(accuracy) {
            return Ok(total);
        }
        let mut refined = Vec::new();
        for p in panels {
            let relevant = (0..3).any(|i| {
                (!accuracy.reached(&total.min[i]) && p.extrema.min[i].lo <= total.min[i].hi)
                    || (!accuracy.reached(&total.max[i]) && p.extrema.max[i].hi >= total.max[i].lo)
            });
            if relevant {
                let mid = (&p.lo + &p.hi).mul_pow2(-1);
                refined.push(panel(pieces, p.piece, p.lo, mid.clone(), budget, dec)?);
                refined.push(panel(pieces, p.piece, mid, p.hi, budget, dec)?);
            } else {
                refined.push(p);
            }
        }
        panels = refined;
    }
}

pub(crate) fn body(
    model: &Model,
    body: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<Option<BoundingBox>> {
    accuracy.validate();
    let faces = model.faces_of_body(body);
    if faces.is_empty() {
        return Ok(Answer::exact(None));
    }
    let mut query = Query::new(model, dec)?;
    let mut total = None;
    for face in faces {
        merge(&mut total, query.face(face, accuracy, dec)?);
    }
    let answer = total.unwrap().finish(accuracy);
    Ok(Answer::with_marks(Some(answer.value), answer.marks))
}
pub(crate) fn face(
    model: &Model,
    face: FaceId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<BoundingBox> {
    accuracy.validate();
    let mut query = Query::new(model, dec)?;
    Ok(query.face(face, accuracy, dec)?.finish(accuracy))
}
