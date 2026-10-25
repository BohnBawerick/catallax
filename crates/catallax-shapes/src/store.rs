//! Handles, write-once arenas and the geometry store.

use im::Vector;
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

use catallax_numbers::{
    Ball, Coincidence, Decider, Fact, FactBook, FactKind, FactOrigin, NotExact, Outcome, Rational,
    Real, RealPoly, Refusal, Sign, StrictSign, ZeroProof,
};

use crate::math;
use crate::{Curve, MeetOf, Point3, PointDef, Position, Surface, Vector3};

/// A typed index into an `Arena`. Meaningful only for the arena that issued it.
pub struct Handle<T> {
    index: u32,
    _kind: PhantomData<fn() -> T>,
}

impl<T> Handle<T> {
    /// Position in the arena.
    pub fn index(&self) -> u32 {
        self.index
    }
}

impl<T> Clone for Handle<T> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<T> Copy for Handle<T> {}

impl<T> PartialEq for Handle<T> {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}

impl<T> Eq for Handle<T> {}

impl<T> PartialOrd for Handle<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Ord for Handle<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.index.cmp(&other.index)
    }
}

impl<T> Hash for Handle<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
    }
}

impl<T> fmt::Debug for Handle<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "#{}", self.index)
    }
}

/// An append-only store. Clones share chunks; appending copies only the changed path.
#[derive(Clone, Debug)]
pub struct Arena<T: Clone> {
    values: Vector<T>,
}

impl<T: Clone> Default for Arena<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Clone> Arena<T> {
    /// An empty arena.
    pub fn new() -> Self {
        Self {
            values: Vector::new(),
        }
    }

    /// Adds an entity and returns its handle; handles are issued in order from 0.
    pub fn push(&mut self, value: T) -> Handle<T> {
        let index = u32::try_from(self.values.len()).expect("arena handle space exhausted");
        self.values.push_back(value);
        Handle {
            index,
            _kind: PhantomData,
        }
    }

    /// The entity behind a handle. Panics on a handle from another arena that is out of range.
    pub fn get(&self, handle: Handle<T>) -> &T {
        &self.values[handle.index as usize]
    }

    /// Number of entities.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether the arena is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// All handles in order.
    pub fn handles(&self) -> Vec<Handle<T>> {
        (0..self.values.len())
            .map(|i| Handle {
                index: u32::try_from(i).expect("valid arena index"),
                _kind: PhantomData,
            })
            .collect()
    }
}

/// Handle of a surface.
pub type SurfaceId = Handle<Surface>;
/// Handle of a curve.
pub type CurveId = Handle<Curve>;
/// Handle of a point definition.
pub type PointId = Handle<PointDef>;

/// Names a geometric object in the fact book.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum GeomKey {
    /// A surface.
    Surface(SurfaceId),
    /// A curve.
    Curve(CurveId),
    /// A point.
    Point(PointId),
}

/// Which side of a surface, relative to its normal.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum Side {
    /// The side the normal points to (implicit positive).
    Front,
    /// The other side (implicit negative).
    Back,
}

/// Proven not to coincide or not to be incident.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Apart;

/// The write-once store of surfaces, curves and points, with its fact book.
#[derive(Clone, Debug, Default)]
pub struct GeomStore {
    surfaces: Arena<Surface>,
    curves: Arena<Curve>,
    points: Arena<PointDef>,
    facts: FactBook<GeomKey>,
    implicits: Vector<Option<RealPoly>>,
    approximate_points: Vector<[f64; 3]>,
}

impl GeomStore {
    /// An empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a surface without searching unrelated geometry for relations.
    pub fn add_surface(&mut self, s: Surface) -> SurfaceId {
        self.implicits.push_back(s.implicit());
        self.surfaces.push(s)
    }

    /// Adds a curve; an intersection curve also records its incidence facts.
    pub fn add_curve(&mut self, c: Curve) -> CurveId {
        if let Curve::Intersection(i) = &c {
            assert_ne!(
                i.surfaces[0], i.surfaces[1],
                "intersection requires distinct surfaces"
            );
        }
        let id = self.curves.push(c);
        if let Curve::Intersection(i) = self.curves.get(id) {
            let surfaces = i.surfaces;
            for s in surfaces {
                self.incident(GeomKey::Curve(id), GeomKey::Surface(s));
            }
        }
        id
    }

    /// Adds a point; a meeting point also records its incidence facts.
    pub fn add_point(&mut self, p: PointDef) -> PointId {
        if let PointDef::Meet(m) = &p {
            if let MeetOf::Surfaces(surfaces) = &m.of {
                assert!(
                    surfaces.len() >= 2,
                    "meeting requires at least two surfaces"
                );
                for (i, s) in surfaces.iter().enumerate() {
                    assert!(
                        !surfaces[..i].contains(s),
                        "meeting requires distinct surfaces"
                    );
                }
            }
            if let Position::Isolated(root) = &m.position {
                assert_eq!(
                    root.root_box().dimension(),
                    3,
                    "point root requires three variables"
                );
            }
        }
        let id = self.points.push(p);
        let approximate = match self.points.get(id) {
            PointDef::Explicit(p) => [&p.x, &p.y, &p.z].map(Real::to_f64),
            PointDef::Meet(_) => {
                let p = self.position(id);
                [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()]
            }
        };
        self.approximate_points.push_back(approximate);
        if let PointDef::Meet(m) = self.points.get(id) {
            let facts = &mut self.facts;
            let mut record = |object| {
                facts.record(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Point(id),
                    b: object,
                    at: None,
                    origin: FactOrigin::Construction("shapes.definition"),
                });
            };
            match &m.of {
                MeetOf::Surfaces(s) => {
                    for surface in s {
                        record(GeomKey::Surface(*surface));
                    }
                }
                MeetOf::CurveSurface(c, s) => {
                    record(GeomKey::Curve(*c));
                    record(GeomKey::Surface(*s));
                }
                MeetOf::Curves(a, b) => {
                    record(GeomKey::Curve(*a));
                    record(GeomKey::Curve(*b));
                }
            }
        }
        id
    }

    /// The surface behind a handle.
    pub fn surface(&self, id: SurfaceId) -> &Surface {
        self.surfaces.get(id)
    }

    /// The curve behind a handle.
    pub fn curve(&self, id: CurveId) -> &Curve {
        self.curves.get(id)
    }

    /// The point definition behind a handle.
    pub fn point(&self, id: PointId) -> &PointDef {
        self.points.get(id)
    }

    /// Number of surfaces.
    pub fn surface_count(&self) -> usize {
        self.surfaces.len()
    }

    /// Number of curves.
    pub fn curve_count(&self) -> usize {
        self.curves.len()
    }

    /// Number of points.
    pub fn point_count(&self) -> usize {
        self.points.len()
    }

    /// The fact book.
    pub fn facts(&self) -> &FactBook<GeomKey> {
        &self.facts
    }

    /// Records a fact and circular-carrier relations it guarantees; `false`
    /// when the supplied fact is already present.
    pub fn record_fact(&mut self, fact: Fact<GeomKey>) -> bool {
        let recorded = self.facts.record(fact);
        if fact.kind == FactKind::Incident
            && fact.at.is_none()
            && let (GeomKey::Curve(curve), GeomKey::Surface(surface)) = (fact.a, fact.b)
            && let Curve::Circle(circle) = self.curve(curve)
        {
            let surface = self.surface(surface);
            // Circular sections of circular cylinders and cones share their
            // carrier's axis. Tori also admit non-axial circular sections.
            let coaxial = match surface {
                Surface::Cylinder(_) | Surface::Cone(_) => true,
                Surface::Torus(torus) => crate::intent::parallel(&circle.normal, &torus.axis),
                _ => false,
            };
            let concentric = crate::intent::center(surface)
                .is_some_and(|center| crate::intent::same_point(center, &circle.center));
            for (kind, holds) in [
                (FactKind::Coaxial, coaxial),
                (FactKind::Concentric, concentric),
            ] {
                if holds {
                    self.facts.record(Fact {
                        kind,
                        a: fact.a,
                        b: fact.b,
                        at: None,
                        origin: fact.origin,
                    });
                }
            }
        }
        recorded
    }

    /// Borrowed global analytic equation, absent for a parametric support.
    pub fn implicit(&self, s: SurfaceId) -> Option<&RealPoly> {
        self.implicits[s.index as usize].as_ref()
    }

    /// Exact position of a point.
    pub fn position(&self, p: PointId) -> Point3 {
        match self.point(p) {
            PointDef::Explicit(p) => p.clone(),
            PointDef::Meet(m) => match &m.position {
                Position::Closed(p) => p.clone(),
                Position::Isolated(root) => {
                    Point3::new(root.coordinate(0), root.coordinate(1), root.coordinate(2))
                }
            },
        }
    }

    /// Cached approximate coordinates, for display and search only.
    pub fn approx_point(&self, p: PointId) -> [f64; 3] {
        self.approximate_points[p.index as usize]
    }

    /// Proven enclosure of a point, as narrow as `max_width` if the budget allows.
    pub fn enclose_point(
        &self,
        p: PointId,
        max_width: &Rational,
        dec: &mut Decider,
    ) -> Result<[Ball; 3], Refusal> {
        let p = self.position(p);
        Ok([
            dec.enclose(&p.x, max_width)?,
            dec.enclose(&p.y, max_width)?,
            dec.enclose(&p.z, max_width)?,
        ])
    }

    /// Outward normal covector at a surface point. A parametric support must
    /// prove its finite chart preimage; it has no global implicit gradient.
    pub fn normal(&self, s: SurfaceId, p: PointId, dec: &mut Decider) -> Result<Vector3, Refusal> {
        self.surface(s).normal_at(&self.position(p), dec)
    }

    /// Which side of `s` the point `p` is on; coincident when on its carrier.
    /// Facts first; a not-exact outcome carries a distance bound.
    pub fn side(&self, s: SurfaceId, p: PointId, dec: &mut Decider) -> Outcome<Side> {
        if self
            .facts
            .holds(FactKind::Incident, GeomKey::Point(p), GeomKey::Surface(s))
        {
            return fact_decision(dec);
        }
        self.side_of(s, &self.position(p), dec)
    }

    /// Which side of `s` an exact point is on; coincident when on its carrier.
    /// A not-exact outcome carries a distance bound.
    pub fn side_of(&self, s: SurfaceId, p: &Point3, dec: &mut Decider) -> Outcome<Side> {
        let surface = self.surface(s);
        if let Surface::BSpline(spline) = surface {
            return match spline.side_at(p, dec) {
                Ok(Sign::Positive) => Outcome::Proven(Side::Front),
                Ok(Sign::Negative) => Outcome::Proven(Side::Back),
                Ok(Sign::Zero) => Outcome::ProvenCoincident(Coincidence {
                    proof: ZeroProof::Algebraic,
                }),
                Err(refusal) => Outcome::Refused(refusal),
            };
        }
        let Some(value) = surface.value(p) else {
            return Outcome::Refused(missing_equation());
        };
        match dec.sign(&value).map(|s| match s {
            StrictSign::Positive => Side::Front,
            StrictSign::Negative => Side::Back,
        }) {
            Outcome::NotExact(mut n) => match surface.distance_bound(&n.bound, dec) {
                Ok(bound) => {
                    n.bound = bound;
                    Outcome::NotExact(n)
                }
                Err(refusal) => Outcome::Refused(refusal),
            },
            other => other,
        }
    }

    /// Whether `p` lies on the common zero set of the curve's equations.
    /// Facts first; a not-exact outcome carries a distance bound.
    pub fn point_on_curve(&self, c: CurveId, p: PointId, dec: &mut Decider) -> Outcome<Apart> {
        let curve = self.curve(c);
        let known_supports = match curve {
            Curve::Intersection(c) => c.surfaces.map(|s| {
                self.facts
                    .holds(FactKind::Incident, GeomKey::Point(p), GeomKey::Surface(s))
            }),
            _ => [false; 2],
        };
        if known_supports.into_iter().all(|known| known)
            || self
                .facts
                .holds(FactKind::Incident, GeomKey::Point(p), GeomKey::Curve(c))
        {
            return fact_decision(dec);
        }
        let position = self.position(p);
        if let Curve::BSpline(curve) = curve {
            return match curve.parameter_of(&position, dec) {
                Ok(Some(_)) => Outcome::ProvenCoincident(Coincidence {
                    proof: ZeroProof::Algebraic,
                }),
                Ok(None) => Outcome::Proven(Apart),
                Err(refusal) => Outcome::Refused(refusal),
            };
        }
        let coordinates = [position.x, position.y, position.z];
        let Some(equations) = curve.equations(self) else {
            return Outcome::Refused(missing_equation());
        };
        // A retained incidence discharges only its own support equation. Keep
        // the remaining equations and their distance-bound indices unchanged.
        let mut indices = [0, 1, 2];
        let mut count = 0;
        for i in 0..equations.len() {
            if known_supports.get(i).copied().unwrap_or(false) {
                dec.note_fact_decision();
            } else {
                indices[count] = i;
                count += 1;
            }
        }
        let indices = &indices[..count];
        decide_zeroes(
            indices
                .iter()
                .map(|&i| equations[i].eval_real(&coordinates)),
            dec,
            |bounds, dec| {
                let mut original = std::array::from_fn::<_, 3, _>(|_| Rational::zero());
                for (&i, bound) in indices.iter().zip(bounds) {
                    original[i] = bound.clone();
                }
                let bounds = &original[..equations.len()];
                match curve {
                    Curve::Line(c) => Ok(math::sqrt_upper(&sum_squares(bounds))
                        * math::inverse_length_upper(&c.dir, dec)?),
                    Curve::Circle(c) => {
                        let h_squared =
                            &Real::from(&bounds[0] * &bounds[0]) / &c.normal.norm_squared();
                        let value =
                            &Real::from(bounds[1].clone()) + &(&Real::from_i64(2) * &h_squared);
                        if let Some(value) = value.as_rational() {
                            return Ok(math::sqrt_upper(&value));
                        }
                        Ok(dec.enclose_at(&value.sqrt(), 64)?.hi())
                    }
                    Curve::Ellipse(c) => {
                        let major = if let Some(value) = c.major_sq.as_rational() {
                            math::sqrt_upper(&value)
                        } else {
                            dec.enclose_at(&c.major_sq.sqrt(), 64)?.hi()
                        };
                        Ok(&major * &bounds[1]
                            + &bounds[0] * &math::inverse_length_upper(&c.normal, dec)?)
                    }
                    Curve::Intersection(c) => Ok(self
                        .surface(c.surfaces[0])
                        .distance_bound(&bounds[0], dec)?
                        .max(
                            self.surface(c.surfaces[1])
                                .distance_bound(&bounds[1], dec)?,
                        )),
                    Curve::BSpline(_) => Err(missing_equation()),
                }
            },
        )
    }

    /// Whether two points coincide. Facts first; a not-exact outcome carries
    /// a distance bound.
    pub fn same_point(&self, a: PointId, b: PointId, dec: &mut Decider) -> Outcome<Apart> {
        if self
            .facts
            .holds(FactKind::Coincident, GeomKey::Point(a), GeomKey::Point(b))
        {
            return fact_decision(dec);
        }
        let d = self.position(a).minus(&self.position(b));
        decide_zeroes([d.x, d.y, d.z], dec, |bounds, _| {
            Ok(math::sqrt_upper(&sum_squares(bounds)))
        })
    }

    /// Whether two surfaces define the same point set. Facts and rational
    /// arithmetic retain their fast paths; algebraic comparisons require proofs.
    pub fn same_surface(&self, a: SurfaceId, b: SurfaceId, dec: &mut Decider) -> Outcome<Apart> {
        if self.facts.holds(
            FactKind::Coincident,
            GeomKey::Surface(a),
            GeomKey::Surface(b),
        ) {
            return fact_decision(dec);
        }
        let known = |kind, dec: &mut Decider| {
            let holds = self
                .facts
                .holds(kind, GeomKey::Surface(a), GeomKey::Surface(b));
            if holds {
                dec.note_fact_decision();
            }
            holds
        };
        let zero =
            |value: Real, dec: &mut Decider| math::sign(&value, dec).map(|s| s == Sign::Zero);
        let equal = (|| -> Result<bool, Refusal> {
            Ok(match (self.surface(a), self.surface(b)) {
                (Surface::Plane(a), Surface::Plane(b)) => {
                    (known(FactKind::Parallel, dec)
                        || zero(a.normal.cross(&b.normal).norm_squared(), dec)?)
                        && zero(b.origin.minus(&a.origin).dot(&a.normal), dec)?
                }
                (Surface::Cylinder(a), Surface::Cylinder(b)) => {
                    zero(&a.radius_sq - &b.radius_sq, dec)?
                        && (known(FactKind::Coaxial, dec)
                            || (zero(a.axis.cross(&b.axis).norm_squared(), dec)?
                                && zero(
                                    b.origin.minus(&a.origin).cross(&a.axis).norm_squared(),
                                    dec,
                                )?))
                }
                (Surface::Cone(a), Surface::Cone(b)) => {
                    zero(a.apex.minus(&b.apex).norm_squared(), dec)?
                        && zero(&a.slope - &b.slope, dec)?
                        && (known(FactKind::Parallel, dec)
                            || zero(a.axis.cross(&b.axis).norm_squared(), dec)?)
                        && math::sign(&a.axis.dot(&b.axis), dec)? == Sign::Positive
                }
                (Surface::Sphere(a), Surface::Sphere(b)) => {
                    (known(FactKind::Concentric, dec)
                        || zero(a.center.minus(&b.center).norm_squared(), dec)?)
                        && zero(&a.radius - &b.radius, dec)?
                }
                (Surface::Torus(a), Surface::Torus(b)) => {
                    (known(FactKind::Concentric, dec)
                        || zero(a.center.minus(&b.center).norm_squared(), dec)?)
                        && zero(&a.major - &b.major, dec)?
                        && zero(&a.minor - &b.minor, dec)?
                        && (known(FactKind::Parallel, dec)
                            || zero(a.axis.cross(&b.axis).norm_squared(), dec)?)
                }
                (Surface::EllipticCylinder(_), Surface::EllipticCylinder(_)) => self
                    .implicit(a)
                    .ok_or_else(missing_equation)?
                    .equals(self.implicit(b).ok_or_else(missing_equation)?, dec)?,
                (Surface::BSpline(a), Surface::BSpline(b)) => {
                    dec.check_deadline()?;
                    if !a.carrier().same_definition(b.carrier()) {
                        return Err(Refusal::new(
                            "shapes.parametric-query",
                            "unrelated spline supports need a certified equality solver",
                        ));
                    }
                    true
                }
                (Surface::BSpline(_), _) | (_, Surface::BSpline(_)) => {
                    return Err(Refusal::new(
                        "shapes.parametric-query",
                        "cross-kind spline support equality is not covered",
                    ));
                }
                _ => false,
            })
        })();
        match equal {
            Ok(true) => Outcome::ProvenCoincident(Coincidence {
                proof: if self.surface(a).has_rational_definition()
                    && self.surface(b).has_rational_definition()
                {
                    ZeroProof::ExactArithmetic
                } else {
                    ZeroProof::Algebraic
                },
            }),
            Ok(false) => Outcome::Proven(Apart),
            Err(refusal) => Outcome::Refused(refusal),
        }
    }

    /// Sign of the dot product of two co-oriented surface normals at `at`.
    /// A not-exact outcome bounds the dot product itself, not a displacement.
    pub fn normals_agree(
        &self,
        a: SurfaceId,
        b: SurfaceId,
        at: PointId,
        dec: &mut Decider,
    ) -> Outcome<StrictSign> {
        let p = self.position(at);
        let first = match self.surface(a).normal_at(&p, dec) {
            Ok(normal) => normal,
            Err(refusal) => return Outcome::Refused(refusal),
        };
        let second = match self.surface(b).normal_at(&p, dec) {
            Ok(normal) => normal,
            Err(refusal) => return Outcome::Refused(refusal),
        };
        dec.sign(&first.dot(&second))
    }

    fn incident(&mut self, a: GeomKey, b: GeomKey) {
        self.facts.record(Fact {
            kind: FactKind::Incident,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction("shapes.definition"),
        });
    }
}

fn missing_equation() -> Refusal {
    Refusal::new(
        "shapes.parametric-query",
        "this query needs a finite parametric solver rather than global equations",
    )
}

fn fact_decision<T>(dec: &mut Decider) -> Outcome<T> {
    dec.note_fact_decision();
    Outcome::ProvenCoincident(Coincidence {
        proof: ZeroProof::ConstructionFact,
    })
}

fn sum_squares(bounds: &[Rational]) -> Rational {
    bounds.iter().fold(Rational::zero(), |s, b| s + b * b)
}

fn decide_zeroes(
    values: impl IntoIterator<Item = Real>,
    dec: &mut Decider,
    distance_bound: impl FnOnce(&[Rational], &Decider) -> Result<Rational, Refusal>,
) -> Outcome<Apart> {
    // Keep equation indices: curve-distance bounds use different units
    // for different equations, even when a later equation is tested first.
    let mut pending = [None, None, None];
    let mut count = 0;
    for (i, value) in values.into_iter().enumerate() {
        pending[i] = Some(value);
        count = i + 1;
    }
    let mut order = [0, 1, 2];
    if let Some(i) = (0..count).find(|&i| {
        let value = pending[i].as_ref().unwrap();
        value.approx().sign().is_some()
            || value.as_rational().is_some_and(|r| r.sign() != Sign::Zero)
    }) {
        order.swap(0, i);
    } else {
        for (i, value) in pending.iter().enumerate().take(count) {
            let value = value.as_ref().unwrap();
            if value.as_rational().is_some() {
                continue;
            }
            let Some(width) = Rational::from_f64_exact(value.approx().width()) else {
                continue;
            };
            if width.sign() != Sign::Positive {
                continue;
            }
            match dec.enclose(value, &width) {
                Ok(ball) if ball.sign().is_some() => {
                    order.swap(0, i);
                    break;
                }
                Err(r) if r.code == "numbers.deadline" => return Outcome::Refused(r),
                _ => {}
            }
        }
    }
    let mut bounds = std::array::from_fn::<_, 3, _>(|_| Rational::zero());
    let mut first_uncertain: Option<NotExact> = None;
    let mut coincidence = None;
    for i in order.into_iter().take(count) {
        let value = pending[i].take().unwrap();
        match dec.sign(&value) {
            Outcome::Proven(_) => return Outcome::Proven(Apart),
            Outcome::Refused(r) => return Outcome::Refused(r),
            Outcome::ProvenCoincident(c) => coincidence = Some(c),
            Outcome::NotExact(n) => {
                bounds[i] = n.bound.clone();
                if first_uncertain.is_none() {
                    first_uncertain = Some(n);
                }
            }
        }
    }
    if let Some(mut n) = first_uncertain {
        match distance_bound(&bounds[..count], dec) {
            Ok(bound) => {
                n.bound = bound;
                Outcome::NotExact(n)
            }
            Err(refusal) => Outcome::Refused(refusal),
        }
    } else {
        Outcome::ProvenCoincident(coincidence.expect("at least one equation"))
    }
}

#[cfg(test)]
mod contact_fact_tests {
    use super::*;
    use crate::{IntersectionCurve, Plane, Sphere};
    use catallax_numbers::BudgetPolicy;

    fn fixture() -> (GeomStore, SurfaceId, CurveId, PointId) {
        let mut store = GeomStore::new();
        let mut dec = Decider::default();
        let sphere = store.add_surface(Surface::Sphere(
            Sphere::new(
                Point3::origin(),
                Real::one(),
                Vector3::new(Real::zero(), Real::zero(), Real::one()),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap(),
        ));
        let plane = store.add_surface(Surface::Plane(
            Plane::new(
                Point3::origin(),
                Vector3::new(Real::zero(), Real::zero(), Real::one()),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap(),
        ));
        // 1/3 + 2/3 = 1 proves this point is on the unit sphere.
        let third = Real::one() / Real::from_i64(3);
        let p = Point3::new(
            third.clone().sqrt(),
            (third * Real::from_i64(2)).sqrt(),
            Real::zero(),
        );
        assert!(
            dec.sign(&store.surface(sphere).value(&p).unwrap())
                .is_coincident()
        );
        let point = store.add_point(PointDef::Explicit(p));
        let curve = store.add_curve(Curve::Intersection(IntersectionCurve {
            surfaces: [sphere, plane],
            witness: point,
        }));
        (store, sphere, curve, point)
    }
    fn small() -> Decider {
        Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 8,
        })
    }
    #[test]
    fn retained_partial_support_regression() {
        let (mut store, sphere, curve, point) = fixture();
        assert!(
            !store
                .point_on_curve(curve, point, &mut small())
                .is_coincident()
        );
        store.record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Surface(sphere),
            at: None,
            origin: FactOrigin::Proof(ZeroProof::Algebraic),
        });
        let facts = store.facts().len();
        let mut dec = small();
        assert!(store.point_on_curve(curve, point, &mut dec).is_coincident());
        assert_eq!(
            store.facts().len(),
            facts,
            "a read-only decision must not append facts"
        );
        assert!(dec.stats().coincident_fact > 0);
    }
    #[test]
    fn one_support_fact_does_not_discharge_the_other_support() {
        let (mut store, sphere, curve, _) = fixture();
        let point = store.add_point(PointDef::Explicit(Point3::new(
            Real::zero(),
            Real::zero(),
            Real::one(),
        )));
        store.record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Surface(sphere),
            at: None,
            origin: FactOrigin::Proof(ZeroProof::ExactArithmetic),
        });
        assert!(matches!(
            store.point_on_curve(curve, point, &mut small()),
            Outcome::Proven(Apart)
        ));
        let (mut store, sphere, curve, point) = fixture();
        store.record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Surface(sphere),
            at: None,
            origin: FactOrigin::Proof(ZeroProof::Algebraic),
        });
        let mut zero = Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 0,
        });
        assert!(
            !store
                .point_on_curve(curve, point, &mut zero)
                .is_coincident()
        );
    }
    #[test]
    fn partial_support_keeps_the_unknown_equations_distance_units() {
        let mut store = GeomStore::new();
        let mut dec = Decider::default();
        let z = Real::one() / Real::from_i64(3);
        let point = Point3::new((Real::from_i64(49) - z.square()).sqrt(), Real::zero(), z);
        let sphere = store.add_surface(Surface::Sphere(
            Sphere::new(
                Point3::origin(),
                Real::from_i64(7),
                Vector3::new(Real::zero(), Real::zero(), Real::one()),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap(),
        ));
        let plane = store.add_surface(Surface::Plane(
            Plane::new(
                Point3::origin(),
                Vector3::new(Real::zero(), Real::zero(), Real::one()),
                Vector3::new(Real::one(), Real::zero(), Real::zero()),
                &mut dec,
            )
            .unwrap(),
        ));
        assert!(
            dec.sign(&store.surface(sphere).value(&point).unwrap())
                .is_coincident()
        );
        let point = store.add_point(PointDef::Explicit(point));
        let witness = store.add_point(PointDef::Explicit(Point3::new(
            Real::from_i64(7),
            Real::zero(),
            Real::zero(),
        )));
        let curve = store.add_curve(Curve::Intersection(IntersectionCurve {
            surfaces: [sphere, plane],
            witness,
        }));
        store.record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Point(point),
            b: GeomKey::Surface(sphere),
            at: None,
            origin: FactOrigin::Proof(ZeroProof::Algebraic),
        });
        let policy = BudgetPolicy {
            max_precision_bits: 64,
            max_work: 0,
        };
        let value = store.surface(plane).value(&store.position(point)).unwrap();
        let Outcome::NotExact(expected) = Decider::new(policy).sign(&value) else {
            panic!("zero budget retains the unresolved plane value");
        };
        let expected = store
            .surface(plane)
            .distance_bound(&expected.bound, &Decider::new(policy))
            .unwrap();
        let Outcome::NotExact(actual) =
            store.point_on_curve(curve, point, &mut Decider::new(policy))
        else {
            panic!("the unknown support remains unresolved");
        };
        assert_eq!(actual.bound, expected);
        assert!(actual.bound >= Rational::one() / Rational::from_i64(3));
    }
}
