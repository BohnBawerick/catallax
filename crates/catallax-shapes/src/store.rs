//! Handles, write-once arenas and the geometry store.

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

/// An append-only store. Entities never change once added.
#[derive(Clone, Debug)]
pub struct Arena<T> {
    values: Vec<T>,
}

impl<T> Default for Arena<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Arena<T> {
    /// An empty arena.
    pub fn new() -> Self {
        Self { values: Vec::new() }
    }

    /// Adds an entity and returns its handle; handles are issued in order from 0.
    pub fn push(&mut self, value: T) -> Handle<T> {
        let index = u32::try_from(self.values.len()).expect("arena handle space exhausted");
        self.values.push(value);
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
    implicits: Vec<RealPoly>,
    approximate_points: Vec<[f64; 3]>,
}

impl GeomStore {
    /// An empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a surface and records identities guaranteed by its definition.
    pub fn add_surface(&mut self, s: Surface) -> SurfaceId {
        self.add_surface_with_origin(s, "shapes.definition")
    }

    /// Adds a surface, retaining the constructing routine's provenance on
    /// structural facts. Supply this before insertion when provenance also
    /// distinguishes auxiliary geometry from intrinsic construction inputs.
    pub fn add_surface_with_origin(&mut self, s: Surface, origin: &'static str) -> SurfaceId {
        self.implicits.push(s.implicit());
        let id = self.surfaces.push(s);
        for index in 0..id.index {
            let other = Handle {
                index,
                _kind: PhantomData,
            };
            let a = self.surface(other);
            let b = self.surface(id);
            let coincident = crate::intent::coincident(a, b);
            let concentric = crate::intent::center(a)
                .zip(crate::intent::center(b))
                .is_some_and(|(a, b)| crate::intent::same_point(a, b));
            let parallel = crate::intent::direction(a)
                .zip(crate::intent::direction(b))
                .is_some_and(|(a, b)| crate::intent::parallel(a, b));
            let coaxial = crate::intent::axis(a)
                .zip(crate::intent::axis(b))
                .is_some_and(|((a, axis), (b, _))| parallel && crate::intent::on_axis(a, axis, b));
            for (kind, holds) in [
                (FactKind::Coincident, coincident),
                (FactKind::Concentric, concentric),
                (FactKind::Parallel, parallel),
                (FactKind::Coaxial, coaxial),
            ] {
                if holds {
                    self.construction(kind, GeomKey::Surface(other), GeomKey::Surface(id), origin);
                }
            }
        }
        if matches!(self.surface(id), Surface::Plane(_)) {
            for index in 0..self.points.len() {
                self.record_plane_incidence(
                    id,
                    Handle {
                        index: u32::try_from(index).expect("valid point index"),
                        _kind: PhantomData,
                    },
                    origin,
                );
            }
        }
        id
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
        if matches!(self.curve(id), Curve::Circle(_)) {
            for index in 0..self.surfaces.len() {
                let surface = Handle {
                    index: u32::try_from(index).expect("valid surface index"),
                    _kind: PhantomData,
                };
                let Curve::Circle(circle) = self.curve(id) else {
                    unreachable!()
                };
                let s = self.surface(surface);
                let coaxial = crate::intent::axis(s).is_some_and(|(origin, axis)| {
                    crate::intent::parallel(&circle.normal, axis)
                        && crate::intent::on_axis(origin, axis, &circle.center)
                });
                let concentric = crate::intent::center(s)
                    .is_some_and(|center| crate::intent::same_point(center, &circle.center));
                for (kind, holds) in [
                    (FactKind::Coaxial, coaxial),
                    (FactKind::Concentric, concentric),
                ] {
                    if holds {
                        self.construction(
                            kind,
                            GeomKey::Curve(id),
                            GeomKey::Surface(surface),
                            "shapes.definition",
                        );
                    }
                }
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
        self.approximate_points.push(approximate);
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
        for index in 0..self.surfaces.len() {
            self.record_plane_incidence(
                Handle {
                    index: u32::try_from(index).expect("valid surface index"),
                    _kind: PhantomData,
                },
                id,
                "shapes.definition",
            );
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
            // Full incidence plus a perpendicular section plane proves that
            // the circle's centre lies on the rotational carrier's axis.
            let coaxial = match surface {
                Surface::Cylinder(_) | Surface::Cone(_) | Surface::Torus(_) => {
                    crate::intent::axis(surface)
                        .is_some_and(|(_, axis)| crate::intent::parallel(&circle.normal, axis))
                }
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

    /// Implicit polynomial of a surface (cached).
    pub fn implicit(&self, s: SurfaceId) -> RealPoly {
        self.implicits[s.index as usize].clone()
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

    /// Gradient of the surface's implicit polynomial at a point.
    pub fn normal(&self, s: SurfaceId, p: PointId) -> Vector3 {
        self.surface(s).gradient(&self.position(p))
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
        match dec.sign(&surface.value(p)).map(|s| match s {
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
        let by_surfaces = match curve {
            Curve::Intersection(c) => c.surfaces.iter().all(|s| {
                self.facts
                    .holds(FactKind::Incident, GeomKey::Point(p), GeomKey::Surface(*s))
            }),
            _ => false,
        };
        if by_surfaces
            || self
                .facts
                .holds(FactKind::Incident, GeomKey::Point(p), GeomKey::Curve(c))
        {
            return fact_decision(dec);
        }
        let position = self.position(p);
        let coordinates = [position.x, position.y, position.z];
        let equations = curve.equations(self);
        decide_zeroes(
            equations.iter().map(|e| e.eval_real(&coordinates)),
            dec,
            |bounds, dec| match curve {
                Curve::Line(c) => Ok(math::sqrt_upper(&sum_squares(bounds))
                    * math::inverse_length_upper(&c.dir, dec)?),
                Curve::Circle(c) => {
                    let h_squared = &Real::from(&bounds[0] * &bounds[0]) / &c.normal.norm_squared();
                    let value = &Real::from(bounds[1].clone()) + &(&Real::from_i64(2) * &h_squared);
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
                (Surface::EllipticCylinder(_), Surface::EllipticCylinder(_)) => self.implicits
                    [a.index as usize]
                    .equals(&self.implicits[b.index as usize], dec)?,
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

    /// Sign of the dot product of the two surfaces' gradients at `at`. A
    /// not-exact outcome bounds the dot product itself.
    pub fn normals_agree(
        &self,
        a: SurfaceId,
        b: SurfaceId,
        at: PointId,
        dec: &mut Decider,
    ) -> Outcome<StrictSign> {
        let p = self.position(at);
        dec.sign(
            &self
                .surface(a)
                .gradient(&p)
                .dot(&self.surface(b).gradient(&p)),
        )
    }

    fn record_plane_incidence(&mut self, surface: SurfaceId, point: PointId, origin: &'static str) {
        let position = match self.point(point) {
            PointDef::Explicit(p) => p,
            PointDef::Meet(crate::Meet {
                position: Position::Closed(p),
                ..
            }) => p,
            _ => return,
        };
        if let Surface::Plane(plane) = self.surface(surface)
            && crate::intent::on_plane(&plane.origin, &plane.normal, position)
        {
            self.construction(
                FactKind::Incident,
                GeomKey::Point(point),
                GeomKey::Surface(surface),
                origin,
            );
        }
    }

    fn construction(&mut self, kind: FactKind, a: GeomKey, b: GeomKey, origin: &'static str) {
        self.facts.record(Fact {
            kind,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction(origin),
        });
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
