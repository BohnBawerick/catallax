//! Handles, write-once arenas and the geometry store.

use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

use catallax_numbers::{
    Ball, Coincidence, Decider, Fact, FactBook, FactKind, FactOrigin, MPoly, NotExact, Outcome,
    Rational, Real, Refusal, Sign, StrictSign, ZeroProof,
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
    implicits: Vec<MPoly>,
    approximate_points: Vec<[f64; 3]>,
}

impl GeomStore {
    /// An empty store.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a surface.
    pub fn add_surface(&mut self, s: Surface) -> SurfaceId {
        self.implicits.push(s.implicit());
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
            PointDef::Explicit(p) => [&p.x, &p.y, &p.z].map(Rational::to_f64_nearest),
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

    /// Records a fact; `false` if already present.
    pub fn record_fact(&mut self, fact: Fact<GeomKey>) -> bool {
        self.facts.record(fact)
    }

    /// Implicit polynomial of a surface (cached).
    pub fn implicit(&self, s: SurfaceId) -> MPoly {
        self.implicits[s.index as usize].clone()
    }

    /// Exact position of a point.
    pub fn position(&self, p: PointId) -> Point3 {
        match self.point(p) {
            PointDef::Explicit(p) => p.to_real(),
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
            Outcome::NotExact(mut n) => {
                n.bound = surface.distance_bound(&n.bound);
                Outcome::NotExact(n)
            }
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
            |bounds| match curve {
                Curve::Line(c) => {
                    math::sqrt_upper(&sum_squares(bounds)) * math::inverse_length_upper(&c.dir)
                }
                Curve::Circle(c) => {
                    let h_squared = &bounds[0] * &bounds[0] / c.normal.norm_squared();
                    math::sqrt_upper(&(&bounds[1] + &(Rational::from_i64(2) * h_squared)))
                }
                Curve::Ellipse(c) => {
                    &math::sqrt_upper(&c.major_sq) * &bounds[1]
                        + &bounds[0] * &math::inverse_length_upper(&c.normal)
                }
                Curve::Intersection(c) => self
                    .surface(c.surfaces[0])
                    .distance_bound(&bounds[0])
                    .max(self.surface(c.surfaces[1]).distance_bound(&bounds[1])),
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
        decide_zeroes([d.x, d.y, d.z], dec, |bounds| {
            math::sqrt_upper(&sum_squares(bounds))
        })
    }

    /// Whether two surfaces are the same point set. Facts first, then exact
    /// rational arithmetic on the definitions; never not exact or refused.
    pub fn same_surface(&self, a: SurfaceId, b: SurfaceId, dec: &mut Decider) -> Outcome<Apart> {
        if self.facts.holds(
            FactKind::Coincident,
            GeomKey::Surface(a),
            GeomKey::Surface(b),
        ) {
            return fact_decision(dec);
        }
        let equal = match (self.surface(a), self.surface(b)) {
            (Surface::Plane(a), Surface::Plane(b)) => {
                a.normal.cross(&b.normal).is_zero()
                    && b.origin.minus(&a.origin).dot(&a.normal).is_zero()
            }
            (Surface::Cylinder(a), Surface::Cylinder(b)) => {
                a.radius == b.radius
                    && a.axis.cross(&b.axis).is_zero()
                    && b.origin.minus(&a.origin).cross(&a.axis).is_zero()
            }
            (Surface::Cone(a), Surface::Cone(b)) => {
                a.apex == b.apex
                    && a.slope == b.slope
                    && a.axis.cross(&b.axis).is_zero()
                    && a.axis.dot(&b.axis).sign() == Sign::Positive
            }
            (Surface::Sphere(a), Surface::Sphere(b)) => {
                a.center == b.center && a.radius == b.radius
            }
            (Surface::Torus(a), Surface::Torus(b)) => {
                a.center == b.center
                    && a.major == b.major
                    && a.minor == b.minor
                    && a.axis.cross(&b.axis).is_zero()
            }
            _ => false,
        };
        if equal {
            Outcome::ProvenCoincident(Coincidence {
                proof: ZeroProof::ExactArithmetic,
            })
        } else {
            Outcome::Proven(Apart)
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
    distance_bound: impl FnOnce(&[Rational]) -> Rational,
) -> Outcome<Apart> {
    let mut bounds = std::array::from_fn::<_, 3, _>(|_| Rational::zero());
    let mut count = 0;
    let mut first_uncertain: Option<NotExact> = None;
    let mut coincidence = None;
    for (i, value) in values.into_iter().enumerate() {
        count = i + 1;
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
        n.bound = distance_bound(&bounds[..count]);
        Outcome::NotExact(n)
    } else {
        Outcome::ProvenCoincident(coincidence.expect("at least one equation"))
    }
}
