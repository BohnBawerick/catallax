//! A rolling-ball path piece compiled as a field leaf (4.2a step 2).
//!
//! The leaf is the piece's region: membership from [`path::classify`], with a
//! membership the cover cannot settle reported as `Indeterminate`. Its value
//! is a residual, never a distance: a lower bound of the distance to the
//! region's boundary, which lies on the two untrimmed supports and the face.
//! The supports' 1-Lipschitz residuals bound the first two. The distance `m`
//! to the centre path is 1-Lipschitz and equals `r` on the face, so
//! `|m - r|` bounds the third; closed-form paths evaluate `m` at their
//! nearest station, ellipses bound `r - m` from the meridian station. Outside
//! the region, the distance to its box is a further lower bound.

use super::generator::{Generator, GeneratorCover};
use super::path::{self, Cover, Frame, Path};
use super::revolved::{Revolved, RevolvedCover};
use super::{Range, Relation};
use crate::Location;
use crate::metric::{Arithmetic, PRECISION};
use crate::support::Support;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_recipe::Origin;
use catallax_shapes::{Isometry, Point3, RPoint3, Surface};

/// The path of a [`Piece`].
#[derive(Clone, Debug)]
pub enum Kind {
    /// A plane with a cylinder, cone or sphere.
    Revolved(Revolved),
    /// Supports parallel to one line.
    Generator(Generator),
}

/// A compiled rolling-ball path piece: the material a convex fillet removes or
/// a concave one adds along one sharp edge run.
#[derive(Clone, Debug)]
pub struct Piece {
    kind: Kind,
    origin: Origin,
    convex: bool,
    supports: [Support; 2],
    // Rational box of the region, when its corners are rational.
    bounds: Option<[[Rational; 2]; 3]>,
}

fn magnitude(x: &Ball) -> Rational {
    let zero = Rational::zero();
    if x.lo() >= zero {
        x.lo()
    } else if x.hi() <= zero {
        -x.hi()
    } else {
        zero
    }
}
fn point(p: &RPoint3) -> Point3 {
    p.to_real()
}

impl Piece {
    /// Wraps a path derived by the field, with the two supports it rolls on.
    pub fn new(
        kind: Kind,
        origin: Origin,
        convex: bool,
        supports: [&Surface; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let hull = match &kind {
            Kind::Revolved(path) => path.hull(dec)?,
            Kind::Generator(path) => path.hull(dec)?,
        };
        Ok(Self {
            kind,
            origin,
            convex,
            supports: supports.map(Support::new),
            bounds: Some(hull.map(|b| [b.lo(), b.hi()])),
        })
    }
    /// The path.
    pub fn kind(&self) -> &Kind {
        &self.kind
    }
    /// The blend's name.
    pub fn origin(&self) -> &Origin {
        &self.origin
    }
    /// Whether the region is removed (convex edge) or added (concave edge).
    pub fn convex(&self) -> bool {
        self.convex
    }
    /// Rational box holding the region, when known.
    pub fn bounds(&self) -> Option<&[[Rational; 2]; 3]> {
        self.bounds.as_ref()
    }
    /// Strict membership in the region.
    pub fn classify(&self, point: &Point3, dec: &mut Decider) -> Result<Relation, Refusal> {
        path::classify(self, point, dec)
    }
    /// Guaranteed range for the distance to the piece's face.
    pub fn distance(
        &self,
        point: &Point3,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<Range, Refusal> {
        path::distance(self, point, width, dec)
    }
    fn radius(&self) -> &Real {
        match &self.kind {
            Kind::Revolved(path) => path.radius(),
            Kind::Generator(path) => path.radius(),
        }
    }

    pub(crate) fn location(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let point = Point3::new(
            p[0].clone().into(),
            p[1].clone().into(),
            p[2].clone().into(),
        );
        self.locate(&point, dec)
    }

    /// Membership at an exact point, with an unsettled cover `Indeterminate`.
    pub(crate) fn locate(&self, point: &Point3, dec: &mut Decider) -> Result<Location, Refusal> {
        match self.classify(point, dec) {
            Ok(Relation::Interior) => Ok(Location::Inside),
            Ok(Relation::Exterior) => Ok(Location::Outside),
            Ok(Relation::Boundary) => Ok(Location::Boundary),
            Err(e) if e.code == "field.rolling-region-unresolved" => Ok(Location::Indeterminate),
            Err(e) => Err(e),
        }
    }

    pub(crate) fn residual(
        &self,
        p: &RPoint3,
        location: Location,
        dec: &mut Decider,
    ) -> Result<Ball, Refusal> {
        if location == Location::Boundary {
            return Ok(Ball::from_rational(&Rational::zero(), PRECISION));
        }
        let balls = [&p.x, &p.y, &p.z].map(|x| Ball::from_rational(x, PRECISION));
        let mut lower = magnitude(&self.supports[0].proven(&balls, dec)?)
            .min(magnitude(&self.supports[1].proven(&balls, dec)?));
        let r = crate::geometry::enclose(self.radius(), PRECISION, dec)?;
        let local = self.local(&point(p));
        let distance = |c: &[Ball; 3]| {
            let d: [Ball; 3] = std::array::from_fn(|i| &balls[i] - &c[i]);
            (&(&d[0].square() + &d[1].square()) + &d[2].square()).sqrt()
        };
        let uniform = self.uniform();
        let face = if uniform && let Some(frame) = self.nearest(&local, dec)? {
            distance(&frame.centre).map(|m| magnitude(&(&m - &r)))
        } else if let Kind::Revolved(path) = &self.kind
            && let Some(station) = path.meridian_station(&local, dec)?
        {
            // The station bounds m from above, so r - m from below.
            distance(&station).map(|m| (&r - &m).lo().max(Rational::zero()))
        } else {
            None
        };
        lower = lower.min(face.unwrap_or_else(Rational::zero));
        if location == Location::Outside
            && let Some(bounds) = &self.bounds
        {
            let mut gap = Rational::zero();
            for (axis, [lo, hi]) in bounds.iter().enumerate() {
                let x = [&p.x, &p.y, &p.z][axis];
                let g = if x < lo {
                    lo - x
                } else if x > hi {
                    x - hi
                } else {
                    Rational::zero()
                };
                gap = gap + &g * &g;
            }
            let boxed = Ball::from_rational(&gap, PRECISION)
                .sqrt()
                .map(|b| b.lo().max(Rational::zero()))
                .unwrap_or_else(Rational::zero);
            lower = lower.max(boxed);
        }
        let lower = Ball::from_rational(&lower, PRECISION);
        Ok(match location {
            Location::Inside => -&lower,
            Location::Indeterminate => lower.hull(&-&lower),
            _ => lower,
        })
    }

    pub(crate) fn fast<A: Arithmetic>(&self, p: [f64; 3]) -> f64 {
        match &self.kind {
            Kind::Revolved(path) => path.fast(p, A::sqrt),
            Kind::Generator(path) => path.fast(p, A::sqrt),
        }
    }

    pub(crate) fn transform(&mut self, pose: &Isometry) {
        match &mut self.kind {
            Kind::Revolved(path) => path.transform(pose),
            Kind::Generator(path) => path.transform(pose),
        }
        for support in &mut self.supports {
            support.transform(pose);
        }
        // A rational map keeps the box's corners rational; otherwise the box
        // is dropped and only the weaker residual remains.
        self.bounds = self.bounds.take().and_then(|bounds| {
            let mut corners = Vec::with_capacity(8);
            for corner in 0..8 {
                let [x, y, z] = std::array::from_fn(|i| bounds[i][(corner >> i) & 1].clone());
                let image = pose.apply_point(&RPoint3::new(x, y, z).to_real());
                corners.push([
                    image.x.rational_ref()?.clone(),
                    image.y.rational_ref()?.clone(),
                    image.z.rational_ref()?.clone(),
                ]);
            }
            let mut moved: [[Rational; 2]; 3] =
                std::array::from_fn(|i| [corners[0][i].clone(), corners[0][i].clone()]);
            for corner in &corners {
                for (axis, value) in corner.iter().enumerate() {
                    moved[axis][0] = moved[axis][0].clone().min(value.clone());
                    moved[axis][1] = moved[axis][1].clone().max(value.clone());
                }
            }
            Some(moved)
        });
    }
}

/// The interval cover of a [`Piece`]. One is built per query and never
/// stored in bulk, so the size gap between the variants needs no box.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug)]
pub enum PieceCover {
    /// A plane with a cylinder, cone or sphere.
    Revolved(RevolvedCover),
    /// Supports parallel to one line.
    Generator(GeneratorCover),
}
impl Cover for PieceCover {
    fn radius(&self) -> &Ball {
        match self {
            Self::Revolved(c) => c.radius(),
            Self::Generator(c) => c.radius(),
        }
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        match self {
            Self::Revolved(c) => c.centre(u),
            Self::Generator(c) => c.centre(u),
        }
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        match self {
            Self::Revolved(c) => c.frame(u),
            Self::Generator(c) => c.frame(u),
        }
    }
}

impl Path for Piece {
    type Cover = PieceCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        match &self.kind {
            Kind::Revolved(path) => path.local(point),
            Kind::Generator(path) => path.local(point),
        }
    }
    fn walls(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        match &self.kind {
            Kind::Revolved(path) => path.walls(local, dec),
            Kind::Generator(path) => path.walls(local, dec),
        }
    }
    fn walls_are_faces(&self) -> bool {
        match &self.kind {
            Kind::Revolved(path) => path.walls_are_faces(),
            Kind::Generator(path) => path.walls_are_faces(),
        }
    }
    fn uniform(&self) -> bool {
        match &self.kind {
            Kind::Revolved(path) => path.uniform(),
            Kind::Generator(path) => path.uniform(),
        }
    }
    fn separation(&self, local: &[Real; 3]) -> Option<Real> {
        match &self.kind {
            Kind::Revolved(path) => path.separation(local),
            Kind::Generator(path) => path.separation(local),
        }
    }
    fn cover(&self, dec: &mut Decider) -> Result<PieceCover, Refusal> {
        Ok(match &self.kind {
            Kind::Revolved(path) => PieceCover::Revolved(path.cover(dec)?),
            Kind::Generator(path) => PieceCover::Generator(path.cover(dec)?),
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        match &self.kind {
            Kind::Revolved(path) => path.nearest(local, dec),
            Kind::Generator(path) => path.nearest(local, dec),
        }
    }
}
