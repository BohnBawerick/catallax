//! A rolling-ball path piece compiled as a field leaf (4.2a steps 2 to 4).
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

use super::coaxial::{Coaxial, CoaxialCover};
use super::generator::{Generator, GeneratorCover};
use super::path::{self, Cover, Frame, Path};
use super::revolved::{Revolved, RevolvedCover, root};
use super::rim::{Rim, RimCover};
use super::saddle::{Saddle, SaddleCover};
use super::{Range, Relation};
use crate::Location;
use crate::metric::{Arithmetic, PRECISION};
use crate::support::Support;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Isometry, Point3, RPoint3, Surface, Vector3};

/// Runs `$body` with `$path` bound to the path of any [`Kind`].
macro_rules! each_kind {
    ($kind:expr, $path:ident => $body:expr) => {
        match $kind {
            Kind::Revolved($path) => $body,
            Kind::Generator($path) => $body,
            Kind::Coaxial($path) => $body,
            Kind::Saddle($path) => $body,
            Kind::Rim($path) => $body,
        }
    };
}
/// Runs `$body` with `$cover` bound to the cover of any [`PieceCover`].
macro_rules! each_cover {
    ($kind:expr, $cover:ident => $body:expr) => {
        match $kind {
            PieceCover::Revolved($cover) => $body,
            PieceCover::Generator($cover) => $body,
            PieceCover::Coaxial($cover) => $body,
            PieceCover::Saddle($cover) => $body,
            PieceCover::Rim($cover) => $body,
        }
    };
}

/// The path of a [`Piece`]. A piece lives boxed in its field leaf, so the
/// size gap between the variants costs nothing worth a second box.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug)]
pub enum Kind {
    /// A plane with a cylinder, cone or sphere.
    Revolved(Revolved),
    /// Supports parallel to one line.
    Generator(Generator),
    /// Two curved supports about one axis.
    Coaxial(Coaxial),
    /// A cylinder with a second cylinder or a sphere off its axis.
    Saddle(Saddle),
    /// A plane and a perpendicular cylinder, with the plane contact crossing
    /// circular mouths.
    Rim(Rim),
}

/// A compiled rolling-ball path piece: the material a convex fillet removes or
/// a concave one adds along one sharp edge run.
#[derive(Clone, Debug)]
pub struct Piece {
    kind: Kind,
    origin: Origin,
    convex: bool,
    supports: [Support; 2],
    // The supports' exact surfaces, for the side of a wall on one.
    carriers: [Surface; 2],
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
fn point_of(p: &[Rational; 3]) -> Point3 {
    Point3::new(
        p[0].clone().into(),
        p[1].clone().into(),
        p[2].clone().into(),
    )
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
        let hull = each_kind!(&kind, path => path.hull(dec)?);
        let mut piece = Self {
            kind,
            origin,
            convex,
            supports: supports.map(Support::new),
            carriers: supports.map(Surface::clone),
            bounds: Some(hull.map(|b| [b.lo(), b.hi()])),
        };
        piece.enclose_cover(dec)?;
        Ok(piece)
    }
    /// Encloses a straight piece's interval cover once, when the piece is
    /// built (4.2a step 4). Other paths enclose theirs per query.
    fn enclose_cover(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        if let Kind::Generator(path) = &mut self.kind {
            path.enclose_cover(dec)?;
        }
        Ok(())
    }
    /// After a rigid motion, which drops a straight piece's cover: the piece
    /// is derived again from its moved supports and encloses its cover.
    pub(crate) fn rebuild_after_transform(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        if let Kind::Generator(path) = &mut self.kind {
            path.rederive(dec)?;
        }
        Ok(())
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
    /// Boxes that together hold the region, tighter than [`Piece::bounds`]
    /// for a closed path.
    pub(crate) fn boxes(&self, dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        each_kind!(&self.kind, path => path.boxes(dec))
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
        each_kind!(&self.kind, path => path.radius())
    }

    /// The exact surface that holds the face, for the face registry (4.2a):
    /// a ring torus about a centre circle wider than the ball, a sphere where
    /// the centre stays put, or a cylinder about a straight centre line. Its
    /// equation has the sign of the path's separation, so it is positive in
    /// the region near the face and negative in every ball. A horn or spindle
    /// torus, an ellipse path's pipe, a saddle and a rim have no such
    /// surface; their faces stay unregistered.
    pub(crate) fn face(&self) -> Option<Surface> {
        match &self.kind {
            Kind::Revolved(path) => path.face(),
            Kind::Generator(path) => Some(path.face()),
            Kind::Coaxial(path) => Some(path.face()),
            Kind::Saddle(_) | Kind::Rim(_) => None,
        }
    }

    /// Whether a boundary point of the region lies on the face. A boundary
    /// point lies on a ball of the family (the face) or, outside every ball,
    /// on a wall; the closed-form separation tells which exactly. A path
    /// without a registered face answers no.
    pub(crate) fn on_face(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        if self.face().is_none() {
            return Ok(false);
        }
        let Some(value) = self.separation(&self.local(point)) else {
            return Ok(false);
        };
        Ok(crate::geometry::sign(&value, dec)? == Sign::Zero)
    }

    /// An exact point of the face near `point`, for a distance witness where
    /// the face's surface has no rational point near it (a ring torus of
    /// irrational radius). It is the surface's nearest point in closed form,
    /// as OpenCASCADE's `Extrema_ExtPElS` finds it for a sphere, cylinder
    /// or torus: along the ray from the sphere's centre, along the
    /// perpendicular from the cylinder's axis, and along the ray from the
    /// torus's core point in the meridian of `point`. On a ring torus that
    /// meridian point stays in the meridian's half-plane, so the point lies
    /// on the tube. The point is on the surface by construction, so the
    /// separation, which is the surface's equation up to a positive factor
    /// (see [`Piece::face`]), is zero there with no exact zero test. Where
    /// every band, wall and junction also encloses strictly above zero, the
    /// point is in the open fan, so it is a boundary point of the region on
    /// its face: [`path::classify`] answers `Boundary` there and
    /// [`Piece::on_face`] holds. `None` on an axis, core or centre, where a
    /// clip does not enclose above zero, or without a registered face.
    pub(crate) fn face_point(
        &self,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<Point3>, Refusal> {
        let positive = |value: &Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(crate::geometry::enclose(value, PRECISION, dec)?.lo() > Rational::zero())
        };
        // The part of `v` normal to `axis`, which need not be a unit.
        let across = |v: &Vector3, axis: &Vector3| {
            v.minus(&axis.scaled(&(&axis.dot(v) / &axis.norm_squared())))
        };
        let foot = match self.face() {
            Some(Surface::Sphere(s)) => {
                let v = point.minus(&s.center);
                let n2 = v.norm_squared();
                if !positive(&n2, dec)? {
                    return Ok(None);
                }
                s.center.plus(&v.scaled(&(&s.radius / &root(&n2))))
            }
            Some(Surface::Cylinder(s)) => {
                let w = across(&point.minus(&s.origin), &s.axis);
                let n2 = w.norm_squared();
                if !positive(&n2, dec)? {
                    return Ok(None);
                }
                point.plus(&w.scaled(&(&root(&(&s.radius_sq / &n2)) - &Real::one())))
            }
            Some(Surface::Torus(s)) => {
                if !positive(&(&s.major - &s.minor), dec)? {
                    return Ok(None);
                }
                let w = across(&point.minus(&s.center), &s.axis);
                let n2 = w.norm_squared();
                if !positive(&n2, dec)? {
                    return Ok(None);
                }
                let core = s.center.plus(&w.scaled(&(&s.major / &root(&n2))));
                let d = point.minus(&core);
                let d2 = d.norm_squared();
                if !positive(&d2, dec)? {
                    return Ok(None);
                }
                core.plus(&d.scaled(&(&s.minor / &root(&d2))))
            }
            _ => return Ok(None),
        };
        let local = self.local(&foot);
        let clips = [
            self.bands(&local, dec)?,
            self.walls(&local, dec)?,
            self.junctions(&local, dec)?,
        ];
        for value in clips.iter().flatten() {
            if !positive(value, dec)? {
                return Ok(None);
            }
        }
        Ok(Some(foot))
    }

    /// The region's side of `p + t d` for small `t > 0`. Off the boundary it
    /// is the membership itself. On the face, away from every wall and
    /// junction, the region near `p` is where the separation is positive,
    /// which is where the face's equation is (see [`Piece::face`]); the first
    /// nonzero coefficient of that equation along the ray decides it, as for
    /// a primitive's constraint. Off every ball on one wall, the region near
    /// `p` is one side of the support that wall lies on ([`Piece::wall_side`]),
    /// and the support's equation decides it the same way. Where a junction
    /// or two walls pass, or a wall meets the face, it stays unresolved.
    pub(crate) fn ray(
        &self,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let location = self.location(p, dec)?;
        if location != Location::Boundary {
            return Ok(location);
        }
        let local = self.local(&point_of(p));
        let walls = self.walls(&local, dec)?;
        let junctions = self.junctions(&local, dec)?;
        let mut clear = true;
        for value in walls.iter().chain(&junctions) {
            if crate::geometry::sign(value, dec)? != Sign::Positive {
                clear = false;
                break;
            }
        }
        let (equation, positive) = if clear {
            let Some(surface) = self.face() else {
                return Ok(Location::Indeterminate);
            };
            // A positive constraint along the ray leaves every ball, into
            // the region; a negative one enters a ball.
            (crate::geometry::constraint(surface)?, true)
        } else {
            match self.wall_side(p, &local, &walls, &junctions, dec)? {
                Some(side) => side,
                None => return Ok(Location::Indeterminate),
            }
        };
        Ok(
            match (
                crate::geometry::ray_constraints(&[equation], p, d, dec)?,
                positive,
            ) {
                (Location::Outside, true) => Location::Inside,
                (Location::Inside, true) => Location::Outside,
                (other, _) => other,
            },
        )
    }

    /// The polynomial of the one smooth surface bounding the region near a
    /// boundary point, where the region is one side of it: the face's
    /// equation away from every wall and junction, or the support under one
    /// wall off every ball ([`Piece::wall_side`]). `None` elsewhere.
    pub(crate) fn smooth_side(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Option<catallax_numbers::RealPoly>, Refusal> {
        let local = self.local(&point_of(p));
        let walls = self.walls(&local, dec)?;
        let junctions = self.junctions(&local, dec)?;
        let mut clear = true;
        for value in walls.iter().chain(&junctions) {
            if crate::geometry::sign(value, dec)? != Sign::Positive {
                clear = false;
                break;
            }
        }
        if clear {
            return match self.face() {
                Some(surface) => Ok(Some(crate::geometry::constraint(surface)?.polynomial)),
                None => Ok(None),
            };
        }
        Ok(self
            .wall_side(p, &local, &walls, &junctions, dec)?
            .map(|(constraint, _)| constraint.polynomial))
    }

    /// The support under a boundary point of the region that lies off every
    /// ball on exactly one wall, as its exact equation and whether the region
    /// near the point is where that equation is positive.
    ///
    /// The point is on the boundary off every ball, so the path's walls are
    /// faces there ([`Path::walls_are_faces`]): it lies on a face of the
    /// region, which off the balls is one of the two supports. With one
    /// support through the point and one wall zero there, that wall's zero
    /// set near the point is the support, so the region near the point is
    /// one side of it. Each ball touches both supports from the region's
    /// side, so that side is the one holding the ball centres, read at a
    /// station's centre. Returns `None` where a junction or more than one
    /// wall passes, on a ball, where both supports or neither pass, or where
    /// the centre's side is not settled by its enclosure.
    fn wall_side(
        &self,
        p: &[Rational; 3],
        local: &[Real; 3],
        walls: &[Real],
        junctions: &[Real],
        dec: &mut Decider,
    ) -> Result<Option<(crate::geometry::Constraint, bool)>, Refusal> {
        if !self.walls_are_faces() {
            return Ok(None);
        }
        for value in junctions {
            if crate::geometry::sign(value, dec)? != Sign::Positive {
                return Ok(None);
            }
        }
        let mut zeros = 0;
        for value in walls {
            match crate::geometry::sign(value, dec)? {
                Sign::Negative => return Ok(None),
                Sign::Zero => zeros += 1,
                Sign::Positive => {}
            }
        }
        if zeros != 1 {
            return Ok(None);
        }
        if let Some(value) = self.separation(local)
            && crate::geometry::sign(&value, dec)? != Sign::Positive
        {
            return Ok(None);
        }
        let mut found = None;
        for surface in &self.carriers {
            let constraint = crate::geometry::constraint(surface.clone())?;
            if crate::polynomial_sign(&constraint.polynomial, p, dec)? == Sign::Zero {
                if found.is_some() {
                    return Ok(None);
                }
                found = Some(constraint);
            }
        }
        let Some(constraint) = found else {
            return Ok(None);
        };
        let Some(centre) = self.centre_near(local, dec)? else {
            return Ok(None);
        };
        let side = constraint.polynomial.eval_ball(&centre, dec)?;
        Ok(if side.lo() > Rational::zero() {
            Some((constraint, true))
        } else if side.hi() < Rational::zero() {
            Some((constraint, false))
        } else {
            None
        })
    }

    /// A world enclosure of a ball centre: the closed-form nearest
    /// station's, or the first station the cover proves to exist among
    /// `u = 0, +-1/2, +-1`. Frames are in the path's local coordinates,
    /// which are world coordinates for every path but the rim's.
    fn centre_near(
        &self,
        local: &[Real; 3],
        dec: &mut Decider,
    ) -> Result<Option<[Ball; 3]>, Refusal> {
        let mut centre = None;
        if let Some(frame) = self.nearest(local, dec)? {
            centre = Some(frame.centre);
        } else {
            let cover = self.cover(dec)?;
            for (n, d) in [(0, 1), (1, 2), (-1, 2), (1, 1), (-1, 1)] {
                let u = Rational::from_ratio(n, d).expect("nonzero denominator");
                let station = Ball::from_rational(&u, PRECISION);
                if cover.exists(&station)? {
                    centre = Some(cover.centre(&station)?);
                    break;
                }
            }
        }
        match (centre, &self.kind) {
            (Some(centre), Kind::Rim(path)) => Ok(Some(path.world_point(&centre, dec)?)),
            (centre, _) => Ok(centre),
        }
    }

    /// The normal of the plane bounding the region near a boundary point,
    /// where the region is one side of a plane support there (see
    /// [`Piece::wall_side`]). It lets the plan settle such a point from the
    /// cells of the planes through it, as for a polyhedral leaf. Points off
    /// every plane support answer `false` before any other exact test.
    pub(crate) fn planes_at(
        &self,
        p: &[Rational; 3],
        planes: &mut Vec<Vector3>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut on_plane = false;
        for surface in &self.carriers {
            if let Surface::Plane(_) = surface {
                let constraint = crate::geometry::constraint(surface.clone())?;
                on_plane |= crate::polynomial_sign(&constraint.polynomial, p, dec)? == Sign::Zero;
            }
        }
        if !on_plane || self.location(p, dec)? != Location::Boundary {
            return Ok(false);
        }
        let local = self.local(&point_of(p));
        let walls = self.walls(&local, dec)?;
        let junctions = self.junctions(&local, dec)?;
        let Some((constraint, _)) = self.wall_side(p, &local, &walls, &junctions, dec)? else {
            return Ok(false);
        };
        let Surface::Plane(plane) = constraint.surface() else {
            return Ok(false);
        };
        let mut duplicate = false;
        for n in planes.iter() {
            if crate::geometry::sign(&plane.normal.cross(n).norm_squared(), dec)? == Sign::Zero {
                duplicate = true;
                break;
            }
        }
        if !duplicate {
            planes.push(plane.normal.clone());
        }
        Ok(true)
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
        each_kind!(&self.kind, path => path.fast(p, A::sqrt))
    }

    pub(crate) fn transform(&mut self, pose: &Isometry) {
        each_kind!(&mut self.kind, path => path.transform(pose));
        for support in &mut self.supports {
            support.transform(pose);
        }
        for carrier in &mut self.carriers {
            *carrier = pose.surface(carrier);
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
    /// Two curved supports about one axis.
    Coaxial(CoaxialCover),
    /// A cylinder with a second cylinder or a sphere off its axis.
    Saddle(SaddleCover),
    /// A ring with mouth contacts.
    Rim(RimCover),
}
impl Cover for PieceCover {
    fn radius(&self) -> &Ball {
        each_cover!(self, c => c.radius())
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        each_cover!(self, c => c.centre(u))
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        each_cover!(self, c => c.frame(u))
    }
    fn excludes(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        each_cover!(self, c => c.excludes(p, u))
    }
    fn below(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        each_cover!(self, c => c.below(p, u))
    }
    fn exists(&self, u: &Ball) -> Result<bool, Refusal> {
        each_cover!(self, c => c.exists(u))
    }
    fn empty(&self, u: &Ball) -> Result<bool, Refusal> {
        each_cover!(self, c => c.empty(u))
    }
}

impl Path for Piece {
    type Cover = PieceCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        each_kind!(&self.kind, path => path.local(point))
    }
    fn walls(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        each_kind!(&self.kind, path => path.walls(local, dec))
    }
    fn walls_are_faces(&self) -> bool {
        each_kind!(&self.kind, path => path.walls_are_faces())
    }
    fn bands(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        each_kind!(&self.kind, path => path.bands(local, dec))
    }
    fn junctions(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        each_kind!(&self.kind, path => path.junctions(local, dec))
    }
    fn uniform(&self) -> bool {
        each_kind!(&self.kind, path => path.uniform())
    }
    fn separation(&self, local: &[Real; 3]) -> Option<Real> {
        each_kind!(&self.kind, path => path.separation(local))
    }
    fn cover(&self, dec: &mut Decider) -> Result<PieceCover, Refusal> {
        Ok(match &self.kind {
            Kind::Revolved(path) => PieceCover::Revolved(path.cover(dec)?),
            Kind::Generator(path) => PieceCover::Generator(path.cover(dec)?),
            Kind::Coaxial(path) => PieceCover::Coaxial(path.cover(dec)?),
            Kind::Saddle(path) => PieceCover::Saddle(path.cover(dec)?),
            Kind::Rim(path) => PieceCover::Rim(path.cover(dec)?),
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        each_kind!(&self.kind, path => path.nearest(local, dec))
    }
    fn guess(&self, local: &[Real; 3]) -> Option<Rational> {
        each_kind!(&self.kind, path => path.guess(local))
    }
}
