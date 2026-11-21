//! One general routine for a rolling-ball path piece (4.2a step 1).
//!
//! A path piece is a finite ball-centre curve `c(u)`, `u` in `[-1, 1]`, a ball
//! radius `r`, and the two contact feet `a(u)`, `b(u)` where the ball touches
//! its supports. The piece's face `F` is the union over `u` of the short arc
//! of the characteristic circle (centre `c(u)`, radius `r`, in the plane normal
//! to `c'(u)`) between the two feet. The piece's region is its fan (exact clip
//! walls supplied by the path) less every ball of the family.
//!
//! [`classify`] proves strict membership in that region, as the M18 rolling
//! region always did. [`distance`] returns a guaranteed range for the
//! Euclidean distance from a point to `F` (Q23 of the Checkpoint 2
//! decisions). Uniform paths with a closed-form nearest station (a straight
//! path with constant feet, a circle whose feet turn with its meridian)
//! evaluate it there, as OpenCASCADE's `Extrema_ExtPElC` does for elementary
//! curves. Other paths use a whole-interval panel cover of `[-1, 1]` with
//! lower bounds per panel and exact witness points on `F` at panel midpoints,
//! the interval branch-and-bound counterpart of OpenCASCADE's sampled
//! `Extrema_GExtPC`. The milestone note is `docs/M42A-FIELD-FILLETS.md`.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use super::{BITS, DEPTH, PANELS, Relation, gap, interval, q};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::Point3;

/// Enclosures at one station, or over a whole parameter panel, in the path's
/// local coordinates.
#[derive(Clone, Debug)]
pub struct Frame {
    /// Ball centre.
    pub centre: [Ball; 3],
    /// Any nonzero multiple of the centre's derivative.
    pub tangent: [Ball; 3],
    /// Contact feet on the two supports; the face arc runs between them.
    pub feet: [[Ball; 3]; 2],
}

/// Interval enclosures of a path over parameter panels inside `[-1, 1]`.
pub trait Cover {
    /// Enclosure of the ball radius.
    fn radius(&self) -> &Ball;
    /// Centre over the panel.
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal>;
    /// Centre, tangent and feet over the panel.
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal>;
    /// Exclusion witness at a station: the point lies in an open ball.
    fn excludes(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = self.centre(u)?;
        Ok(squared(p, &c).hi() < self.radius().square().lo())
    }
    /// Whole-panel clearance: the point lies outside every closed ball.
    fn below(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = self.centre(u)?;
        Ok(squared(p, &c).lo() > self.radius().square().hi())
    }
    /// Whether every parameter of the panel has a station. A path with gaps,
    /// such as a loop over part of its parameter circle, says no where it
    /// cannot prove it; its enclosures then still hold every station the
    /// panel has, but a midpoint frame is no face witness.
    fn exists(&self, _u: &Ball) -> Result<bool, Refusal> {
        Ok(true)
    }
    /// Whether the panel provably has no station at all.
    fn empty(&self, _u: &Ball) -> Result<bool, Refusal> {
        Ok(false)
    }
}

/// A finite rolling-ball path piece with exact clip walls.
pub trait Path {
    /// Interval cover type.
    type Cover: Cover;
    /// Exact local coordinates of a world point.
    fn local(&self, point: &Point3) -> [Real; 3];
    /// Exact clip values. The fan is where every value is non-negative. A
    /// path may decide an exact sign to choose its formula, so it takes the
    /// caller's decider.
    fn walls(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal>;
    /// Whether a zero wall value outside every ball lies on a face of the
    /// region. Exact elementary walls are faces; looser clipping walls are not.
    fn walls_are_faces(&self) -> bool {
        false
    }
    /// Exact clips that are positive everywhere in the region, so a point
    /// where one is zero or negative is outside it. They keep far zones of
    /// the walls out and are never faces.
    fn bands(&self, _local: &[Real; 3], _dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        Ok(Vec::new())
    }
    /// Exact clips where the piece meets a tangent neighbour piece (4.2a
    /// step 4). The region is where each is non-negative. A junction is no
    /// face: the neighbour's section there is the same, so a point on it
    /// inside the fan and outside every ball lies inside the two pieces'
    /// union, and a zero value does not make it a boundary point.
    fn junctions(&self, _local: &[Real; 3], _dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        Ok(Vec::new())
    }
    /// Closed-form value whose sign is that of `min_u |x-c(u)|^2 - r^2`.
    fn separation(&self, _local: &[Real; 3]) -> Option<Real> {
        None
    }
    /// The interval cover for membership, prepared under the caller's decider.
    fn cover(&self, dec: &mut Decider) -> Result<Self::Cover, Refusal>;
    /// The interval cover for distance queries, which also need feet and
    /// tangents. Defaults to [`Path::cover`].
    fn frames(&self, dec: &mut Decider) -> Result<Self::Cover, Refusal> {
        self.cover(dec)
    }
    /// Closed-form frame at the station nearest the point, clamped to the piece.
    fn nearest(&self, _local: &[Real; 3], _dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        Ok(None)
    }
    /// Whether the face is uniform along the path: a straight path whose feet
    /// keep their offsets from the centre, or a circle whose feet turn with
    /// its meridian. Then the face point nearest any point lies on the arc of
    /// the station [`Path::nearest`] returns: along a straight path the squared
    /// distance to a station's arc is the same section distance plus the
    /// squared axial gap, least at the clamped projection; on a circle the
    /// face is a surface of revolution, whose nearest point lies in the
    /// point's meridian half-plane. Only a uniform path's nearest station
    /// bounds the distance to the whole face. Elsewhere the nearest centre
    /// says nothing about the face: with feet that turn along a straight
    /// path, another station's arc can be much nearer.
    fn uniform(&self) -> bool {
        false
    }
    /// A parameter whose station's arc may lie near the point, from a
    /// binary64 search. [`distance`] evaluates the exact face point there
    /// first, so a good guess settles a query at once; a poor one only costs
    /// that evaluation. It never enters a bound unless its exact witness does.
    fn guess(&self, _local: &[Real; 3]) -> Option<Rational> {
        None
    }
}

/// Guaranteed range for the distance from a point to a piece's face, in mm.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Range {
    /// Proven lower bound.
    pub lower: Rational,
    /// Proven upper bound from an exact point on the face; `None` if no
    /// witness was certified.
    pub upper: Option<Rational>,
    /// Whether the requested width was reached before the refinement limits.
    pub converged: bool,
}

/// Membership and distance range of one point (Q23).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Answer {
    /// Strict membership, or an exact boundary on closed-form paths.
    pub relation: Relation,
    /// Guaranteed range for the distance to the piece's face.
    pub distance: Range,
}

fn sub(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}
fn cross(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}
pub(super) fn squared(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    let d = sub(a, b);
    dot(&d, &d)
}
fn root(x: &Ball) -> Result<Ball, Refusal> {
    x.sqrt()
        .ok_or_else(|| gap("negative squared distance enclosure"))
}
fn zero() -> Rational {
    q(0)
}

/// Bounds on the distance from `p` to the face arc of `frame`. For a point
/// station both bounds hold; for a panel only the lower bound is used.
///
/// The arc lies on the circle of radius `r` about the centre in the normal
/// plane. With `v` the offset from the centre, `v_n` its normal-plane part and
/// `v_t` its tangent part, the distance to that circle is
/// `sqrt((|v_n|-r)^2 + |v_t|^2)`. It is the distance to the arc when `v_n`
/// points into the arc's sector; otherwise the nearer foot is the nearest arc
/// point. Every foot lies on the face, so feet always give upper bounds.
fn arc_bounds(p: &[Ball; 3], frame: &Frame, r: &Ball) -> Result<(Rational, Rational), Refusal> {
    let feet = frame
        .feet
        .each_ref()
        .map(|f| root(&squared(p, f)))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let foot_lo = feet[0].lo().min(feet[1].lo());
    let foot_hi = feet[0].hi().min(feet[1].hi());
    let v = sub(p, &frame.centre);
    let tt = dot(&frame.tangent, &frame.tangent);
    if tt.lo() <= zero() {
        // No usable tangent: the sphere bound ||v|-r| is still a lower bound.
        let norm = root(&dot(&v, &v))?;
        let sphere = (&norm - r).lo().max(-(&norm - r).hi()).max(zero());
        return Ok((sphere.min(foot_lo), foot_hi));
    }
    let vt = dot(&v, &frame.tangent);
    let ratio = vt
        .checked_div(&tt)
        .ok_or_else(|| gap("tangent enclosure contains zero"))?;
    let normal: [Ball; 3] = std::array::from_fn(|i| &v[i] - &(&ratio * &frame.tangent[i]));
    // vt^2/tt is non-negative; drop the enclosure's negative part.
    let along = &vt * &ratio;
    let along = interval(&along.lo().max(zero()), &along.hi().max(zero()));
    let circle = root(&(&(&root(&dot(&normal, &normal))? - r).square() + &along))?;
    let a = sub(&frame.feet[0], &frame.centre);
    let b = sub(&frame.feet[1], &frame.centre);
    let n = cross(&a, &b);
    let first = dot(&cross(&a, &normal), &n);
    let second = dot(&cross(&normal, &b), &n);
    let inside = first.lo() > zero() && second.lo() > zero();
    let outside = first.hi() < zero() || second.hi() < zero();
    Ok(if inside {
        (circle.lo(), circle.hi())
    } else if outside {
        (foot_lo, foot_hi)
    } else {
        (circle.lo().min(foot_lo), foot_hi)
    })
}

/// Strict membership in the region of a path piece.
///
/// Every wall is decided exactly. A closed-form separation decides the balls
/// exactly, including the boundary. Otherwise a whole-interval cover of
/// `[-1, 1]` must prove each panel clear, and a midpoint witness inside a ball
/// proves exclusion. Sampling never proves inclusion.
pub fn classify<P: Path>(path: &P, point: &Point3, dec: &mut Decider) -> Result<Relation, Refusal> {
    dec.charge_work(1)?;
    let local = path.local(point);
    for c in path.bands(&local, dec)? {
        if sign(&c, dec)? != Sign::Positive {
            return Ok(Relation::Exterior);
        }
    }
    let mut wall = false;
    for c in path.walls(&local, dec)? {
        match sign(&c, dec)? {
            Sign::Negative => return Ok(Relation::Exterior),
            Sign::Zero => wall = true,
            Sign::Positive => {}
        }
    }
    for c in path.junctions(&local, dec)? {
        if sign(&c, dec)? == Sign::Negative {
            return Ok(Relation::Exterior);
        }
    }
    let unseparated = || gap("membership touches a clipping wall; incidence is not certified");
    if let Some(value) = path.separation(&local) {
        return match sign(&value, dec)? {
            Sign::Negative => Ok(Relation::Exterior),
            Sign::Zero => Ok(Relation::Boundary),
            Sign::Positive if !wall => Ok(Relation::Interior),
            Sign::Positive if path.walls_are_faces() => Ok(Relation::Boundary),
            Sign::Positive => Err(unseparated()),
        };
    }
    let p = [
        enclose(&local[0], BITS, dec)?,
        enclose(&local[1], BITS, dec)?,
        enclose(&local[2], BITS, dec)?,
    ];
    let cover = path.cover(dec)?;
    let mut pending = vec![(q(-1), q(1), 0)];
    let mut panels = 0;
    // A panel at the depth limit stays unresolved, but the rest of the cover
    // may still hold an exclusion witness: a point tangent to one ball can lie
    // inside another. Inclusion is refused if any panel stays unresolved.
    let mut unresolved = false;
    while let Some((lo, hi, depth)) = pending.pop() {
        dec.charge_work(1)?;
        panels += 1;
        if panels > PANELS {
            return Err(gap("finite centre-curve cover did not separate the query"));
        }
        let mid = (&lo + &hi) / q(2);
        let panel = interval(&lo, &hi);
        if cover.empty(&panel)? {
            continue;
        }
        // A whole-ball witness at a station that exists is a sound exclusion
        // from the region.
        let station = interval(&mid, &mid);
        if cover.exists(&station)? && cover.excludes(&p, &station)? {
            return Ok(Relation::Exterior);
        }
        if cover.below(&p, &panel)? {
            continue;
        }
        if depth == DEPTH {
            unresolved = true;
            continue;
        }
        pending.push((mid.clone(), hi, depth + 1));
        pending.push((lo, mid, depth + 1));
    }
    if unresolved {
        return Err(gap("finite centre-curve subdivision depth exhausted"));
    }
    match (wall, path.walls_are_faces()) {
        (false, _) => Ok(Relation::Interior),
        (true, true) => Ok(Relation::Boundary),
        (true, false) => Err(unseparated()),
    }
}

/// Guaranteed range for the distance from `point` to the piece's face.
///
/// Refinement stops once the range is at most `width` wide, or at the cover's
/// panel and depth limits; the range is guaranteed either way and
/// [`Range::converged`] says which. Panels are refined lowest bound first, ties
/// by parameter, so the result and its work do not depend on scheduling.
pub fn distance<P: Path>(
    path: &P,
    point: &Point3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Range, Refusal> {
    dec.charge_work(1)?;
    if width.sign() == Sign::Negative {
        return Err(gap("requested distance width must be non-negative"));
    }
    let local = path.local(point);
    let p = [
        enclose(&local[0], BITS, dec)?,
        enclose(&local[1], BITS, dec)?,
        enclose(&local[2], BITS, dec)?,
    ];
    if path.uniform()
        && let Some(frame) = path.nearest(&local, dec)?
    {
        let cover = path.frames(dec)?;
        let (lower, upper) = arc_bounds(&p, &frame, cover.radius())?;
        let converged = &upper - &lower <= *width;
        return Ok(Range {
            lower,
            upper: Some(upper),
            converged,
        });
    }
    let cover = path.frames(dec)?;
    let r = cover.radius().clone();
    // A panel with no station holds no face point; it leaves the heap.
    let lower_of = |lo: &Rational, hi: &Rational| -> Result<Option<Rational>, Refusal> {
        let u = interval(lo, hi);
        if cover.empty(&u)? {
            return Ok(None);
        }
        Ok(Some(arc_bounds(&p, &cover.frame(&u)?, &r)?.0))
    };
    let mut upper: Option<Rational> = None;
    if let Some(u) = path.guess(&local)
        && u >= q(-1)
        && u <= q(1)
    {
        let station = interval(&u, &u);
        if cover.exists(&station)? {
            upper = Some(arc_bounds(&p, &cover.frame(&station)?, &r)?.1);
        }
    }
    let mut heap = BinaryHeap::new();
    if let Some(bound) = lower_of(&q(-1), &q(1))? {
        heap.push(Reverse((bound, q(-1), q(1), 0usize)));
    }
    let mut panels = 0;
    while let Some(Reverse((lower, lo, hi, depth))) = heap.pop() {
        dec.charge_work(1)?;
        panels += 1;
        if upper.as_ref().is_some_and(|u| &(u - &lower) <= width) {
            return Ok(Range {
                lower,
                upper,
                converged: true,
            });
        }
        if panels > PANELS || depth == DEPTH {
            return Ok(Range {
                lower,
                upper,
                converged: false,
            });
        }
        let mid = (&lo + &hi) / q(2);
        let station = interval(&mid, &mid);
        if cover.exists(&station)? {
            let (_, witness) = arc_bounds(&p, &cover.frame(&station)?, &r)?;
            if upper.as_ref().is_none_or(|u| witness < *u) {
                upper = Some(witness);
            }
        }
        for (a, b) in [(lo, mid.clone()), (mid, hi)] {
            if let Some(bound) = lower_of(&a, &b)? {
                heap.push(Reverse((bound, a, b, depth + 1)));
            }
        }
    }
    // Every panel was empty, or the last ones were: no face point is left.
    Err(gap("the path has no station"))
}

/// Membership plus the guaranteed distance range of one point (Q23).
pub fn query<P: Path>(
    path: &P,
    point: &Point3,
    width: &Rational,
    dec: &mut Decider,
) -> Result<Answer, Refusal> {
    Ok(Answer {
        relation: classify(path, point, dec)?,
        distance: distance(path, point, width, dec)?,
    })
}
