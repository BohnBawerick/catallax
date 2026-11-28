//! Datum targets and the exact proof that one lies on its trimmed planar face.
//!
//! A target's point, line or area must stay inside the closed region its face
//! bounds, including holes cut through that face. The checks are exact: a
//! target that touches a hole rim from inside passes, one that crosses it is
//! refused, and nothing is snapped onto the face.

use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, StrictSign};
use catallax_shapes::distance::{PlanarBoundary, PlanarEdge};
use catallax_shapes::{Circle, Point3, Vector3};

/// A datum target: a point, line or area of one face that establishes part of
/// a declared datum (CAx-IF PMI recommended practices 4.1, section 6.6).
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct DatumTarget<F> {
    /// The letter of the datum the target belongs to, such as `A`.
    pub datum: String,
    /// The target's number within its datum: `1` in `A1`. Positive.
    pub number: u32,
    /// The face the target lies on.
    pub face: F,
    /// The target point, or the centre of a line or area, in millimetres in
    /// the part's coordinates.
    pub point_mm: [Rational; 3],
    /// The target's shape and size.
    pub kind: TargetKind,
}

/// The shape of a datum target. Sizes are full lengths in millimetres.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum TargetKind {
    /// A single point.
    Point,
    /// A line centred on the target point.
    Line {
        /// The full length.
        length_mm: Rational,
        /// The line's direction along the face, in part axes. Not a length.
        direction: [Rational; 3],
    },
    /// A bounded area centred on the target point.
    Area(TargetArea),
}

/// The footprint of an area datum target.
#[derive(Clone, PartialEq, Eq, Debug)]
#[allow(clippy::large_enum_variant)]
pub enum TargetArea {
    /// A disc.
    Circular {
        /// The full diameter.
        diameter_mm: Rational,
    },
    /// A rectangle.
    Rectangular {
        /// The full length, along `direction`.
        length_mm: Rational,
        /// The full width, across `direction` on the face.
        width_mm: Rational,
        /// The length axis along the face, in part axes. Not a length.
        direction: [Rational; 3],
    },
}

impl TargetKind {
    /// The name CAx-IF gives the implicit target type, written in
    /// `placed_datum_target_feature.description` (table 10).
    pub fn description(&self) -> &'static str {
        match self {
            TargetKind::Point => "point",
            TargetKind::Line { .. } => "line",
            TargetKind::Area(TargetArea::Circular { .. }) => "circle",
            TargetKind::Area(TargetArea::Rectangular { .. }) => "rectangle",
        }
    }

    /// The direction a line or rectangle runs in, if it has one.
    pub fn direction(&self) -> Option<&[Rational; 3]> {
        match self {
            TargetKind::Line { direction, .. }
            | TargetKind::Area(TargetArea::Rectangular { direction, .. }) => Some(direction),
            TargetKind::Point | TargetKind::Area(TargetArea::Circular { .. }) => None,
        }
    }

    /// Every length of the target, which must all be positive.
    pub fn sizes(&self) -> Vec<&Rational> {
        match self {
            TargetKind::Point => Vec::new(),
            TargetKind::Line { length_mm, .. } => vec![length_mm],
            TargetKind::Area(TargetArea::Circular { diameter_mm }) => vec![diameter_mm],
            TargetKind::Area(TargetArea::Rectangular {
                length_mm,
                width_mm,
                ..
            }) => vec![length_mm, width_mm],
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Sign {
    Negative,
    Zero,
    Positive,
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(value) = value.as_rational() {
        return Ok(match value.cmp(&Rational::zero()) {
            std::cmp::Ordering::Less => Sign::Negative,
            std::cmp::Ordering::Equal => Sign::Zero,
            std::cmp::Ordering::Greater => Sign::Positive,
        });
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(inexact) => Err(Refusal::from(inexact)),
        Outcome::Refused(refusal) => Err(refusal),
    }
}

/// `a < b`, decided exactly.
fn less(a: &Real, b: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&(b - a), dec)? == Sign::Positive)
}

fn real(value: &Rational) -> Real {
    Real::from_rational(value.clone())
}

fn half(value: &Rational) -> Real {
    real(&(value * &Rational::from_ratio(1, 2).expect("a nonzero denominator")))
}

/// A point of the part from millimetre coordinates.
pub(crate) fn point(at: &[Rational; 3]) -> Point3 {
    Point3::new(real(&at[0]), real(&at[1]), real(&at[2]))
}

/// A vector of the part from its components.
pub(crate) fn vector(at: &[Rational; 3]) -> Vector3 {
    Vector3::new(real(&at[0]), real(&at[1]), real(&at[2]))
}

fn length(vector: &Vector3) -> Real {
    let squared = vector.norm_squared();
    match squared.as_rational() {
        Some(value) => crate::number::sqrt(&value),
        None => squared.sqrt(),
    }
}

/// Whether the direction from an arc's centre, `v`, lies within the arc. The
/// arc runs counter-clockwise about the circle's normal from its first end to
/// its second; a whole circle has no ends.
fn arc_holds(
    circle: &Circle,
    ends: Option<[&Point3; 2]>,
    v: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let Some([a, b]) = ends else {
        return Ok(true);
    };
    let a = a.minus(&circle.center);
    let b = b.minus(&circle.center);
    let below = |v: &Vector3, dec: &mut Decider| -> Result<bool, Refusal> {
        let y = sign(&a.cross(v).dot(&circle.normal), dec)?;
        Ok(y == Sign::Negative || (y == Sign::Zero && sign(&a.dot(v), dec)? == Sign::Negative))
    };
    let b_below = below(&b, dec)?;
    let v_below = below(v, dec)?;
    Ok(if b_below != v_below {
        !v_below
    } else {
        sign(&v.cross(&b).dot(&circle.normal), dec)? != Sign::Negative
    })
}

fn arc_ends(boundary: &PlanarBoundary, ends: Option<[usize; 2]>) -> Option<[&Point3; 2]> {
    ends.map(|[a, b]| [&boundary.vertices[a], &boundary.vertices[b]])
}

fn segment_breaks(
    boundary: &PlanarBoundary,
    edge: &PlanarEdge,
    a: &Point3,
    b: &Point3,
    dec: &mut Decider,
) -> Result<Vec<Real>, Refusal> {
    let normal = &boundary.plane.normal;
    let d = b.minus(a);
    match edge {
        PlanarEdge::Segment([e, f]) => {
            let e = &boundary.vertices[*e];
            let f = &boundary.vertices[*f];
            let along = f.minus(e);
            let from = e.minus(a);
            let denominator = d.cross(&along).dot(normal);
            if sign(&denominator, dec)? != Sign::Zero {
                let t = from.cross(&along).dot(normal) / denominator.clone();
                let s = from.cross(&d).dot(normal) / denominator;
                return Ok(
                    (sign(&s, dec)? != Sign::Negative && !less(&Real::one(), &s, dec)?)
                        .then_some(t)
                        .into_iter()
                        .collect(),
                );
            }
            if sign(&from.cross(&d).dot(normal), dec)? != Sign::Zero {
                return Ok(Vec::new());
            }
            let norm = d.norm_squared();
            Ok(vec![from.dot(&d) / norm.clone(), f.minus(a).dot(&d) / norm])
        }
        PlanarEdge::Arc { circle, ends } => {
            let from = a.minus(&circle.center);
            let quadratic = d.norm_squared();
            let linear = from.dot(&d);
            let constant = &from.norm_squared() - &circle.radius_sq;
            let discriminant = &linear.square() - &(&quadratic * &constant);
            let discriminant_sign = sign(&discriminant, dec)?;
            if discriminant_sign == Sign::Negative {
                return Ok(Vec::new());
            }
            let root = if discriminant_sign == Sign::Zero {
                Real::zero()
            } else {
                discriminant.sqrt()
            };
            let mut parameters = Vec::new();
            for t in [
                (&(-&linear) - &root) / quadratic.clone(),
                (&(-&linear) + &root) / quadratic.clone(),
            ] {
                let v = from.plus(&d.scaled(&t));
                if arc_holds(circle, arc_ends(boundary, *ends), &v, dec)? {
                    parameters.push(t);
                }
            }
            Ok(parameters)
        }
    }
}

/// The squared distance from `c` to a boundary edge is at least `radius_sq`.
fn edge_clear_of_disc(
    boundary: &PlanarBoundary,
    edge: &PlanarEdge,
    c: &Point3,
    radius_sq: &Real,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let far = |p: &Point3, dec: &mut Decider| -> Result<bool, Refusal> {
        Ok(!less(&p.minus(c).norm_squared(), radius_sq, dec)?)
    };
    match edge {
        PlanarEdge::Segment([e, f]) => {
            let e = &boundary.vertices[*e];
            let f = &boundary.vertices[*f];
            let along = f.minus(e);
            let t = c.minus(e).dot(&along) / along.norm_squared();
            let closest = if sign(&t, dec)? != Sign::Positive {
                e.clone()
            } else if !less(&t, &Real::one(), dec)? {
                f.clone()
            } else {
                e.plus(&along.scaled(&t))
            };
            far(&closest, dec)
        }
        PlanarEdge::Arc { circle, ends } => {
            let ends_points = arc_ends(boundary, *ends);
            if let Some(points) = ends_points {
                for point in points {
                    if !far(point, dec)? {
                        return Ok(false);
                    }
                }
            }
            let towards = c.minus(&circle.center);
            let apart_sq = towards.norm_squared();
            if sign(&apart_sq, dec)? == Sign::Zero {
                // Every point of the arc is one radius from its centre.
                return Ok(!less(&circle.radius_sq, radius_sq, dec)?);
            }
            if !arc_holds(circle, ends_points, &towards, dec)? {
                return Ok(true);
            }
            // The arc's nearest point lies on the ray from its centre towards c.
            let gap = &length(&towards) - &circle.radius_sq.sqrt();
            Ok(!less(&gap.square(), radius_sq, dec)?)
        }
    }
}

/// A rectangle on the boundary's plane: centre, length axis and the axis
/// across it, each with its half extent scaled by the axis length.
struct Rectangle {
    centre: Point3,
    axes: [(Vector3, Real); 2],
}

impl Rectangle {
    fn corners(&self) -> [Point3; 4] {
        let [(u, hu), (w, hw)] = &self.axes;
        let u = u.scaled(&(hu / &u.norm_squared()));
        let w = w.scaled(&(hw / &w.norm_squared()));
        [
            self.centre.plus(&u).plus(&w),
            self.centre.plus(&u).plus(&w.scaled(&Real::from_i64(-1))),
            self.centre.plus(&u.scaled(&Real::from_i64(-1))).plus(&w),
            self.centre
                .plus(&u.scaled(&Real::from_i64(-1)))
                .plus(&w.scaled(&Real::from_i64(-1))),
        ]
    }

    /// Whether the closed segment from `a` to `b` meets the open rectangle.
    fn meets_segment(&self, a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let d = b.minus(a);
        let from = a.minus(&self.centre);
        let mut lows = vec![Real::zero()];
        let mut highs = vec![Real::one()];
        for (axis, extent) in &self.axes {
            // -extent < at + rate t < extent, with extent = half size times |axis|.
            let at = from.dot(axis);
            let rate = d.dot(axis);
            match sign(&rate, dec)? {
                Sign::Zero => {
                    if !less(&at.square(), &extent.square(), dec)? {
                        return Ok(false);
                    }
                }
                rate_sign => {
                    let low = (&(-extent) - &at) / rate.clone();
                    let high = (extent - &at) / rate;
                    if rate_sign == Sign::Positive {
                        lows.push(low);
                        highs.push(high);
                    } else {
                        lows.push(high);
                        highs.push(low);
                    }
                }
            }
        }
        for low in &lows {
            for high in &highs {
                if !less(low, high, dec)? {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// The squared distance from `q` to the closed rectangle, for `q` on its plane.
    fn distance_sq(&self, q: &Point3, dec: &mut Decider) -> Result<Real, Refusal> {
        let from = q.minus(&self.centre);
        let mut total = Real::zero();
        for (axis, extent) in &self.axes {
            // Coordinate along the axis, scaled by |axis|; compared unscaled.
            let along = from.dot(axis);
            let magnitude = if sign(&along, dec)? == Sign::Negative {
                -&along
            } else {
                along
            };
            let excess = &magnitude - extent;
            if sign(&excess, dec)? == Sign::Positive {
                total = &total + &(&excess.square() / &axis.norm_squared());
            }
        }
        Ok(total)
    }

    fn meets_edge(
        &self,
        boundary: &PlanarBoundary,
        edge: &PlanarEdge,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let PlanarEdge::Arc { circle, ends } = edge else {
            let PlanarEdge::Segment([a, b]) = edge else {
                unreachable!()
            };
            return self.meets_segment(&boundary.vertices[*a], &boundary.vertices[*b], dec);
        };
        let corners = self.corners();
        let mut inside_disc = true;
        for corner in &corners {
            inside_disc &= !less(
                &circle.radius_sq,
                &corner.minus(&circle.center).norm_squared(),
                dec,
            )?;
        }
        if inside_disc
            || !less(
                &self.distance_sq(&circle.center, dec)?,
                &circle.radius_sq,
                dec,
            )?
        {
            return Ok(false);
        }
        let Some([start, end]) = arc_ends(boundary, *ends) else {
            return Ok(true);
        };
        let mut candidates = vec![start.clone(), end.clone()];
        for (a, b) in [(0, 1), (0, 2), (1, 3), (2, 3)] {
            let from = &corners[a];
            let to = &corners[b];
            let side = to.minus(from);
            for t in segment_breaks(boundary, edge, from, to, dec)? {
                if sign(&t, dec)? != Sign::Negative && !less(&Real::one(), &t, dec)? {
                    candidates.push(from.plus(&side.scaled(&t)));
                }
            }
        }
        for point in candidates {
            let radial = point.minus(&circle.center);
            let tangent = circle.normal.cross(&radial);
            let from = point.minus(&self.centre);
            for forward in [true, false] {
                let at_end = if forward { end } else { start };
                if sign(&point.minus(at_end).norm_squared(), dec)? == Sign::Zero {
                    continue;
                }
                let mut enters = true;
                for (axis, extent) in &self.axes {
                    let at = from.dot(axis);
                    let delta = tangent.dot(axis);
                    let delta = if forward { delta } else { -&delta };
                    for positive in [true, false] {
                        let side = if positive { at.clone() } else { -&at };
                        match sign(&(&side - extent), dec)? {
                            Sign::Positive => enters = false,
                            Sign::Negative => {}
                            Sign::Zero => {
                                let speed = if positive { delta.clone() } else { -&delta };
                                let radial = radial.dot(axis);
                                let curvature = if positive { -&radial } else { radial };
                                enters &= match sign(&speed, dec)? {
                                    Sign::Negative => true,
                                    Sign::Positive => false,
                                    Sign::Zero => sign(&curvature, dec)? == Sign::Negative,
                                };
                            }
                        }
                    }
                }
                if enters {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

/// Why a datum target does not lie on its trimmed planar face, or `None` when
/// it does: its point on the face, a line's ends and the stretch between them,
/// or an area's whole footprint. `boundary` is the complete trim of the face,
/// in the target's millimetres. Refuses only when a decision cannot be made.
pub fn off_face(
    boundary: &PlanarBoundary,
    target: &DatumTarget<impl Sized>,
    dec: &mut Decider,
) -> Result<Option<&'static str>, Refusal> {
    if target
        .kind
        .sizes()
        .into_iter()
        .any(|size| size <= &Rational::zero())
    {
        return Ok(Some("a datum target's sizes must be more than zero"));
    }
    let plane = &boundary.plane;
    let centre = point(&target.point_mm);
    if sign(&centre.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero {
        return Ok(Some("the target point is not on its face's plane"));
    }
    if let Some(direction) = target.kind.direction() {
        let direction = vector(direction);
        if sign(&direction.norm_squared(), dec)? == Sign::Zero
            || sign(&direction.dot(&plane.normal), dec)? != Sign::Zero
        {
            return Ok(Some("the target direction does not run along its face"));
        }
    }
    if !boundary.contains(&centre, dec)? {
        return Ok(Some("the target point is outside its face"));
    }
    match &target.kind {
        TargetKind::Point => {}
        TargetKind::Line {
            length_mm,
            direction,
        } => {
            let direction = vector(direction);
            let offset = direction.scaled(&(&half(length_mm) / &length(&direction)));
            let ends = [
                centre.plus(&offset),
                centre.plus(&offset.scaled(&Real::from_i64(-1))),
            ];
            for end in &ends {
                if !boundary.contains(end, dec)? {
                    return Ok(Some("the target line leaves its face"));
                }
            }
            let direction = ends[1].minus(&ends[0]);
            let mut breaks = vec![Real::zero(), Real::one()];
            for edge in &boundary.edges {
                for t in segment_breaks(boundary, edge, &ends[0], &ends[1], dec)? {
                    for (index, next) in breaks.iter().enumerate().skip(1) {
                        match sign(&(&t - next), dec)? {
                            Sign::Zero => break,
                            Sign::Positive => continue,
                            Sign::Negative => {
                                if less(&breaks[index - 1], &t, dec)? {
                                    breaks.insert(index, t);
                                }
                                break;
                            }
                        }
                    }
                }
            }
            for pair in breaks.windows(2) {
                let mid = (&pair[0] + &pair[1]) / Real::from_i64(2);
                if !boundary.contains(&ends[0].plus(&direction.scaled(&mid)), dec)? {
                    return Ok(Some("the target line leaves its face"));
                }
            }
        }
        TargetKind::Area(TargetArea::Circular { diameter_mm }) => {
            let radius_sq = half(diameter_mm).square();
            for edge in &boundary.edges {
                if !edge_clear_of_disc(boundary, edge, &centre, &radius_sq, dec)? {
                    return Ok(Some("the target area crosses its face's boundary"));
                }
            }
        }
        TargetKind::Area(TargetArea::Rectangular {
            length_mm,
            width_mm,
            direction,
        }) => {
            let u = vector(direction);
            let w = plane.normal.cross(&u);
            let rectangle = Rectangle {
                axes: [
                    (u.clone(), &half(length_mm) * &length(&u)),
                    (w.clone(), &half(width_mm) * &length(&w)),
                ],
                centre,
            };
            for edge in &boundary.edges {
                if rectangle.meets_edge(boundary, edge, dec)? {
                    return Ok(Some("the target area crosses its face's boundary"));
                }
            }
        }
    }
    Ok(None)
}
