// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: only the box operations the intersection chain
// enters; statuses become `Fault` and optional C++ outputs become returned
// values. `Extent2::clip_line` is translated directly from
// SmExtent2d::ClipLine2d.
// Adapted from rust-brep crates/kernel/src/extent1d.rs, extent2d.rs and
// extent3d.rs, commit 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, translations
// of SMLib SmExtent1d.cpp, SmExtent2d.cpp and SmExtent3d.cpp.
//! Intervals and axis-aligned boxes with SMLib's semantics.
//!
//! An extent's default is SMLib's empty box, minimum `1e20` and maximum
//! `-1e20`. The values `±1234567` and `±1111.111` mark unbounded sides.

use crate::Fault;
use crate::scalar::{BIG_DOUBLE, EFF_ZERO, EFF_ZERO_SQRT, PI, are_same, asin, max, min};
use crate::vector::{Vector2, Vector3};

const EMPTY_MIN: f64 = 1e20;
const EMPTY_MAX: f64 = -1e20;

/// `smos_3Max`, in SMLib's comparison order.
pub(crate) fn max3(a: f64, b: f64, c: f64) -> f64 {
    if a > b {
        if a > c { a } else { c }
    } else if b > c {
        b
    } else {
        c
    }
}

/// `smos_3Min`, in SMLib's comparison order.
pub(crate) fn min3(a: f64, b: f64, c: f64) -> f64 {
    if a < b {
        if a < c { a } else { c }
    } else if b < c {
        b
    } else {
        c
    }
}

/// The count rule shared by the line clips: one touch when the ends agree
/// within `SM_EFF_ZERO`, a span when ordered, nothing otherwise.
pub(crate) fn clip_result(enter: f64, exit: f64) -> Option<(f64, f64)> {
    if (enter - exit).abs() < EFF_ZERO {
        let mid = (enter + exit) / 2.0;
        Some((mid, mid))
    } else if enter < exit {
        Some((enter, exit))
    } else {
        None
    }
}

/// One slab of a line clip: the offsets of its two faces from the line's
/// point and the direction's component, all measured along the slab normal
/// (reversed, so `low > high`, for a left-handed pseudo-box basis).
#[derive(Clone, Copy, Debug)]
pub(crate) struct Slab {
    /// Offset of the low face; `-inf` for an unbounded side.
    pub low: f64,
    /// Offset of the high face; `+inf` for an unbounded side.
    pub high: f64,
    /// The direction's component along the normal.
    pub along: f64,
}

impl Slab {
    /// Whether SMLib's clip ignores this slab although the line is not
    /// parallel to it: a nonzero component no larger than `SM_EFF_ZERO`.
    pub(crate) fn is_ignored_but_crossed(self) -> bool {
        self.along != 0.0 && self.along.abs() <= EFF_ZERO
    }
}

/// Catallax's exact clip of a line against three slabs, used where SMLib's
/// clip ignores a slab the line still crosses: every nonzero component
/// bounds the span, and a parallel slab either keeps the whole line (its
/// point between the faces, within the scaled zero) or misses. Ends that
/// cross within the scaled zero touch at their middle, which is computed
/// without overflow when both ends sit near `f64::MAX`.
///
/// The span is the mathematical slab intersection, with no cap at
/// `±SM_BIG_DOUBLE`: a crossing of a huge box along a tiny component stays
/// representable. A side no slab bounds stays unbounded, which a span
/// cannot report, so such a clip is [`Fault::InvalidInput`] instead of a
/// narrowed answer.
pub(crate) fn exact_clip(slabs: [Slab; 3]) -> Result<Option<(f64, f64)>, Fault> {
    let (mut enter, mut exit) = (f64::NEG_INFINITY, f64::INFINITY);
    for slab in slabs {
        if slab.along == 0.0 {
            // A left-handed pseudo-box basis measures its faces in reverse.
            let (low, high) = (slab.low.min(slab.high), slab.low.max(slab.high));
            if low > EFF_ZERO * (1.0 + low.abs()) || high < -EFF_ZERO * (1.0 + high.abs()) {
                return Ok(None);
            }
            continue;
        }
        let (a, b) = (slab.low / slab.along, slab.high / slab.along);
        enter = enter.max(a.min(b));
        exit = exit.min(a.max(b));
    }
    if !(enter.is_finite() && exit.is_finite()) {
        return Err(Fault::InvalidInput);
    }
    if enter <= exit {
        return Ok(Some((enter, exit)));
    }
    if enter - exit <= EFF_ZERO * (1.0 + enter.abs().max(exit.abs())) {
        // Both ends can sit near `f64::MAX`, where the plain sum overflows;
        // half the difference is then safe and still lands between them.
        let mid = (enter + exit) / 2.0;
        let mid = if mid.is_finite() {
            mid
        } else {
            enter + (exit - enter) / 2.0
        };
        if !mid.is_finite() {
            return Err(Fault::InvalidInput);
        }
        return Ok(Some((mid, mid)));
    }
    Ok(None)
}

/// The smallest span covering both clips; a miss adds nothing.
pub(crate) fn clip_hull(a: Option<(f64, f64)>, b: Option<(f64, f64)>) -> Option<(f64, f64)> {
    match (a, b) {
        (Some(a), Some(b)) => Some((a.0.min(b.0), a.1.max(b.1))),
        (a, None) => a,
        (None, b) => b,
    }
}

/// `SM_IS_INFINITE`: SMLib's markers of an unbounded side.
fn infinite(x: f64) -> bool {
    x == 1234567.0 || x == -1234567.0 || x == 1111.111 || x == -1111.111
}

/// Where a value lies against an interval (`SmExtentPointType`).
/// The discriminants are SMLib's numbering.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum PointType {
    /// Within tolerance of the minimum only.
    Start = 0,
    /// Strictly between the ends.
    Inside = 1,
    /// Within tolerance of the maximum only.
    End = 2,
    /// Within tolerance of both ends.
    Both = 3,
    /// Outside the interval.
    Outside = 4,
}

/// A closed scalar interval (`SmExtent1d`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Extent1 {
    /// Lower end.
    pub min: f64,
    /// Upper end.
    pub max: f64,
}

impl Default for Extent1 {
    /// SMLib's empty interval.
    fn default() -> Self {
        Self {
            min: EMPTY_MIN,
            max: EMPTY_MAX,
        }
    }
}

impl Extent1 {
    /// The interval between `a` and `b`, in either order.
    pub fn new(a: f64, b: f64) -> Self {
        if a > b {
            Self { min: b, max: a }
        } else {
            Self { min: a, max: b }
        }
    }

    /// `IsInit`: still SMLib's empty interval.
    pub fn is_init(self) -> bool {
        self.min == EMPTY_MIN && self.max == EMPTY_MAX
    }

    /// `SetMinMax`; refuses a reversed pair.
    pub fn set_min_max(&mut self, min: f64, max: f64) -> Result<(), Fault> {
        if min > max {
            return Err(Fault::InvalidInput);
        }
        self.min = min;
        self.max = max;
        Ok(())
    }

    /// `AddValue`: grows the interval to include `x`.
    pub fn add_value(&mut self, x: f64) {
        if x < self.min {
            self.min = x;
        }
        if x > self.max {
            self.max = x;
        }
    }

    /// `GetMid`.
    pub fn mid(self) -> f64 {
        (self.max + self.min) / 2.0
    }

    /// `GetLength`.
    pub fn length(self) -> f64 {
        self.max - self.min
    }

    /// `GetMaxDimension`: the larger end magnitude.
    pub fn max_dimension(self) -> f64 {
        max(self.max.abs(), self.min.abs())
    }

    /// `ExpandAbsolute`: moves each bounded end outwards by `x`.
    pub fn expand_absolute(&mut self, x: f64) {
        self.min -= if infinite(self.min) { 0.0 } else { x };
        self.max += if infinite(self.max) { 0.0 } else { x };
    }

    /// `ClampValue`.
    pub fn clamp_value(self, x: f64) -> f64 {
        if x < self.min {
            self.min
        } else if x > self.max {
            self.max
        } else {
            x
        }
    }

    /// `Evaluate`: the value at normalised parameter `x`, exact at 0 and 1.
    pub fn evaluate(self, x: f64) -> f64 {
        if x == 0.0 {
            self.min
        } else if x == 1.0 {
            self.max
        } else {
            self.min + x * (self.max - self.min)
        }
    }

    /// `AreDisjoint`: separated by more than `tol`.
    pub fn are_disjoint(self, other: Self, tol: f64) -> bool {
        self.min > other.max - tol || self.max < other.min + tol
    }

    /// `Intersect`: the overlap; refuses disjoint intervals (`SM_ERR`).
    pub fn intersect(self, other: Self) -> Result<Self, Fault> {
        if self.are_disjoint(other, 0.0) {
            return Err(Fault::Degenerate);
        }
        Ok(Self {
            min: max(self.min, other.min),
            max: min(self.max, other.max),
        })
    }

    /// `AreEqual(SmExtent1d, tol)`.
    pub fn are_equal(self, other: Self, tol: f64) -> bool {
        (self.min - other.min).abs() <= tol.abs() && (self.max - other.max).abs() <= tol.abs()
    }

    /// `AreEqual(min, max, tol)`.
    pub fn equals_ends(self, min: f64, max: f64, tol: f64) -> bool {
        (self.min - min).abs() <= tol && (self.max - max).abs() <= tol
    }

    /// `IsContainedBy`.
    pub fn is_contained_by(self, other: Self, tol: f64) -> bool {
        other.min - tol <= self.min && other.max + tol >= self.max
    }

    /// `ContainsValue`.
    pub fn contains_value(self, x: f64, tol: f64) -> bool {
        x >= self.min - tol && x <= self.max + tol
    }

    /// `IsValueOnBoundary`: whether `x` is within `tol` of an end, and the
    /// other end (the farther one when `x` is near both).
    pub fn is_value_on_boundary(self, x: f64, tol: f64) -> (bool, f64) {
        let low = (x - self.min).abs();
        let high = (x - self.max).abs();
        let near_low = low < tol;
        let near_high = high < tol;
        let other = if near_low && near_high {
            if low < high { self.max } else { self.min }
        } else if near_low {
            self.max
        } else {
            self.min
        };
        (near_low || near_high, other)
    }

    /// `ClassifyPoint`.
    pub fn classify_point(self, x: f64, tol: f64) -> PointType {
        let low = (x - self.min).abs() <= tol;
        let high = (x - self.max).abs() <= tol;
        if low && high {
            PointType::Both
        } else if low {
            PointType::Start
        } else if high {
            PointType::End
        } else if x - self.min >= 0.0 && self.max - x >= 0.0 {
            PointType::Inside
        } else {
            PointType::Outside
        }
    }

    /// `IsBounded`: neither end is an unbounded marker.
    pub fn is_bounded(self) -> bool {
        !infinite(self.min) && !infinite(self.max)
    }
}

/// An axis-aligned parameter rectangle (`SmExtent2d`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Extent2 {
    /// Lower corner.
    pub min: Vector2,
    /// Upper corner.
    pub max: Vector2,
}

impl Default for Extent2 {
    /// SMLib's empty rectangle.
    fn default() -> Self {
        Self::from_corners(
            Vector2::new(EMPTY_MIN, EMPTY_MIN),
            Vector2::new(EMPTY_MAX, EMPTY_MAX),
        )
    }
}

impl Extent2 {
    /// The rectangle spanned by two points, each coordinate pair ordered.
    pub fn new(mut a: Vector2, mut b: Vector2) -> Self {
        if a.x > b.x {
            std::mem::swap(&mut a.x, &mut b.x);
        }
        if a.y > b.y {
            std::mem::swap(&mut a.y, &mut b.y);
        }
        Self { min: a, max: b }
    }

    /// The rectangle with these corners, kept in the given order.
    pub const fn from_corners(min: Vector2, max: Vector2) -> Self {
        Self { min, max }
    }

    /// `SetMinMax`; refuses a reversed pair.
    pub fn set_min_max(&mut self, min: Vector2, max: Vector2) -> Result<(), Fault> {
        if min.x > max.x || min.y > max.y {
            return Err(Fault::InvalidInput);
        }
        *self = Self { min, max };
        Ok(())
    }

    /// `GetSize`.
    pub fn size(self) -> Vector2 {
        self.max - self.min
    }

    /// `XLength`.
    pub fn x_length(self) -> f64 {
        self.max.x - self.min.x
    }

    /// `YLength`.
    pub fn y_length(self) -> f64 {
        self.max.y - self.min.y
    }

    /// `GetUInterval`.
    pub fn u_interval(self) -> Extent1 {
        Extent1::new(self.min.x, self.max.x)
    }

    /// `GetVInterval`.
    pub fn v_interval(self) -> Extent1 {
        Extent1::new(self.min.y, self.max.y)
    }

    /// `Union`.
    pub fn union(self, other: Self) -> Self {
        Self::from_corners(
            Vector2::new(min(self.min.x, other.min.x), min(self.min.y, other.min.y)),
            Vector2::new(max(self.max.x, other.max.x), max(self.max.y, other.max.y)),
        )
    }

    /// `Intersect`; refuses disjoint rectangles.
    pub fn intersect(self, other: Self) -> Result<Self, Fault> {
        if self.are_disjoint(other, 0.0) {
            return Err(Fault::InvalidInput);
        }
        Ok(Self::from_corners(
            Vector2::new(max(self.min.x, other.min.x), max(self.min.y, other.min.y)),
            Vector2::new(min(self.max.x, other.max.x), min(self.max.y, other.max.y)),
        ))
    }

    /// `AreDisjoint`: separated by more than `tol` in u or v.
    pub fn are_disjoint(self, other: Self, tol: f64) -> bool {
        self.min.x > other.max.x + tol
            || self.min.y > other.max.y + tol
            || self.max.x < other.min.x - tol
            || self.max.y < other.min.y - tol
    }

    /// `IsContainedBy`.
    pub fn is_contained_by(self, other: Self, tol: f64) -> bool {
        !(self.min.x < other.min.x - tol
            || self.min.y < other.min.y - tol
            || self.max.x > other.max.x + tol
            || self.max.y > other.max.y + tol)
    }

    /// `ContainsPoint2d` without the distance output.
    pub fn contains_point(self, p: Vector2, tol: f64) -> bool {
        !(p.x + tol < self.min.x
            || p.x - tol > self.max.x
            || p.y + tol < self.min.y
            || p.y - tol > self.max.y)
    }

    /// `ContainsPoint2d` with the distance output: whether `p` is inside
    /// within `tol`, and its distance to the rectangle when outside.
    pub fn contains_point_distance(self, p: Vector2, tol: f64) -> (bool, f64) {
        let mut inside = true;
        let mut squared = 0.0;
        if p.x + tol < self.min.x {
            inside = false;
            squared += (self.min.x - p.x) * (self.min.x - p.x);
        } else if p.x - tol > self.max.x {
            inside = false;
            squared += (p.x - self.max.x) * (p.x - self.max.x);
        }
        if p.y + tol < self.min.y {
            inside = false;
            squared += (self.min.y - p.y) * (self.min.y - p.y);
        } else if p.y - tol > self.max.y {
            inside = false;
            squared += (p.y - self.max.y) * (p.y - self.max.y);
        }
        (inside, crate::scalar::sqrt(squared))
    }

    /// `IsPoint2dOnBoundary`: whether `p` is within `tol` of a side, and the
    /// inward normal of the sides it is on.
    pub fn point_on_boundary(self, p: Vector2, tol: f64) -> (bool, Vector2) {
        let mut normal = Vector2::new(0.0, 0.0);
        let mut on = false;
        if (p.x - self.min.x).abs() <= tol {
            on = true;
            normal.x = 1.0;
        } else if (p.x - self.max.x).abs() <= tol {
            on = true;
            normal.x = -1.0;
        }
        if (p.y - self.min.y).abs() <= tol {
            on = true;
            normal.y = 1.0;
        } else if (p.y - self.max.y).abs() <= tol {
            on = true;
            normal.y = -1.0;
        }
        (on, normal)
    }

    /// `GetPoint2dBoundaries`: bits 1, 2, 4 and 8 for the u-min, v-min,
    /// u-max and v-max sides within `tol` of `p`.
    pub fn point_boundaries(self, p: Vector2, tol: f64) -> u32 {
        let mut bits = 0;
        if (p.x - self.min.x).abs() <= tol {
            bits |= 1;
        }
        if (p.y - self.min.y).abs() <= tol {
            bits |= 2;
        }
        if (p.x - self.max.x).abs() <= tol {
            bits |= 4;
        }
        if (p.y - self.max.y).abs() <= tol {
            bits |= 8;
        }
        bits
    }

    /// `ClampPoint2d`.
    pub fn clamp_point(self, p: Vector2) -> Vector2 {
        let mut value = p;
        if p.x < self.min.x {
            value.x = self.min.x;
        }
        if p.y < self.min.y {
            value.y = self.min.y;
        }
        if p.x > self.max.x {
            value.x = self.max.x;
        }
        if p.y > self.max.y {
            value.y = self.max.y;
        }
        value
    }

    /// `Evaluate`: the clamped point at normalised parameters, exact at 0 and 1.
    pub fn evaluate(self, u: f64, v: f64) -> Vector2 {
        let x = if u == 0.0 {
            self.min.x
        } else if u == 1.0 {
            self.max.x
        } else {
            self.min.x + u * (self.max.x - self.min.x)
        };
        let y = if v == 0.0 {
            self.min.y
        } else if v == 1.0 {
            self.max.y
        } else {
            self.min.y + v * (self.max.y - self.min.y)
        };
        self.clamp_point(Vector2::new(x, y))
    }

    /// `ClipLine2d`: the line parameter `param` if `point + param * vector`
    /// lies inside, else the parameter of the first side crossing. Returns 0
    /// when `point` is on a side and the line leaves through it, and
    /// `SM_BIG_DOUBLE` when a crossing cannot be computed. Meaningful only
    /// when `point` is inside the rectangle.
    pub fn clip_line(self, point: Vector2, vector: Vector2, param: f64) -> f64 {
        let sign = if param < 0.0 { -1.0 } else { 1.0 };
        let leave = EFF_ZERO * 100.0;
        if (are_same(self.min.x, point.x) && vector.x * sign < -leave)
            || (are_same(self.min.y, point.y) && vector.y * sign < -leave)
            || (are_same(self.max.x, point.x) && vector.x * sign > leave)
            || (are_same(self.max.y, point.y) && vector.y * sign > leave)
        {
            return 0.0;
        }
        let mut result = param;
        let tol = EFF_ZERO_SQRT / 10000.00;
        let mut end = point + param * vector;
        if self.contains_point(end, tol) {
            return result;
        }
        let side_x = if end.x < self.min.x - tol {
            Some(self.min.x)
        } else if end.x > self.max.x + tol {
            Some(self.max.x)
        } else {
            None
        };
        if let Some(side) = side_x {
            if (end.x - point.x).abs() < EFF_ZERO * (side - point.x).abs() {
                return BIG_DOUBLE;
            }
            result = result * (side - point.x) / (end.x - point.x);
            end = point + result * vector;
        }
        let side_y = if end.y < self.min.y - tol {
            Some(self.min.y)
        } else if end.y > self.max.y + tol {
            Some(self.max.y)
        } else {
            None
        };
        if let Some(side) = side_y {
            if (end.y - point.y).abs() < EFF_ZERO * (side - point.y).abs() {
                return BIG_DOUBLE;
            }
            result = result * (side - point.y) / (end.y - point.y);
        }
        result
    }
}

/// An axis-aligned box (`SmExtent3d`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Extent3 {
    /// Lower corner.
    pub min: Vector3,
    /// Upper corner.
    pub max: Vector3,
}

impl Default for Extent3 {
    /// SMLib's empty box.
    fn default() -> Self {
        Self {
            min: Vector3::new(EMPTY_MIN, EMPTY_MIN, EMPTY_MIN),
            max: Vector3::new(EMPTY_MAX, EMPTY_MAX, EMPTY_MAX),
        }
    }
}

fn coords(v: Vector3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn vector([x, y, z]: [f64; 3]) -> Vector3 {
    Vector3::new(x, y, z)
}

impl Extent3 {
    /// `IntersectLine`: the parameters where `point + t * vector` enters and
    /// leaves the box (one value for a single touch, `None` for a miss).
    /// When the middle of the clipped span is off the box, an end within
    /// the scaled zero of it still counts as a touch. A vanishing direction
    /// is [`Fault::InvalidInput`].
    ///
    /// SMLib ignores a direction component no larger than `SM_EFF_ZERO`,
    /// which along a long box can turn a crossing into a miss or a single
    /// touch. When it ignores a nonzero component, Catallax also clips
    /// exactly and answers the hull of both, so the clip
    /// never prunes a crossing; other lines keep SMLib's answer. An exact
    /// span unbounded on a side is [`Fault::InvalidInput`], since no
    /// parameter can report it.
    pub fn intersect_line(
        self,
        point: Vector3,
        vector: Vector3,
    ) -> Result<Option<(f64, f64)>, Fault> {
        if vector.length_squared() < crate::scalar::EFF_ZERO_SQ {
            return Err(Fault::InvalidInput);
        }
        let smlib = self.intersect_line_smlib(point, vector);
        let slab = |lo: f64, hi: f64, p: f64, v: f64| Slab {
            low: if infinite(lo) {
                f64::NEG_INFINITY
            } else {
                lo - p
            },
            high: if infinite(hi) { f64::INFINITY } else { hi - p },
            along: v,
        };
        let slabs = [
            slab(self.min.x, self.max.x, point.x, vector.x),
            slab(self.min.y, self.max.y, point.y, vector.y),
            slab(self.min.z, self.max.z, point.z, vector.z),
        ];
        if slabs.iter().any(|s| s.is_ignored_but_crossed()) {
            return Ok(clip_hull(smlib, exact_clip(slabs)?));
        }
        Ok(smlib)
    }

    /// `SmExtent3d::IntersectLine` as SMLib computes it.
    fn intersect_line_smlib(self, point: Vector3, vector: Vector3) -> Option<(f64, f64)> {
        let axis = |lo: f64, hi: f64, p: f64, v: f64| {
            if v > EFF_ZERO {
                (
                    if infinite(lo) {
                        -1234567.0
                    } else {
                        (lo - p) / v
                    },
                    if infinite(hi) {
                        1234567.0
                    } else {
                        (hi - p) / v
                    },
                )
            } else if v < -EFF_ZERO {
                (
                    if infinite(hi) {
                        -1234567.0
                    } else {
                        (hi - p) / v
                    },
                    if infinite(lo) {
                        1234567.0
                    } else {
                        (lo - p) / v
                    },
                )
            } else {
                (-BIG_DOUBLE, BIG_DOUBLE)
            }
        };
        let x = axis(self.min.x, self.max.x, point.x, vector.x);
        let y = axis(self.min.y, self.max.y, point.y, vector.y);
        let z = axis(self.min.z, self.max.z, point.z, vector.z);
        let enter = max3(z.0, x.0, y.0);
        let exit = min3(z.1, x.1, y.1);
        let middle = point + ((enter + exit) / 2.0) * vector;
        let tol = EFF_ZERO * (1.0 + middle.max_dimension());
        if self.contains_point(middle, tol) {
            return clip_result(enter, exit);
        }
        if self.contains_point(point + enter * vector, tol) {
            return Some((enter, enter));
        }
        if self.contains_point(point + exit * vector, tol) {
            return Some((exit, exit));
        }
        None
    }

    /// The box spanned by two points, each coordinate pair ordered.
    pub fn new(a: Vector3, b: Vector3) -> Self {
        let (mut low, mut high) = (coords(a), coords(b));
        for i in 0..3 {
            if low[i] > high[i] {
                std::mem::swap(&mut low[i], &mut high[i]);
            }
        }
        Self {
            min: vector(low),
            max: vector(high),
        }
    }

    /// The box with these corners, kept in the given order.
    pub const fn from_corners(min: Vector3, max: Vector3) -> Self {
        Self { min, max }
    }

    /// The box holding exactly `point`.
    pub fn singleton(point: Vector3) -> Self {
        Self {
            min: point,
            max: point,
        }
    }

    /// `AddPoint3d`: grows the box to include `point`.
    pub fn add_point(&mut self, point: Vector3) {
        let p = coords(point);
        let (mut low, mut high) = (coords(self.min), coords(self.max));
        for i in 0..3 {
            if p[i] < low[i] {
                low[i] = p[i];
            }
            if p[i] > high[i] {
                high[i] = p[i];
            }
        }
        (self.min, self.max) = (vector(low), vector(high));
    }

    /// `SetMinMax`; refuses a reversed pair.
    pub fn set_min_max(&mut self, min: Vector3, max: Vector3) -> Result<(), Fault> {
        if coords(min).into_iter().zip(coords(max)).any(|(a, b)| a > b) {
            return Err(Fault::InvalidInput);
        }
        *self = Self { min, max };
        Ok(())
    }

    /// `SetUInterval`: replaces the x range.
    pub fn set_u_interval(&mut self, interval: Extent1) {
        self.min.x = interval.min;
        self.max.x = interval.max;
    }

    /// `SetVInterval`: replaces the y range.
    pub fn set_v_interval(&mut self, interval: Extent1) {
        self.min.y = interval.min;
        self.max.y = interval.max;
    }

    /// `GetUInterval`: the x range.
    pub fn u_interval(self) -> Extent1 {
        Extent1::new(self.min.x, self.max.x)
    }

    /// `GetVInterval`: the y range.
    pub fn v_interval(self) -> Extent1 {
        Extent1::new(self.min.y, self.max.y)
    }

    /// `GetSize`.
    pub fn size(self) -> Vector3 {
        self.max - self.min
    }

    /// `GetMaxDimension`: the largest corner coordinate magnitude.
    pub fn max_dimension(self) -> f64 {
        max(self.min.max_dimension(), self.max.max_dimension())
    }

    /// `ExpandAbsolute`: moves every bounded side outwards by `amount`.
    pub fn expand_absolute(&mut self, amount: f64) {
        let (mut low, mut high) = (coords(self.min), coords(self.max));
        for i in 0..3 {
            low[i] -= if infinite(low[i]) { 0.0 } else { amount };
            high[i] += if infinite(high[i]) { 0.0 } else { amount };
        }
        (self.min, self.max) = (vector(low), vector(high));
    }

    /// `Union`.
    pub fn union(self, other: Self) -> Self {
        let (a, b) = (coords(self.min), coords(self.max));
        let (c, d) = (coords(other.min), coords(other.max));
        Self {
            min: vector(std::array::from_fn(|i| min(a[i], c[i]))),
            max: vector(std::array::from_fn(|i| max(b[i], d[i]))),
        }
    }

    /// `Intersect`; refuses disjoint boxes.
    pub fn intersect(self, other: Self) -> Result<Self, Fault> {
        if self.are_disjoint(other, 0.0) {
            return Err(Fault::InvalidInput);
        }
        let (a, b) = (coords(self.min), coords(self.max));
        let (c, d) = (coords(other.min), coords(other.max));
        Ok(Self {
            min: vector(std::array::from_fn(|i| max(a[i], c[i]))),
            max: vector(std::array::from_fn(|i| min(b[i], d[i]))),
        })
    }

    /// `AreDisjoint`: separated by more than `tol` on some axis.
    pub fn are_disjoint(self, other: Self, tol: f64) -> bool {
        coords(self.min)
            .into_iter()
            .zip(coords(other.max))
            .any(|(a, b)| a > b + tol)
            || coords(self.max)
                .into_iter()
                .zip(coords(other.min))
                .any(|(a, b)| a < b - tol)
    }

    /// `ContainsPoint3d`.
    pub fn contains_point(self, point: Vector3, tol: f64) -> bool {
        let p = coords(point);
        !(p.iter().zip(coords(self.min)).any(|(p, lo)| *p < lo - tol)
            || p.iter().zip(coords(self.max)).any(|(p, hi)| *p > hi + tol))
    }

    /// `IsPointSized`: every side length is at most `tol`.
    pub fn is_point_sized(self, tol: f64) -> bool {
        !coords(self.size()).into_iter().any(|d| d.abs() > tol)
    }

    /// `HasNegativeVolume`: some axis is reversed.
    pub fn has_negative_volume(self) -> bool {
        coords(self.min)
            .into_iter()
            .zip(coords(self.max))
            .any(|(a, b)| a > b)
    }

    /// `Evaluate`: the point at normalised parameters, exact at 0 and 1.
    pub fn evaluate(self, normalized: Vector3) -> Vector3 {
        let (lo, hi, p) = (coords(self.min), coords(self.max), coords(normalized));
        vector(std::array::from_fn(|i| {
            if p[i] == 0.0 {
                lo[i]
            } else if p[i] == 1.0 {
                hi[i]
            } else {
                lo[i] + p[i] * (hi[i] - lo[i])
            }
        }))
    }

    /// `ComputeSphereBound`: the centre, the radius of the enclosing sphere,
    /// and the radius of the largest sphere inside the box about the centre.
    pub fn sphere_bound(self) -> (Vector3, f64, f64) {
        let center = self.evaluate(Vector3::new(0.5, 0.5, 0.5));
        let corner = self.evaluate(Vector3::new(1.0, 1.0, 1.0));
        let radius = corner.distance(center);
        let inner = min(
            min((center.x - corner.x).abs(), (center.y - corner.y).abs()),
            (center.z - corner.z).abs(),
        );
        (center, radius, inner)
    }

    /// `ComputeConeBoundTo`: an axis and half-angle containing every
    /// direction from a point of `self` to a point of `other`. When the
    /// boxes overlap every direction is possible and the half-angle is pi.
    pub fn cone_bound_to(self, other: Self) -> Result<(Vector3, f64), Fault> {
        let half = self.size() / 2.0;
        let grown = Self::from_corners(other.min - half, other.max + half);
        let (center, radius, _) = grown.sphere_bound();
        let this_center = self.evaluate(Vector3::new(0.5, 0.5, 0.5));
        let direction = center - this_center;
        let distance = direction.length();
        let scale = 1.0 + center.max_dimension();
        if distance < radius + EFF_ZERO * scale {
            if direction.length_squared() > EFF_ZERO * scale {
                let mut axis = direction;
                axis.unitize().map_err(|_| Fault::Degenerate)?;
                return Ok((axis, PI));
            }
            return Ok((Vector3::new(0.0, 0.0, 1.0), PI));
        }
        Ok((direction, asin(radius / distance)))
    }
}
