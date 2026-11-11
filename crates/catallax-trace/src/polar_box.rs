// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: a value type with the polar domain held as four
// numbers; statuses become returned values. Only the operations the line
// intersection reaches are carried, and the lazily built boxes of a surface
// cache live beside the cache instead of inside its nodes.
// Translated from SMLib SmPolarBox.cpp (AddVector3d, ComputePolarCoord,
// GetPolarArea, EvaluateNormalized, GetBoundaryVectors,
// HasPerpendicularToVector), SmBSplineSurface.cpp (the polar-box part of
// CalculateBoundingBox), NL_SrfBasic.cpp (N_SrfIsFlatCheap) and
// SmSurfaceCache.cpp (SmBezierAux2d::GetPolarBox, BuildPolarBox); rust-brep
// has no port.
//! Polar boxes: bounds on the directions of a surface's normals, as SMLib
//! keeps them on its subdivision nodes. The line intersection asks whether a
//! patch's normals can be perpendicular to the line.

use crate::extent::Extent2;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, acos, cos, sin};
use crate::solution::NodeId;
use crate::spline::BSplineSurface;
use crate::surface_cache::SurfaceCache;
use crate::surface_eval::evaluate_normal;
use crate::vector::{Vector2, Vector3};

/// `SmPolarBoxState`, from most degenerate to most general.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolarState {
    /// No vector yet.
    Uninitialized,
    /// The basis is set from a guess but no vector is in yet.
    PointWithBasis,
    /// One direction.
    SinglePoint,
    /// Directions on an arc in the `v0`/`v2` plane.
    Arc,
    /// Directions in a region of the sphere.
    Region,
    /// Some direction lies more than 90 degrees from `v0`.
    Unbounded,
}

/// `SmPolarBox`: a basis and the range of two angles from its first axis
/// (`x` about the second axis, `y` about the third).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarBox {
    /// The state.
    pub state: PolarState,
    /// `m_aBasis`.
    pub basis: [Vector3; 3],
    /// `m_sPolarDomain`: the low corner.
    pub min: Vector2,
    /// `m_sPolarDomain`: the high corner.
    pub max: Vector2,
}

impl Default for PolarBox {
    fn default() -> Self {
        Self {
            state: PolarState::Uninitialized,
            basis: [Vector3::default(); 3],
            min: Vector2::new(0.0, 0.0),
            max: Vector2::new(0.0, 0.0),
        }
    }
}

impl PolarBox {
    fn reset_domain(&mut self) {
        self.min = Vector2::new(0.0, 0.0);
        self.max = Vector2::new(0.0, 0.0);
    }

    fn add_angles(&mut self, p: Vector2) {
        if p.x < self.min.x {
            self.min.x = p.x;
        }
        if p.x > self.max.x {
            self.max.x = p.x;
        }
        if p.y < self.min.y {
            self.min.y = p.y;
        }
        if p.y > self.max.y {
            self.max.y = p.y;
        }
    }

    fn domain(&self) -> Extent2 {
        Extent2::new(self.min, self.max)
    }

    /// `SmPolarBox::AddVector3d`: grows the box to hold `vector`. The first
    /// vector sets the first axis, or with `guess` (a unit third axis) the
    /// first axis is the vector's projection across it.
    pub fn add_vector(&mut self, vector: Vector3, guess: Option<Vector3>) {
        if vector.length_squared() < EFF_ZERO_SQ {
            return;
        }
        match self.state {
            PolarState::Uninitialized => {
                self.state = PolarState::SinglePoint;
                self.reset_domain();
                let ortho = |b: &mut Self| {
                    if let Ok(axes) = vector.make_unit_ortho_vectors(None) {
                        b.basis = axes;
                    }
                };
                match guess {
                    None => ortho(self),
                    Some(guess) => {
                        let mut first = guess.cross(vector).cross(guess);
                        if first.length_squared() > EFF_ZERO_SQ {
                            let _ = first.unitize();
                            self.basis[0] = first;
                            self.basis[2] = guess;
                            self.basis[1] = self.basis[2].cross(self.basis[0]);
                            self.state = PolarState::PointWithBasis;
                            self.add_vector(vector, None);
                        } else {
                            ortho(self);
                        }
                    }
                }
            }
            PolarState::PointWithBasis | PolarState::SinglePoint => {
                let mut cross = vector.cross(self.basis[0]);
                let length = cross.length();
                let dot = self.basis[0].dot(vector);
                if length < EFF_ZERO {
                    if dot < EFF_ZERO {
                        self.state = PolarState::Unbounded;
                    }
                } else {
                    self.reset_domain();
                    if self.state == PolarState::SinglePoint {
                        let _ = cross.unitize();
                        self.basis[1] = cross;
                        self.basis[2] = self.basis[0].cross(self.basis[1]);
                        let _ = self.basis[2].unitize();
                    }
                    if dot < EFF_ZERO {
                        self.state = PolarState::Unbounded;
                    } else {
                        self.state = PolarState::Arc;
                        let angles = self.polar_coord(vector);
                        self.add_angles(angles);
                        if angles.y.abs() > EFF_ZERO_SQRT {
                            self.state = PolarState::Region;
                        }
                    }
                }
            }
            PolarState::Arc => {
                if self.basis[0].dot(vector) < EFF_ZERO {
                    self.state = PolarState::Unbounded;
                } else {
                    let angles = self.polar_coord(vector);
                    if angles.y.abs() > EFF_ZERO_SQRT {
                        self.state = PolarState::Region;
                    }
                    self.add_angles(angles);
                }
            }
            PolarState::Region => {
                let scaled_zero = EFF_ZERO * (1.0 + vector.max_dimension());
                if self.basis[0].dot(vector) < scaled_zero {
                    self.state = PolarState::Unbounded;
                } else {
                    let angles = self.polar_coord(vector);
                    self.add_angles(angles);
                    if (self.max.x - self.min.x).abs() < EFF_ZERO_SQRT {
                        let old = -self.basis[1];
                        self.basis[1] = self.basis[2];
                        self.basis[2] = old;
                        // SMLib resets the domain to itself here.
                        self.state = PolarState::Arc;
                    }
                }
            }
            PolarState::Unbounded => {}
        }
    }

    /// `ComputePolarCoord`: the two angles of `vector` from the first axis.
    fn polar_coord(&self, vector: Vector3) -> Vector2 {
        if matches!(
            self.state,
            PolarState::Uninitialized | PolarState::Unbounded
        ) {
            return Vector2::default();
        }
        let angle = |axis: Vector3| {
            let mut projection = axis.cross(vector).cross(axis);
            if projection.length_squared() < EFF_ZERO_SQ {
                return 0.0;
            }
            let _ = projection.unitize();
            let dot = projection.dot(self.basis[0]);
            let mut angle = if dot > 1.0 - EFF_ZERO {
                0.0
            } else if dot < -1.0 + EFF_ZERO {
                PI
            } else {
                acos(dot)
            };
            if self.basis[0].cross(projection).dot(axis) < 0.0 {
                angle = -angle;
            }
            angle
        };
        Vector2::new(angle(self.basis[1]), angle(self.basis[2]))
    }

    /// `GetPolarArea`: the area of the unit sphere the box covers.
    pub fn polar_area(&self) -> f64 {
        match self.state {
            PolarState::Region => {
                let (a1, a2) = (self.min, self.max);
                let [b0, b1, b2] = self.basis;
                let first = |x: f64| -sin(x) * b0 - cos(x) * b2;
                let second = |y: f64| -sin(y) * b0 + cos(y) * b1;
                let (n1min, n1max) = (first(a1.x), first(a2.x));
                let (n2min, n2max) = (second(a1.y), second(a2.y));
                let angle = |a: Vector3, b: Vector3| a.angle_between(b).unwrap_or(0.0);
                2.0 * PI
                    - angle(n1min, n2min)
                    - angle(n1min, -n2max)
                    - angle(n1max, n2max)
                    - angle(n1max, -n2min)
            }
            PolarState::Unbounded => 4.0 * PI,
            _ => 0.0,
        }
    }

    /// `EvaluateNormalized`: the unit direction at normalised angles.
    fn evaluate_normalized(&self, alpha: f64, beta: f64) -> Vector3 {
        match self.state {
            PolarState::PointWithBasis | PolarState::SinglePoint => self.basis[0],
            PolarState::Arc => {
                let p = self.domain().evaluate(alpha, beta);
                cos(p.x) * self.basis[0] - sin(p.x) * self.basis[2]
            }
            PolarState::Region => {
                let p = self.domain().evaluate(alpha, beta);
                let (ca1, sa1, ca2, sa2) = (cos(p.x), sin(p.x), cos(p.y), sin(p.y));
                if ca1 <= EFF_ZERO || ca2 <= EFF_ZERO {
                    return Vector3::new(0.0, 0.0, 0.0);
                }
                let mut vector =
                    ca1 * self.basis[0] - sa1 * self.basis[2] + sa2 * ca1 / ca2 * self.basis[1];
                let _ = vector.unitize();
                vector
            }
            _ => Vector3::new(0.0, 0.0, 0.0),
        }
    }

    /// `GetBoundaryVectors`: the one, two or four corner directions.
    fn boundary_vectors(&self) -> Vec<Vector3> {
        match self.state {
            PolarState::SinglePoint | PolarState::PointWithBasis => vec![self.basis[0]],
            PolarState::Arc => vec![
                self.evaluate_normalized(0.0, 0.0),
                self.evaluate_normalized(1.0, 1.0),
            ],
            PolarState::Region => vec![
                self.evaluate_normalized(0.0, 0.0),
                self.evaluate_normalized(1.0, 0.0),
                self.evaluate_normalized(1.0, 1.0),
                self.evaluate_normalized(0.0, 1.0),
            ],
            _ => Vec::new(),
        }
    }

    /// `HasPerpendicularToVector`: whether some direction in the box may be
    /// perpendicular to `test` (always, for an empty or unbounded box).
    pub fn has_perpendicular_to_vector(&self, test: Vector3) -> bool {
        if matches!(
            self.state,
            PolarState::Uninitialized | PolarState::Unbounded
        ) {
            return true;
        }
        let mut unit = test;
        let _ = unit.unitize();
        let corners = self.boundary_vectors();
        let limit = EFF_ZERO_SQRT * 100.0;
        if corners.len() == 1 {
            return corners[0].dot(unit).abs() < limit;
        }
        let mut sign = 0;
        for corner in corners {
            if corner.dot(unit).abs() < limit {
                return true;
            }
            if corner.dot(test) > limit {
                if sign < 0 {
                    return true;
                } else if sign == 0 {
                    sign = 1;
                }
            }
            if corner.dot(test) < -limit {
                if sign > 0 {
                    return true;
                } else if sign == 0 {
                    sign = -1;
                }
            }
        }
        false
    }
}

/// `N_CPtToPtEuclid`, which multiplies by the weight's reciprocal (unlike
/// `TO_EUCLID`).
fn euclid_by_reciprocal(surface: &BSplineSurface, i: usize, j: usize) -> Vector3 {
    let [x, y, z, w] = surface.raw_point([i, j]);
    let z = if z == crate::spline::NO_Z {
        None
    } else {
        Some(z)
    };
    if w == crate::spline::NO_WEIGHT {
        Vector3::new(x, y, z.unwrap_or(0.0))
    } else {
        let r = 1.0 / w;
        Vector3::new(x * r, y * r, z.map_or(0.0, |z| z * r))
    }
}

/// `N_SrfIsFlatCheap` with its flatness answer only: no closed direction,
/// and every control point within `eps` of the plane through the corners'
/// average with the normal of the corner midpoints' spans.
fn is_flat_cheap(surface: &BSplineSurface, eps: f64) -> bool {
    const MTOL: f64 = 1.0e-7;
    let [nu, nv] = surface.counts();
    let (n, m) = (nu - 1, nv - 1);
    let closed_u = (0..=m).all(|j| {
        euclid_by_reciprocal(surface, 0, j).distance(euclid_by_reciprocal(surface, n, j)) <= MTOL
    });
    let closed_v = (0..=n).all(|i| {
        euclid_by_reciprocal(surface, i, 0).distance(euclid_by_reciprocal(surface, i, m)) <= MTOL
    });
    if closed_u || closed_v {
        return false;
    }
    let p00 = euclid_by_reciprocal(surface, 0, 0);
    let p10 = euclid_by_reciprocal(surface, n, 0);
    let p01 = euclid_by_reciprocal(surface, 0, m);
    let p11 = euclid_by_reciprocal(surface, n, m);
    let combine4 =
        |f: fn(Vector3) -> f64| 0.25 * f(p00) + 0.25 * f(p10) + 0.25 * f(p01) + 0.25 * f(p11);
    let p = Vector3::new(combine4(|v| v.x), combine4(|v| v.y), combine4(|v| v.z));
    let a = Vector3::new(
        0.5 * p00.x + 0.5 * p10.x,
        0.5 * p00.y + 0.5 * p10.y,
        0.5 * p00.z + 0.5 * p10.z,
    );
    let b = Vector3::new(
        0.5 * p10.x + 0.5 * p11.x,
        0.5 * p10.y + 0.5 * p11.y,
        0.5 * p10.z + 0.5 * p11.z,
    );
    let normal = (a - p).cross(b - p);
    if normal.length() < MTOL {
        return false;
    }
    let mut unit = normal;
    if unit.unitize().is_err() {
        return false;
    }
    let d = -p.dot(unit);
    for i in 0..=n {
        for j in 0..=m {
            let q = euclid_by_reciprocal(surface, i, j);
            if (unit.x * q.x + unit.y * q.y + unit.z * q.z + d).abs() > eps {
                return false;
            }
        }
    }
    true
}

/// The polar-box part of `SmBSplineSurface::CalculateBoundingBox` (and of
/// `SmPlane`'s): the box of the normals over `domain`. Low-degree and
/// analytic surfaces sample ten normals in `domain` (one for a flat
/// bilinear patch); others take three normals there and the cross products
/// of neighbouring control-net spans of the whole surface. With `guess`,
/// a second box is built about the guess's axis and kept when smaller.
pub fn surface_polar_box(
    surface: &BSplineSurface,
    domain: Extent2,
    guess: Option<&PolarBox>,
) -> PolarBox {
    let guess = guess.filter(|g| g.state != PolarState::Uninitialized);
    let guess_axis = guess.map(|g| {
        if g.max.x - g.min.x > g.max.y - g.min.y {
            g.basis[2]
        } else {
            g.basis[1]
        }
    });
    let [du, dv] = surface.degrees();
    let [ku, kv] = surface.knots().map(<[f64]>::len);
    let is_plane = surface.plane_frame().is_some();
    let mut planar = du == 1 && ku <= 4 && dv == 1 && kv <= 4;
    let analytic = is_plane
        || !((du == 1 && ku > 4)
            || (du == 2 && ku > 6)
            || du > 2
            || (dv == 1 && kv > 4)
            || (dv == 2 && kv > 6)
            || dv > 2
            || (du == 2 && dv == 2));
    if planar && !is_flat_cheap(surface, 0.0001) {
        planar = false;
    }
    let normal = |u: f64, v: f64| evaluate_normal(surface, domain.evaluate(u, v), [true, true]);
    let mut polar = PolarBox::default();
    if analytic {
        let center = normal(0.5, 0.5);
        if planar {
            polar.add_vector(center, None);
            return polar;
        }
        let corners = [
            normal(0.0, 0.0),
            normal(1.0, 0.0),
            normal(1.0, 1.0),
            normal(0.0, 1.0),
        ];
        let [nu, nv] = surface.counts();
        let first = if nu - 1 <= 2 && nv - 1 <= 2 {
            (corners[0] + corners[1] + corners[2] + corners[3]) / 4.0
        } else {
            center
        };
        if first.length_squared() > EFF_ZERO_SQ {
            let sides = [
                normal(0.0, 0.5),
                normal(1.0, 0.5),
                normal(0.5, 0.0),
                normal(0.5, 1.0),
            ];
            let all = [
                first, sides[0], sides[1], sides[2], sides[3], corners[0], corners[1], corners[2],
                corners[3], center,
            ];
            for v in all {
                polar.add_vector(v, None);
            }
            // SMLib builds a trial box about the guess here but never uses
            // it for analytic surfaces.
        }
        return polar;
    }
    let mut trial = PolarBox::default();
    let add = |polar: &mut PolarBox, trial: &mut PolarBox, v: Vector3| {
        polar.add_vector(v, None);
        if guess_axis.is_some() {
            trial.add_vector(v, guess_axis);
        }
    };
    for (u, v) in [(0.5, 0.5), (1.0, 0.5), (0.5, 1.0)] {
        add(&mut polar, &mut trial, normal(u, v));
    }
    let [nu, nv] = surface.counts();
    let (n, m) = (nu - 1, nv - 1);
    let point = |i: usize, j: usize| {
        surface
            .control_point([i, j], crate::spline::PointForm::Euclidean)
            .map(|p| p.0)
            .unwrap_or_default()
    };
    for i in 0..=n {
        for j in 0..=m {
            let p = point(i, j);
            // SMLib tests `j > 1` and `i > 1` here, not `> 0`.
            let dv0 = (j > 1).then(|| p - point(i, j - 1));
            let dv1 = (j < m).then(|| point(i, j + 1) - p);
            for du in [
                (i > 1).then(|| p - point(i - 1, j)),
                (i < n).then(|| point(i + 1, j) - p),
            ]
            .into_iter()
            .flatten()
            {
                if let Some(dv0) = dv0 {
                    add(&mut polar, &mut trial, du.cross(dv0));
                }
                if let Some(dv1) = dv1 {
                    add(&mut polar, &mut trial, du.cross(dv1));
                }
            }
        }
    }
    if guess_axis.is_some() && trial.polar_area() < polar.polar_area() - EFF_ZERO {
        polar = trial;
    }
    polar
}

/// The polar boxes of a surface cache's nodes, built on demand as
/// `SmBezierAux2d::GetPolarBox` builds them: from the nearest ancestor that
/// has one (or the root) down to the node, each from the original surface
/// over the node's domain with its parent's box as the guess.
pub struct PolarBoxes<'a> {
    surface: &'a BSplineSurface,
    cache: &'a SurfaceCache,
    parents: Vec<Option<usize>>,
    boxes: Vec<Option<PolarBox>>,
}

impl<'a> PolarBoxes<'a> {
    /// Empty boxes for every node of `cache`, built from `surface`.
    pub fn new(surface: &'a BSplineSurface, cache: &'a SurfaceCache) -> Self {
        let count = cache.node_data.len();
        let mut parents = vec![None; count];
        let mut stack = vec![cache.top];
        while let Some(id) = stack.pop() {
            if let Some([a, b]) = cache.trees.get(id).and_then(|n| n.children) {
                parents[a.0 as usize] = Some(id.0 as usize);
                parents[b.0 as usize] = Some(id.0 as usize);
                stack.push(a);
                stack.push(b);
            }
        }
        Self {
            surface,
            cache,
            parents,
            boxes: vec![None; count],
        }
    }

    /// `GetPolarBox` of the node `id` of the cache's own tree.
    pub fn get(&mut self, id: usize) -> PolarBox {
        if let Some(found) = self.boxes[id] {
            return found;
        }
        let mut chain = vec![id];
        let mut start: Option<PolarBox> = None;
        let mut cursor = self.parents[id];
        while let Some(parent) = cursor {
            if let Some(found) = self.boxes[parent] {
                start = Some(found);
                break;
            }
            chain.push(parent);
            cursor = self.parents[parent];
        }
        let mut last = match start {
            Some(found) => found,
            None => {
                let top = chain.pop().expect("chain");
                let built = surface_polar_box(self.surface, self.domain(top), None);
                self.boxes[top] = Some(built);
                built
            }
        };
        while let Some(node) = chain.pop() {
            let guess = (last.state != PolarState::Unbounded).then_some(last);
            let built = surface_polar_box(self.surface, self.domain(node), guess.as_ref());
            self.boxes[node] = Some(built);
            last = built;
        }
        last
    }

    fn domain(&self, node: usize) -> Extent2 {
        self.cache.uv_domain(NodeId(node as u32))
    }
}
