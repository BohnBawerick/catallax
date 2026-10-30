//! Boundary-to-mesh coverage for an intersection of axis planes and one gyroid.
//!
//! The candidate mesh bounds a polyhedral solid `S`, so a segment from a point of
//! `S` to a point outside it meets the mesh. Coverage therefore needs no nearest
//! triangle. Each boundary point `p` receives two nearby points on opposite sides
//! of `S`, and the mesh lies within the longer of the two steps.
//!
//! Which side a point takes follows from how `S` is built. Take a grid tetrahedron
//! and the part of it on one side of its patch. That part lies in the hull of its
//! corners and patch vertices. On the side outside `S`, every corner has `g >= c`
//! when the whole tetrahedron is strictly inside the clip planes, and every patch
//! vertex is within `l/512` of such a point for an edge of length `l`. On the side
//! inside `S` the same holds with `g < c`, near the clip planes too. Linear
//! interpolation over a hull of radius `R` errs by at most `H R² / 2` for a Hessian
//! bound `H`, and `L` bounds the gradient. With `η = H R² / 2 + L l / 512`:
//!
//! - a point of the clean zone with `g < c - η` is inside `S`;
//! - a point with `g > c + η` is outside `S`.
//!
//! The clean zone is the union of grid cells strictly inside the clip planes. Any
//! point beyond a clip plane by more than a patch vertex can reach is also outside.
use super::*;
use crate::{gyroid::Gyroid, metric::Scalar};

type Vector = [f64; 3];

fn interval(lo: f64, hi: f64) -> Interval {
    Interval::new(lo, hi).unwrap()
}
fn length(d: Vector) -> Interval {
    (point(d[0]).square() + point(d[1]).square() + point(d[2]).square())
        .sqrt()
        .expect("sum of squares")
}

/// The clip box and gyroid of a supported plan.
pub(super) struct Clip<'a> {
    gyroid: &'a Gyroid,
    leaf: usize,
    lower: [Interval; 3],
    upper: [Interval; 3],
}

impl<'a> Clip<'a> {
    pub(super) fn new(classifier: &Classifier<'a>) -> Option<Self> {
        let plan = classifier.plan;
        if plan.instructions.iter().any(
            |i| matches!(i, Instruction::Boolean { kind, .. } if *kind != BooleanKind::Intersect),
        ) {
            return None;
        }
        let mut network = None;
        let mut lower: [Option<Interval>; 3] = [None; 3];
        let mut upper: [Option<Interval>; 3] = [None; 3];
        for (id, definition) in classifier.definitions.iter().enumerate() {
            match definition {
                None => (),
                Some(Definition::Gyroid) => {
                    let geometry::Leaf::Gyroid(gyroid) = &plan.leaves[id] else {
                        unreachable!("compiled gyroid");
                    };
                    if network.replace((id, gyroid.as_ref())).is_some() {
                        return None;
                    }
                }
                Some(Definition::Cell(supports)) => {
                    for support in supports {
                        let (axis, coefficient, position) = axis_plane(support)?;
                        // The support is negative inside, so a negative coefficient
                        // bounds the coordinate from below.
                        let side = if coefficient.hi() < 0.0 {
                            &mut lower[axis]
                        } else {
                            &mut upper[axis]
                        };
                        if side.is_some_and(|old| old != position) {
                            return None;
                        }
                        *side = Some(position);
                    }
                }
                Some(Definition::Primitive(block)) if block.kind == Kind::Block => {
                    // An unrotated block is `0 <= x + t <= size` on each axis.
                    let unrotated = (0..3).all(|i| {
                        (0..3).all(|j| block.matrix[i][j] == point(f64::from(u8::from(i == j))))
                    });
                    if !unrotated {
                        return None;
                    }
                    for axis in 0..3 {
                        let low = point(0.0) - block.translation[axis];
                        let high = block.parameters[axis] - block.translation[axis];
                        if lower[axis].replace(low).is_some_and(|old| old != low)
                            || upper[axis].replace(high).is_some_and(|old| old != high)
                        {
                            return None;
                        }
                    }
                }
                _ => return None,
            }
        }
        let (leaf, gyroid) = network?;
        Some(Self {
            gyroid,
            leaf,
            lower: [lower[0]?, lower[1]?, lower[2]?],
            upper: [upper[0]?, upper[1]?, upper[2]?],
        })
    }

    /// The gyroid leaf's index and definition.
    pub(super) fn gyroid(&self) -> (usize, &'a Gyroid) {
        (self.leaf, self.gyroid)
    }
}

/// The separation test over one candidate grid.
pub(super) struct Separation<'a, 'b> {
    clip: &'b Clip<'a>,
    grid: &'b refined::Grid,
    clean: [[f64; 2]; 3],
    curvature: Interval,
    lipschitz: Interval,
    bound: Interval,
}

impl<'a, 'b> Separation<'a, 'b> {
    /// The test for one grid, or `None` when no cell lies strictly inside the clip.
    pub(super) fn new(
        clip: &'b Clip<'a>,
        grid: &'b refined::Grid,
        bound: Interval,
    ) -> Option<Self> {
        Some(Self {
            clip,
            grid,
            clean: grid.clean_zone(&clip.lower, &clip.upper)?,
            curvature: clip.gyroid.curvature(),
            lipschitz: clip.gyroid.lipschitz(),
            bound,
        })
    }

    /// `η` over a box, for the grid cells that can hold one of its points.
    fn margin(&self, box_: &[Interval; 3]) -> Interval {
        let widths = std::array::from_fn::<f64, 3, _>(|i| self.grid.width(i, box_[i]));
        let diagonal = length(widths);
        let half = diagonal * point(0.5);
        point(0.5) * self.curvature * half.square() + self.lipschitz * diagonal * point(1.0 / 512.0)
    }

    /// The range of displacements that land every point of `box_ ∩ clip` in the
    /// clean zone.
    fn clearance(&self, box_: &[Interval; 3]) -> Option<[[f64; 2]; 3]> {
        let mut range = [[0.0; 2]; 3];
        for i in 0..3 {
            let lo =
                (point(self.clean[i][0]) - point(box_[i].lo().max(self.clip.lower[i].lo()))).hi();
            let hi =
                (point(self.clean[i][1]) - point(box_[i].hi().min(self.clip.upper[i].hi()))).lo();
            if lo > hi {
                return None;
            }
            range[i] = [lo, hi];
        }
        Some(range)
    }

    /// The longest step along `direction` within `range`, keeping clamped components.
    fn step(&self, direction: Vector, range: [[f64; 2]; 3]) -> Option<Vector> {
        let limit = self.bound.lo() * 0.999_999;
        let mut d = [0.0; 3];
        let mut free = [true; 3];
        for _ in 0..3 {
            let fixed: f64 = (0..3).filter(|&i| !free[i]).map(|i| d[i] * d[i]).sum();
            let remaining = limit * limit - fixed;
            if remaining < 0.0 {
                return None;
            }
            let wanted: f64 = (0..3)
                .filter(|&i| free[i])
                .map(|i| direction[i] * direction[i])
                .sum();
            if wanted <= 0.0 {
                break;
            }
            let scale = Scalar::sqrt(&(remaining / wanted));
            let mut clamped = false;
            for i in 0..3 {
                if free[i] {
                    d[i] = direction[i] * scale;
                    if d[i] < range[i][0] || d[i] > range[i][1] {
                        d[i] = d[i].clamp(range[i][0], range[i][1]);
                        free[i] = false;
                        clamped = true;
                    }
                }
            }
            if !clamped {
                break;
            }
        }
        Some(d)
    }

    /// Whether the box holds no solid at all: `g` rises away from one corner of the
    /// box's clipped part and that corner is not inside the network. A boundary point
    /// in such a box is a limit of solid points, which lie in a neighbouring box whose
    /// own test covers it. This settles a level set that touches a clip corner from
    /// outside the material, where neither side has a point to separate.
    fn solid_free(&self, box_: &[Interval; 3], dec: &mut Decider) -> Result<bool, Refusal> {
        let gyroid = self.clip.gyroid;
        let centre = box_.map(|b| b.lo() * 0.5 + b.hi() * 0.5);
        let Some(jet) = gyroid.jet(centre.map(point)) else {
            return Ok(false);
        };
        let radius = length(std::array::from_fn(|i| box_[i].hi() - box_[i].lo()));
        let drift = (self.curvature * radius).hi();
        let mut corner = [0.0; 3];
        for i in 0..3 {
            let slope = jet.gradient[i];
            corner[i] = if slope.lo() > drift {
                box_[i].lo().max(self.clip.lower[i].lo())
            } else if slope.hi() < -drift {
                box_[i].hi().min(self.clip.upper[i].hi())
            } else {
                return Ok(false);
            };
        }
        let Some(corner) = corner
            .iter()
            .map(|&v| Rational::from_f64_exact(v))
            .collect::<Option<Vec<_>>>()
        else {
            return Ok(false);
        };
        let corner: [Rational; 3] = corner.try_into().expect("three coordinates");
        Ok(gyroid.location(&corner, dec)? != crate::Location::Inside)
    }

    /// Whether each boundary point in `box_` has both a point of `S` and a point
    /// outside `S` within the bound.
    fn separated(&self, box_: &[Interval; 3]) -> bool {
        let gyroid = self.clip.gyroid;
        let Some(residual) = gyroid.residual(box_) else {
            return false;
        };
        let Some(range) = self.clearance(box_) else {
            return false;
        };
        let moved = |d: Vector| std::array::from_fn(|i| box_[i] + point(d[i]));
        if residual.hi() < 0.0 {
            // Only a clip face can bound this box. Step inward just far enough, and
            // a point beyond the plane by the patch-vertex reach is outside `S`.
            let d = range.map(|r| 0.0_f64.clamp(r[0], r[1]));
            let inner = moved(d);
            let Some(deeper) = gyroid.residual(&inner) else {
                return false;
            };
            return length(d).hi() < self.bound.lo() && (deeper + self.margin(&inner)).hi() < 0.0;
        }
        let centre = box_.map(|b| b.lo() * 0.5 + b.hi() * 0.5);
        let Some(jet) = gyroid.jet(centre.map(point)) else {
            return false;
        };
        let slope = jet.gradient.map(|v| v.mid());
        let size = Scalar::sqrt(&slope.iter().map(|v| v * v).sum::<f64>());
        if size <= 0.0 {
            return false;
        }
        let unit = slope.map(|v| v / size);
        let radius = length(std::array::from_fn(|i| {
            (point(box_[i].hi()) - point(centre[i]))
                .hi()
                .max((point(centre[i]) - point(box_[i].lo())).hi())
        }));
        // Change of `g` along a step `d` from any point of the box, about the
        // centre's gradient. The gradient moves by at most `H` per millimetre.
        let change = |d: Vector| {
            let rate = (0..3).fold(point(0.0), |sum, i| sum + jet.gradient[i] * point(d[i]));
            let steps = length(d);
            let slack = (self.curvature * steps * (radius + steps * point(0.5))).hi();
            rate + interval(-slack, slack)
        };
        let Some(inner) = self.step(unit.map(|v| -v), range) else {
            return false;
        };
        let outer = unit.map(|v| v * self.bound.lo() * 0.999_999);
        length(inner).hi() < self.bound.lo()
            && length(outer).hi() < self.bound.lo()
            && (change(inner) + self.margin(&moved(inner))).hi() < 0.0
            && (change(outer) - self.margin(&moved(outer))).lo() > 0.0
    }

    /// Covers every boundary point of the body within the bound of the mesh, or names
    /// a facet near the first box that cannot be covered.
    pub(super) fn cover(
        &self,
        mesh: &mut Mesh,
        classifier: &mut Classifier<'_>,
        bounds: &[[Rational; 2]; 3],
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        let root: [Interval; 3] = std::array::from_fn(|i| {
            Interval::new(
                Interval::from_rational(&bounds[i][0]).lo(),
                Interval::from_rational(&bounds[i][1]).hi(),
            )
            .unwrap()
        });
        let mut pending = vec![(root, 0u32)];
        while let Some((box_, depth)) = pending.pop() {
            spend(&mut mesh.cells, options, dec)?;
            if classifier.classify(&box_, dec)? != BoxState::Unknown
                || self.separated(&box_)
                || self.solid_free(&box_, dec)?
            {
                continue;
            }
            if depth == 40 {
                let centre = box_.map(|b| b.lo() * 0.5 + b.hi() * 0.5);
                return Ok(Some(Tree::new(mesh).nearest(mesh, centre)));
            }
            let axis = (0..3)
                .max_by(|&a, &b| {
                    (box_[a].hi() - box_[a].lo()).total_cmp(&(box_[b].hi() - box_[b].lo()))
                })
                .unwrap();
            let mid = box_[axis].lo() * 0.5 + box_[axis].hi() * 0.5;
            let mut a = box_;
            let mut b = box_;
            a[axis] = Interval::new(box_[axis].lo(), mid).unwrap();
            b[axis] = Interval::new(mid, box_[axis].hi()).unwrap();
            pending.push((a, depth + 1));
            pending.push((b, depth + 1));
        }
        Ok(None)
    }
}

/// Brackets a facet patch along the gyroid's gradient at its centre: moving the whole
/// patch one way proves it inside the solid, and the other way proves it outside.
/// Mesh facets are tilted against the axes, so this avoids subdividing them into
/// patches small enough for an axis-aligned box to clear the surface.
pub(super) fn bracket_along_gradient(
    clip: &Clip<'_>,
    patch: &Triangle,
    classifier: &mut Classifier<'_>,
    bound: Interval,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (leaf, gyroid) = clip.gyroid();
    let centre: Vector = std::array::from_fn(|i| (patch[0][i] + patch[1][i] + patch[2][i]) / 3.0);
    let Some(jet) = gyroid.jet(centre.map(point)) else {
        return Ok(false);
    };
    let slope = jet.gradient.map(|v| v.mid());
    let size = Scalar::sqrt(&slope.iter().map(|v| v * v).sum::<f64>());
    if size <= 0.0 {
        return Ok(false);
    }
    let reach = bound.lo() * 0.999_999;
    let shift: Vector = slope.map(|v| v / size * reach);
    if length(shift).hi() >= bound.lo() {
        return Ok(false);
    }
    let mut states = [BoxState::Unknown; 2];
    for (state, sign) in states.iter_mut().zip([-1.0, 1.0]) {
        let moved = patch.map(|p| {
            std::array::from_fn::<Interval, 3, _>(|i| point(p[i]) + point(shift[i]) * point(sign))
        });
        let Some(residual) = gyroid.facet_residual(&moved) else {
            return Ok(false);
        };
        let network = if residual.lo() > 0.0 {
            BoxState::Outside
        } else if residual.hi() < 0.0 {
            BoxState::Inside
        } else {
            return Ok(false);
        };
        let box_: [Interval; 3] = std::array::from_fn(|i| {
            interval(
                moved
                    .iter()
                    .map(|p| p[i].lo())
                    .fold(f64::INFINITY, f64::min),
                moved
                    .iter()
                    .map(|p| p[i].hi())
                    .fold(f64::NEG_INFINITY, f64::max),
            )
        });
        *state = classifier.classify_with(&box_, dec, |i| (i == leaf).then_some(network))?;
    }
    Ok(matches!(
        states,
        [BoxState::Inside, BoxState::Outside] | [BoxState::Outside, BoxState::Inside]
    ))
}
