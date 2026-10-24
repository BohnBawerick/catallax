//! Unit-density stand-in bodies for field-only features.
//!
//! Principal second moments determine two box dimensions directly. Volume
//! determines the third. A symmetric pair adds the remaining second moment
//! through the parallel-axis theorem, without a mass-convergence iteration.
//!
//! Rationalised principal axes and certified square-root midpoints propose
//! rational geometry. Its achieved inertia and residual are evaluated exactly;
//! callers must check the residual and the original field enclosures before
//! publishing. Decimal STEP rounding is a separate exchange error.

use catallax_numbers::{Ball, Rational, Refusal, Sign};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model};

/// Unit-density mass properties in model axes, inertia about the centroid.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct MassTarget {
    /// Volume, positive.
    pub volume: Rational,
    /// Centre of mass.
    pub centroid: [Rational; 3],
    /// Symmetric inertia tensor about the centroid, model axes.
    pub inertia: [[Rational; 3]; 3],
}

/// The built stand-in and its exact achieved properties.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct StandInBodies {
    /// One central block or a symmetric pair in deterministic construction order.
    pub bodies: Vec<BodyId>,
    /// Exact unit-density properties of the blocks together.
    pub achieved: MassTarget,
    /// Largest absolute difference between achieved and target inertia components.
    pub inertia_residual: Rational,
    /// Conservative Euclidean support radius about the reported centroid.
    pub support_radius: Rational,
}

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn zero3() -> [Rational; 3] {
    std::array::from_fn(|_| q(0))
}
fn zero33() -> [[Rational; 3]; 3] {
    std::array::from_fn(|_| zero3())
}
fn add33(a: &mut [[Rational; 3]; 3], b: &[[Rational; 3]; 3]) {
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = &a[i][j] + &b[i][j];
        }
    }
}
fn mat_mul(a: &[[Rational; 3]; 3], b: &[[Rational; 3]; 3]) -> [[Rational; 3]; 3] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..3).fold(q(0), |s, k| s + &a[i][k] * &b[k][j]))
    })
}
fn transpose(a: &[[Rational; 3]; 3]) -> [[Rational; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i].clone()))
}
fn mat_vec(a: &[[Rational; 3]; 3], v: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| (0..3).fold(q(0), |s, k| s + &a[i][k] * &v[k]))
}

/// Keeps about 48 significant bits of a positive double as a dyadic rational.
/// Uses only exponent arithmetic, so it is the same on every platform.
fn dyadic(x: f64) -> Rational {
    if !x.is_finite() || x <= 0.0 {
        return q(0);
    }
    let exponent = (x.to_bits() >> 52) as i64 - 1023;
    let shift = 48 - exponent;
    let scaled = if shift >= 0 {
        x * 2f64.powi(shift.min(1000) as i32)
    } else {
        x / 2f64.powi((-shift).min(1000) as i32)
    };
    let n = scaled.round() as i64;
    Rational::from_i64(n).mul_pow2(-shift)
}

/// Cyclic Jacobi eigen-decomposition of a symmetric 3x3 matrix in doubles.
/// Returns eigenvalues and a right-handed matrix whose columns are eigenvectors.
fn jacobi(m: [[f64; 3]; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let mut a = m;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..64 {
        let off = a[0][1].abs() + a[0][2].abs() + a[1][2].abs();
        if off == 0.0 {
            break;
        }
        for (p, r) in [(0, 1), (0, 2), (1, 2)] {
            if a[p][r] == 0.0 {
                continue;
            }
            let theta = (a[r][r] - a[p][p]) / (2.0 * a[p][r]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            for row in &mut a {
                let (akp, akr) = (row[p], row[r]);
                row[p] = c * akp - s * akr;
                row[r] = s * akp + c * akr;
            }
            let (before, after) = a.split_at_mut(r);
            for (apk, ark) in before[p].iter_mut().zip(&mut after[0]) {
                let (old_p, old_r) = (*apk, *ark);
                *apk = c * old_p - s * old_r;
                *ark = s * old_p + c * old_r;
            }
            for row in &mut v {
                let (vkp, vkr) = (row[p], row[r]);
                row[p] = c * vkp - s * vkr;
                row[r] = s * vkp + c * vkr;
            }
        }
    }
    let det = v[0][0] * (v[1][1] * v[2][2] - v[1][2] * v[2][1])
        - v[0][1] * (v[1][0] * v[2][2] - v[1][2] * v[2][0])
        + v[0][2] * (v[1][0] * v[2][1] - v[1][1] * v[2][0]);
    if det < 0.0 {
        for row in &mut v {
            row[2] = -row[2];
        }
    }
    ([a[0][0], a[1][1], a[2][2]], v)
}

/// Quaternion `[w, x, y, z]` of a proper rotation matrix, by Shepperd's method.
fn quaternion(r: &[[f64; 3]; 3]) -> [f64; 4] {
    let trace = r[0][0] + r[1][1] + r[2][2];
    if trace > 0.0 {
        let s = (trace + 1.0).sqrt() * 2.0;
        [
            0.25 * s,
            (r[2][1] - r[1][2]) / s,
            (r[0][2] - r[2][0]) / s,
            (r[1][0] - r[0][1]) / s,
        ]
    } else if r[0][0] > r[1][1] && r[0][0] > r[2][2] {
        let s = (1.0 + r[0][0] - r[1][1] - r[2][2]).sqrt() * 2.0;
        [
            (r[2][1] - r[1][2]) / s,
            0.25 * s,
            (r[0][1] + r[1][0]) / s,
            (r[0][2] + r[2][0]) / s,
        ]
    } else if r[1][1] > r[2][2] {
        let s = (1.0 + r[1][1] - r[0][0] - r[2][2]).sqrt() * 2.0;
        [
            (r[0][2] - r[2][0]) / s,
            (r[0][1] + r[1][0]) / s,
            0.25 * s,
            (r[1][2] + r[2][1]) / s,
        ]
    } else {
        let s = (1.0 + r[2][2] - r[0][0] - r[1][1]).sqrt() * 2.0;
        [
            (r[1][0] - r[0][1]) / s,
            (r[0][2] + r[2][0]) / s,
            (r[1][2] + r[2][1]) / s,
            0.25 * s,
        ]
    }
}

/// Exact dyadic quaternion components near a double quaternion.
fn dyadic_quaternion(v: [f64; 4]) -> [Rational; 4] {
    v.map(|x| {
        let magnitude = dyadic(x.abs());
        if x < 0.0 { -&magnitude } else { magnitude }
    })
}

/// A block in the principal frame: centre offset and full side lengths.
struct Block {
    offset: [Rational; 3],
    sides: [Rational; 3],
}
impl Block {
    fn mass(&self) -> Rational {
        &(&self.sides[0] * &self.sides[1]) * &self.sides[2]
    }
    /// Unit-density inertia about the principal-frame centroid.
    fn inertia(&self) -> [[Rational; 3]; 3] {
        let m = self.mass();
        let twelfth = Rational::from_ratio(1, 12).expect("nonzero denominator");
        let r = &self.offset;
        let rr = (0..3).fold(q(0), |s, i| s + &r[i] * &r[i]);
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                let own = if i == j {
                    let others = (0..3)
                        .filter(|k| *k != i)
                        .fold(q(0), |s, k| s + &self.sides[k] * &self.sides[k]);
                    &(&m * &others) * &twelfth
                } else {
                    q(0)
                };
                let shifted = if i == j {
                    &m * &(&rr - &(&r[i] * &r[i]))
                } else {
                    -&(&m * &(&r[i] * &r[j]))
                };
                own + shifted
            })
        })
    }
}

struct Layout {
    first: Block,
    second: Option<Block>,
}
impl Layout {
    fn blocks(&self) -> impl Iterator<Item = &Block> {
        std::iter::once(&self.first).chain(self.second.as_ref())
    }
}
struct Candidate {
    layout: Layout,
    inertia: [[Rational; 3]; 3],
    residual: Rational,
}
impl Candidate {
    fn new(
        layout: Layout,
        r: &[[Rational; 3]; 3],
        rt: &[[Rational; 3]; 3],
        target: &MassTarget,
    ) -> Self {
        let mut local = zero33();
        for block in layout.blocks() {
            add33(&mut local, &block.inertia());
        }
        let inertia = mat_mul(r, &mat_mul(&local, rt));
        let mut residual = q(0);
        for (i, row) in inertia.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                residual = residual.max((value - &target.inertia[i][j]).abs());
            }
        }
        Self {
            layout,
            inertia,
            residual,
        }
    }
}
fn square_root(value: &Rational) -> Result<Rational, Refusal> {
    Ball::from_rational(value, 192)
        .sqrt()
        .map(|root| root.mid())
        .ok_or_else(|| {
            Refusal::new(
                "stand-in.nonphysical-target",
                "a nonnegative spatial second moment is required",
            )
        })
}

/// Builds a closed-form proposal and reports its exact achieved properties.
/// A proposal outside the caller's property enclosures must not be exported.
pub fn equivalent_blocks(model: &mut Model, target: &MassTarget) -> Result<StandInBodies, Refusal> {
    if target.volume.sign() != Sign::Positive {
        return Err(Refusal::new(
            "stand-in.empty",
            "a mass stand-in needs a positive volume",
        ));
    }
    let f: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| target.inertia[i][j].to_f64_nearest()));
    if f.iter().flatten().any(|value| !value.is_finite()) {
        return Err(Refusal::new(
            "stand-in.numeric-range",
            "principal axes exceed the finite numeric proposal range",
        ));
    }
    let (_, axes) = jacobi(f);
    let [w, x, y, z] = dyadic_quaternion(quaternion(&axes));
    let rotation = Rotation::from_quaternion(w, x, y, z)?;
    let r: [[Rational; 3]; 3] = {
        let m = rotation
            .rational_matrix()
            .expect("rational quaternion gives a rational matrix");
        std::array::from_fn(|i| std::array::from_fn(|j| m[i][j].clone()))
    };
    let rt = transpose(&r);
    let local = mat_mul(&rt, &mat_mul(&target.inertia, &r));
    let variance: [Rational; 3] = std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        (&(&local[j][j] + &local[k][k]) - &local[i][i])
            .mul_pow2(-1)
            .checked_div(&target.volume)
            .expect("positive volume")
    });
    if variance.iter().any(|value| value.sign() != Sign::Positive) {
        return Err(Refusal::new(
            "stand-in.nonphysical-target",
            "positive spatial second moments are required",
        ));
    }
    let full_sides = [
        square_root(&(&variance[0] * &q(12)))?,
        square_root(&(&variance[1] * &q(12)))?,
        square_root(&(&variance[2] * &q(12)))?,
    ];
    let mut best: Option<Candidate> = None;
    for axis in 0..3 {
        let product = &full_sides[(axis + 1) % 3] * &full_sides[(axis + 2) % 3];
        let mut sides = full_sides.clone();
        sides[axis] = target.volume.checked_div(&product).expect("positive sides");
        let one = Candidate::new(
            Layout {
                first: Block {
                    offset: zero3(),
                    sides: sides.clone(),
                },
                second: None,
            },
            &r,
            &rt,
            target,
        );
        if best
            .as_ref()
            .is_none_or(|best| one.residual < best.residual)
        {
            best = Some(one);
        }
        if best.as_ref().is_some_and(|best| best.residual.is_zero()) {
            break;
        }

        // Each box has half the volume. Cross-axis variances are fixed by its
        // sides; the separation supplies the remaining variance on this axis.
        sides[axis] = sides[axis].mul_pow2(-1);
        let own = (&sides[axis] * &sides[axis])
            .checked_div(&q(12))
            .expect("nonzero denominator");
        let separation2 = &variance[axis] - &own;
        if separation2.sign() == Sign::Negative {
            continue;
        }
        let separation = square_root(&separation2)?;
        if separation.mul_pow2(1) < sides[axis] {
            continue;
        }
        let mut offset = zero3();
        offset[axis] = separation;
        let opposite = offset.each_ref().map(|value| -value);
        let pair = Candidate::new(
            Layout {
                first: Block {
                    offset,
                    sides: sides.clone(),
                },
                second: Some(Block {
                    offset: opposite,
                    sides,
                }),
            },
            &r,
            &rt,
            target,
        );
        if best
            .as_ref()
            .is_none_or(|best| pair.residual < best.residual)
        {
            best = Some(pair);
        }
    }
    let Candidate {
        layout,
        inertia,
        residual,
    } = best.expect("three positive single-box proposals");
    let mut bodies = Vec::with_capacity(1 + usize::from(layout.second.is_some()));
    let mut support_radius = q(0);
    for block in layout.blocks() {
        let half: [Rational; 3] = std::array::from_fn(|i| block.sides[i].mul_pow2(-1));
        // A rational rotation preserves Euclidean distance; the local L1 norm
        // bounds the support radius without an additional square root.
        let radius = (0..3).fold(q(0), |sum, i| &(&sum + &block.offset[i].abs()) + &half[i]);
        support_radius = support_radius.max(radius);
        let corner: [Rational; 3] = std::array::from_fn(|i| &block.offset[i] - &half[i]);
        let world = mat_vec(&r, &corner);
        let translation = RVector3::new(
            &target.centroid[0] + &world[0],
            &target.centroid[1] + &world[1],
            &target.centroid[2] + &world[2],
        );
        let pose = Pose::from_rational(rotation.clone(), translation);
        bodies.push(catallax_topology::block(model, block.sides.clone(), &pose)?);
    }
    Ok(StandInBodies {
        bodies,
        achieved: MassTarget {
            volume: target.volume.clone(),
            centroid: target.centroid.clone(),
            inertia,
        },
        inertia_residual: residual,
        support_radius,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }

    #[test]
    fn blocks_reproduce_a_rotated_plate_target() {
        // A 10 x 4 x 1 plate rotated 30 degrees about z, centred at (5, 6, 7).
        let (c, s) = (0.8660254037844386_f64, 0.5_f64);
        let (ix, iy, iz) = (40.0 / 12.0 * 17.0, 40.0 / 12.0 * 101.0, 40.0 / 12.0 * 116.0);
        let rot = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let mut world = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                world[i][j] = (0..3)
                    .map(|k| rot[i][k] * [ix, iy, iz][k] * rot[j][k])
                    .sum::<f64>();
            }
        }
        let target = MassTarget {
            volume: q(40),
            centroid: [q(5), q(6), q(7)],
            inertia: std::array::from_fn(|i| {
                std::array::from_fn(|j| Rational::from_f64_exact(world[i][j]).unwrap())
            }),
        };
        let mut model = Model::new();
        let built = equivalent_blocks(&mut model, &target).unwrap();
        assert_eq!(built.achieved.volume, q(40));
        assert_eq!(built.achieved.centroid, target.centroid);
        assert!(
            built.inertia_residual < r("1e-9"),
            "residual {}",
            built.inertia_residual
        );
        // The exact measurement of the blocks agrees with the reported numbers.
        let accuracy = catallax_measure::Accuracy::default();
        let mut dec = catallax_numbers::Decider::default();
        let mut volume = q(0);
        for body in &built.bodies {
            let v = catallax_measure::volume(&model, *body, &accuracy, &mut dec)
                .unwrap()
                .value;
            assert!(v.is_exact());
            volume = volume + v.lo;
        }
        assert_eq!(volume, q(40));
    }

    #[test]
    fn a_compact_cube_target_is_also_solvable() {
        let target = MassTarget {
            volume: q(8),
            centroid: zero3(),
            inertia: std::array::from_fn(|i| {
                std::array::from_fn(|j| if i == j { r("16/3") } else { q(0) })
            }),
        };
        let mut model = Model::new();
        let built = equivalent_blocks(&mut model, &target).unwrap();
        assert_eq!(built.achieved.volume, q(8));
        assert!(built.inertia_residual < r("1e-9"));
    }
}
