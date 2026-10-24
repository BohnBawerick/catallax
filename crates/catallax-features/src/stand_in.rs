//! Mass stand-in bodies for field-only features.
//!
//! PRD section 7 and PLAN.md section 15 decide that the B-rep carries a
//! field-only lattice as its outer envelope plus bodies matching the mass,
//! centre of gravity and inertia the field measured. This module builds those
//! bodies from a rational mass target: up to seven exact blocks, symmetric pairs
//! on the rationalised principal axes and one block at the centroid.
//!
//! Volume and centroid match the target exactly. The inertia tensor matches up
//! to a reported residual: the principal frame is a rational quaternion near the
//! exact eigenframe, and block self-inertia enters the mass solve through a
//! fixed-point iteration. Every number in the model is rational, so the bodies
//! are deterministic and STEP-exact.

use catallax_numbers::{Rational, Refusal, Sign};
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
    /// The blocks, centre block first, then the pairs in axis order.
    pub bodies: Vec<BodyId>,
    /// Exact unit-density properties of the blocks together.
    pub achieved: MassTarget,
    /// Largest absolute difference between achieved and target inertia components.
    pub inertia_residual: Rational,
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

/// Deterministic cube root by Newton steps on doubles. No platform `cbrt`.
fn cbrt(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let mut y = f64::from_bits(x.to_bits() / 3 + 0x2a9f_7893_0000_0000);
    for _ in 0..40 {
        y = (2.0 * y + x / (y * y)) / 3.0;
    }
    y
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

/// Blocks follow the target's aspect ratio. Two dyadic sides fix the shape;
/// the third side closes its volume exactly.
fn sides_for(
    mass: &Rational,
    along: usize,
    aspect: &[Rational; 3],
) -> Result<[Rational; 3], Refusal> {
    let aspect_volume = &(&aspect[0] * &aspect[1]) * &aspect[2];
    let scale = dyadic(cbrt(mass.to_f64_nearest() / aspect_volume.to_f64_nearest()));
    if scale.is_zero() {
        return Err(Refusal::new(
            "stand-in.numeric-range",
            "surrogate dimensions exceed the finite numeric solve range",
        ));
    }
    let mut sides = std::array::from_fn(|i| {
        if i == along {
            q(0)
        } else {
            &scale * &aspect[i]
        }
    });
    let product = &sides[(along + 1) % 3] * &sides[(along + 2) % 3];
    sides[along] = mass
        .checked_div(&product)
        .expect("positive surrogate sides");
    Ok(sides)
}

/// Builds the stand-in blocks into `model` and returns their exact properties.
pub fn equivalent_blocks(model: &mut Model, target: &MassTarget) -> Result<StandInBodies, Refusal> {
    if target.volume.sign() != Sign::Positive {
        return Err(Refusal::new(
            "stand-in.empty",
            "a mass stand-in needs a positive volume",
        ));
    }
    let f: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| target.inertia[i][j].to_f64_nearest()));
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
    // Target second moments in the principal frame: P_x = integral of x^2.
    let local = mat_mul(&rt, &mat_mul(&target.inertia, &r));
    let d = [
        local[0][0].clone(),
        local[1][1].clone(),
        local[2][2].clone(),
    ];
    let p: [Rational; 3] = std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        (&(&d[j] + &d[k]) - &d[i]).mul_pow2(-1)
    });
    let v = &target.volume;
    let aspect: [Rational; 3] =
        std::array::from_fn(|i| dyadic((p[i].to_f64_nearest() / v.to_f64_nearest()).sqrt()));
    if aspect.iter().any(Rational::is_zero) {
        return Err(Refusal::new(
            "stand-in.numeric-range",
            "a positive spatial mass target is required in every principal direction",
        ));
    }
    // Satellite distances: dyadic numbers near the gyration radius, grown until
    // the masses stay positive and the blocks are apart.
    let mut distances: [Rational; 3] = std::array::from_fn(|i| {
        if p[i].sign() == Sign::Positive {
            dyadic((6.0 * p[i].to_f64_nearest() / v.to_f64_nearest()).sqrt())
        } else {
            q(1)
        }
    });
    let mut blocks = Vec::new();
    for _attempt in 0..64 {
        let mut masses = [q(0), q(0), q(0)];
        let mut centre_mass = v.clone();
        let mut previous_masses: Option<[Rational; 3]> = None;
        let mut candidate: Vec<Block> = Vec::new();
        // Quantise each proposed satellite mass, then close the centre mass
        // exactly. This bounds denominator growth without pretending that the
        // resulting inertia residual is zero.
        for _round in 0..64 {
            // Self inertia of the current blocks, expressed as second moments.
            let mut own = zero33();
            for block in &candidate {
                let shape = Block {
                    offset: zero3(),
                    sides: block.sides.clone(),
                };
                add33(&mut own, &shape.inertia());
            }
            let own_p: [Rational; 3] = std::array::from_fn(|i| {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                (&(&own[j][j] + &own[k][k]) - &own[i][i]).mul_pow2(-1)
            });
            let next: [Rational; 3] = std::array::from_fn(|i| {
                let remaining = &p[i] - &own_p[i];
                if remaining.sign() == Sign::Positive {
                    remaining
                        .checked_div(&(&distances[i] * &distances[i]).mul_pow2(1))
                        .expect("positive distance")
                } else {
                    q(0)
                }
            });
            let mut next = next.map(|mass| dyadic(mass.to_f64_nearest()));
            // The compact-box limit has a slowly contracting fixed point.
            // A safeguarded Aitken step accelerates a decreasing sequence
            // without increasing the iteration or arithmetic budgets.
            if let Some(previous) = &previous_masses {
                for axis in 0..3 {
                    if next[axis] <= masses[axis] && masses[axis] <= previous[axis] {
                        let delta = &next[axis] - &masses[axis];
                        let denominator = &delta - &(&masses[axis] - &previous[axis]);
                        if let Some(correction) = (&delta * &delta).checked_div(&denominator) {
                            let accelerated = &next[axis] - &correction;
                            if accelerated.sign() != Sign::Negative && accelerated < next[axis] {
                                next[axis] = dyadic(accelerated.to_f64_nearest());
                            }
                        }
                    }
                }
            }
            let next_centre = v - &(&(&next[0] + &next[1]) + &next[2]).mul_pow2(1);
            let settled = next == masses && next_centre == centre_mass && !candidate.is_empty();
            previous_masses = Some(masses.clone());
            masses = next;
            centre_mass = next_centre;
            if centre_mass.sign() != Sign::Positive {
                break;
            }
            candidate = vec![Block {
                offset: zero3(),
                sides: sides_for(&centre_mass, 0, &aspect)?,
            }];
            for axis in 0..3 {
                if masses[axis].is_zero() {
                    continue;
                }
                for sign in [-1, 1] {
                    let mut offset = zero3();
                    offset[axis] = &distances[axis] * &q(sign);
                    candidate.push(Block {
                        offset,
                        sides: sides_for(&masses[axis], axis, &aspect)?,
                    });
                }
            }
            if settled {
                break;
            }
        }
        let apart = centre_mass.sign() == Sign::Positive
            && candidate.iter().enumerate().all(|(i, a)| {
                candidate.iter().skip(i + 1).all(|b| {
                    (0..3).any(|axis| {
                        let gap = (&a.offset[axis] - &b.offset[axis]).abs();
                        gap > (&a.sides[axis] + &b.sides[axis]).mul_pow2(-1)
                    })
                })
            });
        if apart {
            blocks = candidate;
            break;
        }
        for distance in &mut distances {
            *distance = distance.mul_pow2(1);
        }
    }
    if blocks.is_empty() {
        return Err(Refusal::new(
            "stand-in.unsolved",
            "no block arrangement reproduced the mass target",
        ));
    }
    // Exact achieved properties in model axes.
    let mut volume = q(0);
    let mut first = zero3();
    let mut local_inertia = zero33();
    for block in &blocks {
        let m = block.mass();
        for (i, value) in first.iter_mut().enumerate() {
            *value = &*value + &(&m * &block.offset[i]);
        }
        volume = volume + m;
        add33(&mut local_inertia, &block.inertia());
    }
    let shift = mat_vec(&r, &first);
    let centroid: [Rational; 3] = std::array::from_fn(|i| {
        &target.centroid[i]
            + &shift[i]
                .checked_div(&volume)
                .expect("positive stand-in volume")
    });
    let inertia = mat_mul(&r, &mat_mul(&local_inertia, &rt));
    let mut residual = q(0);
    for (i, row) in inertia.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            residual = residual.max((value - &target.inertia[i][j]).abs());
        }
    }
    let mut bodies = Vec::with_capacity(blocks.len());
    for block in &blocks {
        let half: [Rational; 3] = std::array::from_fn(|i| block.sides[i].mul_pow2(-1));
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
            volume,
            centroid,
            inertia,
        },
        inertia_residual: residual,
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
