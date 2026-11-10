//! The closed-mesh leaf of RECIPE.md section 4.19: the region whose integer
//! winding number is one, placed by a rational rigid frame.
//!
//! The mesh was proved closed, oriented and free of crossings when it was read.
//! A point is classified by counting signed crossings along a symbolically
//! perturbed ray, an exact integer with no threshold. The value is the true
//! signed Euclidean distance to the surface, from the exact squared distance to
//! the nearest triangle. The leaf neither builds nor queries B-rep topology.
//!
//! Volume, centroid and inertia of the leaf are exact rationals, summed here
//! over the signed tetrahedra the triangles span with the origin.

use std::sync::Arc;

use catallax_import::mesh::ClosedMesh;
use catallax_numbers::{Ball, Decider, Integer, Rational, Refusal};
use catallax_shapes::{Isometry, RPoint3, Rotation};

use crate::Location;
use crate::metric::{PRECISION, Scalar};

#[derive(Clone)]
pub(crate) struct MeshLeaf {
    mesh: Arc<ClosedMesh>,
    /// World point = rotation * local point + translation.
    rotation: [[Rational; 3]; 3],
    translation: [Rational; 3],
    inverse_rotation: [[Rational; 3]; 3],
    fast_inverse: [[f64; 3]; 3],
    fast_translation: [f64; 3],
}

impl std::fmt::Debug for MeshLeaf {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("MeshLeaf")
            .field("mesh", &self.mesh)
            .field("translation", &self.translation)
            .finish()
    }
}

fn frame_error() -> Refusal {
    Refusal::new(
        "import.mesh-pose",
        "a mesh pose must be an exact rational rigid motion",
    )
}

fn dot(a: &[Rational; 3], b: &[Rational; 3]) -> Rational {
    &(&a[0] * &b[0]) + &(&(&a[1] * &b[1]) + &(&a[2] * &b[2]))
}

fn apply(m: &[[Rational; 3]; 3], v: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| dot(&m[i], v))
}

fn transpose(m: &[[Rational; 3]; 3]) -> [[Rational; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| m[j][i].clone()))
}

fn fast(m: &[[Rational; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| m[i][j].to_f64_nearest()))
}

impl MeshLeaf {
    pub(crate) fn new(
        mesh: Arc<ClosedMesh>,
        pose: &catallax_recipe::Pose,
    ) -> Result<Self, Refusal> {
        let [w, x, y, z] = pose.rotation.each_ref().map(|n| n.value().clone());
        let m = Rotation::from_quaternion(w, x, y, z)?;
        let m = m.rational_matrix().ok_or_else(frame_error)?;
        let rotation: [[Rational; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| m[i][j].clone()));
        let translation = pose.translation.0.each_ref().map(|n| n.value().clone());
        Ok(Self::framed(mesh, rotation, translation))
    }

    fn framed(
        mesh: Arc<ClosedMesh>,
        rotation: [[Rational; 3]; 3],
        translation: [Rational; 3],
    ) -> Self {
        let inverse_rotation = transpose(&rotation);
        Self {
            mesh,
            fast_inverse: fast(&inverse_rotation),
            fast_translation: translation.each_ref().map(Rational::to_f64_nearest),
            rotation,
            translation,
            inverse_rotation,
        }
    }

    fn local(&self, p: &[Rational; 3]) -> [Rational; 3] {
        let shifted: [Rational; 3] = std::array::from_fn(|i| &p[i] - &self.translation[i]);
        apply(&self.inverse_rotation, &shifted)
    }

    /// Composes a rational rigid motion, reflections included, in front of
    /// the frame. Classification maps the query back into the mesh's own
    /// coordinates, so a reflected copy keeps its winding numbers.
    pub(crate) fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        let m = map.rational_matrix().ok_or_else(frame_error)?;
        let m: [[Rational; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| m[i][j].clone()));
        let t = map.translation.as_rational().ok_or_else(frame_error)?;
        let t = [t.x, t.y, t.z];
        let rotation = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(Rational::zero(), |s, k| s + &m[i][k] * &self.rotation[k][j])
            })
        });
        let shifted = apply(&m, &self.translation);
        let translation = std::array::from_fn(|i| &shifted[i] + &t[i]);
        *self = Self::framed(self.mesh.clone(), rotation, translation);
        Ok(())
    }

    /// Exact classification of a rational world point by its winding number.
    /// Work is charged to `dec`.
    pub(crate) fn location(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        Ok(match self.mesh.winding(&self.local(p), dec)? {
            None => Location::Boundary,
            Some(0) => Location::Outside,
            Some(_) => Location::Inside,
        })
    }

    /// Enclosure of the true signed distance, negative inside.
    pub(crate) fn proven_distance(
        &self,
        p: &RPoint3,
        location: Location,
        dec: &mut Decider,
    ) -> Result<Ball, Refusal> {
        if location == Location::Boundary {
            return Ok(Ball::from_rational(&Rational::zero(), PRECISION));
        }
        let local = self.local(&[p.x.clone(), p.y.clone(), p.z.clone()]);
        let squared = self
            .mesh
            .squared_distance(&local, dec)?
            .expect("a compiled mesh leaf has triangles");
        let distance = Ball::from_rational(&squared, PRECISION)
            .sqrt()
            .expect("non-negative squared distance");
        Ok(if location == Location::Inside {
            -distance
        } else {
            distance
        })
    }

    /// Tentative signed distance from binary64 searches only: a nearest
    /// triangle for the magnitude and a binary64 winding count for the sign.
    pub(crate) fn fast_distance(&self, p: [f64; 3]) -> f64 {
        let shifted: [f64; 3] = std::array::from_fn(|i| p[i] - self.fast_translation[i]);
        let local: [f64; 3] = std::array::from_fn(|i| {
            (self.fast_inverse[i][0] * shifted[0] + self.fast_inverse[i][1] * shifted[1])
                + self.fast_inverse[i][2] * shifted[2]
        });
        let magnitude = Scalar::sqrt(&self.mesh.fast_squared_distance(local));
        if self.mesh.fast_winding(local) == 0 {
            magnitude
        } else {
            -magnitude
        }
    }

    /// Exact world bounds of the eight corners of the mesh's own bounding box,
    /// an enclosure of the placed mesh.
    pub(crate) fn bounds(&self) -> [[Rational; 2]; 3] {
        let local = self
            .mesh
            .bounds()
            .expect("a compiled mesh leaf has triangles");
        let mut out: Option<[[Rational; 2]; 3]> = None;
        for corner in 0..8 {
            let l: [Rational; 3] = std::array::from_fn(|i| local[i][corner >> i & 1].clone());
            let rotated = apply(&self.rotation, &l);
            let w: [Rational; 3] = std::array::from_fn(|i| &rotated[i] + &self.translation[i]);
            out = Some(match out {
                None => w.map(|x| [x.clone(), x]),
                Some(b) => std::array::from_fn(|i| {
                    [
                        b[i][0].clone().min(w[i].clone()),
                        b[i][1].clone().max(w[i].clone()),
                    ]
                }),
            });
        }
        out.expect("eight corners")
    }

    /// Exact unit-density moments in world coordinates, in the order volume,
    /// first moments x, y, z, then second moments xx, yy, zz, xy, xz, yz.
    /// Work is charged to `dec`.
    pub(crate) fn moments(&self, dec: &mut Decider) -> Result<[Rational; 10], Refusal> {
        let local = local_moments(&self.mesh, dec)?;
        let volume = local[0].clone();
        let first: [Rational; 3] = std::array::from_fn(|i| local[1 + i].clone());
        let second: [[Rational; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| local[4 + crate::mass::pair_index(i, j)].clone())
        });
        // World x = R l + t: M' = R M + V t and S' = R S R^T + (R M) t^T + t (R M)^T + V t t^T.
        let t = &self.translation;
        let rm = apply(&self.rotation, &first);
        let rs: [[Rational; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(Rational::zero(), |s, k| {
                    s + &self.rotation[i][k] * &second[k][j]
                })
            })
        });
        let rsr: [[Rational; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(Rational::zero(), |s, k| {
                    s + &rs[i][k] * &self.rotation[j][k]
                })
            })
        });
        let mut out: [Rational; 10] = Default::default();
        out[0] = volume.clone();
        for i in 0..3 {
            out[1 + i] = &rm[i] + &(&volume * &t[i]);
        }
        for (index, (i, j)) in crate::mass::PAIRS.iter().copied().enumerate() {
            out[4 + index] = &(&(&rsr[i][j] + &(&rm[i] * &t[j])) + &(&t[i] * &rm[j]))
                + &(&volume * &(&t[i] * &t[j]));
        }
        Ok(out)
    }
}

/// Exact moments in the mesh's scaled coordinates. Coordinates are put over
/// one common denominator so the sums run in integers; the sums are divided
/// once at the end.
fn local_moments(mesh: &ClosedMesh, dec: &mut Decider) -> Result<[Rational; 10], Refusal> {
    let count = mesh.vertex_count() as u32;
    // The common denominator can grow with every distinct denominator, for
    // example under OBJ weights, so each step is charged by its operand sizes.
    let limbs = |x: &Integer| x.bit_length() / 64 + 1;
    let mut common = Integer::one();
    for v in 0..count {
        for x in mesh.file_vertex(v) {
            let d = x.denominator();
            dec.charge_work(4u64.saturating_mul(limbs(&common).saturating_mul(limbs(d))))?;
            let g = common.gcd(d);
            let (factor, _) = d.div_rem_floor(&g).expect("nonzero gcd");
            common = &common * &factor;
        }
    }
    let mut integers: Vec<[Integer; 3]> = Vec::with_capacity(count as usize);
    for v in 0..count {
        let point = mesh.file_vertex(v);
        let size = point
            .iter()
            .map(|x| limbs(x.numerator()).max(limbs(x.denominator())))
            .max()
            .unwrap_or(1);
        dec.charge_work(6u64.saturating_mul(limbs(&common).saturating_mul(size)))?;
        integers.push(point.each_ref().map(|x| {
            let (factor, _) = common
                .div_rem_floor(x.denominator())
                .expect("nonzero denominator");
            x.numerator() * &factor
        }));
    }
    // Each triangle costs about sixty integer operations on numbers as large
    // as the largest scaled coordinate; charges go out in batches.
    let size = integers
        .iter()
        .flatten()
        .map(|x| x.bit_length() / 64 + 1)
        .max()
        .unwrap_or(1);
    let per_triangle = 60u64.saturating_mul(size.saturating_mul(size));
    let zero = Integer::zero;
    let mut volume6 = zero();
    let mut first24: [Integer; 3] = std::array::from_fn(|_| zero());
    let mut second120: [Integer; 6] = std::array::from_fn(|_| zero());
    for (n, tri) in mesh.triangles().iter().enumerate() {
        if n % 1024 == 0 {
            let batch = (mesh.triangles().len() - n).min(1024) as u64;
            dec.charge_work(batch.saturating_mul(per_triangle))?;
        }
        let [a, b, c] = tri.map(|v| &integers[v as usize]);
        let det = &(&(&a[0] * &(&(&b[1] * &c[2]) - &(&b[2] * &c[1])))
            - &(&a[1] * &(&(&b[0] * &c[2]) - &(&b[2] * &c[0]))))
            + &(&a[2] * &(&(&b[0] * &c[1]) - &(&b[1] * &c[0])));
        let sum: [Integer; 3] = std::array::from_fn(|i| &(&a[i] + &b[i]) + &c[i]);
        volume6 = &volume6 + &det;
        for i in 0..3 {
            first24[i] = &first24[i] + &(&det * &sum[i]);
        }
        for (index, (i, j)) in crate::mass::PAIRS.iter().copied().enumerate() {
            let products =
                &(&(&(&a[i] * &a[j]) + &(&b[i] * &b[j])) + &(&c[i] * &c[j])) + &(&sum[i] * &sum[j]);
            second120[index] = &second120[index] + &(&det * &products);
        }
    }
    let scale = mesh.scale();
    let unit = Rational::from_integer(common);
    let ratio = scale
        .checked_div(&unit)
        .expect("nonzero common denominator");
    let power = |n: i32| ratio.pow(n).expect("nonzero ratio");
    let divide = |sum: Integer, by: i64, n: i32| {
        (Rational::from_integer(sum) * power(n))
            .checked_div(&Rational::from_i64(by))
            .expect("nonzero constant")
    };
    let mut out: [Rational; 10] = Default::default();
    out[0] = divide(volume6, 6, 3);
    for (i, s) in first24.into_iter().enumerate() {
        out[1 + i] = divide(s, 24, 4);
    }
    for (index, s) in second120.into_iter().enumerate() {
        out[4 + index] = divide(s, 120, 5);
    }
    Ok(out)
}
