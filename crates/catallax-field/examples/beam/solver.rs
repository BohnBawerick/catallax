//! Example-private unproved linear elasticity. No production API or export route.
#![allow(clippy::needless_range_loop)] // Tensor indices make the element algebra explicit.
use catallax_numbers::Refusal;
use std::collections::BTreeMap;

type Matrix = [[f64; 24]; 24];
type Strain = [[f64; 24]; 6];

pub(super) struct Estimate {
    pub tip: f64,
    pub reaction: [f64; 3],
    pub reaction_moment: [f64; 3],
    pub iterations: usize,
    pub residual: f64,
    pub stress: Vec<(f64, f64, f64)>,
}

fn fail(reason: &str) -> Refusal {
    Refusal::new("fea.solver", reason)
}
fn constitutive(e: f64, nu: f64) -> [[f64; 6]; 6] {
    let mu = e / (2. * (1. + nu));
    let lambda = e * nu / ((1. + nu) * (1. - 2. * nu));
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i < 3 && j < 3 {
                lambda + if i == j { 2. * mu } else { 0. }
            } else if i == j {
                mu
            } else {
                0.
            }
        })
    })
}
// Node bits are x,y,z. Engineering shear order: xy,yz,xz.
fn strain(h: [f64; 3], p: [f64; 3]) -> Strain {
    let mut b = [[0.; 24]; 6];
    for a in 0..8 {
        let s: [f64; 3] = std::array::from_fn(|d| if a & (1 << d) == 0 { -1. } else { 1. });
        let g: [f64; 3] = std::array::from_fn(|d| {
            let mut v = s[d] / (4. * h[d]);
            for t in 0..3 {
                if t != d {
                    v *= 1. + s[t] * p[t];
                }
            }
            v
        });
        b[0][3 * a] = g[0];
        b[1][3 * a + 1] = g[1];
        b[2][3 * a + 2] = g[2];
        b[3][3 * a] = g[1];
        b[3][3 * a + 1] = g[0];
        b[4][3 * a + 1] = g[2];
        b[4][3 * a + 2] = g[1];
        b[5][3 * a] = g[2];
        b[5][3 * a + 2] = g[0];
    }
    b
}
fn stiffness(h: [f64; 3], c: &[[f64; 6]; 6]) -> Matrix {
    let mut k = [[0.; 24]; 24];
    let q = 1. / 3_f64.sqrt();
    for x in [-q, q] {
        for y in [-q, q] {
            for z in [-q, q] {
                let b = strain(h, [x, y, z]);
                let cb: Strain = std::array::from_fn(|i| {
                    std::array::from_fn(|j| (0..6).map(|r| c[i][r] * b[r][j]).sum())
                });
                for i in 0..24 {
                    for j in 0..24 {
                        k[i][j] += (0..6).map(|r| b[r][i] * cb[r][j]).sum::<f64>()
                            * h.iter().product::<f64>()
                            / 8.;
                    }
                }
            }
        }
    }
    k
}
fn node(n: [usize; 3], i: usize, j: usize, k: usize) -> usize {
    (k * (n[1] + 1) + j) * (n[0] + 1) + i
}
fn dofs(n: [usize; 3], i: usize, j: usize, k: usize) -> [usize; 24] {
    std::array::from_fn(|a| {
        3 * node(
            n,
            i + ((a / 3) & 1),
            j + (((a / 3) >> 1) & 1),
            k + ((a / 3) >> 2),
        ) + a % 3
    })
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
struct Sparse {
    offsets: Vec<usize>,
    entries: Vec<(usize, f64)>,
    diagonal: Vec<f64>,
}
impl Sparse {
    fn apply(&self, x: &[f64]) -> Vec<f64> {
        (0..x.len())
            .map(|r| {
                self.entries[self.offsets[r]..self.offsets[r + 1]]
                    .iter()
                    .map(|(c, v)| v * x[*c])
                    .sum()
            })
            .collect()
    }
}
fn cg(
    a: &Sparse,
    f: &[f64],
    tolerance: f64,
    limit: usize,
) -> Result<(Vec<f64>, usize, f64), Refusal> {
    let mut x = vec![0.; f.len()];
    let mut r = f.to_vec();
    let mut z: Vec<_> = r.iter().zip(&a.diagonal).map(|(r, d)| r / d).collect();
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    if f.iter().all(|value| *value == 0.) {
        return Ok((x, 0, 0.));
    }
    let norm = dot(f, f).sqrt();
    if !norm.is_finite() || norm <= 0. {
        return Err(fail("nonzero load has an unsupported force norm"));
    }
    for iteration in 1..=limit {
        let ap = a.apply(&p);
        let pap = dot(&p, &ap);
        if !pap.is_finite() || pap <= 0. || !rz.is_finite() {
            return Err(fail("nonpositive or nonfinite CG curvature"));
        }
        let alpha = rz / pap;
        for i in 0..x.len() {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        if dot(&r, &r).sqrt() / norm <= tolerance {
            // Recompute the actual residual; do not trust accumulated recurrence error.
            let ax = a.apply(&x);
            let actual: Vec<_> = f.iter().zip(ax).map(|(f, a)| f - a).collect();
            let relative = dot(&actual, &actual).sqrt() / norm;
            if relative.is_finite() && relative <= tolerance {
                return Ok((x, iteration, relative));
            }
            r = actual;
            z = r.iter().zip(&a.diagonal).map(|(r, d)| r / d).collect();
            p = z.clone();
            rz = dot(&r, &z);
            continue;
        }
        z = r.iter().zip(&a.diagonal).map(|(r, d)| r / d).collect();
        let next = dot(&r, &z);
        let beta = next / rz;
        for i in 0..p.len() {
            p[i] = z[i] + beta * p[i];
        }
        rz = next;
    }
    Err(fail(
        "iteration limit reached without verified residual convergence",
    ))
}

pub(super) fn solve(
    n: [usize; 3],
    size: [f64; 3],
    e: f64,
    nu: f64,
    force: f64,
    tolerance: f64,
    limit: usize,
) -> Result<Estimate, Refusal> {
    if n.iter().any(|v| *v == 0 || *v > 80)
        || n.iter().product::<usize>() > 5120
        || size.iter().any(|v| !v.is_finite() || *v <= 0.)
        || !e.is_finite()
        || e <= 0.
        || !nu.is_finite()
        || !(0. ..0.45).contains(&nu)
        || !force.is_finite()
        || !tolerance.is_finite()
        || tolerance <= 0.
        || tolerance > 1e-9
        || limit > 10_000
    {
        return Err(fail("unsupported grid, material, load or solve policy"));
    }
    let nodal_force = force / (n[0] * n[1] * 4) as f64;
    if force != 0. && nodal_force == 0. {
        return Err(fail("nonzero load underflows during force assembly"));
    }
    let h = std::array::from_fn(|d| size[d] / n[d] as f64);
    let c = constitutive(e, nu);
    let ke = stiffness(h, &c);
    let count = 3 * (n[0] + 1) * (n[1] + 1) * (n[2] + 1);
    let fixed = 3 * (n[0] + 1) * (n[1] + 1);
    let mut rows = vec![BTreeMap::<usize, f64>::new(); count];
    for (i, row) in rows.iter_mut().enumerate().take(fixed) {
        row.insert(i, 1.);
    }
    let mut f = vec![0.; count];
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let d = dofs(n, i, j, k);
                for a in 0..24 {
                    if d[a] >= fixed {
                        for b in 0..24 {
                            if d[b] >= fixed {
                                *rows[d[a]].entry(d[b]).or_default() += ke[a][b];
                            }
                        }
                    }
                }
                if k + 1 == n[2] {
                    for a in 4..8 {
                        f[d[3 * a + 1]] += nodal_force;
                    }
                }
            }
        }
    }
    let mut a = Sparse {
        offsets: vec![0],
        entries: Vec::new(),
        diagonal: Vec::new(),
    };
    for (i, row) in rows.into_iter().enumerate() {
        let diagonal = row[&i];
        if !diagonal.is_finite() || diagonal <= 0. {
            return Err(fail("invalid stiffness diagonal"));
        }
        a.diagonal.push(diagonal);
        a.entries.extend(row);
        a.offsets.push(a.entries.len());
    }
    let (u, iterations, residual) = cg(&a, &f, tolerance, limit)?;
    let mut reaction = [0.; 3];
    let mut reaction_moment = [0.; 3];
    // Recover reactions with the original element stiffness before elimination.
    for j in 0..n[1] {
        for i in 0..n[0] {
            let d = dofs(n, i, j, 0);
            for local in 0..4 {
                let r: [f64; 3] = std::array::from_fn(|axis| {
                    (0..24).map(|b| ke[3 * local + axis][b] * u[d[b]]).sum()
                });
                for axis in 0..3 {
                    reaction[axis] += r[axis];
                }
                let x = (i + (local & 1)) as f64 * h[0] - size[0] / 2.;
                let y = (j + (local >> 1)) as f64 * h[1] - size[1] / 2.;
                reaction_moment[0] += y * r[2];
                reaction_moment[1] -= x * r[2];
                reaction_moment[2] += x * r[1] - y * r[0];
            }
        }
    }
    let mut tip = 0.;
    // Same consistent surface weights as the uniform traction, including edges.
    for j in 0..n[1] {
        for i in 0..n[0] {
            for b in 0..4 {
                tip +=
                    u[3 * node(n, i + (b & 1), j + (b >> 1), n[2]) + 1] / (4 * n[0] * n[1]) as f64;
            }
        }
    }
    let b = strain(h, [0.; 3]);
    let mut stress = Vec::new();
    for k in 0..n[2] {
        let d = dofs(n, n[0] / 2, n[1] - 1, k);
        let eps: [f64; 6] = std::array::from_fn(|r| (0..24).map(|a| b[r][a] * u[d[a]]).sum());
        let s: [f64; 6] = std::array::from_fn(|r| (0..6).map(|a| c[r][a] * eps[a]).sum());
        let vm = (0.5 * ((s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2))
            + 3. * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2)))
        .sqrt();
        if !vm.is_finite() {
            return Err(fail("nonfinite stress"));
        }
        stress.push(((k as f64 + 0.5) * h[2], size[1] - h[1] / 2., vm));
    }
    Ok(Estimate {
        tip,
        reaction,
        reaction_moment,
        iterations,
        residual,
        stress,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn affine_patch_and_rigid_motion() {
        let h = [2., 3., 4.];
        let c = constitutive(200_000., 0.3);
        let k = stiffness(h, &c);
        let u: [f64; 24] = std::array::from_fn(|a| {
            let axis = a % 3;
            let coordinate = ((a / 3 >> axis) & 1) as f64 * h[axis];
            coordinate * if axis == 2 { 0.001 } else { -0.0003 }
        });
        // Exact uniaxial affine patch: energy u.K.u = E*strain^2*volume.
        let energy: f64 = (0..24)
            .map(|i| (0..24).map(|j| u[i] * k[i][j] * u[j]).sum::<f64>())
            .sum();
        assert!((energy - 4.8).abs() < 1e-10);
        for i in 0..24 {
            for j in 0..24 {
                assert!((k[i][j] - k[j][i]).abs() < 1e-9);
            }
            for axis in 0..3 {
                assert!((0..8).map(|a| k[i][3 * a + axis]).sum::<f64>().abs() < 1e-9);
            }
        }
    }
    #[test]
    fn nonzero_underflow_refuses_but_exact_zero_load_converges() {
        let run = |force| {
            solve(
                [2, 2, 20],
                [10., 10., 100.],
                200_000.,
                0.3,
                force,
                1e-9,
                10_000,
            )
        };
        // The first witness has nonzero entries but a zero squared norm;
        // the second loses each surface-cell contribution during assembly.
        for force in [-1e-200, -f64::from_bits(1)] {
            let Err(error) = run(force) else {
                panic!("nonzero underflow must not produce a converged estimate");
            };
            assert_eq!(error.code, "fea.solver");
        }
        for force in [0., -0.] {
            let answer = run(force).unwrap();
            assert_eq!(answer.iterations, 0);
            assert_eq!(answer.tip, 0.);
            assert_eq!(answer.residual, 0.);
            assert_eq!(answer.reaction, [0.; 3]);
            assert_eq!(answer.reaction_moment, [0.; 3]);
            assert!(answer.stress.iter().all(|sample| sample.2 == 0.));
        }
    }
    #[test]
    fn independent_workers_are_bit_identical() {
        let run = || {
            solve(
                [2, 2, 20],
                [10., 10., 100.],
                200_000.,
                0.3,
                -10.,
                1e-9,
                10_000,
            )
            .unwrap()
        };
        let expected = run();
        for workers in [1, 2, 4] {
            std::thread::scope(|scope| {
                let handles: Vec<_> = (0..workers).map(|_| scope.spawn(run)).collect();
                for handle in handles {
                    let answer = handle.join().unwrap();
                    assert_eq!(answer.tip.to_bits(), expected.tip.to_bits());
                    assert_eq!(answer.residual.to_bits(), expected.residual.to_bits());
                    assert_eq!(answer.iterations, expected.iterations);
                    for (a, b) in answer.stress.iter().zip(&expected.stress) {
                        assert_eq!(a.2.to_bits(), b.2.to_bits());
                    }
                }
            });
        }
    }
    #[test]
    fn refinement_sign_linearity_and_exhaustion() {
        let run = |m, f, limit| {
            solve(
                [m, m, 10 * m],
                [10., 10., 100.],
                200_000.,
                0.3,
                f,
                1e-9,
                limit,
            )
        };
        let coarse = run(2, -10., 10_000).unwrap();
        let fine = run(4, -10., 10_000).unwrap();
        let opposite = run(2, 20., 10_000).unwrap();
        let finest = run(8, -10., 10_000).unwrap();
        assert!((fine.reaction[1] - 10.).abs() < 1e-6);
        assert!((fine.reaction_moment[0] + 1000.).abs() < 1e-5);
        let theory = -0.02;
        assert!(fine.tip < 0. && (fine.tip / theory - 1.).abs() < (coarse.tip / theory - 1.).abs());
        assert!((finest.tip / theory - 1.).abs() < (fine.tip / theory - 1.).abs());
        assert!((finest.tip / theory - 1.).abs() < 0.02);
        assert!((opposite.tip + 2. * coarse.tip).abs() < 1e-10);
        for (a, b) in coarse.stress.iter().zip(opposite.stress) {
            assert!((b.2 - 2. * a.2).abs() < 1e-8);
        }
        assert!(run(2, -10., 1).is_err());
        assert!(solve([0, 2, 20], [10., 10., 100.], 200_000., 0.3, -10., 1e-9, 100).is_err());
    }
}
