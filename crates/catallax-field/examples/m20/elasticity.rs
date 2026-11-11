//! Example-private unproved linear elasticity. No production API or export route.
#![allow(clippy::needless_range_loop)] // Tensor indices make the element algebra explicit.
use catallax_numbers::{Decider, Refusal};
use std::collections::BTreeMap;

type Matrix = [[f64; 24]; 24];
type Strain = [[f64; 24]; 6];

pub(super) struct Estimate {
    pub displacement: Vec<f64>,
    pub products: u64,
    pub applications: u64,
    pub tip: f64,
    pub reaction: [f64; 3],
    pub reaction_moment: [f64; 3],
    pub iterations: usize,
    pub residual: f64,
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
fn cg(
    a: &Sparse,
    f: &[f64],
    tolerance: f64,
    limit: usize,
    workers: usize,
    dec: &mut Decider,
) -> Result<(Vec<f64>, usize, f64, u64), Refusal> {
    std::thread::scope(|scope| {
        let mut senders = Vec::new();
        let mut receivers = Vec::new();
        for w in 0..workers.min(f.len()) {
            let workers = workers.min(f.len());
            let lo = w * f.len() / workers;
            let hi = (w + 1) * f.len() / workers;
            let (send, receive) = std::sync::mpsc::channel::<std::sync::Arc<Vec<f64>>>();
            let (answer, result) = std::sync::mpsc::channel();
            std::thread::Builder::new()
                .spawn_scoped(scope, move || {
                    while let Ok(x) = receive.recv() {
                        let rows: Vec<f64> = (lo..hi)
                            .map(|r| {
                                a.entries[a.offsets[r]..a.offsets[r + 1]]
                                    .iter()
                                    .map(|(c, v)| v * x[*c])
                                    .sum()
                            })
                            .collect();
                        if answer.send(rows).is_err() {
                            break;
                        }
                    }
                })
                .map_err(|_| fail("unable to spawn requested row workers"))?;
            senders.push(send);
            receivers.push(result);
        }
        eprintln!("actual_nonempty_row_workers={}", senders.len());
        let mut applications = 0;
        let mut apply = |x: &[f64]| -> Result<Vec<f64>, Refusal> {
            dec.charge_work(f.len() as u64)?;
            applications += 1;
            let x = std::sync::Arc::new(x.to_vec());
            for send in &senders {
                send.send(x.clone()).map_err(|_| fail("worker stopped"))?;
            }
            let mut rows = Vec::with_capacity(f.len());
            for receive in &receivers {
                rows.extend(receive.recv().map_err(|_| fail("worker stopped"))?);
            }
            Ok(rows)
        };
        let answer = cg_inner(a, f, tolerance, limit, &mut apply);
        drop(senders);
        answer.map(|(u, i, r)| (u, i, r, applications))
    })
}
fn cg_inner(
    a: &Sparse,
    f: &[f64],
    tolerance: f64,
    limit: usize,
    apply: &mut impl FnMut(&[f64]) -> Result<Vec<f64>, Refusal>,
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
        let ap = apply(&p)?;
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
            let ax = apply(&x)?;
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

#[allow(clippy::too_many_arguments)] // Explicit frozen physical inputs plus worker/budget policy.
pub(super) fn solve(
    n: [usize; 3],
    size: [f64; 3],
    e: f64,
    nu: f64,
    force: f64,
    tolerance: f64,
    limit: usize,
    workers: usize,
    dec: &mut Decider,
) -> Result<Estimate, Refusal> {
    if workers == 0
        || n.iter().any(|v| *v == 0 || *v > 80)
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
                dec.charge_work(1)?;
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
    let (u, iterations, residual, applications) = cg(&a, &f, tolerance, limit, workers, dec)?;
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
    if u.iter().any(|v| !v.is_finite()) {
        return Err(fail("nonfinite displacement"));
    }
    Ok(Estimate {
        displacement: u,
        products: applications * a.entries.len() as u64,
        applications,
        tip,
        reaction,
        reaction_moment,
        iterations,
        residual,
    })
}

// Internal axes (y,z,x), tensor order yy,zz,xx,yz,zx,yx.
pub(super) fn tensor(n: [usize; 3], u: &[f64], cell: [usize; 3], local: [f64; 3]) -> [f64; 6] {
    let h = [5. / n[0] as f64, 5. / n[1] as f64, 100. / n[2] as f64];
    let b = strain(h, local);
    let c = constitutive(200000., 0.3);
    let d = dofs(n, cell[0], cell[1], cell[2]);
    let eps: [f64; 6] = std::array::from_fn(|r| (0..24).map(|a| b[r][a] * u[d[a]]).sum());
    let t: [f64; 6] = std::array::from_fn(|r| (0..6).map(|a| c[r][a] * eps[a]).sum());
    [t[2], t[0], t[1], t[5], t[3], t[4]]
}
pub(super) fn vm(s: [f64; 6]) -> f64 {
    (0.5 * ((s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2))
        + 3. * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
        .sqrt()
}
pub(super) fn uz(n: [usize; 3], u: &[f64], cell: [usize; 3], local: [f64; 3]) -> f64 {
    let d = dofs(n, cell[0], cell[1], cell[2]);
    (0..8)
        .map(|a| {
            let weight = (0..3)
                .map(|i| (1. + if a >> i & 1 == 0 { -local[i] } else { local[i] }) / 2.)
                .product::<f64>();
            weight * u[d[3 * a + 1]]
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn affine_patch_and_tensor_with_shear() {
        let n = [2, 2, 20];
        let mut u = vec![0.; 3 * 3 * 3 * 21];
        for k in 0..=20 {
            for j in 0..=2 {
                for i in 0..=2 {
                    // Internal long-axis uniaxial strain plus engineering xy shear.
                    let x = i as f64 * 2.5;
                    let y = j as f64 * 2.5;
                    let z = k as f64 * 5.;
                    let at = 3 * node(n, i, j, k);
                    u[at] = -0.0003 * x + 0.002 * y;
                    u[at + 1] = -0.0003 * y;
                    u[at + 2] = 0.001 * z;
                }
            }
        }
        let s = tensor(n, &u, [1, 1, 10], [0.2, -0.4, 0.8]);
        assert!((s[0] - 200.).abs() < 1e-8);
        assert!(s[1].abs() < 1e-8 && s[2].abs() < 1e-8);
        assert!((s[4] - 200000. / 2.6 * 0.002).abs() < 1e-8);
        assert!((vm(s) - (200_f64.powi(2) + 3. * s[4].powi(2)).sqrt()).abs() < 1e-8);
        let c = constitutive(200000., 0.3);
        let ke = stiffness([2., 3., 4.], &c);
        for row in &ke {
            for axis in 0..3 {
                assert!((0..8).map(|a| row[3 * a + axis]).sum::<f64>().abs() < 1e-9);
            }
        }
    }
    #[test]
    fn zero_and_underflow_are_distinct() {
        for force in [0., -1e-200, -f64::from_bits(1)] {
            let answer = solve(
                [2, 2, 20],
                [5., 5., 100.],
                200000.,
                0.3,
                force,
                1e-9,
                10000,
                2,
                &mut Decider::default(),
            );
            if force == 0. {
                let a = answer.unwrap();
                assert_eq!(a.tip, 0.);
                assert_eq!(a.reaction, [0.; 3]);
            } else {
                assert!(answer.is_err());
            }
        }
    }
}

#[cfg(test)]
mod partition_tests {
    use super::*;
    #[test]
    fn overprovisioned_and_uneven_rows() {
        let a = Sparse {
            offsets: vec![0, 1, 2, 3],
            entries: vec![(0, 2.), (1, 3.), (2, 4.)],
            diagonal: vec![2., 3., 4.],
        };
        let mut expected = None;
        for w in [1, 2, 4, 7] {
            let result = cg(&a, &[2., 6., 12.], 1e-9, 10, w, &mut Decider::default()).unwrap();
            if let Some(ref old) = expected {
                assert_eq!(&result, old);
            } else {
                expected = Some(result);
            }
        }
    }
}
