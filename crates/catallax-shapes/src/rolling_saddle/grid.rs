//! A structured display grid of the carrier, sized by second derivatives.
//!
//! A ring has many short spans, so subdividing every patch and measuring it, as
//! [`crate::rational_sheet::sheet_grid`] does, costs seconds. Here each span is
//! bounded once. With `S = P / W` on a patch,
//!
//! ```text
//! S_uu = (P_uu - 2 S_u W_u - S W_uu) / W
//! S_uv = (P_uv - S_u W_v - S_v W_u - S W_uv) / W
//! ```
//!
//! and `S_vv` like `S_uu`, each enclosed from Bernstein differences of the net.
//! A triangulated cell of `du` by `dv` then stays within
//! `du^2 M_uu / 8 + dv^2 M_vv / 8 + du dv M_uv / 4` of the surface, at matching
//! parameters. That bound is looser than measuring each cell, so a grid has more
//! cells, and it is built once per tolerance with no subdivision.

use super::*;

// Cuts of each span across the ring before differencing, so that one patch does
// not hold the whole turn of a section.
const CUTS: usize = 2;
const MAX_CUTS: u32 = 12;

// A grid in the law's frame, with the chart parameter of each column. Columns are
// even only on a uniform carrier.
#[derive(Debug)]
pub(super) struct LocalGrid {
    pub(super) grid: EnvelopeGrid,
    pub(super) parameters: Vec<Rational>,
}

fn hull(values: impl Iterator<Item = Ball>) -> Ball {
    values
        .reduce(|all, value| all.hull(&value))
        .expect("a patch has controls")
}

// One patch's second-derivative bounds, in chart units: `uu`, `vv`, `uv`.
fn turns_of(patch: &Patch) -> Result<[Rational; 3], Refusal> {
    let bits = patch.net[0][3].precision();
    let at = |i: usize, j: usize| &patch.net[i * 3 + j];
    let first = at(0, 0);
    let offset = [0, 1, 2].map(|k| {
        first[k]
            .checked_div(&first[3])
            .map(|value| b(&value.mid(), bits))
    });
    let [Some(mx), Some(my), Some(mz)] = offset else {
        return Err(unresolved("subdivided saddle weight is not positive"));
    };
    let offset = [mx, my, mz];
    let net: Vec<[Ball; 4]> = patch
        .net
        .iter()
        .map(|h| {
            [
                &h[0] - &(&h[3] * &offset[0]),
                &h[1] - &(&h[3] * &offset[1]),
                &h[2] - &(&h[3] * &offset[2]),
                h[3].clone(),
            ]
        })
        .collect();
    let [du, dv] = [0, 1].map(|axis| &patch.cell[axis][1] - &patch.cell[axis][0]);
    let order = patch.degree;
    let degree = q(order as i64);
    let scale = |[a, c]: [u32; 2]| -> Ball {
        let falling = |n: &Rational, order: u32| {
            (0..order).fold(Rational::one(), |p, i| p * (n - &q(i64::from(i))))
        };
        let along_u = falling(&degree, a) / du.pow(a as i32).expect("nonzero patch");
        let along_v = falling(&q(2), c) / dv.pow(c as i32).expect("nonzero patch");
        b(&(along_u * along_v), bits)
    };
    // Bernstein differences of order `(a, c)`, each from the order before it.
    let step = |table: &[[Ball; 4]], (a, c): (usize, usize), axis: usize| -> Vec<[Ball; 4]> {
        let columns = 3 - c;
        let rows = order + 1 - a;
        let (rows_out, columns_out) = if axis == 0 {
            (rows - 1, columns)
        } else {
            (rows, columns - 1)
        };
        let mut next = Vec::with_capacity(rows_out * columns_out);
        for i in 0..rows_out {
            for j in 0..columns_out {
                let far = if axis == 0 {
                    (i + 1) * columns + j
                } else {
                    i * columns + j + 1
                };
                next.push(std::array::from_fn(|k| {
                    &table[far][k] - &table[i * columns + j][k]
                }));
            }
        }
        next
    };
    let first_u = step(&net, (0, 0), 0);
    let first_v = step(&net, (0, 0), 1);
    let tables = [
        ([0, 0], net.clone()),
        ([1, 0], first_u.clone()),
        ([0, 1], first_v),
        ([2, 0], step(&first_u, (1, 0), 0)),
        ([1, 1], step(&first_u, (1, 0), 1)),
        ([0, 2], step(&step(&net, (0, 0), 1), (0, 1), 1)),
    ];
    let difference = |order: [u32; 2], k: usize| -> Ball {
        let (_, table) = tables
            .iter()
            .find(|(key, _)| *key == order)
            .expect("each order has a table");
        hull(table.iter().map(|h| h[k].clone()))
    };
    let derivative = |order: [u32; 2], k: usize| &difference(order, k) * &scale(order);
    let weight = difference([0, 0], 3);
    if weight.lo() <= Rational::zero() {
        return Err(unresolved("subdivided saddle weight is not positive"));
    }
    let over = |value: Ball| {
        value
            .checked_div(&weight)
            .ok_or_else(|| unresolved("subdivided saddle weight is not positive"))
    };
    let [wu, wv, wuu, wvv, wuv] =
        [[1, 0], [0, 1], [2, 0], [0, 2], [1, 1]].map(|order: [u32; 2]| derivative(order, 3));
    let mut sums = [Rational::zero(), Rational::zero(), Rational::zero()];
    for k in 0..3 {
        let s = over(difference([0, 0], k))?;
        let su = over(&derivative([1, 0], k) - &(&s * &wu))?;
        let sv = over(&derivative([0, 1], k) - &(&s * &wv))?;
        let two = b(&q(2), bits);
        let suu = over(&(&derivative([2, 0], k) - &(&(&two * &su) * &wu)) - &(&s * &wuu))?;
        let svv = over(&(&derivative([0, 2], k) - &(&(&two * &sv) * &wv)) - &(&s * &wvv))?;
        let suv = over(&(&(&derivative([1, 1], k) - &(&su * &wv)) - &(&sv * &wu)) - &(&s * &wuv))?;
        for (sum, part) in sums.iter_mut().zip([suu, svv, suv]) {
            *sum = &*sum + &magnitude(&part);
        }
    }
    Ok(sums)
}

impl RollingSaddle {
    // The largest second derivatives of the carrier, kept on the shared net.
    fn turns(&self) -> Result<[Rational; 3], Refusal> {
        let mut cache = self
            .data
            .net
            .turns
            .lock()
            .expect("no panic holds this lock");
        if let Some(turns) = cache.as_ref() {
            return Ok(turns.clone());
        }
        let mut patches = self.patches(self.data.net.bits);
        for _ in 0..CUTS {
            patches = along(patches, 1);
        }
        let mut worst = [Rational::zero(), Rational::zero(), Rational::zero()];
        for patch in &patches {
            for (top, value) in worst.iter_mut().zip(turns_of(patch)?) {
                *top = top.clone().max(value);
            }
        }
        *cache = Some(worst.clone());
        Ok(worst)
    }

    pub(super) fn local_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<LocalGrid, Refusal> {
        let [uu, vv, uv] = self.turns()?;
        let spans = self.span_count();
        // Cuts count from the widest span, and a narrower one takes only the
        // halvings that bring it down to the same cell width.
        let widths: Vec<Rational> = self
            .data
            .nodes
            .windows(2)
            .map(|pair| &pair[1].u - &pair[0].u)
            .collect();
        let span = self.span_widths()[0].clone();
        let halvings = |cuts: u32| -> Vec<u32> {
            let cell = span.mul_pow2(-i64::from(cuts));
            widths
                .iter()
                .map(|width| {
                    let mut own = 0;
                    while width.mul_pow2(-i64::from(own)) > cell {
                        own += 1;
                    }
                    own
                })
                .collect()
        };
        let chord = |cuts: [u32; 2]| {
            let du = span.mul_pow2(-i64::from(cuts[0]));
            let dv = Rational::one().mul_pow2(-i64::from(cuts[1]));
            (&du * &du * uu.clone() + &dv * &dv * vv.clone()).mul_pow2(-3)
                + (&du * &dv * uv.clone()).mul_pow2(-2)
        };
        let (mut cuts_u, mut cuts_v) = (0u32, 2u32);
        loop {
            dec.check_deadline()?;
            if chord([cuts_u, cuts_v]) <= *tolerance {
                break;
            }
            let columns: usize = halvings(cuts_u).iter().map(|own| 1usize << own).sum();
            if cuts_u >= MAX_CUTS
                || cuts_v >= MAX_CUTS
                || (columns + 1).saturating_mul((1 << cuts_v) + 1) > max_points
            {
                return Err(Refusal::new(
                    "mesh.resource-limit",
                    "saddle grid cannot reach the chord tolerance within its limits",
                ));
            }
            if chord([cuts_u + 1, cuts_v]) <= chord([cuts_u, cuts_v + 1]) {
                cuts_u += 1;
            } else {
                cuts_v += 1;
            }
        }
        let cuts = halvings(cuts_u);
        let rows = (1usize << cuts_v) + 1;
        let columns = cuts.iter().map(|own| 1usize << own).sum::<usize>() + 1;
        let mut points = Vec::with_capacity(columns * rows);
        let bits = self.data.net.bits;
        let patches = self.patches(bits);
        let v_step = Rational::one().mul_pow2(-i64::from(cuts_v));
        let basis: Vec<[Ball; 3]> = (0..rows)
            .map(|row| {
                let v = &v_step * &q(row as i64);
                let rest = Rational::one() - v.clone();
                [
                    b(&(&rest * &rest), bits),
                    b(&(&v * &rest * q(2)), bits),
                    b(&(&v * &v), bits),
                ]
            })
            .collect();
        let mut stations = Vec::with_capacity(columns);
        for (index, own) in cuts.iter().enumerate() {
            let per_span = 1usize << own;
            stations.extend((0..per_span).map(|step| (index, step, per_span)));
        }
        stations.push((spans - 1, 1, 1));
        let mut parameters = Vec::with_capacity(columns);
        for (index, step, per_span) in stations {
            dec.check_deadline()?;
            let t = Rational::from_ratio(step as i64, per_span as i64).expect("nonzero");
            parameters.push(&self.data.nodes[index].u + &(&widths[index] * &t));
            let section = casteljau(&patches[index], &t, bits);
            for weights in &basis {
                let sum = |k: usize| {
                    (0..3).fold(b(&Rational::zero(), bits), |total, j| {
                        total + &section[j][k] * &weights[j]
                    })
                };
                let w = sum(3);
                let at = [0, 1, 2].map(|k| sum(k).checked_div(&w));
                let [Some(x), Some(y), Some(z)] = at else {
                    return Err(unresolved("saddle grid weight is not positive"));
                };
                points.push([x, y, z]);
            }
        }
        Ok(LocalGrid {
            grid: EnvelopeGrid {
                columns,
                rows,
                points,
                chord: chord([cuts_u, cuts_v]),
            },
            parameters,
        })
    }
}

// The three homogeneous controls of the section at `t` along one patch.
fn casteljau(patch: &Patch, t: &Rational, bits: u32) -> [[Ball; 4]; 3] {
    let (right, left) = (b(t, bits), b(&(Rational::one() - t.clone()), bits));
    [0, 1, 2].map(|j| {
        let mut row: Vec<[Ball; 4]> = (0..=patch.degree)
            .map(|i| patch.net[i * 3 + j].clone())
            .collect();
        while row.len() > 1 {
            row = row
                .windows(2)
                .map(|pair| {
                    std::array::from_fn(|k| &(&pair[0][k] * &left) + &(&pair[1][k] * &right))
                })
                .collect();
        }
        row.pop().expect("a patch has controls")
    })
}
