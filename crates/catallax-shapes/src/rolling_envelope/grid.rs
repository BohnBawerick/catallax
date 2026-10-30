//! A structured display grid of the carrier with a recomputed chord bound.
//!
//! Each cell is one subdivided rational Bezier patch P / W. Against the bilinear
//! map B through its corners, S - B = (P - B W) / W. Both P and B W are polynomials
//! of one degree higher in each direction, so the Bernstein coefficients of the
//! numerator over the smallest weight bound the distance at matching parameters.
//! That bound follows the patch's real curvature, whereas comparing rational
//! controls with their Greville points adds a term from the weights' own slope.
//! The quad's split then adds its twist.

use super::*;

/// Grid over the chart, u-major. Boundary rows are the carrier's own rows.
#[derive(Clone, Debug)]
pub struct EnvelopeGrid {
    /// Grid lines along `u`.
    pub columns: usize,
    /// Grid lines along `v`.
    pub rows: usize,
    /// Enclosed world points, `columns * rows` of them.
    pub points: Vec<[Ball; 3]>,
    /// Two-way bound between the carrier and the grid's triangles, in mm.
    pub chord: Rational,
}

fn l1(a: &[Ball; 3], c: &[Ball; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| sum + magnitude(&(&a[i] - &c[i])))
}

fn along(cell: Vec<Patch>, axis: usize) -> Vec<Patch> {
    cell.into_iter()
        .flat_map(|patch| {
            let half = (&patch.cell[axis][0] + &patch.cell[axis][1]) / q(2);
            let mut low = patch.cell.clone();
            let mut high = patch.cell.clone();
            low[axis][1] = half.clone();
            high[axis][0] = half;
            let [left, right] = split_net(&patch.net, [DEGREE, 2], axis);
            [
                Patch {
                    cell: low,
                    net: left,
                },
                Patch {
                    cell: high,
                    net: right,
                },
            ]
        })
        .collect()
}

impl Patch {
    fn corners(&self) -> Result<[[Ball; 3]; 4], Refusal> {
        let at = |i: usize, j: usize| -> Result<[Ball; 3], Refusal> {
            let h = &self.net[i * 3 + j];
            let divide = |value: &Ball| {
                value
                    .checked_div(&h[3])
                    .ok_or_else(|| unresolved("subdivided envelope weight is not positive"))
            };
            Ok([divide(&h[0])?, divide(&h[1])?, divide(&h[2])?])
        };
        Ok([at(0, 0)?, at(DEGREE, 0)?, at(0, 2)?, at(DEGREE, 2)?])
    }
    fn chord(&self) -> Result<Rational, Refusal> {
        let [a, c, d, e] = self.corners()?;
        let bits = a[0].precision();
        let zero: [Ball; 3] = std::array::from_fn(|_| b(&Rational::zero(), bits));
        // (1 - t) and t raise a Bernstein coefficient i of degree n to i and i + 1.
        let lifted = |degree: usize, step: usize, i: usize| {
            (if step == 0 { degree + 1 - i } else { i + 1 }) as i64
        };
        let corners = [[&a, &d], [&c, &e]];
        let mut numerator = Rational::zero();
        for k in 0..=DEGREE + 1 {
            for l in 0..=3usize {
                let mut sum = zero.clone();
                for (step_u, along_v) in corners.iter().enumerate() {
                    let Some(i) = k.checked_sub(step_u).filter(|i| *i <= DEGREE) else {
                        continue;
                    };
                    for (step_v, corner) in along_v.iter().enumerate() {
                        let Some(j) = l.checked_sub(step_v).filter(|j| *j <= 2) else {
                            continue;
                        };
                        let h = &self.net[i * 3 + j];
                        let factor = b(
                            &Rational::from_ratio(
                                lifted(DEGREE, step_u, i) * lifted(2, step_v, j),
                                (DEGREE as i64 + 1) * 3,
                            )
                            .expect("degree"),
                            bits,
                        );
                        for (m, value) in sum.iter_mut().enumerate() {
                            *value = &*value + &(&(&h[m] - &(&corner[m] * &h[3])) * &factor);
                        }
                    }
                }
                numerator = numerator.max(l1(&sum, &zero));
            }
        }
        let low = self.net.iter().fold(None::<Rational>, |low, h| {
            let w = h[3].lo();
            Some(low.map_or(w.clone(), |value| value.min(w)))
        });
        let low = low.expect("controls");
        if low <= Rational::zero() {
            return Err(unresolved("subdivided envelope weight is not positive"));
        }
        let twist: [Ball; 3] = std::array::from_fn(|k| {
            &(&(&(&a[k] - &c[k]) - &d[k]) + &e[k]) * &b(&Rational::from_ratio(1, 4).unwrap(), bits)
        });
        Ok(numerator / low + l1(&twist, &zero))
    }
}

impl RollingEnvelope {
    /// Grid with `chord` at most `tolerance`. Doubling keeps every carrier knot
    /// on a grid line. A tolerance out of reach after 2^10 cells per span
    /// direction refuses instead of widening the report.
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<EnvelopeGrid, Refusal> {
        let local = self.local_grid(tolerance, max_points, dec)?;
        let frame = LocalFrame::new(&self.data.definition, self.data.net.bits, dec)?;
        Ok(EnvelopeGrid {
            columns: local.columns,
            rows: local.rows,
            points: local.points.iter().map(|p| frame.world(p)).collect(),
            chord: local.chord.clone(),
        })
    }

    fn local_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<Arc<EnvelopeGrid>, Refusal> {
        let mut cache = self.data.net.grid.lock().expect("no panic holds this lock");
        if let Some((key, grid)) = cache.as_ref()
            && key == tolerance
        {
            return Ok(grid.clone());
        }
        // Two cuts across start each meridian below a sixteenth of a turn per cell.
        let mut patches = along(along(Patch::all(&self.data.carrier, &self.data.net), 1), 1);
        let spans = patches.len() / 4;
        let (mut cuts_u, mut cuts_v) = (0u32, 2u32);
        let (patches, chord) = loop {
            dec.check_deadline()?;
            let mut worst = Rational::zero();
            let mut worst_shape = [Rational::zero(), Rational::zero()];
            for patch in &patches {
                dec.check_deadline()?;
                let error = patch.chord()?;
                if error > worst {
                    worst = error;
                    let [a, c, d, _] = patch.corners()?;
                    worst_shape = [l1(&a, &c), l1(&a, &d)];
                }
            }
            if worst <= *tolerance {
                break (patches, worst);
            }
            let cells = (spans << cuts_u) << cuts_v;
            if cuts_u >= 10 || cuts_v >= 12 || cells.saturating_mul(4) > max_points {
                return Err(Refusal::new(
                    "mesh.resource-limit",
                    "envelope grid cannot reach the chord tolerance within its limits",
                ));
            }
            // Refine across the longer side of the worst cell.
            if worst_shape[0] >= worst_shape[1] {
                patches = along(patches, 0);
                cuts_u += 1;
            } else {
                patches = along(patches, 1);
                cuts_v += 1;
            }
        };
        let per_span = 1usize << cuts_v;
        let rows = per_span + 1;
        let columns = patches.len() / per_span + 1;
        let mut points = Vec::with_capacity(columns * rows);
        let mut ordered: Vec<&Patch> = patches.iter().collect();
        ordered.sort_by(|a, c| {
            a.cell[0][0]
                .cmp(&c.cell[0][0])
                .then_with(|| a.cell[1][0].cmp(&c.cell[1][0]))
        });
        for column in 0..columns {
            let source = column.min(columns - 2);
            for row in 0..rows {
                let patch = ordered[source * per_span + row.min(per_span - 1)];
                let [a, c, d, e] = patch.corners()?;
                let corner = match (column > source, row >= per_span) {
                    (false, false) => a,
                    (true, false) => c,
                    (false, true) => d,
                    (true, true) => e,
                };
                points.push(corner);
            }
        }
        let grid = Arc::new(EnvelopeGrid {
            columns,
            rows,
            points,
            chord,
        });
        *cache = Some((tolerance.clone(), grid.clone()));
        Ok(grid)
    }
}
