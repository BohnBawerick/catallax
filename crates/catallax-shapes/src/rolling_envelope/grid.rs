//! A structured display grid of the carrier with a recomputed chord bound.
//!
//! Each cell is one subdivided Bezier patch. Its controls lie within `e` of the
//! bilinear map through the cell's corners at the matching Greville points, so
//! the rational surface and the cell's two triangles stay within
//! `e + twist + (max w / min w - 1) * diameter` of each other in both directions:
//! `twist` covers the quad's split, and the weight term the difference between
//! rational and plain Bernstein combinations of the same points.

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
        let points = self.cartesian()?;
        let at = |i: usize, j: usize| points[i * 3 + j].clone();
        Ok([at(0, 0), at(DEGREE, 0), at(0, 2), at(DEGREE, 2)])
    }
    fn chord(&self) -> Result<Rational, Refusal> {
        let points = self.cartesian()?;
        let [a, c, d, e] = self.corners()?;
        let bits = a[0].precision();
        let mut error = Rational::zero();
        for i in 0..=DEGREE {
            for j in 0..=2 {
                let x = b(
                    &Rational::from_ratio(i as i64, DEGREE as i64).expect("degree"),
                    bits,
                );
                let y = b(&Rational::from_ratio(j as i64, 2).expect("two"), bits);
                let rest = |t: &Ball| &b(&Rational::one(), bits) - t;
                let bilinear: [Ball; 3] = std::array::from_fn(|k| {
                    &(&(&a[k] * &(&rest(&x) * &rest(&y))) + &(&c[k] * &(&x * &rest(&y))))
                        + &(&(&d[k] * &(&rest(&x) * &y)) + &(&e[k] * &(&x * &y)))
                });
                error = error.max(l1(&points[i * 3 + j], &bilinear));
            }
        }
        let twist: [Ball; 3] = std::array::from_fn(|k| {
            &(&(&(&a[k] - &c[k]) - &d[k]) + &e[k]) * &b(&Rational::from_ratio(1, 4).unwrap(), bits)
        });
        let zero: [Ball; 3] = std::array::from_fn(|_| b(&Rational::zero(), bits));
        let diameter = points[1..]
            .iter()
            .fold(Rational::zero(), |worst, p| worst.max(l1(p, &points[0])));
        let weights = self.net.iter().map(|h| &h[3]);
        let (mut low, mut high) = (None::<Rational>, None::<Rational>);
        for w in weights {
            low = Some(low.map_or(w.lo(), |value| value.min(w.lo())));
            high = Some(high.map_or(w.hi(), |value| value.max(w.hi())));
        }
        let (low, high) = (low.expect("controls"), high.expect("controls"));
        if low <= Rational::zero() {
            return Err(unresolved("subdivided envelope weight is not positive"));
        }
        Ok(error + l1(&twist, &zero) + (high / low - Rational::one()) * diameter)
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
