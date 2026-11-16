// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: exact rational arithmetic, immutable outputs and no knot snapping.
// Adapted from SMLib NL_CrvTool.cpp, N_CrvInsertKnot and N_CrvSplit, and
// NL_Bezier.cpp, N_BezElevateDegree, at c9c5979f04cea8c38213004d4f04df0dcfb40d05.

use super::*;

impl RationalSplineCurve {
    /// Inserts an interior knot `times` times without changing the parameterized curve.
    /// The active domain is clamped first. Internal multiplicity may not exceed degree.
    pub fn insert_knot(
        &self,
        knot: &Rational,
        times: usize,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let mut curve = self.clone().into_clamped(dec)?;
        let [lo, hi] = curve.domain();
        if knot <= lo || knot >= hi {
            return Err(invalid(
                "inserted knot must be strictly inside the active domain",
            ));
        }
        let mult = curve.knots.iter().filter(|k| *k == knot).count();
        if times > curve.degree.saturating_sub(mult) {
            return Err(invalid("inserted knot multiplicity exceeds degree"));
        }
        let mut net = weighted(&curve.controls, &curve.weights, dec)?;
        for _ in 0..times {
            dec.check_deadline()?;
            let p = curve.degree;
            let k = curve.knots.partition_point(|u| u <= knot) - 1;
            let s = curve.knots.iter().filter(|u| *u == knot).count();
            let mut next = Vec::with_capacity(net.len() + 1);
            next.extend_from_slice(&net[..=k - p]);
            for i in k - p + 1..=k - s {
                dec.check_deadline()?;
                let alpha = (knot - &curve.knots[i]) / (&curve.knots[i + p] - &curve.knots[i]);
                next.push(blend(&net[i], &net[i - 1], &alpha));
            }
            next.extend_from_slice(&net[k - s..]);
            curve.knots.insert(k + 1, knot.clone());
            net = next;
        }
        let (controls, weights) = euclidean(net, dec)?;
        Self::new(curve.degree, curve.knots, controls, weights, dec)
    }

    /// Splits strictly inside the domain. Both pieces retain their original parameters.
    pub fn split(&self, knot: &Rational, dec: &Decider) -> Result<[Self; 2], Refusal> {
        let [lo, hi] = self.domain();
        if knot <= lo || knot >= hi {
            return Err(invalid(
                "split parameter must be strictly inside the active domain",
            ));
        }
        Ok([
            self.extracted([lo, knot], dec)?,
            self.extracted([knot, hi], dec)?,
        ])
    }

    /// Raises to `degree` exactly by Bezier decomposition and elevation.
    /// Original parameters and geometry are unchanged. Internal knots retain a C0
    /// representation, even where the represented curve has higher continuity.
    pub fn raise_degree(&self, degree: usize, dec: &Decider) -> Result<Self, Refusal> {
        if degree < self.degree {
            return Err(invalid("degree raise cannot lower degree"));
        }
        let curve = self.clone().into_clamped(dec)?;
        if degree == curve.degree {
            return Ok(curve);
        }
        let mut breaks = curve.knots.clone();
        breaks.dedup();
        let mut net = Vec::new();
        for pair in breaks.windows(2) {
            let piece = curve.extracted([&pair[0], &pair[1]], dec)?;
            let mut row = weighted(&piece.controls, &piece.weights, dec)?;
            while row.len() <= degree {
                let n = row.len();
                let mut next = Vec::with_capacity(n + 1);
                next.push(row[0].clone());
                for i in 1..n {
                    dec.check_deadline()?;
                    let alpha = Rational::from_ratio(i as i64, n as i64).unwrap();
                    next.push(blend(&row[i - 1], &row[i], &alpha));
                }
                next.push(row[n - 1].clone());
                row = next;
            }
            let discontinuous =
                curve.knots.iter().filter(|u| *u == &pair[0]).count() > curve.degree;
            if net.is_empty() || discontinuous {
                net.extend(row);
            } else {
                net.extend(row.into_iter().skip(1));
            }
        }
        let knots = breaks
            .iter()
            .enumerate()
            .flat_map(|(i, u)| {
                let endpoint = i == 0 || i + 1 == breaks.len();
                let discontinuous = curve.knots.iter().filter(|k| *k == u).count() > curve.degree;
                std::iter::repeat_n(u.clone(), degree + usize::from(endpoint || discontinuous))
            })
            .collect();
        let (controls, weights) = euclidean(net, dec)?;
        Self::new(degree, knots, controls, weights, dec)
    }

    /// Reverses traversal, retaining the active interval.
    pub fn reversed(&self, dec: &Decider) -> Result<Self, Refusal> {
        let curve = self.clone().into_clamped(dec)?;
        let [lo, hi] = curve.domain();
        let sum = lo + hi;
        Self::new(
            curve.degree,
            curve.knots.iter().rev().map(|u| &sum - u).collect(),
            curve.controls.iter().rev().cloned().collect(),
            curve.weights.iter().rev().cloned().collect(),
            dec,
        )
    }

    /// Affinely maps the active domain to a nonempty interval.
    pub fn reparameterized(&self, domain: [&Rational; 2], dec: &Decider) -> Result<Self, Refusal> {
        if domain[0] >= domain[1] {
            return Err(invalid(
                "reparameterization requires an increasing interval",
            ));
        }
        let curve = self.clone().into_clamped(dec)?;
        let [lo, hi] = curve.domain();
        let scale = (domain[1] - domain[0]) / (hi - lo);
        let knots = curve
            .knots
            .iter()
            .map(|u| domain[0] + &(&(u - lo) * &scale))
            .collect();
        Self::new(curve.degree, knots, curve.controls, curve.weights, dec)
    }

    /// Joins directed pieces with exactly matching endpoints. Each piece receives
    /// one unit interval; its own knot spacing is preserved by an affine map.
    /// No gap is snapped and no tangent continuity is claimed.
    pub fn joined(pieces: &[(Self, bool)], dec: &Decider) -> Result<Self, Refusal> {
        let degree = pieces
            .iter()
            .map(|(c, _)| c.degree)
            .max()
            .ok_or_else(|| invalid("join requires at least one curve"))?;
        let mut controls: Vec<RationalPoint> = Vec::new();
        let mut weights: Vec<Rational> = Vec::new();
        let mut knots = Vec::new();
        for (index, (piece, forward)) in pieces.iter().enumerate() {
            dec.check_deadline()?;
            let curve = if *forward {
                piece.clone()
            } else {
                piece.reversed(dec)?
            };
            let curve = curve.raise_degree(degree, dec)?.reparameterized(
                [
                    &Rational::from_i64(index as i64),
                    &Rational::from_i64(index as i64 + 1),
                ],
                dec,
            )?;
            let offset = usize::from(!controls.is_empty());
            let scale = if let Some(last) = controls.last() {
                if last != &curve.controls[0] {
                    return Err(invalid("joined endpoints differ exactly"));
                }
                weights.last().unwrap() / &curve.weights[0]
            } else {
                Rational::one()
            };
            if offset == 0 {
                knots.extend_from_slice(&curve.knots);
            } else {
                knots.pop();
                knots.extend_from_slice(&curve.knots[degree + 1..]);
            }
            controls.extend(curve.controls.into_iter().skip(offset));
            weights.extend(curve.weights.into_iter().skip(offset).map(|w| &w * &scale));
        }
        Self::new(degree, knots, controls, weights, dec)
    }

    /// Evaluates the exact rational definition at a rational parameter.
    pub fn point_at(&self, parameter: &Rational, dec: &Decider) -> Result<RationalPoint, Refusal> {
        let [lo, hi] = self.domain();
        if parameter < lo || parameter > hi {
            return Err(invalid("curve parameter outside domain"));
        }
        let k = if parameter == hi {
            self.controls.len() - 1
        } else {
            self.knots.partition_point(|u| u <= parameter) - 1
        };
        let mut row = weighted(
            &self.controls[k - self.degree..=k],
            &self.weights[k - self.degree..=k],
            dec,
        )?;
        for r in 1..=self.degree {
            for j in (r..=self.degree).rev() {
                dec.check_deadline()?;
                let i = k - self.degree + j;
                let alpha = (parameter - &self.knots[i])
                    / (&self.knots[i + self.degree + 1 - r] - &self.knots[i]);
                row[j] = blend(&row[j], &row[j - 1], &alpha);
            }
        }
        Ok(std::array::from_fn(|i| {
            &row[self.degree][i] / &row[self.degree][3]
        }))
    }
}

pub(super) fn bezier_knots(breaks: &[Rational], degree: usize) -> Vec<Rational> {
    breaks
        .iter()
        .enumerate()
        .flat_map(|(i, u)| {
            std::iter::repeat_n(
                u.clone(),
                degree + usize::from(i == 0 || i + 1 == breaks.len()),
            )
        })
        .collect()
}

impl RationalSplineSurface {
    fn edit_axis(
        &self,
        axis: usize,
        dec: &Decider,
        mut edit: impl FnMut(RationalSplineCurve) -> Result<RationalSplineCurve, Refusal>,
    ) -> Result<Self, Refusal> {
        if axis > 1 {
            return Err(invalid("surface axis must be zero or one"));
        }
        let other = 1 - axis;
        let mut curves = Vec::with_capacity(self.count[other]);
        for j in 0..self.count[other] {
            let indices: Vec<_> = (0..self.count[axis])
                .map(|i| {
                    if axis == 0 {
                        i * self.count[1] + j
                    } else {
                        j * self.count[1] + i
                    }
                })
                .collect();
            let curve = RationalSplineCurve::new(
                self.degrees[axis],
                self.knots[axis].clone(),
                indices.iter().map(|i| self.controls[*i].clone()).collect(),
                indices.iter().map(|i| self.weights[*i].clone()).collect(),
                dec,
            )?;
            curves.push(edit(curve)?);
        }
        let mut degrees = self.degrees;
        let mut knots = self.knots.clone();
        let mut count = self.count;
        degrees[axis] = curves[0].degree;
        knots[axis] = curves[0].knots.clone();
        count[axis] = curves[0].controls.len();
        let mut controls = Vec::with_capacity(count[0] * count[1]);
        let mut weights = Vec::with_capacity(count[0] * count[1]);
        for u in 0..count[0] {
            for v in 0..count[1] {
                let (i, j) = if axis == 0 { (u, v) } else { (v, u) };
                controls.push(curves[j].controls[i].clone());
                weights.push(curves[j].weights[i].clone());
            }
        }
        Self::new(degrees, knots, count, controls, weights, dec)
    }

    /// Inserts a knot in the selected tensor direction, preserving the surface exactly.
    pub fn insert_knot(
        &self,
        axis: usize,
        knot: &Rational,
        times: usize,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        self.edit_axis(axis, dec, |c| c.insert_knot(knot, times, dec))
    }

    /// Raises one tensor degree exactly, with a piecewise Bezier representation.
    pub fn raise_degree(&self, axis: usize, degree: usize, dec: &Decider) -> Result<Self, Refusal> {
        self.edit_axis(axis, dec, |c| c.raise_degree(degree, dec))
    }

    /// Splits the active rectangle without changing either piece's parameters.
    pub fn split(&self, axis: usize, knot: &Rational, dec: &Decider) -> Result<[Self; 2], Refusal> {
        if axis > 1 {
            return Err(invalid("surface axis must be zero or one"));
        }
        let domain = self.domain();
        if knot <= domain[axis][0] || knot >= domain[axis][1] {
            return Err(invalid("surface split must be strictly inside the domain"));
        }
        let mut left = domain;
        let mut right = domain;
        left[axis][1] = knot;
        right[axis][0] = knot;
        Ok([self.extracted(left, dec)?, self.extracted(right, dec)?])
    }

    /// Evaluates the exact rational surface at rational parameters.
    pub fn point_at(
        &self,
        u: &Rational,
        v: &Rational,
        dec: &Decider,
    ) -> Result<RationalPoint, Refusal> {
        let mut row = Vec::with_capacity(self.count[1]);
        let mut row_weights = Vec::with_capacity(self.count[1]);
        // Evaluate the homogeneous u curves, retaining their denominator for v.
        let [lo, hi] = self.domain()[0];
        if u < lo || u > hi {
            return Err(invalid("surface parameter outside domain"));
        }
        let k = if u == hi {
            self.count[0] - 1
        } else {
            self.knots[0].partition_point(|t| t <= u) - 1
        };
        for j in 0..self.count[1] {
            let mut net: Vec<_> = (k - self.degrees[0]..=k)
                .map(|i| {
                    let index = i * self.count[1] + j;
                    let w = &self.weights[index];
                    [
                        &self.controls[index][0] * w,
                        &self.controls[index][1] * w,
                        &self.controls[index][2] * w,
                        w.clone(),
                    ]
                })
                .collect();
            for r in 1..=self.degrees[0] {
                for a in (r..=self.degrees[0]).rev() {
                    dec.check_deadline()?;
                    let i = k - self.degrees[0] + a;
                    let alpha = (u - &self.knots[0][i])
                        / (&self.knots[0][i + self.degrees[0] + 1 - r] - &self.knots[0][i]);
                    net[a] = blend(&net[a], &net[a - 1], &alpha);
                }
            }
            let h = &net[self.degrees[0]];
            row.push(std::array::from_fn(|a| &h[a] / &h[3]));
            row_weights.push(h[3].clone());
        }
        RationalSplineCurve::new(
            self.degrees[1],
            self.knots[1].clone(),
            row,
            row_weights,
            dec,
        )?
        .point_at(v, dec)
    }
}
