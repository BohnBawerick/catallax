// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: exact Rational homogeneous arithmetic, validated finite
// domains, typed refusals and caller-owned work/cancellation checkpoints.
// Adapted from rust-brep crates/kernel/src/curve/construction.rs,
// N_CrvExtractCrvSeg, commit 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab.
//! Exact finite-domain extraction of exchange spline definitions. These data
//! carry no regularity, trim, incidence or solid certificate.

use crate::{BSplineCurve, BSplineSurface, Point3};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};

/// Euclidean millimetre coordinates. Weights are always stored separately.
pub type RationalPoint = [Rational; 3];
// Weighted coordinates, never Euclidean coordinates plus a weight.
type Homogeneous = [Rational; 4];

fn invalid(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.spline-exchange-data", detail)
}

fn axis(degree: usize, knots: &[Rational], count: usize, dec: &Decider) -> Result<(), Refusal> {
    dec.check_deadline()?;
    if degree == 0
        || degree >= count
        || count.checked_add(degree).and_then(|n| n.checked_add(1)) != Some(knots.len())
    {
        return Err(invalid(
            "invalid degree, control count or expanded knot count",
        ));
    }
    let mut multiplicity = 1;
    for pair in knots.windows(2) {
        dec.check_deadline()?;
        if pair[0] > pair[1] {
            return Err(invalid("decreasing knots"));
        }
        multiplicity = if pair[0] == pair[1] {
            multiplicity + 1
        } else {
            1
        };
        if multiplicity > degree + 1 {
            return Err(invalid("knot multiplicity exceeds degree plus one"));
        }
    }
    if knots[degree] >= knots[count] {
        return Err(invalid("empty active knot domain"));
    }
    Ok(())
}

fn weights(weights: &[Rational], count: usize, dec: &Decider) -> Result<(), Refusal> {
    if weights.len() != count {
        return Err(invalid("control and weight counts differ"));
    }
    for weight in weights {
        dec.check_deadline()?;
        if weight.sign() != Sign::Positive {
            return Err(invalid("weights must be strictly positive"));
        }
    }
    Ok(())
}

fn weighted(
    points: &[RationalPoint],
    weights: &[Rational],
    dec: &Decider,
) -> Result<Vec<Homogeneous>, Refusal> {
    points
        .iter()
        .zip(weights)
        .map(|(p, w)| {
            dec.check_deadline()?;
            Ok([&p[0] * w, &p[1] * w, &p[2] * w, w.clone()])
        })
        .collect()
}

fn euclidean(
    net: Vec<Homogeneous>,
    dec: &Decider,
) -> Result<(Vec<RationalPoint>, Vec<Rational>), Refusal> {
    let mut points = Vec::with_capacity(net.len());
    let mut weights = Vec::with_capacity(net.len());
    for h in net {
        dec.check_deadline()?;
        if h[3].sign() != Sign::Positive {
            return Err(invalid("extraction produced a nonpositive weight"));
        }
        points.push(std::array::from_fn(|i| &h[i] / &h[3]));
        weights.push(h[3].clone());
    }
    Ok((points, weights))
}

fn blend(a: &Homogeneous, b: &Homogeneous, alpha: &Rational) -> Homogeneous {
    let beta = &Rational::one() - alpha;
    std::array::from_fn(|i| alpha * &a[i] + &beta * &b[i])
}

// Exact version of the donor's two endpoint insertion loops. Equality is exact:
// no legacy parameter tolerance may absorb a small but nonzero knot interval.
fn extract(
    degree: usize,
    knots: &[Rational],
    points: &[Homogeneous],
    domain: [&Rational; 2],
    dec: &Decider,
) -> Result<(Vec<Rational>, Vec<Homogeneous>), Refusal> {
    let [left, right] = domain;
    if left >= right || left < &knots[degree] || right > &knots[points.len()] {
        return Err(invalid("extraction domain is outside the active interval"));
    }
    let span = |t: &Rational| knots.partition_point(|k| k <= t) - 1;
    let spl = span(left);
    let spr = span(right);
    let ml = knots.iter().filter(|k| *k == left).count();
    let mr = knots.iter().filter(|k| *k == right).count();
    let start = spl - degree;
    let end = spr
        .checked_sub(mr)
        .ok_or_else(|| invalid("invalid extraction end"))?;
    if start > end || end >= points.len() {
        return Err(invalid("invalid extraction control range"));
    }
    let mut out = points[start..=end].to_vec();
    let n = out.len() - 1;
    for i in 1..=degree.saturating_sub(ml) {
        for j in 0..=degree - i - ml {
            dec.check_deadline()?;
            let origin = &knots[start + i + j];
            let denominator = &knots[spl + j + 1] - origin;
            if denominator.is_zero() {
                return Err(invalid("zero extraction denominator"));
            }
            let alpha = (left - origin) / denominator;
            out[j] = blend(&out[j + 1], &out[j], &alpha);
        }
    }
    let lr = spr - degree;
    // Use signed arithmetic: an unclamped right endpoint can have multiplicity 0.
    let lk = n as isize - degree as isize + mr as isize;
    for i in 1..=degree.saturating_sub(mr) {
        for j in (0..=degree - i - mr).rev() {
            dec.check_deadline()?;
            let k = usize::try_from(lk + i as isize + j as isize)
                .map_err(|_| invalid("invalid extraction index"))?;
            let origin = knots[lr + i + j].clone().max(left.clone());
            let denominator = &knots[spr + j + 1] - &origin;
            if denominator.is_zero() {
                return Err(invalid("zero extraction denominator"));
            }
            let alpha = (right - &origin) / denominator;
            out[k] = blend(&out[k], &out[k - 1], &alpha);
        }
    }
    let mut out_knots = Vec::with_capacity(out.len() + degree + 1);
    out_knots.extend(std::iter::repeat_n(left.clone(), degree + 1));
    out_knots.extend(knots[spl + 1..=spr - mr].iter().cloned());
    out_knots.extend(std::iter::repeat_n(right.clone(), degree + 1));
    Ok((out_knots, out))
}

/// Exact source curve, including nonclamped knots. Decimal definitions never
/// pass through binary64. Closure and periodicity require separate evidence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RationalSplineCurve {
    degree: usize,
    knots: Vec<Rational>,
    controls: Vec<RationalPoint>,
    weights: Vec<Rational>,
}

impl RationalSplineCurve {
    /// Validates the source storage, active domain and positive weights.
    pub fn new(
        degree: usize,
        knots: Vec<Rational>,
        controls: Vec<RationalPoint>,
        weights_: Vec<Rational>,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        axis(degree, &knots, controls.len(), dec)?;
        weights(&weights_, controls.len(), dec)?;
        Ok(Self {
            degree,
            knots,
            controls,
            weights: weights_,
        })
    }
    /// Degree.
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// Original expanded knots.
    pub fn knots(&self) -> &[Rational] {
        &self.knots
    }
    /// Original Euclidean controls, in millimetres.
    pub fn controls(&self) -> &[RationalPoint] {
        &self.controls
    }
    /// Original positive, dimensionless weights.
    pub fn weights(&self) -> &[Rational] {
        &self.weights
    }
    /// Active finite parameter interval.
    pub fn domain(&self) -> [&Rational; 2] {
        [&self.knots[self.degree], &self.knots[self.controls.len()]]
    }
    /// Clamps an exact subinterval by homogeneous endpoint insertion.
    pub fn extracted(&self, domain: [&Rational; 2], dec: &Decider) -> Result<Self, Refusal> {
        let (knots, net) = extract(
            self.degree,
            &self.knots,
            &weighted(&self.controls, &self.weights, dec)?,
            domain,
            dec,
        )?;
        let (controls, weights) = euclidean(net, dec)?;
        Self::new(self.degree, knots, controls, weights, dec)
    }
    /// Clamps the active interval, reusing storage when it is already clamped.
    pub fn into_clamped(self, dec: &Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if clamped_axis(self.degree, &self.knots, self.controls.len()) {
            Ok(self)
        } else {
            self.extracted(self.domain(), dec)
        }
    }
    /// Catallax's certified evaluator over the exact active interval. This does
    /// not certify the curve's incidence on an imported face.
    pub fn to_curve(&self, dec: &mut Decider) -> Result<BSplineCurve, Refusal> {
        let clamped;
        let source = if clamped_axis(self.degree, &self.knots, self.controls.len()) {
            self
        } else {
            clamped = self.extracted(self.domain(), dec)?;
            &clamped
        };
        BSplineCurve::new(
            source.degree,
            source.knots.clone(),
            real_points(&source.controls),
            source
                .weights
                .iter()
                .cloned()
                .map(Real::from_rational)
                .collect(),
            dec,
        )
    }
}

fn real_points(points: &[RationalPoint]) -> Vec<Point3> {
    points
        .iter()
        .map(|p| {
            let [x, y, z] = p.clone().map(Real::from_rational);
            Point3::new(x, y, z)
        })
        .collect()
}

/// Exact u-major, v-minor tensor-product source surface. This is control data,
/// not a checked `SplineSupport`, and cannot authorize solid decisions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RationalSplineSurface {
    degrees: [usize; 2],
    knots: [Vec<Rational>; 2],
    count: [usize; 2],
    controls: Vec<RationalPoint>,
    weights: Vec<Rational>,
}

impl RationalSplineSurface {
    /// Validates both axes, rectangular storage and positive weights.
    pub fn new(
        degrees: [usize; 2],
        knots: [Vec<Rational>; 2],
        count: [usize; 2],
        controls: Vec<RationalPoint>,
        weights_: Vec<Rational>,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let size = count[0]
            .checked_mul(count[1])
            .ok_or_else(|| invalid("surface size overflow"))?;
        if controls.len() != size {
            return Err(invalid("nonrectangular surface controls"));
        }
        for i in 0..2 {
            axis(degrees[i], &knots[i], count[i], dec)?;
        }
        weights(&weights_, size, dec)?;
        Ok(Self {
            degrees,
            knots,
            count,
            controls,
            weights: weights_,
        })
    }
    /// Tensor degrees.
    pub fn degrees(&self) -> [usize; 2] {
        self.degrees
    }
    /// Original expanded knot vectors.
    pub fn knots(&self) -> [&[Rational]; 2] {
        self.knots.each_ref().map(Vec::as_slice)
    }
    /// U, V control counts.
    pub fn control_count(&self) -> [usize; 2] {
        self.count
    }
    /// Euclidean controls in millimetres, u-major and v-minor.
    pub fn controls(&self) -> &[RationalPoint] {
        &self.controls
    }
    /// Positive dimensionless weights in control order.
    pub fn weights(&self) -> &[Rational] {
        &self.weights
    }
    /// Active parameter rectangle.
    pub fn domain(&self) -> [[&Rational; 2]; 2] {
        std::array::from_fn(|i| {
            [
                &self.knots[i][self.degrees[i]],
                &self.knots[i][self.count[i]],
            ]
        })
    }
    /// Exact rectangular restriction, applying the donor's curve extraction in
    /// each tensor direction. No sampling or fitted controls are involved.
    pub fn extracted(&self, domain: [[&Rational; 2]; 2], dec: &Decider) -> Result<Self, Refusal> {
        let source = weighted(&self.controls, &self.weights, dec)?;
        let mut columns = Vec::with_capacity(self.count[1]);
        let mut ku = Vec::new();
        for v in 0..self.count[1] {
            dec.check_deadline()?;
            let column: Vec<_> = (0..self.count[0])
                .map(|u| source[u * self.count[1] + v].clone())
                .collect();
            let (knots, points) =
                extract(self.degrees[0], &self.knots[0], &column, domain[0], dec)?;
            ku = knots;
            columns.push(points);
        }
        let count_u = columns[0].len();
        let mut net = Vec::new();
        let mut kv = Vec::new();
        let mut count_v = 0;
        for u in 0..count_u {
            let row: Vec<_> = columns.iter().map(|c| c[u].clone()).collect();
            let (knots, points) = extract(self.degrees[1], &self.knots[1], &row, domain[1], dec)?;
            kv = knots;
            count_v = points.len();
            net.extend(points);
        }
        let (controls, weights) = euclidean(net, dec)?;
        Self::new(
            self.degrees,
            [ku, kv],
            [count_u, count_v],
            controls,
            weights,
            dec,
        )
    }
    /// Clamps both active intervals, reusing an already clamped source net.
    pub fn into_clamped(self, dec: &Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if (0..2).all(|i| clamped_axis(self.degrees[i], &self.knots[i], self.count[i])) {
            Ok(self)
        } else {
            self.extracted(self.domain(), dec)
        }
    }
    /// Catallax's existing certified interval evaluator over the exact source
    /// domain. Surface regularity, trims and solid validity remain unproved.
    pub fn to_surface(&self, dec: &mut Decider) -> Result<BSplineSurface, Refusal> {
        let clamped;
        let source = if (0..2).all(|i| clamped_axis(self.degrees[i], &self.knots[i], self.count[i]))
        {
            self
        } else {
            clamped = self.extracted(self.domain(), dec)?;
            &clamped
        };
        BSplineSurface::new(
            source.degrees,
            source.knots.clone(),
            source.count,
            real_points(&source.controls),
            source
                .weights
                .iter()
                .cloned()
                .map(Real::from_rational)
                .collect(),
            dec,
        )
    }
}

fn clamped_axis(degree: usize, knots: &[Rational], count: usize) -> bool {
    knots[..=degree].iter().all(|k| k == &knots[degree])
        && knots[count..].iter().all(|k| k == &knots[count])
}
#[cfg(test)]
mod tests;
