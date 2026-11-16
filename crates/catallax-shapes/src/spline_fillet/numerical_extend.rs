//! Exact polynomial extension of only a fillet support's two end rows.
//! The original charts are kept on their old domains; this supplies room for
//! numerical end-face cuts, not a general surface-extension operation.

use super::*;
use crate::{BSplineSurface, RationalPatch};
use std::collections::BTreeSet;

type Domain = [[Rational; 2]; 2];

fn homogeneous(patch: &RationalPatch, dec: &mut Decider) -> Result<Vec<[Rational; 4]>, Refusal> {
    let carrier = patch.carrier();
    dec.charge_work(carrier.controls().len() as u64 * 8)?;
    carrier.controls().iter().zip(carrier.weights()).map(|(point, weight)| {
        let p = point.as_rational().ok_or_else(|| bad("fillet end extension needs rational chart controls"))?;
        let w = weight.as_rational().ok_or_else(|| bad("fillet end extension needs rational chart weights"))?;
        for value in [&p.x, &p.y, &p.z] { dec.reserve_exact_arithmetic(value, &w, true)?; }
        Ok([&p.x * &w, &p.y * &w, &p.z * &w, w])
    }).collect()
}
fn knot_vector(cuts: &[Rational], degree: usize) -> Vec<Rational> {
    let mut knots = vec![cuts[0].clone(); degree + 1];
    for value in &cuts[1..cuts.len() - 1] { knots.extend(std::iter::repeat_n(value.clone(), degree)); }
    knots.extend(std::iter::repeat_n(cuts.last().expect("two cuts").clone(), degree + 1));
    knots
}

impl NumericalFilletSpine {
    /// Extends the outer charts along the selected edge parameter by `fraction`
    /// of the support's width at each end. Exact homogeneous seam agreement is
    /// required before joining the extended chart nets. The resulting atlas
    /// constructor rechecks regularity and G1 joins. Physical face trims do not
    /// change; a numerical fillet must still be cut back to its end faces.
    pub fn extended_support(
        atlas: &RationalAtlas,
        axis: usize,
        fraction: &Rational,
        dec: &mut Decider,
    ) -> Result<RationalAtlas, Refusal> {
        if axis > 1 || fraction.sign() != Sign::Positive { return Err(bad("fillet end extension needs a valid axis and positive fraction")); }
        let mut cuts: [BTreeSet<Rational>; 2] = std::array::from_fn(|_| BTreeSet::new());
        for chart in atlas.charts() {
            dec.charge_work(4)?;
            for (cut, bounds) in cuts.iter_mut().zip(chart.patch().domain()) { cut.extend(bounds.iter().cloned()); }
        }
        let mut cuts = cuts.map(|cut| cut.into_iter().collect::<Vec<_>>());
        let count = [cuts[0].len() - 1, cuts[1].len() - 1];
        if count[0].checked_mul(count[1]) != Some(atlas.charts().len()) {
            return Err(bad("fillet end extension requires a complete rectangular chart cover"));
        }
        let last = cuts[axis].len() - 1;
        let delta = &(&cuts[axis][last] - &cuts[axis][0]) * fraction;
        let mut nets = Vec::with_capacity(atlas.charts().len());
        let mut degree = None;
        for i in 0..count[0] {
            for j in 0..count[1] {
                let old: Domain = [
                    [cuts[0][i].clone(), cuts[0][i + 1].clone()],
                    [cuts[1][j].clone(), cuts[1][j + 1].clone()],
                ];
                let chart = atlas.charts().iter().find(|chart| chart.patch().domain() == &old)
                    .ok_or_else(|| bad("fillet end extension has a missing or overlapping chart cell"))?;
                let mut domain = old.clone();
                if old[axis][0] == cuts[axis][0] { domain[axis][0] = &domain[axis][0] - &delta; }
                if old[axis][1] == cuts[axis][last] { domain[axis][1] = &domain[axis][1] + &delta; }
                let extended = chart.patch().extended(domain, dec)?;
                let this_degree = extended.carrier().degrees();
                if degree.is_some_and(|degree| degree != this_degree) { return Err(bad("fillet end charts have different tensor degrees")); }
                degree = Some(this_degree);
                nets.push(homogeneous(&extended, dec)?);
            }
        }
        let [p, q] = degree.ok_or_else(|| bad("fillet support has no charts"))?;
        if p == 0 || q == 0 {
            return Err(bad("fillet end extension needs two nonconstant parameters"));
        }
        for i in 0..count[0] {
            for j in 0..count[1] {
                let here = &nets[i * count[1] + j];
                if i + 1 < count[0] {
                    let next = &nets[(i + 1) * count[1] + j];
                    for k in 0..=q {
                        dec.charge_work(4)?;
                        if here[p * (q + 1) + k] != next[k] { return Err(bad("extended fillet charts disagree on a homogeneous u seam")); }
                    }
                }
                if j + 1 < count[1] {
                    let next = &nets[i * count[1] + j + 1];
                    for k in 0..=p {
                        dec.charge_work(4)?;
                        if here[k * (q + 1) + q] != next[k * (q + 1)] { return Err(bad("extended fillet charts disagree on a homogeneous v seam")); }
                    }
                }
            }
        }
        cuts[axis][0] = &cuts[axis][0] - &delta;
        cuts[axis][last] = &cuts[axis][last] + &delta;
        let sizes = [count[0] * p + 1, count[1] * q + 1];
        let mut controls = Vec::with_capacity(sizes[0] * sizes[1]);
        let mut weights = Vec::with_capacity(sizes[0] * sizes[1]);
        for i in 0..sizes[0] {
            let chart_u = (i / p).min(count[0] - 1);
            let local_u = i - chart_u * p;
            for j in 0..sizes[1] {
                let chart_v = (j / q).min(count[1] - 1);
                let local_v = j - chart_v * q;
                let h = &nets[chart_u * count[1] + chart_v][local_u * (q + 1) + local_v];
                dec.charge_work(8)?;
                controls.push(Point3::new(
                    Real::from_rational(&h[0] / &h[3]), Real::from_rational(&h[1] / &h[3]), Real::from_rational(&h[2] / &h[3]),
                ));
                weights.push(Real::from_rational(h[3].clone()));
            }
        }
        let source = BSplineSurface::new(
            [p, q], [knot_vector(&cuts[0], p), knot_vector(&cuts[1], q)], sizes, controls, weights, dec,
        )?;
        let [extended] = <[_; 1]>::try_from(RationalAtlas::partition(source, [false, false], dec)?)
            .map_err(|_| bad("extended fillet support is not one regular G1 logical face"))?;
        Ok(extended)
    }
}
