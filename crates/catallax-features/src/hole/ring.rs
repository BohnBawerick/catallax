//! Native parameter correspondence for shared circular stock boundaries.

use catallax_boolean::numerical::NumericalCylinderSource;
use catallax_numbers::{
    Ball, Decider, NumericalCheck, NumericalEvidence, NumericalSample, Rational, Real, Refusal,
    Sign, operation_budget_active, with_operation_budget,
};
use catallax_shapes::{BSplineCurve, Point3, Surface};
use catallax_topology::numerical::{DEFAULT_WORK, NumericalOptions};

const BITS: u32 = 128;

/// A numerical pcurve at the unchanged shared 3D ring's native parameter.
#[derive(Clone, Debug)]
pub struct CylinderRingCurve {
    key: String,
    quarter: usize,
    curve: BSplineCurve,
    evidence: NumericalEvidence,
}

impl CylinderRingCurve {
    /// Caller identity for this shared ring fragment.
    pub fn key(&self) -> &str {
        &self.key
    }
    /// The actual canonical cylinder quarter owning this pcurve.
    pub fn quarter(&self) -> usize {
        self.quarter
    }
    /// Native cylinder UV at the original 3D edge parameter, without refitting it.
    pub fn curve(&self) -> &BSplineCurve {
        &self.curve
    }
    /// Positive parameter hull and actual endpoint/composition checks.
    pub fn evidence(&self) -> &NumericalEvidence {
        &self.evidence
    }
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("hole.ring-parameters", detail)
}

fn distance(a: &[Ball; 3], b: &[Ball; 3]) -> Result<Rational, Refusal> {
    let delta: [Ball; 3] = std::array::from_fn(|i| &a[i] - &b[i]);
    Ok((delta[0].square() + delta[1].square() + delta[2].square())
        .sqrt()
        .ok_or_else(|| bad("shared ring distance is unresolved"))?
        .hi())
}

fn point(net: &[Real; 4], dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_budgeted_at(&(&net[0] / &net[3]), BITS)?,
        dec.enclose_budgeted_at(&(&net[1] / &net[3]), BITS)?,
        dec.enclose_budgeted_at(&(&net[2] / &net[3]), BITS)?,
    ])
}

/// Generates a pcurve for one native quadratic circular ring fragment on an
/// actual canonical cylinder quarter. The caller selects the shared original
/// ring and its cylinder boundary row, zero or one. The finite interval remains
/// in the original 3D carrier's parameter, including descending edge ranges.
/// `None` means no single corresponding quarter was found for this source span.
/// Checks remain numerical under the caller's preset tolerance. Complete shell
/// admission and approximate export remain separate obligations.
pub fn cylinder_ring_curve(
    key: &str,
    source: &BSplineCurve,
    parameters: [Rational; 2],
    cylinder: &NumericalCylinderSource,
    row: usize,
    options: &NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<CylinderRingCurve>, Refusal> {
    let run = |dec: &mut Decider| build(key, source, parameters, cylinder, row, options, dec);
    if operation_budget_active() {
        return run(dec);
    }
    let (result, work) = with_operation_budget(DEFAULT_WORK, || run(dec));
    if work.exhausted {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "shared cylinder ring parameters exhausted their operation allowance",
        ));
    }
    result
}

fn build(
    key: &str,
    source: &BSplineCurve,
    parameters: [Rational; 2],
    cylinder: &NumericalCylinderSource,
    row: usize,
    options: &NumericalOptions,
    dec: &mut Decider,
) -> Result<Option<CylinderRingCurve>, Refusal> {
    dec.charge_work(1)?;
    if key.trim().is_empty() || row > 1 || options.tolerance_mm.sign() != Sign::Positive {
        return Err(bad(
            "ring identity, cylinder boundary row and positive tolerance are required",
        ));
    }
    if source.degree() != 2 || parameters[0] == parameters[1] {
        return Ok(None);
    }
    let lo = parameters[0].clone().min(parameters[1].clone());
    let hi = parameters[0].clone().max(parameters[1].clone());
    let mut selected = None;
    for (index, bounds) in source.spans() {
        dec.charge_work(1)?;
        if lo >= *bounds[0] && hi <= *bounds[1] {
            selected = Some((index, bounds));
            break;
        }
    }
    let Some((index, bounds)) = selected else {
        return Ok(None);
    };
    let net = source.span_bezier_controls(index, dec)?;
    let ends = [point(&net[0], dec)?, point(&net[2], dec)?];
    let v = Rational::from_i64(row as i64);
    let mut owners = Vec::new();
    for quarter in cylinder.quarters() {
        dec.charge_work(1)?;
        let Surface::BSpline(support) = quarter.surface() else {
            return Err(bad("canonical cylinder lost its native support"));
        };
        let target = support.carrier();
        let [u, _] = target.domain();
        let corners = [
            target.evaluate_box([u[0], u[0]], [&v, &v], BITS, dec)?,
            target.evaluate_box([u[1], u[1]], [&v, &v], BITS, dec)?,
        ];
        for reverse in [false, true] {
            let order = if reverse { [1, 0] } else { [0, 1] };
            let errors = [
                distance(&ends[0], &corners[order[0]])?,
                distance(&ends[1], &corners[order[1]])?,
            ];
            if errors.iter().all(|error| *error <= options.tolerance_mm) {
                owners.push((quarter, reverse, errors));
            }
        }
    }
    let (quarter, reverse, errors) = match owners.as_slice() {
        [] => return Ok(None),
        [owner] => owner,
        _ => return Err(bad("shared ring has ambiguous cylinder-quarter ownership")),
    };
    let Surface::BSpline(support) = quarter.surface() else {
        unreachable!()
    };
    let target = support.carrier();
    if target.degrees() != [2, 1] || target.control_count() != [3, 2] {
        return Err(bad("shared ring needs a quadratic/linear cylinder quarter"));
    }
    let source_ratio = &net[1][3] / &net[0][3];
    let target_ratio = &target.weights()[2 + row] / &target.weights()[row];
    let weight = if *reverse {
        &(&Real::from_i64(2) * &source_ratio) * &target_ratio
    } else {
        &source_ratio / &target_ratio
    };
    if dec.enclose_budgeted_at(&weight, BITS)?.lo().sign() != Sign::Positive {
        return Err(bad("ring parameter weights must be strictly positive"));
    }
    let [mut u, _] = target.domain();
    if *reverse {
        u.reverse();
    }
    let curve = BSplineCurve::new(
        1,
        vec![
            bounds[0].clone(),
            bounds[0].clone(),
            bounds[1].clone(),
            bounds[1].clone(),
        ],
        u.map(|u| Point3::new(Real::from(u.clone()), Real::from(v.clone()), Real::zero()))
            .to_vec(),
        vec![Real::one(), weight],
        dec,
    )?;
    let mut samples = Vec::new();
    for i in 0..=16 {
        dec.charge_work(1)?;
        let t = &parameters[0]
            + &((&parameters[1] - &parameters[0]) * Rational::from_ratio(i, 16).unwrap());
        let at = if &t == bounds[0] {
            curve.control(0).x.clone()
        } else if &t == bounds[1] {
            curve.control(1).x.clone()
        } else {
            curve.point_at(&Real::from(t.clone()), dec)?.x
        };
        let u = dec.enclose_budgeted_at(&at, BITS)?;
        let actual = target.evaluate_box([&u.lo(), &u.hi()], [&v, &v], BITS, dec)?;
        let original = source.evaluate_box([&t, &t], BITS, dec)?;
        samples.push(NumericalSample {
            parameters: vec![t, u.mid(), v.clone()],
            distance_upper_mm: distance(&original, &actual)?,
        });
    }
    let evidence = NumericalEvidence::new(
        options.tolerance_mm.clone(),
        vec![
            NumericalCheck::Passed {
                name: "hole.ring.positive-native-parameter-hull".into(),
                scope: key.into(),
            },
            NumericalCheck::Samples {
                name: "hole.ring.original-quarter-endpoints".into(),
                scope: key.into(),
                samples: bounds
                    .iter()
                    .zip(errors)
                    .map(|(t, error)| NumericalSample {
                        parameters: vec![(*t).clone()],
                        distance_upper_mm: error.clone(),
                    })
                    .collect(),
            },
            NumericalCheck::Samples {
                name: "hole.ring.same-native-parameter-composition".into(),
                scope: key.into(),
                samples,
            },
        ],
        vec![],
    )?;
    Ok(Some(CylinderRingCurve {
        key: key.into(),
        quarter: quarter.id(),
        curve,
        evidence,
    }))
}
