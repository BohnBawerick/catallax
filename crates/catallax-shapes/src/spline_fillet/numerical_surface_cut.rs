//! Selected transverse neighbour cuts using the existing surface marcher.
//! As in OCCT IntWalk and SMLib's Hermite 3D/UV checks, retain both native
//! curves and check their compositions. This is not a complete branch census.

use super::*;
use crate::{SplineSupport, Surface};
use catallax_numbers::{NumericalCheck, NumericalSample};

/// One numerical graph across a fillet, with a native curve on both surfaces.
/// Face-domain clipping is separate from these untrimmed carrier equations.
#[derive(Clone, Debug)]
pub struct NumericalSurfaceCut {
    curve: BSplineCurve,
    parameters: [BSplineCurve; 2],
    tau_ends: [Rational; 2],
    sampled_error_mm: Rational,
    check: NumericalCheck,
}

impl NumericalSurfaceCut {
    /// World curve from the first fillet contact rail to the second.
    pub fn curve(&self) -> &BSplineCurve {
        &self.curve
    }
    /// Native curves at the world curve's parameter: `[fillet, neighbour]`.
    pub fn parameters(&self) -> &[BSplineCurve; 2] {
        &self.parameters
    }
    /// Longitudinal fillet parameters at the first and second contact rails.
    pub fn tau_ends(&self) -> &[Rational; 2] {
        &self.tau_ends
    }
    /// Largest checked 3D/native composition residual, not a Hausdorff bound.
    pub fn sampled_error_mm(&self) -> &Rational {
        &self.sampled_error_mm
    }
    /// The finite sample locations and their enclosed composition residuals.
    pub fn check(&self) -> &NumericalCheck {
        &self.check
    }
}

fn native_curve(
    nodes: &[NumericalSpineNode],
    side: usize,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let (knots, controls) = carrier::hermite(nodes, |i| &nodes[i].parameters[side]);
    let weights = vec![Real::one(); controls.len()];
    let points = controls
        .into_iter()
        .map(|[u, v]| Point3::new(u.into(), v.into(), Real::zero()))
        .collect();
    BSplineCurve::new(3, knots, points, weights, dec)
}

fn composition_error(a: &[Ball; 3], b: &[Ball; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        sum + (&a[i].lo() - &b[i].hi())
            .abs()
            .max((&a[i].hi() - &b[i].lo()).abs())
    })
}

impl NumericalSplineFillet {
    /// Intersects a neighbour as one selected graph `tau(s)` across the strip.
    /// `seed` gives longitudinal fillet tau and both neighbour parameters at
    /// the first rail. All regular support adapters use the same zero-offset
    /// intersection equations. Tangencies, missing charts and failed solves
    /// retain their typed causes; no exact intersection certificate is required.
    pub fn surface_cut(
        &self,
        neighbour: Surface,
        seed: [Rational; 3],
        dec: &mut Decider,
    ) -> Result<NumericalSurfaceCut, Refusal> {
        let fillet = Surface::BSpline(SplineSupport::RationalAtlas(self.surface().clone()));
        let surfaces = [fillet, neighbour];
        eprintln!("cut attribution start work={}", catallax_numbers::current_operation_work().used);
        let supports = [
            NumericalFilletSupport::new(surfaces[0].clone(), Rational::zero(), dec)?,
            NumericalFilletSupport::new(surfaces[1].clone(), Rational::zero(), dec)?,
        ];
        eprintln!("cut attribution supports work={}", catallax_numbers::current_operation_work().used);
        let limit = self.spine().tolerance_mm().mul_pow2(-2);
        let mut spine = NumericalFilletSpine::trace(
            supports,
            1,
            [Rational::zero(), Rational::one()],
            seed,
            self.spine().tolerance_mm().mul_pow2(-1),
            dec,
        )?;
        eprintln!("cut attribution trace work={} nodes={}", catallax_numbers::current_operation_work().used, spine.nodes.len());
        for depth in 0..=MAX_DEPTH {
            let parameters = [
                native_curve(&spine.nodes, 0, dec)?,
                native_curve(&spine.nodes, 1, dec)?,
            ];
            eprintln!("cut attribution depth={depth} native curves work={} nodes={}", catallax_numbers::current_operation_work().used, spine.nodes.len());
            let mut upper = Rational::zero();
            let mut samples = Vec::with_capacity((spine.nodes.len() - 1) * 8 + 1);
            let mut refine = Vec::new();
            for (index, pair) in spine.nodes.windows(2).enumerate() {
                let mut misses = false;
                for j in 0..=8 {
                    if index != 0 && j == 0 {
                        continue;
                    }
                    dec.charge_work(1)?;
                    let share = Rational::from_ratio(j, 8).expect("nonzero denominator");
                    let across = &pair[0].tau + &((&pair[1].tau - &pair[0].tau) * share);
                    let point = spine
                        .centre_curve
                        .evaluate_box([&across, &across], 128, dec)?;
                    let mut error = Rational::zero();
                    for side in 0..2 {
                        let uv = parameters[side].evaluate_box([&across, &across], 128, dec)?;
                        let composed = surfaces[side].evaluate(&uv[0], &uv[1], dec)?;
                        error = error.max(composition_error(&point, &composed));
                    }
                    misses |= error > limit;
                    upper = upper.max(error.clone());
                    samples.push(NumericalSample {
                        parameters: vec![across],
                        distance_upper_mm: error,
                    });
                }
                if misses {
                    refine.push(index);
                }
            }
            eprintln!("cut attribution depth={depth} samples work={} refine={}", catallax_numbers::current_operation_work().used, refine.len());
            if upper <= limit {
                return Ok(NumericalSurfaceCut {
                    tau_ends: [
                        spine.nodes[0].parameters[0][0].clone(),
                        spine.nodes.last().expect("marched nodes").parameters[0][0].clone(),
                    ],
                    curve: spine.centre_curve,
                    parameters,
                    sampled_error_mm: upper,
                    check: NumericalCheck::Samples {
                        name: "fillet.neighbour.world-and-both-native-curves".into(),
                        scope: "selected transverse neighbour cut".into(),
                        samples,
                    },
                });
            }
            if depth == MAX_DEPTH
                || spine.nodes.len() + refine.len() > MAX_NODES
                || refine.is_empty()
            {
                return Err(bad(format!(
                    "neighbour cut reaches sampled world/native residual {upper} mm against {limit} mm"
                )));
            }
            for index in refine.into_iter().rev() {
                let a = &spine.nodes[index];
                let b = &spine.nodes[index + 1];
                let across = (&a.tau + &b.tau).mul_pow2(-1);
                let seed = std::array::from_fn(|i| {
                    (&a.parameters[i / 2][i % 2] + &b.parameters[i / 2][i % 2])
                        .mul_pow2(-1)
                        .to_f64_nearest()
                });
                let node = spine.problem.solve(&across, seed, dec)?;
                spine.nodes.insert(index + 1, node);
            }
            spine.centre_curve = centre_fit(&spine.nodes, dec)?;
        }
        unreachable!("last refinement level returns its measured refusal")
    }
}
