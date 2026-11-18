//! Guide-normal sections, following OCCT BRepBlend_Walking and SMLib's
//! constant-radius offset equality plus a section-plane equation.

use super::*;

struct SectionPlane {
    point: [Ball; 3],
    tangent: [Ball; 3],
    inverse_length: Ball,
    point64: [f64; 3],
    normal64: [f64; 3],
}

impl SectionPlane {
    fn new(guide: &BSplineCurve, tau: &Rational, dec: &mut Decider) -> Result<Self, Refusal> {
        // Continue the endpoint two-jet, not a tangent line: a curvature jump
        // in the guide changes the offset walk's speed at the join. This only
        // selects sections outside the edge; it does not extend support faces.
        let domain = guide.domain();
        let anchor = if tau < domain[0] {
            domain[0]
        } else if tau > domain[1] {
            domain[1]
        } else {
            tau
        };
        let [point, tangent] = if anchor == tau {
            guide.first_box_reusing([anchor, anchor], 128, dec, &mut Vec::new())?
        } else {
            let jet = guide.jet_box([anchor, anchor], 128, dec)?;
            let delta = tau - anchor;
            let offset = Ball::from_rational(&delta, 128);
            let quadratic = Ball::from_rational(&(&delta * &delta).mul_pow2(-1), 128);
            let point = std::array::from_fn(|i| {
                let linear = &jet.value[i] + &(&jet.first[i] * &offset);
                &linear + &(&jet.second[i] * &quadratic)
            });
            let tangent = std::array::from_fn(|i| &jet.first[i] + &(&jet.second[i] * &offset));
            [point, tangent]
        };
        dec.charge_work(16)?;
        let square = &(&tangent[0].square() + &tangent[1].square()) + &tangent[2].square();
        if square.lo().sign() != Sign::Positive {
            return Err(bad(
                "rolling-ball guide has a singular or unresolved tangent",
            ));
        }
        let inverse_length = square
            .sqrt()
            .and_then(|length| length.recip())
            .ok_or_else(|| bad("rolling-ball guide tangent cannot be normalized"))?;
        let point64 = point.each_ref().map(|value| value.mid().to_f64_nearest());
        let tangent64 = tangent.each_ref().map(|value| value.mid().to_f64_nearest());
        let length64 = norm(tangent64);
        let normal64 = tangent64.map(|value| value / length64);
        if !point64
            .iter()
            .chain(&normal64)
            .all(|value| value.is_finite())
        {
            return Err(bad("rolling-ball guide is outside binary64 range"));
        }
        Ok(Self {
            point,
            tangent,
            inverse_length,
            point64,
            normal64,
        })
    }

    fn residual(&self, centre: [f64; 3]) -> f64 {
        (0..3)
            .map(|i| (centre[i] - self.point64[i]) * self.normal64[i])
            .sum()
    }

    fn checked_error(
        &self,
        centre: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Rational, Refusal> {
        dec.charge_work(10)?;
        // Keep the common positive inverse length outside the dot product:
        // distributing its interval loses correlation between cancelling terms.
        let residual = (0..3).fold(Ball::from_rational(&Rational::zero(), 128), |sum, i| {
            let delta = &Ball::from_rational(&centre[i], 128) - &self.point[i];
            sum + (&delta * &self.tangent[i])
        });
        let residual = &residual * &self.inverse_length;
        Ok(residual.lo().abs().max(residual.hi().abs()))
    }
}

fn residual_norm<const N: usize>(values: [f64; N]) -> f64 {
    values.iter().map(|value| value * value).sum::<f64>().sqrt()
}

/// The iterate is only a proposal. Callers independently enclose each accepted
/// station, including the constraint used to select its branch.
pub(super) fn newton<const N: usize>(
    mut x: [f64; N],
    tolerance: f64,
    mut residual: impl FnMut([f64; N], &mut Decider) -> Result<[f64; N], Refusal>,
    dec: &mut Decider,
) -> Result<[f64; N], Refusal> {
    for _ in 0..40 {
        dec.check_deadline()?;
        let value = residual(x, dec)?;
        let error = residual_norm(value);
        if error <= tolerance / 64.0 {
            return Ok(x);
        }
        let mut jacobian = [[0.0; N]; N];
        for j in 0..N {
            let h = 1e-6 * x[j].abs().max(1.0);
            let mut up = x;
            let mut down = x;
            up[j] += h;
            down[j] -= h;
            let (difference, width): ([f64; N], f64) =
                match (residual(up, dec), residual(down, dec)) {
                    (Ok(a), Ok(b)) => (std::array::from_fn(|i| a[i] - b[i]), 2.0 * h),
                    (Ok(a), Err(e)) if e.code == "fillet.numerical-spine" => {
                        (std::array::from_fn(|i| a[i] - value[i]), h)
                    }
                    (Err(e), Ok(b)) if e.code == "fillet.numerical-spine" => {
                        (std::array::from_fn(|i| value[i] - b[i]), h)
                    }
                    (Err(e), _) | (_, Err(e)) => return Err(e),
                };
            for i in 0..N {
                jacobian[i][j] = difference[i] / width;
            }
        }
        dec.charge_work((4 * N * N * N) as u64)?;
        let step = inverse(jacobian, value)
            .ok_or_else(|| bad("singular numerical rolling-ball Jacobian"))?;
        let mut accepted = None;
        let mut trial_failure = None;
        for power in 0..16 {
            let share = 2.0f64.powi(-power);
            let next = std::array::from_fn(|i| x[i] - share * step[i]);
            match residual(next, dec) {
                Ok(value) if residual_norm(value) < error => {
                    accepted = Some(next);
                    break;
                }
                Ok(_) => {}
                Err(e) if e.code == "fillet.numerical-spine" => trial_failure = Some(e),
                Err(e) => return Err(e),
            }
        }
        x = accepted.ok_or_else(|| {
            let mut refusal = bad(format!(
                "rolling-ball Newton step did not reduce the residual; residual approx {error:.9e} mm, iterate {x:?}, full step {step:?}",
            ));
            if let Some(failed) = trial_failure {
                refusal.detail.push_str(&format!("; rejected trial: {}", failed.detail));
            }
            refusal
        })?;
    }
    Err(bad("rolling-ball Newton iteration limit reached"))
}

impl Problem {
    pub(super) fn solve_guided(
        &self,
        guide: &BSplineCurve,
        tau: &Rational,
        seed: [f64; 4],
        dec: &mut Decider,
    ) -> Result<NumericalSpineNode, Refusal> {
        let plane = SectionPlane::new(guide, tau, dec)?;
        let x = newton(
            seed,
            self.tolerance,
            |x, dec| {
                let a = self.supports[0].centre([x[0], x[1]], dec)?;
                let b = self.supports[1].centre([x[2], x[3]], dec)?;
                Ok([a[0] - b[0], a[1] - b[1], a[2] - b[2], plane.residual(a)])
            },
            dec,
        )?;
        let parameters = [
            [rational(x[0])?, rational(x[1])?],
            [rational(x[2])?, rational(x[3])?],
        ];
        let (ca, a) = self.supports[0].checked(&parameters[0], dec)?;
        let (cb, b) = self.supports[1].checked(&parameters[1], dec)?;
        let centre = std::array::from_fn(|i| (&a[i].mid() + &b[i].mid()).mul_pow2(-1));
        let offset = offset_error(&centre, [&a, &b]);
        let section = plane.checked_error(&centre, dec)?;
        let offset_upper = offset.to_f64_nearest();
        let section_upper = section.to_f64_nearest();
        let limit = self.tolerance / 8.0;
        if offset_upper > limit || section_upper > limit {
            let midpoint_residual = plane
                .residual(centre.each_ref().map(Rational::to_f64_nearest))
                .abs();
            let diameter = |values: &[Ball; 3]| {
                values
                    .iter()
                    .map(|value| (value.hi() - value.lo()).to_f64_nearest())
                    .fold(0.0, f64::max)
            };
            let point_diameter = diameter(&plane.point);
            let normal = plane
                .tangent
                .each_ref()
                .map(|value| value * &plane.inverse_length);
            let normal_diameter = diameter(&normal);
            return Err(bad(format!(
                "rolling-ball station at tau {tau}: offset residual upper {offset_upper:.9e} mm (exceeds limit: {}); guide-plane residual upper {section_upper:.9e} mm (exceeds limit: {}); limit {limit:.9e} mm; checked-centre midpoint plane residual {midpoint_residual:.9e} mm; guide point diameter {point_diameter:.9e} mm; guide unit-normal diameter {normal_diameter:.9e}",
                offset_upper > limit,
                section_upper > limit,
            )));
        }
        Ok(NumericalSpineNode {
            tau: tau.clone(),
            parameters,
            charts: [ca, cb],
            centre,
        })
    }
}
