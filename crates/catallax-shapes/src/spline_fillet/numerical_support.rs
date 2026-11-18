//! Signed offsets used by numerical rolling-ball continuation.

use super::*;
use crate::spline_fillet::{SUPPORT_SPANS, definition};
use crate::{RationalPatch, SplineSupport, Surface};
use std::sync::OnceLock;

#[derive(Debug)]
enum Geometry {
    Atlas(Vec<AtlasSupport>),
    Native {
        surface: Surface,
        proposal: Option<super::analytic::Analytic64>,
    },
}

#[derive(Debug)]
struct AtlasSupport {
    patch: RationalPatch,
    domain: [[Rational; 2]; 2],
    side: Sign,
    expanded: OnceLock<Support>,
}

impl AtlasSupport {
    fn new(patch: &RationalPatch, side: Sign) -> Result<Self, Refusal> {
        // Keep the eager adapter's admission checks for every chart. Only
        // the exact tangent/normal expansion waits for a chart to be visited.
        let spans = patch.rational_span_polynomials().ok_or_else(|| {
            definition("a spline fillet needs supports with rational coefficients")
        })?;
        let [(domain, _)] = <[_; 1]>::try_from(spans).map_err(|_| {
            Refusal::new(
                SUPPORT_SPANS,
                "a spline fillet on a support of more than one knot span is not covered yet",
            )
        })?;
        Ok(Self {
            patch: patch.clone(),
            domain,
            side,
            expanded: OnceLock::new(),
        })
    }

    fn get(&self, dec: &mut Decider) -> Result<&Support, Refusal> {
        if let Some(support) = self.expanded.get() {
            return Ok(support);
        }
        let support = Support::new(&self.patch, self.side, dec)?;
        let _ = self.expanded.set(support);
        Ok(self.expanded.get().expect("initialized support"))
    }
}

/// A support's native parameterization and signed rolling-ball distance.
/// Binary64 values propose contacts. Interval evaluations check accepted ones.
#[derive(Debug)]
pub struct NumericalFilletSupport {
    geometry: Geometry,
    pub(super) distance: Rational,
    distance64: f64,
}

pub(super) struct Contact {
    pub(super) chart: usize,
    pub(super) point: [Ball; 3],
    pub(super) normal: [Ball; 3],
}

impl NumericalFilletSupport {
    /// Retains the surface's native parameters. A negative distance places the
    /// ball opposite `S_u cross S_v`; zero traces the unoffset sharp guide.
    pub fn new(
        surface: Surface,
        signed_distance: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let distance64 = signed_distance.to_f64_nearest();
        if !distance64.is_finite() || (!signed_distance.is_zero() && distance64 == 0.0) {
            return Err(bad("support offset distance is outside binary64 range"));
        }
        if let Surface::BSpline(SplineSupport::RationalAtlas(atlas)) = &surface {
            return Self::from_atlas(atlas, signed_distance, dec);
        }
        Ok(Self {
            geometry: Geometry::Native {
                proposal: super::analytic::Analytic64::new(&surface, dec)?,
                surface,
            },
            distance: signed_distance,
            distance64,
        })
    }

    pub(super) fn from_atlas(
        atlas: &RationalAtlas,
        distance: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let side = if distance.sign() == Sign::Negative {
            Sign::Negative
        } else {
            Sign::Positive
        };
        let charts = atlas
            .charts()
            .iter()
            .map(|chart| {
                dec.check_deadline()?;
                AtlasSupport::new(chart.patch(), side)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            geometry: Geometry::Atlas(charts),
            distance64: distance.to_f64_nearest(),
            distance,
        })
    }

    fn chart<'a>(
        charts: &'a [AtlasSupport],
        at: [f64; 2],
        dec: &mut Decider,
    ) -> Result<(usize, &'a Support), Refusal> {
        let (index, chart) = charts
            .iter()
            .enumerate()
            .find(|(_, chart)| {
                (0..2).all(|i| {
                    at[i] >= chart.domain[i][0].to_f64_nearest()
                        && at[i] <= chart.domain[i][1].to_f64_nearest()
                })
            })
            .ok_or_else(|| bad("contact parameters leave the retained atlas cover"))?;
        Ok((index, chart.get(dec)?))
    }

    pub(super) fn contact_f64(
        &self,
        at: [f64; 2],
        dec: &mut Decider,
    ) -> Result<([f64; 3], [f64; 3]), Refusal> {
        match &self.geometry {
            Geometry::Atlas(charts) => {
                let (_, chart) = Self::chart(charts, at, dec)?;
                dec.charge_work(chart.evaluation_work())?;
                Ok(chart.contact_f64(at))
            }
            Geometry::Native {
                proposal: Some(proposal),
                ..
            } => {
                dec.charge_work(64)?;
                proposal.contact(at)
            }
            Geometry::Native {
                surface,
                proposal: None,
            } => {
                let parameters = [rational(at[0])?, rational(at[1])?];
                let (point, normal) = native_contact(surface, &parameters, 64, dec)?;
                Ok((
                    point.each_ref().map(|value| value.mid().to_f64_nearest()),
                    normal.each_ref().map(|value| value.mid().to_f64_nearest()),
                ))
            }
        }
    }

    pub(super) fn centre(&self, at: [f64; 2], dec: &mut Decider) -> Result<[f64; 3], Refusal> {
        let (point, normal) = self.contact_f64(at, dec)?;
        let d = self.distance64;
        let centre = std::array::from_fn(|i| point[i] + d * normal[i]);
        if centre.iter().all(|v| v.is_finite()) {
            Ok(centre)
        } else {
            Err(bad("non-finite offset evaluation"))
        }
    }

    pub(super) fn contact(
        &self,
        at: &[Rational; 2],
        dec: &mut Decider,
    ) -> Result<Contact, Refusal> {
        match &self.geometry {
            Geometry::Atlas(charts) => {
                let (index, chart) =
                    Self::chart(charts, at.each_ref().map(Rational::to_f64_nearest), dec)?;
                if !(0..2).all(|i| at[i] >= chart.domain[i][0] && at[i] <= chart.domain[i][1]) {
                    return Err(bad("rounded contact leaves its selected chart"));
                }
                let mut meter = Meter::default();
                let contact = chart.contact(
                    &at.each_ref().map(|t| Jet::constant(ball(t), 0)),
                    &mut meter,
                )?;
                meter.charge(dec)?;
                Ok(Contact {
                    chart: index,
                    point: contact
                        .point
                        .map(|jet| jet.0.into_iter().next().expect("order zero")),
                    normal: contact
                        .unit
                        .map(|jet| jet.0.into_iter().next().expect("order zero")),
                })
            }
            Geometry::Native { surface, .. } => {
                let (point, normal) = native_contact(surface, at, 128, dec)?;
                Ok(Contact {
                    chart: 0,
                    point,
                    normal,
                })
            }
        }
    }

    pub(super) fn checked(
        &self,
        at: &[Rational; 2],
        dec: &mut Decider,
    ) -> Result<(usize, [Ball; 3]), Refusal> {
        let contact = self.contact(at, dec)?;
        dec.charge_work(6)?;
        let distance = Ball::from_rational(&self.distance, contact.point[0].precision());
        Ok((
            contact.chart,
            std::array::from_fn(|i| &contact.point[i] + &(&distance * &contact.normal[i])),
        ))
    }
}

fn native_contact(
    surface: &Surface,
    at: &[Rational; 2],
    precision: u32,
    dec: &mut Decider,
) -> Result<([Ball; 3], [Ball; 3]), Refusal> {
    if matches!(surface, Surface::Cone(_)) && at[1].sign() != Sign::Positive {
        return Err(bad(
            "rolling-ball contact leaves the cone's positive-height chart",
        ));
    }
    let parameters = at
        .each_ref()
        .map(|value| Ball::from_rational(value, precision));
    let point = surface.evaluate(&parameters[0], &parameters[1], dec)?;
    let [du, dv] = surface.derivatives(&parameters[0], &parameters[1], dec)?;
    dec.charge_work(32)?;
    let cross: [Ball; 3] = std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &(&du[j] * &dv[k]) - &(&du[k] * &dv[j])
    });
    let square = &(&cross[0].square() + &cross[1].square()) + &cross[2].square();
    if square.lo().sign() != Sign::Positive {
        return Err(bad(
            "rolling-ball contact has an unresolved or singular native normal",
        ));
    }
    let inverse_length = square
        .sqrt()
        .and_then(|length| length.recip())
        .ok_or_else(|| bad("rolling-ball contact normal cannot be normalized"))?;
    Ok((point, cross.map(|value| &value * &inverse_length)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BSplineSurface;
    use catallax_numbers::{BudgetPolicy, with_operation_budget};

    #[test]
    fn lazy_atlas_preserves_eager_points_normals_and_seam_choice() {
        let ((), work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::new(BudgetPolicy {
                max_work: 100_000_000,
                ..Default::default()
            });
            let half = Rational::from_ratio(1, 2).unwrap();
            let mut controls = Vec::new();
            for u in [Rational::zero(), half.clone(), Rational::one()] {
                for v in [Rational::zero(), Rational::one()] {
                    controls.push(Point3::new(
                        u.clone().into(),
                        v.clone().into(),
                        (&u * &v).into(),
                    ));
                }
            }
            let source = BSplineSurface::new(
                [1, 1],
                [
                    vec![
                        Rational::zero(),
                        Rational::zero(),
                        half,
                        Rational::one(),
                        Rational::one(),
                    ],
                    [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
                ],
                [3, 2],
                controls,
                vec![Real::one(); 6],
                &mut dec,
            )
            .unwrap();
            let [atlas] = <[_; 1]>::try_from(
                RationalAtlas::partition(source, [false, false], &mut dec).unwrap(),
            )
            .unwrap();
            for distance in [-1, 0, 1] {
                let side = if distance < 0 {
                    Sign::Negative
                } else {
                    Sign::Positive
                };
                let eager = atlas
                    .charts()
                    .iter()
                    .map(|chart| Support::new(chart.patch(), side, &mut dec).unwrap())
                    .collect::<Vec<_>>();
                let lazy = NumericalFilletSupport::from_atlas(
                    &atlas,
                    Rational::from_i64(distance),
                    &mut dec,
                )
                .unwrap();
                for numerator in [0, 1, 2, 3, 4] {
                    let at = [
                        Rational::from_ratio(numerator, 4).unwrap(),
                        Rational::from_ratio(1, 3).unwrap(),
                    ];
                    let proposed = at.each_ref().map(Rational::to_f64_nearest);
                    let (index, expected) = eager
                        .iter()
                        .enumerate()
                        .find(|(_, chart)| {
                            (0..2).all(|i| {
                                proposed[i] >= chart.domain[i][0].to_f64_nearest()
                                    && proposed[i] <= chart.domain[i][1].to_f64_nearest()
                            })
                        })
                        .unwrap();
                    assert_eq!(
                        lazy.contact_f64(proposed, &mut dec).unwrap(),
                        expected.contact_f64(proposed)
                    );
                    let actual = lazy.contact(&at, &mut dec).unwrap();
                    assert_eq!(actual.chart, index);
                    let mut meter = Meter::default();
                    let expected = expected
                        .contact(
                            &at.each_ref().map(|value| Jet::constant(ball(value), 0)),
                            &mut meter,
                        )
                        .unwrap();
                    meter.charge(&mut dec).unwrap();
                    for axis in 0..3 {
                        assert_eq!(
                            (actual.point[axis].lo(), actual.point[axis].hi()),
                            (
                                expected.point[axis].0[0].lo(),
                                expected.point[axis].0[0].hi()
                            ),
                        );
                        assert_eq!(
                            (actual.normal[axis].lo(), actual.normal[axis].hi()),
                            (expected.unit[axis].0[0].lo(), expected.unit[axis].0[0].hi()),
                        );
                    }
                }
            }
        });
        assert!(!work.exhausted);
    }
}
