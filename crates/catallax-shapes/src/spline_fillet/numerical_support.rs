//! Signed offsets used by numerical rolling-ball continuation.

use super::*;
use crate::{SplineSupport, Surface};

#[derive(Debug)]
enum Geometry {
    Atlas(Vec<Support>),
    Native(Surface),
}

/// A support's native parameterization and signed rolling-ball distance.
/// Binary64 values propose contacts. Interval evaluations check accepted ones.
#[derive(Debug)]
pub struct NumericalFilletSupport {
    geometry: Geometry,
    pub(super) distance: Rational,
}

pub(super) struct Contact {
    pub(super) chart: usize,
    pub(super) point: [Ball; 3],
    pub(super) normal: [Ball; 3],
}

impl NumericalFilletSupport {
    /// Retains the surface's native parameters. A negative distance places the
    /// ball opposite `S_u cross S_v`; it does not reverse the parameterization.
    pub fn new(
        surface: Surface,
        signed_distance: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if signed_distance.sign() == Sign::Zero
            || !signed_distance.to_f64_nearest().is_finite()
        {
            return Err(bad("rolling-ball support needs a finite nonzero distance"));
        }
        if let Surface::BSpline(SplineSupport::RationalAtlas(atlas)) = &surface {
            return Self::from_atlas(atlas, signed_distance, dec);
        }
        Ok(Self {
            geometry: Geometry::Native(surface),
            distance: signed_distance,
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
            .map(|chart| Support::new(chart.patch(), side, dec))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            geometry: Geometry::Atlas(charts),
            distance,
        })
    }

    fn chart(charts: &[Support], at: [f64; 2]) -> Result<(usize, &Support), Refusal> {
        charts
            .iter()
            .enumerate()
            .find(|(_, chart)| {
                (0..2).all(|i| {
                    at[i] >= chart.domain[i][0].to_f64_nearest()
                        && at[i] <= chart.domain[i][1].to_f64_nearest()
                })
            })
            .ok_or_else(|| bad("contact parameters leave the retained atlas cover"))
    }

    pub(super) fn contact_f64(
        &self,
        at: [f64; 2],
        dec: &mut Decider,
    ) -> Result<([f64; 3], [f64; 3]), Refusal> {
        match &self.geometry {
            Geometry::Atlas(charts) => {
                let (_, chart) = Self::chart(charts, at)?;
                dec.charge_work(chart.evaluation_work())?;
                Ok(chart.contact_f64(at))
            }
            Geometry::Native(surface) => {
                let parameters = [rational(at[0])?, rational(at[1])?];
                let (point, normal) = native_contact(surface, &parameters, 64, dec)?;
                Ok((
                    point.each_ref().map(|value| value.mid().to_f64_nearest()),
                    normal.each_ref().map(|value| value.mid().to_f64_nearest()),
                ))
            }
        }
    }

    pub(super) fn centre(
        &self,
        at: [f64; 2],
        dec: &mut Decider,
    ) -> Result<[f64; 3], Refusal> {
        let (point, normal) = self.contact_f64(at, dec)?;
        let d = self.distance.to_f64_nearest();
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
                    Self::chart(charts, at.each_ref().map(Rational::to_f64_nearest))?;
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
                    point: contact.point.map(|jet| jet.0.into_iter().next().expect("order zero")),
                    normal: contact.unit.map(|jet| jet.0.into_iter().next().expect("order zero")),
                })
            }
            Geometry::Native(surface) => {
                let (point, normal) = native_contact(surface, at, 128, dec)?;
                Ok(Contact { chart: 0, point, normal })
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
    let parameters = at.each_ref().map(|value| Ball::from_rational(value, precision));
    let point = surface.evaluate(&parameters[0], &parameters[1], dec)?;
    let [du, dv] = surface.derivatives(&parameters[0], &parameters[1], dec)?;
    dec.charge_work(32)?;
    let cross: [Ball; 3] = std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &(&du[j] * &dv[k]) - &(&du[k] * &dv[j])
    });
    let square = &(&cross[0].square() + &cross[1].square()) + &cross[2].square();
    if square.lo().sign() != Sign::Positive {
        return Err(bad("rolling-ball contact has an unresolved or singular native normal"));
    }
    let inverse_length = square
        .sqrt()
        .and_then(|length| length.recip())
        .ok_or_else(|| bad("rolling-ball contact normal cannot be normalized"))?;
    Ok((point, cross.map(|value| &value * &inverse_length)))
}
