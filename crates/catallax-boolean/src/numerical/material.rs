//! Numerical material-side samples on retained native curves.

use super::*;

pub(super) struct CylinderSide {
    origin: [Ball; 3],
    axis: [Ball; 3],
    radius_squared: Ball,
}
impl CylinderSide {
    pub fn new(tool: &Cylinder, dec: &mut Decider) -> Result<Self, Refusal> {
        let axis = unit(&tool.axis, dec)?;
        let mut enclosed = |v: [&Real; 3]| -> Result<[Ball; 3], Refusal> {
            Ok([
                interval(v[0], dec)?,
                interval(v[1], dec)?,
                interval(v[2], dec)?,
            ])
        };
        let origin = enclosed([&tool.origin.x, &tool.origin.y, &tool.origin.z])?;
        let axis = enclosed([&axis.x, &axis.y, &axis.z])?;
        Ok(Self {
            origin,
            axis,
            radius_squared: interval(&tool.radius_sq, dec)?,
        })
    }
    pub fn covers(
        &self,
        source: &BSplineSurface,
        parameters: &[Rational; 2],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        for control in source.controls() {
            dec.charge_work(1)?;
            let point = [
                interval(&control.x, dec)?,
                interval(&control.y, dec)?,
                interval(&control.z, dec)?,
            ];
            let offset: [Ball; 3] = std::array::from_fn(|i| &point[i] - &self.origin[i]);
            let height = &offset[0] * &self.axis[0]
                + &offset[1] * &self.axis[1]
                + &offset[2] * &self.axis[2];
            if height.lo() < parameters[0] || height.hi() > parameters[1] {
                return Err(bad(
                    "this through-hole consumer requires a tool interval covering the retained cone control hull",
                ));
            }
        }
        Ok(())
    }
    fn radial(&self, point: &[Ball; 3]) -> [Ball; 3] {
        let d: [Ball; 3] = std::array::from_fn(|i| &point[i] - &self.origin[i]);
        let height = &d[0] * &self.axis[0] + &d[1] * &self.axis[1] + &d[2] * &self.axis[2];
        std::array::from_fn(|i| &d[i] - &(&height * &self.axis[i]))
    }
    pub fn outside(
        &self,
        source: &BSplineSurface,
        uv: &[Rational; 2],
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        dec.charge_work(1)?;
        let point = source.evaluate_box([&uv[0], &uv[0]], [&uv[1], &uv[1]], BITS, dec)?;
        let radial = self.radial(&point);
        let signed =
            radial[0].square() + radial[1].square() + radial[2].square() - &self.radius_squared;
        if signed.lo().sign() == Sign::Positive {
            Ok(true)
        } else if signed.hi().sign() == Sign::Negative {
            Ok(false)
        } else {
            Err(bad(
                "finite boundary material sample needs curve-contact interval ownership",
            ))
        }
    }
    pub fn outside_left(
        &self,
        source: &BSplineSurface,
        uv: &BSplineCurve,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut orientation = None;
        for (_, span) in uv.spans() {
            dec.charge_work(1)?;
            let t = (span[0] + span[1]) / q(2);
            let at = uv.point_at(&Real::from(t.clone()), dec)?;
            let [u, v, _] = events::rational_point(at)?;
            let jet = source.jet_box([&u, &u], [&v, &v], BITS, dec)?;
            let derivative = uv.first_at(&Real::from(t), dec)?;
            let dx = interval(&derivative.x, dec)?;
            let dy = interval(&derivative.y, dec)?;
            let radial = self.radial(&jet.value);
            let gradient = jet
                .first
                .each_ref()
                .map(|d| &radial[0] * &d[0] + &radial[1] * &d[1] + &radial[2] * &d[2]);
            let left = &gradient[1] * &dx - &gradient[0] * &dy;
            let sign = if left.lo().sign() == Sign::Positive {
                true
            } else if left.hi().sign() == Sign::Negative {
                false
            } else {
                return Err(bad(
                    "a fitted native material-side sample is tangent or unresolved",
                ));
            };
            if orientation.is_some_and(|previous| previous != sign) {
                return Err(bad(
                    "fitted material-side samples change orientation inside a retained piece",
                ));
            }
            orientation = Some(sign);
        }
        orientation.ok_or_else(|| bad("a fitted native curve has no finite nonzero span"))
    }
}
