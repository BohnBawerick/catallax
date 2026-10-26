//! Exact angular domains for definition-native revolutions.

use crate::{Location, geometry, profile};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign, sin_cos_degrees};
use catallax_shapes::{RVector3, Vector3};

#[derive(Clone, Debug)]
pub(super) struct Angle {
    degrees: Rational,
    pub(super) cosine: Real,
    pub(super) sine: Real,
}

/// Preserve rational dimensions and retain other square roots exactly.
pub(super) fn radius(squared: &Rational) -> Real {
    profile::sqrt(squared)
        .map(Real::from)
        .unwrap_or_else(|| Real::from(squared.clone()).sqrt())
}

impl Angle {
    pub(super) fn new(degrees: &Rational, dec: &mut Decider) -> Result<Self, Refusal> {
        let (cosine, sine) = if *degrees == Rational::from_i64(360) {
            (Real::one(), Real::zero())
        } else if *degrees == Rational::from_i64(90) {
            (Real::zero(), Real::one())
        } else if *degrees == Rational::from_i64(180) {
            (-Real::one(), Real::zero())
        } else if *degrees == Rational::from_i64(270) {
            (Real::zero(), -Real::one())
        } else {
            let [sine, cosine] = sin_cos_degrees(degrees, dec)?;
            (cosine, sine)
        };
        Ok(Self {
            degrees: degrees.clone(),
            cosine,
            sine,
        })
    }

    pub(super) fn is(&self, degrees: i64) -> bool {
        self.degrees == Rational::from_i64(degrees)
    }

    pub(super) fn fast(&self) -> (u16, [f64; 2]) {
        let rotation = [self.cosine.to_f64(), self.sine.to_f64()];
        for cardinal in [90, 180, 270, 360] {
            if self.is(cardinal) {
                return (cardinal as u16, rotation);
            }
        }
        // Keep the existing cardinal fast paths; 1 and 2 select the new
        // non-cardinal convex and reflex angular predicates respectively.
        let sector = if self.degrees <= Rational::from_i64(180) {
            1
        } else {
            2
        };
        (sector, rotation)
    }

    // b is the coefficient along z cross axis, so its physical length is
    // b*sqrt(norm). Each flip is a distinct preimage in the sketch plane.
    // Reuse settled coordinate signs to preserve the cardinal query work path.
    pub(super) fn location(
        &self,
        flip: i64,
        coordinates: [&Real; 2],
        norm: &Rational,
        signs: [Sign; 2],
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let [b, z] = coordinates;
        if self.is(360) {
            return Ok(Location::Inside);
        }
        for cardinal in [90, 180, 270] {
            if self.is(cardinal) {
                return Ok(super::angular(cardinal as u16, flip, signs[0], signs[1]));
            }
        }
        let oriented = |sign| {
            if flip == 1 {
                sign
            } else {
                match sign {
                    Sign::Positive => Sign::Negative,
                    Sign::Negative => Sign::Positive,
                    Sign::Zero => Sign::Zero,
                }
            }
        };
        let start = oriented(signs[1]);
        let end = oriented(geometry::sign(
            &((b * &radius(norm)) * self.sine.clone() - z * &self.cosine),
            dec,
        )?);
        let outside = if self.degrees <= Rational::from_i64(180) {
            start == Sign::Negative || end == Sign::Negative
        } else {
            start == Sign::Negative && end == Sign::Negative
        };
        Ok(if outside {
            Location::Outside
        } else if start == Sign::Zero || end == Sign::Zero {
            // For a reflex sweep, the opposite half of each cap plane is
            // interior. Only its directed ray is part of the angular boundary.
            if self.degrees > Rational::from_i64(180)
                && (start == Sign::Positive || end == Sign::Positive)
            {
                Location::Inside
            } else {
                Location::Boundary
            }
        } else {
            Location::Inside
        })
    }

    pub(super) fn end_normal(&self, axis: &RVector3) -> Vector3 {
        let z = RVector3::new(Rational::zero(), Rational::zero(), Rational::one());
        if self.sine.rational_ref().is_some_and(Rational::is_zero) {
            return z.to_real().scaled(&self.cosine);
        }
        z.to_real().scaled(&self.cosine).plus(
            &axis
                .cross(&z)
                .to_real()
                .scaled(&(&self.sine / &radius(&axis.norm_squared()))),
        )
    }

    pub(super) fn rational_rotation(&self, norm: &Rational) -> Option<(Rational, Rational)> {
        let cosine = self.cosine.rational_ref()?.clone();
        let sine = self.sine.rational_ref()?;
        let transverse = if sine.is_zero() {
            Rational::zero()
        } else {
            sine / &profile::sqrt(norm)?
        };
        Some((cosine, transverse))
    }

    pub(super) fn unrotate_radial(
        &self,
        radial: &Vector3,
        axis: &Vector3,
        norm: &Rational,
    ) -> Vector3 {
        if self.sine.rational_ref().is_some_and(Rational::is_zero) {
            return radial.scaled(&self.cosine);
        }
        radial
            .scaled(&self.cosine)
            .plus(&axis.cross(radial).scaled(&(-&self.sine / radius(norm))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn algebraic_cap_native_display_keeps_its_chart_coverage_refusal() {
        let recipe = catallax_recipe::Recipe::parse(include_str!(
            "../../../catallax-features/tests/fixtures/algebraic-revolutions/line-30.json"
        ))
        .unwrap();
        let plan = crate::Plan::compile(&recipe).unwrap();
        let refusal =
            crate::display::mesh_plan(&plan, &Default::default(), &mut Decider::default())
                .unwrap_err();
        assert_eq!(refusal.code, "field.display-unsupported");
    }

    #[test]
    fn directed_angular_domain_keeps_both_sides_and_reflex_caps_exact() {
        let norm = Rational::from_i64(2);
        for degrees in [30, 45, 60, 120, 210, 225, 300, 330] {
            let mut dec = Decider::default();
            let angle = Angle::new(&Rational::from_i64(degrees), &mut dec).unwrap();
            let half = Angle::new(&Rational::from_ratio(degrees, 2).unwrap(), &mut dec).unwrap();
            for flip in [-1, 1] {
                let f = Real::from_i64(flip);
                for (cosine, sine, want) in [
                    (Real::one(), Real::zero(), Location::Boundary),
                    (angle.cosine.clone(), angle.sine.clone(), Location::Boundary),
                    (half.cosine.clone(), half.sine.clone(), Location::Inside),
                    (
                        Real::one(),
                        -Real::from(Rational::from_ratio(1, 100).unwrap()),
                        Location::Outside,
                    ),
                ] {
                    let b = (&f * &cosine) / radius(&norm);
                    let z = &f * &sine;
                    let signs = [
                        geometry::sign(&b, &mut dec).unwrap(),
                        geometry::sign(&z, &mut dec).unwrap(),
                    ];
                    assert_eq!(
                        angle
                            .location(flip, [&b, &z], &norm, signs, &mut dec)
                            .unwrap(),
                        want,
                        "{degrees} {flip}"
                    );
                }
            }
        }
    }
}
