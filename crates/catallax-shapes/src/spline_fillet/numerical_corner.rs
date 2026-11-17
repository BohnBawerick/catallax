//! A common rolling ball for three incident fillets. Like SMLib's spherical
//! corner path and OCCT's special-case corners, the ball is retained only
//! after its contacts meet the declared numerical tolerance.

use super::*;

/// A numerically located common ball, not a completed trimmed corner face.
#[derive(Clone, Debug)]
pub struct NumericalFilletCorner {
    centre: [Rational; 3],
    radius: Rational,
    parameters: [[Rational; 2]; 3],
    contacts: [[Rational; 3]; 3],
    tolerance_mm: Rational,
    sampled_error_mm: Rational,
}

impl NumericalFilletCorner {
    /// Intersects three signed offsets near the supplied native contacts.
    /// Distances must have one common positive magnitude. The Newton proposal
    /// is checked against interval contacts and a nonsingular normal triple.
    pub fn solve(
        supports: &[NumericalFilletSupport; 3],
        seed: [[Rational; 2]; 3],
        tolerance_mm: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let radius = supports[0].distance.abs();
        let tolerance = tolerance_mm.to_f64_nearest();
        if radius.sign() != Sign::Positive
            || supports
                .iter()
                .any(|support| support.distance.abs() != radius)
            || tolerance_mm.sign() != Sign::Positive
            || !tolerance.is_finite()
            || tolerance <= 0.0
        {
            return Err(bad(
                "corner needs equal nonzero ball radii and a finite positive tolerance",
            ));
        }
        let x = walk::newton(
            std::array::from_fn::<_, 6, _>(|i| seed[i / 2][i % 2].to_f64_nearest()),
            tolerance,
            |x, dec| {
                let a = supports[0].centre([x[0], x[1]], dec)?;
                let b = supports[1].centre([x[2], x[3]], dec)?;
                let c = supports[2].centre([x[4], x[5]], dec)?;
                Ok([
                    a[0] - b[0],
                    a[1] - b[1],
                    a[2] - b[2],
                    a[0] - c[0],
                    a[1] - c[1],
                    a[2] - c[2],
                ])
            },
            dec,
        )?;
        let parameters = [
            [rational(x[0])?, rational(x[1])?],
            [rational(x[2])?, rational(x[3])?],
            [rational(x[4])?, rational(x[5])?],
        ];
        let contacts = [
            supports[0].contact(&parameters[0], dec)?,
            supports[1].contact(&parameters[1], dec)?,
            supports[2].contact(&parameters[2], dec)?,
        ];
        let offsets: [[Ball; 3]; 3] = std::array::from_fn(|i| {
            let contact = &contacts[i];
            let distance = Ball::from_rational(&supports[i].distance, contact.point[0].precision());
            std::array::from_fn(|j| &contact.point[j] + &(&distance * &contact.normal[j]))
        });
        dec.charge_work(32)?;
        let normals = contacts.each_ref().map(|contact| &contact.normal);
        let determinant = (0..3).fold(Ball::from_rational(&Rational::zero(), 128), |sum, i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            let minor = &(&normals[1][j] * &normals[2][k]) - &(&normals[1][k] * &normals[2][j]);
            sum + (&normals[0][i] * &minor)
        });
        if determinant.contains_zero() {
            return Err(bad(
                "corner contacts do not isolate a nonsingular common ball",
            ));
        }
        let centre = offsets[0].each_ref().map(Ball::mid);
        let error = offset_error(&centre, [&offsets[0], &offsets[1]])
            .max(offset_error(&centre, [&offsets[0], &offsets[2]]));
        if error > tolerance_mm.mul_pow2(-3) {
            return Err(bad(
                "corner's checked offset residual exceeds its numerical tolerance",
            ));
        }
        Ok(Self {
            centre,
            radius,
            parameters,
            contacts: std::array::from_fn(|i| contacts[i].point.each_ref().map(Ball::mid)),
            tolerance_mm,
            sampled_error_mm: error,
        })
    }

    /// Common centre checked at the three retained contact stations.
    pub fn centre(&self) -> &[Rational; 3] {
        &self.centre
    }
    /// Requested rolling-ball radius in millimetres.
    pub fn radius(&self) -> &Rational {
        &self.radius
    }
    /// Native contact parameters on the original three supports.
    pub fn parameters(&self) -> &[[Rational; 2]; 3] {
        &self.parameters
    }
    /// Numerical source contacts. They are not exact points of a new B-rep.
    pub fn contacts(&self) -> &[[Rational; 3]; 3] {
        &self.contacts
    }
    /// Declared numerical tolerance in millimetres.
    pub fn tolerance_mm(&self) -> &Rational {
        &self.tolerance_mm
    }
    /// Enclosed centre-to-offset disagreement at the three checked contacts.
    pub fn sampled_error_mm(&self) -> &Rational {
        &self.sampled_error_mm
    }
}
