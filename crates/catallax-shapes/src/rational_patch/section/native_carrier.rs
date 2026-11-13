//! Exact analytic incidence of a complete retained atlas chart cover.

use super::*;
use crate::RationalAtlas;
use std::collections::BTreeSet;

impl RationalAtlas {
    /// Prove that every retained chart lies on the analytic carrier. Positive
    /// spline weights make its homogeneous cut polynomial equivalent to the
    /// native equation. This proves inclusion, not a cover of the whole native
    /// surface or a bounded intersection window.
    pub fn certify_analytic_carrier(
        &self,
        surface: &SectionSurface,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        surface.validate(dec)?;
        let mut spans = BTreeSet::new();
        for chart in self.charts() {
            dec.check_deadline()?;
            let patch = chart.patch();
            if !patch.carrier().same_definition(self.source()) {
                return Err(Refusal::new(
                    "shapes.rational-atlas-carrier",
                    "the chart needs an exact binding to the retained tensor source",
                ));
            }
            for span in &patch.data.spans {
                if !spans.insert(span.indices) {
                    continue;
                }
                let equation = surface.equation(span);
                for (_, coefficient) in equation.terms() {
                    if exact_sign(&coefficient, dec)? != Sign::Zero {
                        return Err(Refusal::new(
                            "shapes.rational-atlas-carrier",
                            "a retained span's exact analytic cutting polynomial is nonzero",
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BSplineSurface;

    #[test]
    fn the_complete_span_identity_rejects_an_interior_bulge_on_a_planar_boundary() {
        let plane = SectionSurface::Plane(Plane {
            origin: Point3::new(Real::zero(), Real::zero(), Real::zero()),
            normal: Vector3::new(Real::zero(), Real::zero(), Real::one()),
            ref_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        });
        for bulged in [false, true] {
            let mut dec = Decider::default();
            let controls = (0..3)
                .flat_map(|u| {
                    (0..3).map(move |v| {
                        Point3::new(
                            Real::from_i64(u),
                            Real::from_i64(v),
                            Real::from_i64(i64::from(bulged && u == 1 && v == 1)),
                        )
                    })
                })
                .collect();
            let knots = [0, 0, 0, 1, 1, 1].map(Rational::from_i64).to_vec();
            let source = BSplineSurface::new(
                [2, 2],
                [knots.clone(), knots],
                [3, 3],
                controls,
                vec![Real::one(); 9],
                &mut dec,
            )
            .unwrap();
            let atlases =
                RationalAtlas::partition_with_polynomial_bounds(source, [false, false], &mut dec)
                    .unwrap();
            assert_eq!(atlases.len(), 1);
            let proof = atlases[0].certify_analytic_carrier(&plane, &mut dec);
            if bulged {
                assert_eq!(proof.unwrap_err().code, "shapes.rational-atlas-carrier");
            } else {
                proof.unwrap();
            }
        }
    }
}
