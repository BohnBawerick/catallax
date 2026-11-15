//! Exact finite tensor geometry for STEP extrusions of rational spline curves.

use super::*;
use crate::exact::dot;

fn range(values: impl IntoIterator<Item = Rational>) -> Option<[Rational; 2]> {
    let mut values = values.into_iter();
    let first = values.next()?;
    Some(values.fold([first.clone(), first], |[lo, hi], value| {
        [lo.min(value.clone()), hi.max(value)]
    }))
}

fn curve_range(edge: &OrientedEdge, direction: &Vec3) -> [Rational; 2] {
    match &edge.curve {
        Curve::Line(..) => range([
            dot(&edge.start.point, direction),
            dot(&edge.end.point, direction),
        ])
        .expect("two endpoints"),
        Curve::BSpline(curve) => range(curve.geometry.controls().iter().map(|p| dot(p, direction)))
            .expect("validated source controls"),
        Curve::Circle(placement, radius) | Curve::Ellipse(placement, radius, _) => {
            let radius = if let Curve::Ellipse(_, a, b) = &edge.curve {
                a.max(b)
            } else {
                radius
            };
            // The L1 norm bounds the Euclidean support radius in every frame.
            let spread = radius
                * direction
                    .iter()
                    .fold(Rational::zero(), |sum, component| sum + component.abs());
            let center = dot(&placement.origin, direction);
            [&center - &spread, center + spread]
        }
    }
}

impl Reader<'_> {
    pub(super) fn linear_extrusion(
        &self,
        id: u64,
        loops: &[FaceLoop],
    ) -> Result<Surface, ImportRefusal> {
        let args = self.args(id, "SURFACE_OF_LINEAR_EXTRUSION")?;
        let profile = self.curve(reference(args, 1)?)?;
        let Curve::BSpline(profile) = profile else {
            return Err(unsupported(id, self.entity(id)?));
        };
        let axis = self.args(reference(args, 2)?, "VECTOR")?;
        let mut direction = self.direction(reference(axis, 1)?)?;
        let magnitude = self.length(axis, 2)?;
        match magnitude.sign() {
            catallax_numbers::Sign::Zero => {
                return Err(malformed("extrusion vector must not be zero"));
            }
            catallax_numbers::Sign::Negative => direction = direction.map(|v| -v),
            catallax_numbers::Sign::Positive => {}
        }
        let closed = profile.closed;
        let profile = profile.geometry;
        let [profile_lo, profile_hi] = range(
            profile
                .controls()
                .iter()
                .map(|point| dot(point, &direction)),
        )
        .expect("validated source controls");
        let boundary = loops.iter().flat_map(|bound| {
            bound
                .edges
                .iter()
                .flat_map(|edge| curve_range(edge, &direction))
                .chain(bound.vertex.iter().map(|v| dot(&v.point, &direction)))
        });
        let [boundary_lo, boundary_hi] =
            range(boundary).ok_or_else(|| ImportRefusal::NotSolid {
                detail: format!("extrusion #{id} has no finite face boundary"),
            })?;
        // Direction ratios rescale v positively without changing the carrier or its normal.
        let norm_sq = dot(&direction, &direction);
        let lo = (boundary_lo - profile_hi) / &norm_sq;
        let hi = (boundary_hi - profile_lo) / norm_sq;
        if lo >= hi {
            return Err(malformed(format!(
                "extrusion #{id} has no finite surface interval"
            )));
        }
        // Positive weights bound the whole profile and each whole trim curve.
        // Projection along the sweep is monotone, so its finite face lies in this strip.
        let controls = profile
            .controls()
            .iter()
            .flat_map(|point| {
                [&lo, &hi].map(|t| std::array::from_fn(|i| &point[i] + &(t * &direction[i])))
            })
            .collect();
        let geometry = RationalSplineSurface::new(
            [profile.degree(), 1],
            [
                profile.knots().to_vec(),
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                ],
            ],
            [profile.controls().len(), 2],
            controls,
            profile
                .weights()
                .iter()
                .flat_map(|weight| [weight.clone(), weight.clone()])
                .collect(),
            self.dec,
        )
        .map_err(spline::source_error)?;
        Ok(Surface::BSpline(StepSplineSurface {
            geometry,
            closed: [closed, Some(false)],
            self_intersect: None,
        }))
    }
}
