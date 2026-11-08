//! Exact regular normal extrusions of monotone polynomial Bezier spans.

use crate::{BSplineSurface, RPoint3, RVector3};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use std::sync::Arc;

#[derive(Debug)]
struct Data {
    controls: Vec<RPoint3>,
    direction: RVector3,
    carrier: BSplineSurface,
}

/// A finite quadratic or cubic Bezier span extruded normally to its plane.
/// Monotone projection proves injectivity; endpoint tangents prove regularity
/// where the interior projection derivative may tend to zero.
#[derive(Clone, Debug)]
pub struct BezierExtrusion {
    data: Arc<Data>,
}

impl BezierExtrusion {
    /// Checks coplanarity, normal extrusion, regularity and embedding exactly.
    /// A span whose chord projection cannot prove embedding is refused.
    pub fn new(
        controls: Vec<RPoint3>,
        direction: RVector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        if !(3..=4).contains(&controls.len()) {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-controls",
                "a quadratic or cubic span needs three or four controls",
            ));
        }
        if direction.is_zero() {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-direction",
                "the extrusion vector must be nonzero",
            ));
        }
        let first = &controls[0];
        for point in &controls[1..] {
            if !point.minus(first).dot(&direction).is_zero() {
                return Err(Refusal::new(
                    "shapes.bezier-extrusion-plane",
                    "the controls must lie in a plane normal to the extrusion vector",
                ));
            }
        }
        let last = controls.last().expect("checked control count");
        let chord = last.minus(first);
        if chord.is_zero()
            || controls[1].minus(first).is_zero()
            || last.minus(&controls[controls.len() - 2]).is_zero()
        {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-singular",
                "span endpoints and their first tangents must be distinct and nonzero",
            ));
        }
        // Nonnegative Bernstein coefficients and one positive coefficient make
        // the chord projection derivative positive throughout the open span.
        let slopes: Vec<_> = controls
            .windows(2)
            .map(|pair| pair[1].minus(&pair[0]).dot(&chord))
            .collect();
        if slopes.iter().any(|s| s.sign() == Sign::Negative) || slopes.iter().all(Rational::is_zero)
        {
            return Err(Refusal::new(
                "shapes.bezier-extrusion-embedding",
                "the span has no proved monotone chord projection",
            ));
        }
        let count = controls.len();
        let knots = std::iter::repeat_n(Rational::zero(), count)
            .chain(std::iter::repeat_n(Rational::one(), count))
            .collect();
        let net = controls
            .iter()
            .flat_map(|point| [point.to_real(), point.plus(&direction).to_real()])
            .collect();
        let carrier = BSplineSurface::new(
            [count - 1, 1],
            [
                knots,
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                ],
            ],
            [count, 2],
            net,
            vec![Real::one(); count * 2],
            dec,
        )?;
        Ok(Self {
            data: Arc::new(Data {
                controls,
                direction,
                carrier,
            }),
        })
    }

    /// The exact polynomial ruled carrier on the unit parameter square.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }

    /// The retained lower-boundary controls in parameter order.
    pub fn controls(&self) -> &[RPoint3] {
        &self.data.controls
    }

    /// The displacement from the lower boundary to the upper boundary.
    pub fn direction(&self) -> &RVector3 {
        &self.data.direction
    }
}
