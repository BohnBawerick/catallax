//! Exact offset intent, independent of the eventual B-rep carrier.

use catallax_numbers::{Decider, Real, Refusal, Sign};

use crate::{EllipticCylinder, Pose, Surface, math};

/// Chosen side of the outward-oriented elliptic-cylinder support.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OffsetSide {
    /// Positive distance along the original outward normal.
    Outward,
    /// Negative distance along the original outward normal.
    Inward,
}

/// Immutable definition of an elliptic-cylinder support offset.
///
/// This is intent, not a fitted carrier, valid trimmed face or distance
/// certificate. Inward offsets still need a regular-domain and branch proof.
/// The original support and unsigned distance match the field measuring stick's
/// `FaceDefinition::Offset` fields; the side must also survive in the producer.
#[derive(Clone, Debug)]
pub struct EllipticOffsetDefinition {
    original: Surface,
    distance: Real,
    side: OffsetSide,
}

impl EllipticOffsetDefinition {
    /// Validates the original analytic definition and signed offset expression.
    /// Negative distance selects the inward branch; zero retains the exact support.
    pub fn new(
        original: EllipticCylinder,
        signed_distance: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let original = EllipticCylinder::new(
            original.origin,
            original.axis,
            original.major_dir,
            original.major_sq,
            original.minor_sq,
            dec,
        )?;
        math::validate(&signed_distance, dec)?;
        let (side, distance) = match math::sign(&signed_distance, dec)? {
            Sign::Negative => (OffsetSide::Inward, -signed_distance),
            Sign::Zero | Sign::Positive => (OffsetSide::Outward, signed_distance),
        };
        dec.check_deadline()?;
        Ok(Self {
            original: Surface::EllipticCylinder(original),
            distance,
            side,
        })
    }

    /// Snapshot of the exact original support, not a rounded or fitted surface.
    pub fn original(&self) -> &Surface {
        &self.original
    }

    /// Nonnegative distance in millimetres, matching the field checker input.
    /// Unsigned distance alone does not identify the selected side.
    pub fn distance(&self) -> &Real {
        &self.distance
    }

    /// Branch chosen relative to the original support's outward normal.
    pub fn side(&self) -> OffsetSide {
        self.side
    }

    /// Rigidly transports the original law and selected side without refitting.
    pub fn transformed(&self, pose: &Pose) -> Self {
        Self {
            original: self.original.transformed(pose),
            distance: self.distance.clone(),
            side: self.side,
        }
    }
}
