//! Volume from retained finite field cells, independently of source topology.

use super::{Enclosure, QuantityAccuracy, reached, valid};
use crate::{Instruction, Plan, draft::Cell, geometry::Leaf, metric::Kind};
use catallax_numbers::{Ball, Decider, Real, Refusal, Sign};
use catallax_shapes::Point3;

/// Encloses the volume of one retained planar or axial analytic field cell, or
/// of one sphere or torus primitive.
/// Other definitions keep the general `properties` query and its cell controls.
pub fn volume(
    plan: &Plan,
    accuracy: &QuantityAccuracy,
    dec: &mut Decider,
) -> Result<Enclosure, Refusal> {
    if !valid(accuracy) {
        return Err(Refusal::new(
            "field.mass-options",
            "nonnegative volume Accuracy is required",
        ));
    }
    let unsupported = || {
        Refusal::new(
            "field.volume-definition",
            "analytic volume needs one retained planar or axial field cell, sphere or torus; use mass::properties for general field integration",
        )
    };
    let Instruction::Leaf(id) = plan.instructions[plan.root] else {
        return Err(unsupported());
    };
    let region = match &plan.leaves[id] {
        Leaf::Draft(region) => region,
        Leaf::Primitive(primitive) => {
            let [a, b, _] = primitive.parameters.clone().map(Real::from);
            let (factor, pi) = match primitive.kind {
                Kind::Sphere => (Real::from_i64(4) * a.pow_u32(3) / Real::from_i64(3), 1),
                Kind::Torus => (Real::from_i64(2) * a * b.square(), 2),
                _ => return Err(unsupported()),
            };
            return enclose(factor, pi, accuracy, dec);
        }
        _ => return Err(unsupported()),
    };
    let (factor, pi) = match &region.cell {
        Cell::Round(cell) => {
            let a = &cell.radii[0];
            let b = &cell.radii[1];
            (
                (&cell.heights[1] - &cell.heights[0]) * (a.square() + a * b + b.square())
                    / Real::from_i64(3),
                1,
            )
        }
        Cell::Polyhedral(cell) => {
            let mut six = Real::zero();
            for face in &cell.faces {
                for cycle in &face.loops {
                    if cycle.len() < 3 {
                        continue;
                    }
                    let a = cell.vertices[cycle[0]].minus(&Point3::origin());
                    for pair in cycle[1..].windows(2) {
                        let b = cell.vertices[pair[0]].minus(&Point3::origin());
                        let c = cell.vertices[pair[1]].minus(&Point3::origin());
                        six = six + a.dot(&b.cross(&c));
                    }
                }
            }
            (six / Real::from_i64(6), 0)
        }
        _ => return Err(unsupported()),
    };
    enclose(factor, pi, accuracy, dec)
}

/// Encloses `factor * pi^pi_power` to the requested accuracy.
fn enclose(
    factor: Real,
    pi_power: u32,
    accuracy: &QuantityAccuracy,
    dec: &mut Decider,
) -> Result<Enclosure, Refusal> {
    if crate::geometry::sign(&factor, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "field.volume-orientation",
            "retained field cell has no positive oriented volume",
        ));
    }
    let maximum = dec.policy().max_precision_bits;
    if maximum == 0 {
        return Err(Refusal::new(
            "numbers.precision-limit",
            "no volume enclosure precision available",
        ));
    }
    let mut bits = maximum.min(128);
    loop {
        dec.check_deadline()?;
        let mut value = crate::geometry::enclose(&factor, bits, dec)?;
        for _ in 0..pi_power {
            value = value * Ball::pi_checked(bits, dec)?;
        }
        let value = Enclosure {
            lo: value.lo(),
            hi: value.hi(),
        };
        if reached(&value, accuracy) {
            return Ok(value);
        }
        if bits == maximum {
            return Err(Refusal::new(
                "field.volume-accuracy",
                format!(
                    "volume width {} exceeds the requested accuracy at {bits} bits",
                    value.width()
                ),
            ));
        }
        bits = bits.saturating_mul(2).min(maximum);
    }
}
