//! Canonical feature definitions and their stock-origin correspondence.

use std::collections::BTreeMap;

use crate::{
    Axis, ElementId, Loop, Number, Operation, Origin, Point2, Point3, Recipe, Role, Segment, Step,
    StepId, Turn, Vector3,
};

/// A stock solid expressed as a sketch and one sweep, without changing its identity.
/// Generic authored revolutions still name radial segments as sides. Only this
/// explicit expansion maps a cone's radial meridians to its stock cap names.
#[derive(Clone, Debug)]
pub struct CanonicalStock {
    recipe: Recipe,
    origins: BTreeMap<Origin, Origin>,
}

impl CanonicalStock {
    /// Expands a stock operation. Non-stock operations return `None`.
    /// Dimension and placement validity remain the evaluators' responsibility.
    pub fn new(operation: &Operation, step: &StepId) -> Option<Self> {
        let id = |text: &str| ElementId::new(text).expect("fixed canonical segment id");
        let p = |x: Number, y: Number| Point2([x, y]);
        let line = |name: &str, x: Number, y: Number| Segment::Line {
            id: id(name),
            to: p(x, y),
        };
        let circle = |x: Number, radius: Number| Segment::Circle {
            id: id("wall"),
            centre: p(x, 0.into()),
            radius,
        };
        let zero = Number::from(0);
        let profile = StepId::new(if step.as_str() == "canonical_profile" {
            "canonical_section"
        } else {
            "canonical_profile"
        })
        .unwrap();
        let (pose, loops, distance) = match operation {
            Operation::Block { size, pose } => (
                pose,
                vec![Loop(vec![
                    line("y_min", size[0].clone(), zero.clone()),
                    line("x_max", size[0].clone(), size[1].clone()),
                    line("y_max", zero.clone(), size[1].clone()),
                    line("x_min", zero.clone(), zero.clone()),
                ])],
                Some(size[2].clone()),
            ),
            Operation::Cylinder {
                radius,
                height,
                pose,
            } => (
                pose,
                vec![Loop(vec![circle(zero.clone(), radius.clone())])],
                Some(height.clone()),
            ),
            Operation::Sphere { radius, pose } => (
                pose,
                vec![Loop(vec![
                    Segment::Arc {
                        id: id("wall"),
                        centre: p(zero.clone(), zero.clone()),
                        to: p(zero.clone(), radius.clone()),
                        turn: Turn::Ccw,
                    },
                    line("axis", zero.clone(), Number::new(-radius.value())),
                ])],
                None,
            ),
            Operation::Torus { major, minor, pose } => (
                pose,
                vec![Loop(vec![circle(major.clone(), minor.clone())])],
                None,
            ),
            Operation::Cone {
                base_radius,
                top_radius,
                height,
                pose,
            } => {
                let mut segments = Vec::new();
                if !base_radius.value().is_zero() {
                    segments.push(line("start", base_radius.clone(), zero.clone()));
                }
                segments.push(line("wall", top_radius.clone(), height.clone()));
                if !top_radius.value().is_zero() {
                    segments.push(line("end", zero.clone(), height.clone()));
                }
                segments.push(line("axis", zero.clone(), zero.clone()));
                (pose, vec![Loop(segments)], None)
            }
            _ => return None,
        };
        let mut plane = pose.clone();
        let sweep = if let Some(distance) = distance {
            Operation::Extrude {
                profile: profile.clone(),
                distance,
                direction: None,
            }
        } else {
            // Compose the stock frame with +90 degrees about local x. The
            // meridian's local y is the stock z axis. Quaternions are homogeneous.
            let [w, x, y, z] = pose.rotation.each_ref().map(|n| n.value());
            plane.rotation = [w - x, w + x, y + z, z - y].map(Number::new);
            let two = catallax_numbers::Rational::from_i64(2);
            let norm = w * w + x * x + y * y + z * z;
            if norm.is_zero() {
                // Leave the invalid quaternion to the normal placement check.
                plane.rotation = pose.rotation.clone();
            }
            let axis = if norm.is_zero() {
                Vector3([zero.clone(), zero.clone(), zero.clone()])
            } else {
                Vector3([
                    Number::new(&(&two * &(x * z + w * y)) / &norm),
                    Number::new(&(&two * &(y * z - w * x)) / &norm),
                    Number::new(&(w * w - x * x - y * y + z * z) / &norm),
                ])
            };
            Operation::Revolve {
                profile: profile.clone(),
                axis: Axis {
                    point: Point3(pose.translation.0.clone()),
                    direction: axis,
                },
                angle: 360.into(),
            }
        };
        let mut origins = BTreeMap::new();
        if matches!(operation, Operation::Cone { .. }) {
            for (segment, role) in [("start", Role::Start), ("end", Role::End)] {
                origins.insert(
                    Origin {
                        step: step.clone(),
                        role: Role::Side {
                            segment: id(segment),
                            path: None,
                        },
                    },
                    Origin {
                        step: step.clone(),
                        role,
                    },
                );
            }
        }
        Some(Self {
            recipe: Recipe::new(
                "canonical stock feature",
                vec![
                    Step {
                        id: profile,
                        op: Operation::Sketch { plane, loops },
                    },
                    Step {
                        id: step.clone(),
                        op: sweep,
                    },
                ],
                step.clone(),
            ),
            origins,
        })
    }

    /// The immutable sketch and sweep definition. Its construction step is never
    /// the owner of an output face.
    pub fn recipe(&self) -> &Recipe {
        &self.recipe
    }

    /// Output origin corresponding to a generated feature origin.
    pub fn origin<'a>(&'a self, generated: &'a Origin) -> &'a Origin {
        self.origins.get(generated).unwrap_or(generated)
    }
}
