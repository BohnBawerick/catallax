//! A bounded quadric target, retained independently of a B-rep evaluation.

use super::*;
use catallax_recipe::{Operation, Origin, Recipe, Role, StepId};
use catallax_shapes::{Cone, Cylinder, Pose, RVector3, Rotation, Sphere};

/// A stock quadric and its closed planar trims.
#[derive(Clone, Debug)]
pub struct Scene {
    pub(super) surface: Surface,
    pub(super) trims: Vec<Plane>,
    pub(super) radius: Rational,
    pub(super) owner: StepId,
    pub(super) operation: Operation,
}

pub(super) fn placement(pose: &catallax_recipe::Pose) -> Result<Pose, Refusal> {
    let [w, x, y, z] = pose.rotation.each_ref().map(|n| n.value().clone());
    let [tx, ty, tz] = pose.translation.0.each_ref().map(|n| n.value().clone());
    Ok(Pose::from_rational(
        Rotation::from_quaternion(w, x, y, z)?,
        RVector3::new(tx, ty, tz),
    ))
}

impl Scene {
    pub(super) fn resolve_wall(&self, reference: &FaceRef) -> Result<(), Refusal> {
        let wall = self.wall_origin();
        let cap = |role| Origin {
            step: self.owner.clone(),
            role,
        };
        let mut caps = vec![];
        match &self.operation {
            Operation::Cylinder { .. } => caps.extend([cap(Role::Start), cap(Role::End)]),
            Operation::Cone { top_radius, .. } => {
                caps.push(cap(Role::Start));
                if !top_radius.value().is_zero() {
                    caps.push(cap(Role::End));
                }
            }
            _ => (),
        }
        if !reference.face.matches(&wall) {
            return Err(if caps.iter().any(|o| reference.face.matches(o)) {
                coverage("selected planar stock caps need a separate first-hit reduction")
            } else {
                Refusal::new("recipe.missing-reference", "limit matches no target face")
            });
        }
        if reference
            .touching
            .iter()
            .any(|p| !caps.iter().any(|o| p.matches(o)))
        {
            return Err(Refusal::new(
                "recipe.missing-reference",
                "selected wall has no shared edge with a touching pattern",
            ));
        }
        if reference.near.is_some() && caps.iter().any(|o| reference.face.matches(o)) {
            return Err(coverage(
                "near resolution between curved and planar faces needs a trimmed-distance proof",
            ));
        }
        Ok(())
    }

    /// Compiles a sphere, cylinder or conical frustum, without building topology.
    /// Other target graphs remain explicit coverage gaps.
    pub fn compile(recipe: &Recipe, target: &StepId, dec: &mut Decider) -> Result<Self, Refusal> {
        recipe
            .check()
            .map_err(|e| Refusal::new(e.code(), e.to_string()))?;
        dec.check_deadline()?;
        let op = &recipe
            .step(target)
            .ok_or_else(|| Refusal::new("recipe.input-failed", "curved limit target is missing"))?
            .op;
        let (pose, local_radius, surface, height) = match op {
            Operation::Sphere { radius, pose } => (
                pose,
                radius.value().clone(),
                Surface::Sphere(Sphere::new(
                    Point3::origin(),
                    radius.value().clone().into(),
                    axis(2),
                    axis(0),
                    dec,
                )?),
                None,
            ),
            Operation::Cylinder {
                radius,
                height,
                pose,
            } => {
                positive(height.value(), "topology.non-positive-size", dec)?;
                (
                    pose,
                    radius.value() + height.value(),
                    Surface::Cylinder(Cylinder::new(
                        Point3::origin(),
                        axis(2),
                        radius.value().clone().into(),
                        axis(0),
                        dec,
                    )?),
                    Some(height.value().clone()),
                )
            }
            Operation::Cone {
                base_radius,
                top_radius,
                height,
                pose,
            } => {
                positive(base_radius.value(), "topology.non-positive-radius", dec)?;
                positive(height.value(), "topology.non-positive-size", dec)?;
                if top_radius.value().sign() == Sign::Negative {
                    return Err(Refusal::new(
                        "topology.negative-radius",
                        "cone top radius is negative",
                    ));
                }
                let slope = (top_radius.value() - base_radius.value()) / height.value().clone();
                let surface = if slope.is_zero() {
                    Surface::Cylinder(Cylinder::new(
                        Point3::origin(),
                        axis(2),
                        base_radius.value().clone().into(),
                        axis(0),
                        dec,
                    )?)
                } else {
                    let apex = Point3::new(
                        Real::zero(),
                        Real::zero(),
                        Real::from(-base_radius.value() / slope.clone()),
                    );
                    Surface::Cone(Cone::new(
                        apex,
                        axis(2).scaled(&Real::from_i64(if slope.sign() == Sign::Positive {
                            1
                        } else {
                            -1
                        })),
                        Real::from(slope.abs()),
                        axis(0),
                        dec,
                    )?)
                };
                (
                    pose,
                    base_radius.value().clone().max(top_radius.value().clone())
                        + height.value().clone(),
                    surface,
                    Some(height.value().clone()),
                )
            }
            _ => {
                return Err(coverage(
                    "curved limits currently require one stock sphere, cylinder or cone",
                ));
            }
        };
        let map = placement(pose)?;
        let trims = height
            .into_iter()
            .flat_map(|h| {
                [
                    Plane {
                        origin: Point3::origin(),
                        normal: axis(2).scaled(&Real::from_i64(-1)),
                        ref_dir: axis(0),
                    },
                    Plane {
                        origin: Point3::new(Real::zero(), Real::zero(), h.into()),
                        normal: axis(2),
                        ref_dir: axis(0),
                    },
                ]
            })
            .map(|p| match Surface::Plane(p).transformed(&map) {
                Surface::Plane(p) => p,
                _ => unreachable!("rigid plane image"),
            })
            .collect();
        let radius = local_radius
            + pose
                .translation
                .0
                .iter()
                .fold(Rational::zero(), |s, n| s + n.value().abs());
        Ok(Self {
            surface: surface.transformed(&map),
            trims,
            radius,
            owner: target.clone(),
            operation: op.clone(),
        })
    }

    /// The unchanged analytic limiting support.
    pub fn surface(&self) -> &Surface {
        &self.surface
    }

    /// The exact closed trim halfspaces; their material side is nonpositive.
    pub fn trims(&self) -> &[Plane] {
        &self.trims
    }

    /// A certified sphere about world zero containing the entire target.
    pub fn radius(&self) -> &Rational {
        &self.radius
    }

    /// The stock wall's persistent origin.
    pub fn wall_origin(&self) -> Origin {
        Origin {
            step: self.owner.clone(),
            role: Role::Side {
                segment: catallax_recipe::ElementId::new("wall").expect("fixed wall id"),
                path: None,
            },
        }
    }
}

fn positive(value: &Rational, code: &'static str, dec: &mut Decider) -> Result<(), Refusal> {
    dec.check_deadline()?;
    if value.sign() != Sign::Positive {
        return Err(Refusal::new(code, "stock dimension must be positive"));
    }
    Ok(())
}
