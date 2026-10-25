//! Coaxial annular first-hit caps using exact circular torus sections.

use super::curved::Branch;
use crate::{Profile, profile::sign};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{FaceRef, Operation, Origin, Recipe, Role, Segment, StepId};
use catallax_shapes::{Point3, RPoint3, RVector3, Surface, Torus};
use catallax_topology::{BodyId, Model};

/// A whole-annulus branch proof. The finite construction cap lies strictly
/// within material for entry, or strictly beyond material for exit.
#[derive(Clone, Debug)]
pub struct AnnularCap {
    distance: Rational,
    branch: Branch,
    surface: Surface,
    radius: Rational,
    direction: RVector3,
}
fn coverage(detail: &str) -> Refusal {
    Refusal::new("features.unsupported-curved-limit", detail)
}

impl AnnularCap {
    /// Certifies every radial coordinate, rather than accepting sampled rays.
    pub fn compile(
        recipe: &Recipe,
        target: &StepId,
        profile: &StepId,
        direction: Option<&catallax_recipe::Vector3>,
        limit: Option<&FaceRef>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        recipe
            .check()
            .map_err(|e| Refusal::new(e.code(), e.to_string()))?;
        let Some(Operation::Torus { major, minor, pose }) = recipe.step(target).map(|s| &s.op)
        else {
            return Err(coverage(
                "annular first-hit construction requires a stock torus",
            ));
        };
        let placement = crate::profile::placement(pose, dec)?;
        let z = RVector3::new(Rational::zero(), Rational::zero(), Rational::one());
        let axis = crate::profile::rational_vector(&placement, &z);
        let center = placement.apply_point(&Point3::origin());
        let surface = Surface::Torus(Torus::new(
            center.clone(),
            axis.to_real(),
            major.value().clone().into(),
            minor.value().clone().into(),
            placement.apply_vector(
                &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
            ),
            dec,
        )?);
        if let Some(reference) = limit {
            let wall = Origin {
                step: target.clone(),
                role: Role::Side {
                    segment: catallax_recipe::ElementId::new("wall").expect("stock wall id"),
                    path: None,
                },
            };
            if !reference.face.matches(&wall) || !reference.touching.is_empty() {
                return Err(Refusal::new(
                    "recipe.missing-reference",
                    "torus limit has no matching shared-edge face",
                ));
            }
        }
        let Some(Operation::Sketch { plane, loops }) = recipe.step(profile).map(|s| &s.op) else {
            return Err(Refusal::new(
                "recipe.input-failed",
                "torus cap input is not a sketch",
            ));
        };
        if loops.len() != 2 || loops.iter().any(|l| l.0.len() != 1) {
            return Err(coverage(
                "torus limits require two coaxial full circles; other sections need certified torus branch tracing",
            ));
        }
        fn circle(
            boundary: &catallax_recipe::Loop,
        ) -> Result<(&catallax_recipe::Point2, &Rational), Refusal> {
            let Segment::Circle { centre, radius, .. } = &boundary.0[0] else {
                return Err(coverage(
                    "non-circular torus cap walls need certified torus branch tracing",
                ));
            };
            Ok((centre, radius.value()))
        }
        let circles = [circle(&loops[0])?, circle(&loops[1])?];
        if circles[0].0 != circles[1].0 {
            return Err(coverage("torus cap circles are not concentric"));
        }
        let sketch = crate::profile::placement(plane, dec)?;
        let normal = crate::profile::rational_vector(&sketch, &z);
        let direction = direction
            .map(|v| {
                let [x, y, z] = v.0.each_ref().map(|n| n.value().clone());
                RVector3::new(x, y, z)
            })
            .unwrap_or_else(|| normal.clone());
        if direction.is_zero() {
            return Err(Refusal::new(
                "extrude.zero-direction",
                "extrusion direction is zero",
            ));
        }
        if direction.dot(&normal).is_zero() {
            return Err(Refusal::new(
                "extrude.zero-thickness",
                "ray lies in the sketch plane",
            ));
        }
        if !axis.cross(&normal).is_zero() || !axis.cross(&direction).is_zero() {
            return Err(coverage(
                "oblique torus limits need certified torus branch tracing",
            ));
        }
        let centre = &circles[0].0.0;
        let origin = sketch.apply_point(
            &RPoint3::new(
                centre[0].value().clone(),
                centre[1].value().clone(),
                Rational::zero(),
            )
            .to_real(),
        );
        let displacement = origin.minus(&center);
        if sign(&displacement.cross(&axis.to_real()).norm_squared(), dec)? != Sign::Zero {
            return Err(coverage(
                "offset torus cap axes need certified torus branch tracing",
            ));
        }
        let height = displacement.dot(&axis.to_real());
        let height = height
            .rational_ref()
            .ok_or_else(|| coverage("torus axial placement is not rational"))?;
        let forward = direction.dot(&axis).sign();
        let start = if forward == Sign::Positive {
            height.clone()
        } else {
            -height
        };
        let inner = circles[0].1.min(circles[1].1);
        let outer = circles[0].1.max(circles[1].1);
        let major = major.value();
        let minor = minor.value();
        if inner <= &(major - minor) || outer >= &(major + minor) || inner == outer {
            return Err(coverage(
                "torus profile reaches a silhouette or has rays without finite hits",
            ));
        }
        let deviation = (inner - major).abs().max((outer - major).abs());
        let (branch, distance) = if start < -minor {
            (Branch::Entry, -start)
        } else if &start * &start + &deviation * &deviation < minor * minor {
            (Branch::Exit, minor - &start + Rational::one())
        } else {
            return Err(coverage(
                "torus start state is not uniformly outside-before-entry or inside material",
            ));
        };
        let radius = major
            + minor
            + pose
                .translation
                .0
                .iter()
                .fold(Rational::zero(), |s, n| s + n.value().abs());
        let profile_bound = outer
            + &plane
                .translation
                .0
                .iter()
                .fold(Rational::zero(), |s, n| s + n.value().abs())
            + centre
                .iter()
                .fold(Rational::zero(), |s, n| s + n.value().abs());
        Ok(Self {
            distance,
            branch,
            surface,
            radius: radius.max(profile_bound),
            direction,
        })
    }

    /// Physical construction length, independent of authored direction scale.
    pub fn distance(&self) -> &Rational {
        &self.distance
    }
    /// The uniform first crossing over the complete annulus.
    pub fn branch(&self) -> Branch {
        self.branch
    }
    /// Exact limiting support, also consumed by the independent field compiler.
    pub fn surface(&self) -> &Surface {
        &self.surface
    }
    /// Certified enclosing sphere for the result about world zero.
    pub fn radius(&self) -> &Rational {
        &self.radius
    }

    /// Builds only exact coaxial circle sections; no general torus tracing.
    pub(crate) fn build(
        &self,
        model: &mut Model,
        source: BodyId,
        profile: &Profile,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        dec.check_deadline()?;
        let column = crate::extrude(model, profile, step, &self.distance, Some(&self.direction))?;
        let face = model.faces_of_body(source)[0];
        let face = model.add_named_face(
            model.face(face).clone(),
            vec![Origin {
                step: step.clone(),
                role: Role::End,
            }],
        );
        let shell = model.add_shell(catallax_topology::Shell { faces: vec![face] });
        let lump = model.add_lump(catallax_topology::Lump {
            shells: vec![shell],
        });
        let source = model.add_body(catallax_topology::Body {
            lumps: vec![lump],
            ..model.body(source).clone()
        });
        let kind = if self.branch == Branch::Entry {
            catallax_boolean::BoolOp::Subtract
        } else {
            catallax_boolean::BoolOp::Intersect
        };
        let result = catallax_boolean::boolean(model, column, source, kind, dec)
            .map_err(|error| error.refusal)?
            .body;
        if !model.check(result).is_empty() {
            return Err(Refusal::new(
                "features.invalid-topology",
                "toroidal cap failed topology validation",
            ));
        }
        Ok(result)
    }
}
