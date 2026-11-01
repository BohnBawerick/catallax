//! Analytic one-line sweeps with the authored displacement and path provenance.

use catallax_numbers::{Rational, Refusal};
use catallax_recipe::{ElementId, PathSegment, Point3, Role, StepId};
use catallax_shapes::RVector3;
use catallax_topology::{BodyId, Model};

use crate::Profile;

#[derive(Clone, Debug)]
pub(crate) enum Path {
    Empty,
    Line { id: ElementId, vector: RVector3 },
}

impl Path {
    pub(crate) fn new(start: &Point3, segments: &[PathSegment]) -> Result<Self, Refusal> {
        match segments {
            [] => Ok(Self::Empty),
            [PathSegment::Line { id, to }] => {
                let vector = RVector3::new(
                    to.0[0].value() - start.0[0].value(),
                    to.0[1].value() - start.0[1].value(),
                    to.0[2].value() - start.0[2].value(),
                );
                if vector.is_zero() {
                    return Err(Refusal::new(
                        "path.zero-length",
                        "a line path segment must have distinct endpoints",
                    ));
                }
                Ok(Self::Line {
                    id: id.clone(),
                    vector,
                })
            }
            _ => Err(Refusal::new(
                "features.unsupported-path",
                "B-rep sweeps currently need an empty path or one straight line segment",
            )),
        }
    }

    /// Translation never exceeds the L1 length of the endpoint displacement.
    pub(crate) fn extent(&self) -> Rational {
        match self {
            Self::Empty => Rational::zero(),
            Self::Line { vector, .. } => vector.x.abs() + vector.y.abs() + vector.z.abs(),
        }
    }

    pub(crate) fn build(
        &self,
        model: &mut Model,
        profile: &Profile,
        step: &StepId,
    ) -> Result<BodyId, Refusal> {
        let Self::Line { id, vector } = self else {
            return Ok(catallax_topology::empty(model));
        };
        let normal = crate::profile::rational_vector(
            &profile.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
        if profile.is_empty() || vector.dot(&normal).is_zero() {
            // A planar union has no 3D interior, even when it has positive area.
            return Ok(catallax_topology::empty(model));
        }
        // On one line T'=0, so both fixed and rotation-minimising modes have Q=I.
        // The nonzero normal component proves injectivity over the full profile:
        // normal height determines t, then q = x - t*vector. No folds or envelope.
        crate::extrude::extrude_vector(
            model,
            profile,
            step,
            vector,
            [Some(Role::Start), Some(Role::End)],
            Some(id),
        )
    }
}
