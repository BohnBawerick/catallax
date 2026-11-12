//! Trim routing for patch faces that share a physical rail with a fillet.

use catallax_shapes::{SplineSupport, Surface};
use catallax_topology::{FaceId, Loop, Model};

/// Selects general certified winding instead of the complete-four-row test.
/// This recognises an adjacent face, and grants no membership or incidence.
pub(crate) fn trimmed_patch(model: &Model, face: FaceId) -> bool {
    model.face(face).loops.iter().any(|lp| {
        let Loop::Coedges(coedges) = model.loop_(*lp) else {
            return false;
        };
        coedges.iter().any(|coedge| {
            model
                .coedges_of_edge(model.coedge(*coedge).edge)
                .into_iter()
                .any(|other| {
                    model.owner_of_coedge(other).is_some_and(|(_, neighbor)| {
                        neighbor != face
                            && matches!(
                                model.geom().surface(model.face(neighbor).surface),
                                Surface::BSpline(SplineSupport::SplineFillet(_))
                            )
                    })
                })
        })
    })
}
