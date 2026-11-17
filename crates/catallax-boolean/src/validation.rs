//! Structural validity and opposite-sheet overlaps on matching exact carriers.

use catallax_numbers::{Decider, Outcome, Refusal, StrictSign};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Model, TopologyFault};

/// Checks topology and rejects positive-area overlaps of oppositely oriented faces.
///
/// Carrier matching covers analytic supports and definition-identical spline
/// carriers. This is not a general surface self-intersection test. Shared edges,
/// isolated contacts and disjoint trims on a common carrier are permitted.
/// Proof queries use a lazy model fork, leaving the input and its handles unchanged.
pub fn check_body(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<Vec<TopologyFault>, Refusal> {
    let mut faults = model.check(body);
    if !faults.is_empty() {
        return Ok(faults);
    }
    let faces = model.faces_of_body(body);
    let mut scratch = None;
    let mut bounds = crate::bounds::FaceBounds::new(true);
    for (index, &first) in faces.iter().enumerate() {
        for &second in &faces[index + 1..] {
            dec.charge_work(1)?;
            let a = model.face(first);
            let b = model.face(second);
            if a.surface == b.surface && a.sense == b.sense {
                continue;
            }
            // Distinct parametric definitions need a separate equality solver;
            // do not mistake equal bounding boxes for equal finite carriers.
            match (
                model.geom().surface(a.surface),
                model.geom().surface(b.surface),
            ) {
                (Surface::BSpline(a), Surface::BSpline(b))
                    if !a.carrier().same_definition(b.carrier()) =>
                {
                    continue;
                }
                (Surface::BSpline(_), Surface::BSpline(_)) => {}
                (Surface::BSpline(_), _) | (_, Surface::BSpline(_)) => continue,
                _ => {}
            }
            if !crate::geometry::coincident(
                model.geom().same_surface(a.surface, b.surface, dec),
                "boolean.body-check-ambiguous",
            )? || !opposite(model, first, second, &mut scratch, dec)?
            {
                continue;
            }
            let work = scratch.get_or_insert_with(|| model.clone());
            if bounds.separated(work, first, second, dec) {
                continue;
            }
            let graph = crate::graph::arrangement(work, body, &[first, second], dec)?.value;
            let (_, owners) =
                crate::imprinting::arrangement(work, &graph, &[first, second], dec)?.value;
            for (piece, owner) in owners {
                dec.charge_work(1)?;
                if owner != first {
                    continue;
                }
                let overlaps = crate::points::find(work, piece, dec, |model, point, _, dec| {
                    Ok(match crate::region::locate(model, second, point, dec)? {
                        crate::region::Location::Inside => Some(true),
                        crate::region::Location::Outside => Some(false),
                        _ => None,
                    })
                })?;
                if overlaps {
                    faults.push(TopologyFault {
                        code: "topology.coincident-opposite-faces",
                        detail: format!(
                            "faces {first:?} and {second:?} have opposite overlapping interiors"
                        ),
                    });
                    break;
                }
            }
        }
    }
    Ok(faults)
}

fn agreement(outcome: Outcome<StrictSign>) -> Result<bool, Refusal> {
    match outcome {
        Outcome::Proven(StrictSign::Positive) => Ok(true),
        Outcome::Proven(StrictSign::Negative) => Ok(false),
        Outcome::Refused(error) => Err(error),
        other => Err(crate::geometry::refusal(
            "boolean.body-check-ambiguous",
            "coincident face normals could not be ordered",
        )
        .keeping_undecided(other)),
    }
}

fn opposite(
    model: &Model,
    first: FaceId,
    second: FaceId,
    scratch: &mut Option<Model>,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let a = model.face(first);
    let b = model.face(second);
    let agrees = if a.surface == b.surface {
        true
    } else if let (Surface::Plane(a), Surface::Plane(b)) = (
        model.geom().surface(a.surface),
        model.geom().surface(b.surface),
    ) {
        agreement(dec.sign(&a.normal.dot(&b.normal)))?
    } else {
        let work = scratch.get_or_insert_with(|| model.clone());
        crate::points::find(work, first, dec, |model, point, _, dec| {
            agreement(model.geom().normals_agree(a.surface, b.surface, point, dec)).map(Some)
        })?
    };
    Ok(agrees != (a.sense == b.sense))
}
