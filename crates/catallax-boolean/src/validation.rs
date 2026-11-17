//! Structural validity and opposite-sheet overlaps on matching exact carriers.

use catallax_numbers::{Decider, Outcome, Refusal, StrictSign};
use catallax_shapes::Surface;
use catallax_topology::{BodyId, FaceId, Model, TopologyFault};

/// Checks topology and rejects positive-area overlaps of oppositely oriented faces.
///
/// Carrier matching covers analytic supports and definition-identical spline
/// carriers. This is not a general surface self-intersection test. Shared edges,
/// isolated contacts and disjoint trims on a common carrier are permitted.
///
/// Finite boxes are screened before carrier and normal proofs. Index construction
/// and candidate traversal consume the caller's allowance; unknown bounds retain
/// the full proof path.
/// Bound and proof queries use a model fork, leaving the input and its handles unchanged.
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
    if faces.len() < 2 {
        return Ok(faults);
    }
    dec.charge_work(faces.len() as u64)?;
    let mut work = model.clone();
    let mut boxes = Vec::with_capacity(faces.len());
    for &face in &faces {
        dec.check_deadline()?;
        dec.charge_work(1)?;
        // Optional bounds never replace the full proof for unsupported or
        // unsettled trims. The broad phase continues to charge the caller.
        boxes.push(
            crate::bounds::face_box(&mut work, face, true, dec)
                .ok()
                .flatten()
                .map(|bounds| bounds.map(|axis| [axis.lo(), axis.hi()])),
        );
    }
    crate::broad_phase::pairs(&boxes, dec, |first, second, dec| {
        dec.charge_work(1)?;
        let first = faces[first];
        let second = faces[second];
        let a = work.face(first);
        let b = work.face(second);
        if a.surface == b.surface && a.sense == b.sense {
            return Ok(());
        }
        // Distinct parametric definitions need a separate equality solver;
        // do not mistake equal bounding boxes for equal finite carriers.
        match (
            work.geom().surface(a.surface),
            work.geom().surface(b.surface),
        ) {
            (Surface::BSpline(a), Surface::BSpline(b))
                if !a.carrier().same_definition(b.carrier()) =>
            {
                return Ok(());
            }
            (Surface::BSpline(_), Surface::BSpline(_)) => {}
            (Surface::BSpline(_), _) | (_, Surface::BSpline(_)) => return Ok(()),
            _ => {}
        }
        if !crate::geometry::coincident(
            work.geom().same_surface(a.surface, b.surface, dec),
            "boolean.body-check-ambiguous",
        )? || !opposite(&mut work, first, second, dec)?
        {
            return Ok(());
        }
        let graph = crate::graph::arrangement(&mut work, body, &[first, second], dec)?.value;
        let (_, owners) =
            crate::imprinting::arrangement(&mut work, &graph, &[first, second], dec)?.value;
        for (piece, owner) in owners {
            dec.charge_work(1)?;
            if owner != first {
                continue;
            }
            let overlaps = crate::points::find(&mut work, piece, dec, |model, point, _, dec| {
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
        Ok(())
    })?;
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
    model: &mut Model,
    first: FaceId,
    second: FaceId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let a = model.face(first);
    let b = model.face(second);
    let (a_surface, b_surface, same_sense) = (a.surface, b.surface, a.sense == b.sense);
    let agrees = if a_surface == b_surface {
        true
    } else if let (Surface::Plane(a), Surface::Plane(b)) = (
        model.geom().surface(a_surface),
        model.geom().surface(b_surface),
    ) {
        agreement(dec.sign(&a.normal.dot(&b.normal)))?
    } else {
        crate::points::find(model, first, dec, |model, point, _, dec| {
            agreement(model.geom().normals_agree(a_surface, b_surface, point, dec)).map(Some)
        })?
    };
    Ok(agrees != same_sense)
}
