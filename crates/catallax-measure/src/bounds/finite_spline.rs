//! Bounded spline extrema and checked planar-cap boundary ownership.

use super::{Accuracy, Budget, Extrema, spline};
use catallax_numbers::{Decider, Refusal};
use catallax_shapes::{Curve, SplineSupport, Surface};
use catallax_topology::{
    BezierExtrusionBand, BodyId, EdgeId, FaceId, Loop, Model, RationalPatchFace,
};
use std::collections::BTreeSet;

pub(super) fn edge(
    model: &Model,
    id: EdgeId,
    accuracy: &Accuracy,
    budget: &mut Budget,
    dec: &Decider,
) -> Result<Extrema, Refusal> {
    let Curve::BSpline(curve) = model.geom().curve(model.edge(id).curve) else {
        return Err(Refusal::new(
            "measure.unsupported-surface",
            "a checked Bezier rail lost its polynomial curve",
        ));
    };
    spline::curve_extrema(curve, curve.domain(), accuracy, budget, dec)
}

pub(super) fn band(
    model: &Model,
    face: FaceId,
    accuracy: &Accuracy,
    budget: &mut Budget,
    dec: &mut Decider,
) -> Result<Option<Extrema>, Refusal> {
    let Some(band) = BezierExtrusionBand::of_face(model, face, dec)? else {
        return Ok(None);
    };
    // Every coordinate is affine along a generator. Its minimum and maximum
    // on the complete rectangular patch are attained on one of these rails.
    let [a, b] = band.rails();
    let mut extrema = edge(model, a, accuracy, budget, dec)?;
    extrema.merge(&edge(model, b, accuracy, budget, dec)?);
    Ok(Some(extrema))
}

pub(super) fn rims(
    model: &Model,
    body: Option<BodyId>,
    face: FaceId,
    dec: &mut Decider,
) -> Result<BTreeSet<EdgeId>, Refusal> {
    let faces = if let Some(body) = body {
        model.faces_of_body(body)
    } else {
        let mut faces = BTreeSet::new();
        for lp in &model.face(face).loops {
            if let Loop::Coedges(coedges) = model.loop_(*lp) {
                for c in coedges {
                    for adjacent in model.coedges_of_edge(model.coedge(*c).edge) {
                        if let Some((_, owner)) = model.owner_of_coedge(adjacent) {
                            faces.insert(owner);
                        }
                    }
                }
            }
        }
        faces.into_iter().collect()
    };
    let mut result = BTreeSet::new();
    for face in faces {
        match model.geom().surface(model.face(face).surface) {
            Surface::BSpline(SplineSupport::BezierExtrusion(_)) => {
                let band = BezierExtrusionBand::of_face(model, face, dec)?.ok_or_else(|| {
                    Refusal::new(
                        "measure.unsupported-surface",
                        "a Bezier rim has no checked band",
                    )
                })?;
                result.extend(band.rails());
            }
            Surface::BSpline(SplineSupport::RationalPatch(_)) => {
                let checked = RationalPatchFace::of_face(model, face, dec)?.ok_or_else(|| {
                    Refusal::new(
                        "measure.rational-patch-domain",
                        "a spline cap row has no checked rational face",
                    )
                })?;
                for (edge, _, _, _) in checked.rows() {
                    if matches!(
                        model.geom().curve(model.edge(edge).curve),
                        Curve::BSpline(_)
                    ) {
                        result.insert(edge);
                    }
                }
            }
            _ => {}
        }
    }
    Ok(result)
}

pub(super) fn cap(model: &Model, face: FaceId, rims: &BTreeSet<EdgeId>) -> bool {
    let mut found = false;
    for lp in &model.face(face).loops {
        if let Loop::Coedges(coedges) = model.loop_(*lp) {
            for c in coedges {
                let edge = model.coedge(*c).edge;
                if matches!(
                    model.geom().curve(model.edge(edge).curve),
                    Curve::BSpline(_)
                ) {
                    if !rims.contains(&edge) {
                        return false;
                    }
                    found = true;
                }
            }
        }
    }
    found
}
