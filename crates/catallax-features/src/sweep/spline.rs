//! Exact translated polygon sections with a whole-path embedding proof.

use super::*;
use crate::loft::topology::{LoopInfo, assemble};
use catallax_numbers::Real;
use catallax_shapes::{BSplineCurve, BSplineSurface, RationalAtlas};

fn unproven() -> Refusal {
    Refusal::new(
        "sweep.spline-admissibility-unproven",
        "the spline must have strictly monotone Bezier controls along the profile normal; self-crossing paths and returning or tangent sections need a separate envelope proof",
    )
}

pub(super) fn build(
    model: &mut Model,
    profile: &Profile,
    path: &PathSegment,
    step: &StepId,
    id: &ElementId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let PathSegment::Spline {
        degree,
        control_points,
        weights,
        knots,
        ..
    } = path
    else {
        unreachable!("spline dispatch");
    };
    if profile.is_empty() {
        return Ok(catallax_topology::empty(model));
    }
    if profile.ellipse_frame.is_some()
        || profile
            .loops
            .iter()
            .flat_map(|b| &b.segments)
            .any(|s| s.circle.is_some() || s.bezier.is_some())
    {
        return Err(unsupported(
            "fixed spline sweeps currently need polygon profiles",
        ));
    }
    let points: Vec<_> = control_points.iter().map(point).collect();
    let knots: Vec<_> = knots.iter().map(|k| k.value().clone()).collect();
    let weights: Vec<_> = weights
        .iter()
        .map(|w| Real::from(w.value().clone()))
        .collect();
    let curve = BSplineCurve::new(
        *degree,
        knots.clone(),
        points.iter().map(|p| p.to_real()).collect(),
        weights.clone(),
        dec,
    )?;
    let normal = crate::profile::rational_vector(
        &profile.pose,
        &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
    );
    let mut direction = None;
    for (span, _) in curve.spans() {
        let controls = curve.span_bezier_controls(span, dec)?;
        let heights = controls
            .iter()
            .map(|p| {
                let h = &(&p[0] * &Real::from(normal.x.clone())
                    + &p[1] * &Real::from(normal.y.clone())
                    + &p[2] * &Real::from(normal.z.clone()))
                    / &p[3];
                h.rational_ref().cloned().ok_or_else(unproven)
            })
            .collect::<Result<Vec<_>, _>>()?;
        for pair in heights.windows(2) {
            let sign = (&pair[1] - &pair[0]).sign();
            if sign == Sign::Zero || direction.is_some_and(|d| d != sign) {
                return Err(unproven());
            }
            direction = Some(sign);
        }
    }
    let positive = direction.ok_or_else(unproven)? == Sign::Positive;
    // Positive weights and ordered projected Bezier controls make n.gamma(t)
    // strictly monotone. Each material section occupies a distinct plane, so
    // the complete volume map is injective and every side has full rank.
    let start = points.first().ok_or_else(unproven)?;
    let displacement: Vec<_> = points.iter().map(|p| p.minus(start)).collect();
    let mut atlases = Vec::new();
    let mut infos = Vec::new();
    for boundary in &profile.loops {
        let mut vertices = Vec::new();
        let mut ids = Vec::new();
        for segment in &boundary.segments {
            for (index, id) in segment.ids.iter().enumerate() {
                let fraction = Rational::from_ratio(index as i64, segment.ids.len() as i64)
                    .expect("segment origins");
                let xy: [Rational; 2] = std::array::from_fn(|axis| {
                    &segment.start[axis]
                        + &(&fraction * &(&segment.end[axis] - &segment.start[axis]))
                });
                vertices.push(profile.pose.apply_point(
                    &RPoint3::new(xy[0].clone(), xy[1].clone(), Rational::zero()).to_real(),
                ));
                ids.push(id.clone());
            }
        }
        vertices.push(vertices[0].clone());
        let count = vertices.len();
        let mut u_knots = vec![Rational::zero()];
        u_knots.extend((0..count).map(|i| {
            Rational::from_ratio(i as i64, (count - 1) as i64).expect("polygon boundary")
        }));
        u_knots.push(Rational::one());
        let mut net = Vec::new();
        let mut tensor_weights = Vec::new();
        for vertex in vertices {
            for (delta, weight) in displacement.iter().zip(&weights) {
                net.push(vertex.plus(&delta.to_real()));
                tensor_weights.push(weight.clone());
            }
        }
        // Tensor partition of unity gives P(u) + gamma(v) - gamma(0) exactly.
        let source = BSplineSurface::new(
            [1, *degree],
            [u_knots, knots.clone()],
            [count, points.len()],
            net,
            tensor_weights,
            dec,
        )?;
        atlases.push(RationalAtlas::partition(source, [true, false], dec)?);
        infos.push(LoopInfo {
            ids,
            depth: boundary.depth,
            parent: boundary.parent,
            material_forward: true,
        });
    }
    let delta = displacement.last().expect("validated path controls");
    let mut end = profile.clone();
    end.pose = end
        .pose
        .then(&Pose::from_rational(Rotation::identity(), delta.clone()));
    let body = assemble(
        model,
        atlases,
        &infos,
        &[profile, &end],
        [positive, positive],
        step,
        dec,
    )?;
    model.map_face_origins(body, |origin| Origin {
        step: step.clone(),
        role: match &origin.role {
            Role::Side { segment, .. } => Role::Side {
                segment: segment.clone(),
                path: Some(id.clone()),
            },
            role => role.clone(),
        },
    });
    Ok(body)
}
