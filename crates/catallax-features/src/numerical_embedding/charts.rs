//! Lift periodic parameters along the retained boundary, not vertex enumeration.

use super::*;

#[allow(clippy::too_many_arguments)]
pub(super) fn lift(
    model: &Model,
    body: BodyId,
    carriers: &[Carrier<'_>],
    edge_ends: &BTreeMap<EdgeId, [usize; 2]>,
    points: &[[f64; 3]],
    uvs: &mut BTreeMap<(usize, usize), Uv>,
    tolerances: &Tolerances,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for (index, id) in model.faces_of_body(body).iter().enumerate() {
        let Carrier::Native(surface) = &carriers[index] else { continue };
        let periodic = surface.domain().map(|kind| kind == catallax_shapes::ParamKind::Periodic);
        if !periodic.iter().any(|axis| *axis) {
            continue;
        }
        let face = ParametricFace::Native(surface);
        let mut anchor = None;
        for lp in &model.face(*id).loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else { unreachable!("prechecked") };
            let Some(first) = uses.first() else { continue };
            let first = model.coedge(*first);
            let edge = first.edge;
            if model.edge(edge).start == model.edge(edge).end
                && matches!(model.geom().curve(model.edge(edge).curve), Curve::Circle(_))
            {
                // Complete native rings get their explicit seam cell separately.
                continue;
            }
            let first_vertex = edge_ends[&edge][usize::from(first.sense == Sense::Reversed)];
            let initial = match anchor {
                Some(seed) => face.invert(points[first_vertex], Some(seed), tolerances, dec)?,
                None => uvs[&(index, first_vertex)].approx,
            };
            uvs.insert((index, first_vertex), Uv {
                exact: [exact(initial[0])?, exact(initial[1])?], approx: initial,
            });
            let mut previous = initial;
            let mut bounds = initial.map(|value| [value, value]);
            for coedge in uses {
                let coedge = model.coedge(*coedge);
                let edge = model.edge(coedge.edge);
                let natural_reverse = edge.sense == Sense::Reversed;
                let reversed = natural_reverse ^ (coedge.sense == Sense::Reversed);
                let mut ends = [edge.start, edge.end].map(|v| v.map(|v| model.vertex(v).point));
                if natural_reverse { ends.reverse(); }
                let trace = catallax_intersect::trace(model.geom(), edge.curve, ends[0], ends[1], dec)?.value;
                let pieces = trace.pieces();
                for ordinal in 0..pieces.len() {
                    let piece = &pieces[if reversed { pieces.len() - 1 - ordinal } else { ordinal }];
                    for sample in 1..=4 {
                        dec.charge_work(1)?;
                        let share = Rational::from_ratio(sample, 4).unwrap();
                        let parameter = if reversed { Rational::one() - share } else { share };
                        let point = piece.evaluate(&Ball::from_rational(&parameter, SAMPLE_BITS), dec)?
                            .each_ref().map(|v| v.mid().to_f64_nearest());
                        previous = face.invert(point, Some(previous), tolerances, dec)?;
                    }
                }
                let ends = edge_ends[&coedge.edge];
                let end = ends[usize::from(coedge.sense == Sense::Same)];
                let at = face.invert(points[end], Some(previous), tolerances, dec)?;
                for axis in 0..2 {
                    bounds[axis][0] = bounds[axis][0].min(at[axis]);
                    bounds[axis][1] = bounds[axis][1].max(at[axis]);
                    if end == first_vertex && periodic[axis]
                        && (at[axis] - initial[axis]).abs() > std::f64::consts::PI
                    {
                        return Err(unsupported(&format!(
                            "finite boundary on {id:?} winds in periodic axis {axis} and needs an explicit seam"
                        )));
                    }
                }
                if end != first_vertex {
                    uvs.insert((index, end), Uv {
                        exact: [exact(at[0])?, exact(at[1])?], approx: at,
                    });
                }
                previous = at;
            }
            anchor.get_or_insert_with(|| bounds.map(|[lo, hi]| (lo + hi) / 2.0));
        }
    }
    Ok(())
}
