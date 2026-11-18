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
    endpoints: &mut BTreeMap<(usize, EdgeId, usize), Uv>,
    tolerances: &Tolerances,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for (index, id) in model.faces_of_body(body).iter().enumerate() {
        let Carrier::Native(surface) = &carriers[index] else {
            continue;
        };
        let periodic = surface
            .domain()
            .map(|kind| kind == catallax_shapes::ParamKind::Periodic);
        if !periodic.iter().any(|axis| *axis) {
            continue;
        }
        let face = ParametricFace::Native(surface);
        let mut anchor = None;
        for lp in &model.face(*id).loops {
            let Loop::Coedges(uses) = model.loop_(*lp) else {
                unreachable!("prechecked")
            };
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
            uvs.insert(
                (index, first_vertex),
                Uv {
                    exact: [exact(initial[0])?, exact(initial[1])?],
                    approx: initial,
                },
            );
            let mut previous = initial;
            let mut bounds = initial.map(|value| [value, value]);
            for coedge in uses {
                let coedge = model.coedge(*coedge);
                let edge = model.edge(coedge.edge);
                let natural_reverse = edge.sense == Sense::Reversed;
                let reversed = natural_reverse ^ (coedge.sense == Sense::Reversed);
                let mut ends = [edge.start, edge.end].map(|v| v.map(|v| model.vertex(v).point));
                if natural_reverse {
                    ends.reverse();
                }
                let trace =
                    catallax_intersect::trace(model.geom(), edge.curve, ends[0], ends[1], dec)?
                        .value;
                let pieces = trace.pieces();
                for ordinal in 0..pieces.len() {
                    let piece = &pieces[if reversed {
                        pieces.len() - 1 - ordinal
                    } else {
                        ordinal
                    }];
                    for sample in 1..=4 {
                        dec.charge_work(1)?;
                        let share = Rational::from_ratio(sample, 4).unwrap();
                        let parameter = if reversed {
                            Rational::one() - share
                        } else {
                            share
                        };
                        let point = piece
                            .evaluate(&Ball::from_rational(&parameter, SAMPLE_BITS), dec)?
                            .each_ref()
                            .map(|v| v.mid().to_f64_nearest());
                        previous = face.invert(point, Some(previous), tolerances, dec)?;
                    }
                }
                let ends = edge_ends[&coedge.edge];
                let end = ends[usize::from(coedge.sense == Sense::Same)];
                let at = face.invert(points[end], Some(previous), tolerances, dec)?;
                let mut exact_at = [exact(at[0])?, exact(at[1])?];
                for axis in 0..2 {
                    bounds[axis][0] = bounds[axis][0].min(at[axis]);
                    bounds[axis][1] = bounds[axis][1].max(at[axis]);
                    if end == first_vertex && periodic[axis] {
                        let turns = ((at[axis] - initial[axis]) / std::f64::consts::TAU).round();
                        if turns.abs() > 1.0 {
                            return Err(unsupported(
                                "a boundary winds more than once around its carrier",
                            ));
                        }
                        exact_at[axis] = &uvs[&(index, first_vertex)].exact[axis]
                            + &(Ball::pi(BITS).mid().mul_pow2(1) * exact(turns)?);
                    }
                }
                if end != first_vertex {
                    uvs.insert(
                        (index, end),
                        Uv {
                            exact: [exact(at[0])?, exact(at[1])?],
                            approx: at,
                        },
                    );
                } else {
                    endpoints.insert(
                        (index, coedge.edge, end),
                        Uv {
                            exact: exact_at,
                            approx: at,
                        },
                    );
                }
                previous = at;
            }
            anchor.get_or_insert_with(|| bounds.map(|[lo, hi]| (lo + hi) / 2.0));
        }
    }
    Ok(())
}

pub(super) fn parameter_cuts(
    face: ParametricFace<'_>,
    edge: &Native,
    range: [f64; 2],
    ends: [&Uv; 2],
    tolerances: &Tolerances,
    dec: &mut Decider,
) -> Result<Vec<f64>, Refusal> {
    let mut cuts = Vec::new();
    if let Native::Spline(curve) = edge {
        cuts.extend_from_slice(curve.internal_knots());
    }
    if let Some(spline) = face.spline()
        && (0..2).any(|axis| !spline.surface.internal_knots(axis).is_empty())
    {
        const PROBES: usize = DENSITIES[0];
        let parameters: [f64; PROBES + 1] =
            std::array::from_fn(|i| lerp(range, i as f64 / PROBES as f64));
        let mut uv = [[0.0; 2]; PROBES + 1];
        uv[0] = ends[0].approx;
        uv[PROBES] = ends[1].approx;
        for i in 1..PROBES {
            uv[i] = face.invert(edge.point(parameters[i]), Some(uv[i - 1]), tolerances, dec)?;
        }
        for axis in 0..2 {
            let mut previous = None;
            for knot in spline.surface.internal_knots(axis) {
                if previous == Some(*knot) {
                    continue;
                }
                previous = Some(*knot);
                for i in 0..PROBES {
                    dec.charge_work(1)?;
                    let a = uv[i][axis] - knot;
                    let b = uv[i + 1][axis] - knot;
                    if a == 0.0 {
                        cuts.push(parameters[i]);
                    }
                    if b == 0.0 {
                        cuts.push(parameters[i + 1]);
                    }
                    if a * b >= 0.0 {
                        continue;
                    }
                    let mut bracket = [parameters[i], parameters[i + 1]];
                    let mut seeds = [uv[i], uv[i + 1]];
                    let increasing = b > a;
                    for _ in 0..64 {
                        dec.charge_work(1)?;
                        let parameter = lerp(bracket, 0.5);
                        // The bracket has closed to adjacent binary64 values.
                        #[allow(clippy::float_cmp)]
                        if parameter == bracket[0] || parameter == bracket[1] {
                            break;
                        }
                        let seed =
                            std::array::from_fn(|axis| lerp([seeds[0][axis], seeds[1][axis]], 0.5));
                        let at = face.invert(edge.point(parameter), Some(seed), tolerances, dec)?;
                        let side = usize::from((at[axis] < *knot) != increasing);
                        bracket[side] = parameter;
                        seeds[side] = at;
                    }
                    cuts.push(lerp(bracket, 0.5));
                }
            }
        }
    }
    let noise = 16.0 * f64::EPSILON * range[0].abs().max(range[1].abs()).max(range[1] - range[0]);
    cuts.retain(|t| *t > range[0] + noise && *t < range[1] - noise);
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|a, b| (*a - *b).abs() <= noise);
    cuts.insert(0, range[0]);
    cuts.push(range[1]);
    Ok(cuts)
}
