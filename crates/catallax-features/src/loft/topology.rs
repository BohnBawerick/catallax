//! Sew exact chart cells, then publish only maximal logical faces and edges.

use super::*;
use catallax_numbers::MarkLog;
use catallax_recipe::{ElementId, Origin, Role, StepId};
use catallax_shapes::{Curve, Plane, PointDef, RationalAtlas, SplineBoundary, Surface, Vector3};
use catallax_topology::{
    AtlasFaceGroup, Body, BodyId, BodyKind, Edge, EdgeId, Lump, Model, Sense, Shell, Vertex,
    VertexId,
};
use std::collections::BTreeMap;

pub(crate) struct LoopInfo {
    pub ids: Vec<ElementId>,
    pub depth: usize,
    pub parent: Option<usize>,
    pub material_forward: bool,
}

type PointKey = (usize, Rational, Rational);
type RowKey = (usize, usize, Rational, Rational, Rational);
type CapRow = (Rational, EdgeId, Sense);

fn normalized_u(value: &Rational) -> Rational {
    if *value == Rational::one() {
        Rational::zero()
    } else {
        value.clone()
    }
}
fn opposite(sense: Sense) -> Sense {
    crate::extrude::reverse(sense)
}

fn vertex(
    model: &mut Model,
    vertices: &mut BTreeMap<PointKey, VertexId>,
    key: PointKey,
    point: &catallax_shapes::Point3,
) -> VertexId {
    *vertices.entry(key).or_insert_with(|| {
        let point = model
            .geom_mut()
            .add_point(PointDef::Explicit(point.clone()));
        model.add_vertex(Vertex { point })
    })
}

fn plane(profile: &Profile) -> Plane {
    Plane {
        origin: profile
            .pose
            .apply_point(&catallax_shapes::RPoint3::origin().to_real()),
        normal: profile
            .pose
            .rotation
            .apply(&Vector3::new(Real::zero(), Real::zero(), Real::one())),
        ref_dir: profile.pose.rotation.apply(&Vector3::new(
            Real::one(),
            Real::zero(),
            Real::zero(),
        )),
    }
}

fn cap_path(mut rows: Vec<CapRow>) -> Result<Vec<(EdgeId, Sense)>, Refusal> {
    rows.sort_by(|a, b| a.0.cmp(&b.0));
    let direction = rows
        .first()
        .ok_or_else(|| unproven("the loft cap has no complete boundary cover"))?
        .2;
    if rows.iter().any(|r| r.2 != direction) {
        return Err(unproven("cap boundary orientation changes within a loop"));
    }
    if direction == Sense::Reversed {
        rows.reverse();
    }
    Ok(rows
        .into_iter()
        .map(|(_, edge, sense)| (edge, sense))
        .collect())
}

pub(crate) fn assemble(
    model: &mut Model,
    atlases: Vec<Vec<RationalAtlas>>,
    infos: &[LoopInfo],
    profiles: &[&Profile],
    orientation: [bool; 2],
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let [positive, end_positive] = orientation;
    let mut vertices = BTreeMap::new();
    let mut rows: BTreeMap<RowKey, EdgeId> = BTreeMap::new();
    let mut groups = Vec::new();
    let mut side_faces = vec![Vec::new(); infos.len()];
    let mut cap_rows: Vec<[Vec<CapRow>; 2]> =
        (0..infos.len()).map(|_| [Vec::new(), Vec::new()]).collect();
    let end_station = Rational::from_i64((profiles.len() - 1) as i64);
    for (loop_index, (components, info)) in atlases.into_iter().zip(infos).enumerate() {
        let sense = if info.material_forward == positive {
            Sense::Same
        } else {
            Sense::Reversed
        };
        for atlas in components {
            let mut faces = Vec::new();
            for chart in atlas.charts() {
                let patch = chart.patch();
                let bounds = patch.domain();
                let [nu, nv] = patch.carrier().control_count();
                let controls = patch.carrier().controls();
                let keys = [
                    (
                        loop_index,
                        normalized_u(&bounds[0][0]),
                        bounds[1][0].clone(),
                    ),
                    (
                        loop_index,
                        normalized_u(&bounds[0][1]),
                        bounds[1][0].clone(),
                    ),
                    (
                        loop_index,
                        normalized_u(&bounds[0][1]),
                        bounds[1][1].clone(),
                    ),
                    (
                        loop_index,
                        normalized_u(&bounds[0][0]),
                        bounds[1][1].clone(),
                    ),
                ];
                let corner_indices = [0, (nu - 1) * nv, nu * nv - 1, nv - 1];
                let corners: [VertexId; 4] = std::array::from_fn(|i| {
                    vertex(
                        model,
                        &mut vertices,
                        keys[i].clone(),
                        &controls[corner_indices[i]],
                    )
                });
                let descriptions = [
                    (
                        SplineBoundary::VMin,
                        (
                            loop_index,
                            0,
                            bounds[1][0].clone(),
                            bounds[0][0].clone(),
                            bounds[0][1].clone(),
                        ),
                        [0, 1],
                        Sense::Same,
                    ),
                    (
                        SplineBoundary::UMax,
                        (
                            loop_index,
                            1,
                            normalized_u(&bounds[0][1]),
                            bounds[1][0].clone(),
                            bounds[1][1].clone(),
                        ),
                        [1, 2],
                        Sense::Same,
                    ),
                    (
                        SplineBoundary::VMax,
                        (
                            loop_index,
                            0,
                            bounds[1][1].clone(),
                            bounds[0][0].clone(),
                            bounds[0][1].clone(),
                        ),
                        [3, 2],
                        Sense::Reversed,
                    ),
                    (
                        SplineBoundary::UMin,
                        (
                            loop_index,
                            1,
                            normalized_u(&bounds[0][0]),
                            bounds[1][0].clone(),
                            bounds[1][1].clone(),
                        ),
                        [0, 3],
                        Sense::Reversed,
                    ),
                ];
                let mut path = Vec::new();
                for (side, key, ends, direction) in descriptions {
                    let edge = *rows.entry(key).or_insert_with(|| {
                        let curve = model
                            .geom_mut()
                            .add_curve(Curve::BSpline(patch.boundary(side)));
                        model.add_edge(Edge {
                            curve,
                            sense: Sense::Same,
                            start: Some(corners[ends[0]]),
                            end: Some(corners[ends[1]]),
                        })
                    });
                    let used = if sense == Sense::Same {
                        direction
                    } else {
                        opposite(direction)
                    };
                    if side == SplineBoundary::VMin && bounds[1][0] == Rational::zero() {
                        cap_rows[loop_index][0].push((bounds[0][0].clone(), edge, opposite(used)));
                    }
                    if side == SplineBoundary::VMax && bounds[1][1] == end_station {
                        cap_rows[loop_index][1].push((bounds[0][0].clone(), edge, opposite(used)));
                    }
                    path.push((edge, used));
                }
                if sense == Sense::Reversed {
                    path.reverse();
                }
                let middle = (&bounds[0][0] + &bounds[0][1]).mul_pow2(-1);
                let segment = info
                    .ids
                    .iter()
                    .enumerate()
                    .find(|(index, _)| {
                        middle
                            < Rational::from_ratio((index + 1) as i64, info.ids.len() as i64)
                                .expect("authored segment count")
                    })
                    .map(|(_, id)| id)
                    .expect("common chart stays within an authored segment");
                let face = crate::extrude::face(
                    model,
                    Surface::BSpline(patch.clone().into()),
                    sense,
                    vec![path],
                    None,
                    vec![Origin {
                        step: step.clone(),
                        role: Role::Side {
                            segment: segment.clone(),
                            path: None,
                        },
                    }],
                );
                faces.push(face);
                side_faces[loop_index].push(face);
            }
            groups.push(AtlasFaceGroup { atlas, faces });
        }
    }
    let cap_planes = [
        plane(profiles[0]),
        plane(profiles.last().expect("loft sections")),
    ];
    let mut lumps = Vec::new();
    for (outer, info) in infos
        .iter()
        .enumerate()
        .filter(|(_, info)| info.depth % 2 == 0)
    {
        let indices: Vec<_> = std::iter::once(outer)
            .chain(infos.iter().enumerate().filter_map(|(index, candidate)| {
                (candidate.parent == Some(outer) && candidate.depth == info.depth + 1)
                    .then_some(index)
            }))
            .collect();
        let mut faces = Vec::new();
        for cap in 0..2 {
            let paths = indices
                .iter()
                .map(|index| cap_path(std::mem::take(&mut cap_rows[*index][cap])))
                .collect::<Result<Vec<_>, Refusal>>()?;
            let sense = if if cap == 0 { !positive } else { end_positive } {
                Sense::Same
            } else {
                Sense::Reversed
            };
            faces.push(crate::extrude::face(
                model,
                Surface::Plane(cap_planes[cap].clone()),
                sense,
                paths,
                None,
                vec![Origin {
                    step: step.clone(),
                    role: if cap == 0 { Role::Start } else { Role::End },
                }],
            ));
        }
        for index in indices {
            faces.append(&mut side_faces[index]);
        }
        let shell = model.add_shell(Shell { faces });
        lumps.push(model.add_lump(Lump {
            shells: vec![shell],
        }));
    }
    let cells = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
        marks: MarkLog::new(),
    });
    catallax_topology::coalesce_rational_charts(model, cells, groups, dec)
}
