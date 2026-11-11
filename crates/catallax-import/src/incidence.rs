//! Bounded planar source repair. Exact predicates are never relaxed.

use std::collections::{BTreeMap, BTreeSet};

use catallax_numbers::{Decider, Rational};

use crate::exact::{cross, dot, is_zero, sub};
use crate::gap::{auto_cap, decimal, sqrt_enclosure};
use crate::{
    Curve, Disposition, Face, GapKind, ImportRefusal, Placement, PlanarRepairDiagnostic,
    PlanarRepairMethod, Proposal, Solid, SourceIncidence, StepFile, Surface, Vec3,
};

pub(crate) struct Plan {
    pub solid: usize,
    pub proposal: Proposal,
    pub record: SourceIncidence,
}

pub(crate) struct Prepared {
    pub file: StepFile,
    pub methods: BTreeMap<usize, PlanarRepairMethod>,
    pub plans: Vec<Plan>,
}

pub(crate) fn refused(
    solid: &Solid,
    face: &Face,
    stage: &str,
    measurement: &str,
    distance2: Option<&Rational>,
    reason: impl Into<String>,
) -> ImportRefusal {
    let residual = if distance2.is_none() && Decider::default().check_deadline().is_ok() {
        if let Surface::Plane(p) = &face.surface {
            let normal2 = dot(&p.axis, &p.axis);
            if normal2.is_zero() {
                None
            } else {
                face.loops
                    .iter()
                    .flat_map(|b| {
                        b.edges
                            .iter()
                            .flat_map(|e| [&e.start.point, &e.end.point])
                            .chain(b.vertex.iter().map(|v| &v.point))
                    })
                    .map(|v| {
                        let value = dot(&p.axis, &sub(v, &p.origin));
                        &(&value * &value) / &normal2
                    })
                    .max()
            }
        } else {
            None
        }
    } else {
        None
    };
    let value = distance2.or(residual.as_ref());
    let measurement = if distance2.is_none() && residual.is_some() {
        "vertex_plane_residual"
    } else {
        measurement
    };
    let range = value.map(sqrt_enclosure);
    ImportRefusal::PlanarRepairNotVerified {
        detail: Box::new(PlanarRepairDiagnostic {
            solid: solid.name.clone(),
            face: face.name.clone(),
            stage: stage.to_string(),
            measurement: measurement.to_string(),
            deviation_lower_mm: range.as_ref().map(|r| decimal(&r.0)),
            deviation_upper_mm: range.as_ref().map(|r| decimal(&r.1)),
            limit_mm: crate::AUTO_CAP_MM.to_string(),
            feature_size_squared_mm2: None,
            reason: reason.into(),
        }),
    }
}

fn budget() -> Result<(), ImportRefusal> {
    Decider::default()
        .check_deadline()
        .map_err(|_| ImportRefusal::ResourceLimit)
}

fn failure_text(error: &ImportRefusal) -> String {
    match error {
        ImportRefusal::PlanarRepairNotVerified { detail } => {
            let bound = detail
                .deviation_upper_mm
                .as_ref()
                .map(|value| format!("; {} upper bound {value} mm", detail.measurement))
                .unwrap_or_default();
            format!("{} ({}{bound})", detail.reason, detail.stage)
        }
        ImportRefusal::ResourceLimit => "the deterministic work allowance ended".into(),
        ImportRefusal::NotSolid { detail } => detail.clone(),
        _ => "the candidate construction was not verified".into(),
    }
}

fn failure_information(error: &ImportRefusal) -> u8 {
    match error {
        ImportRefusal::PlanarRepairNotVerified { detail } => match detail.measurement.as_str() {
            "whole_face_source_bound" => 3,
            "vertex_correspondence" => 2,
            "source_boundary_correspondence" => 1,
            _ => 0,
        },
        _ => 0,
    }
}

fn method_name(method: PlanarRepairMethod) -> &'static str {
    match method {
        PlanarRepairMethod::RefitPlanes => "plane refitting",
        PlanarRepairMethod::SnapSharedVertices => "shared-vertex projection",
    }
}

fn linear(face: &Face) -> bool {
    matches!(face.surface, Surface::Plane(_))
        && !face.loops.is_empty()
        && face
            .loops
            .iter()
            .flat_map(|b| &b.edges)
            .all(|e| matches!(e.curve, Curve::Line(..)))
}

fn refit_lines(solid: &Solid, result: &mut Solid) -> Result<(), ImportRefusal> {
    for (before, face) in solid.faces.iter().zip(&mut result.faces) {
        for edge in face.loops.iter_mut().flat_map(|b| &mut b.edges) {
            let Curve::Line(origin, direction) = &mut edge.curve else {
                continue;
            };
            if is_zero(&sub(&edge.end.point, &edge.start.point)) {
                return Err(refused(
                    solid,
                    before,
                    "collapsed_edge",
                    "unavailable",
                    None,
                    "the repair would collapse an edge",
                ));
            }
            if !is_zero(&cross(&sub(&edge.start.point, origin), direction))
                || !is_zero(&cross(&sub(&edge.end.point, origin), direction))
            {
                if !linear(before) {
                    continue;
                }
                *origin = edge.start.point.clone();
                let d = sub(&edge.end.point, &edge.start.point);
                *direction = if edge.curve_same_sense {
                    d
                } else {
                    d.map(|v| -v)
                };
            }
        }
    }
    let mut shared = BTreeMap::new();
    for (before, face) in solid.faces.iter().zip(&result.faces) {
        for edge in face.loops.iter().flat_map(|b| &b.edges) {
            let geometry = (&edge.curve, &edge.start.point, &edge.end.point);
            if let Some(previous) = shared.insert(edge.edge.as_str(), geometry)
                && previous != geometry
            {
                return Err(refused(
                    solid,
                    before,
                    "shared_edge",
                    "unavailable",
                    None,
                    "a shared edge would have inconsistent geometry on its faces",
                ));
            }
        }
    }
    Ok(())
}

fn refit(solid: &Solid) -> Result<Solid, ImportRefusal> {
    let mut result = solid.clone();
    for (before, face) in solid.faces.iter().zip(&mut result.faces) {
        budget()?;
        if crate::intersect::planar_incidence_error(before).is_some() {
            crate::repair::refit(face).map_err(|reason| {
                refused(solid, before, "plane_refit", "unavailable", None, reason)
            })?;
        }
    }
    refit_lines(solid, &mut result)?;
    Ok(result)
}

fn reduce(rows: &mut [Vec<Rational>], variables: usize) -> Option<usize> {
    let mut rank = 0;
    for column in 0..variables {
        let Some(pivot) = (rank..rows.len()).find(|&i| !rows[i][column].is_zero()) else {
            continue;
        };
        rows.swap(rank, pivot);
        let divisor = rows[rank][column].clone();
        for value in &mut rows[rank] {
            *value = &*value / &divisor;
        }
        let pivot = rows[rank].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i == rank || row[column].is_zero() {
                continue;
            }
            let factor = row[column].clone();
            for (value, basis) in row.iter_mut().zip(&pivot) {
                *value = &*value - &(&factor * basis);
            }
        }
        rank += 1;
    }
    (!rows
        .iter()
        .any(|row| row[..variables].iter().all(Rational::is_zero) && !row[variables].is_zero()))
    .then_some(rank)
}

/// Orthogonal projection onto a consistent affine intersection of any rank.
fn project(point: &Vec3, planes: &[&Placement]) -> Option<Vec3> {
    let mut rows: Vec<Vec<Rational>> = planes
        .iter()
        .map(|p| {
            let mut row = p.axis.to_vec();
            row.push(dot(&p.axis, &p.origin));
            row
        })
        .collect();
    let rank = reduce(&mut rows, 3)?;
    let normals: Vec<Vec3> = rows[..rank]
        .iter()
        .map(|r| std::array::from_fn(|i| r[i].clone()))
        .collect();
    let mut gram: Vec<Vec<Rational>> = (0..rank)
        .map(|i| {
            let mut row: Vec<_> = normals.iter().map(|n| dot(&normals[i], n)).collect();
            row.push(&dot(&normals[i], point) - &rows[i][3]);
            row
        })
        .collect();
    if reduce(&mut gram, rank)? != rank {
        return None;
    }
    let result = std::array::from_fn(|i| {
        normals
            .iter()
            .enumerate()
            .fold(point[i].clone(), |value, (j, n)| {
                value - &gram[j][rank] * &n[i]
            })
    });
    planes
        .iter()
        .all(|p| dot(&p.axis, &sub(&result, &p.origin)).is_zero())
        .then_some(result)
}

fn snap(solid: &Solid) -> Result<Solid, ImportRefusal> {
    let mut uses: BTreeMap<Vec3, Vec<usize>> = BTreeMap::new();
    let mut fixed = BTreeSet::new();
    for (fi, face) in solid.faces.iter().enumerate() {
        for bound in &face.loops {
            for vertex in bound
                .edges
                .iter()
                .flat_map(|e| [&e.start, &e.end])
                .chain(bound.vertex.iter())
            {
                if linear(face) {
                    let faces = uses.entry(vertex.point.clone()).or_default();
                    if !faces.contains(&fi) {
                        faces.push(fi);
                    }
                } else {
                    fixed.insert(vertex.point.clone());
                }
            }
        }
    }
    let mut moves = BTreeMap::new();
    let mut destinations: BTreeMap<_, _> = fixed.iter().map(|p| (p.clone(), p.clone())).collect();
    for (point, faces) in &uses {
        budget()?;
        let planes: Vec<_> = faces
            .iter()
            .map(|&i| {
                let Surface::Plane(p) = &solid.faces[i].surface else {
                    unreachable!()
                };
                p
            })
            .collect();
        let anchor = &solid.faces[faces
            .iter()
            .copied()
            .find(|&i| {
                let Surface::Plane(p) = &solid.faces[i].surface else {
                    unreachable!()
                };
                !dot(&p.axis, &sub(point, &p.origin)).is_zero()
            })
            .unwrap_or(faces[0])];
        let result = if planes
            .iter()
            .all(|p| dot(&p.axis, &sub(point, &p.origin)).is_zero())
        {
            point.clone()
        } else if fixed.contains(point) {
            return Err(refused(
                solid,
                anchor,
                "curved_dependency",
                "unavailable",
                None,
                "a shared vertex on a curved boundary cannot move in the planar path",
            ));
        } else {
            project(point, &planes).ok_or_else(|| {
                refused(
                    solid,
                    anchor,
                    "vertex_constraints",
                    "unavailable",
                    None,
                    "incident planes have no exact shared-vertex solution",
                )
            })?
        };
        let delta = sub(point, &result);
        let movement2 = dot(&delta, &delta);
        let cap = auto_cap();
        if movement2 > &cap * &cap {
            return Err(refused(
                solid,
                anchor,
                "vertex_movement",
                "vertex_correspondence",
                Some(&movement2),
                "the proposed shared-vertex movement exceeds the import cap",
            ));
        }
        if let Some(other) = destinations.insert(result.clone(), point.clone())
            && other != *point
        {
            return Err(refused(
                solid,
                anchor,
                "vertex_identity",
                "vertex_correspondence",
                Some(&movement2),
                "the repair would join distinct source vertices",
            ));
        }
        moves.insert(point.clone(), result);
    }
    let mut result = solid.clone();
    for face in &mut result.faces {
        for bound in &mut face.loops {
            for vertex in bound
                .edges
                .iter_mut()
                .flat_map(|e| [&mut e.start, &mut e.end])
                .chain(bound.vertex.iter_mut())
            {
                if let Some(point) = moves.get(&vertex.point) {
                    vertex.point = point.clone();
                }
            }
        }
    }
    refit_lines(solid, &mut result)?;
    Ok(result)
}

pub(crate) fn prepare(file: &StepFile, sha: &str) -> Result<Prepared, ImportRefusal> {
    let mut prepared = Prepared {
        file: file.clone(),
        methods: BTreeMap::new(),
        plans: Vec::new(),
    };
    if !file
        .solids
        .iter()
        .flat_map(|s| &s.faces)
        .any(|face| crate::intersect::planar_incidence_error(face).is_some())
    {
        return Ok(prepared);
    }
    for (si, solid) in file.solids.iter().enumerate() {
        if !solid
            .faces
            .iter()
            .any(|f| crate::intersect::planar_incidence_error(f).is_some())
        {
            continue;
        }
        let mut failed = Vec::new();
        let mut selected = None;
        for method in [
            PlanarRepairMethod::RefitPlanes,
            PlanarRepairMethod::SnapSharedVertices,
        ] {
            budget()?;
            let candidate = match method {
                PlanarRepairMethod::RefitPlanes => refit(solid),
                PlanarRepairMethod::SnapSharedVertices => snap(solid),
            };
            let result = candidate.and_then(|candidate| {
                let (_, records) = crate::deviation::certify_incidence(solid, &candidate, method)?;
                Ok((candidate, records))
            });
            match result {
                Ok(result) => {
                    selected = Some((method, result));
                    break;
                }
                Err(ImportRefusal::ResourceLimit) => return Err(ImportRefusal::ResourceLimit),
                Err(error) => failed.push((method, error)),
            }
        }
        let Some((method, (candidate, records))) = selected else {
            let index = failed
                .iter()
                .enumerate()
                .max_by_key(|(_, (_, e))| failure_information(e))
                .map(|(i, _)| i)
                .expect("both routes refused");
            let (_, mut error) = failed.remove(index);
            if let ImportRefusal::PlanarRepairNotVerified { detail } = &mut error {
                for (method, other) in &failed {
                    detail.reason.push_str(&format!(
                        "; {} was also refused because {}",
                        method_name(*method),
                        failure_text(other)
                    ));
                }
            }
            return Err(error);
        };
        prepared.methods.insert(si, method);
        for record in records {
            let id = format!("p{}", prepared.plans.len());
            let construction = match method {
                PlanarRepairMethod::RefitPlanes => {
                    "Planar incidence repaired without moving shared vertices."
                }
                PlanarRepairMethod::SnapSharedVertices => {
                    "Planar incidence repaired using shared-vertex projection."
                }
            };
            let mut reason = format!(
                "{construction} Certified whole-face source bound {} mm within {} mm.",
                record.source_to_result_upper_mm,
                crate::AUTO_CAP_MM
            );
            for (method, error) in &failed {
                reason.push_str(&format!(
                    " {} was refused because {}.",
                    method_name(*method),
                    failure_text(error)
                ));
            }
            prepared.plans.push(Plan {
                solid: si,
                proposal: Proposal {
                    question: format!("repair:{}:{id}", &sha[..12]),
                    id,
                    entities: vec![record.face.clone()],
                    kind: GapKind::PlanarIncidence,
                    distance_min_mm: record.geometric_deviation_lower_mm.clone(),
                    distance_max_mm: record.source_to_result_upper_mm.clone(),
                    tolerance_mm: crate::AUTO_CAP_MM.to_string(),
                    disposition: Disposition::Auto,
                    reason,
                },
                record,
            });
        }
        prepared.file.solids[si] = candidate;
    }
    Ok(prepared)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn vector(values: [i64; 3]) -> Vec3 {
        values.map(Rational::from_i64)
    }

    fn plane(origin: [i64; 3], axis: [i64; 3]) -> Placement {
        Placement {
            origin: vector(origin),
            axis: vector(axis),
            ref_dir: vector([0, 0, 0]),
        }
    }

    #[test]
    fn projection_handles_every_constraint_rank_and_redundant_planes() {
        let point = vector([3, 4, 5]);
        let x = plane([1, 0, 0], [1, 0, 0]);
        let y = plane([0, 2, 0], [0, 1, 0]);
        let z = plane([0, 0, 3], [0, 0, 1]);
        let redundant = plane([1, 2, 3], [2, 3, 4]);
        assert_eq!(project(&point, &[]), Some(point.clone()));
        assert_eq!(project(&point, &[&x]), Some(vector([1, 4, 5])));
        assert_eq!(project(&point, &[&x, &y]), Some(vector([1, 2, 5])));
        assert_eq!(
            project(&point, &[&x, &y, &z, &redundant]),
            Some(vector([1, 2, 3]))
        );
    }

    #[test]
    fn projection_uses_euclidean_distance_for_oblique_constraints() {
        let point = vector([3, 4, 5]);
        let diagonal = plane([1, 0, 0], [1, 1, 0]);
        assert_eq!(project(&point, &[&diagonal]), Some(vector([0, 1, 5])));
        let parallel = plane([2, 0, 0], [2, 2, 0]);
        assert_eq!(project(&point, &[&diagonal, &parallel]), None);
    }
}
