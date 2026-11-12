//! Bounded, charged planar source repair. Exact predicates are never relaxed.

use crate::exact::is_zero;
use crate::gap::{auto_cap, decimal};
use crate::planar_arithmetic as a;
use crate::{
    Curve, Disposition, Face, GapKind, ImportRefusal, Placement, PlanarRepairDiagnostic,
    PlanarRepairMethod, Proposal, Solid, SourceIncidence, StepFile, Surface, Vec3,
};
use catallax_numbers::{Decider, Rational, Sign};
use std::collections::{BTreeMap, BTreeSet};

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
    let reason = reason.into();
    let measured = (|| -> Result<_, ImportRefusal> {
        let mut d = Decider::default();
        a::checkpoint(&mut d)?;
        let mut residual = None;
        if distance2.is_none()
            && let Surface::Plane(p) = &face.surface
        {
            let normal2 = a::dot(&p.axis, &p.axis, &mut d)?;
            if !normal2.is_zero() {
                for v in face.loops.iter().flat_map(|b| {
                    b.edges
                        .iter()
                        .flat_map(|e| [&e.start.point, &e.end.point])
                        .chain(b.vertex.iter().map(|v| &v.point))
                }) {
                    let value = a::dot(&p.axis, &a::sub(v, &p.origin, &mut d)?, &mut d)?;
                    let value2 = a::arithmetic(
                        &a::arithmetic(&value, &value, '*', &mut d)?,
                        &normal2,
                        '/',
                        &mut d,
                    )?;
                    residual = Some(if let Some(previous) = residual {
                        a::maximum(previous, value2, &mut d)?
                    } else {
                        value2
                    });
                }
            }
        }
        let value = distance2.or(residual.as_ref());
        let measurement = if distance2.is_none() && residual.is_some() {
            "vertex_plane_residual"
        } else {
            measurement
        };
        let range = value.map(|v| a::sqrt_enclosure(v, &mut d)).transpose()?;
        Ok(PlanarRepairDiagnostic {
            solid: solid.name.clone(),
            face: face.name.clone(),
            stage: stage.into(),
            measurement: measurement.into(),
            deviation_lower_mm: range.as_ref().map(|r| decimal(&r.0)),
            deviation_upper_mm: range.as_ref().map(|r| decimal(&r.1)),
            limit_mm: crate::AUTO_CAP_MM.into(),
            feature_size_squared_mm2: None,
            reason,
        })
    })();
    match measured {
        Ok(detail) => ImportRefusal::PlanarRepairNotVerified {
            detail: Box::new(detail),
        },
        Err(error) => error,
    }
}

fn budget() -> Result<(), ImportRefusal> {
    a::checkpoint(&mut Decider::default())
}

fn off_plane(face: &Face) -> Result<bool, ImportRefusal> {
    let mut d = Decider::default();
    a::checkpoint(&mut d)?;
    if !linear(face) {
        return Ok(false);
    }
    let Surface::Plane(p) = &face.surface else {
        unreachable!()
    };
    for v in face
        .loops
        .iter()
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start.point, &e.end.point])
    {
        if !a::dot(&p.axis, &a::sub(v, &p.origin, &mut d)?, &mut d)?.is_zero() {
            return Ok(true);
        }
    }
    Ok(false)
}

fn refit_lines(solid: &Solid, result: &mut Solid) -> Result<(), ImportRefusal> {
    let mut d = Decider::default();
    for (before, face) in solid.faces.iter().zip(&mut result.faces) {
        for edge in face.loops.iter_mut().flat_map(|b| &mut b.edges) {
            a::checkpoint(&mut d)?;
            let Curve::Line(origin, direction) = &mut edge.curve else {
                continue;
            };
            let delta = a::sub(&edge.end.point, &edge.start.point, &mut d)?;
            if is_zero(&delta) {
                return Err(refused(
                    solid,
                    before,
                    "collapsed_edge",
                    "unavailable",
                    None,
                    "the repair would collapse an edge",
                ));
            }
            if !is_zero(&a::cross(
                &a::sub(&edge.start.point, origin, &mut d)?,
                direction,
                &mut d,
            )?) || !is_zero(&a::cross(
                &a::sub(&edge.end.point, origin, &mut d)?,
                direction,
                &mut d,
            )?) {
                if !linear(before) {
                    continue;
                }
                *origin = edge.start.point.clone();
                *direction = if edge.curve_same_sense {
                    delta
                } else {
                    delta.map(|v| -v)
                };
            }
        }
    }
    let mut shared = BTreeMap::new();
    for (before, face) in solid.faces.iter().zip(&result.faces) {
        for edge in face.loops.iter().flat_map(|b| &b.edges) {
            a::checkpoint(&mut d)?;
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

fn refit_face(face: &mut Face, d: &mut Decider) -> Result<bool, ImportRefusal> {
    let Surface::Plane(old) = &face.surface else {
        return Ok(false);
    };
    let points: Vec<_> = face
        .loops
        .iter()
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start.point, &e.end.point])
        .collect();
    let Some(p0) = points.first() else {
        return Ok(false);
    };
    let mut p1 = None;
    for point in &points {
        if !a::same(point, p0, d)? {
            p1 = Some(*point);
            break;
        }
    }
    let Some(p1) = p1 else { return Ok(false) };
    let edge = a::sub(p1, p0, d)?;
    let mut normal = None;
    for point in &points {
        let n = a::cross(&edge, &a::sub(point, p0, d)?, d)?;
        if !is_zero(&n) {
            normal = Some(n);
            break;
        }
    }
    let Some(normal) = normal else {
        return Ok(false);
    };
    for point in &points {
        if !a::dot(&a::sub(point, p0, d)?, &normal, d)?.is_zero() {
            return Ok(false);
        }
    }
    let sign = a::dot(&normal, &old.axis, d)?;
    let axis = if sign.sign() == Sign::Negative {
        normal.map(|v| -v)
    } else if sign.is_zero() {
        return Ok(false);
    } else {
        normal
    };
    let ref_dir = if a::dot(&old.ref_dir, &axis, d)?.is_zero() {
        old.ref_dir.clone()
    } else {
        edge
    };
    face.surface = Surface::Plane(Placement {
        origin: (*p0).clone(),
        axis,
        ref_dir,
    });
    Ok(true)
}

fn refit(solid: &Solid) -> Result<Solid, ImportRefusal> {
    let mut d = Decider::default();
    let mut result = solid.clone();
    for (before, face) in solid.faces.iter().zip(&mut result.faces) {
        a::checkpoint(&mut d)?;
        if off_plane(before)? && !refit_face(face, &mut d)? {
            return Err(refused(
                solid,
                before,
                "plane_refit",
                "unavailable",
                None,
                "the vertices do not define one consistently oriented planar face",
            ));
        }
    }
    refit_lines(solid, &mut result)?;
    Ok(result)
}

fn project(point: &Vec3, planes: &[&Placement]) -> Result<Option<Vec3>, ImportRefusal> {
    let mut d = Decider::default();
    let planes: Vec<_> = planes.iter().map(|p| (*p).clone()).collect();
    let result = match crate::incidence_repair::project(point, &planes, &mut d) {
        Ok(p) => p,
        Err(_) if catallax_numbers::current_operation_work().exhausted => {
            return Err(ImportRefusal::ResourceLimit);
        }
        Err(_) => return Ok(None),
    };
    for p in &planes {
        if !a::dot(&p.axis, &a::sub(&result, &p.origin, &mut d)?, &mut d)?.is_zero() {
            return Ok(None);
        }
    }
    Ok(Some(result))
}

fn snap(solid: &Solid) -> Result<Solid, ImportRefusal> {
    let mut d = Decider::default();
    let mut uses: BTreeMap<Vec3, Vec<usize>> = BTreeMap::new();
    let mut fixed = BTreeSet::new();
    for (fi, face) in solid.faces.iter().enumerate() {
        for vertex in face.loops.iter().flat_map(|b| {
            b.edges
                .iter()
                .flat_map(|e| [&e.start, &e.end])
                .chain(b.vertex.iter())
        }) {
            a::checkpoint(&mut d)?;
            if linear(face) {
                let faces = uses.entry(vertex.point.clone()).or_default();
                if !faces.contains(&fi) {
                    faces.push(fi)
                }
            } else {
                fixed.insert(vertex.point.clone());
            }
        }
    }
    let mut moves = BTreeMap::new();
    let mut destinations: BTreeMap<_, _> = fixed.iter().map(|p| (p.clone(), p.clone())).collect();
    let cap = auto_cap();
    let cap2 = a::arithmetic(&cap, &cap, '*', &mut d)?;
    for (point, faces) in &uses {
        a::checkpoint(&mut d)?;
        let planes: Vec<_> = faces
            .iter()
            .map(|&i| {
                let Surface::Plane(p) = &solid.faces[i].surface else {
                    unreachable!()
                };
                p
            })
            .collect();
        let mut anchor = &solid.faces[faces[0]];
        let mut incident = true;
        for &index in faces {
            let Surface::Plane(p) = &solid.faces[index].surface else {
                unreachable!()
            };
            if !a::dot(&p.axis, &a::sub(point, &p.origin, &mut d)?, &mut d)?.is_zero() {
                anchor = &solid.faces[index];
                incident = false;
                break;
            }
        }
        let result = if incident {
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
            project(point, &planes)?.ok_or_else(|| {
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
        let movement2 = a::squared(point, &result, &mut d)?;
        if a::compare(&movement2, &cap2, &mut d)?.is_gt() {
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
            && !a::same(&other, point, &mut d)?
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
        for vertex in face.loops.iter_mut().flat_map(|b| {
            b.edges
                .iter_mut()
                .flat_map(|e| [&mut e.start, &mut e.end])
                .chain(b.vertex.iter_mut())
        }) {
            a::checkpoint(&mut d)?;
            if let Some(point) = moves.get(&vertex.point) {
                vertex.point = point.clone()
            }
        }
    }
    refit_lines(solid, &mut result)?;
    Ok(result)
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

pub(crate) fn prepare(file: &StepFile, sha: &str) -> Result<Prepared, ImportRefusal> {
    let mut prepared = Prepared {
        file: file.clone(),
        methods: BTreeMap::new(),
        plans: Vec::new(),
    };
    let mut affected = Vec::new();
    for solid in &file.solids {
        let mut bad = false;
        for face in &solid.faces {
            if off_plane(face)? {
                bad = true;
                break;
            }
        }
        affected.push(bad);
    }
    for (si, solid) in file.solids.iter().enumerate() {
        if !affected[si] {
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
        assert_eq!(project(&point, &[]).unwrap(), Some(point.clone()));
        assert_eq!(project(&point, &[&x]).unwrap(), Some(vector([1, 4, 5])));
        assert_eq!(project(&point, &[&x, &y]).unwrap(), Some(vector([1, 2, 5])));
        assert_eq!(
            project(&point, &[&x, &y, &z, &redundant]).unwrap(),
            Some(vector([1, 2, 3]))
        );
    }

    #[test]
    fn projection_uses_euclidean_distance_for_oblique_constraints() {
        let point = vector([3, 4, 5]);
        let diagonal = plane([1, 0, 0], [1, 1, 0]);
        assert_eq!(
            project(&point, &[&diagonal]).unwrap(),
            Some(vector([0, 1, 5]))
        );
        let parallel = plane([2, 0, 0], [2, 2, 0]);
        assert_eq!(project(&point, &[&diagonal, &parallel]).unwrap(), None);
    }
}
