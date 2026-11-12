//! A full cylindrical band refitted to two unchanged, complete circular rims.
use catallax_numbers::{Decider, Rational, Sign};

use crate::exact::is_zero;
use crate::planar_arithmetic as a;
use crate::{Curve, CylinderRadiusRepair, FaceTolerance, ImportRefusal, StepFile, Surface};

pub(crate) struct Prepared {
    pub file: StepFile,
    pub record: CylinderRadiusRepair,
    pub tolerance: FaceTolerance,
}

pub(crate) fn prepare(file: &StepFile) -> Result<Option<Prepared>, ImportRefusal> {
    let [solid] = file.solids.as_slice() else {
        return Ok(None);
    };
    if solid.faces.len() != 3 {
        return Ok(None);
    }
    let cylinders = solid
        .faces
        .iter()
        .enumerate()
        .filter_map(|(i, face)| {
            if let Surface::Cylinder(frame, radius) = &face.surface {
                Some((i, frame, radius))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let [(index, frame, radius)] = cylinders.as_slice() else {
        return Ok(None);
    };
    let mut dec = Decider::default();
    let normal2 = a::dot(&frame.axis, &frame.axis, &mut dec)?;
    if normal2.is_zero() {
        return Ok(None);
    }
    let mut rims = Vec::new();
    let mut target_radius = None;
    let mut centres = Vec::new();
    for (i, face) in solid.faces.iter().enumerate() {
        a::checkpoint(&mut dec)?;
        if i == *index {
            continue;
        }
        let Surface::Plane(plane) = &face.surface else {
            return Ok(None);
        };
        let [bound] = face.loops.as_slice() else {
            return Ok(None);
        };
        let [edge] = bound.edges.as_slice() else {
            return Ok(None);
        };
        let Curve::Circle(circle, r) = &edge.curve else {
            return Ok(None);
        };
        if bound.vertex.is_some()
            || edge.start.name != edge.end.name
            || !a::same(&edge.start.point, &edge.end.point, &mut dec)?
            || r.sign() != Sign::Positive
            || is_zero(&plane.axis)
            || is_zero(&circle.axis)
            || !is_zero(&a::cross(&plane.axis, &frame.axis, &mut dec)?)
            || !is_zero(&a::cross(&circle.axis, &frame.axis, &mut dec)?)
            || !is_zero(&a::cross(
                &a::sub(&circle.origin, &frame.origin, &mut dec)?,
                &frame.axis,
                &mut dec,
            )?)
            || !a::dot(
                &a::sub(&circle.origin, &plane.origin, &mut dec)?,
                &plane.axis,
                &mut dec,
            )?
            .is_zero()
        {
            return Ok(None);
        }
        if let Some(previous) = &target_radius {
            if !a::compare(previous, r, &mut dec)?.is_eq() {
                return Ok(None);
            }
        } else {
            target_radius = Some(r.clone());
        }
        centres.push(circle.origin.clone());
        rims.push(edge.edge.clone());
    }
    let side = &solid.faces[*index];
    if side.loops.iter().any(|bound| bound.vertex.is_some()) {
        return Ok(None);
    }
    let mut physical = Vec::new();
    for bound in &side.loops {
        for edge in &bound.edges {
            let seam = bound
                .edges
                .iter()
                .any(|other| other.edge == edge.edge && other.forward != edge.forward);
            if seam {
                if !matches!(edge.curve, Curve::Line(..)) {
                    return Ok(None);
                }
            } else {
                physical.push(edge);
            }
        }
    }
    if physical.len() != 2
        || physical[0].edge == physical[1].edge
        || physical
            .iter()
            .any(|edge| !matches!(edge.curve, Curve::Circle(..)) || !rims.contains(&edge.edge))
    {
        return Ok(None);
    }
    let Some(target_radius) = target_radius else {
        return Ok(None);
    };
    if a::compare(radius, &target_radius, &mut dec)?.is_eq() {
        return Ok(None);
    }
    let height2 = a::squared(&centres[0], &centres[1], &mut dec)?;
    if height2.is_zero() || rims[0] == rims[1] {
        return Ok(None);
    }
    let mut feature2 = a::minimum(
        height2,
        a::arithmetic(radius, radius, '*', &mut dec)?,
        &mut dec,
    )?;
    feature2 = a::minimum(
        feature2,
        a::arithmetic(&target_radius, &target_radius, '*', &mut dec)?,
        &mut dec,
    )?;
    let delta = a::arithmetic(radius, &target_radius, '-', &mut dec)?.abs();
    let delta2 = a::arithmetic(&delta, &delta, '*', &mut dec)?;
    let tenth = a::arithmetic(&delta2, &Rational::from_i64(100), '*', &mut dec)?;
    if a::compare(&delta, &crate::gap::auto_cap(), &mut dec)?.is_gt()
        || a::compare(&tenth, &feature2, &mut dec)?.is_gt()
    {
        return Err(ImportRefusal::ValidationFailed {
            id: solid.faces[*index].name.clone(),
            detail: format!(
                "radius discrepancy {delta} mm exceeds the 0.01 mm import cap or a tenth of the nearest feature; feature squared {feature2} mm^2"
            ),
        });
    }
    let mut result = file.clone();
    let Surface::Cylinder(_, repaired_radius) = &mut result.solids[0].faces[*index].surface else {
        unreachable!()
    };
    *repaired_radius = target_radius.clone();
    // Complete rings and the candidate's exact topology establish full angular
    // coverage. Equal-angle, equal-height pairing attains |r_source-r_result|
    // in both Hausdorff directions over the entire band, including its rims.
    let (lower, upper) = delta.to_decimal_bounds(20);
    Ok(Some(Prepared {
        file: result,
        record: CylinderRadiusRepair {
            solid: solid.name.clone(),
            face: solid.faces[*index].name.clone(),
            source_definition: "cylinder_between_cap_planes/1".into(),
            source_radius_mm: radius.to_string(),
            repaired_radius_mm: target_radius.to_string(),
            unchanged_rim_edges: rims,
        },
        tolerance: FaceTolerance {
            solid: solid.name.clone(),
            face: solid.faces[*index].name.clone(),
            deviation_lower_mm: lower,
            deviation_upper_mm: upper,
            approved_tolerance_mm: crate::AUTO_CAP_MM.into(),
            feature_size_squared_mm2: feature2.to_string(),
        },
    }))
}
