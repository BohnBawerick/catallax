//! Numerical continuation of certified rim-defined cone carriers.

use super::*;
use catallax_shapes::Cone;

pub(super) fn eligible(
    bytes: &[u8],
    error: &ImportRefusal,
    dec: &mut Decider,
) -> Result<bool, ImportError> {
    let entity = match error {
        ImportRefusal::Unsupported { entity, type_name }
            if type_name.starts_with("CONICAL_SURFACE") || type_name.starts_with("EDGE_CURVE") =>
        {
            Some(entity.as_str())
        }
        _ => return Ok(false),
    };
    dec.charge_work(bytes.len() as u64)?;
    let file = crate::read_step(bytes)?;
    in_file(&file, entity, dec)
}

pub(super) fn in_file(
    file: &crate::StepFile,
    entity: Option<&str>,
    dec: &mut Decider,
) -> Result<bool, ImportError> {
    for (_, face) in file.named_faces() {
        dec.charge_work(1)?;
        if !matches!(face.surface, crate::Surface::Cone(..)) {
            continue;
        }
        if entity.is_some_and(|key| {
            key != face.name
                && !face
                    .loops
                    .iter()
                    .flat_map(|b| &b.edges)
                    .any(|e| e.edge == key)
        }) {
            continue;
        }
        // A failed cap or neighbouring-feature certificate cannot be bypassed.
        if crate::cone::rim_cone(face).is_ok_and(|cone| cone.is_some()) {
            return Ok(true);
        }
    }
    Ok(false)
}

fn fitted(face: &crate::Face) -> Result<crate::cone::RimCone, ImportError> {
    crate::cone::rim_cone(face)?.ok_or_else(|| unsupported("expected a fitted cone"))
}

pub(super) fn surface(face: &crate::Face) -> Result<Surface, ImportError> {
    let cone = fitted(face)?;
    let crate::Surface::Cone(frame, ..) = &face.surface else {
        unreachable!()
    };
    Ok(Surface::Cone(Cone {
        apex: point(&cone.apex),
        axis: vector(&cone.axis),
        slope: cone.physical_slope(),
        ref_dir: vector(&frame.ref_dir),
    }))
}

pub(super) fn parameters(
    p: &Vec3,
    face: &crate::Face,
    dec: &mut Decider,
) -> Result<[Rational; 2], ImportError> {
    let cone = fitted(face)?;
    let crate::Surface::Cone(source, ..) = &face.surface else {
        unreachable!()
    };
    let frame = crate::Placement {
        origin: cone.apex,
        axis: cone.axis,
        ref_dir: source.ref_dir.clone(),
    };
    let [x, y, z] = super::analytic::coordinates(p, &frame, dec)?;
    Ok([super::analytic::angle(&x, &y)?, z.mid()])
}

pub(super) fn rim_delta(
    edge: &crate::OrientedEdge,
    face: &crate::Face,
    forward: bool,
    dec: &mut Decider,
) -> Result<Rational, ImportError> {
    let crate::Curve::Circle(circle, _) = &edge.curve else {
        return Err(unsupported("a periodic cone rim must be a circle"));
    };
    let cone = fitted(face)?;
    let aligned =
        crate::planar_arithmetic::dot(&circle.axis, &cone.axis, dec)?.sign() == Sign::Positive;
    Ok(if (edge.curve_same_sense == forward) == aligned {
        super::analytic::turn()
    } else {
        -super::analytic::turn()
    })
}

pub(super) fn same_generator(
    face: &crate::Face,
    start: &Vec3,
    end: &Vec3,
) -> Result<bool, ImportError> {
    let cone = fitted(face)?;
    let delta = crate::exact::sub(end, start);
    Ok(cone.contains_line(start, &delta))
}
