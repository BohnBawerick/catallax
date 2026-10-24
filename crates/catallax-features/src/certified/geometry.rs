//! Complete three-face prism validation; construction facts and marks are not proof.

use super::{sign, unproven};
use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, MarkLog, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{
    BSplineCurve, Curve, EllipticCylinder, Plane, Point3, SplineBoundary, SplineSurface, Surface,
    Vector3,
};
use catallax_topology::{
    Body, BodyId, BodyKind, Edge, EdgeId, FaceId, Loop, Lump, Model, Sense, Shell,
};

#[derive(Clone, Debug)]
pub(super) struct PrismTopology {
    /// Lower cap, upper cap, outward side.
    pub faces: [FaceId; 3],
    pub rings: [EdgeId; 2],
}

#[derive(Clone, Copy)]
pub(super) enum ProfileSupport<'a> {
    Analytic,
    Spline(&'a SplineSurface),
}

fn sense(sense: Sense) -> i64 {
    match sense {
        Sense::Same => 1,
        Sense::Reversed => -1,
    }
}

fn equal(a: &Real, b: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    if sign(&(a - b), dec)? != Sign::Zero {
        return Err(unproven(
            "prism geometry does not match its retained definition",
        ));
    }
    Ok(())
}

fn same_point(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<(), Refusal> {
    equal(&a.x, &b.x, dec)?;
    equal(&a.y, &b.y, dec)?;
    equal(&a.z, &b.z, dec)
}

fn parallel(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<i64, Refusal> {
    if sign(&a.norm_squared(), dec)? != Sign::Positive
        || sign(&b.norm_squared(), dec)? != Sign::Positive
        || sign(&a.cross(b).norm_squared(), dec)? != Sign::Zero
    {
        return Err(unproven(
            "prism boundary direction is not the claimed nonzero parallel",
        ));
    }
    match sign(&a.dot(b), dec)? {
        Sign::Positive => Ok(1),
        Sign::Negative => Ok(-1),
        Sign::Zero => Err(unproven("parallel prism directions have zero pairing")),
    }
}

fn ring_use(model: &Model, face: FaceId, loop_index: usize) -> Result<(EdgeId, Sense), Refusal> {
    let lp = model
        .face(face)
        .loops
        .get(loop_index)
        .ok_or_else(|| unproven("missing prism trim"))?;
    let Loop::Coedges(coedges) = model.loop_(*lp) else {
        return Err(unproven("prism trim is not a whole closed ring"));
    };
    let [coedge] = coedges.as_slice() else {
        return Err(unproven(
            "prism trim contains a split or extra physical edge",
        ));
    };
    let coedge = model.coedge(*coedge);
    let edge = model.edge(coedge.edge);
    if edge.start.is_some() || edge.end.is_some() {
        return Err(unproven("a prism ring must not invent seam vertices"));
    }
    Ok((coedge.edge, coedge.sense))
}

fn analytic_ring(
    model: &Model,
    edge: EdgeId,
    original: &EllipticCylinder,
    axis: &Vector3,
    level: &Real,
    dec: &mut Decider,
) -> Result<i64, Refusal> {
    let CurveForm::Ellipse {
        center,
        normal,
        major_dir,
        major_sq,
        minor_sq,
    } = curve_form(model.geom(), model.edge(edge).curve, dec)?
    else {
        return Err(unproven(
            "source prism trim is not the complete analytic ellipse",
        ));
    };
    same_point(&center, &original.origin.plus(&axis.scaled(level)), dec)?;
    equal(&major_sq, &original.major_sq, dec)?;
    equal(&minor_sq, &original.minor_sq, dec)?;
    parallel(&major_dir, &original.major_dir, dec)?;
    if sign(&normal.dot(&major_dir), dec)? != Sign::Zero {
        return Err(unproven("source ellipse frame is not planar"));
    }
    parallel(&normal, axis, dec)
}

/// Checks actual topology, complete trims, material sense and naming, not a
/// producer's entity list or a zero-valued witness on each supporting plane.
pub(super) fn check_prism(
    model: &Model,
    body: BodyId,
    original: &EllipticCylinder,
    axial: [&Real; 2],
    profile: ProfileSupport<'_>,
    origins: &[Vec<Origin>; 3],
    dec: &mut Decider,
) -> Result<PrismTopology, Refusal> {
    dec.check_deadline()?;
    if let ProfileSupport::Spline(expected) = profile {
        let checked = catallax_topology::SplinePrism::check(model, body, dec)?;
        if !checked
            .spline()
            .carrier()
            .same_definition(expected.carrier())
        {
            return Err(unproven(
                "represented side is detached from its checked control net",
            ));
        }
        for (face, expected) in checked.faces().iter().zip(origins) {
            if expected.is_empty() || model.face_origins(*face) != expected.as_slice() {
                return Err(unproven(
                    "spline face origins do not match the immutable recipe",
                ));
            }
        }
        return Ok(PrismTopology {
            faces: checked.faces(),
            rings: checked.rings(),
        });
    }
    let entity = model.body(body);
    let [lump] = entity.lumps.as_slice() else {
        return Err(unproven(
            "prism certificate requires exactly one material lump",
        ));
    };
    let [shell] = model.lump(*lump).shells.as_slice() else {
        return Err(unproven(
            "prism certificate cannot omit a void or extra shell",
        ));
    };
    let faces = &model.shell(*shell).faces;
    if entity.kind != BodyKind::Solid
        || faces.len() != 3
        || faces[0] == faces[1]
        || faces[0] == faces[2]
        || faces[1] == faces[2]
        || !model.check(body).is_empty()
    {
        return Err(unproven(
            "actual prism topology is not a valid three-face solid",
        ));
    }
    let axis = original
        .axis
        .scaled(&(Real::one() / original.axis.norm_squared().sqrt()));
    let mut caps = [None, None];
    let mut side = None;
    for &face in faces {
        match model.geom().surface(model.face(face).surface) {
            Surface::Plane(plane) => {
                Plane::new(
                    plane.origin.clone(),
                    plane.normal.clone(),
                    plane.ref_dir.clone(),
                    dec,
                )?;
                let normal_sense = parallel(&plane.normal, &axis, dec)?;
                let level = plane.origin.minus(&original.origin).dot(&axis);
                let index = if sign(&(&level - axial[0]), dec)? == Sign::Zero {
                    0
                } else if sign(&(&level - axial[1]), dec)? == Sign::Zero {
                    1
                } else {
                    return Err(unproven("cap plane is outside its exact target interval"));
                };
                if caps[index].replace(face).is_some()
                    || model.face(face).loops.len() != 1
                    || normal_sense * sense(model.face(face).sense) != [-1, 1][index]
                {
                    return Err(unproven(
                        "cap inventory, trims or outward material sense disagree",
                    ));
                }
            }
            surface => {
                if side.replace(face).is_some()
                    || model.face(face).sense != Sense::Same
                    || model.face(face).loops.len() != 2
                {
                    return Err(unproven("prism side inventory or material sense disagrees"));
                }
                let Surface::EllipticCylinder(actual) = surface else {
                    return Err(unproven(
                        "source prism side has the wrong genuine support kind",
                    ));
                };
                if !Surface::EllipticCylinder(actual.clone())
                    .same_definition(&Surface::EllipticCylinder(original.clone()))
                {
                    return Err(unproven("source side lost its original exact support"));
                }
            }
        }
    }
    let [Some(lower), Some(upper)] = caps else {
        return Err(unproven("prism certificate is missing a cap"));
    };
    let side = side.ok_or_else(|| unproven("prism certificate is missing its side"))?;
    let ordered = [lower, upper, side];
    for (face, expected) in ordered.iter().zip(origins) {
        if expected.is_empty() || model.face_origins(*face) != expected.as_slice() {
            return Err(unproven(
                "prism face origins do not match the immutable recipe",
            ));
        }
    }
    let mut rings = [None, None];
    for index in 0..2 {
        let (edge, cap_use) = ring_use(model, ordered[index], 0)?;
        let natural = analytic_ring(model, edge, original, &axis, axial[index], dec)?;
        if natural * sense(model.edge(edge).sense) * sense(cap_use) != [-1, 1][index] {
            return Err(unproven("cap loop runs against its material orientation"));
        }
        rings[index] = Some((edge, natural));
    }
    let [
        Some((lower_edge, lower_natural)),
        Some((upper_edge, upper_natural)),
    ] = rings
    else {
        return Err(unproven("prism ring inventory is incomplete"));
    };
    if lower_edge == upper_edge {
        return Err(unproven("distinct prism caps share one physical ring"));
    }
    let mut used = [false, false];
    for index in 0..2 {
        let (edge, coedge) = ring_use(model, side, index)?;
        let slot = if edge == lower_edge {
            0
        } else if edge == upper_edge {
            1
        } else {
            return Err(unproven("side trim is detached from both cap rings"));
        };
        if std::mem::replace(&mut used[slot], true)
            || [lower_natural, upper_natural][slot] * sense(model.edge(edge).sense) * sense(coedge)
                != [1, -1][slot]
        {
            return Err(unproven(
                "side/cap correspondence has duplicate or misoriented rings",
            ));
        }
    }
    Ok(PrismTopology {
        faces: ordered,
        rings: [lower_edge, upper_edge],
    })
}

pub(super) fn build_prism(
    model: &mut Model,
    spline: &SplineSurface,
    origins: &[Vec<Origin>; 3],
    marks: MarkLog,
) -> BodyId {
    let rings = [SplineBoundary::VMin, SplineBoundary::VMax].map(|side| {
        let curve = model
            .geom_mut()
            .add_curve(Curve::BSpline(BSplineCurve::boundary(
                spline.carrier(),
                side,
            )));
        model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        })
    });
    let Surface::EllipticCylinder(original) = spline.definition().original() else {
        unreachable!("checked elliptic sweep retains an elliptic source");
    };
    let caps = std::array::from_fn::<_, 2, _>(|index| {
        crate::extrude::face(
            model,
            Surface::Plane(Plane {
                origin: original
                    .origin
                    .plus(&spline.axis().scaled(spline.axial()[index])),
                normal: spline.axis().clone(),
                ref_dir: original.major_dir.clone(),
            }),
            [Sense::Reversed, Sense::Same][index],
            vec![vec![(rings[index], [Sense::Reversed, Sense::Same][index])]],
            None,
            origins[index].clone(),
        )
    });
    let side = crate::extrude::face(
        model,
        Surface::BSpline(spline.clone().into()),
        Sense::Same,
        vec![
            vec![(rings[0], Sense::Same)],
            vec![(rings[1], Sense::Reversed)],
        ],
        None,
        origins[2].clone(),
    );
    let shell = model.add_shell(Shell {
        faces: vec![caps[0], caps[1], side],
    });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks,
    })
}
