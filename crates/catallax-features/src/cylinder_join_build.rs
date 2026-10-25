//! Finite trim witnesses and B-rep assembly for bounded cylinder joins.
use crate::{
    cylinder_join::{CylinderJoin, same_support, unsupported},
    cylinder_skin::CylinderSkin,
    polyhedron::{ConvexCell, PlaneFace},
    profile::sign,
};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, StepId, ThickenSide};
use catallax_shapes::{Circle, Curve, Surface};
use catallax_topology::{Body, BodyId, BodyKind, EdgeId, FaceId, Loop, Lump, Model, Sense, Shell};
use std::collections::BTreeMap;

fn equal(a: &Real, b: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&(a - b), dec)? == Sign::Zero)
}
fn same_ring(a: &Circle, b: &Circle, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(equal(&a.radius_sq, &b.radius_sq, dec)?
        && equal(
            &a.center.minus(&b.center).norm_squared(),
            &Real::zero(),
            dec,
        )?
        && equal(
            &a.normal.cross(&b.normal).norm_squared(),
            &Real::zero(),
            dec,
        )?)
}
/// A recognition succeeds only after every active source trim has matched.
pub(crate) fn skin(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<Option<CylinderJoin>, Refusal> {
    let source = model.faces_of_body(body);
    let cylinders: Vec<_> = source
        .iter()
        .copied()
        .filter(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Cylinder(_)
            )
        })
        .collect();
    if cylinders.is_empty()
        || source.len() == 3
        || source.iter().any(|f| {
            !matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::Cylinder(_) | Surface::Plane(_)
            )
        })
    {
        return Ok(None);
    }
    if !matches!((cylinders.len(), source.len()), (1, 8) | (2, 5))
        || model.body(body).lumps.len() != 1
        || model.lump(model.body(body).lumps[0]).shells.len() != 1
        || source.iter().any(|f| model.face(*f).sense != Sense::Same)
    {
        return Err(unsupported(
            "cylinder joins require one bounded complete boss or two-band shaft skin",
        ));
    }
    let mut incidence: BTreeMap<EdgeId, Vec<(FaceId, Sense)>> = BTreeMap::new();
    for f in &source {
        for lp in &model.face(*f).loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                return Err(unsupported(
                    "singular join trims are outside the bounded lane",
                ));
            };
            for c in cs {
                let c = model.coedge(*c);
                incidence.entry(c.edge).or_default().push((*f, c.sense));
            }
        }
    }
    if incidence
        .values()
        .any(|v| v.len() != 2 || v[0].0 == v[1].0 || v[0].1 == v[1].1)
    {
        return Err(unsupported(
            "join source trims must have two opposite incident uses",
        ));
    }
    let mut bands = Vec::new();
    for f in &cylinders {
        let Surface::Cylinder(cylinder) = model.geom().surface(model.face(*f).surface) else {
            unreachable!()
        };
        if model.face(*f).loops.len() != 2 {
            return Err(unsupported(
                "a cylindrical source patch needs two complete circular rims",
            ));
        }
        let mut ends = Vec::new();
        for lp in &model.face(*f).loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                unreachable!()
            };
            let [coedge] = cs.as_slice() else {
                return Err(unsupported(
                    "partial circular joins need a finite arc arrangement",
                ));
            };
            let edge = model.edge(model.coedge(*coedge).edge);
            let Curve::Circle(c) = model.geom().curve(edge.curve) else {
                return Err(unsupported(
                    "join bands require analytic circular trim witnesses",
                ));
            };
            if edge.start.is_some()
                || edge.end.is_some()
                || !equal(&c.radius_sq, &cylinder.radius_sq, dec)?
                || !equal(
                    &c.normal.cross(&cylinder.axis).norm_squared(),
                    &Real::zero(),
                    dec,
                )?
                || !equal(
                    &c.center
                        .minus(&cylinder.origin)
                        .cross(&cylinder.axis)
                        .norm_squared(),
                    &Real::zero(),
                    dec,
                )?
            {
                return Err(unsupported(
                    "a circular trim is not a complete section of its cylinder",
                ));
            }
            let other = incidence[&model.coedge(*coedge).edge]
                .iter()
                .find(|(other, _)| other != f)
                .unwrap()
                .0;
            if !matches!(
                model.geom().surface(model.face(other).surface),
                Surface::Plane(_)
            ) {
                return Err(unsupported(
                    "a circular rim needs its perpendicular adjacent plane",
                ));
            }
            ends.push((c.center.clone(), model.face_origins(other).to_vec(), other));
        }
        if sign(&ends[1].0.minus(&ends[0].0).dot(&cylinder.axis), dec)? == Sign::Negative {
            ends.swap(0, 1);
        }
        let band = CylinderSkin {
            cylinder: cylinder.clone(),
            ends: [ends[0].0.clone(), ends[1].0.clone()],
            origins: [
                ends[0].1.clone(),
                ends[1].1.clone(),
                model.face_origins(*f).to_vec(),
            ],
        };
        if sign(
            &band.ends[1].minus(&band.ends[0]).dot(&band.cylinder.axis),
            dec,
        )? != Sign::Positive
        {
            return Err(unsupported("a cylinder band needs positive axial extent"));
        }
        for end in &ends {
            let Surface::Plane(p) = model.geom().surface(model.face(end.2).surface) else {
                unreachable!()
            };
            if !equal(
                &p.normal.cross(&cylinder.axis).norm_squared(),
                &Real::zero(),
                dec,
            )? || !equal(&end.0.minus(&p.origin).dot(&p.normal), &Real::zero(), dec)?
            {
                return Err(unsupported(
                    "a join cap must contain the complete perpendicular section",
                ));
            }
            // An internal shoulder may face either way. The end names are
            // retained separately from the synthetic primitive cap orientation.
        }
        bands.push((band, ends[1].2));
    }
    let join = if bands.len() == 2 {
        CylinderJoin::stepped(bands.remove(0).0, bands.remove(0).0, dec)?
    } else {
        let (band, cap) = bands.remove(0);
        let planes = source
            .iter()
            .filter(|f| **f != cap && !cylinders.contains(f))
            .map(|f| {
                let Surface::Plane(p) = model.geom().surface(model.face(*f).surface) else {
                    unreachable!()
                };
                PlaneFace {
                    plane: p.clone(),
                    origins: model.face_origins(*f).to_vec(),
                }
            })
            .collect();
        CylinderJoin::boss(ConvexCell::new(planes, dec)?, band, dec)?
    };
    validate(model, &source, &incidence, &join, dec)?;
    Ok(Some(join))
}
fn validate(
    model: &Model,
    source: &[FaceId],
    incidence: &BTreeMap<EdgeId, Vec<(FaceId, Sense)>>,
    join: &CylinderJoin,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut mapped = Vec::new();
    for f in &join.faces {
        let mut found = Vec::new();
        for source in source {
            if same_support(
                &f.surface,
                model.geom().surface(model.face(*source).surface),
                dec,
            )? {
                found.push(*source);
            }
        }
        let [f] = found.as_slice() else {
            return Err(unsupported(
                "physical join patches need unique oriented supports",
            ));
        };
        mapped.push(*f);
    }
    let mut used = Vec::new();
    for boundary in &join.boundaries {
        let Curve::Circle(expected) = &boundary.curve else {
            unreachable!()
        };
        let mut matched = None;
        for (edge, uses) in incidence {
            if used.contains(edge)
                || !boundary
                    .faces
                    .iter()
                    .all(|i| uses.iter().any(|(f, _)| *f == mapped[*i]))
            {
                continue;
            }
            let e = model.edge(*edge);
            let Curve::Circle(actual) = model.geom().curve(e.curve) else {
                continue;
            };
            if e.start.is_some() || e.end.is_some() || !same_ring(expected, actual, dec)? {
                continue;
            }
            let same = sign(&expected.normal.dot(&actual.normal), dec)? == Sign::Positive;
            for (i, index) in boundary.faces.iter().enumerate() {
                let sense = uses.iter().find(|(f, _)| *f == mapped[*index]).unwrap().1;
                if (sense == e.sense) != (same == (i == 1)) {
                    return Err(unsupported(
                        "a circular trim has the wrong material-side orientation",
                    ));
                }
            }
            matched = Some(*edge);
            break;
        }
        used.push(matched.ok_or_else(|| {
            unsupported("source circular trim adjacency differs from the finite join witness")
        })?);
    }
    if let Some((cell, _)) = &join.plate {
        for (i, plane) in cell.faces.iter().enumerate() {
            let mut found = None;
            for f in source {
                if same_support(
                    &Surface::Plane(plane.plane.clone()),
                    model.geom().surface(model.face(*f).surface),
                    dec,
                )? {
                    found = Some(*f);
                }
            }
            let f = found.ok_or_else(|| unsupported("a plate face is missing"))?;
            let polygon = &cell.polygons[i];
            let mut lines = Vec::new();
            for lp in &model.face(f).loops {
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    unreachable!()
                };
                if cs.len() == 1 && used.contains(&model.coedge(cs[0]).edge) {
                    continue;
                }
                if cs.len() != polygon.len() {
                    return Err(unsupported(
                        "a plate trim is not the complete source rectangle",
                    ));
                }
                for c in cs {
                    let c = model.coedge(*c);
                    let e = model.edge(c.edge);
                    if !matches!(model.geom().curve(e.curve), Curve::Line(_)) {
                        return Err(unsupported(
                            "plate outer trims require exact straight segments",
                        ));
                    }
                    let (Some(a), Some(b)) = (e.start, e.end) else {
                        return Err(unsupported("a plate segment needs finite endpoints"));
                    };
                    let (a, b) = if c.sense == Sense::Same {
                        (a, b)
                    } else {
                        (b, a)
                    };
                    lines.push((
                        model.geom().position(model.vertex(a).point),
                        model.geom().position(model.vertex(b).point),
                    ));
                    used.push(c.edge);
                }
            }
            if lines.len() != polygon.len() {
                return Err(unsupported("a plate face has extra or missing outer loops"));
            }
            for j in 0..polygon.len() {
                let a = &cell.vertices()[polygon[j]];
                let b = &cell.vertices()[polygon[(j + 1) % polygon.len()]];
                let mut found = false;
                for (x, y) in &lines {
                    if equal(&a.minus(x).norm_squared(), &Real::zero(), dec)?
                        && equal(&b.minus(y).norm_squared(), &Real::zero(), dec)?
                    {
                        found = true;
                        break;
                    }
                }
                if !found {
                    return Err(unsupported(
                        "a plate's directed finite segment differs from its rectangle witness",
                    ));
                }
            }
        }
    }
    used.sort();
    used.dedup();
    if used.len() != incidence.len() {
        return Err(unsupported("join source contains an unmatched active trim"));
    }
    Ok(())
}
pub(crate) fn outward(
    model: &mut Model,
    join: &CylinderJoin,
    d: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let (faces, boundaries) = join.outward(d, step, dec)?;
    let circular = crate::cone_build::sections(model, &faces, &boundaries, None);
    let Some((cell, top)) = &join.plate else {
        return Ok(circular);
    };
    let rounded = crate::rounded_polyhedron::build(model, cell, d, step, Default::default(), dec)?;
    let shifted = Surface::Plane(cell.faces[*top].plane.clone())
        .parallel_support(d, dec)?
        .unwrap();
    let mut combined = Vec::new();
    let circular_faces = model.faces_of_body(circular);
    for f in model.faces_of_body(rounded) {
        if same_support(&shifted, model.geom().surface(model.face(f).surface), dec)? {
            let paths = model
                .face(f)
                .loops
                .iter()
                .chain(model.face(circular_faces[0]).loops.iter())
                .map(|lp| {
                    let Loop::Coedges(cs) = model.loop_(*lp) else {
                        unreachable!()
                    };
                    cs.iter()
                        .map(|c| {
                            let c = model.coedge(*c);
                            (c.edge, c.sense)
                        })
                        .collect()
                })
                .collect();
            combined.push(crate::extrude::face(
                model,
                shifted.clone(),
                Sense::Same,
                paths,
                None,
                model.face_origins(f).to_vec(),
            ));
        } else {
            combined.push(f);
        }
    }
    combined.extend(circular_faces.into_iter().skip(1));
    let shell = model.add_shell(Shell { faces: combined });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    }))
}
pub(crate) fn cylinder(model: &mut Model, skin: &CylinderSkin) -> BodyId {
    let faces = skin
        .surfaces()
        .into_iter()
        .zip(skin.origins.clone())
        .map(|(surface, origins)| crate::cone_skin::ConeFace { surface, origins })
        .collect::<Vec<_>>();
    let boundaries =
        [(0, [0, 2]), (1, [2, 1])].map(|(end, faces)| crate::cone_skin::ConeBoundary {
            curve: Curve::Circle(Circle {
                center: skin.ends[end].clone(),
                normal: skin.cylinder.axis.clone(),
                radius_sq: skin.cylinder.radius_sq.clone(),
                ref_dir: skin.cylinder.ref_dir.clone(),
            }),
            faces,
        });
    crate::cone_build::sections(model, &faces, &boundaries, None)
}
pub(crate) fn thicken(
    model: &mut Model,
    join: &CylinderJoin,
    refs: &[FaceRef],
    t: &Real,
    side: ThickenSide,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let index = join.select(refs)?;
    if let Some((prism, hole)) = join.plate_sweep(index, t, side, step, dec)? {
        let prism = prism.build(model);
        let hole = cylinder(model, &hole);
        return Ok(catallax_boolean::boolean(
            model,
            prism,
            hole,
            catallax_boolean::BoolOp::Subtract,
            dec,
        )
        .map_err(|e| e.refusal)?
        .body);
    }
    let band = join.face_band(index, t, side, step, dec)?;
    let outer = cylinder(model, &band.outer);
    if let Some(inner) = band.inner {
        let inner = cylinder(model, &inner);
        Ok(
            catallax_boolean::boolean(model, outer, inner, catallax_boolean::BoolOp::Subtract, dec)
                .map_err(|e| e.refusal)?
                .body,
        )
    } else {
        Ok(outer)
    }
}
