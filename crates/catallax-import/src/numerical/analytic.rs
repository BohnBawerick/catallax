//! Retained analytic source geometry with rebuilt, checked finite parameters.

use super::*;
use catallax_numbers::Ball;
use catallax_shapes::{Circle, Cylinder};
pub(super) mod constructed;
mod toroids;

const BITS: u32 = 128;
const OPTIONAL_EXACT_WORK: u64 = 1_000_000;

#[cfg(test)]
use crate::incidence_features::common as fixtures;

pub(super) fn eligible(error: &ImportRefusal) -> bool {
    let ImportRefusal::PlanarRepairNotVerified { detail } = error else {
        return false;
    };
    match detail.stage.as_str() {
        "curved_dependency" | "source_cover" | "vertex_constraints" | "shared_edge"
        | "plane_refit" => true,
        // A missing certificate can take the numerical route. A measured,
        // violated neighbouring-feature limit remains an honest refusal.
        "nearest_feature" => detail.feature_size_squared_mm2.is_none(),
        _ => false,
    }
}

pub(super) fn coordinates(
    p: &Vec3,
    frame: &crate::Placement,
    dec: &mut Decider,
) -> Result<[Ball; 3], ImportError> {
    use crate::planar_arithmetic as a;
    let delta = a::sub(p, &frame.origin, dec)?;
    let y = a::cross(&frame.axis, &frame.ref_dir, dec)?;
    let mut result = Vec::new();
    for axis in [&frame.ref_dir, &y, &frame.axis] {
        let norm = Ball::from_rational(&a::dot(axis, axis, dec)?, BITS)
            .sqrt()
            .ok_or_else(|| unsupported("analytic source frame has no positive norm"))?;
        let value = Ball::from_rational(&a::dot(&delta, axis, dec)?, BITS);
        result.push(
            value
                .checked_div(&norm)
                .ok_or_else(|| unsupported("analytic source frame is degenerate"))?,
        );
    }
    Ok(result.try_into().unwrap())
}

pub(super) fn angle(x: &Ball, y: &Ball) -> Result<Rational, ImportError> {
    if y.lo().is_zero() && y.hi().is_zero() {
        if x.lo().sign() == Sign::Positive {
            return Ok(Rational::zero());
        }
        if x.hi().sign() == Sign::Negative {
            return Ok(Ball::pi(BITS).mid());
        }
    }
    let value =
        Ball::atan2(y, x).ok_or_else(|| unsupported("analytic parameter reaches its axis"))?;
    if value.hi() - value.lo() > Rational::from_ratio(1, 1 << 30).unwrap() {
        return Err(unsupported(
            "analytic angle crosses an unresolved inverse-chart seam",
        ));
    }
    Ok(value.mid())
}

pub(super) fn turn() -> Rational {
    Ball::pi(BITS).mid() * Rational::from_i64(2)
}

fn near(
    mut value: Rational,
    previous: &Rational,
    dec: &mut Decider,
) -> Result<Rational, ImportError> {
    let period = turn();
    let half = Ball::pi(BITS).mid();
    for _ in 0..4 {
        let delta = crate::planar_arithmetic::arithmetic(&value, previous, '-', dec)?;
        if delta > half {
            value = crate::planar_arithmetic::arithmetic(&value, &period, '-', dec)?;
        } else if delta < -&half {
            value = crate::planar_arithmetic::arithmetic(&value, &period, '+', dec)?;
        } else {
            return Ok(value);
        }
    }
    Err(unsupported(
        "analytic boundary requires a distant periodic chart",
    ))
}

fn surface(face: &crate::Face) -> Result<(Surface, &crate::Placement), ImportError> {
    Ok(match &face.surface {
        crate::Surface::Plane(frame) => (
            Surface::Plane(Plane {
                origin: point(&frame.origin),
                normal: vector(&frame.axis),
                ref_dir: vector(&frame.ref_dir),
            }),
            frame,
        ),
        crate::Surface::Cylinder(frame, radius) => (
            Surface::Cylinder(Cylinder {
                origin: point(&frame.origin),
                axis: vector(&frame.axis),
                radius_sq: Real::from(radius * radius),
                ref_dir: vector(&frame.ref_dir),
            }),
            frame,
        ),
        crate::Surface::Cone(frame, ..) => (cones::surface(face)?, frame),
        crate::Surface::Torus(frame, major, minor) => (
            Surface::Torus(catallax_shapes::Torus {
                center: point(&frame.origin),
                axis: vector(&frame.axis),
                major: major.clone().into(),
                minor: minor.clone().into(),
                ref_dir: vector(&frame.ref_dir),
            }),
            frame,
        ),
        _ => {
            return Err(unsupported(
                "the retained analytic source family requires planes, cylinders, fitted cones or ring tori",
            ));
        }
    })
}

fn parameters(
    p: &Vec3,
    face: &crate::Face,
    dec: &mut Decider,
) -> Result<[Rational; 2], ImportError> {
    if let crate::Surface::Torus(frame, major, _) = &face.surface {
        return toroids::parameters(p, frame, major, dec);
    }
    if matches!(face.surface, crate::Surface::Cone(..)) {
        return cones::parameters(p, face, dec);
    }
    let (_, frame) = surface(face)?;
    let [x, y, z] = coordinates(p, frame, dec)?;
    Ok(if matches!(face.surface, crate::Surface::Plane(_)) {
        [x.mid(), y.mid()]
    } else {
        [angle(&x, &y)?, z.mid()]
    })
}

fn add_vertex(
    definition: &mut n::NumericalDefinition,
    indices: &mut BTreeMap<String, usize>,
    source: &crate::Vertex,
    dec: &mut Decider,
) -> Result<usize, ImportError> {
    dec.charge_work(1)?;
    if let Some(index) = indices.get(&source.name) {
        if definition.vertices[*index].position != source.point {
            return Err(unsupported(
                "one source vertex has inconsistent definitions",
            ));
        }
        return Ok(*index);
    }
    let index = definition.vertices.len();
    indices.insert(source.name.clone(), index);
    definition.vertices.push(n::NumericalVertex {
        key: source.name.clone(),
        position: source.point.clone(),
    });
    Ok(index)
}

fn curve(
    edge: &crate::OrientedEdge,
    dec: &mut Decider,
) -> Result<(Curve, [Rational; 2]), ImportError> {
    use crate::planar_arithmetic as a;
    match &edge.curve {
        crate::Curve::Line(origin, dir) => {
            let norm = a::dot(dir, dir, dec)?;
            if norm.sign() != Sign::Positive {
                return Err(unsupported("source line has a zero direction"));
            }
            let mut ts = Vec::new();
            for p in [&edge.start.point, &edge.end.point] {
                let delta = a::sub(p, origin, dec)?;
                ts.push(a::arithmetic(&a::dot(&delta, dir, dec)?, &norm, '/', dec)?);
            }
            Ok((
                Curve::Line(Line {
                    origin: point(origin),
                    dir: vector(dir),
                }),
                ts.try_into().unwrap(),
            ))
        }
        crate::Curve::Circle(frame, radius) => {
            if edge.start.point == edge.end.point && edge.start.name != edge.end.name {
                return Err(unsupported(
                    "coincident circle endpoints need the same explicit source vertex",
                ));
            }
            let mut ts = Vec::new();
            for p in [&edge.start.point, &edge.end.point] {
                let [x, y, _] = coordinates(p, frame, dec)?;
                ts.push(angle(&x, &y)?);
            }
            if edge.curve_same_sense {
                if ts[1] <= ts[0] {
                    ts[1] = &ts[1] + &turn();
                }
            } else if ts[1] >= ts[0] {
                ts[1] = &ts[1] - &turn();
            }
            Ok((
                Curve::Circle(Circle {
                    center: point(&frame.origin),
                    normal: vector(&frame.axis),
                    radius_sq: Real::from(radius * radius),
                    ref_dir: vector(&frame.ref_dir),
                }),
                ts.try_into().unwrap(),
            ))
        }
        _ => Err(unsupported(
            "the retained analytic source family requires finite lines and circular arcs",
        )),
    }
}

fn definition(
    solid: &Solid,
    repairs: &mut Vec<NumericalTopologyRepair>,
    dec: &mut Decider,
) -> Result<n::NumericalDefinition, ImportError> {
    let mut result = n::NumericalDefinition {
        name: solid.name.clone(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::<String, usize>::new();
    let mut source_curves = BTreeMap::new();
    for face in &solid.faces {
        dec.charge_work(1)?;
        if let crate::Surface::Retained(source) = &face.surface {
            toroids::retained_refusal(source, dec)?;
        }
        let (carrier, _) = surface(face)?;
        let mut loops = Vec::new();
        let mut bounds: Vec<_> = face.loops.iter().collect();
        bounds.sort_by_key(|bound| !bound.outer);
        for bound in bounds {
            if bound.vertex.is_some() {
                return Err(unsupported(
                    "analytic vertex-loop faces need separate singular admission",
                ));
            }
            for chain in crate::build::chains(&face.name, bound.orientation, &bound.edges)? {
                let mut trims = Vec::<n::NumericalTrim>::new();
                for (source, forward) in chain {
                    let ends = [
                        add_vertex(&mut result, &mut vertices, &source.start, dec)?,
                        add_vertex(&mut result, &mut vertices, &source.end, dec)?,
                    ];
                    let (curve, ts) = curve(source, dec)?;
                    let index = if let Some(index) = edges.get(&source.edge) {
                        let prior = &result.edges[*index];
                        if prior.vertices != ends
                            || prior.parameters != ts
                            || source_curves.get(&source.edge) != Some(&source.curve)
                        {
                            return Err(unsupported(
                                "one shared source edge has inconsistent geometry",
                            ));
                        }
                        *index
                    } else {
                        let index = result.edges.len();
                        result.edges.push(n::NumericalEdge {
                            key: source.edge.clone(),
                            curve,
                            vertices: ends,
                            parameters: ts,
                        });
                        edges.insert(source.edge.clone(), index);
                        source_curves.insert(source.edge.clone(), source.curve.clone());
                        index
                    };
                    let (a, b) = if forward {
                        (&source.start.point, &source.end.point)
                    } else {
                        (&source.end.point, &source.start.point)
                    };
                    let mut start = parameters(a, face, dec)?;
                    let mut end = parameters(b, face, dec)?;
                    if matches!(face.surface, crate::Surface::Torus(..)) {
                        if let Some(previous) = trims.last() {
                            for axis in 0..2 {
                                start[axis] =
                                    near(start[axis].clone(), &previous.uv[1][axis], dec)?;
                            }
                        }
                        end = toroids::unwrap_end(
                            &result.edges[index],
                            face,
                            forward,
                            &start,
                            end,
                            dec,
                        )?;
                    } else if matches!(
                        face.surface,
                        crate::Surface::Cylinder(..) | crate::Surface::Cone(..)
                    ) {
                        if let Some(previous) = trims.last() {
                            start[0] = near(start[0].clone(), &previous.uv[1][0], dec)?;
                        }
                        let edge = &result.edges[index];
                        if edge.vertices[0] == edge.vertices[1] {
                            end[0] = &start[0] + &rim_delta(source, face, forward, dec)?;
                        } else {
                            let t =
                                (&edge.parameters[0] + &edge.parameters[1]) / Rational::from_i64(2);
                            let p = edge
                                .curve
                                .evaluate(&Ball::from_rational(&t, BITS), dec)?
                                .ok_or_else(|| unsupported("analytic edge has no evaluator"))?;
                            let p = p.each_ref().map(Ball::mid);
                            let middle =
                                near(parameters(&p, face, dec)?[0].clone(), &start[0], dec)?;
                            end[0] = near(end[0].clone(), &middle, dec)?;
                        }
                    }
                    trims.push(n::NumericalTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        edge: index,
                        sense: if forward {
                            Sense::Same
                        } else {
                            Sense::Reversed
                        },
                        uv: [start, end],
                    });
                }
                loops.push(trims);
            }
        }
        if matches!(
            face.surface,
            crate::Surface::Cylinder(..) | crate::Surface::Cone(..)
        ) {
            join_rims(&mut result, &mut loops, solid, face, repairs, dec)?;
        }
        if matches!(face.surface, crate::Surface::Torus(..)) {
            toroids::join_rims(&mut result, &mut loops, solid, face, repairs, dec)?;
        }
        result.faces.push(n::NumericalFace {
            key: face.name.clone(),
            surface: carrier,
            sense: if face.same_sense {
                Sense::Same
            } else {
                Sense::Reversed
            },
            loops,
        });
    }
    Ok(result)
}

// This first periodic family has coaxial complete circles and existing anchors
// on one cylinder generatrix. It needs neither projection nor vertex splitting.
fn rim_delta(
    edge: &crate::OrientedEdge,
    face: &crate::Face,
    forward: bool,
    dec: &mut Decider,
) -> Result<Rational, ImportError> {
    use crate::planar_arithmetic as a;
    if matches!(face.surface, crate::Surface::Cone(..)) {
        return cones::rim_delta(edge, face, forward, dec);
    }
    let (crate::Curve::Circle(circle, r), crate::Surface::Cylinder(cylinder, radius)) =
        (&edge.curve, &face.surface)
    else {
        return Err(unsupported(
            "a periodic cylindrical rim must be a complete circle",
        ));
    };
    let offset = a::sub(&circle.origin, &cylinder.origin, dec)?;
    let parallel = a::cross(&circle.axis, &cylinder.axis, dec)?;
    let centered = a::cross(&offset, &cylinder.axis, dec)?;
    if r != radius || parallel.iter().chain(&centered).any(|q| !q.is_zero()) {
        return Err(unsupported(
            "periodic rims require retained coaxial equal-radius circles",
        ));
    }
    let aligned = a::dot(&circle.axis, &cylinder.axis, dec)?.sign() == Sign::Positive;
    Ok(if (edge.curve_same_sense == forward) == aligned {
        turn()
    } else {
        -turn()
    })
}

fn join_rims(
    definition: &mut n::NumericalDefinition,
    loops: &mut Vec<Vec<n::NumericalTrim>>,
    solid: &Solid,
    face: &crate::Face,
    repairs: &mut Vec<NumericalTopologyRepair>,
    dec: &mut Decider,
) -> Result<(), ImportError> {
    use crate::planar_arithmetic as a;
    if loops.len() != 2 || loops.iter().any(|boundary| boundary.len() != 1) {
        return Ok(());
    }
    dec.charge_work(1)?;
    let first = loops[0][0].clone();
    let mut second = loops[1][0].clone();
    let rims = [
        &definition.edges[first.edge],
        &definition.edges[second.edge],
    ];
    if rims
        .iter()
        .any(|e| e.vertices[0] != e.vertices[1] || !matches!(e.curve, Curve::Circle(_)))
    {
        return Ok(());
    }
    let vertices = [rims[0].vertices[0], rims[1].vertices[0]];
    let start = &definition.vertices[vertices[0]].position;
    let end = &definition.vertices[vertices[1]].position;
    let delta = a::sub(end, start, dec)?;
    let same_generator = match &face.surface {
        crate::Surface::Cylinder(frame, _) => a::cross(&delta, &frame.axis, dec)?
            .iter()
            .all(Rational::is_zero),
        crate::Surface::Cone(..) => cones::same_generator(face, start, end)?,
        _ => unreachable!(),
    };
    if delta.iter().all(Rational::is_zero)
        || !same_generator
        || &first.uv[1][0] - &first.uv[0][0] != &second.uv[0][0] - &second.uv[1][0]
    {
        return Err(unsupported(
            "periodic rim anchors must share one generatrix with opposite winding",
        ));
    }
    let shift = &first.uv[1][0] - &second.uv[0][0];
    for uv in &mut second.uv {
        uv[0] = &uv[0] + &shift;
    }
    let seam = definition.edges.len();
    let key = format!("derived:{}:periodic-seam", face.name);
    repairs.push(NumericalTopologyRepair {
        solid: solid.name.clone(),
        face: face.name.clone(),
        source_rim_edges: rims.map(|edge| edge.key.clone()),
        source_vertices: vertices.map(|v| definition.vertices[v].key.clone()),
        derived_seam_edge: key.clone(),
        method: if matches!(face.surface, crate::Surface::Cone(..)) {
            "paired_cone_generatrix/1"
        } else {
            "paired_cylinder_generatrix/1"
        }
        .into(),
    });
    definition.edges.push(n::NumericalEdge {
        key,
        curve: Curve::Line(Line {
            origin: point(start),
            dir: vector(&delta),
        }),
        parameters: [Rational::zero(), Rational::one()],
        vertices,
    });
    *loops = vec![vec![
        first.clone(),
        n::NumericalTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            edge: seam,
            sense: Sense::Same,
            uv: [first.uv[1].clone(), second.uv[0].clone()],
        },
        second.clone(),
        n::NumericalTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            edge: seam,
            sense: Sense::Reversed,
            uv: [second.uv[1].clone(), first.uv[0].clone()],
        },
    ]];
    Ok(())
}

fn charge_face(face: &crate::Face, dec: &mut Decider) -> Result<(), ImportError> {
    dec.charge_work(1)?;
    for boundary in &face.loops {
        dec.charge_work(1 + boundary.edges.len() as u64)?;
    }
    Ok(())
}

#[derive(Debug)]
struct SourceChecks {
    pairs: Vec<NumericalSourcePairCheck>,
    optional_work_reserved: u64,
    optional_work_used: u64,
}

fn polygonal(face: &crate::Face) -> bool {
    matches!(face.surface, crate::Surface::Plane(_))
        && !face.loops.is_empty()
        && face
            .loops
            .iter()
            .flat_map(|boundary| &boundary.edges)
            .all(|edge| matches!(edge.curve, crate::Curve::Line(..)))
}

fn unresolved_pair(
    solid: &str,
    faces: [String; 2],
    refusal: ImportRefusal,
    work: u64,
) -> NumericalSourcePairCheck {
    NumericalSourcePairCheck::Unresolved {
        solid: solid.to_owned(),
        faces,
        refusal: Box::new(refusal),
        exact_attempt_work_used: work,
    }
}

fn optional_unresolved(
    solid: &str,
    faces: [String; 2],
    allowance: u64,
    work: u64,
) -> NumericalSourcePairCheck {
    let refusal = Refusal::new(
        "import.source-proof-work-limit",
        format!(
            "optional curved-pair proof was not completed within the collective {allowance}-work allowance"
        ),
    );
    unresolved_pair(
        solid,
        faces.clone(),
        ImportRefusal::SelfIntersectionNotVerified {
            solid: solid.to_owned(),
            faces,
            reason: refusal.to_string(),
            refusal: Some(Box::new(refusal)),
        },
        work,
    )
}

fn source_pair(
    solid: &str,
    a: &crate::Face,
    b: &crate::Face,
    dec: &mut Decider,
) -> Result<(), ImportError> {
    // Charge retained tables before scratch setup. Existing exact arithmetic
    // remains on the active caller/optional-proof counter.
    charge_face(a, dec)?;
    charge_face(b, dec)?;
    crate::intersect::check_solid(&Solid {
        name: solid.to_owned(),
        faces: vec![a.clone(), b.clone()],
    })
    .map_err(Into::into)
}

/// All cheap checks run on the shared request before optional proofs. One gap
/// or expensive proof therefore cannot suppress decidable pairs in later solids.
fn intersection_checks(
    file: &crate::StepFile,
    optional_allowance: u64,
    dec: &mut Decider,
) -> Result<SourceChecks, ImportError> {
    let mut checks = Vec::new();
    let mut optional = Vec::new();
    for solid in &file.solids {
        let mut incidence = Vec::new();
        let mut bounds = Vec::new();
        for face in &solid.faces {
            charge_face(face, dec)?;
            incidence.push(crate::intersect::planar_incidence_error(face));
            bounds.push(crate::curved::face_bounds(face));
            dec.check_deadline()?;
        }
        for (i, a) in solid.faces.iter().enumerate() {
            for (j, b) in solid.faces.iter().enumerate().skip(i + 1) {
                dec.charge_work(1)?;
                let faces = [a.name.clone(), b.name.clone()];
                if let Some(refusal) = incidence[i].as_ref().or(incidence[j].as_ref()) {
                    checks.push(unresolved_pair(&solid.name, faces, refusal.clone(), 0));
                } else if let (Some(a), Some(b)) = (&bounds[i], &bounds[j])
                    && crate::curved::apart(a, b)
                {
                    checks.push(NumericalSourcePairCheck::CheckedNoCrossing {
                        solid: solid.name.clone(),
                        faces,
                        exact_attempt_work_used: 0,
                    });
                } else if polygonal(a) && polygonal(b) {
                    let before = catallax_numbers::current_operation_work().used;
                    source_pair(&solid.name, a, b, dec)?;
                    dec.check_deadline()?;
                    checks.push(NumericalSourcePairCheck::CheckedNoCrossing {
                        solid: solid.name.clone(),
                        faces,
                        exact_attempt_work_used: catallax_numbers::current_operation_work().used
                            - before,
                    });
                } else {
                    optional.push((checks.len(), solid.name.as_str(), a, b));
                    checks.push(optional_unresolved(
                        &solid.name,
                        faces,
                        optional_allowance,
                        0,
                    ));
                }
            }
        }
    }
    if optional.is_empty() || optional_allowance == 0 {
        return Ok(SourceChecks {
            pairs: checks,
            optional_work_reserved: 0,
            optional_work_used: 0,
        });
    }
    // A single collective prepayment bounds costly curved proofs. Cheap pairs
    // above pay only their actual shared work, regardless of their count.
    let allowance = catallax_numbers::OperationAllowance::reserve(optional_allowance)?;
    let (result, work) = allowance.run(|| -> Result<(), ImportError> {
        for (slot, solid, a, b) in optional {
            if catallax_numbers::current_operation_work().exhausted {
                break;
            }
            dec.check_deadline()?;
            let before = catallax_numbers::current_operation_work().used;
            let result = source_pair(solid, a, b, dec);
            let after = catallax_numbers::current_operation_work();
            let used = after.used - before;
            let faces = [a.name.clone(), b.name.clone()];
            if after.exhausted {
                checks[slot] = optional_unresolved(solid, faces, optional_allowance, used);
                break;
            }
            checks[slot] = match result {
                Ok(()) => NumericalSourcePairCheck::CheckedNoCrossing {
                    solid: solid.to_owned(),
                    faces,
                    exact_attempt_work_used: used,
                },
                Err(ImportError::Source(
                    refusal @ ImportRefusal::SelfIntersectionNotVerified { .. },
                )) => unresolved_pair(solid, faces, refusal, used),
                Err(error) => return Err(error),
            };
        }
        Ok(())
    });
    dec.check_deadline()?;
    result?;
    Ok(SourceChecks {
        pairs: checks,
        optional_work_reserved: optional_allowance,
        optional_work_used: work.used,
    })
}

#[cfg(test)]
#[test]
fn exhausted_optional_proof_is_unresolved_and_paid_without_exhausting_the_request() {
    use catallax_numbers::{current_operation_work, with_operation_budget};
    let mut source = fixtures::Soup::new();
    source.cylinder(["0", "0", "0"], 2, "2", "1", fixtures::Caps::Both);
    let source = source.finish();
    let file = crate::read_step(source.as_bytes()).unwrap();
    let ((checks, paid), work) = with_operation_budget(100_000, || {
        let checks = intersection_checks(&file, 1, &mut Decider::default()).unwrap();
        (checks, current_operation_work().used)
    });
    assert!(!work.exhausted && paid >= 1 && paid == work.used);
    assert_eq!(checks.optional_work_reserved, 1);
    assert!(
        checks.pairs.iter().any(|check| matches!(check,
        NumericalSourcePairCheck::Unresolved { refusal, .. }
        if matches!(refusal.as_ref(), ImportRefusal::SelfIntersectionNotVerified {
            refusal: Some(reason), ..
        } if reason.code == "import.source-proof-work-limit"))),
        "{checks:?}"
    );
    let (result, work) =
        with_operation_budget(1, || intersection_checks(&file, 1, &mut Decider::default()));
    assert!(result.is_err() && work.exhausted);
}

#[cfg(test)]
#[test]
fn cheap_crossings_are_checked_before_any_optional_proof_allowance() {
    let mut source = fixtures::Soup::new();
    source.cylinder(["10", "10", "10"], 2, "2", "1", fixtures::Caps::Both);
    source.box_(["0", "0", "0"], ["2", "2", "2"]);
    source.box_(["1", "1", "1"], ["3", "3", "3"]);
    let file = crate::read_step(source.finish().as_bytes()).unwrap();
    let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || {
        intersection_checks(&file, 0, &mut Decider::default())
    });
    assert!(!work.exhausted);
    assert!(matches!(
        result,
        Err(ImportError::Source(ImportRefusal::SelfIntersecting { .. }))
    ));
}

#[cfg(test)]
#[test]
fn many_planar_pairs_do_not_prepay_optional_work() {
    let points = [
        ["0", "0", "0"],
        ["1", "1", "0"],
        ["2", "4", "0"],
        ["3", "9", "0"],
        ["4", "16", "0"],
        ["5", "25", "0"],
        ["6", "36", "0"],
        ["7", "49", "0"],
        ["8", "64", "0"],
        ["9", "81", "0"],
        ["10", "100", "0"],
        ["11", "121", "0"],
        ["12", "144", "0"],
        ["13", "169", "0"],
        ["14", "196", "0"],
        ["15", "225", "0"],
        ["16", "256", "0"],
        ["17", "289", "0"],
        ["18", "324", "0"],
        ["19", "361", "0"],
        ["20", "400", "0"],
        ["21", "441", "0"],
        ["22", "484", "0"],
        ["23", "529", "0"],
        ["24", "576", "0"],
        ["25", "625", "0"],
        ["26", "676", "0"],
        ["27", "729", "0"],
    ];
    let mut faces = vec![points.into_iter().rev().collect()];
    for i in 0..points.len() {
        faces.push(vec![
            points[i],
            points[(i + 1) % points.len()],
            ["13", "260", "1"],
        ]);
    }
    let source = fixtures::polyhedra_step(&[faces]);
    let file = crate::read_step(source.as_bytes()).unwrap();
    assert_eq!(file.solids[0].faces.len(), 29);
    let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || {
        intersection_checks(&file, 50_000, &mut Decider::default())
    });
    assert!(result.is_ok() && !work.exhausted, "{result:?}, {work:?}");
    let checks = result.unwrap();
    assert_eq!(checks.pairs.len(), 406);
    assert_eq!(checks.optional_work_reserved, 0);
}

pub(super) fn import(
    bytes: &[u8],
    verified_sha256: &str,
    refusal: ImportRefusal,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<StepImport, ImportError> {
    if options.tolerance_mm.sign() != Sign::Positive
        || options.tolerance_mm > crate::gap::auto_cap()
    {
        return Err(Refusal::new(
            "import.numerical-tolerance",
            "source tolerance must be positive and at most 0.01 mm",
        )
        .into());
    }
    dec.charge_work(bytes.len() as u64)?;
    let file = crate::read_step(bytes)?;
    admit_definition(&file, verified_sha256, refusal, options, false, dec)
}

pub(super) fn admit_definition(
    file: &crate::StepFile,
    verified_sha256: &str,
    refusal: ImportRefusal,
    options: &n::NumericalOptions,
    permit_planar: bool,
    dec: &mut Decider,
) -> Result<StepImport, ImportError> {
    if !file.surface_models.is_empty() {
        return Err(unsupported(
            "numerical analytic continuation must retain every surface model and its surface body kind",
        ));
    }
    if !permit_planar
        && !file
            .solids
            .iter()
            .flat_map(|solid| &solid.faces)
            .any(|face| {
                matches!(
                    face.surface,
                    crate::Surface::Cylinder(..)
                        | crate::Surface::Cone(..)
                        | crate::Surface::Torus(..)
                ) || face
                    .loops
                    .iter()
                    .flat_map(|bound| &bound.edges)
                    .any(|edge| matches!(edge.curve, crate::Curve::Circle(..)))
            })
    {
        // Keep the previous planar-only policy and its export controls until
        // that separate family receives its own complete qualification.
        return Err(refusal.into());
    }
    let face_tolerances = file
        .solids
        .iter()
        .map(crate::cone::tolerances)
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    if face_tolerances.iter().any(|bound| {
        Rational::parse(&bound.deviation_upper_mm).expect("producer bound") > options.tolerance_mm
    }) {
        return Err(Refusal::new(
            "import.numerical-tolerance",
            "fitted cone deviation exceeds the requested tolerance",
        )
        .into());
    }
    let continuation_options = cones::continuation_options(file, options, dec)?;
    let intersection_checks = intersection_checks(file, OPTIONAL_EXACT_WORK, dec)?;
    let mut bodies = Vec::new();
    let mut topology_repairs = Vec::new();
    for solid in &file.solids {
        bodies.push((
            solid.name.clone(),
            n::build(
                definition(solid, &mut topology_repairs, dec)?,
                continuation_options.clone(),
                dec,
            )?,
        ));
    }
    if bodies.is_empty() {
        return Err(unsupported("source contains no numerical solids"));
    }
    let record = RepairRecord {
        format: crate::REPAIR_FORMAT.into(),
        source: crate::SourceRef {
            name: "source.step".into(),
            sha256: verified_sha256.to_owned(),
        },
        derivative: None,
        changes: Vec::new(),
        face_tolerances,
        source_incidence: Vec::new(),
        cylinder_radius: Vec::new(),
    };
    Ok(StepImport::Numerical(Box::new(NumericalImport {
        bodies,
        record,
        source_admission: Some(NumericalSourceAdmission {
            strict_refusal: refusal,
            tolerance_mm: options.tolerance_mm.clone(),
            method: "retained_shared_analytic_boundaries/3".into(),
            topology_repairs,
            intersection_checks: intersection_checks.pairs,
            optional_exact_work_reserved: intersection_checks.optional_work_reserved,
            optional_exact_work_used: intersection_checks.optional_work_used,
        }),
    })))
}
