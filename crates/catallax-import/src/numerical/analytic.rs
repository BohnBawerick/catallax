//! Retained analytic source geometry with rebuilt, checked finite parameters.

use super::*;
use catallax_numbers::Ball;
use catallax_shapes::{Circle, Cylinder};

const BITS: u32 = 128;

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

fn coordinates(
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

fn angle(x: &Ball, y: &Ball) -> Result<Rational, ImportError> {
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

fn turn() -> Rational {
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
        _ => {
            return Err(unsupported(
                "the retained analytic source family requires planes and cylinders",
            ));
        }
    })
}

fn parameters(
    p: &Vec3,
    face: &crate::Face,
    dec: &mut Decider,
) -> Result<[Rational; 2], ImportError> {
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
            if edge.start.point == edge.end.point {
                return Err(unsupported(
                    "closed circle edges require an explicit periodic seam construction",
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

fn definition(solid: &Solid, dec: &mut Decider) -> Result<n::NumericalDefinition, ImportError> {
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
                    if matches!(face.surface, crate::Surface::Cylinder(..)) {
                        if let Some(previous) = trims.last() {
                            start[0] = near(start[0].clone(), &previous.uv[1][0], dec)?;
                        }
                        let edge = &result.edges[index];
                        let t = (&edge.parameters[0] + &edge.parameters[1]) / Rational::from_i64(2);
                        let p = edge
                            .curve
                            .evaluate(&Ball::from_rational(&t, BITS), dec)?
                            .ok_or_else(|| unsupported("analytic edge has no evaluator"))?;
                        let p = p.each_ref().map(Ball::mid);
                        let middle = near(parameters(&p, face, dec)?[0].clone(), &start[0], dec)?;
                        end[0] = near(end[0].clone(), &middle, dec)?;
                    }
                    trims.push(n::NumericalTrim {
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

pub(super) fn import(
    bytes: &[u8],
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
    if !file
        .solids
        .iter()
        .flat_map(|solid| &solid.faces)
        .any(|face| {
            matches!(face.surface, crate::Surface::Cylinder(..))
                || face
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
    match crate::check_self_intersection(&file) {
        Ok(()) => {}
        // This is the failed exact source-incidence prerequisite being replaced
        // by numerical checks. As in N01, samples establish no global body proof.
        Err(ImportRefusal::NotSolid { detail }) if detail.contains("exactly") => {}
        Err(error) => return Err(error.into()),
    }
    let mut bodies = Vec::new();
    for solid in &file.solids {
        bodies.push((
            solid.name.clone(),
            n::build(definition(solid, dec)?, options.clone(), dec)?,
        ));
    }
    if bodies.is_empty() {
        return Err(unsupported("source contains no numerical solids"));
    }
    let record = RepairRecord {
        format: crate::REPAIR_FORMAT.into(),
        source: crate::SourceRef {
            name: "source.step".into(),
            sha256: crate::sha256_hex(bytes),
        },
        derivative: None,
        changes: Vec::new(),
        face_tolerances: Vec::new(),
        source_incidence: Vec::new(),
        cylinder_radius: Vec::new(),
    };
    Ok(StepImport::Numerical(Box::new(NumericalImport {
        bodies,
        record,
        source_admission: Some(NumericalSourceAdmission {
            strict_refusal: refusal,
            tolerance_mm: options.tolerance_mm.clone(),
            method: "retained_shared_analytic_boundaries/1".into(),
        }),
    })))
}
