//! A bounded cut chart for a spline blend and exact ribbons for its sides.

use super::rational_patch::{Definition, Support};
use super::*;
use catallax_numbers::{Outcome, ZeroProof};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, CutCurve, Point3, RationalPatch, SplineBoundary, SplineFillet,
};
use catallax_topology::{Loop, RationalPatchFace};

mod grid;
pub(super) use grid::Mapping;

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("mesh.spline-fillet-unproven", detail)
}

fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(&a.minus(b).norm_squared()) {
        Outcome::Proven(_) => Ok(false),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => Ok(true),
        other => Err(bad("boundary identity needs an exact zero").keeping_undecided(other)),
    }
}

fn coedges(model: &Model, face: FaceId) -> Result<Vec<catallax_topology::CoedgeId>, Refusal> {
    let [lp] = model.face(face).loops.as_slice() else {
        return Err(bad("a spline fillet chart requires one loop"));
    };
    let Loop::Coedges(edges) = model.loop_(*lp) else {
        return Err(bad("a spline fillet chart requires physical edges"));
    };
    if edges.len() != 4 {
        return Err(bad(format!(
            "a spline fillet chart requires 4 edges, got {}",
            edges.len()
        )));
    }
    let mut seen = std::collections::BTreeSet::new();
    let mut cycle = Vec::new();
    for &c in edges {
        let c = model.coedge(c);
        let e = model.edge(c.edge);
        let (Some(a), Some(b)) = (e.start, e.end) else {
            return Err(bad("a spline fillet chart needs finite edge ends"));
        };
        if !seen.insert(c.edge) {
            return Err(bad("a spline fillet chart repeats an edge"));
        }
        cycle.push(if c.sense == Sense::Same {
            (a, b)
        } else {
            (b, a)
        });
    }
    if (0..4).any(|i| cycle[i].1 != cycle[(i + 1) % 4].0) {
        return Err(bad("a spline fillet chart has an unsewn vertex cycle"));
    }
    Ok(edges.clone())
}

fn patch_of(model: &Model, face: FaceId) -> Option<&RationalPatch> {
    match model.geom().surface(model.face(face).surface) {
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => Some(patch),
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) if atlas.charts().len() == 1 => {
            Some(atlas.charts()[0].patch())
        }
        _ => None,
    }
}

fn plane_of(model: &Model, face: FaceId) -> Result<([Rational; 3], Rational), Refusal> {
    let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
        return Err(bad("a spline fillet cut needs a planar neighbour"));
    };
    let [nx, ny, nz] = [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(Real::as_rational);
    let normal =
        [nx, ny, nz].map(|v| v.ok_or_else(|| bad("a spline fillet cap needs a rational normal")));
    let [nx, ny, nz] = normal;
    let normal = [nx?, ny?, nz?];
    let offset = plane
        .normal
        .dot(&catallax_shapes::Vector3::new(
            plane.origin.x.clone(),
            plane.origin.y.clone(),
            plane.origin.z.clone(),
        ))
        .as_rational()
        .ok_or_else(|| bad("a spline fillet cap needs a rational offset"))?;
    Ok((normal, offset))
}

pub(super) fn prepare(
    model: &Model,
    faces: &[FaceId],
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<Junction, Refusal> {
    catallax_numbers::reuse_operation_enclosures();
    let blends: Vec<_> = faces
        .iter()
        .filter_map(
            |&face| match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(SplineSupport::SplineFillet(fillet)) => Some((face, fillet)),
                _ => None,
            },
        )
        .collect();
    let [(blend, fillet)] = blends.as_slice() else {
        return Err(bad(format!(
            "one spline blend is supported, got {}",
            blends.len()
        )));
    };
    let (blend, fillet) = (*blend, *fillet);
    let borders = fillet.borders();
    if borders.len() != 4 || (borders[0] < borders[1]) != (borders[2] < borders[3]) {
        return Err(bad(format!(
            "four consistently ordered rail borders required, got {borders:?}"
        )));
    }
    let order = if borders[0] < borders[1] {
        [0, 1]
    } else {
        [1, 0]
    };
    let ranges: [[Rational; 2]; 2] =
        std::array::from_fn(|side| order.map(|end| borders[side * 2 + end].clone()));
    let rails = [
        fillet.rail_curve(0, &ranges[0][0], &ranges[0][1], dec)?,
        fillet.rail_curve(1, &ranges[1][0], &ranges[1][1], dec)?,
    ];
    let mut uses = BTreeMap::<EdgeId, Vec<FaceId>>::new();
    for &face in faces {
        for lp in &model.face(face).loops {
            if let Loop::Coedges(cs) = model.loop_(*lp) {
                for c in cs {
                    uses.entry(model.coedge(*c).edge).or_default().push(face);
                }
            }
        }
    }
    let mut rows = [None; 4];
    let mut side_faces = [None; 2];
    let mut cuts: [Option<CutCurve>; 2] = [None, None];
    let mut planes = [None, None];
    for c in coedges(model, blend)? {
        let c = model.coedge(c);
        let e = model.edge(c.edge);
        let neighbours = uses[&c.edge]
            .iter()
            .copied()
            .filter(|&f| f != blend)
            .collect::<Vec<_>>();
        let [other] = neighbours.as_slice() else {
            return Err(bad("a blend edge needs exactly one neighbour"));
        };
        let start = model.geom().position(model.vertex(e.start.unwrap()).point);
        let end = model.geom().position(model.vertex(e.end.unwrap()).point);
        let (row, forward) = match model.geom().curve(e.curve) {
            Curve::BSpline(curve) => {
                let side = (0..2)
                    .find(|&side| rails[side].same_definition(curve))
                    .ok_or_else(|| bad("a blend rail is not its complete fitted row"))?;
                side_faces[side] = Some(*other);
                let a = fillet
                    .rail_point(side, &ranges[side][0])
                    .expect("rail node");
                let b = fillet
                    .rail_point(side, &ranges[side][1])
                    .expect("rail node");
                let forward = same(&start, &a, dec)? && same(&end, &b, dec)?;
                if !forward && !(same(&start, &b, dec)? && same(&end, &a, dec)?) {
                    return Err(bad("a rail's vertices do not match its nodes"));
                }
                if forward != (e.sense == Sense::Same) {
                    return Err(bad("a rail reverses its carrier"));
                }
                (side, forward)
            }
            Curve::Intersection(curve) => {
                let surfaces = [model.face(blend).surface, model.face(*other).surface];
                let mut agrees = |a, b| -> Result<bool, Refusal> {
                    let a = model.geom().surface(a);
                    let b = model.geom().surface(b);
                    if a.same_definition(b) {
                        return Ok(true);
                    }
                    let (Surface::Plane(a), Surface::Plane(b)) = (a, b) else {
                        return Ok(false);
                    };
                    for value in [
                        a.normal.cross(&b.normal).norm_squared(),
                        a.normal.dot(&b.origin.minus(&a.origin)),
                    ] {
                        match dec.sign(&value) {
                            Outcome::Proven(_) => return Ok(false),
                            Outcome::ProvenCoincident(proof)
                                if proof.proof != ZeroProof::SeparationBound => {}
                            other => {
                                return Err(bad(
                                    "cap support identity needs exact plane coincidence",
                                )
                                .keeping_undecided(other));
                            }
                        }
                    }
                    Ok(true)
                };
                if !(agrees(curve.surfaces[0], surfaces[0])?
                    && agrees(curve.surfaces[1], surfaces[1])?)
                    && !(agrees(curve.surfaces[1], surfaces[0])?
                        && agrees(curve.surfaces[0], surfaces[1])?)
                {
                    return Err(bad(
                        "a cap edge is not the blend's intersection with its neighbour",
                    ));
                }
                let (normal, offset) = plane_of(model, *other)?;
                let cut = fillet
                    .plane_cut(normal.clone(), offset.clone(), dec)?
                    .ok_or_else(|| bad("the cap has no complete cut"))?;
                let end_index = (0..2)
                    .find(|&i| cut.ends == [ranges[0][i].clone(), ranges[1][i].clone()])
                    .ok_or_else(|| {
                        bad(format!(
                            "cap cut {:?} does not meet the rail borders {ranges:?}",
                            cut.ends
                        ))
                    })?;
                let a = fillet.rail_point(0, &cut.ends[0]).expect("cut end");
                let b = fillet.rail_point(1, &cut.ends[1]).expect("cut end");
                let forward = same(&start, &a, dec)? && same(&end, &b, dec)?;
                if !forward && !(same(&start, &b, dec)? && same(&end, &a, dec)?) {
                    return Err(bad("a cut's vertices do not match its rails"));
                }
                cuts[end_index] =
                    Some(fillet.cut_curve(&cut, &options.chord_tolerance.mul_pow2(-8), dec)?);
                planes[end_index] = Some((normal, offset));
                (2 + end_index, forward)
            }
            _ => return Err(bad("a blend boundary is neither a rail nor a cap cut")),
        };
        let positive = [true, false, false, true][row];
        if rows[row].is_some()
            || ((c.sense == Sense::Same) == forward)
                != (positive == (model.face(blend).sense == Sense::Same))
        {
            return Err(bad("a blend loop duplicates or reverses a chart row"));
        }
        rows[row] = Some((
            c.edge,
            [
                BoundaryRow::Rail,
                BoundaryRow::Mouth,
                BoundaryRow::Negative,
                BoundaryRow::Positive,
            ][row],
            c.sense,
            forward,
        ));
    }
    let [Some(a), Some(b), Some(c), Some(d)] = rows else {
        return Err(bad("the blend is missing a row"));
    };
    let [Some(first), Some(second)] = side_faces else {
        return Err(bad("the blend is missing a support face"));
    };
    if first == second {
        return Err(bad("the blend needs two distinct support faces"));
    }
    let [Some(low), Some(high)] = cuts else {
        return Err(bad("the blend needs two cap cuts"));
    };
    let [Some(low_plane), Some(high_plane)] = planes else {
        return Err(bad("the blend needs two planes"));
    };
    let mapping = Mapping::new(
        fillet.clone(),
        rails.clone(),
        [low, high],
        [low_plane, high_plane],
        dec,
    )?;
    let mut definitions = Vec::new();
    for &face in faces {
        if face == blend {
            definitions.push(Definition {
                face,
                support: Support::Fillet(Box::new(mapping.clone())),
                rows: [a, b, c, d],
            });
        } else if face == first || face == second {
            let side = usize::from(face == second);
            definitions.push(ribbon(model, face, fillet, side, &rails[side], dec)?);
        } else if matches!(
            model.geom().surface(model.face(face).surface),
            Surface::BSpline(_)
        ) {
            let checked = RationalPatchFace::of_chart_face(model, face, dec)?
                .ok_or_else(|| bad("an untrimmed spline neighbour needs a checked patch"))?;
            definitions.push(Definition {
                face,
                support: Support::Patch(checked.support().clone()),
                rows: checked.rows(),
            });
        }
    }
    super::rational_patch::prepare_definitions(model, definitions, options, true, dec)
}

fn polynomial_bilinear(patch: &RationalPatch) -> bool {
    patch.carrier().degrees() == [1, 1]
        && patch.carrier().control_count() == [2, 2]
        && patch
            .carrier()
            .weights()
            .iter()
            .all(|w| w.as_rational() == Some(Rational::one()))
}

fn ribbon(
    model: &Model,
    face: FaceId,
    fillet: &SplineFillet,
    side: usize,
    rail: &BSplineCurve,
    dec: &mut Decider,
) -> Result<Definition, Refusal> {
    let patch = patch_of(model, face).ok_or_else(|| bad("a trimmed side is not a single patch"))?;
    let extended = &fillet.supports()[side];
    if !polynomial_bilinear(patch)
        || !polynomial_bilinear(extended)
        || rail
            .weights()
            .any(|w| w.as_rational() != Some(Rational::one()))
    {
        return Err(bad(
            "trimmed sides currently require polynomial bilinear supports and unit-weight rails",
        ));
    }
    // Four exact corner identities prove equality of the two bilinear maps.
    for u in &patch.domain()[0] {
        for v in &patch.domain()[1] {
            let at = [Real::from(u.clone()), Real::from(v.clone())];
            if !same(
                &patch.point_at(&at, dec)?,
                &extended.point_at(&at, dec)?,
                dec,
            )? {
                return Err(bad(
                    "a fitted support differs from its neighbour's bilinear map",
                ));
            }
        }
    }
    let range = rail.domain();
    let spans: Vec<_> = fillet
        .rail(side)
        .iter()
        .filter(|s| &s.tau[0] >= range[0] && &s.tau[1] <= range[1])
        .collect();
    if spans.is_empty()
        || &spans[0].tau[0] != range[0]
        || &spans.last().unwrap().tau[1] != range[1]
        || spans
            .windows(2)
            .any(|pair| pair[0].tau[1] != pair[1].tau[0])
    {
        return Err(bad(
            "a rail's parameter cubics do not cover its entire domain",
        ));
    }
    let start = &spans[0].parameters;
    let end = &spans.last().unwrap().parameters;
    let along = (0..2)
        .find(|&axis| {
            let domain = &patch.domain()[axis];
            (start[axis][0] == domain[0] && end[axis][3] == domain[1])
                || (start[axis][0] == domain[1] && end[axis][3] == domain[0])
        })
        .ok_or_else(|| bad("a rail does not span a support axis"))?;
    let across = 1 - along;
    let increasing = start[along][0] < end[along][3];
    for span in &spans {
        if span.parameters[along].windows(2).any(|p| {
            if increasing {
                p[0] >= p[1]
            } else {
                p[0] <= p[1]
            }
        }) || span.parameters[across]
            .iter()
            .any(|p| p <= &patch.domain()[across][0] || p >= &patch.domain()[across][1])
        {
            return Err(bad(
                "a rail needs a strictly monotone graph strictly inside its side",
            ));
        }
    }
    let cs = coedges(model, face)?;
    let vertices: Vec<_> = cs
        .iter()
        .flat_map(|c| {
            let e = model.edge(model.coedge(*c).edge);
            [e.start.unwrap(), e.end.unwrap()]
        })
        .collect();
    let mut opposite = None;
    for bound in &patch.domain()[across] {
        let mut both = true;
        for value in &patch.domain()[along] {
            let mut at = [Real::zero(), Real::zero()];
            at[across] = Real::from(bound.clone());
            at[along] = Real::from(value.clone());
            let corner = patch.point_at(&at, dec)?;
            let mut found = false;
            for &v in &vertices {
                if same(&corner, &model.geom().position(model.vertex(v).point), dec)? {
                    found = true;
                    break;
                }
            }
            both &= found;
        }
        if both {
            if opposite.is_some() {
                return Err(bad("a trimmed side has two opposite borders"));
            }
            opposite = Some(bound.clone());
        }
    }
    let opposite = opposite.ok_or_else(|| bad("a trimmed side has no complete opposite border"))?;
    let upper = opposite == patch.domain()[across][1];
    let same_orientation = if along == 0 {
        increasing == upper
    } else {
        increasing != upper
    };
    let degree = rail.degree();
    if degree < 3 || rail.control_count() != spans.len() * degree + 1 {
        return Err(bad("a rail does not retain its Bezier span controls"));
    }
    let mut outer = Vec::new();
    for (i, span) in spans.iter().enumerate() {
        let mut parameters = span.parameters[along].to_vec();
        while parameters.len() <= degree {
            dec.charge_work(parameters.len() as u64 * 4)?;
            let n = parameters.len();
            let mut elevated = vec![parameters[0].clone()];
            for j in 1..n {
                let fraction = Rational::from_ratio(j as i64, n as i64).unwrap();
                elevated
                    .push(&parameters[j] + &(&fraction * &(&parameters[j - 1] - &parameters[j])));
            }
            elevated.push(parameters[n - 1].clone());
            parameters = elevated;
        }
        for parameter in parameters.into_iter().skip(usize::from(i != 0)) {
            let mut at = [Real::zero(), Real::zero()];
            at[across] = Real::from(opposite.clone());
            at[along] = Real::from(parameter);
            outer.push(patch.point_at(&at, dec)?);
        }
    }
    let controls = rail
        .controls()
        .zip(outer)
        .flat_map(|(a, b)| [a.clone(), b])
        .collect::<Vec<_>>();
    let carrier = BSplineSurface::new(
        [degree, 1],
        [
            rail.knots().to_vec(),
            vec![
                Rational::zero(),
                Rational::zero(),
                Rational::one(),
                Rational::one(),
            ],
        ],
        [rail.control_count(), 2],
        controls.clone(),
        vec![Real::one(); controls.len()],
        dec,
    )?;
    let rows = checked_rows(model, face, &carrier, same_orientation, dec)?;
    Ok(Definition {
        face,
        support: Support::Measured(carrier),
        rows,
    })
}

fn checked_rows(
    model: &Model,
    face: FaceId,
    patch: &BSplineSurface,
    same_orientation: bool,
    dec: &mut Decider,
) -> Result<[(EdgeId, BoundaryRow, Sense, bool); 4], Refusal> {
    let sides = [
        SplineBoundary::VMin,
        SplineBoundary::VMax,
        SplineBoundary::UMin,
        SplineBoundary::UMax,
    ];
    let mut rows = [None; 4];
    for id in coedges(model, face)? {
        let c = model.coedge(id);
        let e = model.edge(c.edge);
        let start = model.geom().position(model.vertex(e.start.unwrap()).point);
        let end = model.geom().position(model.vertex(e.end.unwrap()).point);
        let mut matched = false;
        for (i, side) in sides.iter().enumerate() {
            let curve = BSplineCurve::boundary(patch, *side);
            let domain = curve.domain();
            let a = curve.point_at(&Real::from(domain[0].clone()), dec)?;
            let b = curve.point_at(&Real::from(domain[1].clone()), dec)?;
            let forward = if same(&start, &a, dec)? && same(&end, &b, dec)? {
                true
            } else if same(&start, &b, dec)? && same(&end, &a, dec)? {
                false
            } else {
                continue;
            };
            let direction = patch
                .boundary_orientation(*side, model.geom().curve(e.curve), dec)?
                .ok_or_else(|| bad("a ribbon row has no whole-curve identity"))?;
            let natural_same = (model.face(face).sense == Sense::Same) == same_orientation;
            if rows[i].is_some()
                || ((e.sense == Sense::Same) != (direction == forward))
                || (((c.sense == Sense::Same) == forward)
                    != ([true, false, false, true][i] == natural_same))
            {
                return Err(bad(
                    "a ribbon row duplicates or reverses its finite boundary",
                ));
            }
            rows[i] = Some((
                c.edge,
                [
                    BoundaryRow::Rail,
                    BoundaryRow::Mouth,
                    BoundaryRow::Negative,
                    BoundaryRow::Positive,
                ][i],
                c.sense,
                forward,
            ));
            matched = true;
            break;
        }
        if !matched {
            return Err(bad("a trimmed side has an additional or unmatched trim"));
        }
    }
    let [Some(a), Some(b), Some(c), Some(d)] = rows else {
        return Err(bad("a ribbon is missing a row"));
    };
    Ok([a, b, c, d])
}
