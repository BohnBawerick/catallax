//! Charged whole-face certificates for planar source repair.

use crate::exact::is_zero;
use crate::gap::decimal;
use crate::planar_arithmetic as a;
use crate::{
    Curve, Face, FaceTolerance, ImportRefusal, PlanarRepairMethod, Solid, SourceIncidence, Surface,
    Vec3,
};
use catallax_numbers::{Decider, Integer, Rational, Sign};
use std::collections::{BTreeMap, BTreeSet};

type Proof<T> = Result<T, ImportRefusal>;

fn polygon<'a>(face: &'a Face, d: &mut Decider) -> Proof<Option<Vec<&'a Vec3>>> {
    a::checkpoint(d)?;
    let Surface::Plane(plane) = &face.surface else {
        return Ok(None);
    };
    let [bound] = face.loops.as_slice() else {
        return Ok(None);
    };
    if bound.edges.len() < 3 || is_zero(&plane.axis) {
        return Ok(None);
    }
    let mut points = Vec::with_capacity(bound.edges.len());
    for (i, edge) in bound.edges.iter().enumerate() {
        a::checkpoint(d)?;
        let Curve::Line(origin, direction) = &edge.curve else {
            return Ok(None);
        };
        let (from, to) = if edge.forward {
            (&edge.start.point, &edge.end.point)
        } else {
            (&edge.end.point, &edge.start.point)
        };
        let next = &bound.edges[(i + 1) % bound.edges.len()];
        let next_start = if next.forward {
            &next.start.point
        } else {
            &next.end.point
        };
        if !a::same(to, next_start, d)?
            || a::same(from, to, d)?
            || is_zero(direction)
            || !is_zero(&a::cross(&a::sub(from, origin, d)?, direction, d)?)
            || !is_zero(&a::cross(&a::sub(to, origin, d)?, direction, d)?)
            || !a::dot(&a::sub(from, &plane.origin, d)?, &plane.axis, d)?.is_zero()
        {
            return Ok(None);
        }
        points.push(from);
    }
    let mut side = None;
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let edge = a::sub(points[j], points[i], d)?;
        for k in 0..points.len() {
            a::checkpoint(d)?;
            if k == i || k == j {
                continue;
            }
            let sign = a::dot(
                &a::cross(&edge, &a::sub(points[k], points[i], d)?, d)?,
                &plane.axis,
                d,
            )?
            .sign();
            if sign == Sign::Zero || side.is_some_and(|s| s != sign) {
                return Ok(None);
            }
            side = Some(sign);
        }
    }
    Ok(Some(points))
}

fn feature_size2(face: &Face, d: &mut Decider) -> Proof<Option<Rational>> {
    let Some(points) = polygon(face, d)? else {
        return Ok(None);
    };
    let mut size: Option<Rational> = None;
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let edge = a::sub(points[j], points[i], d)?;
        let length2 = a::dot(&edge, &edge, d)?;
        size = Some(if let Some(s) = size {
            a::minimum(s, length2.clone(), d)?
        } else {
            length2.clone()
        });
        for k in 0..points.len() {
            a::checkpoint(d)?;
            if k == i || k == j {
                continue;
            }
            let area = a::cross(&edge, &a::sub(points[k], points[i], d)?, d)?;
            let altitude2 = a::arithmetic(&a::dot(&area, &area, d)?, &length2, '/', d)?;
            size = Some(a::minimum(size.unwrap(), altitude2, d)?);
        }
    }
    Ok(size)
}

fn distance2(point: &Vec3, polygon: &[&Vec3], d: &mut Decider) -> Proof<Rational> {
    let normal = a::cross(
        &a::sub(polygon[1], polygon[0], d)?,
        &a::sub(polygon[2], polygon[0], d)?,
        d,
    )?;
    let normal2 = a::dot(&normal, &normal, d)?;
    let mut inside = true;
    for i in 0..polygon.len() {
        a::checkpoint(d)?;
        let edge = a::sub(polygon[(i + 1) % polygon.len()], polygon[i], d)?;
        let side = a::dot(
            &a::cross(&edge, &a::sub(point, polygon[i], d)?, d)?,
            &normal,
            d,
        )?
        .sign();
        if side == Sign::Negative {
            inside = false;
            break;
        }
    }
    if inside {
        let height = a::dot(&a::sub(point, polygon[0], d)?, &normal, d)?;
        return a::arithmetic(&a::arithmetic(&height, &height, '*', d)?, &normal2, '/', d);
    }
    let mut best = None;
    for i in 0..polygon.len() {
        a::checkpoint(d)?;
        let edge = a::sub(polygon[(i + 1) % polygon.len()], polygon[i], d)?;
        let delta = a::sub(point, polygon[i], d)?;
        let t = a::arithmetic(
            &a::dot(&delta, &edge, d)?,
            &a::dot(&edge, &edge, d)?,
            '/',
            d,
        )?;
        let t = a::minimum(a::maximum(t, Rational::zero(), d)?, Rational::one(), d)?;
        let residual = a::offset(&delta, &edge, &t, '-', d)?;
        let value = a::dot(&residual, &residual, d)?;
        best = Some(if let Some(previous) = best {
            a::minimum(previous, value, d)?
        } else {
            value
        });
    }
    Ok(best.expect("covered polygon"))
}

fn projected(face: &Face, d: &mut Decider) -> Proof<Option<Face>> {
    let Surface::Plane(plane) = &face.surface else {
        return Ok(None);
    };
    let normal2 = a::dot(&plane.axis, &plane.axis, d)?;
    if normal2.is_zero() {
        return Ok(None);
    }
    let mut result = face.clone();
    for edge in result.loops.iter_mut().flat_map(|b| &mut b.edges) {
        a::checkpoint(d)?;
        let Curve::Line(origin, direction) = &mut edge.curve else {
            return Ok(None);
        };
        for point in [&mut edge.start.point, &mut edge.end.point] {
            let t = a::arithmetic(
                &a::dot(&plane.axis, &a::sub(point, &plane.origin, d)?, d)?,
                &normal2,
                '/',
                d,
            )?;
            *point = a::offset(point, &plane.axis, &t, '-', d)?;
        }
        let delta = a::sub(&edge.end.point, &edge.start.point, d)?;
        if is_zero(&delta) {
            return Ok(None);
        }
        *origin = edge.start.point.clone();
        *direction = if edge.curve_same_sense {
            delta
        } else {
            delta.map(|v| -v)
        };
    }
    Ok(Some(result))
}

fn source_correspondence(
    before: &Face,
    defined: &Face,
    after: &Face,
    d: &mut Decider,
) -> Proof<Option<(Rational, Rational, Rational)>> {
    if before.name != after.name
        || before.same_sense != after.same_sense
        || before.loops.len() != after.loops.len()
    {
        return Ok(None);
    }
    let (mut discrepancy2, mut movement2, mut total2) =
        (Rational::zero(), Rational::zero(), Rational::zero());
    for ((x, source), y) in before.loops.iter().zip(&defined.loops).zip(&after.loops) {
        if x.edges.len() != y.edges.len() || x.outer != y.outer || x.orientation != y.orientation {
            return Ok(None);
        }
        for ((x, source), y) in x.edges.iter().zip(&source.edges).zip(&y.edges) {
            a::checkpoint(d)?;
            if x.edge != y.edge
                || x.start.name != y.start.name
                || x.end.name != y.end.name
                || x.forward != y.forward
                || x.curve_same_sense != y.curve_same_sense
            {
                return Ok(None);
            }
            let Curve::Line(origin, direction) = &x.curve else {
                return Ok(None);
            };
            let length2 = a::dot(direction, direction, d)?;
            if length2.is_zero() {
                return Ok(None);
            }
            for ((v, q), target) in [&x.start.point, &x.end.point]
                .into_iter()
                .zip([&source.start.point, &source.end.point])
                .zip([&y.start.point, &y.end.point])
            {
                let t = a::arithmetic(
                    &a::dot(&a::sub(v, origin, d)?, direction, d)?,
                    &length2,
                    '/',
                    d,
                )?;
                let on_line = a::offset(origin, direction, &t, '+', d)?;
                discrepancy2 = a::maximum(discrepancy2, a::squared(v, q, d)?, d)?;
                discrepancy2 = a::maximum(discrepancy2, a::squared(&on_line, q, d)?, d)?;
                let vertex2 = a::squared(v, target, d)?;
                movement2 = a::maximum(movement2, vertex2.clone(), d)?;
                total2 = a::maximum(total2, vertex2, d)?;
                total2 = a::maximum(total2, a::squared(&on_line, target, d)?, d)?;
            }
        }
    }
    Ok(Some((discrepancy2, movement2, total2)))
}

// Canonical rational storage identifies source vertices without a numeric sort.
fn key(point: &Vec3) -> [(Integer, Integer); 3] {
    std::array::from_fn(|i| (point[i].numerator().clone(), point[i].denominator().clone()))
}

pub(crate) fn certify_incidence(
    source: &Solid,
    repaired: &Solid,
    method: PlanarRepairMethod,
) -> Proof<(Vec<FaceTolerance>, Vec<SourceIncidence>)> {
    if source.faces.len() != repaired.faces.len() {
        return Err(ImportRefusal::NotSolid {
            detail: "the planar repair changed the source face count".into(),
        });
    }
    let mut d = Decider::default();
    let mut uses: BTreeMap<_, BTreeSet<usize>> = BTreeMap::new();
    for (index, face) in source.faces.iter().enumerate() {
        for point in face
            .loops
            .iter()
            .flat_map(|b| &b.edges)
            .flat_map(|e| [&e.start.point, &e.end.point])
        {
            a::checkpoint(&mut d)?;
            uses.entry(key(point)).or_default().insert(index);
        }
    }
    let (mut tolerances, mut records) = (Vec::new(), Vec::new());
    let cap = crate::gap::auto_cap();
    let cap2 = a::arithmetic(&cap, &cap, '*', &mut d)?;
    for (before, after) in source.faces.iter().zip(&repaired.faces) {
        a::checkpoint(&mut d)?;
        if before == after {
            continue;
        }
        let fail = |stage: &str, measurement: &str, value: Option<&Rational>, reason: String| {
            crate::planar_incidence::refused(source, before, stage, measurement, value, reason)
        };
        let defined = projected(before, &mut d)?.ok_or_else(|| {
            fail(
                "source_cover",
                "unavailable",
                None,
                "the source needs a straight planar boundary interpretation".into(),
            )
        })?;
        let (discrepancy2, movement2, correspondence2) =
            source_correspondence(before, &defined, after, &mut d)?.ok_or_else(|| {
                fail(
                    "source_curve",
                    "unavailable",
                    None,
                    "the source boundary correspondence could not be bounded".into(),
                )
            })?;
        let x = polygon(&defined, &mut d)?.ok_or_else(|| {
            fail(
                "source_cover",
                "source_boundary_correspondence",
                Some(&discrepancy2),
                "the projected source needs a complete convex single-bound cover".into(),
            )
        })?;
        let y = polygon(after, &mut d)?.ok_or_else(|| {
            fail(
                "result_cover",
                "source_boundary_correspondence",
                Some(&discrepancy2),
                "the repaired face needs a complete convex single-bound cover".into(),
            )
        })?;
        if x.len() != y.len() {
            return Err(fail(
                "face_correspondence",
                "unavailable",
                None,
                "the repair changed the face correspondence".into(),
            ));
        }
        let mut geometry2 = Rational::zero();
        for point in &x {
            geometry2 = a::maximum(geometry2, distance2(point, &y, &mut d)?, &mut d)?;
        }
        for point in &y {
            geometry2 = a::maximum(geometry2, distance2(point, &x, &mut d)?, &mut d)?;
        }
        // Complete convex covers bound geometric deviation in both directions.
        // Corresponding affine maps also cover the raw fan and finite line pieces.
        let total2 = a::maximum(correspondence2, geometry2.clone(), &mut d)?;
        if a::compare(&total2, &cap2, &mut d)?.is_gt() {
            return Err(fail(
                "import_cap",
                "whole_face_source_bound",
                Some(&total2),
                "the proposed source bound exceeds 0.01 mm".into(),
            ));
        }
        let mut feature2 = a::minimum(
            feature_size2(&defined, &mut d)?.expect("checked cover"),
            feature_size2(after, &mut d)?.expect("checked cover"),
            &mut d,
        )?;
        let mut neighbours = BTreeSet::new();
        for point in before
            .loops
            .iter()
            .flat_map(|b| &b.edges)
            .flat_map(|e| [&e.start.point, &e.end.point])
        {
            a::checkpoint(&mut d)?;
            neighbours.extend(uses[&key(point)].iter().copied());
        }
        for index in neighbours {
            a::checkpoint(&mut d)?;
            let neighbour = &source.faces[index];
            if neighbour.name == before.name {
                continue;
            }
            let projected_neighbour = projected(neighbour, &mut d)?.ok_or_else(|| {
                fail(
                    "nearest_feature",
                    "whole_face_source_bound",
                    Some(&total2),
                    format!(
                        "source neighbour {} needs a feature-size certificate",
                        neighbour.name
                    ),
                )
            })?;
            for face in [&projected_neighbour, &repaired.faces[index]] {
                let bound = feature_size2(face, &mut d)?.ok_or_else(|| {
                    fail(
                        "nearest_feature",
                        "whole_face_source_bound",
                        Some(&total2),
                        format!(
                            "neighbour {} needs a complete convex feature bound",
                            face.name
                        ),
                    )
                })?;
                feature2 = a::minimum(feature2, bound, &mut d)?;
            }
        }
        let (_, total_upper) = a::sqrt_enclosure(&total2, &mut d)?;
        let size = a::arithmetic(
            &a::arithmetic(&total_upper, &total_upper, '*', &mut d)?,
            &Rational::from_i64(100),
            '*',
            &mut d,
        )?;
        if a::compare(&size, &feature2, &mut d)?.is_gt() {
            let mut error = fail(
                "nearest_feature",
                "whole_face_source_bound",
                Some(&total2),
                "the certified source bound exceeds a tenth of the nearest feature".into(),
            );
            if let ImportRefusal::PlanarRepairNotVerified { detail } = &mut error {
                detail.feature_size_squared_mm2 = Some(feature2.to_string());
            }
            return Err(error);
        }
        let (geometry_lower, geometry_upper) = a::sqrt_enclosure(&geometry2, &mut d)?;
        let (discrepancy_lower, discrepancy_upper) = a::sqrt_enclosure(&discrepancy2, &mut d)?;
        let (_, movement_upper) = a::sqrt_enclosure(&movement2, &mut d)?;
        tolerances.push(FaceTolerance {
            solid: source.name.clone(),
            face: before.name.clone(),
            deviation_lower_mm: decimal(&geometry_lower),
            deviation_upper_mm: decimal(&total_upper),
            approved_tolerance_mm: decimal(&cap),
            feature_size_squared_mm2: feature2.to_string(),
        });
        records.push(SourceIncidence {
            solid: source.name.clone(),
            face: before.name.clone(),
            method,
            source_definition: "projected_vertex_boundary/1".into(),
            source_discrepancy_lower_mm: decimal(&discrepancy_lower),
            source_discrepancy_upper_mm: decimal(&discrepancy_upper),
            geometric_deviation_lower_mm: decimal(&geometry_lower),
            geometric_deviation_upper_mm: decimal(&geometry_upper),
            source_to_result_upper_mm: decimal(&total_upper),
            vertex_movement_upper_mm: decimal(&movement_upper),
        });
    }
    Ok((tolerances, records))
}
