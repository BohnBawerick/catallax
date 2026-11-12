//! Exact planar derivatives with bounds against the immutable source carriers.

use std::collections::BTreeMap;

use catallax_numbers::{BoundKind, Decider, Mark, MarkLog, Rational};
use catallax_shapes::spline_tools::RationalSplineCurve;
use serde::Serialize;

use crate::exact::is_zero;
use crate::incidence_arithmetic::{
    arithmetic, checkpoint, compare, cross, dot, equal, maximum, minimum, same, same_sequence,
    sqrt_enclosure, sub,
};
use crate::{Curve, Placement, Solid, Surface, Vec3};

#[derive(Serialize)]
struct LocalTolerance {
    solid: String,
    entity: String,
    entity_kind: &'static str,
    source_geometry_sha256: String,
    definition: &'static str,
    deviation_lower_mm: String,
    deviation_upper_mm: String,
    approved_tolerance_mm: String,
    feature_size_squared_mm2: String,
}

pub(crate) struct Prepared {
    pub solid: Solid,
    pub marks: MarkLog,
}

fn length2(v: &Vec3, dec: &mut Decider) -> Result<Rational, String> {
    dot(v, v, dec)
}
fn distance2(a: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<Rational, String> {
    let difference = sub(a, b, dec)?;
    length2(&difference, dec)
}

pub(crate) fn project(
    point: &Vec3,
    planes: &[Placement],
    dec: &mut Decider,
) -> Result<Vec3, String> {
    project_inner(point, planes, dec)
        .map_err(|e| format!("planar projection/RREF certificate: {e}"))
}

// RREF rejects inconsistent dependent planes before the minimum-norm solve.
fn project_inner(point: &Vec3, planes: &[Placement], dec: &mut Decider) -> Result<Vec3, String> {
    checkpoint(dec)?;
    let mut rows = Vec::new();
    for plane in planes {
        checkpoint(dec)?;
        rows.push(vec![
            plane.axis[0].clone(),
            plane.axis[1].clone(),
            plane.axis[2].clone(),
            dot(&plane.axis, &plane.origin, dec)?,
        ]);
    }
    let mut rank = 0;
    for column in 0..3 {
        checkpoint(dec)?;
        let mut pivot = None;
        for (r, row) in rows.iter().enumerate().skip(rank) {
            checkpoint(dec)?;
            if !row[column].is_zero() {
                pivot = Some(r);
                break;
            }
        }
        let Some(pivot) = pivot else { continue };
        rows.swap(rank, pivot);
        let divisor = rows[rank][column].clone();
        for x in &mut rows[rank] {
            *x = arithmetic(x, &divisor, '/', dec)?;
        }
        let pivot_row = rows[rank].clone();
        for (row, target) in rows.iter_mut().enumerate() {
            checkpoint(dec)?;
            if row == rank {
                continue;
            }
            let multiplier = target[column].clone();
            for (value, pivot) in target.iter_mut().zip(&pivot_row) {
                let product = arithmetic(&multiplier, pivot, '*', dec)?;
                *value = arithmetic(value, &product, '-', dec)?;
            }
        }
        rank += 1;
    }
    for row in rows.iter().skip(rank) {
        checkpoint(dec)?;
        if !row[3].is_zero() {
            return Err("adjacent source planes have no exact common intersection".into());
        }
    }
    rows.truncate(rank);
    let normals = rows
        .iter()
        .map(|r| [r[0].clone(), r[1].clone(), r[2].clone()])
        .collect::<Vec<Vec3>>();
    let mut system = Vec::new();
    for i in 0..rank {
        checkpoint(dec)?;
        let mut row = Vec::new();
        for normal in &normals {
            row.push(dot(&normals[i], normal, dec)?);
        }
        let height = dot(&normals[i], point, dec)?;
        row.push(arithmetic(&height, &rows[i][3], '-', dec)?);
        system.push(row);
    }
    for i in 0..rank {
        checkpoint(dec)?;
        let mut pivot = None;
        for (r, row) in system.iter().enumerate().skip(i) {
            checkpoint(dec)?;
            if !row[i].is_zero() {
                pivot = Some(r);
                break;
            }
        }
        let pivot = pivot.ok_or("dependent plane projection")?;
        system.swap(i, pivot);
        let divisor = system[i][i].clone();
        for x in &mut system[i] {
            *x = arithmetic(x, &divisor, '/', dec)?;
        }
        let pivot_row = system[i].clone();
        for (row, target) in system.iter_mut().enumerate() {
            checkpoint(dec)?;
            if row == i {
                continue;
            }
            let multiplier = target[i].clone();
            for (value, pivot) in target.iter_mut().zip(&pivot_row) {
                let product = arithmetic(&multiplier, pivot, '*', dec)?;
                *value = arithmetic(value, &product, '-', dec)?;
            }
        }
    }
    let mut result = point.clone();
    for j in 0..3 {
        let mut sum = Rational::zero();
        for i in 0..rank {
            let product = arithmetic(&system[i][rank], &normals[i][j], '*', dec)?;
            sum = arithmetic(&sum, &product, '+', dec)?;
        }
        result[j] = arithmetic(&point[j], &sum, '-', dec)?;
    }
    Ok(result)
}

pub(crate) type Net = Vec<[Rational; 4]>;
pub(crate) fn split(net: &Net, dec: &mut Decider) -> Result<(Net, Net), String> {
    checkpoint(dec)?;
    if net.is_empty() {
        return Err("empty displacement net".into());
    }
    let half = Rational::from_ratio(1, 2).unwrap();
    let mut work = net.clone();
    let mut left = vec![work[0].clone()];
    let mut right = vec![work[work.len() - 1].clone()];
    while work.len() > 1 {
        checkpoint(dec)?;
        let mut next = Vec::with_capacity(work.len() - 1);
        for pair in work.windows(2) {
            checkpoint(dec)?;
            let mut point = std::array::from_fn(|_| Rational::zero());
            for i in 0..4 {
                let sum = arithmetic(&pair[0][i], &pair[1][i], '+', dec)?;
                point[i] = arithmetic(&sum, &half, '*', dec)?;
            }
            next.push(point);
        }
        work = next;
        left.push(work[0].clone());
        right.push(work[work.len() - 1].clone());
    }
    right.reverse();
    Ok((left, right))
}
fn net_norm(p: &[Rational; 4], dec: &mut Decider) -> Result<Rational, String> {
    checkpoint(dec)?;
    if p[3].sign() != catallax_numbers::Sign::Positive {
        return Err("nonpositive displacement net weight".into());
    }
    let v = [
        arithmetic(&p[0], &p[3], '/', dec)?,
        arithmetic(&p[1], &p[3], '/', dec)?,
        arithmetic(&p[2], &p[3], '/', dec)?,
    ];
    length2(&v, dec)
}
pub(crate) fn net_bounds(net: &Net, dec: &mut Decider) -> Result<(Rational, Rational), String> {
    checkpoint(dec)?;
    if net.is_empty() {
        return Err("empty displacement net".into());
    }
    let mut upper = Rational::zero();
    for point in net {
        upper = maximum(upper, net_norm(point, dec)?, dec)?;
    }
    let (left, _) = split(net, dec)?;
    let a = net_norm(&net[0], dec)?;
    let b = net_norm(&net[net.len() - 1], dec)?;
    let middle = net_norm(&left[left.len() - 1], dec)?;
    let lower = maximum(a, b, dec)?;
    let lower = maximum(lower, middle, dec)?;
    Ok((lower, upper))
}
fn spline_movement(
    before: &RationalSplineCurve,
    after: &RationalSplineCurve,
    dec: &mut Decider,
) -> Result<(Rational, Rational), String> {
    movement_inner(before, after, dec)
        .map_err(|e| format!("retained spline displacement certificate: {e}"))
}
fn movement_inner(
    before: &RationalSplineCurve,
    after: &RationalSplineCurve,
    dec: &mut Decider,
) -> Result<(Rational, Rational), String> {
    checkpoint(dec)?;
    if before.degree() != after.degree()
        || !same_sequence(before.knots(), after.knots(), dec)?
        || !same_sequence(before.weights(), after.weights(), dec)?
        || before.controls().len() != after.controls().len()
    {
        return Err("source parameter correspondence changed".into());
    }
    let mut difference = Vec::new();
    let mut unchanged = true;
    for (a, b) in before.controls().iter().zip(after.controls()) {
        checkpoint(dec)?;
        unchanged &= same(a, b, dec)?;
        difference.push(sub(a, b, dec)?);
    }
    if unchanged {
        return Ok((Rational::zero(), Rational::zero()));
    }
    let spline = RationalSplineCurve::new(
        before.degree(),
        before.knots().to_vec(),
        difference,
        before.weights().to_vec(),
        dec,
    )
    .map_err(|e| e.detail)?;
    let curve = spline.to_curve(dec).map_err(|e| e.detail)?;
    let mut nets = Vec::new();
    for (span, _) in curve.spans() {
        checkpoint(dec)?;
        let mut net = Vec::new();
        for point in curve
            .span_bezier_controls(span, dec)
            .map_err(|e| e.detail)?
        {
            checkpoint(dec)?;
            let mut result = std::array::from_fn(|_| Rational::zero());
            for i in 0..4 {
                result[i] = point[i]
                    .as_rational()
                    .ok_or("nonrational displacement control")?;
            }
            net.push(result);
        }
        nets.push(net);
    }
    if nets.is_empty() {
        return Err("empty displacement span set".into());
    }
    let mut lower = Rational::zero();
    for _ in 0..4096 {
        checkpoint(dec)?;
        let mut upper = Rational::zero();
        let mut index = 0;
        for (i, net) in nets.iter().enumerate() {
            checkpoint(dec)?;
            let bounds = net_bounds(net, dec)?;
            lower = maximum(lower, bounds.0, dec)?;
            if compare(&bounds.1, &upper, dec)?.is_ge() {
                index = i;
                upper = bounds.1;
            }
        }
        let lo = sqrt_enclosure(&lower, dec)?.0;
        let hi = sqrt_enclosure(&upper, dec)?.1;
        let width = arithmetic(&hi, &lo, '-', dec)?;
        if equal(&upper, &lower, dec)?
            || compare(
                &width,
                &Rational::from_ratio(1, 1_000_000_000_000).unwrap(),
                dec,
            )?
            .is_le()
        {
            return Ok((lo, hi));
        }
        let (a, b) = split(&nets.swap_remove(index), dec)?;
        nets.push(a);
        nets.push(b);
    }
    let mut upper = Rational::zero();
    for net in &nets {
        upper = maximum(upper, net_bounds(net, dec)?.1, dec)?;
    }
    Ok((
        sqrt_enclosure(&lower, dec)?.0,
        sqrt_enclosure(&upper, dec)?.1,
    ))
}

pub(crate) fn prepare(source: &Solid, measured: &crate::IncidenceGap) -> Result<Prepared, String> {
    if source
        .faces
        .iter()
        .any(|f| !matches!(f.surface, Surface::Plane(_)))
    {
        return crate::incidence_curved::prepare(source);
    }
    let mut dec = Decider::default();
    checkpoint(&mut dec)?;
    let mut vertices: BTreeMap<Vec3, Vec<Placement>> = BTreeMap::new();
    let mut planes: BTreeMap<String, Vec<Placement>> = BTreeMap::new();
    let mut pairs: BTreeMap<(Vec3, Vec3), (String, bool)> = BTreeMap::new();
    for face in &source.faces {
        checkpoint(&mut dec)?;
        let Surface::Plane(p) = &face.surface else {
            return Err(
                "exact incidence repair currently requires planar support on every face".into(),
            );
        };
        for edge in face.loops.iter().flat_map(|b| &b.edges) {
            checkpoint(&mut dec)?;
            if matches!(edge.curve, Curve::Circle(..) | Curve::Ellipse(..)) {
                return Err(
                    "exact planar incidence repair requires line or spline boundaries".into(),
                );
            }
            let (a, b) = if crate::incidence_arithmetic::vector_compare(
                &edge.start.point,
                &edge.end.point,
                &mut dec,
            )?
            .is_le()
            {
                (&edge.start.point, &edge.end.point)
            } else {
                (&edge.end.point, &edge.start.point)
            };
            let exact_line = if let Curve::Line(o, d) = &edge.curve {
                let ao = sub(a, o, &mut dec)?;
                let bo = sub(b, o, &mut dec)?;
                is_zero(&cross(&ao, d, &mut dec)?) && is_zero(&cross(&bo, d, &mut dec)?)
            } else {
                false
            };
            if let Some((name, exact)) = pairs.get(&(a.clone(), b.clone())) {
                if name != &edge.edge && !(*exact && exact_line) {
                    return Err("the planar derivative needs a seam-identity certificate for distinct source edges".into());
                }
            } else {
                pairs.insert((a.clone(), b.clone()), (edge.edge.clone(), exact_line));
            }
            planes.entry(edge.edge.clone()).or_default().push(p.clone());
            for v in [&edge.start, &edge.end] {
                vertices.entry(v.point.clone()).or_default().push(p.clone());
            }
        }
    }
    let cap = crate::gap::auto_cap();
    let cap2 = arithmetic(&cap, &cap, '*', &mut dec)?;
    let mut moved = BTreeMap::new();
    let mut displacements = BTreeMap::new();
    let mut movement = Rational::zero();
    for (point, planes) in vertices {
        checkpoint(&mut dec)?;
        let result = project(&point, &planes, &mut dec)?;
        let d2 = distance2(&point, &result, &mut dec)?;
        if compare(&d2, &cap2, &mut dec)?.is_gt() {
            return Err("an adjacent vertex needs movement above the 10 um import cap".into());
        }
        movement = maximum(movement, sqrt_enclosure(&d2, &mut dec)?.1, &mut dec)?;
        displacements.insert(point.clone(), d2);
        moved.insert(point, result);
    }
    let mut images = BTreeMap::new();
    for (original, image) in &moved {
        checkpoint(&mut dec)?;
        if images.insert(image.clone(), original).is_some() {
            return Err("the planar derivative would merge distinct source vertices; seam identity is unproved".into());
        }
    }
    let mut result = source.clone();
    let mut records = BTreeMap::new();
    for edge in result
        .faces
        .iter_mut()
        .flat_map(|f| &mut f.loops)
        .flat_map(|b| &mut b.edges)
    {
        checkpoint(&mut dec)?;
        let before = edge.clone();
        edge.start.point = moved[&before.start.point].clone();
        edge.end.point = moved[&before.end.point].clone();
        if let Some(curve) = records
            .get(&edge.edge)
            .map(|(_, _, curve): &(Rational, Rational, Curve)| curve.clone())
        {
            edge.curve = curve;
            continue;
        }
        let (lo, hi) = match &before.curve {
            Curve::Line(o, d) => {
                let d2 = dot(d, d, &mut dec)?;
                let mut maximum_distance = Rational::zero();
                for (old, new) in [
                    (&before.start.point, &edge.start.point),
                    (&before.end.point, &edge.end.point),
                ] {
                    checkpoint(&mut dec)?;
                    let delta = sub(old, o, &mut dec)?;
                    let numerator = dot(&delta, d, &mut dec)?;
                    let t = arithmetic(&numerator, &d2, '/', &mut dec)?;
                    let mut source_point = o.clone();
                    for i in 0..3 {
                        let offset = arithmetic(&t, &d[i], '*', &mut dec)?;
                        source_point[i] = arithmetic(&o[i], &offset, '+', &mut dec)?;
                    }
                    let source_gap2 = distance2(old, &source_point, &mut dec)?;
                    if compare(&source_gap2, &cap2, &mut dec)?.is_gt() {
                        let (lo, hi) = sqrt_enclosure(&source_gap2, &mut dec)?;
                        return Err(format!(
                            "source endpoint cap: edge {} has a line distance in [{}, {}] mm above 10 um",
                            edge.edge,
                            crate::gap::decimal(&lo),
                            crate::gap::decimal(&hi)
                        ));
                    }
                    maximum_distance = maximum(
                        maximum_distance,
                        distance2(&source_point, new, &mut dec)?,
                        &mut dec,
                    )?;
                }
                edge.curve = Curve::Line(
                    edge.start.point.clone(),
                    if edge.curve_same_sense {
                        sub(&edge.end.point, &edge.start.point, &mut dec)?
                    } else {
                        sub(&edge.start.point, &edge.end.point, &mut dec)?
                    },
                );
                sqrt_enclosure(&maximum_distance, &mut dec)?
            }
            Curve::BSpline(s) => {
                let original = s
                    .geometry
                    .clone()
                    .into_clamped(&dec)
                    .map_err(|e| e.detail)?;
                let mut controls = Vec::new();
                for point in original.controls() {
                    controls.push(project(point, &planes[&edge.edge], &mut dec)?);
                }
                let last = controls.len() - 1;
                let ends = if edge.curve_same_sense {
                    [&edge.start.point, &edge.end.point]
                } else {
                    [&edge.end.point, &edge.start.point]
                };
                let source_ends = if edge.curve_same_sense {
                    [&before.start.point, &before.end.point]
                } else {
                    [&before.end.point, &before.start.point]
                };
                let measured_cap = measured.edge == edge.edge
                    && measured.face.is_none()
                    && compare(
                        &Rational::parse(&measured.distance_upper_mm).expect("measured decimal"),
                        &cap,
                        &mut dec,
                    )?
                    .is_le();
                for (index, end) in [(0, ends[0]), (last, ends[1])] {
                    checkpoint(&mut dec)?;
                    let source_end = if index == 0 {
                        source_ends[0]
                    } else {
                        source_ends[1]
                    };
                    if compare(
                        &distance2(&original.controls()[index], source_end, &mut dec)?,
                        &cap2,
                        &mut dec,
                    )?
                    .is_gt()
                        && !measured_cap
                    {
                        return Err(format!(
                            "source spline endpoint cap is unproved for edge {}; a nearest interior correspondence is required",
                            edge.edge
                        ));
                    }
                    if compare(
                        &distance2(&original.controls()[index], end, &mut dec)?,
                        &cap2,
                        &mut dec,
                    )?
                    .is_gt()
                    {
                        return Err("the finite spline endpoints need a certified interior trim correspondence".into());
                    }
                    controls[index] = end.clone();
                }
                let geometry = RationalSplineCurve::new(
                    original.degree(),
                    original.knots().to_vec(),
                    controls,
                    original.weights().to_vec(),
                    &dec,
                )
                .map_err(|e| e.detail)?;
                let bound = spline_movement(&original, &geometry, &mut dec)?;
                let mut curve = s.clone();
                curve.geometry = geometry;
                edge.curve = Curve::BSpline(curve);
                bound
            }
            _ => unreachable!(),
        };
        if compare(&hi, &cap, &mut dec)?.is_gt() {
            return Err(format!(
                "edge {} needs a retained-curve deviation bound of {} mm, above the 10 um import cap",
                edge.edge,
                crate::gap::decimal(&hi)
            ));
        }
        movement = maximum(movement, hi.clone(), &mut dec)?;
        records.insert(edge.edge.clone(), (lo, hi, edge.curve.clone()));
    }
    if &result == source {
        return Err("the planar derivative cannot change the failed incidence".into());
    }
    let mut cover = result.clone();
    let mut reduced: BTreeMap<String, Curve> = BTreeMap::new();
    for edge in cover
        .faces
        .iter_mut()
        .flat_map(|f| &mut f.loops)
        .flat_map(|b| &mut b.edges)
    {
        checkpoint(&mut dec)?;
        if let Some(curve) = reduced.get(&edge.edge) {
            edge.curve = curve.clone();
            continue;
        }
        if let Curve::BSpline(s) = &edge.curve {
            let chord = sub(&edge.end.point, &edge.start.point, &mut dec)?;
            for p in s.geometry.controls() {
                checkpoint(&mut dec)?;
                let offset = sub(p, &edge.start.point, &mut dec)?;
                if !is_zero(&cross(&offset, &chord, &mut dec)?) {
                    return Err("nearest-feature proof needs straight repaired boundaries".into());
                }
            }
            if let Some(curve) = crate::incidence_straight::reduce(edge, &mut dec)? {
                edge.curve = curve;
                reduced.insert(edge.edge.clone(), edge.curve.clone());
                continue;
            }
            let catallax_shapes::Curve::Line(line) = crate::spline_build::curve(edge, &mut dec)
                .map_err(|e| format!("straight reduction certificate: {e:?}"))?
            else {
                return Err("the straight repaired boundary has no exact line reduction".into());
            };
            let origin = line
                .origin
                .as_rational()
                .ok_or("nonrational line origin after straight reduction")?;
            let direction = line
                .dir
                .as_rational()
                .ok_or("nonrational line direction after straight reduction")?;
            edge.curve = Curve::Line(
                [origin.x, origin.y, origin.z],
                [direction.x, direction.y, direction.z],
            );
        }
        reduced.insert(edge.edge.clone(), edge.curve.clone());
    }
    let mut feature2: Option<Rational> = None;
    for face in &cover.faces {
        checkpoint(&mut dec)?;
        let value = crate::incidence_features::feature_size2(face)?;
        feature2 = Some(if let Some(previous) = feature2 {
            minimum(previous, value, &mut dec)?
        } else {
            value
        });
    }
    let feature2 = feature2.ok_or("empty planar solid")?;
    let margin = arithmetic(&movement, &Rational::from_i64(2), '*', &mut dec)?;
    let feature = arithmetic(
        &sqrt_enclosure(&feature2, &mut dec)?.0,
        &margin,
        '-',
        &mut dec,
    )?;
    let required = arithmetic(&movement, &Rational::from_i64(10), '*', &mut dec)?;
    if feature.sign() != catallax_numbers::Sign::Positive
        || compare(&required, &feature, &mut dec)?.is_gt()
    {
        return Err("deviation exceeds a tenth of the certified neighbouring feature size".into());
    }
    let feature_squared = arithmetic(&feature, &feature, '*', &mut dec)?.to_string();
    let source_geometry_sha256 = crate::sha256_hex(&crate::canonical_bytes(&crate::StepFile {
        solids: vec![source.clone()],
        mm_per_unit: Rational::one(),
    }));
    let mut marks = MarkLog::new();
    for (edge, (lo, hi, curve)) in records {
        checkpoint(&mut dec)?;
        if hi.is_zero() {
            continue;
        }
        let record = LocalTolerance {
            solid: source.name.clone(),
            entity: edge,
            entity_kind: "edge",
            source_geometry_sha256: source_geometry_sha256.clone(),
            definition: if matches!(curve, Curve::Line(..)) {
                "pointwise displacement of the projected finite source line segment"
            } else {
                "pointwise retained-spline displacement; unchanged knots and weights"
            },
            deviation_lower_mm: crate::gap::decimal(&lo),
            deviation_upper_mm: crate::gap::decimal(&hi),
            approved_tolerance_mm: crate::gap::decimal(&cap),
            feature_size_squared_mm2: feature_squared.clone(),
        };
        marks.push(Mark::new(
            "import.incidence-repair",
            BoundKind::Distance,
            hi,
            serde_json::to_string(&record).expect("serializable edge tolerance"),
        ));
    }
    let mut vertices = BTreeMap::new();
    for v in source
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start, &e.end])
    {
        checkpoint(&mut dec)?;
        vertices.insert(v.name.clone(), displacements[&v.point].clone());
    }
    for (entity, d2) in vertices {
        checkpoint(&mut dec)?;
        if d2.is_zero() {
            continue;
        }
        let (lo, hi) = sqrt_enclosure(&d2, &mut dec)?;
        let record = LocalTolerance {
            solid: source.name.clone(),
            entity,
            entity_kind: "vertex",
            source_geometry_sha256: source_geometry_sha256.clone(),
            definition: "Euclidean source vertex displacement",
            deviation_lower_mm: crate::gap::decimal(&lo),
            deviation_upper_mm: crate::gap::decimal(&hi),
            approved_tolerance_mm: crate::gap::decimal(&cap),
            feature_size_squared_mm2: feature_squared.clone(),
        };
        marks.push(Mark::new(
            "import.incidence-vertex-repair",
            BoundKind::Distance,
            hi,
            serde_json::to_string(&record).expect("serializable vertex tolerance"),
        ));
    }
    checkpoint(&mut dec)?;
    Ok(Prepared {
        solid: cover,
        marks,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::with_operation_budget;

    fn point(a: i64, b: i64, c: i64) -> Vec3 {
        [a.into(), b.into(), c.into()]
    }
    fn plane() -> Placement {
        Placement {
            origin: point(0, 0, 0),
            axis: point(0, 0, 1),
            ref_dir: point(1, 0, 0),
        }
    }
    fn run_project(
        point: &Vec3,
        planes: &[Placement],
        limit: u64,
    ) -> (Result<Vec3, String>, catallax_numbers::OperationWork) {
        with_operation_budget(limit, || project(point, planes, &mut Decider::default()))
    }

    #[test]
    fn projection_charges_before_returning_and_refuses_an_exhausted_allowance() {
        let p = point(1, 2, 3);
        let planes = [plane()];
        let (result, work) = run_project(&p, &planes, 20_000_000);
        assert_eq!(result.unwrap(), point(1, 2, 0));
        assert!(work.used > 0);
        assert!(!work.exhausted);
        for limit in [0, work.used / 2] {
            let (result, limited) = run_project(&p, &planes, limit);
            assert!(result.unwrap_err().contains("numbers.operation-budget"));
            assert!(limited.exhausted);
            assert!(limited.used <= limit);
        }
    }

    #[test]
    fn projection_accounts_for_large_exact_operands() {
        let mut support = plane();
        support.axis = point(1, 1, 1);
        let p = point(1, 2, 3);
        let (a, small) = run_project(&p, &[support.clone()], 20_000_000);
        assert_eq!(a.unwrap(), point(-1, 0, 1));
        let shift = Rational::parse(&format!("0.{}1", "0".repeat(300))).unwrap();
        support.origin = [shift.clone(), shift.clone(), shift.clone()];
        let translated = p.map(|x| x + shift.clone());
        let (b, large) = run_project(&translated, &[support.clone()], 20_000_000);
        assert_eq!(b.unwrap(), point(-1, 0, 1).map(|x| x + shift.clone()));
        assert!(!large.exhausted);
        assert!(large.used > small.used);
        let (refused, limited) = run_project(&translated, &[support], small.used);
        assert!(refused.is_err());
        assert!(limited.exhausted);
    }

    #[test]
    fn projection_preserves_dependent_plane_checks() {
        let mut x = plane();
        x.axis = point(1, 0, 0);
        assert_eq!(
            run_project(&point(1, 2, 3), &[plane(), x], 20_000_000)
                .0
                .unwrap(),
            point(0, 2, 0)
        );
        let mut offset = plane();
        offset.origin = point(0, 0, 1);
        let (result, work) = run_project(&point(1, 2, 3), &[plane(), offset], 20_000_000);
        assert!(result.unwrap_err().contains("no exact common intersection"));
        assert!(!work.exhausted);
    }

    #[test]
    fn displacement_net_bounds_charge_subdivision_and_refuse_on_budget() {
        let net = vec![[1.into(), 0.into(), 0.into(), 1.into()]; 3];
        let (result, work) =
            with_operation_budget(20_000_000, || net_bounds(&net, &mut Decider::default()));
        assert_eq!(result.unwrap(), (Rational::one(), Rational::one()));
        assert!(work.used > 0);
        assert!(!work.exhausted);
        let (result, limited) =
            with_operation_budget(work.used / 2, || net_bounds(&net, &mut Decider::default()));
        assert!(result.unwrap_err().contains("numbers.operation-budget"));
        assert!(limited.exhausted);
        let (result, zero) = with_operation_budget(0, || split(&net, &mut Decider::default()));
        assert!(result.is_err());
        assert!(zero.exhausted);
    }

    fn movement(
        bow: Rational,
        limit: u64,
    ) -> (
        Result<(Rational, Rational), String>,
        catallax_numbers::OperationWork,
    ) {
        let knots = vec![0.into(), 0.into(), 0.into(), 1.into(), 1.into(), 1.into()];
        let before = RationalSplineCurve::new(
            2,
            knots.clone(),
            vec![
                point(0, 0, 0),
                [Rational::from_ratio(1, 2).unwrap(), 0.into(), 0.into()],
                point(1, 0, 0),
            ],
            vec![1.into(); 3],
            &Decider::default(),
        )
        .unwrap();
        let mut controls = before.controls().to_vec();
        controls[1][1] = bow;
        let after =
            RationalSplineCurve::new(2, knots, controls, vec![1.into(); 3], &Decider::default())
                .unwrap();
        with_operation_budget(limit, || {
            spline_movement(&before, &after, &mut Decider::default())
        })
    }

    #[test]
    fn spline_movement_charges_work_and_refuses_when_the_allowance_is_spent() {
        let bow = Rational::parse("0.001").unwrap();
        let (result, work) = movement(bow.clone(), 20_000_000);
        assert_eq!(
            result.unwrap(),
            (
                Rational::parse("0.0005").unwrap(),
                Rational::parse("0.0005").unwrap()
            )
        );
        assert!(work.used > 0);
        assert!(!work.exhausted);
        for limit in [0, work.used / 2] {
            let (refused, limited) = movement(bow.clone(), limit);
            assert!(refused.unwrap_err().contains("numbers.operation-budget"));
            assert!(limited.exhausted);
            assert!(limited.used <= limit);
        }
    }

    #[test]
    fn spline_movement_accounts_for_large_subdivision_operands() {
        let bow = Rational::parse("0.001").unwrap();
        let (_, small) = movement(bow.clone(), 20_000_000);
        let extra = Rational::parse(&format!("0.{}1", "0".repeat(300))).unwrap();
        let large_bow = bow + extra;
        let expected = large_bow.mul_pow2(-1);
        let (result, large) = movement(large_bow.clone(), 20_000_000);
        let (lo, hi) = result.unwrap();
        assert!(lo <= expected && hi >= expected);
        assert!(large.used > small.used);
        assert!(!large.exhausted);
        let (result, limited) = movement(large_bow, small.used);
        assert!(result.is_err());
        assert!(limited.exhausted);
    }

    #[test]
    fn displacement_subdivision_accounts_for_operand_size_without_curve_conversion() {
        let small_net = vec![[1.into(), 0.into(), 0.into(), 1.into()]; 3];
        let value = Rational::one() + Rational::parse(&format!("0.{}1", "0".repeat(300))).unwrap();
        let large_net = vec![[value.clone(), 0.into(), 0.into(), 1.into()]; 3];
        let (small, small_work) =
            with_operation_budget(20_000_000, || split(&small_net, &mut Decider::default()));
        let (large, large_work) =
            with_operation_budget(20_000_000, || split(&large_net, &mut Decider::default()));
        assert_eq!(small.unwrap(), (small_net.clone(), small_net.clone()));
        assert_eq!(large.unwrap(), (large_net.clone(), large_net.clone()));
        assert!(large_work.used > small_work.used);
        let (refused, limited) = with_operation_budget(small_work.used, || {
            split(&large_net, &mut Decider::default())
        });
        assert!(refused.is_err());
        assert!(limited.exhausted);
        let (_, small_work) = with_operation_budget(20_000_000, || {
            net_bounds(&small_net, &mut Decider::default())
        });
        let (large, large_work) = with_operation_budget(20_000_000, || {
            net_bounds(&large_net, &mut Decider::default())
        });
        assert_eq!(large.unwrap(), (&value * &value, &value * &value));
        assert!(large_work.used > small_work.used);
        let (refused, limited) = with_operation_budget(small_work.used, || {
            net_bounds(&large_net, &mut Decider::default())
        });
        assert!(refused.is_err());
        assert!(limited.exhausted);
    }
}
