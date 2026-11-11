//! Exact planar derivatives with bounds against the immutable source carriers.

use std::collections::BTreeMap;

use catallax_numbers::{BoundKind, Decider, Mark, MarkLog, Rational};
use catallax_shapes::spline_tools::RationalSplineCurve;
use serde::Serialize;

use crate::exact::{cross, dot, is_zero, sub};
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

fn length2(v: &Vec3) -> Rational {
    dot(v, v)
}

fn available() -> Result<(), String> {
    if catallax_numbers::current_operation_work().exhausted {
        Err("the operation work allowance is exhausted".into())
    } else {
        Ok(())
    }
}

// RREF rejects inconsistent dependent planes before the minimum-norm solve.
fn project(point: &Vec3, planes: &[Placement]) -> Result<Vec3, String> {
    available()?;
    let mut rows: Vec<Vec<Rational>> = planes
        .iter()
        .map(|p| {
            vec![
                p.axis[0].clone(),
                p.axis[1].clone(),
                p.axis[2].clone(),
                dot(&p.axis, &p.origin),
            ]
        })
        .collect();
    let mut rank = 0;
    for column in 0..3 {
        available()?;
        let Some(pivot) = (rank..rows.len()).find(|&r| !rows[r][column].is_zero()) else {
            continue;
        };
        rows.swap(rank, pivot);
        let divisor = rows[rank][column].clone();
        for x in &mut rows[rank] {
            *x = &*x / &divisor;
        }
        let pivot_row = rows[rank].clone();
        for (row, target) in rows.iter_mut().enumerate() {
            available()?;
            if row == rank {
                continue;
            }
            let multiplier = target[column].clone();
            for (value, pivot) in target.iter_mut().zip(&pivot_row) {
                *value = &*value - &(&multiplier * pivot);
            }
        }
        rank += 1;
    }
    if rows[rank..].iter().any(|r| !r[3].is_zero()) {
        return Err("adjacent source planes have no exact common intersection".into());
    }
    rows.truncate(rank);
    let normals: Vec<Vec3> = rows
        .iter()
        .map(|r| [r[0].clone(), r[1].clone(), r[2].clone()])
        .collect();
    let mut system: Vec<Vec<Rational>> = (0..rank)
        .map(|i| {
            let mut row = (0..rank)
                .map(|j| dot(&normals[i], &normals[j]))
                .collect::<Vec<_>>();
            row.push(&dot(&normals[i], point) - &rows[i][3]);
            row
        })
        .collect();
    for i in 0..rank {
        available()?;
        let pivot = (i..rank)
            .find(|&r| !system[r][i].is_zero())
            .ok_or("dependent plane projection")?;
        system.swap(i, pivot);
        let divisor = system[i][i].clone();
        for x in &mut system[i] {
            *x = &*x / &divisor;
        }
        let pivot_row = system[i].clone();
        for (row, target) in system.iter_mut().enumerate() {
            if row == i {
                continue;
            }
            let multiplier = target[i].clone();
            for (value, pivot) in target.iter_mut().zip(&pivot_row) {
                *value = &*value - &(&multiplier * pivot);
            }
        }
    }
    Ok(std::array::from_fn(|j| {
        point[j].clone()
            - (0..rank).fold(Rational::zero(), |sum, i| {
                sum + &system[i][rank] * &normals[i][j]
            })
    }))
}

type Net = Vec<[Rational; 4]>;
fn split(net: &Net) -> (Net, Net) {
    let mut work = net.clone();
    let mut left = vec![work[0].clone()];
    let mut right = vec![work[work.len() - 1].clone()];
    while work.len() > 1 {
        work = work
            .windows(2)
            .map(|pair| std::array::from_fn(|i| (&pair[0][i] + &pair[1][i]).mul_pow2(-1)))
            .collect();
        left.push(work[0].clone());
        right.push(work[work.len() - 1].clone());
    }
    right.reverse();
    (left, right)
}
fn net_bounds(net: &Net) -> (Rational, Rational) {
    let norm = |p: &[Rational; 4]| length2(&std::array::from_fn(|i| &p[i] / &p[3]));
    let upper = net.iter().map(norm).max().unwrap();
    let (left, _) = split(net);
    let lower = norm(&net[0])
        .max(norm(&net[net.len() - 1]))
        .max(norm(&left[left.len() - 1]));
    (lower, upper)
}

fn spline_movement(
    before: &RationalSplineCurve,
    after: &RationalSplineCurve,
    dec: &mut Decider,
) -> Result<(Rational, Rational), String> {
    if before == after {
        return Ok((Rational::zero(), Rational::zero()));
    }
    if before.knots() != after.knots() || before.weights() != after.weights() {
        return Err("source parameter correspondence changed".into());
    }
    let difference = before
        .controls()
        .iter()
        .zip(after.controls())
        .map(|(a, b)| sub(a, b))
        .collect();
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
        let net = curve
            .span_bezier_controls(span, dec)
            .map_err(|e| e.detail)?
            .into_iter()
            .map(|p| {
                let mut result = std::array::from_fn(|_| Rational::zero());
                for i in 0..4 {
                    result[i] = p[i]
                        .as_rational()
                        .ok_or("nonrational displacement control")?;
                }
                Ok(result)
            })
            .collect::<Result<Net, &str>>()
            .map_err(str::to_owned)?;
        nets.push(net);
    }
    let mut lower = Rational::zero();
    for _ in 0..4096 {
        dec.check_deadline().map_err(|e| e.detail)?;
        let bounds = nets.iter().map(net_bounds).collect::<Vec<_>>();
        lower = lower.max(bounds.iter().map(|b| b.0.clone()).max().unwrap());
        let (index, (_, upper)) = bounds
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.1.cmp(&b.1))
            .unwrap();
        let lo = crate::gap::sqrt_enclosure(&lower).0;
        let hi = crate::gap::sqrt_enclosure(upper).1;
        if upper == &lower || &hi - &lo <= Rational::from_ratio(1, 1_000_000_000_000).unwrap() {
            return Ok((lo, hi));
        }
        let (a, b) = split(&nets.swap_remove(index));
        nets.push(a);
        nets.push(b);
    }
    let upper = nets.iter().map(|n| net_bounds(n).1).max().unwrap();
    Ok((
        crate::gap::sqrt_enclosure(&lower).0,
        crate::gap::sqrt_enclosure(&upper).1,
    ))
}

pub(crate) fn prepare(source: &Solid) -> Result<Prepared, String> {
    let mut vertices: BTreeMap<Vec3, Vec<Placement>> = BTreeMap::new();
    let mut planes: BTreeMap<String, Vec<Placement>> = BTreeMap::new();
    let mut pairs: BTreeMap<(Vec3, Vec3), (String, bool)> = BTreeMap::new();
    for face in &source.faces {
        let Surface::Plane(p) = &face.surface else {
            return Err(
                "exact incidence repair currently requires planar support on every face".into(),
            );
        };
        for edge in face.loops.iter().flat_map(|b| &b.edges) {
            available()?;
            if matches!(edge.curve, Curve::Circle(..) | Curve::Ellipse(..)) {
                return Err(
                    "exact planar incidence repair requires line or spline boundaries".into(),
                );
            }
            let (a, b) = if edge.start.point <= edge.end.point {
                (&edge.start.point, &edge.end.point)
            } else {
                (&edge.end.point, &edge.start.point)
            };
            let exact_line = if let Curve::Line(o, d) = &edge.curve {
                is_zero(&cross(&sub(a, o), d)) && is_zero(&cross(&sub(b, o), d))
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
    let cap2 = &cap * &cap;
    let mut moved = BTreeMap::new();
    let mut movement = Rational::zero();
    for (point, planes) in vertices {
        available()?;
        let result = project(&point, &planes)?;
        let d2 = length2(&sub(&point, &result));
        if d2 > cap2 {
            return Err("an adjacent vertex needs movement above the 10 um import cap".into());
        }
        movement = movement.max(crate::gap::sqrt_enclosure(&d2).1);
        moved.insert(point, result);
    }
    let mut images = BTreeMap::new();
    for (original, image) in &moved {
        if images.insert(image.clone(), original).is_some() {
            return Err("the planar derivative would merge distinct source vertices; seam identity is unproved".into());
        }
    }
    let mut result = source.clone();
    let mut records = BTreeMap::new();
    let mut dec = Decider::default();
    for edge in result
        .faces
        .iter_mut()
        .flat_map(|f| &mut f.loops)
        .flat_map(|b| &mut b.edges)
    {
        available()?;
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
                let d2 = dot(d, d);
                let mut maximum = Rational::zero();
                for (old, new) in [
                    (&before.start.point, &edge.start.point),
                    (&before.end.point, &edge.end.point),
                ] {
                    let t = &dot(&sub(old, o), d) / &d2;
                    let source_point = std::array::from_fn(|i| &o[i] + &(&t * &d[i]));
                    maximum = maximum.max(length2(&sub(&source_point, new)));
                }
                edge.curve = Curve::Line(
                    edge.start.point.clone(),
                    sub(&edge.end.point, &edge.start.point),
                );
                crate::gap::sqrt_enclosure(&maximum)
            }
            Curve::BSpline(s) => {
                let original = s
                    .geometry
                    .clone()
                    .into_clamped(&dec)
                    .map_err(|e| e.detail)?;
                let mut controls = original
                    .controls()
                    .iter()
                    .map(|p| project(p, &planes[&edge.edge]))
                    .collect::<Result<Vec<_>, _>>()?;
                let last = controls.len() - 1;
                let ends = if edge.curve_same_sense {
                    [&edge.start.point, &edge.end.point]
                } else {
                    [&edge.end.point, &edge.start.point]
                };
                for (index, end) in [(0, ends[0]), (last, ends[1])] {
                    if length2(&sub(&original.controls()[index], end)) > cap2 {
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
        if hi > cap {
            return Err(format!(
                "edge {} needs a retained-curve deviation bound of {} mm, above the 10 um import cap",
                edge.edge,
                crate::gap::decimal(&hi)
            ));
        }
        movement = movement.max(hi.clone());
        records.insert(edge.edge.clone(), (lo, hi, edge.curve.clone()));
    }
    if &result == source {
        return Err("the planar derivative cannot change the failed incidence".into());
    }
    // The feature proof uses complete convex polygon covers. A nonstraight
    // boundary, hole or nonconvex face cannot acquire this certificate.
    let mut cover = result.clone();
    for edge in cover
        .faces
        .iter_mut()
        .flat_map(|f| &mut f.loops)
        .flat_map(|b| &mut b.edges)
    {
        if let Curve::BSpline(s) = &edge.curve {
            let chord = sub(&edge.end.point, &edge.start.point);
            if s.geometry
                .controls()
                .iter()
                .any(|p| !is_zero(&cross(&sub(p, &edge.start.point), &chord)))
            {
                return Err("nearest-feature proof needs straight repaired boundaries".into());
            }
            edge.curve = Curve::Line(edge.start.point.clone(), chord);
        }
    }
    let feature2 = cover
        .faces
        .iter()
        .map(|f| {
            crate::deviation::feature_size2(f)
                .ok_or("nearest-feature proof needs complete convex planar faces")
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .min()
        .ok_or("empty planar solid")?;
    let feature = crate::gap::sqrt_enclosure(&feature2).0 - &movement * &Rational::from_i64(2);
    if feature <= Rational::zero() || &movement * &Rational::from_i64(10) > feature {
        return Err("deviation exceeds a tenth of the certified neighbouring feature size".into());
    }
    let source_geometry_sha256 = crate::sha256_hex(&crate::canonical_bytes(&crate::StepFile {
        solids: vec![source.clone()],
        mm_per_unit: Rational::one(),
    }));
    let mut marks = MarkLog::new();
    for (edge, (lo, hi, curve)) in records {
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
            feature_size_squared_mm2: (&feature * &feature).to_string(),
        };
        marks.push(Mark::new(
            "import.incidence-repair",
            BoundKind::Distance,
            hi,
            serde_json::to_string(&record).expect("serializable edge tolerance"),
        ));
    }
    let vertices = source
        .faces
        .iter()
        .flat_map(|f| &f.loops)
        .flat_map(|b| &b.edges)
        .flat_map(|e| [&e.start, &e.end])
        .map(|v| (v.name.clone(), length2(&sub(&v.point, &moved[&v.point]))))
        .collect::<BTreeMap<_, _>>();
    for (entity, d2) in vertices {
        if d2.is_zero() {
            continue;
        }
        let (lo, hi) = crate::gap::sqrt_enclosure(&d2);
        let record = LocalTolerance {
            solid: source.name.clone(),
            entity,
            entity_kind: "vertex",
            source_geometry_sha256: source_geometry_sha256.clone(),
            definition: "Euclidean source vertex displacement",
            deviation_lower_mm: crate::gap::decimal(&lo),
            deviation_upper_mm: crate::gap::decimal(&hi),
            approved_tolerance_mm: crate::gap::decimal(&cap),
            feature_size_squared_mm2: (&feature * &feature).to_string(),
        };
        marks.push(Mark::new(
            "import.incidence-vertex-repair",
            BoundKind::Distance,
            hi,
            serde_json::to_string(&record).expect("serializable vertex tolerance"),
        ));
    }
    Ok(Prepared {
        solid: result,
        marks,
    })
}
