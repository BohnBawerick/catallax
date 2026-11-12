//! Whole trimmed-face bounds for repaired convex planar polygons.
use catallax_numbers::{Rational, Sign};

use crate::exact::{cross, dot, is_zero, sub};
use crate::gap::{decimal, sqrt_enclosure};
use crate::{
    Curve, Face, FaceTolerance, ImportRefusal, PlanarRepairMethod, Solid, SourceIncidence, Surface,
    Vec3,
};

/// Checks a complete convex polygon cover, including straight incidence and order.
fn polygon(face: &Face) -> Option<Vec<&Vec3>> {
    let Surface::Plane(plane) = &face.surface else {
        return None;
    };
    let [bound] = face.loops.as_slice() else {
        return None;
    };
    if bound.edges.len() < 3 || is_zero(&plane.axis) {
        return None;
    }
    let mut points = Vec::with_capacity(bound.edges.len());
    for (i, edge) in bound.edges.iter().enumerate() {
        let Curve::Line(origin, direction) = &edge.curve else {
            return None;
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
        if to != next_start
            || from == to
            || is_zero(direction)
            || !is_zero(&cross(&sub(from, origin), direction))
            || !is_zero(&cross(&sub(to, origin), direction))
            || !dot(&sub(from, &plane.origin), &plane.axis).is_zero()
        {
            return None;
        }
        points.push(from);
    }
    // Every nonincident vertex must lie strictly on the same side of every
    // oriented edge. This proves simplicity, strict convexity and fan coverage.
    let mut side = None;
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let edge = sub(points[j], points[i]);
        for k in 0..points.len() {
            if k == i || k == j {
                continue;
            }
            let sign = dot(&cross(&edge, &sub(points[k], points[i])), &plane.axis).sign();
            if sign == Sign::Zero || side.is_some_and(|s| s != sign) {
                return None;
            }
            side = Some(sign);
        }
    }
    Some(points)
}

/// Conservative squared feature size: every edge length and every nonincident
/// vertex's distance to an edge's supporting line. Narrow skew faces count too.
pub(crate) fn feature_size2(face: &Face) -> Option<Rational> {
    let points = polygon(face)?;
    let mut size: Option<Rational> = None;
    for i in 0..points.len() {
        let j = (i + 1) % points.len();
        let edge = sub(points[j], points[i]);
        let length2 = dot(&edge, &edge);
        size = Some(size.map_or_else(|| length2.clone(), |s| s.min(length2.clone())));
        for k in 0..points.len() {
            if k == i || k == j {
                continue;
            }
            let area = cross(&edge, &sub(points[k], points[i]));
            let altitude2 = &dot(&area, &area) / &length2;
            size = Some(size.unwrap().min(altitude2));
        }
    }
    size
}

/// Exact squared distance to a covered convex polygon, including its trim.
fn distance2(point: &Vec3, polygon: &[&Vec3]) -> Rational {
    let normal = cross(&sub(polygon[1], polygon[0]), &sub(polygon[2], polygon[0]));
    let normal2 = dot(&normal, &normal);
    let inside = (0..polygon.len()).all(|i| {
        let edge = sub(polygon[(i + 1) % polygon.len()], polygon[i]);
        let side = dot(&cross(&edge, &sub(point, polygon[i])), &normal).sign();
        side != Sign::Negative
    });
    if inside {
        let height = dot(&sub(point, polygon[0]), &normal);
        return &height * &height / normal2;
    }
    (0..polygon.len())
        .map(|i| {
            let edge = sub(polygon[(i + 1) % polygon.len()], polygon[i]);
            let delta = sub(point, polygon[i]);
            let t = (dot(&delta, &edge) / dot(&edge, &edge))
                .max(Rational::zero())
                .min(Rational::one());
            let residual = std::array::from_fn(|j| &delta[j] - &(&t * &edge[j]));
            dot(&residual, &residual)
        })
        .min()
        .expect("polygon has edges")
}

pub(crate) fn certify(
    source: &Solid,
    repaired: &Solid,
    cap: &Rational,
) -> Result<Vec<FaceTolerance>, String> {
    let mut records = Vec::new();
    for (before, after) in source.faces.iter().zip(&repaired.faces) {
        if before == after {
            continue;
        }
        let a = polygon(before)
            .ok_or_else(|| format!("face {} needs a convex planar cover", before.name))?;
        let b = polygon(after)
            .ok_or_else(|| format!("face {} needs a convex planar cover", after.name))?;
        if a.len() != b.len() || before.name != after.name {
            return Err("the repair changed the face correspondence".into());
        }
        // Distance to a convex set is convex. For every fan triangle, combine
        // the three vertices' closest target points with the same barycentric
        // weights. The combination lies in the target polygon and its distance
        // is no larger than the largest vertex distance. Vertices also supply
        // exact lower witnesses. Apply this in both directions over both covers.
        let distance2 = a
            .iter()
            .map(|p| distance2(p, &b))
            .chain(b.iter().map(|p| distance2(p, &a)))
            .max()
            .expect("covered polygons");
        if distance2.is_zero() {
            continue;
        }
        let (lower, upper) = sqrt_enclosure(&distance2);
        let mut feature2 = feature_size2(before).expect("checked polygon");
        feature2 = feature2.min(feature_size2(after).expect("checked polygon"));
        // Include neighbouring faces, identified by shared source vertices.
        for (neighbour, result_neighbour) in source.faces.iter().zip(&repaired.faces) {
            if neighbour.name == before.name {
                continue;
            }
            let adjacent = neighbour
                .loops
                .iter()
                .flat_map(|l| &l.edges)
                .any(|e| a.iter().any(|p| **p == e.start.point || **p == e.end.point));
            if adjacent {
                feature2 = feature2.min(feature_size2(neighbour).ok_or_else(|| {
                    format!(
                        "nearest feature of face {} could not be bounded",
                        neighbour.name
                    )
                })?);
                feature2 = feature2.min(feature_size2(result_neighbour).ok_or_else(|| {
                    format!(
                        "nearest feature of repaired face {} could not be bounded",
                        neighbour.name
                    )
                })?);
            }
        }
        if upper > *cap {
            return Err(format!(
                "face {} exceeds its approved tolerance",
                before.name
            ));
        }
        if &upper * &upper * Rational::from_i64(100) > feature2 {
            return Err(format!(
                "face {} exceeds a tenth of the nearest feature",
                before.name
            ));
        }
        records.push(FaceTolerance {
            solid: source.name.clone(),
            face: before.name.clone(),
            deviation_lower_mm: decimal(&lower),
            deviation_upper_mm: decimal(&upper),
            approved_tolerance_mm: decimal(cap),
            feature_size_squared_mm2: feature2.to_string(),
        });
    }
    Ok(records)
}

fn squared(a: &Vec3, b: &Vec3) -> Rational {
    let delta = sub(a, b);
    dot(&delta, &delta)
}

/// Defines the finite source patch explicitly, without changing the retained bytes.
fn projected(face: &Face) -> Option<Face> {
    let Surface::Plane(plane) = &face.surface else {
        return None;
    };
    let normal2 = dot(&plane.axis, &plane.axis);
    if normal2.is_zero() {
        return None;
    }
    let project = |point: &Vec3| {
        let t = &dot(&plane.axis, &sub(point, &plane.origin)) / &normal2;
        std::array::from_fn(|i| &point[i] - &(&t * &plane.axis[i]))
    };
    let mut result = face.clone();
    for edge in result.loops.iter_mut().flat_map(|b| &mut b.edges) {
        let Curve::Line(origin, direction) = &mut edge.curve else {
            return None;
        };
        edge.start.point = project(&edge.start.point);
        edge.end.point = project(&edge.end.point);
        let delta = sub(&edge.end.point, &edge.start.point);
        if is_zero(&delta) {
            return None;
        }
        *origin = edge.start.point.clone();
        *direction = if edge.curve_same_sense {
            delta
        } else {
            delta.map(|c| -c)
        };
    }
    Some(result)
}

/// Bounds the complete vertex fan and finite line pieces by their affine maps.
fn source_correspondence(
    before: &Face,
    defined: &Face,
    after: &Face,
) -> Option<(Rational, Rational, Rational)> {
    if before.name != after.name
        || before.same_sense != after.same_sense
        || before.loops.len() != after.loops.len()
    {
        return None;
    }
    let (mut uncertainty2, mut movement2, mut total2) =
        (Rational::zero(), Rational::zero(), Rational::zero());
    for ((a, d), b) in before.loops.iter().zip(&defined.loops).zip(&after.loops) {
        if a.edges.len() != b.edges.len() || a.outer != b.outer || a.orientation != b.orientation {
            return None;
        }
        for ((a, d), b) in a.edges.iter().zip(&d.edges).zip(&b.edges) {
            if a.edge != b.edge
                || a.start.name != b.start.name
                || a.end.name != b.end.name
                || a.forward != b.forward
                || a.curve_same_sense != b.curve_same_sense
            {
                return None;
            }
            let Curve::Line(origin, direction) = &a.curve else {
                return None;
            };
            let length2 = dot(direction, direction);
            if length2.is_zero() {
                return None;
            }
            for ((v, q), target) in [&a.start.point, &a.end.point]
                .into_iter()
                .zip([&d.start.point, &d.end.point])
                .zip([&b.start.point, &b.end.point])
            {
                let t = &dot(&sub(v, origin), direction) / &length2;
                let on_line = std::array::from_fn(|i| &origin[i] + &(&t * &direction[i]));
                uncertainty2 = uncertainty2.max(squared(v, q)).max(squared(&on_line, q));
                movement2 = movement2.max(squared(v, target));
                total2 = total2
                    .max(squared(v, target))
                    .max(squared(&on_line, target));
            }
        }
    }
    Some((uncertainty2, movement2, total2))
}

pub(crate) fn certify_incidence(
    source: &Solid,
    repaired: &Solid,
    method: PlanarRepairMethod,
) -> Result<(Vec<FaceTolerance>, Vec<SourceIncidence>), ImportRefusal> {
    if source.faces.len() != repaired.faces.len() {
        return Err(ImportRefusal::NotSolid {
            detail: "the planar repair changed the source face count".into(),
        });
    }
    let mut tolerances = Vec::new();
    let mut records = Vec::new();
    let cap = crate::gap::auto_cap();
    for (before, after) in source.faces.iter().zip(&repaired.faces) {
        catallax_numbers::Decider::default()
            .check_deadline()
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        if before == after {
            continue;
        }
        let fail = |stage: &str, measurement: &str, value: Option<&Rational>, reason: String| {
            crate::planar_incidence::refused(source, before, stage, measurement, value, reason)
        };
        let defined = projected(before).ok_or_else(|| {
            fail(
                "source_cover",
                "unavailable",
                None,
                "the source needs a straight planar boundary interpretation".into(),
            )
        })?;
        let (uncertainty2, movement2, correspondence2) =
            source_correspondence(before, &defined, after).ok_or_else(|| {
                fail(
                    "source_curve",
                    "unavailable",
                    None,
                    "the source boundary correspondence could not be bounded".into(),
                )
            })?;
        let a = polygon(&defined).ok_or_else(|| {
            fail(
                "source_cover",
                "source_boundary_correspondence",
                Some(&uncertainty2),
                "the projected source needs a complete convex single-bound cover".into(),
            )
        })?;
        let b = polygon(after).ok_or_else(|| {
            fail(
                "result_cover",
                "source_boundary_correspondence",
                Some(&uncertainty2),
                "the repaired face needs a complete convex single-bound cover".into(),
            )
        })?;
        if a.len() != b.len() {
            return Err(fail(
                "face_correspondence",
                "unavailable",
                None,
                "the repair changed the face correspondence".into(),
            ));
        }
        let geometry2 = a
            .iter()
            .map(|p| distance2(p, &b))
            .chain(b.iter().map(|p| distance2(p, &a)))
            .max()
            .expect("covered faces");
        // The two convex patches use complete covers. The raw vertex fan and
        // finite source line pieces use corresponding affine maps. Taking the
        // maximum bounds every source constraint directly, without adding errors.
        let total2 = correspondence2.max(geometry2.clone());
        if total2 > &cap * &cap {
            return Err(fail(
                "import_cap",
                "whole_face_source_bound",
                Some(&total2),
                "the proposed source bound exceeds 0.01 mm".into(),
            ));
        }
        let mut feature2 = feature_size2(&defined)
            .expect("checked source cover")
            .min(feature_size2(after).expect("checked result cover"));
        let vertices: Vec<_> = before
            .loops
            .iter()
            .flat_map(|l| &l.edges)
            .flat_map(|e| [&e.start.point, &e.end.point])
            .collect();
        for (neighbour, result_neighbour) in source.faces.iter().zip(&repaired.faces) {
            if neighbour.name == before.name {
                continue;
            }
            let adjacent = neighbour.loops.iter().flat_map(|l| &l.edges).any(|e| {
                vertices
                    .iter()
                    .any(|p| **p == e.start.point || **p == e.end.point)
            });
            if !adjacent {
                continue;
            }
            let projected_neighbour = projected(neighbour).ok_or_else(|| {
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
            for face in [&projected_neighbour, result_neighbour] {
                feature2 = feature2.min(feature_size2(face).ok_or_else(|| {
                    fail(
                        "nearest_feature",
                        "whole_face_source_bound",
                        Some(&total2),
                        format!(
                            "neighbour {} needs a complete convex feature bound",
                            face.name
                        ),
                    )
                })?);
            }
        }
        let (_, total_upper) = sqrt_enclosure(&total2);
        if &total_upper * &total_upper * Rational::from_i64(100) > feature2 {
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
        let (geometry_lower, geometry_upper) = sqrt_enclosure(&geometry2);
        let (uncertainty_lower, uncertainty_upper) = sqrt_enclosure(&uncertainty2);
        let (_, movement_upper) = sqrt_enclosure(&movement2);
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
            source_discrepancy_lower_mm: decimal(&uncertainty_lower),
            source_discrepancy_upper_mm: decimal(&uncertainty_upper),
            geometric_deviation_lower_mm: decimal(&geometry_lower),
            geometric_deviation_upper_mm: decimal(&geometry_upper),
            source_to_result_upper_mm: decimal(&total_upper),
            vertex_movement_upper_mm: decimal(&movement_upper),
        });
    }
    Ok((tolerances, records))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FaceLoop, OrientedEdge, Placement, Vertex};

    fn face(points: &[[&str; 3]]) -> Face {
        let points: Vec<Vec3> = points
            .iter()
            .map(|p| p.map(|x| Rational::parse(x).unwrap()))
            .collect();
        let normal = cross(&sub(&points[1], &points[0]), &sub(&points[2], &points[0]));
        Face {
            name: "#face".into(),
            surface: Surface::Plane(Placement {
                origin: points[0].clone(),
                axis: normal,
                ref_dir: sub(&points[1], &points[0]),
            }),
            same_sense: true,
            loops: vec![FaceLoop {
                outer: true,
                orientation: true,
                edges: (0..points.len())
                    .map(|i| {
                        let j = (i + 1) % points.len();
                        OrientedEdge {
                            edge: format!("#{i}"),
                            start: Vertex {
                                name: format!("v{i}"),
                                point: points[i].clone(),
                            },
                            end: Vertex {
                                name: format!("v{j}"),
                                point: points[j].clone(),
                            },
                            curve_same_sense: true,
                            forward: true,
                            curve: Curve::Line(points[i].clone(), sub(&points[j], &points[i])),
                        }
                    })
                    .collect(),
                vertex: None,
            }],
        }
    }

    #[test]
    fn source_constraints_are_bounded_directly_without_adding_discrepancies() {
        let after = face(&[
            ["0", "0", "0.006"],
            ["1", "0", "0.006"],
            ["1", "1", "0.006"],
            ["0", "1", "0.006"],
        ]);
        let mut before = after.clone();
        let Surface::Plane(plane) = &mut before.surface else {
            unreachable!()
        };
        plane.origin[2] = Rational::zero();
        for edge in &mut before.loops[0].edges {
            let Curve::Line(origin, _) = &mut edge.curve else {
                unreachable!()
            };
            origin[2] = Rational::parse("0.012").unwrap();
        }
        let source = Solid {
            name: "solid".into(),
            faces: vec![before],
        };
        let result = Solid {
            name: "solid".into(),
            faces: vec![after],
        };
        let (bounds, records) =
            certify_incidence(&source, &result, PlanarRepairMethod::RefitPlanes).unwrap();
        assert_eq!(records[0].source_discrepancy_upper_mm, "0.012");
        assert_eq!(records[0].geometric_deviation_upper_mm, "0.006");
        assert_eq!(records[0].source_to_result_upper_mm, "0.006");
        assert_eq!(bounds[0].deviation_upper_mm, "0.006");
    }

    #[test]
    fn a_boundary_carrier_change_is_not_lost_when_geometric_deviation_is_zero() {
        let after = face(&[
            ["0", "0", "0"],
            ["1", "0", "0"],
            ["1", "1", "0"],
            ["0", "1", "0"],
        ]);
        let mut before = after.clone();
        for edge in &mut before.loops[0].edges {
            let Curve::Line(origin, _) = &mut edge.curve else {
                unreachable!()
            };
            origin[2] = Rational::parse("0.004").unwrap();
        }
        let source = Solid {
            name: "solid".into(),
            faces: vec![before],
        };
        let result = Solid {
            name: "solid".into(),
            faces: vec![after],
        };
        let (bounds, records) =
            certify_incidence(&source, &result, PlanarRepairMethod::RefitPlanes).unwrap();
        assert_eq!(records[0].geometric_deviation_upper_mm, "0");
        assert_eq!(records[0].vertex_movement_upper_mm, "0");
        assert_eq!(records[0].source_to_result_upper_mm, "0.004");
        assert_eq!(bounds[0].deviation_upper_mm, "0.004");
    }

    #[test]
    fn trimmed_polygon_distance_is_not_vertex_correspondence_or_carrier_distance() {
        let a = face(&[
            ["0", "0", "0"],
            ["10", "0", "0"],
            ["10", "10", "0"],
            ["0", "10", "0"],
        ]);
        // Cyclically different vertex correspondence, same complete face.
        let reordered = face(&[
            ["10", "0", "0"],
            ["10", "10", "0"],
            ["0", "10", "0"],
            ["0", "0", "0"],
        ]);
        let solid = |f| Solid {
            name: "#solid".into(),
            faces: vec![f],
        };
        assert!(
            certify(&solid(a.clone()), &solid(reordered), &Rational::one())
                .unwrap()
                .is_empty()
        );
        // Same carrier plane, nonzero deviation of the trimmed faces.
        let b = face(&[
            ["0.004", "0", "0"],
            ["10.004", "0", "0"],
            ["10.004", "10", "0"],
            ["0.004", "10", "0"],
        ]);
        let records = certify(&solid(a), &solid(b), &Rational::one()).unwrap();
        assert_eq!(records[0].deviation_lower_mm, "0.004");
        assert_eq!(records[0].deviation_upper_mm, "0.004");
    }

    #[test]
    fn spatial_distance_and_narrow_skew_features_are_bounded_exactly() {
        let a = face(&[
            ["0", "0", "0"],
            ["10", "0", "0"],
            ["10", "10", "0"],
            ["0", "10", "0"],
        ]);
        let points = polygon(&a).unwrap();
        assert_eq!(
            distance2(
                &["10.003", "10", "0.004"].map(|s| s.parse().unwrap()),
                &points
            ),
            Rational::parse("0.000025").unwrap()
        );
        let skew = face(&[
            ["0", "0", "0"],
            ["10", "0", "0"],
            ["20", "0.001", "0"],
            ["10", "0.001", "0"],
        ]);
        assert!(feature_size2(&skew).unwrap() < Rational::parse("0.000001").unwrap());
        let concave = face(&[
            ["0", "0", "0"],
            ["10", "0", "0"],
            ["1", "1", "0"],
            ["0", "10", "0"],
        ]);
        assert!(polygon(&concave).is_none());
    }
}
