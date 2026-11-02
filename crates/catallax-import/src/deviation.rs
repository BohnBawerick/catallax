//! Whole trimmed-face bounds for repaired convex planar polygons.
use catallax_numbers::{Rational, Sign};

use crate::exact::{cross, dot, is_zero, sub};
use crate::gap::{decimal, sqrt_enclosure};
use crate::{Curve, Face, FaceTolerance, Solid, Surface, Vec3};

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
            }],
        }
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
