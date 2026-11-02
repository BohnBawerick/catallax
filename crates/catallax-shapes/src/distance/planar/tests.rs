use super::*;

fn p(x: i64, y: i64, z: i64) -> Point3 {
    Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
}
fn plane() -> Plane {
    Plane {
        origin: p(0, 0, 0),
        normal: Vector3::new(0.into(), 0.into(), 7.into()),
        ref_dir: Vector3::new(3.into(), 0.into(), 0.into()),
    }
}
fn polygon(loops: &[&[[i64; 2]]]) -> PlanarBoundary {
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    for cycle in loops {
        let first = vertices.len();
        vertices.extend(cycle.iter().map(|[x, y]| p(*x, *y, 0)));
        edges.extend(
            (0..cycle.len())
                .map(|i| PlanarEdge::Segment([first + i, first + (i + 1) % cycle.len()])),
        );
    }
    PlanarBoundary {
        plane: plane(),
        vertices,
        edges,
    }
}
fn circle(radius: i64) -> Circle {
    Circle {
        center: p(0, 0, 0),
        normal: Vector3::new(0.into(), 0.into(), 1.into()),
        radius_sq: Real::from_i64(radius * radius),
        ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
    }
}
fn distance(boundary: &PlanarBoundary, point: Point3, squared: i64) {
    let mut dec = Decider::default();
    let controls = Options {
        width: Rational::parse("1e-12").unwrap(),
        max_cells: 1,
    };
    let answer = boundary
        .point_distance(&point, &controls, &mut dec)
        .unwrap();
    let expected = dec
        .enclose_at(&Real::from_i64(squared).sqrt(), 256)
        .unwrap();
    let upper = answer.upper.as_ref().unwrap();
    assert!(
        answer.lower <= expected.lo() && upper >= &expected.hi(),
        "{answer:?}"
    );
    assert!(
        answer.converged && upper - &answer.lower <= controls.width,
        "{answer:?}"
    );
    let witness = &answer.witnesses.as_ref().unwrap()[1];
    assert!(boundary.contains(witness, &mut dec).unwrap());
    assert_eq!(
        math::sign(
            &(point.minus(witness).norm_squared() - Real::from_i64(squared)),
            &mut dec
        )
        .unwrap(),
        Sign::Zero
    );
}

#[test]
fn nonconvex_polygon_and_holes_reject_the_infinite_plane_projection() {
    let l = polygon(&[&[[0, 0], [4, 0], [4, 1], [1, 1], [1, 4], [0, 4]]]);
    distance(&l, p(3, 3, 2), 8);
    distance(&l, p(0, 4, 0), 0);
    distance(&l, p(0, 2, 3), 9);
    let hole = polygon(&[
        &[[0, 0], [6, 0], [6, 6], [0, 6]],
        &[[1, 1], [1, 3], [3, 3], [3, 1]],
    ]);
    distance(&hole, p(2, 2, 1), 2);
    distance(&hole, p(4, 4, 2), 4);
    distance(&hole, p(1, 2, 0), 0);
}

#[test]
fn circular_arcs_use_the_closed_sweep_not_the_whole_circle() {
    let upper = PlanarBoundary {
        plane: plane(),
        vertices: vec![p(2, 0, 0), p(-2, 0, 0)],
        edges: vec![
            PlanarEdge::Arc {
                circle: circle(2),
                ends: Some([0, 1]),
            },
            PlanarEdge::Segment([1, 0]),
        ],
    };
    distance(&upper, p(0, -1, 3), 10);
    distance(&upper, p(0, 3, 4), 17);
    distance(&upper, p(-3, -1, 0), 2);
    let major = PlanarBoundary {
        plane: plane(),
        vertices: vec![p(2, 0, 0), p(0, -2, 0)],
        edges: vec![
            PlanarEdge::Arc {
                circle: circle(2),
                ends: Some([0, 1]),
            },
            PlanarEdge::Segment([1, 0]),
        ],
    };
    distance(&major, p(3, -3, 1), 9);
    distance(&major, p(-2, 0, 1), 1);
}

#[test]
fn complete_circular_holes_have_axis_and_rim_witnesses() {
    let annulus = PlanarBoundary {
        plane: plane(),
        vertices: vec![],
        edges: vec![
            PlanarEdge::Arc {
                circle: circle(4),
                ends: None,
            },
            PlanarEdge::Arc {
                circle: circle(2),
                ends: None,
            },
        ],
    };
    distance(&annulus, p(0, 0, 3), 13);
    distance(&annulus, p(3, 0, 3), 9);
    distance(&annulus, p(5, 0, 0), 1);
}

#[test]
fn polygon_pair_minima_cover_interior_overlaps_skew_crossings_and_holes() {
    let controls = Options {
        width: Rational::parse("1e-10").unwrap(),
        max_cells: 1,
    };
    let horizontal = polygon(&[&[[-3, -1], [3, -1], [3, 1], [-3, 1]]]);
    let mut across = polygon(&[&[[-1, -3], [1, -3], [1, 3], [-1, 3]]]);
    across.plane.origin = p(0, 0, 2);
    for vertex in &mut across.vertices {
        vertex.z = 2.into();
    }
    let mut dec = Decider::default();
    let answer = horizontal
        .face_distance(&across, &controls, &mut dec)
        .unwrap()
        .unwrap();
    assert_eq!(
        (answer.lower, answer.upper),
        (Rational::from_i64(2), Some(Rational::from_i64(2)))
    );
    let vertical = PlanarBoundary {
        plane: Plane {
            origin: p(0, 0, 0),
            normal: Vector3::new(1.into(), 0.into(), 0.into()),
            ref_dir: Vector3::new(0.into(), 1.into(), 0.into()),
        },
        vertices: vec![p(0, -1, -1), p(0, 1, -1), p(0, 1, 1), p(0, -1, 1)],
        edges: vec![
            PlanarEdge::Segment([0, 1]),
            PlanarEdge::Segment([1, 2]),
            PlanarEdge::Segment([2, 3]),
            PlanarEdge::Segment([3, 0]),
        ],
    };
    let answer = horizontal
        .face_distance(&vertical, &controls, &mut dec)
        .unwrap()
        .unwrap();
    assert_eq!(
        (answer.lower, answer.upper),
        (Rational::zero(), Some(Rational::zero()))
    );
    let hole = polygon(&[
        &[[0, 0], [6, 0], [6, 6], [0, 6]],
        &[[1, 1], [1, 3], [3, 3], [3, 1]],
    ]);
    let mut inside_hole = polygon(&[&[[0, 0], [1, 0], [1, 1], [0, 1]]]);
    let shift = Vector3::new(
        Real::from(Rational::parse("1.5").unwrap()),
        Real::from(Rational::parse("1.5").unwrap()),
        2.into(),
    );
    inside_hole.plane.origin = inside_hole.plane.origin.plus(&shift);
    for vertex in &mut inside_hole.vertices {
        *vertex = vertex.plus(&shift);
    }
    let answer = hole
        .face_distance(&inside_hole, &controls, &mut dec)
        .unwrap()
        .unwrap();
    let expected = dec
        .enclose_at(&(&Real::from_i64(17) / &Real::from_i64(4)).sqrt(), 256)
        .unwrap();
    assert!(answer.lower <= expected.lo() && answer.upper.as_ref().unwrap() >= &expected.hi());
    assert!(answer.converged);
}
