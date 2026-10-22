//! Closed-form profiles, naming and impossible sketch regressions.

use catallax_features::{Profile, extrude};
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{ElementId, Loop, Point2, Pose, Segment, StepId, Turn};
use catallax_topology::{BodyId, Model};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn p(x: i64, y: i64) -> Point2 {
    Point2([x.into(), y.into()])
}
fn id(s: &str) -> ElementId {
    ElementId::new(s).unwrap()
}
fn rectangle(x: i64, y: i64, width: i64, height: i64) -> Loop {
    Loop(
        [
            (x + width, y),
            (x + width, y + height),
            (x, y + height),
            (x, y),
        ]
        .into_iter()
        .enumerate()
        .map(|(i, (x, y))| Segment::Line {
            id: id(&format!("s{i}")),
            to: p(x, y),
        })
        .collect(),
    )
}
fn circle(name: &str, x: i64, y: i64, radius: i64) -> Loop {
    Loop(vec![Segment::Circle {
        id: id(name),
        centre: p(x, y),
        radius: radius.into(),
    }])
}
fn build(loops: &[Loop], height: i64) -> (Model, BodyId) {
    let mut dec = Decider::default();
    let profile = Profile::new(&Pose::identity(), loops, &mut dec).unwrap();
    let mut model = Model::new();
    let body = extrude(
        &mut model,
        &profile,
        &StepId::new("sweep").unwrap(),
        &r(height),
        None,
    )
    .unwrap();
    assert!(model.check(body).is_empty(), "{:?}", model.check(body));
    for face in model.faces_of_body(body) {
        assert!(!model.face_origins(face).is_empty());
    }
    (model, body)
}
fn measure(model: &Model, body: BodyId) -> catallax_measure::MassProperties {
    let answer = catallax_measure::mass_properties(
        model,
        body,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    answer.value
}
fn contains(e: &catallax_measure::Enclosure, expected: Ball) {
    assert!(
        e.lo <= expected.lo() && e.hi >= expected.hi(),
        "{e:?} versus {expected:?}"
    );
    assert!(e.width() <= Rational::parse("1e-9").unwrap());
}
fn pi() -> Ball {
    Ball::pi(192)
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&r(n), 192)
}

#[test]
fn rectangle_signed_distance_keeps_start_names_and_exact_properties() {
    for height in [3, -3] {
        let (model, body) = build(&[rectangle(0, 0, 4, 2)], height);
        let mass = measure(&model, body);
        assert_eq!(mass.volume.lo, r(24));
        assert_eq!(mass.volume.hi, r(24));
        assert_eq!(mass.area.lo, r(52));
        let centroid = mass.centroid.unwrap();
        for (value, expected) in
            centroid
                .iter()
                .zip([r(2), r(1), Rational::from_ratio(height, 2).unwrap()])
        {
            assert!(value.contains(&expected));
        }
        for edge in model.edges_of_body(body) {
            assert!(model.edge_name(body, edge).is_some());
        }
        for vertex in model.vertices_of_body(body) {
            assert!(model.vertex_name(body, vertex).is_some());
        }
    }
}

#[test]
fn circular_sweep_boundaries_keep_their_axis_under_pose() {
    for rotation in [[1i64, 0, 0, 0], [1, 2, 3, 4]] {
        let pose = Pose {
            rotation: rotation.map(Into::into),
            translation: catallax_recipe::Vector3([5.into(), 6.into(), 7.into()]),
        };
        let profile =
            Profile::new(&pose, &[circle("outer", 0, 0, 3)], &mut Decider::default()).unwrap();
        let mut model = Model::new();
        let body = extrude(
            &mut model,
            &profile,
            &StepId::new("sweep").unwrap(),
            &r(2),
            None,
        )
        .unwrap();
        let surface = model
            .faces_of_body(body)
            .into_iter()
            .map(|face| model.face(face).surface)
            .find(|surface| {
                matches!(
                    model.geom().surface(*surface),
                    catallax_shapes::Surface::Cylinder(_)
                )
            })
            .unwrap();
        for edge in model.edges_of_body(body) {
            assert!(model.geom().facts().holds(
                catallax_numbers::FactKind::Coaxial,
                catallax_shapes::GeomKey::Curve(model.edge(edge).curve),
                catallax_shapes::GeomKey::Surface(surface)
            ));
        }
    }
}

#[test]
fn circle_and_ring_have_closed_form_mass_properties() {
    let (model, body) = build(&[circle("outer", 0, 0, 3)], 2);
    let mass = measure(&model, body);
    contains(&mass.volume, b(18) * pi());
    contains(&mass.area, b(30) * pi());
    for (value, expected) in mass.centroid.unwrap().iter().zip([0, 0, 1]) {
        assert!(value.contains(&r(expected)));
    }
    let (model, body) = build(&[circle("hole", 0, 0, 1), circle("outer", 0, 0, 3)], 2);
    let mass = measure(&model, body);
    contains(&mass.volume, b(16) * pi());
    contains(&mass.area, b(32) * pi());
    for (value, expected) in mass.centroid.unwrap().iter().zip([0, 0, 1]) {
        assert!(value.contains(&r(expected)));
    }
}

#[test]
fn slot_uses_exact_arcs_not_polygonal_chords() {
    let slot = Loop(vec![
        Segment::Line {
            id: id("bottom"),
            to: p(4, -1),
        },
        Segment::Arc {
            id: id("right"),
            centre: p(4, 0),
            to: p(4, 1),
            turn: Turn::Ccw,
        },
        Segment::Line {
            id: id("top"),
            to: p(0, 1),
        },
        Segment::Arc {
            id: id("left"),
            centre: p(0, 0),
            to: p(0, -1),
            turn: Turn::Ccw,
        },
    ]);
    let (model, body) = build(&[slot], 3);
    let mass = measure(&model, body);
    contains(&mass.volume, b(24) + b(3) * pi());
    contains(&mass.area, b(40) + b(8) * pi());
    for (value, expected) in
        mass.centroid
            .unwrap()
            .iter()
            .zip([r(2), r(0), Rational::from_ratio(3, 2).unwrap()])
    {
        assert!(value.contains(&expected));
    }
}

#[test]
fn impossible_sketches_have_distinct_stable_reasons() {
    let bow = Loop(vec![
        Segment::Line {
            id: id("a"),
            to: p(2, 2),
        },
        Segment::Line {
            id: id("b"),
            to: p(0, 2),
        },
        Segment::Line {
            id: id("c"),
            to: p(2, 0),
        },
        Segment::Line {
            id: id("d"),
            to: p(0, 0),
        },
    ]);
    let backtrack = Loop(vec![
        Segment::Line {
            id: id("a"),
            to: p(1, 0),
        },
        Segment::Line {
            id: id("b"),
            to: p(0, 0),
        },
    ]);
    for (loops, code) in [
        (vec![bow], "sketch.self-intersection"),
        (vec![backtrack], "sketch.self-intersection"),
        (
            vec![circle("a", 0, 0, 2), circle("b", 4, 0, 2)],
            "sketch.loop-intersection",
        ),
        (vec![circle("a", 0, 0, 0)], "sketch.non-positive-radius"),
        (
            vec![Loop(vec![Segment::Line {
                id: id("a"),
                to: p(0, 0),
            }])],
            "sketch.zero-segment",
        ),
        (
            vec![Loop(vec![
                Segment::Arc {
                    id: id("a"),
                    centre: p(0, 0),
                    to: p(1, 0),
                    turn: Turn::Ccw,
                },
                Segment::Line {
                    id: id("b"),
                    to: p(0, 2),
                },
            ])],
            "sketch.arc-radius",
        ),
    ] {
        assert_eq!(
            Profile::new(&Pose::identity(), &loops, &mut Decider::default())
                .unwrap_err()
                .code,
            code
        );
    }
}

#[test]
fn empty_profile_is_empty_but_zero_distance_is_still_impossible() {
    let profile = Profile::new(&Pose::identity(), &[], &mut Decider::default()).unwrap();
    let mut model = Model::new();
    let step = StepId::new("empty").unwrap();
    let body = extrude(&mut model, &profile, &step, &r(1), None).unwrap();
    assert_eq!(model.body(body).kind, catallax_topology::BodyKind::Empty);
    assert_eq!(
        extrude(&mut model, &profile, &step, &r(0), None)
            .unwrap_err()
            .code,
        "extrude.zero-distance"
    );
}

#[test]
fn pythagorean_oblique_polygon_extrusion_has_exact_properties() {
    let profile = Profile::new(
        &Pose::identity(),
        &[rectangle(0, 0, 4, 2)],
        &mut Decider::default(),
    )
    .unwrap();
    let mut model = Model::new();
    let direction = catallax_shapes::RVector3::new(r(0), r(3), r(4));
    let body = extrude(
        &mut model,
        &profile,
        &StepId::new("slant").unwrap(),
        &r(5),
        Some(&direction),
    )
    .unwrap();
    assert!(model.check(body).is_empty());
    let properties = measure(&model, body);
    assert!(properties.volume.contains(&r(32)));
    assert!(properties.area.contains(&r(72)));
    for (value, expected) in
        properties
            .centroid
            .unwrap()
            .iter()
            .zip([r(2), Rational::from_ratio(5, 2).unwrap(), r(2)])
    {
        assert!(value.contains(&expected));
    }
}

#[test]
fn odd_nesting_preserves_islands_inside_holes() {
    let (model, body) = build(
        &[
            circle("hole", 0, 0, 3),
            circle("island", 0, 0, 1),
            circle("outside", 0, 0, 5),
        ],
        2,
    );
    assert_eq!(model.counts(body).lumps, 2);
    let properties = measure(&model, body);
    contains(&properties.volume, b(34) * pi());
    contains(&properties.area, b(70) * pi());
}

#[test]
fn extrusion_rejects_zero_and_in_plane_directions_before_construction() {
    let profile = Profile::new(
        &Pose::identity(),
        &[rectangle(0, 0, 4, 2)],
        &mut Decider::default(),
    )
    .unwrap();
    let step = StepId::new("sweep").unwrap();
    for (direction, code) in [
        (catallax_shapes::RVector3::zero(), "extrude.zero-direction"),
        (
            catallax_shapes::RVector3::new(r(3), r(4), r(0)),
            "extrude.zero-thickness",
        ),
    ] {
        assert_eq!(
            extrude(&mut Model::new(), &profile, &step, &r(1), Some(&direction))
                .unwrap_err()
                .code,
            code
        );
    }
}
