use super::*;
use catallax_topology::{Coedge, Edge, Face, Loop};

fn face(upper: bool, sense: Sense) -> (Model, FaceId) {
    let mut model = Model::new();
    let axis = Vector3::new(
        0.into(),
        if upper { (-1).into() } else { 1.into() },
        0.into(),
    );
    let apex = Point3::new(
        0.into(),
        if upper { 23.into() } else { (-11).into() },
        0.into(),
    );
    let reference = Vector3::new(1.into(), 0.into(), 0.into());
    let cone = Cone {
        apex: apex.clone(),
        axis: axis.clone(),
        slope: 1.into(),
        ref_dir: reference.clone(),
    };
    let mut loops = Vec::new();
    for (i, h) in [11i64, 12].into_iter().enumerate() {
        let circle = Circle {
            center: apex.plus(&axis.scaled(&Real::from(h))),
            normal: axis.clone(),
            ref_dir: reference.clone(),
            radius_sq: Real::from(h * h),
        };
        let curve = model.geom_mut().add_curve(Curve::Circle(circle));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        let coedge = model.add_coedge(Coedge {
            edge,
            sense: if (i == 0) == (sense == Sense::Same) {
                Sense::Same
            } else {
                Sense::Reversed
            },
        });
        loops.push(model.add_loop(Loop::Coedges(vec![coedge])));
    }
    let surface = model.geom_mut().add_surface(Surface::Cone(cone));
    let face = model.add_face(Face {
        surface,
        sense,
        loops,
    });
    (model, face)
}
fn tool(offset: &str, skew: bool) -> Cylinder {
    Cylinder {
        origin: Point3::new(Rational::parse(offset).unwrap().into(), 12.into(), 0.into()),
        axis: if skew {
            Vector3::new(1.into(), (-1).into(), 1.into())
        } else {
            Vector3::new(0.into(), (-1).into(), 0.into())
        },
        radius_sq: 1.into(),
        ref_dir: Vector3::new(1.into(), 0.into(), 0.into()),
    }
}
fn run(upper: bool, sense: Sense, offset: &str, skew: bool) -> NumericalFacePartition {
    let (model, face) = face(upper, sense);
    let (answer, work) = catallax_numbers::with_operation_budget(20_000_000, || {
        clip_cone_cylinder_face(
            &model,
            face,
            "cone-test",
            &tool(offset, skew),
            if offset == "10" {
                [q(-49), q(49)]
            } else {
                [
                    Rational::from_ratio(-99, 2).unwrap(),
                    Rational::from_ratio(99, 2).unwrap(),
                ]
            },
            NumericalOptions {
                tolerance_mm: Rational::from_ratio(1, 10_000_000).unwrap(),
            },
            &mut Decider::default(),
        )
    });
    let answer = answer.unwrap().unwrap();
    assert!(!work.exhausted, "{work:?}");
    let stages = answer.work;
    let pair = stages.intersection;
    assert_eq!(
        stages.boundaries
            + stages.supports
            + stages.topology
            + pair.windows
            + pair.setup
            + pair.tracing
            + pair.fit_checks,
        work.used
    );
    for event in &answer.events {
        if let Some(vertex) = event.vertex {
            assert!(vertex < answer.definition.vertices.len())
        }
    }
    for link in &answer.ring_links {
        assert!(link.edge < answer.definition.edges.len());
        assert!(link.ring < 2)
    }
    answer
}
#[test]
fn separated_and_point_only_faces_are_supported_without_short_edges() {
    let separated = run(false, Sense::Same, "9.999", false);
    assert_eq!(separated.definition.faces.len(), 4);
    assert_eq!(separated.definition.edges.len(), 12);
    assert_eq!(separated.definition.vertices.len(), 8);
    assert!(separated.intersection.pieces().is_empty());
    assert!(separated.events.is_empty());
    let contact = run(false, Sense::Same, "10", false);
    assert_eq!(contact.definition.faces.len(), 4);
    assert_eq!(contact.definition.edges.len(), 12);
    assert_eq!(contact.definition.vertices.len(), 8);
    assert_eq!(contact.intersection.pieces().len(), 4);
    assert_eq!(contact.events.len(), 1);
    assert!(contact.events[0].point_only);
    assert!(contact.events[0].vertex.is_some());
    assert_eq!(contact.events[0].hits.len(), 4);
    assert_eq!(contact.events[0].key, "cone-test/corner-0-0/meeting");
    assert!(contact.cylinder_trims.is_empty());
    assert!(contact.native_trims.is_empty());
}
#[test]
fn axial_crossing_keeps_both_native_arcs_and_original_ring_links() {
    for sense in [Sense::Same, Sense::Reversed] {
        let answer = run(false, sense, "10.5", false);
        assert_eq!(answer.definition.faces.len(), 4);
        assert_eq!(answer.definition.edges.len(), 14);
        assert_eq!(answer.definition.vertices.len(), 10);
        assert_eq!(answer.native_trims.len(), 2);
        assert_eq!(answer.cylinder_trims.len(), 2);
        assert_eq!(answer.events.len(), 3);
        assert_eq!(answer.ring_links.len(), 8);
        for binding in &answer.native_trims {
            let original = &answer.intersection.pieces()[answer
                .cylinder_trims
                .iter()
                .find(|t| answer.definition.edges[t.edge].key == binding.edge)
                .unwrap()
                .piece];
            assert!(binding.curve.same_definition(&original.uv()[0]));
        }
        assert!(answer.definition.faces.iter().all(|f| f.sense == sense));
    }
}
#[test]
fn skew_retains_the_unchanged_lower_face_and_all_upper_material_regions() {
    let lower = run(false, Sense::Same, "10.5", true);
    assert!(lower.intersection.pieces().is_empty());
    assert_eq!(lower.definition.faces.len(), 4);
    let upper = run(true, Sense::Same, "10.5", true);
    assert_eq!(upper.intersection.pieces().len(), 3);
    assert_eq!(upper.definition.faces.len(), 5);
    assert_eq!(upper.cylinder_trims.len(), 3);
    assert_eq!(upper.native_trims.len(), 3);
    assert_eq!(upper.events.len(), 5);
    for i in 0..4 {
        let Surface::BSpline(SplineSupport::RationalPatch(a)) =
            &lower.cylinder_source.quarters[i].surface
        else {
            unreachable!()
        };
        let Surface::BSpline(SplineSupport::RationalPatch(b)) =
            &upper.cylinder_source.quarters[i].surface
        else {
            unreachable!()
        };
        assert!(a.carrier().same_definition(b.carrier()));
        assert_eq!(upper.cylinder_source.quarters[i].id, i);
    }
}
