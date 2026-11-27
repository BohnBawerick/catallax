//! A certified drawing section cuts the exact solid and views the retained half.

use std::collections::BTreeSet;
use std::num::NonZeroUsize;
use std::time::Instant;

use catallax_app::queries::{SectionPlane, SectionSource, section_solid};
use catallax_app::views::{StrokeSource, hidden_line_view, stroke_drawing};
use catallax_boolean::{BoolOp, boolean};
use catallax_numbers::{Decider, Rational, Sign, set_threads, with_operation_budget};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation, Surface};
use catallax_topology::{BodyId, BodyKind, Loop, Model, Sense, block, cone, cylinder};

fn vector(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}

fn at(x: i64, y: i64, z: i64) -> Pose {
    Pose::from_rational(Rotation::identity(), vector(x, y, z))
}

fn flange(model: &mut Model, dec: &mut Decider) -> BodyId {
    let disk = cylinder(model, 30.into(), 6.into(), &at(0, 0, 0)).unwrap();
    let hub = cylinder(model, 12.into(), 12.into(), &at(0, 0, 4)).unwrap();
    let mut part = boolean(model, disk, hub, BoolOp::Unite, dec).unwrap().body;
    let bore = cylinder(model, 5.into(), 20.into(), &at(0, 0, -2)).unwrap();
    part = boolean(model, part, bore, BoolOp::Subtract, dec)
        .unwrap()
        .body;
    // Rational drawing coordinates on a five-hole bolt circle, including a
    // hole centered exactly on the 144-degree leg.
    for (x, y) in [
        (22_000_000, 0),
        (6_798_374, 20_923_243),
        (-17_798_374, 12_931_276),
        (-17_798_374, -12_931_276),
        (6_798_374, -20_923_243),
    ] {
        let ratio = |v| Rational::from_i64(v) / Rational::from_i64(1_000_000);
        let pose = Pose::from_rational(
            Rotation::identity(),
            RVector3::new(ratio(x), ratio(y), (-2).into()),
        );
        let drill = cylinder(model, 2.into(), 10.into(), &pose).unwrap();
        part = boolean(model, part, drill, BoolOp::Subtract, dec)
            .unwrap()
            .body;
    }
    part
}

fn check_cap(
    cut: &catallax_app::queries::SectionSolid,
    plane: &SectionPlane,
    view: &catallax_app::views::HiddenLineView,
    expected_faces: usize,
) {
    let model = &cut.model;
    let n = plane.normal.to_real();
    let o = plane.origin.to_real();
    let mut faces = 0;
    for face in model.faces_of_body(cut.body) {
        let face = model.face(face);
        let Surface::Plane(support) = model.geom().surface(face.surface) else {
            continue;
        };
        if !support
            .normal
            .cross(&n)
            .norm_squared()
            .rational_ref()
            .is_some_and(Rational::is_zero)
            || !support
                .origin
                .minus(&o)
                .dot(&n)
                .rational_ref()
                .is_some_and(Rational::is_zero)
        {
            continue;
        }
        let facing = support.normal.dot(&n);
        let facing = facing.rational_ref().unwrap().sign();
        assert_eq!(
            facing,
            if face.sense == Sense::Same {
                Sign::Positive
            } else {
                Sign::Negative
            }
        );
        faces += 1;
        for &lp in &face.loops {
            let Loop::Coedges(edges) = model.loop_(lp) else {
                continue;
            };
            for &coedge in edges {
                let edge = model.coedge(coedge).edge;
                assert!(
                    view.strokes
                        .iter()
                        .any(|s| s.source == StrokeSource::Edge(edge) && s.visible),
                    "section edge {edge:?} must be visible from the cut side"
                );
            }
        }
    }
    assert_eq!(faces, expected_faces);
}

#[test]
fn inclined_sections_certify() {
    for (name, normals) in [
        (
            "flange",
            &[
                ("straight", [1, 0, 0], 2),
                ("135", [1, 1, 0], 2),
                ("143.13", [3, 4, 0], 3),
                ("144", [12_931_276, 17_798_374, 0], 3),
            ][..],
        ),
        ("cylinder", &[("ellipse", [2, 3, 1], 1)][..]),
        (
            "frustum",
            &[
                ("ellipse", [2, -3, 2], 1),
                ("hyperbola", [2, -3, 1], 1),
                ("parabola", [10, 0, 3], 1),
            ][..],
        ),
    ] {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let body = match name {
            "flange" => flange(&mut model, &mut dec),
            "cylinder" => cylinder(&mut model, 8.into(), 10.into(), &at(0, 0, -5)).unwrap(),
            "frustum" => cone(&mut model, 8.into(), 5.into(), 10.into(), &at(0, 0, -5)).unwrap(),
            _ => unreachable!(),
        };
        for &(angle, [x, y, z], cap_faces) in normals {
            let plane = SectionPlane {
                origin: RPoint3::origin(),
                normal: vector(x, y, z),
            };
            let mut cancelled = Decider::default();
            cancelled.set_deadline(Some(Instant::now()));
            assert_eq!(
                section_solid(&model, body, &plane, &mut cancelled)
                    .unwrap_err()
                    .refusal
                    .code,
                "numbers.deadline"
            );
            assert_eq!(
                hidden_line_view(&model, body, &plane.normal, &mut cancelled)
                    .unwrap_err()
                    .refusal
                    .code,
                "numbers.deadline"
            );
            let mut first = None;
            for workers in [1, 2, 4, 24] {
                set_threads(NonZeroUsize::new(workers));
                let mut dec = Decider::default();
                let start = Instant::now();
                let (cut, cut_work) = with_operation_budget(u64::MAX, || {
                    section_solid(&model, body, &plane, &mut dec)
                });
                let cut_time = start.elapsed();
                let cut = cut
                    .unwrap_or_else(|e| panic!("{name}/{angle} cut: {e:?}"))
                    .value;
                let start = Instant::now();
                let (view, view_work) = with_operation_budget(u64::MAX, || {
                    hidden_line_view(&cut.model, cut.body, &plane.normal, &mut dec)
                });
                let view_time = start.elapsed();
                eprintln!(
                    "{name}/{angle} {workers} workers: cut {} work {:.3} s; view {} work {:.3} s",
                    cut_work.used,
                    cut_time.as_secs_f64(),
                    view_work.used,
                    view_time.as_secs_f64(),
                );
                let view = view.unwrap_or_else(|e| panic!("{name}/{angle} view: {e:?}"));
                assert!(
                    !cut_work.exhausted && !view_work.exhausted,
                    "{name}/{angle}"
                );
                check_cap(&cut, &plane, &view.value, cap_faces);
                let strokes: Vec<_> = view
                    .value
                    .strokes
                    .iter()
                    .map(|stroke| {
                        let drawn = stroke_drawing(&view.value, stroke, &mut dec).unwrap();
                        assert!(drawn.value.fit.is_none());
                        (
                            stroke.source,
                            stroke.visible,
                            stroke.smooth,
                            drawn.value.curve,
                        )
                    })
                    .collect();
                let mut step = Vec::new();
                catallax_step::write_step(
                    &cut.model,
                    &[cut.body],
                    &Default::default(),
                    &mut dec,
                    &mut step,
                )
                .unwrap_or_else(|e| panic!("{name}/{angle} STEP: {e:?}"));
                let record = (strokes, view.marks, step, cut_work.used, view_work.used);
                if let Some(base) = &first {
                    assert_eq!(base, &record, "{name}/{angle} with {workers} workers");
                } else {
                    first = Some(record);
                }
            }
        }
    }
    set_threads(None);
}

#[test]
fn inclined_flange_view_retains_original_ids_and_cap_ids() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let body = flange(&mut model, &mut dec);
    let input_faces: BTreeSet<_> = model.faces_of_body(body).into_iter().collect();
    let input_edges: BTreeSet<_> = model.edges_of_body(body).into_iter().collect();

    for (angle, normal, cap_count) in [
        (
            "144",
            vector(20_572_483_830_236_564, 28_315_594_803_123_158, 0),
            3,
        ),
        ("135", vector(1, 1, 0), 2),
    ] {
        let plane = SectionPlane {
            origin: RPoint3::new(0.into(), 0.into(), 10.into()),
            normal,
        };
        let cut = section_solid(&model, body, &plane, &mut dec)
            .unwrap_or_else(|e| panic!("{angle} section: {e}"))
            .value;
        let view = hidden_line_view(&cut.model, cut.body, &plane.normal, &mut dec)
            .unwrap_or_else(|e| panic!("{angle} view: {e}"));
        assert!(view.marks.is_empty(), "{angle}: {:?}", view.marks);
        check_cap(&cut, &plane, &view.value, cap_count);

        let faces = cut.model.faces_of_body(cut.body);
        let edges = cut.model.edges_of_body(cut.body);
        assert_eq!(cut.face_sources.len(), faces.len());
        assert_eq!(cut.edge_sources.len(), edges.len());
        let mut caps = BTreeSet::new();
        let mut cap_edges = BTreeSet::new();
        for face in faces {
            let sources = &cut.face_sources[&face];
            assert!(!sources.is_empty(), "{angle} face {face:?}");
            for &source in sources {
                match source {
                    SectionSource::Original(original) => {
                        assert!(input_faces.contains(&original));
                        assert_eq!(cut.model.face(face).surface, model.face(original).surface);
                    }
                    SectionSource::Cap(cap) => {
                        assert_eq!(face, cap);
                        caps.insert(cap);
                        for &lp in &cut.model.face(cap).loops {
                            if let Loop::Coedges(coedges) = cut.model.loop_(lp) {
                                cap_edges.extend(
                                    coedges.iter().map(|&c| (cut.model.coedge(c).edge, cap)),
                                );
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(caps.len(), cap_count);
        let mut inherited = false;
        let mut generated = false;
        for edge in edges {
            let sources = &cut.edge_sources[&edge];
            assert!(!sources.is_empty(), "{angle} edge {edge:?}");
            generated |= sources.iter().all(|s| matches!(s, SectionSource::Cap(_)));
            for &source in sources {
                match source {
                    SectionSource::Original(original) => {
                        inherited = true;
                        assert!(input_edges.contains(&original));
                        assert!(
                            catallax_intersect::same_curve(
                                cut.model.geom(),
                                cut.model.edge(edge).curve,
                                model.edge(original).curve,
                                &mut dec,
                            )
                            .is_coincident()
                        );
                    }
                    SectionSource::Cap(cap) => assert!(cap_edges.remove(&(edge, cap))),
                }
            }
        }
        assert!(inherited && generated && cap_edges.is_empty());
    }
}

#[test]
fn section_outside_body_reports_only_retained_topology() {
    let mut model = Model::new();
    let mut dec = Decider::default();
    let body = block(&mut model, [4.into(), 4.into(), 4.into()], &at(0, 0, 0)).unwrap();
    let mut plane = SectionPlane {
        origin: RPoint3::new((-1).into(), 0.into(), 0.into()),
        normal: vector(1, 1, 0),
    };
    let removed = section_solid(&model, body, &plane, &mut dec).unwrap().value;
    assert_eq!(removed.model.body(removed.body).kind, BodyKind::Empty);
    assert!(removed.face_sources.is_empty() && removed.edge_sources.is_empty());

    plane.origin.x = 9.into();
    let retained = section_solid(&model, body, &plane, &mut dec).unwrap().value;
    assert_eq!(retained.model.body(retained.body).kind, BodyKind::Solid);
    let mut faces: BTreeSet<_> = model.faces_of_body(body).into_iter().collect();
    for face in retained.model.faces_of_body(retained.body) {
        let [SectionSource::Original(original)] = retained.face_sources[&face].as_slice() else {
            panic!("uncut face {face:?} must have one original");
        };
        assert!(faces.remove(original));
    }
    let mut edges: BTreeSet<_> = model.edges_of_body(body).into_iter().collect();
    for edge in retained.model.edges_of_body(retained.body) {
        let [SectionSource::Original(original)] = retained.edge_sources[&edge].as_slice() else {
            panic!("uncut edge {edge:?} must have one original");
        };
        assert!(edges.remove(original));
    }
    assert!(faces.is_empty() && edges.is_empty());
}
