//! A certified drawing section cuts the exact solid and views the retained half.

use std::num::NonZeroUsize;
use std::time::Instant;

use catallax_app::queries::{SectionPlane, section_solid};
use catallax_app::views::{StrokeSource, hidden_line_view, stroke_drawing};
use catallax_boolean::{BoolOp, boolean};
use catallax_numbers::{Decider, Rational, set_threads, with_operation_budget};
use catallax_shapes::{Pose, RPoint3, RVector3, Rotation, Surface};
use catallax_topology::{BodyId, Loop, Model, cone, cylinder};

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
        assert!(
            support
                .normal
                .dot(&n)
                .rational_ref()
                .is_some_and(|q| q > &Rational::zero())
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
