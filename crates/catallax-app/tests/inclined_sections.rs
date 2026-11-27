//! A certified drawing section cuts the exact solid and views the retained half.

use std::num::NonZeroUsize;
use std::time::Instant;

use catallax_app::views::hidden_line_view;
use catallax_boolean::{BoolOp, boolean};
use catallax_numbers::{Decider, Rational, set_threads, with_operation_budget};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model, block, cone, cylinder};

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

fn section_box(model: &mut Model, normal: &RVector3, dec: &mut Decider) -> BodyId {
    let rotation =
        Rotation::from_frame(&normal.to_real(), &vector(0, 0, 1).to_real(), dec).unwrap();
    let translation = rotation.apply_real(&vector(-100, -100, 0).to_real());
    let pose = Pose::new(rotation, translation, dec).unwrap();
    block(model, [200.into(), 200.into(), 100.into()], &pose).unwrap()
}

#[test]
fn inclined_sections_certify() {
    for (name, normals) in [
        (
            "flange",
            &[
                ("straight", [1, 0, 0]),
                ("135", [1, 1, 0]),
                ("143.13", [3, 4, 0]),
                ("144", [12_931_276, 17_798_374, 0]),
            ][..],
        ),
        ("cylinder", &[("compound", [2, 3, 1])][..]),
        ("frustum", &[("compound", [2, -3, 1])][..]),
    ] {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let body = match name {
            "flange" => flange(&mut model, &mut dec),
            "cylinder" => cylinder(&mut model, 8.into(), 10.into(), &at(0, 0, -5)).unwrap(),
            "frustum" => cone(&mut model, 8.into(), 5.into(), 10.into(), &at(0, 0, -5)).unwrap(),
            _ => unreachable!(),
        };
        for &(angle, [x, y, z]) in normals {
            set_threads(NonZeroUsize::new(1));
            let normal = vector(x, y, z);
            let mut model = model.clone();
            let tool = section_box(&mut model, &normal, &mut dec);
            let start = Instant::now();
            let (cut, work) = with_operation_budget(100_000_000, || {
                boolean(&mut model, body, tool, BoolOp::Subtract, &mut dec)
            });
            eprintln!(
                "{name}/{angle}: cut {} work, {:?}",
                work.used,
                start.elapsed()
            );
            let cut = match cut {
                Ok(answer) => answer.body,
                Err(err) => {
                    eprintln!("{name}/{angle}: cut refused: {err:?}");
                    continue;
                }
            };
            let start = Instant::now();
            let (view, work) = with_operation_budget(100_000_000, || {
                hidden_line_view(&model, cut, &normal, &mut dec)
            });
            eprintln!(
                "{name}/{angle}: view {} work, {:?}",
                work.used,
                start.elapsed()
            );
            match view {
                Ok(answer) => assert!(answer.value.strokes.iter().any(|stroke| stroke.visible)),
                Err(err) => eprintln!("{name}/{angle}: view refused: {err:?}"),
            }
        }
    }
    set_threads(None);
}
