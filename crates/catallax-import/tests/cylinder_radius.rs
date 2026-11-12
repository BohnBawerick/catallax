//! Radius discrepancies against complete, unchanged circular rims.
use catallax_import::{ImportRefusal, apply, build_file, read_step, rebuild};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::Pose;
use catallax_topology::{Model, cylinder};

fn source(radius: &str, height: &str) -> Vec<u8> {
    source_at(radius, height, &Pose::identity())
}

fn source_at(radius: &str, height: &str, pose: &Pose) -> Vec<u8> {
    let mut model = Model::new();
    let body = cylinder(
        &mut model,
        Rational::from_i64(2),
        Rational::parse(height).unwrap(),
        pose,
    )
    .unwrap();
    let mut bytes = Vec::new();
    catallax_step::write_step(
        &model,
        &[body],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let result = text
        .lines()
        .map(|line| {
            if line.contains("CYLINDRICAL_SURFACE(") {
                format!("{},{radius});", &line[..line.rfind(',').unwrap()])
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert_ne!(result, text.trim_end());
    result.into_bytes()
}

#[test]
fn a_radius_discrepancy_at_the_ten_micrometre_cap_has_a_source_bound() {
    for radius in ["1.99", "2.01"] {
        let bytes = source(radius, "5");
        assert!(build_file(&mut Model::new(), &read_step(&bytes).unwrap(), None).is_err());
        let mut model = Model::new();
        let imported = apply(
            &mut model,
            &bytes,
            "radius.step",
            &Default::default(),
            &[],
            None,
        )
        .unwrap();
        assert!(model.check(imported.bodies[0].1).is_empty());
        let tolerance = &imported.record.face_tolerances;
        assert_eq!(tolerance.len(), 1);
        assert_eq!(tolerance[0].deviation_lower_mm, "0.01");
        assert_eq!(tolerance[0].deviation_upper_mm, "0.01");
        assert_eq!(tolerance[0].approved_tolerance_mm, "0.01");
        let rebuilt = rebuild(
            &mut Model::new(),
            &bytes,
            &Default::default(),
            &imported.record,
            None,
        )
        .unwrap();
        assert_eq!(rebuilt.record, imported.record);
        assert_eq!(rebuilt.derivative, imported.derivative);
        assert!(imported.derivative.is_some());
    }
}

#[test]
fn an_over_cap_or_thin_cylinder_keeps_a_measured_refusal_and_no_body() {
    for (radius, height) in [("1.989999", "5"), ("1.994", "0.05")] {
        let mut model = Model::new();
        let error = apply(
            &mut model,
            &source(radius, height),
            "radius.step",
            &Default::default(),
            &[],
            None,
        )
        .unwrap_err();
        assert!(
            matches!(error, ImportRefusal::ValidationFailed { ref detail, .. } if detail.contains("radius discrepancy") && detail.contains("mm")),
            "{error:?}"
        );
        let body = cylinder(
            &mut model,
            Rational::one(),
            Rational::one(),
            &Pose::identity(),
        )
        .unwrap();
        assert_eq!(body.index(), 0);
    }
}

#[test]
fn a_forged_cylinder_source_interpretation_rolls_back() {
    let bytes = source("1.996", "5");
    let mut record = apply(
        &mut Model::new(),
        &bytes,
        "radius.step",
        &Default::default(),
        &[],
        None,
    )
    .unwrap()
    .record;
    assert_eq!(record.format, catallax_import::CYLINDER_REPAIR_FORMAT);
    record.cylinder_radius[0].source_definition = "arbitrary_patch".into();
    let mut model = Model::new();
    let before = format!("{model:?}");
    assert!(rebuild(&mut model, &bytes, &Default::default(), &record, None).is_err());
    assert_eq!(format!("{model:?}"), before);
}

#[test]
fn posed_radius_repairs_are_identical_at_one_two_and_four_workers() {
    use catallax_shapes::{RVector3, Rotation};
    let pose = Pose::from_rational(
        Rotation::from_quaternion(2.into(), 1.into(), 0.into(), 0.into()).unwrap(),
        RVector3::new(7.into(), (-3).into(), 11.into()),
    );
    let bytes = source_at("1.996", "5", &pose);
    let run = || {
        let mut model = Model::new();
        let imported = apply(
            &mut model,
            &bytes,
            "radius.step",
            &Default::default(),
            &[],
            None,
        )
        .unwrap();
        assert_eq!(
            imported.record.face_tolerances[0].deviation_upper_mm,
            "0.004"
        );
        let mut step = Vec::new();
        let mut dec = Decider::default();
        catallax_step::write_step(
            &model,
            &[imported.bodies[0].1],
            &Default::default(),
            &mut dec,
            &mut step,
        )
        .unwrap();
        (
            serde_json::to_string(&imported.record).unwrap(),
            imported.derivative,
            step,
            format!("{:?}", dec.stats()),
        )
    };
    let expected = run();
    for kernel_threads in [1, 2, 4] {
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(kernel_threads));
        for workers in [1, 2, 4] {
            std::thread::scope(|scope| {
                let threads: Vec<_> = (0..workers).map(|_| scope.spawn(run)).collect();
                for thread in threads {
                    assert_eq!(thread.join().unwrap(), expected);
                }
            });
        }
    }
    catallax_numbers::set_threads(None);
}

#[test]
fn iges_radius_repairs_replay_through_the_same_source_format() {
    use catallax_import::{SourceFormat, apply_source, rebuild_source};
    let mut model = Model::new();
    let body = cylinder(
        &mut model,
        Rational::from_i64(2),
        Rational::from_i64(5),
        &Pose::identity(),
    )
    .unwrap();
    let mut bytes = Vec::new();
    catallax_iges::write_iges(
        &model,
        &[body],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let mut changed = 0;
    let bytes = text
        .lines()
        .map(|line| {
            if line.as_bytes().get(72) == Some(&b'P') && line[..64].starts_with("192,") {
                let mut parameters = line[..64].trim_end().split(',').collect::<Vec<_>>();
                parameters[3] = "1.996";
                let data = parameters.join(",");
                assert!(data.len() <= 64);
                changed += 1;
                format!("{data:<64}{}", &line[64..])
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
        .into_bytes();
    assert_eq!(changed, 1);
    let imported = apply_source(
        &mut Model::new(),
        SourceFormat::Iges,
        &bytes,
        "radius.igs",
        &Default::default(),
        &[],
        None,
    )
    .unwrap();
    assert_eq!(
        imported.record.face_tolerances[0].deviation_upper_mm,
        "0.004"
    );
    let again = rebuild_source(
        &mut Model::new(),
        SourceFormat::Iges,
        &bytes,
        &Default::default(),
        &imported.record,
        None,
    )
    .unwrap();
    assert_eq!(again.record, imported.record);
    assert_eq!(again.derivative, imported.derivative);
}

#[test]
fn binary_noise_keeps_its_established_record_before_radius_refitting() {
    let bytes = source("1.9999999999999996", "5");
    let imported = apply(
        &mut Model::new(),
        &bytes,
        "noise.step",
        &Default::default(),
        &[],
        None,
    )
    .unwrap();
    assert_eq!(imported.record.format, catallax_import::REPAIR_FORMAT);
    assert!(imported.record.cylinder_radius.is_empty());
    assert_eq!(imported.record.face_tolerances.len(), 3);
    assert!(
        !serde_json::to_string(&imported.record)
            .unwrap()
            .contains("cylinder_radius")
    );
}
