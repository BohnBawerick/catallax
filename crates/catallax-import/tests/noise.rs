//! Binary-noise cleanup through the repair entry, for STEP files.

use catallax_import::{ImportRefusal, apply, build_file, read_step};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_step::{StepOptions, write_step};
use catallax_topology::{Model, cylinder};

/// A cylinder of radius 2 and height 5, with its ring vertices moved by `y`.
fn noisy(y: &str) -> Vec<u8> {
    let mut m = Model::new();
    let body = cylinder(
        &mut m,
        Rational::from_i64(2),
        Rational::from_i64(5),
        &Pose::identity(),
    )
    .unwrap();
    let mut bytes = Vec::new();
    write_step(
        &m,
        &[body],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let moved = text
        .replace("(2.,0.,0.)", &format!("(2.,{y},0.)"))
        .replace("(2.,0.,5.)", &format!("(2.,{y},5.)"));
    assert_ne!(moved, text);
    moved.into_bytes()
}

#[test]
fn a_noisy_seam_vertex_imports_with_a_recorded_bound() {
    let bytes = noisy("-7.3478807948841188E-16");
    let mut exact = Model::new();
    assert!(matches!(
        build_file(&mut exact, &read_step(&bytes).unwrap(), None),
        Err(ImportRefusal::NotSolid { .. })
    ));
    let mut m = Model::new();
    let imported = apply(&mut m, &bytes, "noisy.step", &Default::default(), &[], None).unwrap();
    assert_eq!(imported.bodies.len(), 1);
    assert!(imported.derivative.is_some());
    assert_eq!(imported.record.face_tolerances.len(), 3);
    let limit = Rational::from_ratio(1, 1_000_000_000).unwrap();
    for tolerance in &imported.record.face_tolerances {
        assert!(Rational::parse(&tolerance.deviation_upper_mm).unwrap() < limit);
    }
}

#[test]
fn written_offsets_are_not_noise() {
    // A deliberate 1 um offset written with ten digits is not snapped, and a
    // short decimal never is: the exact refusal stands.
    for y in ["0.0000010000000000", "1E-15"] {
        let mut m = Model::new();
        let error = apply(
            &mut m,
            &noisy(y),
            "offset.step",
            &Default::default(),
            &[],
            None,
        )
        .unwrap_err();
        assert!(
            matches!(error, ImportRefusal::NotSolid { .. }),
            "{y}: {error:?}"
        );
    }
}

/// Two cylinders written as one file: the first exact at the origin, the
/// second translated by (10, 0, 0) with its ring vertices moved by `y`.
fn two_noisy(y: &str) -> Vec<u8> {
    let mut m = Model::new();
    let first = cylinder(
        &mut m,
        Rational::from_i64(2),
        Rational::from_i64(5),
        &Pose::identity(),
    )
    .unwrap();
    let second = cylinder(
        &mut m,
        Rational::from_i64(2),
        Rational::from_i64(5),
        &Pose::from_rational(
            Rotation::identity(),
            RVector3::new(Rational::from_i64(10), Rational::zero(), Rational::zero()),
        ),
    )
    .unwrap();
    let mut bytes = Vec::new();
    write_step(
        &m,
        &[first, second],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let moved = text
        .replace("(12.,0.,0.)", &format!("(12.,{y},0.)"))
        .replace("(12.,0.,5.)", &format!("(12.,{y},5.)"));
    assert_ne!(moved, text);
    moved.into_bytes()
}

/// A caller body plus an exact first solid and a noisy second solid: the
/// failed exact attempt may leave nothing behind, and the successful cleaned
/// retry commits exactly the returned imported bodies.
#[test]
fn a_failed_exact_attempt_leaves_no_partial_geometry() {
    let bytes = two_noisy("-7.3478807948841188E-16");
    let mut m = Model::new();
    let caller = cylinder(
        &mut m,
        Rational::from_i64(2),
        Rational::from_i64(5),
        &Pose::identity(),
    )
    .unwrap();
    let imported = apply(&mut m, &bytes, "two.step", &Default::default(), &[], None).unwrap();
    assert_eq!(imported.bodies.len(), 2);
    assert!(m.check(caller).is_empty());
    for (i, (_, body)) in imported.bodies.iter().enumerate() {
        assert_eq!(body.index(), caller.index() + 1 + i as u32, "{i}");
        assert!(m.check(*body).is_empty());
    }
}

/// An overall refusal after a retry-eligible failure: no solid of the failed
/// exact attempt may remain in the caller's model.
#[test]
fn a_refused_retry_leaves_the_caller_model_unchanged() {
    // The second solid's seam sits one micrometre off, written with ten
    // digits: not noise, so the exact refusal stands.
    let bytes = two_noisy("0.0000010000000000");
    let mut m = Model::new();
    let caller = cylinder(
        &mut m,
        Rational::from_i64(2),
        Rational::from_i64(5),
        &Pose::identity(),
    )
    .unwrap();
    let error = apply(&mut m, &bytes, "two.step", &Default::default(), &[], None).unwrap_err();
    assert!(matches!(error, ImportRefusal::NotSolid { .. }), "{error:?}");
    assert!(m.check(caller).is_empty());
    // The caller's body is still the last one in the model: the failed exact
    // attempt committed no body of its own.
    let probe = cylinder(
        &mut m,
        Rational::from_i64(2),
        Rational::from_i64(5),
        &Pose::identity(),
    )
    .unwrap();
    assert_eq!(probe.index(), 1);
}
