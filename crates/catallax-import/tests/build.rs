//! Topology construction from recognised STEP faces, through the public entries.

mod common;

use catallax_import::{ImportRefusal, build_file, read_step};
use catallax_measure::{Accuracy, volume};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Origin, Role, StepId};
use catallax_shapes::Pose;
use catallax_topology::{BodyId, Model, cylinder};
use common::*;

fn import(text: &str, step: Option<&StepId>) -> Result<(Model, Vec<BodyId>), ImportRefusal> {
    let file = read_step(text.as_bytes())?;
    let mut model = Model::new();
    let bodies = build_file(&mut model, &file, step)?;
    Ok((model, bodies))
}

fn vol(model: &Model, body: BodyId) -> (Rational, Rational) {
    let v = volume(model, body, &Accuracy::default(), &mut Decider::default())
        .expect("volume")
        .value;
    (v.lo, v.hi)
}

#[test]
fn a_clean_box_builds_a_sound_body_with_exact_volume() {
    let (model, bodies) = import(&box_step(BoxFixture::new(["10", "20", "0.5"])), None).unwrap();
    assert!(model.check(bodies[0]).is_empty());
    assert_eq!(model.faces_of_body(bodies[0]).len(), 6);
    assert_eq!(model.edges_of_body(bodies[0]).len(), 12);
    assert_eq!(vol(&model, bodies[0]), (r("100"), r("100")));
}

#[test]
fn inch_files_build_in_millimetres() {
    let mut b = BoxFixture::new(["1", "1", "1"]);
    b.unit = Unit::Inch;
    let (model, bodies) = import(&box_step(b), None).unwrap();
    assert_eq!(vol(&model, bodies[0]), (r("16387.064"), r("16387.064")));
}

#[test]
fn faces_carry_imported_origins_keyed_by_entity_name() {
    let step = StepId::new("imp").unwrap();
    let (model, bodies) = import(&box_step(BoxFixture::new(["1", "1", "1"])), Some(&step)).unwrap();
    for f in model.faces_of_body(bodies[0]) {
        let [origin] = model.face_origins(f) else {
            panic!("one origin")
        };
        assert_eq!(origin.step, step);
        assert!(matches!(&origin.role, Role::Imported { key } if key.starts_with('#')));
        let _: &Origin = origin;
    }
}

#[test]
fn two_boxes_build_two_bodies() {
    let mut b = BoxFixture::new(["1", "1", "1"]);
    b.origin = [r("1.005"), r("0"), r("0")];
    let (model, bodies) =
        import(&boxes_step(&[BoxFixture::new(["1", "1", "1"]), b]), None).unwrap();
    assert_eq!(bodies.len(), 2);
    for body in bodies {
        assert!(model.check(body).is_empty());
    }
}

#[test]
fn a_face_off_by_any_amount_is_not_solid_and_nothing_is_moved() {
    let mut b = BoxFixture::new(["10", "10", "10"]);
    b.shift = Some((1, [r("0.004"), r("0"), r("0")]));
    let refusal = import(&box_step(b), None).unwrap_err();
    let ImportRefusal::NotSolid { detail } = refusal else {
        panic!("expected not_solid, got {refusal:?}");
    };
    assert!(detail.contains("#"), "{detail}");
}

#[test]
fn a_cylinder_written_by_the_exporter_builds_the_same_volume() {
    use catallax_step::{StepOptions, write_step};
    let mut original = Model::new();
    let c = cylinder(
        &mut original,
        Rational::from_i64(3),
        Rational::from_i64(7),
        &Pose::identity(),
    )
    .unwrap();
    let mut bytes = Vec::new();
    write_step(
        &original,
        &[c],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    let (model, bodies) = import(std::str::from_utf8(&bytes).unwrap(), None).unwrap();
    assert!(model.check(bodies[0]).is_empty());
    let (a_lo, a_hi) = vol(&original, c);
    let (b_lo, b_hi) = vol(&model, bodies[0]);
    // 63 pi, written through decimals: the two enclosures must overlap closely.
    assert!(a_lo <= b_hi && b_lo <= a_hi, "enclosures disjoint");
    let width = Rational::parse("1e-6").unwrap();
    assert!(&b_hi - &b_lo < width);
}

#[test]
fn exporter_cones_and_spheres_are_refused_by_entity() {
    use catallax_step::{StepOptions, write_step};
    use catallax_topology::{cone, sphere};
    for (name, expected) in [("cone", "CONICAL_SURFACE"), ("sphere", "VERTEX_LOOP")] {
        let mut m = Model::new();
        let b = if name == "cone" {
            cone(&mut m, 3.into(), 1.into(), 2.into(), &Pose::identity()).unwrap()
        } else {
            sphere(&mut m, 3.into(), &Pose::identity()).unwrap()
        };
        let mut bytes = Vec::new();
        write_step(
            &m,
            &[b],
            &StepOptions::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        let refusal = import(std::str::from_utf8(&bytes).unwrap(), None).unwrap_err();
        let ImportRefusal::Unsupported { entity, type_name } = refusal else {
            panic!("{name}: expected unsupported, got {refusal:?}");
        };
        assert!(
            entity.starts_with('#') && type_name.starts_with(expected),
            "{name}: {type_name}"
        );
    }
}

#[test]
fn the_consumer_entry_checks_the_saved_hash() {
    use catallax_import::{import_step_bodies, sha256_hex};
    let text = box_step(BoxFixture::new(["1", "2", "3"]));
    let step = StepId::new("imp").unwrap();
    let hash = sha256_hex(text.as_bytes());
    assert_eq!(hash.len(), 64);
    let mut model = Model::new();
    let bodies = import_step_bodies(&mut model, text.as_bytes(), &hash, &step).unwrap();
    assert_eq!(vol(&model, bodies[0]), (r("6"), r("6")));
    let mut model = Model::new();
    let refusal = import_step_bodies(
        &mut model,
        text.replace("'box'", "'bos'").as_bytes(),
        &hash,
        &step,
    )
    .unwrap_err();
    assert!(matches!(refusal, ImportRefusal::HashMismatch { expected, .. } if expected == hash));
}

#[test]
fn sha256_matches_the_known_empty_digest() {
    assert_eq!(
        catallax_import::sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
}
