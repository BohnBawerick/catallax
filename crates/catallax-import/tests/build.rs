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
fn placement_projection_and_defaults_preserve_exported_box_bytes() {
    let export = |text: &str| {
        let (model, bodies) = import(text, None).unwrap();
        assert!(model.check(bodies[0]).is_empty());
        assert_eq!(vol(&model, bodies[0]), (r("6000"), r("6000")));
        let mut bytes = Vec::new();
        catallax_step::write_step(
            &model,
            &bodies,
            &Default::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        bytes
    };
    let original = box_step(BoxFixture::new(["10", "20", "30"]));
    let expected_file = read_step(original.as_bytes()).unwrap();
    let expected_bytes = export(&original);
    for text in [
        box_with_placement(0, Some(["1", "0", "0"]), Some(["2", "1", "0"])),
        box_with_placement(0, Some(["1", "0", "0"]), None),
        box_with_placement(4, None, Some(["1", "0", "2"])),
        box_with_placement(4, None, None),
    ] {
        assert_eq!(read_step(text.as_bytes()).unwrap(), expected_file);
        assert_eq!(export(&text), expected_bytes);
    }
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

fn exported(build: impl FnOnce(&mut Model) -> BodyId) -> (Model, BodyId, String) {
    use catallax_step::{StepOptions, write_step};
    let mut original = Model::new();
    let body = build(&mut original);
    let mut bytes = Vec::new();
    write_step(
        &original,
        &[body],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap();
    (original, body, String::from_utf8(bytes).unwrap())
}

fn pi_times(factor: &str) -> (Rational, Rational) {
    let pi = catallax_numbers::Ball::pi(128);
    let f = r(factor);
    (&pi.lo() * &f, &pi.hi() * &f)
}

fn pi2_times(factor: &str) -> (Rational, Rational) {
    let pi = catallax_numbers::Ball::pi(128);
    let f = r(factor);
    (&(&pi.lo() * &pi.lo()) * &f, &(&pi.hi() * &pi.hi()) * &f)
}

fn assert_volume(text: &str, expected: (Rational, Rational), faces: usize) {
    let (model, bodies) = import(text, None).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(model.check(bodies[0]).is_empty());
    assert_eq!(model.faces_of_body(bodies[0]).len(), faces);
    let (lo, hi) = vol(&model, bodies[0]);
    assert!(lo <= expected.1 && expected.0 <= hi, "{lo} {hi}");
    assert!(&hi - &lo < r("1e-6"));
}

fn assert_near(text: &str, expected: (Rational, Rational), faces: usize, slack: &str) {
    let (model, bodies) = import(text, None).unwrap_or_else(|e| panic!("{e:?}"));
    assert!(model.check(bodies[0]).is_empty());
    assert_eq!(model.faces_of_body(bodies[0]).len(), faces);
    let (lo, hi) = vol(&model, bodies[0]);
    let slack = r(slack);
    assert!(
        &lo - &slack <= expected.1 && expected.0 <= &hi + &slack,
        "{lo} {hi}"
    );
}

fn cone_tolerances(text: &str) -> Vec<catallax_import::FaceTolerance> {
    let imported = catallax_import::apply(
        &mut Model::new(),
        text.as_bytes(),
        "part.step",
        &Default::default(),
        &[],
        None,
    )
    .unwrap_or_else(|e| panic!("{e:?}"));
    imported.record.face_tolerances
}

#[test]
fn exporter_cones_build_through_their_rims_with_a_recorded_tolerance() {
    use catallax_topology::cone;
    // The exporter writes the semi-angle in radians, rounded, and the rims as
    // decimals. Neither lies exactly on the other; the rims decide.
    for (top, faces, factor) in [("1", 3, "26/3"), ("0", 2, "6")] {
        let (_, _, text) =
            exported(|m| cone(m, 3.into(), r(top), 2.into(), &Pose::identity()).unwrap());
        assert!(text.contains("CONICAL_SURFACE"));
        assert_near(&text, pi_times(factor), faces, "1e-6");
        let tolerances = cone_tolerances(&text);
        assert_eq!(tolerances.len(), 1, "{tolerances:?}");
        let t = &tolerances[0];
        assert_eq!(t.approved_tolerance_mm, "0.01");
        assert!(r(&t.deviation_lower_mm) <= r(&t.deviation_upper_mm));
        assert!(r(&t.deviation_upper_mm) <= r("1e-9"), "{t:?}");
    }
}

#[test]
fn cones_at_any_degree_angle_build_with_a_frustum_or_an_apex() {
    // 30 degrees: tan = 0.5773502691896...; rims at 1 and 1 + 3 tan 30.
    let mut soup = Soup::new();
    soup.cone(["0", "0", "0"], "1", "30.", Ok(("3", "2.732050808")));
    let text = in_degrees(&soup.finish());
    // pi h (a^2 + ab + b^2) / 3 for the written rims.
    let (a, b) = (r("1"), r("2.732050808"));
    let factor = (&(&(&a * &a) + &(&a * &b)) + &(&b * &b)).to_string();
    assert_near(&text, pi_times(&factor), 3, "1e-6");
    let tolerances = cone_tolerances(&text);
    assert_eq!(tolerances.len(), 1);
    // The saved record rebuilds, cone tolerance included.
    let imported = catallax_import::apply(
        &mut Model::new(),
        text.as_bytes(),
        "part.step",
        &Default::default(),
        &[],
        None,
    )
    .unwrap();
    let rebuilt = catallax_import::rebuild(
        &mut Model::new(),
        text.as_bytes(),
        &Default::default(),
        &imported.record,
        None,
    )
    .unwrap();
    assert_eq!(rebuilt.record, imported.record);
    let mut tampered = imported.record.clone();
    tampered.face_tolerances[0].deviation_upper_mm = "0.000000001".into();
    assert!(
        catallax_import::rebuild(
            &mut Model::new(),
            text.as_bytes(),
            &Default::default(),
            &tampered,
            None
        )
        .is_err()
    );
    // 1 + 3 tan 30 = 2.7320508075688772...: the top rim misses by about 4.3e-10.
    assert!(r(&tolerances[0].deviation_upper_mm) < r("5e-10"));
    assert!(r(&tolerances[0].deviation_lower_mm) > r("3e-10"));
    // 60 degrees to an apex: the apex is the vertex loop, 1/tan 60 below the base.
    let mut soup = Soup::new();
    soup.cone(["0", "0", "0"], "1", "60.", Err("-0.577350269"));
    let text = in_degrees(&soup.finish());
    assert_near(&text, pi_times("0.192450089666666666"), 2, "1e-6");
}

#[test]
fn cones_beyond_the_automatic_cap_are_refused_with_their_numbers() {
    let mut soup = Soup::new();
    soup.cone(["0", "0", "0"], "1", "30.", Ok(("3", "2.75")));
    let text = in_degrees(&soup.finish());
    for refusal in [
        import(&text, None).unwrap_err(),
        catallax_import::apply(
            &mut Model::new(),
            text.as_bytes(),
            "part.step",
            &Default::default(),
            &[],
            None,
        )
        .unwrap_err(),
    ] {
        let ImportRefusal::NotSolid { detail } = refusal else {
            panic!("expected not_solid, got {refusal:?}");
        };
        assert!(detail.contains("0.01 mm automatic cap"), "{detail}");
        assert!(detail.contains("0.0179"), "{detail}");
    }
}

#[test]
fn exporter_spheres_and_tori_build_the_same_volume() {
    use catallax_topology::{sphere, torus};
    // The exporter writes a whole sphere as one vertex loop at a pole and a
    // whole torus as two seams, each used once each way.
    let (_, _, text) = exported(|m| sphere(m, 3.into(), &Pose::identity()).unwrap());
    assert!(text.contains("VERTEX_LOOP"));
    assert_volume(&text, pi_times("36"), 1);
    let (_, _, text) = exported(|m| torus(m, 4.into(), 1.into(), &Pose::identity()).unwrap());
    assert_volume(&text, pi2_times("8"), 1);
}

#[test]
fn whole_spheres_from_a_vertex_loop_or_a_seam_build() {
    for seam in [false, true] {
        let mut soup = Soup::new();
        if seam {
            soup.seam_sphere(["1", "2", "3"], "2");
        } else {
            soup.vertex_loop_sphere(["1", "2", "3"], "2");
        }
        assert_volume(&soup.finish(), pi_times("32/3"), 1);
    }
}

#[test]
fn hemispheres_and_half_tori_check_their_circles_exactly() {
    let mut soup = Soup::new();
    soup.hemisphere(["0", "0", "0"], 2, "3");
    assert_volume(&soup.finish(), pi_times("18"), 2);
    let mut soup = Soup::new();
    soup.half_torus(["0", "0", "1"], 0, "4", "1");
    assert_volume(&soup.finish(), pi2_times("4"), 2);
    let mut soup = Soup::new();
    soup.meridian_half_torus(["0", "0", "0"], "5", "2");
    assert_volume(&soup.finish(), pi2_times("20"), 3);
}

/// Replaces `from` by `to` on every line that holds `entity`.
fn edit(text: &str, entity: &str, from: &str, to: &str) -> String {
    let out: Vec<String> = text
        .lines()
        .map(|line| {
            if line.contains(entity) {
                line.replace(from, to)
            } else {
                line.to_string()
            }
        })
        .collect();
    let out = out.join("\n");
    assert_ne!(out, text.trim_end());
    out
}

#[test]
fn spherical_and_toroidal_faces_off_their_circles_are_refused() {
    let mut soup = Soup::new();
    soup.hemisphere(["0", "0", "0"], 2, "3");
    let text = edit(&soup.finish(), "SPHERICAL_SURFACE", ",3.0)", ",3.001)");
    let refusal = import(&text, None).unwrap_err();
    assert!(
        matches!(&refusal, ImportRefusal::NotSolid { detail } if detail.contains("does not lie exactly")),
        "{refusal:?}"
    );
    // A torus whose circles miss it is refused too, though the self-intersection
    // check, which runs first, may already refuse it as not verified.
    let mut soup = Soup::new();
    soup.meridian_half_torus(["0", "0", "0"], "5", "2");
    let text = edit(&soup.finish(), "TOROIDAL_SURFACE", ",2.0)", ",2.001)");
    assert!(import(&text, None).is_err());
}

#[test]
fn vertex_loops_off_the_surface_or_beside_other_bounds_are_refused() {
    let mut soup = Soup::new();
    soup.vertex_loop_sphere(["0", "0", "0"], "2");
    let text = edit(&soup.finish(), "CARTESIAN_POINT", "-2.0)", "-2.001)");
    let refusal = import(&text, None).unwrap_err();
    assert!(
        matches!(&refusal, ImportRefusal::NotSolid { detail } if detail.contains("vertex loop")),
        "{refusal:?}"
    );
    // The hemisphere's spherical face is written first; give it a second bound,
    // a vertex loop at its north pole.
    let mut soup = Soup::new();
    soup.hemisphere(["0", "0", "0"], 2, "3");
    let mut text = String::new();
    let mut first = true;
    for line in soup.finish().lines() {
        if line == "ENDSEC;" && text.contains("DATA;") {
            text.push_str("#9001=CARTESIAN_POINT('',(0.0,0.0,3.0));\n#9002=VERTEX_POINT('',#9001);\n#9003=VERTEX_LOOP('',#9002);\n#9004=FACE_BOUND('',#9003,.T.);\n");
        }
        if first && line.contains("ADVANCED_FACE") {
            first = false;
            text.push_str(&line.replacen("),#", ",#9004),#", 1));
        } else {
            text.push_str(line);
        }
        text.push('\n');
    }
    let refusal = import(&text, None).unwrap_err();
    assert!(
        matches!(&refusal, ImportRefusal::Unsupported { type_name, .. } if type_name.starts_with("VERTEX_LOOP")),
        "{refusal:?}"
    );
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
