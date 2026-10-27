//! The final lattice stage: request building, refusals and rebuild from an edited body.
use catallax_field::{
    Location, Workspace,
    lattice_stage::{LatticeSpec, Settings, Stage, check_no_cad_after_lattice},
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Grade, LatticeCell, Number, Pose, Recipe};
use catallax_shapes::RPoint3;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn body(id: &str, size: &str) -> Recipe {
    Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/3","name":"body","steps":[{{"id":"{id}","op":"block","size":["{size}","4","4"],"pose":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}}}],"result":"{id}"}}"#
    ))
    .unwrap()
}
fn spec() -> LatticeSpec {
    LatticeSpec {
        cell: LatticeCell::Cubic,
        size: Number::parse("2").unwrap(),
        frame: serde_json::from_str::<Pose>(
            r#"{"rotation":["1","0","0","0"],"translation":["0","0","0"]}"#,
        )
        .unwrap(),
        radius: serde_json::from_str::<Grade>(r#"{"constant":"0.25"}"#).unwrap(),
    }
}
fn settings() -> Settings {
    Settings {
        accuracy: q("0.25"),
        minimum_detail: q("0.5"),
        max_cells: 2_000_000,
        max_vertices: 1_000_000,
    }
}
fn location(stage: &Stage, x: &str, y: &str, z: &str) -> Location {
    let point = RPoint3::new(q(x), q(y), q(z));
    stage
        .plan()
        .proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
        .unwrap()
        .location
}

#[test]
fn the_stage_makes_the_lattice_the_result_of_the_body() {
    let stage = Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(stage.lattice_step().as_str(), "lattice");
    assert_eq!(stage.recipe().result, *stage.lattice_step());
    assert_eq!(location(&stage, "2", "1", "2"), Location::Inside);
    assert_eq!(location(&stage, "1", "1", "1"), Location::Outside);
}

fn first_again() -> Stage {
    Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn an_edited_body_rebuilds_the_lattice_from_the_edited_body() {
    let first = Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    let (earlier, kept_spec, kept_settings) = first.return_before_lattice();
    assert_eq!(kept_settings, settings());
    assert_eq!(kept_spec, spec());
    assert_eq!(earlier, body("block", "4"));
    assert_eq!(
        location(&first_again(), "3", "0.125", "2"),
        Location::Inside
    );
    // The block is shortened along x; the same lattice parameters now stop at x = 2.
    let edited = Stage::new(
        &body("block", "2"),
        kept_spec,
        kept_settings,
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(location(&edited, "3", "0.125", "2"), Location::Outside);
    assert_eq!(location(&edited, "1", "0.125", "2"), Location::Inside);
}

#[test]
fn an_existing_lattice_in_the_body_is_cad_after_lattice() {
    let first = Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    let error =
        Stage::new(first.recipe(), spec(), settings(), &mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "lattice.cad-after-lattice");
}

#[test]
fn a_step_that_reads_a_lattice_is_cad_after_lattice() {
    let text = r#"{"format":"catallax-recipe/3","name":"x","steps":[{"id":"region","op":"block","size":["4","4","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"radius":{"constant":"0.25"}},{"id":"after","op":"offset","body":"beams","distance":"0.1"}],"result":"after"}"#;
    let recipe = Recipe::parse(text).unwrap();
    assert_eq!(
        check_no_cad_after_lattice(&recipe).unwrap_err().code,
        "lattice.cad-after-lattice"
    );
    assert!(check_no_cad_after_lattice(&body("block", "4")).is_ok());
}

#[test]
fn impossible_settings_and_lattices_refuse_with_their_own_codes() {
    let mut dec = Decider::default();
    for bad in [
        Settings {
            accuracy: q("0"),
            ..settings()
        },
        Settings {
            minimum_detail: q("-1"),
            ..settings()
        },
        Settings {
            max_cells: 0,
            ..settings()
        },
        Settings {
            max_vertices: 3,
            ..settings()
        },
    ] {
        let error = Stage::new(&body("block", "4"), spec(), bad, &mut dec).unwrap_err();
        assert_eq!(error.code, "lattice.settings");
    }
    let thin = LatticeSpec {
        radius: serde_json::from_str::<Grade>(r#"{"constant":"0"}"#).unwrap(),
        ..spec()
    };
    let error = Stage::new(&body("block", "4"), thin, settings(), &mut dec).unwrap_err();
    assert_eq!(error.code, "lattice.impossible");
}

#[test]
fn the_lattice_step_id_avoids_a_taken_name() {
    let stage = Stage::new(
        &body("lattice", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(stage.lattice_step().as_str(), "lattice-1");
}
