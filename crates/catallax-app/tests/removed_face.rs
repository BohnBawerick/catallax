//! A Boolean-removed face cannot supply a later extrusion's support.
use catallax_numbers::Decider;
use catallax_recipe::Recipe;
use serde_json::json;

fn recipe(cut_x: &str, width: &str) -> Recipe {
    Recipe::parse(&json!({"format":"catallax-recipe/5","name":"removed support","steps":[
        {"id":"base","op":"block","size":["10","10","10"],
            "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"upper","op":"block","size":[width,"10","5"],
            "pose":{"rotation":["1","0","0","0"],"translation":[cut_x,"0","5"]}},
        {"id":"lower","op":"boolean","kind":"subtract","target":"base","tool":"upper"},
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
            "loops":[[
                {"kind":"line","id":"a","to":["2","0"]},
                {"kind":"line","id":"b","to":["2","2"]},
                {"kind":"line","id":"c","to":["0","2"]},
                {"kind":"line","id":"d","to":["0","0"]}]]},
        {"id":"ghost","op":"extrude_on_face","profile":"profile","body":"lower",
            "face":{"step":"base","role":"end"},"normal_direction":["0","0","1"],
            "x_direction":["1","0","0"],"offset":"0","distance":"1"}
    ],"result":"ghost"}).to_string()).unwrap()
}

#[test]
fn both_kernels_refuse_a_removed_end_face() {
    let recipe = recipe("0", "10");
    let catallax_features::BuildError::Geometry { refusal: error, .. } =
        catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap_err()
    else {
        panic!("expected a geometric refusal")
    };
    assert_eq!(error.code, "reference.lost-face");
    let error = catallax_field::Plan::compile(&recipe).unwrap_err();
    assert_eq!(error.code, "reference.lost-face");
}

#[test]
fn partial_and_split_faces_still_supply_their_supporting_plane() {
    for (cut_x, width) in [("0", "5"), ("4", "2")] {
        let recipe = recipe(cut_x, width);
        let mut dec = Decider::default();
        let mut built = catallax_features::evaluate(&recipe, &mut dec).unwrap();
        let point = catallax_shapes::RPoint3::new(1.into(), 1.into(), "10.5".parse().unwrap());
        assert!(matches!(
            catallax_boolean::classify_point(&mut built.model, built.result, &point, &mut dec),
            catallax_numbers::Outcome::Proven(catallax_boolean::Containment::Inside)
        ));
        let plan = catallax_field::Plan::compile_checked(&recipe, &mut dec).unwrap();
        let query = plan
            .proven_checked(&point, &mut Default::default(), &mut dec)
            .unwrap();
        assert_eq!(query.location, catallax_field::Location::Inside);
    }
}
