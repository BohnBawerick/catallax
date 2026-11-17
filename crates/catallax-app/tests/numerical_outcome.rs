//! Numerical application retention survives source release and downstream export.

use catallax_app::{NumericalOptions, build_numerical};
use catallax_numbers::{Ball, Decider, NumericalRepair};
use catallax_shapes::Surface;

mod numerical_cases;
use numerical_cases::{from_model, q};

#[test]
fn application_outcome_retains_repair_evidence_and_export_provenance() {
    use catallax_app::{EvaluatedApplicationOutcome, NumericalApplicationBody};
    use std::sync::Arc;

    let outcome = {
        let recipe = Arc::new(
            catallax_recipe::Recipe::parse(
                r#"{
            "format":"catallax-recipe/4","name":"retained numerical cube","result":"stock",
            "steps":[
                {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],
                 "translation":["0","0","0"]},"loops":[[
                    {"id":"a","kind":"line","to":["1","0"]},
                    {"id":"b","kind":"line","to":["1","1"]},
                    {"id":"c","kind":"line","to":["0","1"]},
                    {"id":"d","kind":"line","to":["0","0"]}
                ]]},
                {"id":"stock","op":"extrude","profile":"profile",
                 "direction":["0","0","1"],"distance":"1"}
            ]
        }"#,
            )
            .unwrap(),
        );
        let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
        let mut definition = from_model(&built.model, built.result, "0");
        let origins = built
            .model
            .faces_of_body(built.result)
            .into_iter()
            .zip(&definition.faces)
            .map(|(face, numerical)| {
                (
                    numerical.key.clone(),
                    built.model.face_origins(face).to_vec(),
                )
            })
            .collect();
        let position = definition.vertices[0].position.clone();
        let vertex = definition.vertices[0].key.clone();
        definition.vertices[0].position[0] = &position[0] + &q("0.000001");
        let admitted = build_numerical(
            definition,
            NumericalOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let repaired = admitted
            .repair_vertex(&vertex, position, &mut Decider::default())
            .unwrap();
        EvaluatedApplicationOutcome::Numerical(Arc::new(NumericalApplicationBody::new(
            repaired, recipe, origins,
        )))
    };
    let EvaluatedApplicationOutcome::Numerical(retained) = outcome else {
        panic!("numerical geometry must retain its authority");
    };
    let artifact = catallax_step::numerical::build(
        retained.body(),
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let [
        NumericalRepair::VertexMoved {
            distance_upper_mm, ..
        },
    ] = artifact.report().source_evidence.repairs()
    else {
        panic!("retention must not rebuild away the actual repair");
    };
    assert!(*distance_upper_mm >= q("0.000001"));
    assert!(distance_upper_mm - &q("0.000001") < q("1e-25"));
    for face in &artifact.report().source_definition.faces {
        let origins = &retained.face_origins()[&face.key];
        assert!(
            !origins.is_empty(),
            "exported face {} lost its origins",
            face.key
        );
        for origin in origins {
            assert!(
                retained
                    .recipe()
                    .steps
                    .iter()
                    .any(|step| step.id == origin.step)
            );
        }
    }
    let Surface::BSpline(roof) = &retained
        .body()
        .definition()
        .faces
        .iter()
        .find(|face| face.key == "roof")
        .unwrap()
        .surface
    else {
        panic!("retention must preserve the admitted carrier");
    };
    let at = Ball::from_rational(&q("1/2"), 128);
    let point = roof
        .carrier()
        .evaluate(&at, &at, &Decider::default())
        .unwrap();
    assert!(point[2].lo() <= q("1") && point[2].hi() >= q("1"));
    assert!(point[2].hi() - point[2].lo() < q("1e-25"));
}
