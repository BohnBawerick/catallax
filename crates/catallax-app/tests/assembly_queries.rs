//! Independent assembly queries on complete exact definitions.
use catallax_app::{
    assembly,
    queries::{BodyRelation, body_relation},
};
use catallax_features::evaluate;
use catallax_field::{
    Plan,
    assembly::{self as field, PlacedPart, Relation},
    mass,
};
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{Answer, Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{Model, block, sphere};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn exact(value: &str) -> Enclosure {
    Enclosure {
        lo: q(value),
        hi: q(value),
    }
}
fn pose(xyz: [i64; 3]) -> Pose {
    Pose::from_rational(
        Rotation::identity(),
        RVector3::new(xyz[0].into(), xyz[1].into(), xyz[2].into()),
    )
}
fn cube(model: &mut Model, xyz: [i64; 3], side: i64) -> catallax_topology::BodyId {
    block(model, [side.into(), side.into(), side.into()], &pose(xyz)).unwrap()
}
fn recipe() -> Recipe {
    Recipe::parse(r#"{"format":"catallax-recipe/3","name":"placed blocks","steps":[{"id":"a","op":"block","size":["1","2","3"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"b","op":"block","size":["1","2","3"],"pose":{"rotation":["1","0","0","0"],"translation":["4","6","0"]}}],"result":"a"}"#).unwrap()
}

#[test]
fn brep_relation_proves_global_gap_clash_containment_and_contact() {
    let mut model = Model::new();
    let a = cube(&mut model, [0, 0, 0], 2);
    let counts = model.counts(a);
    for (xyz, expected) in [([5, 6, 0], q("5")), ([0, 0, 5], q("3"))] {
        let b = cube(&mut model, xyz, 2);
        for (a, b) in [(a, b), (b, a)] {
            let answer =
                body_relation(&model, a, b, &Accuracy::default(), &mut Decider::default()).unwrap();
            let BodyRelation::Clear { gap } = answer.value else {
                panic!("{answer:?}")
            };
            assert!(gap.is_exact());
            assert_eq!(gap.lo, expected);
        }
    }
    for xyz in [[1, 0, 0], [0, 0, 0]] {
        let b = cube(&mut model, xyz, 2);
        assert!(
            matches!(body_relation(&model, a, b, &Accuracy::default(), &mut Decider::default()).unwrap().value, BodyRelation::Clash { overlap } if overlap.lo > q("0"))
        );
    }
    let outer = cube(&mut model, [-1, -1, -1], 4);
    assert!(matches!(
        body_relation(
            &model,
            a,
            outer,
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap()
        .value,
        BodyRelation::Clash { .. }
    ));
    for xyz in [[2, 0, 0], [2, 2, 0], [2, 2, 2]] {
        let b = cube(&mut model, xyz, 2);
        assert_eq!(
            body_relation(&model, a, b, &Accuracy::default(), &mut Decider::default())
                .unwrap_err()
                .refusal
                .code,
            "assembly.contact"
        );
    }
    assert_eq!(model.counts(a), counts);
}

#[test]
fn spheres_and_irrational_clearance_keep_certified_enclosures() {
    let mut model = Model::new();
    let a = sphere(&mut model, 1.into(), &pose([0, 0, 0])).unwrap();
    let b = sphere(&mut model, 1.into(), &pose([3, 0, 0])).unwrap();
    let BodyRelation::Clear { gap } =
        body_relation(&model, a, b, &Accuracy::default(), &mut Decider::default())
            .unwrap()
            .value
    else {
        panic!()
    };
    assert_eq!(gap, exact("1"));
    let a = cube(&mut model, [0, 0, 0], 1);
    let b = cube(&mut model, [2, 2, 0], 1);
    let BodyRelation::Clear { gap } =
        body_relation(&model, a, b, &Accuracy::default(), &mut Decider::default())
            .unwrap()
            .value
    else {
        panic!()
    };
    assert!(&gap.lo * &gap.lo <= q("2") && &gap.hi * &gap.hi >= q("2"));
}

#[test]
fn independent_recipe_queries_report_one_mass_and_one_gap() {
    let recipe = recipe();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let a = PlacedPart::new(&plan, &Pose::identity(), &mut Decider::default()).unwrap();
    let b = PlacedPart::new(&plan, &pose([4, 6, 0]), &mut Decider::default()).unwrap();
    let brep = body_relation(
        &built.model,
        built.result,
        built.bodies[&catallax_recipe::StepId::new("b").unwrap()],
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let field = field::relation(
        &a,
        &b,
        &field::Options {
            gap_width: q("1e-12"),
            max_pairs: 64,
        },
        &mut Decider::default(),
    );
    let agreed = assembly::relation_agreement(&brep, &field, &Accuracy::default()).unwrap();
    assert!(
        matches!(agreed.value, BodyRelation::Clear { gap } if gap.is_exact() && gap.lo == q("5"))
    );
    let density = q("0.00785");
    let accuracy = mass::Accuracy::default();
    let brep =
        assembly::uniform_mass(&built, &density, &accuracy, &mut Decider::default()).unwrap();
    let field = field::uniform_mass(
        &a,
        &density,
        &accuracy,
        &mass::Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let agreed = assembly::mass_agreement(&brep, &field, &density, &accuracy).unwrap();
    assert_eq!(agreed.value.grams, exact("0.0471"));
    assert_eq!(
        agreed.value.centre_of_gravity.unwrap().map(|v| v.lo),
        [q("0.5"), q("1"), q("1.5")]
    );
}

#[test]
fn agreement_uses_a_hull_and_refuses_mismatch_wide_bounds_and_missing_centres() {
    let make = |lo: &str, hi: &str| assembly::Mass {
        grams: Enclosure {
            lo: q(lo),
            hi: q(hi),
        },
        centre_of_gravity: Some(std::array::from_fn(|_| exact("0"))),
    };
    let accuracy = mass::Accuracy {
        volume: Accuracy {
            relative: q("0"),
            absolute: q("0.2"),
        },
        ..Default::default()
    };
    let brep = Answer::exact(make("1", "1.1"));
    let field = make("1.05", "1.15");
    let result = assembly::mass_agreement(&brep, &field, &q("1"), &accuracy).unwrap();
    assert_eq!(
        result.value.grams,
        Enclosure {
            lo: q("1"),
            hi: q("1.15")
        }
    );
    assert_eq!(
        assembly::mass_agreement(&brep, &make("1.2", "1.3"), &q("1"), &accuracy)
            .unwrap_err()
            .code,
        "assembly.disagreement"
    );
    assert_eq!(
        assembly::mass_agreement(&brep, &make("1.05", "2"), &q("1"), &accuracy)
            .unwrap_err()
            .code,
        "assembly.agreement-width"
    );
    let mut field = field;
    field.centre_of_gravity = None;
    assert_eq!(
        assembly::mass_agreement(&brep, &field, &q("1"), &accuracy)
            .unwrap_err()
            .code,
        "assembly.mass-centroid"
    );
    assert_eq!(
        assembly::mass_agreement(&brep, &make("1", "1"), &q("0"), &accuracy)
            .unwrap_err()
            .code,
        "assembly.density"
    );
    let brep = Answer::exact(BodyRelation::Clear { gap: exact("1") });
    assert_eq!(
        assembly::relation_agreement(
            &brep,
            &Relation::Clash {
                witness: [q("0"), q("0"), q("0")]
            },
            &Accuracy::default()
        )
        .unwrap_err()
        .code,
        "assembly.relation-disagreement"
    );
}

#[test]
fn envelope_mass_and_exhausted_relation_never_become_answers() {
    let recipe = Recipe::parse(include_str!("../../../parts/P6-lattice-bracket.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    assert_eq!(
        assembly::uniform_mass(
            &built,
            &q("1"),
            &mass::Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .refusal
        .code,
        "assembly.mass-envelope"
    );
    let mut model = Model::new();
    let a = cube(&mut model, [0, 0, 0], 2);
    let b = cube(&mut model, [1, 0, 0], 2);
    let (answer, work) = with_operation_budget(0, || {
        body_relation(&model, a, b, &Accuracy::default(), &mut Decider::default())
    });
    assert!(answer.is_err());
    assert!(work.exhausted);
}

#[test]
fn cavity_gap_tiny_positive_gap_empty_and_invalid_requests() {
    let mut model = Model::new();
    let outer = cube(&mut model, [0, 0, 0], 10);
    let inner = cube(&mut model, [1, 1, 1], 8);
    let shell = catallax_boolean::boolean(
        &mut model,
        outer,
        inner,
        catallax_boolean::BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let inside = cube(&mut model, [3, 3, 3], 1);
    let BodyRelation::Clear { gap } = body_relation(
        &model,
        shell,
        inside,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value
    else {
        panic!()
    };
    assert_eq!(gap, exact("2"));
    let a = cube(&mut model, [0, 0, 0], 1);
    let b = block(
        &mut model,
        [1.into(), 1.into(), 1.into()],
        &Pose::from_rational(
            Rotation::identity(),
            RVector3::new(q("1.000000000001"), 0.into(), 0.into()),
        ),
    )
    .unwrap();
    let BodyRelation::Clear { gap } =
        body_relation(&model, a, b, &Accuracy::default(), &mut Decider::default())
            .unwrap()
            .value
    else {
        panic!()
    };
    assert!(gap.contains(&q("1e-12")));
    assert!(gap.lo > q("0") && gap.width() <= q("1e-24"));
    let empty = catallax_boolean::boolean(
        &mut model,
        a,
        a,
        catallax_boolean::BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    assert_eq!(
        body_relation(
            &model,
            a,
            empty,
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .refusal
        .code,
        "app.relation-empty"
    );
    let bad = Accuracy {
        absolute: q("-1"),
        relative: q("0"),
    };
    assert_eq!(
        body_relation(&model, a, b, &bad, &mut Decider::default())
            .unwrap_err()
            .refusal
            .code,
        "app.relation-accuracy"
    );
}
