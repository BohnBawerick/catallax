//! Selected-edge material changes and fillet provenance from the independent field.

use catallax_field::{FaceState, Location, Plan, Workspace};
use catallax_numbers::{BudgetPolicy, Decider, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::RPoint3;

fn point(p: [&str; 3]) -> RPoint3 {
    let [x, y, z] = p.map(|s| Rational::parse(s).unwrap());
    RPoint3::new(x, y, z)
}

#[test]
fn p3_convex_removal_and_concave_addition_are_independent_and_named() {
    let recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    for (p, want) in [
        (["1/10", "15", "99/10"], Location::Outside),
        (["3/5", "15", "49/5"], Location::Inside),
        (["2/5", "15", "49/5"], Location::Boundary),
        (["147/5", "15", "51/5"], Location::Outside),
        (["142/5", "15", "51/5"], Location::Inside),
        (["144/5", "15", "52/5"], Location::Boundary),
    ] {
        let p = point(p);
        let query = plan.proven_checked(&p, &mut workspace, &mut dec).unwrap();
        assert_eq!(query.location, want, "at {p:?}: {query:?}");
    }
    let junction = plan
        .faces()
        .iter()
        .find(|f| {
            f.origin.step.as_str() == "junction" && matches!(f.origin.role, Role::Blend { .. })
        })
        .unwrap();
    let face = plan.select(&junction.origin).unwrap();
    let query = plan
        .face_query_checked(
            face,
            &point(["144/5", "15", "52/5"]),
            &mut workspace,
            &mut dec,
        )
        .unwrap();
    assert_eq!(query.state, FaceState::OnFace);
}

#[test]
fn convex_ring_keeps_the_unselected_end_sharp() {
    let recipe = Recipe::parse(r#"{
        "format":"catallax-recipe/3","name":"convex cylinder",
        "steps":[
            {"id":"base","op":"cylinder","radius":"8","height":"12",
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"base","radius":"2","edges":[
                {"between":[{"step":"base","role":"end"},{"step":"base","role":{"side":{"segment":"wall"}}}]}
            ]}
        ],"result":"round"
    }"#).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    for (p, expected) in [
        (["38/5", "0", "58/5"], Location::Outside),
        (["36/5", "0", "58/5"], Location::Boundary),
        (["7", "0", "11"], Location::Inside),
        (["8", "0", "0"], Location::Boundary),
        (["39/5", "0", "1/5"], Location::Inside),
    ] {
        assert_eq!(
            plan.proven_checked(&point(p), &mut workspace, &mut dec)
                .unwrap()
                .location,
            expected,
            "at {p:?}"
        );
    }
}

#[test]
fn reflected_blend_keeps_its_named_boundary_and_material_sector() {
    let recipe = Recipe::parse(r#"{
        "format":"catallax-recipe/3","name":"placed fillet",
        "steps":[
            {"id":"base","op":"cylinder","radius":"8","height":"12",
             "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
            {"id":"round","op":"fillet","body":"base","radius":"2","edges":[
                {"between":[{"step":"base","role":"end"},{"step":"base","role":{"side":{"segment":"wall"}}}]}
            ]},
            {"id":"copies","op":"pattern","body":"round",
             "layout":{"linear":{"vector":["30","0","0"],"count":2}}},
            {"id":"reflected","op":"mirror","body":"copies",
             "plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}
        ],"result":"reflected"
    }"#).unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    for (p, expected) in [
        (["-37", "0", "11"], Location::Inside),
        (["-188/5", "0", "58/5"], Location::Outside),
        (["-186/5", "0", "58/5"], Location::Boundary),
        (["7", "0", "11"], Location::Inside),
        (["37", "0", "11"], Location::Outside),
    ] {
        assert_eq!(
            plan.proven_checked(&point(p), &mut workspace, &mut dec)
                .unwrap()
                .location,
            expected,
            "at {p:?}"
        );
    }
    let mut named_boundary = 0;
    for face in plan.faces().iter().filter(|face| {
        face.origin.step.as_str() == "reflected"
            && matches!(face.surface, catallax_shapes::Surface::Torus(_))
    }) {
        let selected = plan.select(&face.origin).unwrap();
        if plan
            .face_query_checked(
                selected,
                &point(["-186/5", "0", "58/5"]),
                &mut workspace,
                &mut dec,
            )
            .unwrap()
            .state
            == FaceState::OnFace
        {
            named_boundary += 1;
        }
    }
    assert_eq!(named_boundary, 1);
}

#[test]
fn trihedral_corner_removes_the_cube_beyond_its_ball_and_names_the_sector() {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/trihedral-corner.json"
    ))
    .unwrap();
    let mut dec = Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    });
    let plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let mut workspace = Workspace::default();
    // The ball is centred at 8,8,8 with radius 2. The 93/10 diagonal point
    // lies in all three strip cylinders but outside the ball.
    for (p, expected) in [
        (["91/10", "91/10", "91/10"], Location::Inside),
        (["93/10", "93/10", "93/10"], Location::Outside),
        (["93/10", "93/10", "5"], Location::Inside),
        (["8", "8", "6"], Location::Inside),
        (["26/3", "28/3", "28/3"], Location::Boundary),
        (["10", "10", "10"], Location::Outside),
    ] {
        assert_eq!(
            plan.proven_checked(&point(p), &mut workspace, &mut dec)
                .unwrap()
                .location,
            expected,
            "at {p:?}"
        );
    }
    let corners: Vec<_> = plan
        .faces()
        .iter()
        .filter(|f| matches!(f.origin.role, Role::Corner { .. }))
        .collect();
    let [corner] = corners.as_slice() else {
        panic!("one corner origin expected: {corners:?}");
    };
    let Role::Corner { vertex } = &corner.origin.role else {
        unreachable!()
    };
    assert_eq!(vertex.faces().len(), 3);
    let face = plan.select(&corner.origin).unwrap();
    for (p, expected) in [
        (["26/3", "28/3", "28/3"], FaceState::OnFace),
        (["8", "8", "6"], FaceState::NotOnFace),
    ] {
        assert_eq!(
            plan.face_query_checked(face, &point(p), &mut workspace, &mut dec)
                .unwrap()
                .state,
            expected,
            "at {p:?}"
        );
    }
}
