//! Boolean history must not restore a removed face when its carrier returns.
use catallax_field::{FaceState, Plan, Workspace};
use catallax_numbers::Decider;
use catallax_recipe::{Origin, Recipe};
use catallax_shapes::RPoint3;
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn block(id: &str, x: &str, width: &str) -> Value {
    json!({"id":id,"op":"block","size":[width,"4","4"],
        "pose":{"rotation":["1","0","0","0"],"translation":[x,"0","0"]}})
}
fn boolean(id: &str, kind: &str, target: &str, tool: &str) -> Value {
    json!({"id":id,"op":"boolean","kind":kind,"target":target,"tool":tool})
}
fn cap(step: &str) -> Origin {
    serde_json::from_value(json!({"step":step,"role":{"side":{"segment":"x_max"}}})).unwrap()
}
fn verify(steps: Vec<Value>, x: i64, present: &[Origin], absent: &[Origin]) {
    let result = steps.last().unwrap()["id"].clone();
    let recipe = Recipe::parse(
        &json!({"format":catallax_recipe::FORMAT,
        "name":"boundary history","steps":steps,"result":result})
        .to_string(),
    )
    .unwrap();
    let field = Plan::compile(&recipe).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    let p = RPoint3::new(x.into(), 2.into(), 2.into());
    let mut workspace = Workspace::default();
    let mut actual = BTreeSet::new();
    for face in field.faces() {
        let id = field.select(&face.origin).unwrap();
        let q = field
            .closed_face_query_checked(id, &p, &mut workspace, &mut Decider::default())
            .unwrap();
        assert_ne!(q.state, FaceState::Indeterminate, "{:?}", face.origin);
        if q.state == FaceState::OnFace {
            actual.insert(face.origin.clone());
        }
    }
    assert_eq!(actual, present.iter().cloned().collect());
    for origin in present {
        assert!(
            !built
                .model
                .faces_with_origin(built.result, origin)
                .is_empty(),
            "{origin:?}"
        );
    }
    for origin in absent {
        assert!(
            built
                .model
                .faces_with_origin(built.result, origin)
                .is_empty(),
            "{origin:?}"
        );
    }
}

#[test]
fn union_then_intersection_does_not_restore_a_buried_creator() {
    verify(
        vec![
            block("old", "0", "2"),
            block("cover", "-1", "4"),
            block("fresh", "0", "2"),
            boolean("hidden", "union", "old", "cover"),
            boolean("result", "intersect", "hidden", "fresh"),
        ],
        2,
        &[cap("fresh")],
        &[cap("old")],
    );
}

#[test]
fn union_cannot_restore_a_clip_face_outside_the_intermediate_body() {
    verify(
        vec![
            block("seed", "0", "2"),
            block("clip", "0", "3"),
            block("fresh", "1", "2"),
            boolean("clipped", "intersect", "seed", "clip"),
            boolean("result", "union", "clipped", "fresh"),
        ],
        3,
        &[cap("fresh")],
        &[cap("clip")],
    );
}

#[test]
fn subtraction_cannot_reverse_a_previously_buried_name() {
    verify(
        vec![
            block("old", "0", "2"),
            block("cover", "0", "3"),
            block("cutter", "0", "2"),
            boolean("hidden", "union", "old", "cover"),
            boolean("result", "subtract", "hidden", "cutter"),
        ],
        2,
        &[cap("cutter")],
        &[cap("old")],
    );
}

#[test]
fn folded_empty_branch_does_not_borrow_a_shared_geometry_leaf() {
    verify(
        vec![
            block("old", "0", "2"),
            block("fresh", "0", "2"),
            boolean("empty", "subtract", "old", "old"),
            boolean("result", "union", "empty", "fresh"),
        ],
        2,
        &[cap("fresh")],
        &[cap("old")],
    );
}

#[test]
fn a_separate_surviving_path_keeps_its_coincident_origin() {
    verify(
        vec![
            block("old", "0", "2"),
            block("cover", "-1", "4"),
            block("fresh", "0", "2"),
            boolean("hidden", "union", "old", "cover"),
            boolean("restored", "intersect", "hidden", "fresh"),
            boolean("result", "union", "restored", "old"),
        ],
        2,
        &[cap("old"), cap("fresh")],
        &[],
    );
}

#[test]
fn a_coincident_copy_preserves_the_seed_snapshots_trim_history() {
    let copy = |origin| {
        catallax_features::copy_origin(&catallax_recipe::StepId::new("image").unwrap(), 1, &origin)
    };
    verify(
        vec![
            block("old", "0", "2"),
            block("cover", "-1", "4"),
            block("fresh", "0", "2"),
            boolean("hidden", "union", "old", "cover"),
            boolean("trimmed", "intersect", "hidden", "fresh"),
            json!({"id":"image","op":"mirror","body":"trimmed",
            "plane":{"point":["0","2","0"],"normal":["0","1","0"]},"keep_original":false}),
        ],
        2,
        &[copy(cap("fresh"))],
        &[copy(cap("old"))],
    );
}

#[test]
fn lattice_clipping_keeps_the_regions_intermediate_ownership() {
    let steps = vec![
        block("old", "0", "2"),
        block("cover", "-1", "4"),
        block("fresh", "0", "2"),
        boolean("hidden", "union", "old", "cover"),
        boolean("region", "intersect", "hidden", "fresh"),
        json!({"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"4",
            "frame":{"rotation":["1","0","0","0"],"translation":["0","2","2"]},
            "radius":{"constant":"1/2"}}),
    ];
    let recipe = Recipe::parse(
        &json!({"format":catallax_recipe::FORMAT,"name":"clipped ownership",
            "steps":steps,"result":"beams"})
        .to_string(),
    )
    .unwrap();
    let field = Plan::compile(&recipe).unwrap();
    let p = RPoint3::new(2.into(), 2.into(), 2.into());
    let mut workspace = Workspace::default();
    for (origin, expected) in [
        (cap("old"), FaceState::NotOnFace),
        (cap("fresh"), FaceState::OnFace),
    ] {
        let query = field
            .closed_face_query_checked(
                field.select(&origin).unwrap(),
                &p,
                &mut workspace,
                &mut Decider::default(),
            )
            .unwrap();
        assert_eq!(query.state, expected, "{origin:?}");
    }
}
