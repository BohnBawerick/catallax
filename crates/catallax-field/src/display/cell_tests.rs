use super::*;
use catallax_recipe::{Recipe, Role};
use serde_json::json;

#[test]
fn disjoint_cell_does_not_split_a_planar_face_by_its_other_constraints() {
    let pose = |p: [&str; 3]| json!({"rotation":["1","0","0","0"],"translation":p});
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"disjoint cells","steps":[
        {"id":"a","op":"block","size":["2","2","1"],"pose":pose(["-1","-1","-1"])},
        {"id":"b","op":"block","size":["1","1","1"],"pose":pose(["-1","-1","1"])},
        {"id":"result","op":"boolean","kind":"union","target":"a","tool":"b"}
    ],"result":"result"})
        .to_string(),
    )
    .unwrap();
    let mut dec = Decider::default();
    let mut plan = Plan::compile_checked(&recipe, &mut dec).unwrap();
    let leaf = plan
        .faces
        .iter()
        .find(|f| f.origin.step.as_str() == "b")
        .unwrap()
        .leaf;
    let radius = plan.leaves[leaf].display_radius(&mut dec).unwrap();
    let surfaces = plan
        .faces
        .iter()
        .filter(|f| f.leaf == leaf)
        .map(|f| f.surface.clone())
        .collect();
    // The same finite box represented by an analytic inequality cell. Its
    // z>=1 constraint excludes the entire z=0 cap of a, even though its x<=0
    // and y<=0 constraints change sign across that cap.
    plan.leaves[leaf] = crate::geometry::Leaf::Cell(crate::geometry::Cell::new(surfaces, radius));
    let mesh = mesh_plan(&plan, &Options::default(), &mut dec).unwrap();
    assert!(mesh.max_deviation <= Options::default().tolerance);
    let count = mesh
        .triangles
        .iter()
        .filter(|triangle| {
            let origin = &plan.faces[triangle.face.0].origin;
            origin.step.as_str() == "a" && origin.role == Role::End
        })
        .count();
    assert_eq!(
        count, 2,
        "a disjoint operand must not introduce trim seams into this plane"
    );
}
