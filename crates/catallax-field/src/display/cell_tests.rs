use super::*;
use catallax_recipe::{Recipe, Role};
use serde_json::json;

fn default_mesh(plan: &Plan) -> Mesh {
    let options = Options::default();
    let (result, work) = catallax_numbers::with_operation_budget(100_000_000, || {
        mesh_plan(plan, &options, &mut Decider::default())
    });
    assert!(!work.exhausted);
    let mesh = result.unwrap();
    assert!(mesh.max_deviation <= options.tolerance);
    assert!(mesh.cells <= options.max_cells && mesh.positions.len() <= options.max_vertices);
    mesh
}

fn as_cell(plan: &mut Plan, step: &str) {
    let leaf = plan
        .faces
        .iter()
        .find(|f| f.origin.step.as_str() == step)
        .unwrap()
        .leaf;
    let radius = plan.leaves[leaf]
        .display_radius(&mut Decider::default())
        .unwrap();
    let surfaces = plan
        .faces
        .iter()
        .filter(|f| f.leaf == leaf)
        .map(|f| f.surface.clone())
        .collect();
    plan.leaves[leaf] = crate::geometry::Leaf::Cell(crate::geometry::Cell::new(surfaces, radius));
}

fn covers(
    mesh: &Mesh,
    plan: &Plan,
    step: &str,
    role: &Role,
    point: [Rational; 3],
    axis: usize,
) -> bool {
    let axes: Vec<_> = (0..3).filter(|&i| i != axis).collect();
    let cross = |a: &[Rational; 2], b: &[Rational; 2], c: &[Rational; 2]| {
        &(&b[0] - &a[0]) * &(&c[1] - &a[1]) - &(&b[1] - &a[1]) * &(&c[0] - &a[0])
    };
    let p = [point[axes[0]].clone(), point[axes[1]].clone()];
    mesh.triangles.iter().any(|triangle| {
        let origin = &plan.faces[triangle.face.0].origin;
        if origin.step.as_str() != step || &origin.role != role {
            return false;
        }
        let xyz = triangle
            .vertices
            .map(|i| mesh.positions[i].map(|x| Rational::from_f64_exact(x).unwrap()));
        if xyz.iter().any(|q| q[axis] != point[axis]) {
            return false;
        }
        let q = xyz.map(|q| [q[axes[0]].clone(), q[axes[1]].clone()]);
        if cross(&q[0], &q[1], &q[2]).is_zero() {
            return false;
        }
        let sides = [
            cross(&q[0], &q[1], &p),
            cross(&q[1], &q[2], &p),
            cross(&q[2], &q[0], &p),
        ];
        sides.iter().all(|x| x >= &Rational::zero()) || sides.iter().all(|x| x <= &Rational::zero())
    })
}

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
    let mesh = default_mesh(&plan);
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

#[test]
fn crossing_cell_constraints_preserve_the_cut_and_its_closed_rails() {
    let pose = |p: [&str; 3]| json!({"rotation":["1","0","0","0"],"translation":p});
    let recipe = Recipe::parse(
        &json!({"format":"catallax-recipe/3","name":"cell notch","steps":[
        {"id":"a","op":"block","size":["2","2","1"],"pose":pose(["-1","-1","-1"])},
        {"id":"b","op":"block","size":["1","1","2"],"pose":pose(["0","0","-1"])},
        {"id":"result","op":"boolean","kind":"subtract","target":"a","tool":"b"}
    ],"result":"result"})
        .to_string(),
    )
    .unwrap();
    let mut plan = Plan::compile(&recipe).unwrap();
    as_cell(&mut plan, "b");
    let mesh = default_mesh(&plan);
    let h = Rational::from_ratio(1, 2).unwrap();
    assert!(covers(
        &mesh,
        &plan,
        "a",
        &Role::End,
        [-&h, h.clone(), q(0)],
        2
    ));
    assert!(covers(
        &mesh,
        &plan,
        "a",
        &Role::End,
        [h.clone(), -&h, q(0)],
        2
    ));
    assert!(!covers(
        &mesh,
        &plan,
        "a",
        &Role::End,
        [h.clone(), h.clone(), q(0)],
        2
    ));
    let x = Role::Side {
        segment: catallax_recipe::ElementId::new("x_min").unwrap(),
        path: None,
    };
    let y = Role::Side {
        segment: catallax_recipe::ElementId::new("y_min").unwrap(),
        path: None,
    };
    assert!(covers(&mesh, &plan, "b", &x, [q(0), h.clone(), -&h], 0));
    assert!(covers(&mesh, &plan, "b", &y, [h.clone(), q(0), -&h], 1));
    assert!(covers(&mesh, &plan, "b", &x, [q(0), q(0), -&h], 0));
    assert!(covers(&mesh, &plan, "b", &y, [q(0), q(0), -&h], 1));
}

#[test]
fn concave_feature_keeps_winding_membership_despite_positive_supports() {
    let recipe=Recipe::parse(&json!({"format":"catallax-recipe/3","name":"concave winding","steps":[
        {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
            {"kind":"line","id":"a","to":["3","0"]},{"kind":"line","id":"b","to":["3","1"]},{"kind":"line","id":"c","to":["1","1"]},
            {"kind":"line","id":"d","to":["1","3"]},{"kind":"line","id":"e","to":["0","3"]},{"kind":"line","id":"f","to":["0","0"]}
        ]]}, {"id":"solid","op":"extrude","profile":"profile","distance":"2"}
    ],"result":"solid"}).to_string()).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    let mesh = default_mesh(&plan);
    let h = Rational::from_ratio(1, 2).unwrap();
    let far = Rational::from_ratio(5, 2).unwrap();
    assert!(covers(
        &mesh,
        &plan,
        "solid",
        &Role::End,
        [far.clone(), h.clone(), q(2)],
        2
    ));
    assert!(covers(&mesh, &plan, "solid", &Role::End, [h, far, q(2)], 2));
    assert!(!covers(
        &mesh,
        &plan,
        "solid",
        &Role::End,
        [q(2), q(2), q(2)],
        2
    ));
}

#[test]
fn analytic_cylinder_cell_displays_all_faces_with_default_controls() {
    let recipe=Recipe::parse(&json!({"format":"catallax-recipe/3","name":"cell cylinder","steps":[
        {"id":"solid","op":"cylinder","radius":"1","height":"2","pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}
    ],"result":"solid"}).to_string()).unwrap();
    let mut plan = Plan::compile(&recipe).unwrap();
    as_cell(&mut plan, "solid");
    let mesh = default_mesh(&plan);
    let faces = mesh
        .triangles
        .iter()
        .map(|t| t.face.0)
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(faces.len(), 3);
    assert!(!mesh.positions.is_empty());
}
