//! Stock solids share the exact recipe number and naming contracts.

use catallax_recipe::{EdgeName, ElementId, Operation, Origin, Recipe, Role, StepId, VertexName};

const SOLIDS: &str = r#"{
  "format":"catallax-recipe/3", "name":"stock solids",
  "steps":[
    {"id":"box", "op":"block", "size":["2.00","3/2","4"], "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"tube", "op":"cylinder", "radius":"1/3", "height":"7", "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"ball", "op":"sphere", "radius":"5", "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"point", "op":"cone", "base_radius":"3", "top_radius":"0", "height":"4", "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"ring", "op":"torus", "major":"4", "minor":"1", "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
    {"id":"part", "op":"boolean", "kind":"subtract", "target":"box", "tool":"tube"}
  ], "result":"part"
}"#;

fn origin(step: &str, role: Role) -> Origin {
    Origin {
        step: StepId::new(step).unwrap(),
        role,
    }
}
fn side(step: &str, segment: &str) -> Origin {
    origin(
        step,
        Role::Side {
            segment: ElementId::new(segment).unwrap(),
            path: None,
        },
    )
}

#[test]
fn primitive_definitions_preserve_exact_parameters_and_placement() {
    let recipe = Recipe::parse(SOLIDS).unwrap();
    let canonical = recipe.to_text();
    let restored = Recipe::parse(&canonical).unwrap();
    assert_eq!(restored, recipe);
    let Operation::Block { size, .. } = &restored.steps[0].op else {
        panic!("block lost");
    };
    assert_eq!(size[1].value().to_string(), "3/2");
    let Operation::Cone { top_radius, .. } = &restored.steps[3].op else {
        panic!("cone lost");
    };
    assert!(top_radius.value().is_zero());
    for threads in [1, 2, 4] {
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..threads)
                .map(|_| scope.spawn(|| Recipe::parse(SOLIDS).unwrap().to_text()))
                .collect();
            for handle in handles {
                assert_eq!(handle.join().unwrap(), canonical);
            }
        });
    }
}

#[test]
fn primitive_faces_edges_and_singular_apex_have_recipe_names() {
    let x = side("box", "x_min");
    let y = side("box", "y_min");
    let start = origin("box", Role::Start);
    let edge = EdgeName::new(x.clone(), y.clone()).unwrap();
    assert_eq!(edge, EdgeName::new(y.clone(), x.clone()).unwrap());
    let corner = VertexName::new(vec![start.clone(), y.clone(), x.clone()]).unwrap();
    let apex = VertexName::new(vec![side("point", "wall")]).unwrap();
    assert_eq!(apex.faces(), [side("point", "wall")]);
    assert!(VertexName::new(vec![]).is_err());
    let mut recipe = Recipe::parse(SOLIDS).unwrap();
    recipe.steps.push(catallax_recipe::Step {
        id: StepId::new("round").unwrap(),
        op: Operation::Fillet {
            body: recipe.result.clone(),
            radius: 1.into(),
            edges: catallax_recipe::EdgeSelection::Refs(vec![catallax_recipe::EdgeRef {
                between: [
                    catallax_recipe::OriginPattern {
                        step: x.step.clone(),
                        role: Some(x.role.clone()),
                    },
                    catallax_recipe::OriginPattern {
                        step: y.step.clone(),
                        role: Some(y.role.clone()),
                    },
                ],
                touching: vec![],
                near: None,
            }]),
        },
    });
    assert!(recipe.check().is_ok());
    assert_eq!(Recipe::parse(&recipe.to_text()).unwrap(), recipe);
    assert_eq!(corner, VertexName::new(vec![x, start, y]).unwrap());
}

#[test]
fn nonexistent_primitive_roles_are_rejected_before_geometry() {
    let mut recipe = Recipe::parse(SOLIDS).unwrap();
    let reference = |step: &str, role| catallax_recipe::FaceRef {
        face: catallax_recipe::OriginPattern {
            step: StepId::new(step).unwrap(),
            role: Some(role),
        },
        touching: vec![],
        near: None,
    };
    recipe.steps.push(catallax_recipe::Step {
        id: StepId::new("shell").unwrap(),
        op: Operation::Shell {
            body: recipe.result.clone(),
            open: vec![reference("point", Role::End)],
            thickness: 1.into(),
        },
    });
    assert_eq!(recipe.check().unwrap_err().code(), "recipe.wrong-role");
    let Operation::Shell { open, .. } = &mut recipe.steps.last_mut().unwrap().op else {
        unreachable!()
    };
    open[0] = reference(
        "box",
        Role::Side {
            segment: ElementId::new("not_a_face").unwrap(),
            path: None,
        },
    );
    assert_eq!(recipe.check().unwrap_err().code(), "recipe.wrong-role");
}
