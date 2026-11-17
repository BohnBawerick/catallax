use catallax_boolean::{BoolOp, classify, imprint, intersect, keep, sew};
use catallax_measure::{Accuracy, face_area};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::Surface;
use catallax_topology::{FaceId, Model};

fn show(model: &Model, label: &str, faces: &[FaceId], dec: &mut Decider) {
    for &face in faces {
        let entity = model.face(face);
        let Surface::Plane(plane) = model.geom().surface(entity.surface) else { continue };
        if plane.normal.z.to_f64().abs() < 0.9 { continue; }
        let area = face_area(model, face, &Accuracy {
            relative: Rational::zero(), absolute: Rational::parse("0.0001").unwrap(),
        }, dec);
        println!("{label}: {face:?} surface={:?} z={} normal_z={} sense={:?} area={area:?} origins={:?}", entity.surface, plane.origin.z.to_f64(), plane.normal.z.to_f64(), entity.sense, model.face_origins(face));
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (answer, work) = with_operation_budget(100_000_000, || -> Result<(), Box<dyn std::error::Error>> {
        let mut input: serde_json::Value = serde_json::from_str(include_str!("../tests/fixtures/coplanar-start-cap.json"))?;
        input["steps"].as_array_mut().unwrap().pop();
        input["result"] = "base".into();
        let recipe = Recipe::parse(&input.to_string())?;
        let mut built = catallax_features::evaluate(&recipe, &mut Decider::default())?;
        let a = built.bodies[&StepId::new("base")?];
        let b = built.bodies[&StepId::new("pocket")?];
        let model = &mut built.model;
        let mut dec = Decider::default();
        show(model, "input a", &model.faces_of_body(a), &mut dec);
        show(model, "input b", &model.faces_of_body(b), &mut dec);
        let graph = intersect(model, a, b, &mut dec)?.value;
        println!("sections: {}", graph.edges().len());
        for edge in graph.edges() {
            println!("section {:?} {:?} {:?} a={:?} b={:?}", edge.curve, edge.start, edge.end, edge.faces_of_first, edge.faces_of_second);
        }
        let split = imprint(model, &graph, &mut dec)?.value;
        let first = split.pieces_of_first();
        let second = split.pieces_of_second();
        println!("imprint a={} b={}", first.len(), second.len());
        show(model, "split a", &first, &mut dec);
        show(model, "split b", &second, &mut dec);
        let classes = classify(model, &split, &mut dec)?.value;
        for &face in first.iter().chain(&second) {
            println!("class {face:?} {:?} {:?}", classes.class_of(face), model.face_origins(face));
        }
        let selected = keep(&classes, &split, BoolOp::Subtract);
        let result = sew(model, &split, &selected, &mut dec)?.value;
        let faces = model.faces_of_body(result);
        println!("result faces={} faults={:?}", faces.len(), model.check(result));
        show(model, "result", &faces, &mut dec);
        Ok(())
    });
    println!("work={work:?}");
    answer
}
