//! Merged fitted prisms must display within the unchanged application allowance.
use catallax_app::build_certified;
use catallax_mesh::{Mesh, MeshOptions};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_topology::{BodyId, Model, SplinePrism};
use serde_json::json;
use std::collections::BTreeMap;

fn recipe() -> Recipe {
    Recipe::parse(&json!({
        "format":"catallax-recipe/3", "name":"Equivalent offset prism display",
        "steps":[
            {"id":"profile","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"2"}]]},
            {"id":"sweep","op":"extrude","profile":"profile","distance":"10","direction":["3","0","4"]},
            {"id":"window","op":"block","size":["8","8","4"],"pose":{"rotation":["3","0","1","0"],"translation":["-2","-4","4"]}},
            {"id":"source","op":"boolean","kind":"intersect","target":"sweep","tool":"window"},
            {"id":"offset","op":"offset","body":"source","distance":"-1/4"},
            {"id":"half","op":"offset","body":"source","distance":"-1/8"},
            {"id":"other","op":"offset","body":"half","distance":"-1/8"},
            {"id":"combined","op":"boolean","kind":"union","target":"offset","tool":"other"}
        ], "result":"combined"
    }).to_string()).unwrap()
}

fn check_mesh(mesh: &Mesh, options: &MeshOptions) {
    let mut incidence = BTreeMap::<_, (usize, i32)>::new();
    let mut volume = 0.0;
    for triangle in &mesh.triangles {
        let [a, b, c] = triangle.vertices.map(|v| mesh.positions[v]);
        let ab: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
        let ac: [f64; 3] = std::array::from_fn(|i| c[i] - a[i]);
        let normal = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        assert!(normal.iter().any(|x| x.abs() > 1e-20));
        volume += a[0] * (b[1] * c[2] - b[2] * c[1])
            + a[1] * (b[2] * c[0] - b[0] * c[2])
            + a[2] * (b[0] * c[1] - b[1] * c[0]);
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let entry = incidence.entry((a.min(b), a.max(b))).or_default();
            entry.0 += 1;
            entry.1 += if a < b { 1 } else { -1 };
        }
    }
    assert!(volume > 0.0);
    assert!(
        incidence
            .values()
            .all(|&(count, winding)| count == 2 && winding == 0)
    );
    for face in &mesh.faces {
        assert!(face.max_chord_error <= options.chord_tolerance);
        if let Some(angle) = &options.angular_deviation {
            assert!(
                face.max_angular_deviation
                    .as_ref()
                    .is_some_and(|bound| bound <= angle)
            );
        }
    }
    for edge in &mesh.edges {
        assert_eq!(edge.vertices.first(), edge.vertices.last());
        assert!(edge.max_chord_error <= options.chord_tolerance);
    }
    assert!(mesh.unavailable_controls.is_empty());
}

fn check_sampled_angles(model: &Model, body: BodyId, mesh: &Mesh, limit: &Rational) {
    let mut dec = Decider::default();
    let prism = SplinePrism::check(model, body, &mut dec).unwrap();
    let trace = catallax_intersect::trace(
        model.geom(),
        model.edge(prism.rings()[0]).curve,
        None,
        None,
        &mut dec,
    )
    .unwrap()
    .value;
    let side = prism.faces()[2];
    // The fixture's section frame is perpendicular to its (3, 0, 4) axis.
    let polar = |p: [f64; 3]| p[1].atan2((4.0 * p[0] - 3.0 * p[2]) / 5.0);
    let wrap = |angle: f64| {
        (angle + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
    };
    for piece in trace.pieces() {
        for parameter in ["0", "1/2", "1"] {
            let t = catallax_numbers::Ball::from_rational(&parameter.parse().unwrap(), 128);
            let p = piece
                .evaluate(&t, &dec)
                .unwrap()
                .map(|x| x.mid().to_f64_nearest());
            let d = piece
                .derivative(&t, &dec)
                .unwrap()
                .map(|x| x.mid().to_f64_nearest());
            let normal = [4.0 * d[1], 3.0 * d[2] - 4.0 * d[0], -3.0 * d[1]];
            let mut checked = false;
            for triangle in mesh
                .triangles
                .iter()
                .filter(|triangle| triangle.face == side)
            {
                let points = triangle.vertices.map(|index| mesh.positions[index]);
                let angles = points.map(|point| wrap(polar(point) - polar(p)));
                let lo = angles.into_iter().fold(f64::INFINITY, f64::min);
                let hi = angles.into_iter().fold(f64::NEG_INFINITY, f64::max);
                if hi - lo > std::f64::consts::PI || lo > 1e-10 || hi < -1e-10 {
                    continue;
                }
                let a: [f64; 3] = std::array::from_fn(|i| points[1][i] - points[0][i]);
                let b: [f64; 3] = std::array::from_fn(|i| points[2][i] - points[0][i]);
                let facet = [
                    a[1] * b[2] - a[2] * b[1],
                    a[2] * b[0] - a[0] * b[2],
                    a[0] * b[1] - a[1] * b[0],
                ];
                let dot = (0..3).map(|i| facet[i] * normal[i]).sum::<f64>();
                let norm = |v: [f64; 3]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
                let angle = (dot / (norm(facet) * norm(normal))).clamp(-1.0, 1.0).acos();
                assert!(
                    angle <= limit.to_f64_nearest() + 1e-8,
                    "sampled angle {angle}"
                );
                checked = true;
            }
            assert!(checked, "spline sample has no covering facet");
        }
    }
}

fn display(options: &MeshOptions) {
    let (answer, work) = with_operation_budget(20_000_000, || {
        build_certified(&recipe(), Some(options), &mut Decider::default())
    });
    eprintln!("merged spline display work: {work:?}");
    let app = answer.unwrap();
    let mesh = &app.display().unwrap().as_ref().unwrap().value;
    assert!(!work.exhausted);
    assert_eq!(mesh.faces.len(), 3);
    check_mesh(mesh, options);
    if let Some(limit) = &options.angular_deviation {
        let built = app.built().built();
        check_sampled_angles(&built.model, built.result, mesh, limit);
    }
}

#[test]
fn merged_spline_default_display_fits_application_work() {
    display(&MeshOptions::default());
}

#[test]
fn merged_spline_soliddraft_angle_fits_application_work() {
    display(&MeshOptions {
        angular_deviation: Some(Rational::parse("0.03490658503988659").unwrap()),
        ..MeshOptions::default()
    });
}
