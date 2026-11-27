//! The GPU input retains shared geometry and names without doing proof work.

use std::collections::BTreeMap;

use catallax_app::preview::display_description;
use catallax_features::evaluate;
use catallax_mesh::float::{DISPLAY_DESCRIPTION_VERSION, DisplayOptions, DisplaySurface};
use catallax_numbers::{Decider, with_operation_budget};
use catallax_recipe::Recipe;

#[test]
fn analytic_trims_share_samples_and_project_back_to_the_edge() {
    let recipe = Recipe::parse(include_str!("../benches/fine_mesh/countersink.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let (result, work) = with_operation_budget(0, || {
        display_description(&built, &DisplayOptions::default())
    });
    let display = result.unwrap();
    assert!(!work.exhausted);
    assert_eq!((work.used, work.uncharged), (0, 0));
    let desc = display.description;
    assert_eq!(desc.version, DISPLAY_DESCRIPTION_VERSION);
    let mut uses = BTreeMap::<usize, Vec<Vec<usize>>>::new();
    for face in &desc.faces {
        assert_eq!(
            display.names.faces[&face.face],
            built.model.face_origins(face.face)
        );
        for edge in face.loops.iter().flat_map(|l| &l.uses) {
            uses.entry(edge.edge)
                .or_default()
                .push(edge.samples.iter().map(|s| s.point).collect());
            for sample in &edge.samples {
                let [u, v] = sample.uv;
                let (f, x, y, z) = match face.surface {
                    DisplaySurface::Plane(f) => (f, u, v, 0.0),
                    DisplaySurface::Cylinder { frame, radius } => {
                        let (s, c) = (u * std::f64::consts::TAU).sin_cos();
                        (frame, radius * c, radius * s, v)
                    }
                    DisplaySurface::Cone { frame, slope } => {
                        let (s, c) = (u * std::f64::consts::TAU).sin_cos();
                        (frame, slope * v * c, slope * v * s, v)
                    }
                    ref other => panic!("unexpected countersink carrier {other:?}"),
                };
                let p: [f64; 3] =
                    std::array::from_fn(|i| f.origin[i] + x * f.x[i] + y * f.y[i] + z * f.z[i]);
                let q = desc.edges[edge.edge].points[sample.point];
                assert!(
                    p.into_iter().zip(q).all(|(a, b)| (a - b).abs() < 1e-8),
                    "{p:?} != {q:?}"
                );
            }
        }
    }
    assert_eq!(uses.len(), desc.edges.len());
    for (i, indices) in uses {
        let [a, b] = indices.as_slice() else {
            panic!("edge {i} lost an incident face")
        };
        assert!(a.iter().eq(b.iter().rev()));
        assert_eq!(a.len(), desc.edges[i].points.len());
    }
}

#[test]
fn spline_preview_keeps_the_retained_control_nets() {
    let recipe = Recipe::parse(include_str!(
        "../../catallax-features/tests/fixtures/spline-sweep-planar.json"
    ))
    .unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let (result, work) = with_operation_budget(0, || {
        display_description(&built, &DisplayOptions::default())
    });
    let result = result.unwrap().description;
    assert_eq!((work.used, work.uncharged), (0, 0));
    let mut count = 0;
    for face in &result.faces {
        let DisplaySurface::BSpline(net) = &face.surface else {
            continue;
        };
        let catallax_shapes::Surface::BSpline(source) = built
            .model
            .geom()
            .surface(built.model.face(face.face).surface)
        else {
            panic!("unexpected spline face {}", face.face.index());
        };
        let source = source.carrier();
        let [nu, nv] = source.control_count();
        assert_eq!(net.size, [nu, nv]);
        for v in 0..nv {
            for u in 0..nu {
                let p = &source.controls()[u * nv + v];
                let q = [&p.x, &p.y, &p.z].map(|r| r.to_f64());
                let q = std::array::from_fn(|i| q[i] - result.origin[i]);
                assert_eq!(net.points[v * nu + u], q);
            }
        }
        count += 1;
    }
    assert!(count > 0);
}
