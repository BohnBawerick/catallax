//! Construct six separate native faces, sew them, and export checked material.
use catallax_import::{
    Curve, Face, FaceLoop, OrientedEdge, Placement, Surface, Vertex,
    numerical::sew_constructed_surfaces,
    tolerant::{ToleranceOptions, healing::HealingChange},
};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_topology::{BodyKind, numerical::NumericalOptions};

fn faces() -> Vec<Face> {
    let corners: [[Rational; 3]; 8] = std::array::from_fn(|i| {
        std::array::from_fn(|axis| {
            Rational::from_i64(if i & (4 >> axis) == 0 {
                0
            } else {
                [2, 3, 4][axis]
            })
        })
    });
    [
        [0, 1, 3, 2],
        [4, 6, 7, 5],
        [0, 4, 5, 1],
        [2, 3, 7, 6],
        [0, 2, 6, 4],
        [1, 5, 7, 3],
    ]
    .into_iter()
    .enumerate()
    .map(|(i, cycle)| {
        let axis = i / 2;
        let normal = std::array::from_fn(|j| {
            Rational::from_i64(if j == axis {
                if i % 2 == 0 { -1 } else { 1 }
            } else {
                0
            })
        });
        let ref_dir = std::array::from_fn(|j| Rational::from_i64(i64::from(j == (axis + 1) % 3)));
        let vertex = |j: usize| Vertex {
            name: format!("side-{i}/vertex-{j}"),
            point: corners[cycle[j]].clone(),
        };
        let edges = (0..4)
            .map(|j| {
                let start = vertex(j);
                let end = vertex((j + 1) % 4);
                let direction = std::array::from_fn(|axis| &end.point[axis] - &start.point[axis]);
                OrientedEdge {
                    edge: format!("side-{i}/edge-{j}"),
                    curve: Curve::Line(start.point.clone(), direction),
                    start,
                    end,
                    curve_same_sense: true,
                    forward: true,
                }
            })
            .collect();
        Face {
            name: format!("side-{i}"),
            surface: Surface::Plane(Placement {
                origin: corners[cycle[0]].clone(),
                axis: normal,
                ref_dir,
            }),
            same_sense: true,
            loops: vec![FaceLoop {
                outer: true,
                orientation: true,
                edges,
                vertex: None,
            }],
        }
    })
    .collect()
}

fn main() {
    let folder = std::path::PathBuf::from(std::env::args().nth(1).expect("output directory"));
    std::fs::create_dir_all(&folder).unwrap();
    let (result, work) = with_operation_budget(catallax_topology::numerical::DEFAULT_WORK, || {
        let mut dec = Decider::default();
        let faces = faces();
        let tolerances = ToleranceOptions {
            initial_mm: Rational::zero(),
            maximum_mm: Rational::parse("0.000001").unwrap(),
        };
        let options = NumericalOptions {
            tolerance_mm: Rational::parse("0.000001").unwrap(),
            ..Default::default()
        };
        assert!(
            sew_constructed_surfaces(
                "open-box",
                &faces[..5],
                tolerances.clone(),
                options.clone(),
                &mut dec
            )
            .is_err()
        );
        let sewn =
            sew_constructed_surfaces("constructed-box", &faces, tolerances, options, &mut dec)
                .unwrap();
        let source = &sewn.healing().source().bodies()[0];
        assert!(!source.declared_solid());
        assert_eq!(source.vertices().len(), 24);
        assert!(source.tolerance_changes().is_empty());
        assert_eq!(sewn.body().kind(), BodyKind::Solid);
        let counts = sewn.body().counts();
        assert_eq!((counts.faces, counts.edges, counts.vertices), (6, 12, 8));
        let seam_count = sewn
            .healing()
            .report()
            .changes
            .iter()
            .filter(|change| matches!(change, HealingChange::EdgesSewn { .. }))
            .count();
        assert_eq!(seam_count, 12);
        let exported =
            catallax_step::numerical::build(sewn.body(), &Default::default(), &mut dec).unwrap();
        exported.check_written(exported.bytes(), &mut dec).unwrap();
        std::fs::write(folder.join("sewn-solid.step"), exported.bytes()).unwrap();
        serde_json::json!({"exports":[{"name":"sewn-solid","kind":"Solid","faces":counts.faces,"edges":counts.edges,"vertices":counts.vertices,"bytes":exported.bytes().len()}],"seams":seam_count,"healing":sewn.healing().report(),"source_pair_checks":sewn.source_pair_checks(),"topology_repairs":sewn.topology_repairs()})
    });
    let mut report = result;
    report["work_units"] = serde_json::json!(work.used);
    std::fs::write(
        folder.join("surface-body-report.json"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    println!(
        "{}",
        serde_json::json!({"work_units":work.used,"seams":report["seams"],"exports":report["exports"]})
    );
}
