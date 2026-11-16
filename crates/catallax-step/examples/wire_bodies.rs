//! Exact open/closed spline wire bodies, written for healing-off STEP checks.
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::spline_tools::RationalSplineCurve;
use catallax_topology::{BodyKind, spline_wire};

fn curve(points: [[i64; 3]; 3], dec: &Decider) -> RationalSplineCurve {
    RationalSplineCurve::new(
        2,
        [0, 0, 0, 1, 1, 1]
            .into_iter()
            .map(Rational::from_i64)
            .collect(),
        points
            .into_iter()
            .map(|point| point.map(Rational::from_i64))
            .collect(),
        vec![Rational::one(); 3],
        dec,
    )
    .unwrap()
}

fn main() {
    let folder = std::path::PathBuf::from(std::env::args().nth(1).expect("output directory"));
    std::fs::create_dir_all(&folder).unwrap();
    let (result, work) = with_operation_budget(catallax_topology::numerical::DEFAULT_WORK, || {
        let mut dec = Decider::default();
        let a = curve([[0, 0, 0], [1, 1, 0], [2, 0, 0]], &dec);
        let b = curve([[2, 0, 0], [3, -1, 0], [4, 0, 0]], &dec);
        let closing = curve([[4, 0, 0], [2, -2, 1], [0, 0, 0]], &dec);
        let ring = RationalSplineCurve::new(
            3,
            [0, 0, 0, 0, 1, 1, 1, 1]
                .into_iter()
                .map(Rational::from_i64)
                .collect(),
            [[0, 0, 0], [2, 0, 0], [1, 2, 0], [0, 0, 0]]
                .into_iter()
                .map(|p| p.map(Rational::from_i64))
                .collect(),
            vec![Rational::one(); 4],
            &dec,
        )
        .unwrap();
        assert!(spline_wire(&[a.clone(), a.clone()], &mut dec).is_err());
        let mut report = Vec::new();
        for (name, pieces, expected_vertices) in [
            ("open-wire", vec![a.clone(), b.clone()], 3),
            ("closed-wire", vec![a, b, closing], 3),
            ("single-closed-wire", vec![ring], 1),
        ] {
            let (model, body) = spline_wire(&pieces, &mut dec)?;
            assert_eq!(model.body(body).kind, BodyKind::Wire);
            assert!(model.check(body).is_empty());
            let counts = model.counts(body);
            assert_eq!(counts.faces, 0);
            assert_eq!(counts.edges, pieces.len());
            assert_eq!(counts.vertices, expected_vertices);
            let mut bytes = Vec::new();
            let exported = catallax_step::write_step(
                &model,
                &[body],
                &Default::default(),
                &mut dec,
                &mut bytes,
            )?;
            assert_eq!(exported.value.solids, 0);
            assert!(String::from_utf8_lossy(&bytes).contains("EDGE_BASED_WIREFRAME_MODEL"));
            std::fs::write(folder.join(format!("{name}.step")), &bytes).unwrap();
            report.push(serde_json::json!({"name":name,"edges":counts.edges,"vertices":counts.vertices,"bytes":bytes.len(),"step_max_deviation_mm":exported.value.max_deviation.to_string()}));
        }
        Ok::<_, catallax_numbers::Refusal>(report)
    });
    let report = serde_json::json!({"exports":result.unwrap(),"work_units":work.used});
    std::fs::write(
        folder.join("wire-body-report.json"),
        serde_json::to_vec_pretty(&report).unwrap(),
    )
    .unwrap();
    println!("{report}");
}
