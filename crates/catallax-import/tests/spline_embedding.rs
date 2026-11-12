//! Exact embedding admission through the production STEP entry.
use catallax_import::{ImportRefusal, import_step_bodies, sha256_hex};
use catallax_numbers::with_operation_budget;
use catallax_recipe::StepId;
use catallax_topology::Model;

fn source(self_intersection: bool) -> Vec<u8> {
    fn add(lines: &mut Vec<String>, entity: String) -> usize {
        lines.push(entity);
        lines.len()
    }
    fn refs(ids: &[usize]) -> String {
        ids.iter()
            .map(|i| format!("#{i}"))
            .collect::<Vec<_>>()
            .join(",")
    }
    fn bound(lines: &mut Vec<String>, edges: &[usize], forward: bool, outer: bool) -> usize {
        let mut oriented = Vec::new();
        let mut edges = edges.to_vec();
        if !forward {
            edges.reverse();
        }
        for e in edges {
            oriented.push(add(
                lines,
                format!(
                    "ORIENTED_EDGE('',*,*,#{e},.{sense}.)",
                    sense = if forward { "T" } else { "F" }
                ),
            ));
        }
        let wire = add(lines, format!("EDGE_LOOP('',({}))", refs(&oriented)));
        add(
            lines,
            format!(
                "{}('',#{wire},.T.)",
                if outer {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                }
            ),
        )
    }
    let mut lines = Vec::new();
    let controls = if self_intersection {
        [
            ["12", "-6", "0"],
            ["-6", "7", "0"],
            ["-6", "-7", "0"],
            ["12", "6", "0"],
        ]
    } else {
        [
            ["1", "0", "0"],
            ["1", "0.55228475", "0"],
            ["0.55228475", "1", "0"],
            ["0", "1", "0"],
        ]
    };
    let mut control_ids = Vec::new();
    for p in controls {
        control_ids.push(add(
            &mut lines,
            format!("CARTESIAN_POINT('',({},{},{}))", p[0], p[1], p[2]),
        ));
    }
    let spline = add(
        &mut lines,
        format!(
            "B_SPLINE_CURVE_WITH_KNOTS('',3,({}),.UNSPECIFIED.,.F.,.F.,(4,4),(0.,1.),.UNSPECIFIED.)",
            refs(&control_ids)
        ),
    );
    let axis = add(&mut lines, "DIRECTION('',(0.,0.,1.))".into());
    let reference = add(&mut lines, "DIRECTION('',(1.,0.,0.))".into());
    let mut placements = Vec::new();
    let mut rings = Vec::new();
    for h in [0, 2] {
        let origin = add(&mut lines, format!("CARTESIAN_POINT('',(0.,0.,{h}.))"));
        let placement = add(
            &mut lines,
            format!("AXIS2_PLACEMENT_3D('',#{origin},#{axis},#{reference})"),
        );
        placements.push(placement);
        let circle = add(&mut lines, format!("CIRCLE('',#{placement},1.)"));
        let mut vertices = Vec::new();
        for (i, [x, y]) in [[1, 0], [0, 1], [-1, 0], [0, -1]].into_iter().enumerate() {
            let point = if h == 0 && self_intersection && i < 2 {
                control_ids[if i == 0 { 0 } else { 3 }]
            } else {
                add(&mut lines, format!("CARTESIAN_POINT('',({x}.,{y}.,{h}.))"))
            };
            vertices.push(add(&mut lines, format!("VERTEX_POINT('',#{point})")));
        }
        let mut edges = Vec::new();
        for i in 0..4 {
            let curve = if h == 0 && i == 0 { spline } else { circle };
            edges.push(add(
                &mut lines,
                format!(
                    "EDGE_CURVE('',#{},#{},#{curve},.T.)",
                    vertices[i],
                    vertices[(i + 1) % 4]
                ),
            ));
        }
        rings.push(edges);
    }
    let wall = add(
        &mut lines,
        format!("CYLINDRICAL_SURFACE('',#{},1.)", placements[0]),
    );
    let bottom = bound(&mut lines, &rings[0], true, true);
    let top = bound(&mut lines, &rings[1], false, false);
    let mut faces = vec![add(
        &mut lines,
        format!("ADVANCED_FACE('',(#{bottom},#{top}),#{wall},.T.)"),
    )];
    for i in 0..2 {
        let plane = add(&mut lines, format!("PLANE('',#{})", placements[i]));
        let outer = bound(&mut lines, &rings[i], i == 1, true);
        faces.push(add(
            &mut lines,
            format!(
                "ADVANCED_FACE('',(#{outer}),#{plane},.{sense}.)",
                sense = if i == 1 { "T" } else { "F" }
            ),
        ));
    }
    let shell = add(&mut lines, format!("CLOSED_SHELL('',({}))", refs(&faces)));
    add(
        &mut lines,
        format!("MANIFOLD_SOLID_BREP('synthetic',#{shell})"),
    );
    let unit = add(
        &mut lines,
        "(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))".into(),
    );
    add(
        &mut lines,
        format!(
            "(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{unit})) REPRESENTATION_CONTEXT('',''))"
        ),
    );
    let mut out = String::from("ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n");
    for (i, line) in lines.iter().enumerate() {
        out.push_str(&format!("#{}={line};\n", i + 1));
    }
    out.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    out.into_bytes()
}
fn refusal(bytes: &[u8]) -> ImportRefusal {
    let (result, work) = with_operation_budget(20_000_000, || {
        import_step_bodies(
            &mut Model::new(),
            bytes,
            &sha256_hex(bytes),
            &StepId::new("embedding").unwrap(),
        )
    });
    assert!(!work.exhausted);
    let error = result.unwrap_err();
    println!(
        "embedding-step {}",
        serde_json::json!({"input_sha256":sha256_hex(bytes),"work":work.used,"refusal":error})
    );
    error
}
#[test]
fn a_quarter_circle_reaches_recorded_curved_incidence_repair() {
    let bytes = source(false);
    let mut model = Model::new();
    let (result, work) = with_operation_budget(20_000_000, || {
        import_step_bodies(
            &mut model,
            &bytes,
            &sha256_hex(&bytes),
            &StepId::new("embedding").unwrap(),
        )
    });
    let bodies = result.unwrap();
    assert_eq!(bodies.len(), 1);
    assert!(!work.exhausted);
    let marks = model
        .body(bodies[0])
        .marks
        .iter()
        .filter(|m| m.code == "import.incidence-curved-repair")
        .collect::<Vec<_>>();
    assert_eq!(marks.len(), 1);
    let record: serde_json::Value = serde_json::from_str(&marks[0].detail).unwrap();
    assert_eq!(record["approved_tolerance_mm"], "0.01");
    assert!(marks[0].bound > catallax_numbers::Rational::zero());
    assert!(marks[0].bound <= catallax_numbers::Rational::parse("0.01").unwrap());
    assert!(
        record["definition"]
            .as_str()
            .unwrap()
            .contains("parameter-correspondence")
    );
    println!(
        "curved-step {}",
        serde_json::json!({"input_sha256":sha256_hex(&bytes),"work":work.used,"mark":record})
    );
    for limit in [0, work.used / 2, work.used - 1] {
        let (result, work) = with_operation_budget(limit, || {
            import_step_bodies(
                &mut Model::new(),
                &bytes,
                &sha256_hex(&bytes),
                &StepId::new("embedding").unwrap(),
            )
        });
        assert!(result.is_err());
        assert!(work.exhausted);
    }
}
#[test]
fn a_self_intersecting_source_keeps_the_typed_embedding_refusal() {
    let error = refusal(&source(true));
    assert!(
        matches!(error,ImportRefusal::Unsupported{ref type_name,..} if type_name.contains("spline curve embedding: import.spline-embedding-unproven")),
        "{error:?}"
    );
}

#[test]
fn curved_repair_refuses_a_small_neighbour_feature_and_an_excessive_curve_move() {
    let original = String::from_utf8(source(false)).unwrap();
    for (bytes, reason) in [
        (
            original.replace(",2.))", ",0.001))").into_bytes(),
            "permitted bound",
        ),
        (
            original
                .replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.DECI.,.METRE.)")
                .into_bytes(),
            "proven source gap exceeds the 10 um",
        ),
    ] {
        let error = refusal(&bytes);
        assert!(
            matches!(error, ImportRefusal::IncidenceGap { reason: ref detail, .. } if detail.contains(reason)),
            "{error:?}"
        );
    }
}

#[test]
fn curved_repair_requires_exact_endpoint_and_cap_plane_correspondence() {
    let original = String::from_utf8(source(false)).unwrap();
    for (bytes, reason) in [
        (
            original
                .replacen(
                    "CARTESIAN_POINT('',(1,0,0))",
                    "CARTESIAN_POINT('',(1.000001,0,0))",
                    1,
                )
                .into_bytes(),
            "endpoint correspondence",
        ),
        (
            original
                .replacen(
                    "CARTESIAN_POINT('',(1,0.55228475,0))",
                    "CARTESIAN_POINT('',(1,0.55228475,0.000001))",
                    1,
                )
                .into_bytes(),
            "source cap certificate",
        ),
    ] {
        assert_ne!(bytes, original.as_bytes());
        let error = refusal(&bytes);
        assert!(
            matches!(error, ImportRefusal::IncidenceGap { reason: ref detail, .. } if detail.contains(reason)),
            "{error:?}"
        );
    }
}
