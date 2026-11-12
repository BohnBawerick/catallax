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
fn a_quarter_circle_reaches_measured_incidence_instead_of_an_embedding_refusal() {
    let ImportRefusal::IncidenceGap { gap, .. } = refusal(&source(false)) else {
        panic!("expected the next source incidence refusal")
    };
    assert_eq!(gap.step, "spline_curve_surface_incidence");
    assert_eq!(gap.coverage, "finite_source_feature");
    assert!(gap.converged);
}
#[test]
fn a_self_intersecting_source_keeps_the_typed_embedding_refusal() {
    let error = refusal(&source(true));
    assert!(
        matches!(error,ImportRefusal::Unsupported{ref type_name,..} if type_name.contains("spline curve embedding: import.spline-embedding-unproven")),
        "{error:?}"
    );
}
