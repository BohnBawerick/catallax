//! Imported source preview, shared boundaries, trimming and resource limits.

use catallax_import::{Surface, read_step};
use catallax_mesh::{PreviewOptions, preview_step};
use catallax_numbers::{Decider, with_operation_budget};
use std::collections::{BTreeMap, BTreeSet};

const SOURCE: &[u8] = include_bytes!("../../catallax-import/tests/fixtures/spline-roof.step");

#[test]
fn previews_a_weighted_imported_face_with_shared_edges_and_source_picking() {
    let options = PreviewOptions::default();
    let first = preview_step(SOURCE, &options, &mut Decider::default()).unwrap();
    let second = preview_step(SOURCE, &options, &mut Decider::default()).unwrap();
    assert_eq!(first.mesh, second.mesh);
    let roof = first.source.solids[0]
        .faces
        .iter()
        .find(|f| matches!(f.surface, Surface::BSpline(_)))
        .unwrap();
    let region = first
        .mesh
        .faces
        .iter()
        .find(|f| f.face == roof.name)
        .unwrap();
    let center = first
        .mesh
        .positions
        .iter()
        .position(|p| {
            (p[0] - 0.5).abs() < 1e-9 && (p[1] - 0.5).abs() < 1e-9 && (p[2] - 1.1).abs() < 1e-9
        })
        .unwrap();
    assert!(
        first.mesh.triangles[region.triangles.clone()]
            .iter()
            .any(|t| t.vertices.contains(&center))
    );
    let mut uses = BTreeMap::<_, BTreeSet<_>>::new();
    let mut counts = BTreeMap::new();
    for triangle in &first.mesh.triangles {
        let p = triangle.vertices.map(|i| first.mesh.positions[i]);
        let a: [f64; 3] = std::array::from_fn(|i| p[1][i] - p[0][i]);
        let b: [f64; 3] = std::array::from_fn(|i| p[2][i] - p[0][i]);
        let normal: [f64; 3] = std::array::from_fn(|i| {
            a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3]
        });
        let outward = (0..3)
            .map(|i| normal[i] * ((p[0][i] + p[1][i] + p[2][i]) / 3.0 - 0.5))
            .sum::<f64>();
        assert!(outward > 0.0);
        for i in 0..3 {
            let (a, b) = (triangle.vertices[i], triangle.vertices[(i + 1) % 3]);
            let key = (a.min(b), a.max(b));
            *counts.entry(key).or_insert(0usize) += 1;
            uses.entry(key).or_default().insert(triangle.face.as_str());
        }
    }
    assert!(counts.values().all(|&n| n == 2));
    for edge in &first.mesh.edges {
        let expected: BTreeSet<_> = first.source.solids[0]
            .faces
            .iter()
            .filter(|f| {
                f.loops
                    .iter()
                    .flat_map(|lp| &lp.edges)
                    .any(|e| e.edge == edge.edge)
            })
            .map(|f| f.name.as_str())
            .collect();
        for pair in edge.vertices.windows(2) {
            assert_eq!(
                uses[&(pair[0].min(pair[1]), pair[0].max(pair[1]))],
                expected
            );
        }
    }
}

fn flat_roof_with_hole() -> Vec<u8> {
    let source = std::str::from_utf8(SOURCE)
        .unwrap()
        .replace("1.25", "1.0")
        .replace("(1.,2.,1.)", "(1.,1.,1.)");
    let parsed = read_step(source.as_bytes()).unwrap();
    let roof = parsed.solids[0]
        .faces
        .iter()
        .find(|f| matches!(f.surface, Surface::BSpline(_)))
        .unwrap();
    let mut records = Vec::new();
    let mut add = |s: String| {
        let id = 10000 + records.len();
        records.push(format!("#{id}={s};\n"));
        id
    };
    let coordinates = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]];
    let mut points = Vec::new();
    let mut vertices = Vec::new();
    for [x, y] in coordinates {
        let p = add(format!("CARTESIAN_POINT('',({x},{y},1.))"));
        points.push(p);
        vertices.push(add(format!("VERTEX_POINT('',#{p})")));
    }
    let mut edges = Vec::new();
    for i in 0..4 {
        let j = (i + 1) % 4;
        let direction = add(format!(
            "DIRECTION('',({},{},0.))",
            coordinates[j][0] - coordinates[i][0],
            coordinates[j][1] - coordinates[i][1]
        ));
        let vector = add(format!("VECTOR('',#{direction},1.)"));
        let line = add(format!("LINE('',#{},#{vector})", points[i]));
        let edge = add(format!(
            "EDGE_CURVE('',#{},#{},#{line},.T.)",
            vertices[i], vertices[j]
        ));
        edges.push(add(format!("ORIENTED_EDGE('',*,*,#{edge},.T.)")));
    }
    let lp = add(format!(
        "EDGE_LOOP('',({}))",
        edges
            .iter()
            .map(|i| format!("#{i}"))
            .collect::<Vec<_>>()
            .join(",")
    ));
    let hole = add(format!("FACE_BOUND('',#{lp},.F.)"));
    let line = source
        .lines()
        .find(|l| l.starts_with(&format!("{}=", roof.name)))
        .unwrap();
    let changed = line.replacen("),#", &format!(",#{hole}),#"), 1);
    let mut source = source.replace(line, &changed);
    source.insert_str(source.rfind("ENDSEC;").unwrap(), &records.concat());
    source.into_bytes()
}

#[test]
fn trims_a_hole_and_respects_reversed_bound_orientation() {
    let source = flat_roof_with_hole();
    let preview =
        preview_step(&source, &PreviewOptions::default(), &mut Decider::default()).unwrap();
    let roof = preview.source.solids[0]
        .faces
        .iter()
        .find(|f| matches!(f.surface, Surface::BSpline(_)))
        .unwrap();
    let mut area = 0.0;
    for triangle in preview
        .mesh
        .triangles
        .iter()
        .filter(|t| t.face == roof.name)
    {
        let p = triangle.vertices.map(|i| preview.mesh.positions[i]);
        let center: [f64; 2] = std::array::from_fn(|i| p.iter().map(|p| p[i] / 3.0).sum());
        assert!(!(center[0] > 0.25 && center[0] < 0.75 && center[1] > 0.25 && center[1] < 0.75));
        area += ((p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
            - (p[1][1] - p[0][1]) * (p[2][0] - p[0][0]))
            * 0.5;
    }
    assert!((area - 0.75).abs() < 1e-10);
}

#[test]
fn honors_work_and_vertex_limits_without_a_partial_preview() {
    let options = PreviewOptions {
        max_vertices: 8,
        ..PreviewOptions::default()
    };
    assert!(preview_step(SOURCE, &options, &mut Decider::default()).is_err());
    let (result, work) = with_operation_budget(10, || {
        preview_step(SOURCE, &PreviewOptions::default(), &mut Decider::default())
    });
    assert!(result.is_err());
    assert!(work.exhausted && work.used <= 10);
}

#[test]
fn preview_does_not_replace_a_declared_line_by_its_endpoint_chord() {
    // #111 is above the roof. The declared line then misses its original
    // endpoints, while the old endpoint-only chord would still close the box.
    let source = std::str::from_utf8(SOURCE)
        .unwrap()
        .replace("#21=LINE('',#1,#20);", "#21=LINE('',#111,#20);");
    read_step(source.as_bytes()).unwrap();
    let refusal = preview_step(
        source.as_bytes(),
        &PreviewOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "mesh.import-preview-unresolved");
}
