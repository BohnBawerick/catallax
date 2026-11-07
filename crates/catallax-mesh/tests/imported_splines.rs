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

fn concave_prism(notch_y: f64, reversed: bool) -> Vec<u8> {
    let mut records = Vec::new();
    let mut add = |s: String| {
        let id = records.len() + 1;
        records.push(format!("#{id}={s};\n"));
        id
    };
    let mut outline = [
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.5, 1.0],
        [0.5, notch_y],
        [0.0, notch_y],
    ];
    if reversed {
        outline.reverse();
    }
    let coordinates: Vec<_> = [0.0, 1.0]
        .into_iter()
        .flat_map(|z| outline.iter().map(move |&[x, y]| [x, y, z]))
        .collect();
    let points: Vec<_> = coordinates
        .iter()
        .map(|[x, y, z]| add(format!("CARTESIAN_POINT('',({x},{y},{z}))")))
        .collect();
    let vertices: Vec<_> = points
        .iter()
        .map(|p| add(format!("VERTEX_POINT('',#{p})")))
        .collect();
    let mut shared = BTreeMap::new();
    let mut faces = Vec::new();
    let caps = [0usize, 6].map(|offset| (offset..offset + 6).collect::<Vec<_>>());
    let sides = (0..6).map(|i| vec![i, (i + 1) % 6, (i + 1) % 6 + 6, i + 6]);
    for polygon in caps.into_iter().chain(sides) {
        let mut oriented = Vec::new();
        for i in 0..polygon.len() {
            let (a, b) = (polygon[i], polygon[(i + 1) % polygon.len()]);
            let key = (a.min(b), a.max(b));
            let edge = *shared.entry(key).or_insert_with(|| {
                let (a, b) = key;
                let [dx, dy, dz] = std::array::from_fn(|i| coordinates[b][i] - coordinates[a][i]);
                let direction = add(format!("DIRECTION('',({dx},{dy},{dz}))"));
                let vector = add(format!("VECTOR('',#{direction},1.)"));
                let line = add(format!("LINE('',#{},#{vector})", points[a]));
                add(format!(
                    "EDGE_CURVE('',#{},#{},#{line},.T.)",
                    vertices[a], vertices[b]
                ))
            });
            let sense = if a < b { "T" } else { "F" };
            oriented.push(add(format!("ORIENTED_EDGE('',*,*,#{edge},.{sense}.)")));
        }
        let refs = |ids: &[usize]| {
            ids.iter()
                .map(|id| format!("#{id}"))
                .collect::<Vec<_>>()
                .join(",")
        };
        let lp = add(format!("EDGE_LOOP('',({}))", refs(&oriented)));
        let bound = add(format!("FACE_OUTER_BOUND('',#{lp},.T.)"));
        let a = coordinates[polygon[0]];
        let b = coordinates[polygon[1]];
        let c = coordinates[*polygon.last().unwrap()];
        let u: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
        let v: [f64; 3] = std::array::from_fn(|i| c[i] - a[i]);
        let n: [f64; 3] = std::array::from_fn(|i| {
            u[(i + 1) % 3] * v[(i + 2) % 3] - u[(i + 2) % 3] * v[(i + 1) % 3]
        });
        // Cap charts keep the witness midpoint at (0.5, 0.5).
        let (n, u) = if polygon.len() == 6 {
            ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
        } else {
            (n, u)
        };
        let normal = add(format!("DIRECTION('',({},{},{}))", n[0], n[1], n[2]));
        let reference = add(format!("DIRECTION('',({},{},{}))", u[0], u[1], u[2]));
        let origin = if polygon.len() == 6 {
            add(format!("CARTESIAN_POINT('',(0.,0.,{}))", a[2]))
        } else {
            points[polygon[0]]
        };
        let frame = add(format!(
            "AXIS2_PLACEMENT_3D('',#{origin},#{normal},#{reference})"
        ));
        let plane = add(format!("PLANE('',#{frame})"));
        faces.push(add(format!("ADVANCED_FACE('',(#{bound}),#{plane},.T.)")));
    }
    let face_refs = faces
        .iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",");
    let shell = add(format!("CLOSED_SHELL('',({face_refs}))"));
    let solid = add(format!("MANIFOLD_SOLID_BREP('concave prism',#{shell})"));
    let unit = add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))".into());
    let context = add(format!(
        "(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{unit})) REPRESENTATION_CONTEXT('',''))"
    ));
    add(format!(
        "ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid}),#{context})"
    ));
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n{}ENDSEC;\nEND-ISO-10303-21;\n",
        records.concat()
    )
    .into_bytes()
}

fn assert_concave_prism_shared_edges(notch_y: f64) {
    for reversed in [false, true] {
        let source = concave_prism(notch_y, reversed);
        let preview =
            preview_step(&source, &PreviewOptions::default(), &mut Decider::default()).unwrap();
        for edge in &preview.mesh.edges {
            let incident: BTreeSet<_> = preview.source.solids[0]
                .faces
                .iter()
                .filter(|face| {
                    face.loops
                        .iter()
                        .flat_map(|lp| &lp.edges)
                        .any(|e| e.edge == edge.edge)
                })
                .map(|face| face.name.as_str())
                .collect();
            assert_eq!(incident.len(), 2);
            for pair in edge.vertices.windows(2) {
                for face in &incident {
                    let uses = preview
                        .mesh
                        .triangles
                        .iter()
                        .filter(|t| t.face == *face)
                        .filter(|t| t.vertices.contains(&pair[0]) && t.vertices.contains(&pair[1]))
                        .count();
                    assert_eq!(
                        uses, 1,
                        "edge {} segment {pair:?}, face {face}, notch {notch_y}, reversed {reversed}",
                        edge.edge
                    );
                }
            }
        }
    }
}

#[test]
fn rectangle_seed_at_concave_trim_vertex_preserves_shared_indices() {
    assert_concave_prism_shared_edges(0.5);
}

#[test]
fn rectangle_seed_on_concave_trim_segment_preserves_shared_indices() {
    assert_concave_prism_shared_edges(0.25);
}
