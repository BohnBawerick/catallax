//! One-way USD BRep export of analytic-square-through-hole/1 only.
//! Run with OUT HEIGHT PLACEMENT SOURCE_COMMIT. HEIGHT is 2 or 3; PLACEMENT is
//! identity or placed. See docs/archive/PLAN-2026-11.md, Optional USD BRep and Omniverse exchange.
use std::{collections::BTreeMap, error::Error, fmt::Write as _, fs, path::Path};

use catallax_features::evaluate;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::{Curve, PointDef, Surface};
use catallax_topology::{Loop, Sense};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

type Point = [i64; 3];
type Result<T, E = Box<dyn Error>> = std::result::Result<T, E>;
fn require(ok: bool, message: &str) -> Result<()> {
    if ok { Ok(()) } else { Err(message.into()) }
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: Point, b: Point) -> i64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn axis(a: Point) -> Result<Point> {
    let length: i64 = a.iter().map(|x| x.abs()).sum();
    require(
        length > 0 && a.iter().filter(|x| **x != 0).count() == 1,
        "only coordinate-axis geometry is supported",
    )?;
    Ok(a.map(|x| x / length))
}
fn exact(p: &catallax_shapes::Point3) -> Result<Point> {
    let mut result = [0; 3];
    for (i, value) in [&p.x, &p.y, &p.z].into_iter().enumerate() {
        let rational = value.as_rational().ok_or("expected rational coordinate")?;
        result[i] = rational.to_string().parse()?;
        require(
            rational == Rational::from_i64(result[i]),
            "coordinate is not an exact integer",
        )?;
    }
    Ok(result)
}
fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}
fn recipe(height: i64) -> Value {
    let ring = |prefix: &str, points: [[i64; 2]; 4]| {
        points
            .into_iter()
            .enumerate()
            .map(|(i, p)| {
                json!({
                    "kind":"line", "id":format!("{prefix}{i}"), "to":p.map(|x| x.to_string())
                })
            })
            .collect::<Vec<_>>()
    };
    json!({"format":"catallax-recipe/3", "name":"analytic-square-through-hole/1",
        "steps":[{"id":"profile", "op":"sketch",
            "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
            "loops":[ring("outer",[[5,-4],[5,4],[-5,4],[-5,-4]]),
                     ring("inner",[[1,-1],[1,1],[-1,1],[-1,-1]])]},
            {"id":"solid","op":"extrude","profile":"profile","distance":height.to_string()}],
        "result":"solid"})
}
fn canonical_loops() -> Vec<Vec<Vec<usize>>> {
    let mut faces = vec![
        vec![vec![3, 2, 1, 0], vec![4, 5, 6, 7]],
        vec![vec![8, 9, 10, 11], vec![12, 15, 14, 13]],
    ];
    for i in 0..4 {
        let j = (i + 1) % 4;
        faces.push(vec![vec![i, j, j + 8, i + 8]]);
    }
    for i in 4..8 {
        let j = 4 + (i - 3) % 4;
        faces.push(vec![vec![i, i + 8, j + 8, j]]);
    }
    faces
}
fn cycle_equal(a: &[usize], b: &[usize]) -> bool {
    a.len() == b.len() && (0..a.len()).any(|s| (0..a.len()).all(|i| a[(s + i) % a.len()] == b[i]))
}
fn attr(out: &mut String, ty: &str, name: &str, values: Value) {
    write!(out, "        uniform {ty} {name} = [").unwrap();
    for (i, value) in values
        .as_array()
        .expect("array attribute")
        .iter()
        .enumerate()
    {
        if i != 0 {
            out.push_str(", ");
        }
        if let Some(tuple) = value.as_array() {
            out.push('(');
            for (j, component) in tuple.iter().enumerate() {
                if j != 0 {
                    out.push_str(", ");
                }
                write!(out, "{component}").unwrap();
            }
            out.push(')');
        } else {
            write!(out, "{value}").unwrap();
        }
    }
    out.push_str("]\n");
}

fn export(out: &Path, height: i64, placed: bool, source_commit: &str) -> Result<()> {
    require(matches!(height, 2 | 3), "height must be 2 or 3 mm")?;
    require(
        source_commit.len() == 40 && source_commit.bytes().all(|c| c.is_ascii_hexdigit()),
        "source commit must be a full SHA-1",
    )?;
    let recipe_bytes = serde_json::to_vec_pretty(&recipe(height))?;
    let built = evaluate(
        &Recipe::parse(std::str::from_utf8(&recipe_bytes)?)?,
        &mut Decider::default(),
    )?;
    let model = &built.model;
    let body = built.result;
    require(model.check(body).is_empty(), "source topology check failed")?;
    require(
        model.body(body).marks.is_empty(),
        "source carries non-exact decisions",
    )?;
    let counts = model.counts(body);
    require(
        (
            counts.lumps,
            counts.shells,
            counts.faces,
            counts.loops,
            counts.edges,
            counts.vertices,
        ) == (1, 1, 10, 12, 24, 16),
        "unexpected source incidence",
    )?;
    let mut points: Vec<Point> = vec![
        [-5, -4, 0],
        [5, -4, 0],
        [5, 4, 0],
        [-5, 4, 0],
        [-1, -1, 0],
        [1, -1, 0],
        [1, 1, 0],
        [-1, 1, 0],
    ];
    for i in 0..8 {
        let [x, y, _] = points[i];
        points.push([x, y, height]);
    }
    let mut vertices = BTreeMap::new();
    let mut vertex_names = vec![Value::Null; 16];
    for vertex in model.vertices_of_body(body) {
        let PointDef::Explicit(p) = model.geom().point(model.vertex(vertex).point) else {
            return Err("non-explicit source point".into());
        };
        let p = exact(p)?;
        let i = points
            .iter()
            .position(|q| *q == p)
            .ok_or("unexpected source vertex")?;
        require(vertex_names[i].is_null(), "coincident source vertices")?;
        vertices.insert(vertex, i);
        vertex_names[i] = serde_json::to_value(
            model
                .vertex_name(body, vertex)
                .ok_or("unnamed source vertex")?,
        )?;
    }
    let loops = canonical_loops();
    let mut face_names = vec![Value::Null; 10];
    let mut source_charts = vec![Value::Null; 10];
    for face in model.faces_of_body(body) {
        let source = model.face(face);
        let Surface::Plane(plane) = model.geom().surface(source.surface) else {
            return Err("nonplanar source face".into());
        };
        let paths: Vec<Vec<usize>> = source
            .loops
            .iter()
            .map(|lp| {
                let Loop::Coedges(cs) = model.loop_(*lp) else {
                    return Err("point loop".into());
                };
                cs.iter()
                    .map(|c| {
                        let c = model.coedge(*c);
                        let e = model.edge(c.edge);
                        let v = if c.sense == Sense::Same {
                            e.start
                        } else {
                            e.end
                        };
                        Ok(vertices[&v.ok_or("ring edge")?])
                    })
                    .collect::<Result<_>>()
            })
            .collect::<Result<_>>()?;
        // Loop order is not an outer/inner tag. Match verified oriented cycles.
        let index = loops
            .iter()
            .position(|candidate| {
                candidate.len() == paths.len()
                    && candidate
                        .iter()
                        .all(|p| paths.iter().any(|q| cycle_equal(p, q)))
            })
            .ok_or("source loops do not match the fixture boundary")?;
        require(face_names[index].is_null(), "duplicate source face")?;
        let p = &loops[index][0];
        let n = axis(cross(
            sub(points[p[1]], points[p[0]]),
            sub(points[p[2]], points[p[1]]),
        ))?;
        let pn = [&plane.normal.x, &plane.normal.y, &plane.normal.z];
        for (component, expected) in pn.into_iter().zip(n) {
            let r = component.as_rational().ok_or("irrational plane normal")?;
            let r = if source.sense == Sense::Same { r } else { -r };
            require(
                r.sign() == Rational::from_i64(expected).sign(),
                "source normal opposes exported chart",
            )?;
        }
        let source_origin = exact(&plane.origin)?;
        require(
            paths
                .iter()
                .flatten()
                .all(|i| dot(n, sub(points[*i], source_origin)) == 0),
            "source plane misses a boundary vertex",
        )?;
        face_names[index] = serde_json::to_value(model.face_origins(face))?;
        require(
            face_names[index].as_array().is_some_and(|a| !a.is_empty()),
            "unnamed source face",
        )?;
        source_charts[index] = json!({
            "source_origin":source_origin,
            "source_normal":sub(exact(&plane.origin.plus(&plane.normal))?, source_origin),
            "source_ref_direction":sub(exact(&plane.origin.plus(&plane.ref_dir))?, source_origin),
            "source_sense":format!("{:?}",source.sense),
            "export_chart":"outer first vertex, U first segment, V=N cross U"
        });
    }
    let mut edge_map = BTreeMap::new();
    for edge in model.edges_of_body(body) {
        let e = model.edge(edge);
        let Curve::Line(line) = model.geom().curve(e.curve) else {
            return Err("nonlinear source edge".into());
        };
        let a = vertices[&e.start.ok_or("ring edge")?];
        let b = vertices[&e.end.ok_or("ring edge")?];
        let origin = exact(&line.origin)?;
        let direction = sub(exact(&line.origin.plus(&line.dir))?, origin);
        require(
            direction != [0; 3]
                && [a, b]
                    .into_iter()
                    .all(|i| cross(sub(points[i], origin), direction) == [0; 3]),
            "source line misses an endpoint",
        )?;
        require(
            edge_map.insert((a.min(b), a.max(b)), edge).is_none(),
            "duplicate source edge",
        )?;
    }
    let edges: Vec<_> = edge_map.keys().map(|&(a, b)| [a, b]).collect();
    let edge_names: Vec<_> = edge_map
        .values()
        .map(|e| {
            serde_json::to_value(model.edge_name(body, *e).ok_or("unnamed source edge")?)
                .map_err(Into::into)
        })
        .collect::<Result<_>>()?;
    let mut origins = Vec::new();
    let mut normals = Vec::new();
    let mut refs = Vec::new();
    let mut ranges = Vec::new();
    let mut eu_edges = Vec::new();
    let mut eu_senses = Vec::new();
    let mut uv = Vec::new();
    let mut knots = Vec::new();
    let mut users = vec![Vec::new(); edges.len()];
    for face in &loops {
        let outer = &face[0];
        let o = points[outer[0]];
        let u = axis(sub(points[outer[1]], o))?;
        let n = axis(cross(u, sub(points[outer[2]], points[outer[1]])))?;
        let v = cross(n, u);
        origins.push(o);
        refs.push(u);
        normals.push(n);
        let project = |p: Point| [dot(sub(p, o), u), dot(sub(p, o), v)];
        let all: Vec<_> = face.iter().flatten().map(|i| project(points[*i])).collect();
        ranges.push([
            all.iter().map(|p| p[0]).min().unwrap(),
            all.iter().map(|p| p[1]).min().unwrap(),
        ]);
        ranges.push([
            all.iter().map(|p| p[0]).max().unwrap(),
            all.iter().map(|p| p[1]).max().unwrap(),
        ]);
        for lp in face {
            for i in 0..lp.len() {
                let a = lp[i];
                let b = lp[(i + 1) % lp.len()];
                let e = edges
                    .binary_search(&[a.min(b), a.max(b)])
                    .map_err(|_| "missing shared source edge")?;
                users[e].push(eu_edges.len());
                eu_edges.push(e);
                eu_senses.push(if a < b { "same" } else { "opposite" });
                uv.extend([project(points[a]), project(points[b])]);
                let length: i64 = sub(points[b], points[a]).iter().map(|x| x.abs()).sum();
                knots.extend([0, 0, length, length]);
            }
        }
    }
    let mut radial = vec![0; 48];
    for uses in users {
        require(
            uses.len() == 2 && eu_senses[uses[0]] != eu_senses[uses[1]],
            "nonmanifold source edge",
        )?;
        radial[uses[0]] = uses[1];
        radial[uses[1]] = uses[0];
    }
    let line_origins: Vec<_> = edges.iter().map(|e| points[e[0]]).collect();
    let line_dirs: Vec<_> = edges
        .iter()
        .map(|e| axis(sub(points[e[1]], points[e[0]])))
        .collect::<Result<_>>()?;
    let edge_ranges: Vec<_> = edges
        .iter()
        .flat_map(|e| {
            [
                0,
                sub(points[e[1]], points[e[0]])
                    .iter()
                    .map(|x| x.abs())
                    .sum::<i64>(),
            ]
        })
        .collect();
    let matrix = if placed {
        "( (0, 1, 0, 0), (-1, 0, 0, 0), (0, 0, 1, 0), (20, 30, 40, 1) )"
    } else {
        "( (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1) )"
    };
    let mut usd = format!(
        "#usda 1.0\n(\n    defaultPrim = \"Part\"\n    metersPerUnit = 0.001\n    upAxis = \"Z\"\n)\ndef Xform \"Part\" {{\n    matrix4d xformOp:transform = {matrix}\n    uniform token[] xformOpOrder = [\"xformOp:transform\"]\n    def BrepArray \"Solid\" (\n        prepend apiSchemas = [\"BrepSurfacePlaneAPI\", \"BrepPointAPI:vertexPoint\", \"BrepCurve3dLineAPI:edge3dLine\", \"BrepCurveUvNurbAPI\"]\n        customData = {{ string source = \"Catallax\" }}\n    ) {{\n        uniform token purpose = \"guide\"\n"
    );
    attr(
        &mut usd,
        "double[]",
        "brep:intersectTol3d",
        json!([0.000001]),
    );
    attr(
        &mut usd,
        "double3[]",
        "brep:extent",
        json!([[-5, -4, 0], [5, 4, height]]),
    );
    for (name, value) in [
        ("brep:regionCount", json!([2])),
        ("region:shellCount", json!([1, 1])),
        ("shell:faceuseCount", json!([10, 10])),
        ("shell:wireEdgeCount", json!([0, 0])),
        (
            "faceuse:faceIndex",
            json!((0..10).chain(0..10).collect::<Vec<_>>()),
        ),
        (
            "face:loopCount",
            json!(loops.iter().map(Vec::len).collect::<Vec<_>>()),
        ),
        ("loop:edgeuseCount", json!(vec![4; 12])),
        ("loop:vertexIndex", json!(vec![0; 12])),
        ("edgeuse:edgeIndex", json!(eu_edges)),
        ("edgeuse:nextRadialEUIndex", json!(radial)),
    ] {
        attr(&mut usd, "uint[]", name, value);
    }
    for (name, value) in [
        ("region:type", json!(["voidRegion", "solidRegion"])),
        ("shell:pointType", json!(["none", "none"])),
        (
            "faceuse:orientationType",
            json!(
                vec!["same"; 10]
                    .into_iter()
                    .chain(vec!["opposite"; 10])
                    .collect::<Vec<_>>()
            ),
        ),
        ("face:surfaceType", json!(vec!["BrepSurfacePlaneAPI"; 10])),
        ("face:trimType", json!(vec!["rectangular"; 10])),
        ("edgeuse:orientationType", json!(eu_senses)),
        (
            "edgeuse:thisRadialEntryType",
            json!(
                eu_senses
                    .iter()
                    .map(|s| if *s == "same" {
                        "bottomEntry"
                    } else {
                        "topEntry"
                    })
                    .collect::<Vec<_>>()
            ),
        ),
        ("edge:curveType", json!(vec!["BrepCurve3dLineAPI"; 24])),
        ("vertex:pointType", json!(vec!["BrepPointAPI"; 16])),
    ] {
        attr(&mut usd, "token[]", name, value);
    }
    attr(&mut usd, "double2[]", "face:range", json!(ranges));
    attr(&mut usd, "double[]", "edge:range", json!(edge_ranges));
    attr(&mut usd, "int2[]", "edge:vertexIndices", json!(edges));
    for (name, value) in [
        ("brep:surface:plane:origin", json!(origins)),
        ("brep:edge3dLine:curve3d:line:origin", json!(line_origins)),
        ("brep:vertexPoint:point:position", json!(points)),
    ] {
        attr(&mut usd, "point3d[]", name, value);
    }
    for (name, value) in [
        ("brep:surface:plane:axis", json!(normals)),
        ("brep:surface:plane:refDirection", json!(refs)),
        ("brep:edge3dLine:curve3d:line:direction", json!(line_dirs)),
    ] {
        attr(&mut usd, "vector3d[]", name, value);
    }
    attr(
        &mut usd,
        "double2[]",
        "brep:curveUv:nurb:controlVertices",
        json!(uv),
    );
    for suffix in ["vertexCount", "order"] {
        attr(
            &mut usd,
            "uint[]",
            &format!("brep:curveUv:nurb:{suffix}"),
            json!(vec![2; 48]),
        );
    }
    attr(
        &mut usd,
        "double[]",
        "brep:curveUv:nurb:knots",
        json!(knots),
    );
    attr(
        &mut usd,
        "double[]",
        "brep:curveUv:nurb:weights",
        json!(vec![1; 96]),
    );
    usd.push_str("    }\n    def Mesh \"Display\" {\n        uniform token subdivisionScheme = \"none\"\n        uniform token orientation = \"rightHanded\"\n");
    attr(&mut usd, "point3f[]", "points", json!(points));
    // Four nonoverlapping quads tile each annular cap without filling its hole.
    let mut mesh = Vec::new();
    for i in 0..4 {
        let j = (i + 1) % 4;
        mesh.push(vec![j, i, i + 4, j + 4]);
        mesh.push(vec![i + 8, j + 8, j + 12, i + 12]);
    }
    mesh.extend(loops.iter().skip(2).map(|f| f[0].clone()));
    attr(
        &mut usd,
        "int[]",
        "faceVertexCounts",
        json!(vec![4; mesh.len()]),
    );
    attr(
        &mut usd,
        "int[]",
        "faceVertexIndices",
        json!(mesh.iter().flatten().collect::<Vec<_>>()),
    );
    usd.push_str("    }\n}\n");
    let identities = json!({"faces":face_names,"edges":edge_names,"vertices":vertex_names});
    let bounds = |count: usize| {
        (0..count).map(|index|json!({"index":index,"approximation_mm":"0","encoding_mm":"0","basis":"exact integer endpoints and coordinate-axis carriers"})).collect::<Vec<_>>()
    };
    let unit_rounding =
        Rational::from_f64_exact(0.001).unwrap() - Rational::from_ratio(1, 1000).unwrap();
    let physical_scale_bound =
        unit_rounding.abs() * Rational::from_i64(if placed { 40 + height } else { 5 });
    let tolerance_rounding =
        Rational::from_f64_exact(0.000001).unwrap() - Rational::from_ratio(1, 1000000).unwrap();
    let metadata = json!({"format":"catallax-usd-golden/1","fixture":"analytic-square-through-hole/1","height_mm":height,"placement":if placed {"placed"} else {"identity"},
        "source":{"commit":source_commit,"commit_authority":"caller-supplied","producer_version":env!("CARGO_PKG_VERSION"),"schema_commit":"c9c5979f04cea8c38213004d4f04df0dcfb40d05"},
        "recipe":{"file":"recipe.json","sha256":digest(&recipe_bytes),"revision":format!("height-{height}")},
        "files":{"part.usda":digest(usd.as_bytes())},"identities":identities,"charts":source_charts,
        "units":{"nominal_meters_per_unit":"1/1000","encoded_meters_per_unit":0.001,"rounding_m_per_unit":unit_rounding.to_string(),"max_world_coordinate_error_m":physical_scale_bound.to_string()},
        "intersect_tolerance":{"nominal_stage_units":"1/1000000","encoded_stage_units":0.000001,"role":"interchange setting, not measured error or healing permission","rounding_stage_units":tolerance_rounding.to_string()},
        "local_bounds":{"faces":bounds(10),"edges":bounds(24),"vertices":bounds(16),"uv":bounds(48)},
        "proof":{"source_topology":"checked","source_geometry":"exact analytic construction","acceptance":"requires independent checker","consumer_interoperability":"unproved","healing":"none"},
        "preservation":{"geometry":"exact for this fixture","topology":"shared incidence retained","trims":"explicit linear UV records retained","names":"serialized source identities retained","recipe":"recipe.json retained","local_certificates":"derived zero bounds, no general certificate transport","units":"nominal rational retained; binary64 scale rounding separate","proof_provenance":"source checks recorded; consumer execution unproved"},
        "information_loss":["BrepArray alone omits recipe, persistent identities and local proof detail; retain metadata.json and recipe.json","No USD import, editing, periodic seams, curved carriers or lossless round trip is supported"]});
    // A new directory prevents stale companions from an earlier export.
    fs::create_dir(out)?;
    fs::write(out.join("part.usda"), usd)?;
    fs::write(out.join("recipe.json"), recipe_bytes)?;
    fs::write(
        out.join("metadata.json"),
        serde_json::to_vec_pretty(&metadata)?,
    )?;
    Ok(())
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    require(
        args.len() == 4,
        "usage: usd_brep OUT HEIGHT identity|placed SOURCE_COMMIT",
    )?;
    require(
        matches!(args[2].as_str(), "identity" | "placed"),
        "placement must be identity or placed",
    )?;
    export(
        Path::new(&args[0]),
        args[1].parse()?,
        args[2] == "placed",
        &args[3],
    )?;
    println!("Exported {}", args[0]);
    Ok(())
}

#[test]
fn golden_export_edit_preserves_persistent_identities() -> Result<()> {
    let retained = std::env::var_os("CATALLAX_USD_PROOF");
    let root = retained
        .as_ref()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::env::temp_dir().join(format!(
                "catallax-usd-{}-{}",
                std::process::id(),
                std::time::SystemTime::now()
                    .duration_since(std::time::UNIX_EPOCH)
                    .unwrap()
                    .as_nanos()
            ))
        });
    fs::create_dir_all(&root)?;
    let commit = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .output()?;
    require(commit.status.success(), "cannot resolve source commit")?;
    let commit = std::str::from_utf8(&commit.stdout)?.trim();
    export(&root.join("baseline"), 2, false, commit)?;
    export(&root.join("edited"), 3, false, commit)?;
    export(&root.join("placed"), 2, true, commit)?;
    let read_metadata = |name: &str| -> Result<Value> {
        Ok(serde_json::from_slice(&fs::read(
            root.join(name).join("metadata.json"),
        )?)?)
    };
    let before = read_metadata("baseline")?;
    let after = read_metadata("edited")?;
    let placed = read_metadata("placed")?;
    assert_eq!(before["identities"], after["identities"]);
    assert_eq!(before["identities"], placed["identities"]);
    assert_ne!(before["recipe"]["sha256"], after["recipe"]["sha256"]);
    assert!(export(&root.join("baseline"), 2, false, commit).is_err());
    if retained.is_none() {
        fs::remove_dir_all(&root)?;
    } else {
        println!("Retained export packages at {}", root.display());
    }
    Ok(())
}
