//! N02 public polygonal tube and explicit endpoint-sewing controls.
use catallax_numbers::{Decider, NumericalCheck, NumericalRepair, Rational, with_operation_budget};
use catallax_shapes::{Curve, Line, Plane, Point3, Surface, Vector3};
use catallax_topology::{
    Sense,
    numerical::{
        self, NumericalDefinition, NumericalEdge, NumericalFace, NumericalOptions, NumericalTrim,
        NumericalVertex,
    },
};
use std::collections::BTreeMap;
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn point(p: [i64; 3]) -> Point3 {
    Point3::new(p[0].into(), p[1].into(), p[2].into())
}
fn vector(p: [i64; 3]) -> Vector3 {
    Vector3::new(p[0].into(), p[1].into(), p[2].into())
}
fn sub(a: [i64; 3], b: [i64; 3]) -> [i64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [i64; 3], b: [i64; 3]) -> [i64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: [i64; 3]) -> [i64; 3] {
    a.map(|x| x.signum())
}
fn dot(a: [i64; 3], b: [i64; 3]) -> i64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

fn tube() -> NumericalDefinition {
    let outer = [[0, 0], [4, 0], [4, 4], [0, 4]];
    let hole = [[1, 1], [1, 3], [3, 3], [3, 1]];
    let at = |p: [i64; 2], z| [p[0], p[1], z];
    let mut polygons = vec![
        vec![
            outer.map(|p| at(p, 1)).to_vec(),
            hole.map(|p| at(p, 1)).to_vec(),
        ],
        vec![
            outer.into_iter().rev().map(|p| at(p, 0)).collect(),
            hole.into_iter().rev().map(|p| at(p, 0)).collect(),
        ],
    ];
    for ring in [outer, hole] {
        for i in 0..4 {
            let (a, b) = (ring[i], ring[(i + 1) % 4]);
            polygons.push(vec![vec![at(a, 0), at(b, 0), at(b, 1), at(a, 1)]]);
        }
    }
    let mut result = NumericalDefinition {
        name: "square-through-hole".into(),
        vertices: vec![],
        edges: vec![],
        faces: vec![],
    };
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    let mut dec = Decider::default();
    for (f, loops) in polygons.into_iter().enumerate() {
        let p = &loops[0];
        let origin = p[0];
        let x = unit(sub(p[1], origin));
        let normal = unit(cross(x, sub(p[2], p[1])));
        let y = cross(normal, x);
        let surface =
            Surface::Plane(Plane::new(point(origin), vector(normal), vector(x), &mut dec).unwrap());
        let mut trims = vec![];
        for polygon in loops {
            let mut boundary = vec![];
            for i in 0..polygon.len() {
                let ends = [polygon[i], polygon[(i + 1) % polygon.len()]];
                let ids = ends.map(|p| {
                    *vertices.entry(p).or_insert_with(|| {
                        let index = result.vertices.len();
                        result.vertices.push(NumericalVertex {
                            key: format!("v{index}"),
                            position: p.map(Rational::from),
                        });
                        index
                    })
                });
                let mut key = ids;
                key.sort();
                let edge = *edges.entry(key).or_insert_with(|| {
                    let index = result.edges.len();
                    result.edges.push(NumericalEdge {
                        key: format!("e{index}"),
                        curve: Curve::Line(
                            Line::new(point(ends[0]), vector(sub(ends[1], ends[0])), &mut dec)
                                .unwrap(),
                        ),
                        vertices: ids,
                        parameters: [q("0"), q("1")],
                    });
                    index
                });
                boundary.push(NumericalTrim {
                    edge,
                    sense: if result.edges[edge].vertices == ids {
                        Sense::Same
                    } else {
                        Sense::Reversed
                    },
                    uv: ends.map(|p| {
                        [
                            Rational::from(dot(sub(p, origin), x)),
                            Rational::from(dot(sub(p, origin), y)),
                        ]
                    }),
                });
            }
            trims.push(boundary);
        }
        result.faces.push(NumericalFace {
            key: format!("f{f}"),
            surface,
            sense: Sense::Same,
            loops: trims,
        });
    }
    result
}
fn record(name: &str, text: String) {
    if let Ok(directory) = std::env::var("CATALLAX_N02_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&directory).join(format!("{name}.txt")),
            text,
        )
        .unwrap();
    }
}
fn options() -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q("0.01"),
    }
}
fn split(input: &mut NumericalDefinition, gap: &str) -> [usize; 2] {
    let keeper = input.edges[0].vertices[0];
    let removed = input.vertices.len();
    let mut vertex = input.vertices[keeper].clone();
    vertex.key = "split-endpoint".into();
    vertex.position[0] = vertex.position[0].clone() + q(gap);
    input.vertices.push(vertex);
    input.edges[0].vertices[0] = removed;
    [keeper, removed]
}
#[test]
fn two_loop_material_and_near_gap_sewing_stay_numerical() {
    let original = tube();
    let body = numerical::build(original.clone(), options(), &mut Decider::default()).unwrap();
    let c = body.counts();
    assert_eq!((c.vertices, c.edges, c.faces, c.loops), (16, 24, 10, 12));
    // Independent footprint: outer area 16 minus hole area 4, unit height.
    let signed_area = |boundary: &[NumericalTrim]| {
        boundary.iter().fold(q("0"), |sum, t| {
            sum + (&t.uv[0][0] * &t.uv[1][1] - &t.uv[0][1] * &t.uv[1][0])
        }) / q("2")
    };
    assert_eq!(signed_area(&body.definition().faces[0].loops[0]), q("16"));
    assert_eq!(signed_area(&body.definition().faces[0].loops[1]), q("-4"));
    let Surface::Plane(top) = &body.definition().faces[0].surface else {
        panic!("plane")
    };
    assert_eq!(top.normal.z.as_rational(), Some(q("1")));
    for face in &body.definition().faces[6..] {
        let Surface::Plane(plane) = &face.surface else {
            panic!("plane")
        };
        let toward_hole = Point3::new(2.into(), 2.into(), q("1/2").into()).minus(&plane.origin);
        assert!(toward_hole.dot(&plane.normal).as_rational().unwrap() > q("0"));
    }
    assert!(
        body.exact_topology_faults()
            .iter()
            .all(|f| f.code == "topology.missing-incidence")
    );
    assert!(!body.exact_topology_faults().is_empty());
    for gap in ["0", "0.005"] {
        let mut input = original.clone();
        let pair = split(&mut input, gap);
        let before = format!("{input:?}");
        assert!(numerical::build(input.clone(), options(), &mut Decider::default()).is_err());
        let (answer, work) = with_operation_budget(20_000_000, || {
            numerical::sew(&input, &[pair], options(), &mut Decider::default())
        });
        let sewn = answer.unwrap();
        assert!(!work.exhausted);
        assert_eq!(format!("{input:?}"), before);
        assert_eq!(sewn.counts(), body.counts());
        assert_eq!(sewn.evidence().repairs().len(), usize::from(gap != "0"));
        if gap != "0" {
            let NumericalRepair::VertexMoved {
                vertex,
                before,
                after,
                distance_upper_mm,
            } = &sewn.evidence().repairs()[0];
            assert_eq!(vertex, "split-endpoint");
            assert_ne!(before, after);
            assert!(*distance_upper_mm >= q(gap) && *distance_upper_mm <= q("0.01"));
        }
        assert!(sewn.evidence().checks().iter().any(
            |c| matches!(c,NumericalCheck::Passed{name,..} if name=="sewing.explicit-vertex-pair")
        ));
        if gap != "0" {
            record(
                "near-gap",
                format!(
                    "N02 retained NumericalBody definition and NumericalEvidence (Rust debug, not a wire schema).\nResult category: numerical; no complete body bound or STEP authority.\n{:?}\n{:?}\n",
                    sewn.definition(),
                    sewn.evidence()
                ),
            );
        }
        println!(
            "N02 gap={gap} result=numerical work={} exhausted={} repairs={}",
            work.used,
            work.exhausted,
            sewn.evidence().repairs().len()
        );
    }
}
#[test]
fn excessive_gaps_bad_pairs_and_crossing_or_outside_holes_preserve_input() {
    let mut input = tube();
    let pair = split(&mut input, "0.02");
    let before = format!("{input:?}");
    let e = numerical::sew(&input, &[pair], options(), &mut Decider::default()).unwrap_err();
    assert_eq!(e.code, "numerical.repair-tolerance-exceeded");
    assert_eq!(format!("{input:?}"), before);
    record("excessive-gap-refusal", format!("{e:#?}\n"));
    let mut valid = tube();
    let valid_pair = split(&mut valid, "0");
    for pairs in [
        vec![valid_pair, valid_pair],
        vec![[0, 0]],
        vec![[0, usize::MAX]],
    ] {
        let error = numerical::sew(&valid, &pairs, options(), &mut Decider::default()).unwrap_err();
        assert_eq!(error.code, "numerical.invalid-definition");
    }
    for mode in [
        "crossing",
        "touching",
        "outside",
        "reversed",
        "self-crossing",
    ] {
        let mut input = tube();
        let hole = &mut input.faces[0].loops[1];
        if mode == "reversed" {
            for trim in hole {
                for point in &mut trim.uv {
                    point[0] = q("4") - point[0].clone();
                }
            }
        } else {
            for trim in hole {
                for p in &mut trim.uv {
                    p[0] = p[0].clone()
                        + if mode == "crossing" {
                            q("2")
                        } else if mode == "touching" {
                            q("1")
                        } else if mode == "outside" {
                            q("5")
                        } else {
                            q("0")
                        };
                }
            }
            if mode == "self-crossing" {
                let loop_ = &mut input.faces[0].loops[1];
                let points = [["1", "1"], ["3", "3"], ["1", "3"], ["3", "1"]].map(|p| p.map(q));
                for (i, t) in loop_.iter_mut().enumerate() {
                    t.uv = [points[i].clone(), points[(i + 1) % 4].clone()];
                }
            }
        }
        let before = format!("{input:?}");
        let e = numerical::sew(&input, &[], options(), &mut Decider::default()).unwrap_err();
        assert_eq!(e.code, "numerical.invalid-trim", "{mode}: {e:?}");
        assert_eq!(format!("{input:?}"), before);
        record(&format!("{mode}-refusal"), format!("{e:#?}\n"));
        println!("N02 control={mode} refusal={}", e.code);
    }
}
#[test]
fn sewing_and_admission_share_the_callers_allowance() {
    let mut input = tube();
    let pair = split(&mut input, "0.005");
    let (answer, work) = with_operation_budget(1, || {
        numerical::sew(&input, &[pair], options(), &mut Decider::default())
    });
    assert!(answer.is_err());
    assert!(work.exhausted);
}

#[test]
fn sewing_outer_and_hole_endpoints_refuses_without_mutating_input() {
    let mut input = tube();
    let scale = q("0.001");
    let real_scale = scale.clone().into();
    for vertex in &mut input.vertices {
        for coordinate in &mut vertex.position {
            *coordinate = &*coordinate * &scale;
        }
    }
    for edge in &mut input.edges {
        let Curve::Line(line) = &mut edge.curve else {
            panic!("line")
        };
        for coordinate in [
            &mut line.origin.x,
            &mut line.origin.y,
            &mut line.origin.z,
            &mut line.dir.x,
            &mut line.dir.y,
            &mut line.dir.z,
        ] {
            *coordinate = &*coordinate * &real_scale;
        }
    }
    for face in &mut input.faces {
        let Surface::Plane(plane) = &mut face.surface else {
            panic!("plane")
        };
        for coordinate in [
            &mut plane.origin.x,
            &mut plane.origin.y,
            &mut plane.origin.z,
        ] {
            *coordinate = &*coordinate * &real_scale;
        }
        for trim in face.loops.iter_mut().flatten() {
            for point in &mut trim.uv {
                for coordinate in point {
                    *coordinate = &*coordinate * &scale;
                }
            }
        }
    }
    let find = |position: [&str; 3]| {
        input
            .vertices
            .iter()
            .position(|v| v.position == position.map(q))
            .unwrap()
    };
    let pairs = [
        [find(["0", "0", "0.001"]), find(["0.001", "0.001", "0.001"])],
        [
            find(["0.004", "0", "0.001"]),
            find(["0.003", "0.001", "0.001"]),
        ],
    ];
    for [keeper, removed] in pairs {
        let squared = (0..3).fold(q("0"), |sum, i| {
            let delta = &input.vertices[keeper].position[i] - &input.vertices[removed].position[i];
            sum + &delta * &delta
        });
        assert_eq!(squared, q("0.000002"));
        assert!(squared < (&options().tolerance_mm * &options().tolerance_mm));
    }
    let original = numerical::build(input.clone(), options(), &mut Decider::default()).unwrap();
    let count = original.counts();
    assert_eq!(
        count.vertices as i64 - count.edges as i64 + count.faces as i64 - 2,
        0
    );
    let remap = |v| {
        pairs
            .iter()
            .find(|pair| pair[1] == v)
            .map_or(v, |pair| pair[0])
    };
    assert!(
        input
            .edges
            .iter()
            .all(|e| remap(e.vertices[0]) != remap(e.vertices[1]))
    );
    let before = format!("{input:?}");
    let (answer, work) = with_operation_budget(20_000_000, || {
        numerical::sew(&input, &pairs, options(), &mut Decider::default())
    });
    assert_eq!(format!("{input:?}"), before);
    assert!(!work.exhausted);
    if let Ok(body) = &answer {
        let count = body.counts();
        println!(
            "N02 P2 observed=accepted vertices={} edges={} faces={} euler={} work={}",
            count.vertices,
            count.edges,
            count.faces,
            count.vertices as i64 - count.edges as i64 + count.faces as i64 - 2,
            work.used
        );
    }
    assert!(
        answer.is_err(),
        "outer and hole loops must not share walked endpoint IDs after sewing"
    );
    let refusal = answer.unwrap_err();
    assert_eq!(refusal.code, "numerical.invalid-trim");
    assert!(refusal.detail.contains("share an endpoint"));
    println!(
        "N02 P2 observed=refused code={} work={} input=unchanged",
        refusal.code, work.used
    );
    record("shared-endpoints-refusal", format!("{refusal:#?}\n"));
}
