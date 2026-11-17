//! N02 public polygonal tube and explicit endpoint-sewing controls.
use catallax_numbers::{
    Decider, NumericalCheck, NumericalRepair, Rational, Real, with_operation_budget,
};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Curve, Line, Plane, Point3, RationalPatch, SplineSupport,
    Surface, Vector3,
};
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
    perforated(
        [[0, 0], [4, 0], [4, 4], [0, 4]],
        &[[[1, 1], [1, 3], [3, 3], [3, 1]]],
    )
}

fn perforated(outer: [[i64; 2]; 4], holes: &[[[i64; 2]; 4]]) -> NumericalDefinition {
    let at = |p: [i64; 2], z| [p[0], p[1], z];
    let rings = || std::iter::once(&outer).chain(holes.iter());
    let mut polygons = vec![
        rings()
            .map(|ring| ring.map(|p| at(p, 1)).to_vec())
            .collect(),
        rings()
            .map(|ring| ring.iter().rev().map(|p| at(*p, 0)).collect())
            .collect(),
    ];
    for ring in rings() {
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

#[test]
fn non_binary_patch_boundary_is_sampled_as_a_finite_point() {
    let mut input = tube();
    let face = &input.faces[2];
    let points: Vec<_> = face.loops[0]
        .iter()
        .map(|trim| {
            let edge = &input.edges[trim.edge];
            let vertex = edge.vertices[usize::from(trim.sense == Sense::Reversed)];
            let [x, y, z] = input.vertices[vertex].position.clone().map(Real::from);
            Point3::new(x, y, z)
        })
        .collect();
    let hi = q("894427191/400000000");
    let mut dec = Decider::default();
    let carrier = BSplineSurface::new(
        [1, 1],
        [
            vec![q("0"), q("0"), hi.clone(), hi.clone()],
            vec![q("0"), q("0"), q("1"), q("1")],
        ],
        [2, 2],
        [0, 3, 1, 2].map(|i| points[i].clone()).to_vec(),
        vec![Real::one(); 4],
        &mut dec,
    )
    .unwrap();
    input.faces[2].surface = Surface::BSpline(SplineSupport::RationalPatch(
        RationalPatch::new(carrier, &mut dec).unwrap(),
    ));
    let corners = [
        [q("0"), q("0")],
        [hi.clone(), q("0")],
        [hi.clone(), q("1")],
        [q("0"), q("1")],
    ];
    for (i, trim) in input.faces[2].loops[0].iter_mut().enumerate() {
        trim.uv = [corners[i].clone(), corners[(i + 1) % 4].clone()];
    }
    numerical::build(input.clone(), options(), &mut dec).unwrap();
    let outside = &hi + &q("1/1000000");
    input.faces[2].loops[0][0].uv[1][0] = outside.clone();
    input.faces[2].loops[0][1].uv[0][0] = outside;
    let refusal = numerical::build(input, options(), &mut dec).unwrap_err();
    assert_eq!(refusal.code, "numerical.invalid-trim");
    assert_eq!(
        refusal.detail,
        "parameter trim leaves the retained patch domain"
    );
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

#[test]
fn continuation_retains_source_sewing_checks() {
    let mut input = tube();
    let pair = split(&mut input, "0.005");
    let source = numerical::sew(&input, &[pair], options(), &mut Decider::default()).unwrap();
    let edited = source
        .replace_surface(
            "f0",
            source.definition().faces[0].surface.clone(),
            &mut Decider::default(),
        )
        .unwrap();
    assert_eq!(edited.history().len(), 1);
    assert_eq!(edited.history()[0].evidence, *source.evidence());
    assert!(edited.history()[0].evidence.checks().iter().any(
        |c| matches!(c, NumericalCheck::Passed { name, .. } if name == "sewing.explicit-vertex-pair")
    ));
    assert_eq!(
        edited.source_definition().vertices.len(),
        input.vertices.len()
    );
    assert!(!edited.evidence().checks().iter().any(
        |c| matches!(c, NumericalCheck::Passed { name, .. } if name == "sewing.explicit-vertex-pair")
    ));
}

#[test]
fn curved_hole_sampling_does_not_claim_containment_between_nodes() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut dec = Decider::default();
        let mut input = tube();
        let edge_index = input.faces[0].loops[1][0].edge;
        let edge = &input.edges[edge_index];
        let ends = edge
            .vertices
            .map(|index| input.vertices[index].position.clone());

        // The excursion vanishes at all nine checked nodes but leaves the
        // outer rectangle at 1/16. Controls come from its full polynomial.
        let mut polynomial = vec![q("1")];
        for node in 0..=8 {
            let root = Rational::from_ratio(node, 8).unwrap();
            let mut next = vec![q("0"); polynomial.len() + 1];
            for (power, coefficient) in polynomial.iter().enumerate() {
                next[power] = &next[power] - &(coefficient * &root);
                next[power + 1] = &next[power + 1] + coefficient;
            }
            polynomial = next;
        }
        let at = q("1/16");
        let value = polynomial
            .iter()
            .rev()
            .fold(q("0"), |sum, coefficient| &(&sum * &at) + coefficient);
        let scale = -(&q("2") / &value);
        let degree = polynomial.len() - 1;
        let mut controls = Vec::new();
        for k in 0..=degree {
            let share = Rational::from_ratio(k as i64, degree as i64).unwrap();
            let mut position: [Rational; 3] = std::array::from_fn(|axis| {
                &ends[0][axis] + &(&share * &(&ends[1][axis] - &ends[0][axis]))
            });
            let mut factor = q("1");
            let mut bump = q("0");
            for (power, coefficient) in polynomial.iter().enumerate().take(k + 1) {
                if power > 0 {
                    factor = &factor
                        * &Rational::from_ratio(
                            (k - power + 1) as i64,
                            (degree - power + 1) as i64,
                        )
                        .unwrap();
                }
                bump = &bump + &(coefficient * &factor);
            }
            position[0] = &position[0] + &(&scale * &bump);
            let [x, y, z] = position;
            controls.push(Point3::new(x.into(), y.into(), z.into()));
        }
        let curve = BSplineCurve::new(
            degree,
            std::iter::repeat_n(q("0"), degree + 1)
                .chain(std::iter::repeat_n(q("1"), degree + 1))
                .collect(),
            controls,
            vec![Real::one(); degree + 1],
            &mut dec,
        )
        .unwrap();
        let witness = curve.point_at(&Real::from_rational(at), &mut dec).unwrap();
        assert_eq!(witness.x.as_rational(), Some(q("-1")));
        input.edges[edge_index].curve = Curve::BSpline(curve);
        let body = numerical::build(input, options(), &mut dec).unwrap();
        let checks = body.evidence().checks();
        assert!(
            !checks.iter().any(|check| matches!(check,
                NumericalCheck::Passed { name, scope }
                    if name == "trim.hole-contained-disjoint" && scope == "f0"
            )),
            "the actual curved hole leaves the outer loop between sample nodes"
        );
        assert!(checks.iter().any(|check| matches!(check,
            NumericalCheck::Passed { name, scope }
                if name == "trim.sampled-hole-contained-disjoint" && scope == "f0"
        )));
        assert!(checks.iter().any(|check| matches!(check,
            NumericalCheck::Passed { name, scope }
                if name == "trim.hole-contained-disjoint" && scope == "f1"
        )));
    });
    assert!(!work.exhausted, "{} work units", work.used);
}

#[path = "numerical_trim/native.rs"]
mod native;

#[test]
fn multiple_polygonal_holes_must_be_disjoint_and_not_nested() {
    let outer = [[0, 0], [8, 0], [8, 8], [0, 8]];
    let first = [[1, 1], [1, 3], [3, 3], [3, 1]];
    let second = [[5, 1], [5, 3], [7, 3], [7, 1]];
    let body = numerical::build(
        perforated(outer, &[first, second]),
        options(),
        &mut Decider::default(),
    )
    .unwrap();
    let counts = body.counts();
    assert_eq!((counts.faces, counts.edges, counts.vertices), (14, 36, 24));
    for holes in [
        [first, [[3, 1], [3, 3], [5, 3], [5, 1]]],
        [first, [[2, 1], [2, 3], [4, 3], [4, 1]]],
        [
            [[1, 1], [1, 5], [5, 5], [5, 1]],
            [[2, 2], [2, 3], [3, 3], [3, 2]],
        ],
    ] {
        let error = numerical::build(
            perforated(outer, &holes),
            options(),
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(error.code, "numerical.invalid-trim");
    }
}
