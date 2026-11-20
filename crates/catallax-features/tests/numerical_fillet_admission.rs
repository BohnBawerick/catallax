//! Fillet admission preserves requested tolerance and untouched periodic seam uses.

use catallax_features::fillet::numerical as fillet;
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Curve, Line, Point3, RationalAtlas, SplineSupport, Surface,
    Vector3,
};
use catallax_topology::{
    Sense,
    numerical::{
        self, NumericalBody, NumericalDefinition, NumericalEdge, NumericalFace,
        NumericalNativeTrim, NumericalOptions, NumericalTrim, NumericalVertex,
    },
};
use serde_json::json;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn linear_uv(
    edge: &numerical::NumericalEdge,
    trim: &numerical::NumericalTrim,
    dec: &mut Decider,
) -> BSplineCurve {
    let mut uv = trim.uv.clone();
    if trim.sense == Sense::Reversed {
        uv.reverse();
    }
    let mut range = edge.parameters.clone();
    if range[0] > range[1] {
        range.reverse();
        uv.reverse();
    }
    BSplineCurve::new(
        1,
        vec![
            range[0].clone(),
            range[0].clone(),
            range[1].clone(),
            range[1].clone(),
        ],
        uv.map(|p| Point3::new(p[0].clone().into(), p[1].clone().into(), Real::zero()))
            .to_vec(),
        vec![Real::one(); 2],
        dec,
    )
    .unwrap()
}

fn plane_parameters(surface: &Surface, point: &Point3, dec: &Decider) -> [Rational; 2] {
    let point = [&point.x, &point.y, &point.z].map(|value| dec.enclose_at(value, 128).unwrap());
    surface
        .parameters_of(&point, dec)
        .unwrap()
        .unwrap()
        .map(|value| value.mid())
}

// A closed C1 cubic extrusion gives a real periodic support without requiring
// the plane/spline-only exact-body embedding to convert an analytic cylinder.
fn add_periodic_bore(
    definition: &mut NumericalDefinition,
    bindings: &mut Vec<NumericalNativeTrim>,
    dec: &mut Decider,
) {
    let xy = [
        [2, 0],
        [2, 1],
        [1, 2],
        [0, 2],
        [-1, 2],
        [-2, 1],
        [-2, 0],
        [-2, -1],
        [-1, -2],
        [0, -2],
        [1, -2],
        [2, -1],
        [2, 0],
    ];
    let point = |i: usize, z: i64| {
        Point3::new(
            (q("3/2") + Rational::from_ratio(xy[i][0], 4).unwrap()).into(),
            (q("3/2") + Rational::from_ratio(xy[i][1], 4).unwrap()).into(),
            z.into(),
        )
    };
    let knots: Vec<_> = [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4]
        .map(Rational::from_i64)
        .to_vec();
    let carrier = BSplineSurface::new(
        [3, 1],
        [knots.clone(), vec![q("0"), q("0"), q("3"), q("3")]],
        [13, 2],
        (0..13).flat_map(|i| [point(i, 0), point(i, 3)]).collect(),
        vec![Real::one(); 26],
        dec,
    )
    .unwrap();
    let [atlas]: [_; 1] = RationalAtlas::partition(carrier, [true, false], dec)
        .unwrap()
        .try_into()
        .unwrap();
    let mut arcs = [[0; 4]; 2];
    for (end, z) in [0, 3].into_iter().enumerate() {
        let vertices = definition.vertices.len();
        for i in 0..4 {
            let p = point(3 * i, z);
            definition.vertices.push(NumericalVertex {
                key: format!("bore/v/{end}/{i}"),
                position: [p.x, p.y, p.z].map(|value| value.as_rational().unwrap()),
            });
        }
        let controls: Vec<_> = (0..13).map(|i| point(i, z)).collect();
        let world = BSplineCurve::new(
            3,
            knots.clone(),
            controls.clone(),
            vec![Real::one(); 13],
            dec,
        )
        .unwrap();
        for i in 0..4 {
            arcs[end][i] = definition.edges.len();
            definition.edges.push(NumericalEdge {
                key: format!("bore/arc/{end}/{i}"),
                curve: Curve::BSpline(world.clone()),
                vertices: [vertices + i, vertices + (i + 1) % 4],
                parameters: [
                    Rational::from_i64(i as i64),
                    Rational::from_i64(i as i64 + 1),
                ],
            });
        }
        let cap = definition
            .faces
            .iter()
            .position(|face| {
                matches!(&face.surface, Surface::Plane(plane)
                if plane.normal.z.as_rational().is_some_and(|n| n != q("0"))
                    && plane.origin.z.as_rational() == Some(Rational::from_i64(z)))
            })
            .unwrap();
        let surface = &definition.faces[cap].surface;
        let uv_curve = BSplineCurve::new(
            3,
            knots.clone(),
            controls
                .iter()
                .map(|p| {
                    let [u, v] = plane_parameters(surface, p, dec);
                    Point3::new(u.into(), v.into(), Real::zero())
                })
                .collect(),
            vec![Real::one(); 13],
            dec,
        )
        .unwrap();
        let boundary = definition.faces[cap].loops.len();
        let mut row = Vec::new();
        for i in if end == 0 { [0, 1, 2, 3] } else { [3, 2, 1, 0] } {
            let sense = if end == 0 {
                Sense::Same
            } else {
                Sense::Reversed
            };
            let mut uv = [i, i + 1].map(|t| {
                let p = uv_curve
                    .point_at(&Real::from(Rational::from_i64(t as i64)), dec)
                    .unwrap();
                [p.x.as_rational().unwrap(), p.y.as_rational().unwrap()]
            });
            if sense == Sense::Reversed {
                uv.reverse();
            }
            let edge = arcs[end][i];
            bindings.push(NumericalNativeTrim {
                tag: catallax_topology::numerical::CopyTag::Independent,
                face: definition.faces[cap].key.clone(),
                boundary,
                use_index: row.len(),
                edge: definition.edges[edge].key.clone(),
                curve: uv_curve.clone(),
            });
            row.push(NumericalTrim {
                tag: catallax_topology::numerical::CopyTag::Independent,
                edge,
                sense,
                uv,
            });
        }
        definition.faces[cap].loops.push(row);
    }
    let seam = definition.edges.len();
    definition.edges.push(NumericalEdge {
        key: "bore/seam".into(),
        curve: Curve::Line(Line {
            origin: point(0, 0),
            dir: Vector3::new(0.into(), 0.into(), 1.into()),
        }),
        vertices: [
            definition.edges[arcs[0][0]].vertices[0],
            definition.edges[arcs[1][0]].vertices[0],
        ],
        parameters: [q("0"), q("3")],
    });
    let mut row = vec![NumericalTrim {
        tag: catallax_topology::numerical::CopyTag::Independent,
        edge: seam,
        sense: Sense::Same,
        uv: [[q("0"), q("0")], [q("0"), q("3")]],
    }];
    for i in 0..4 {
        row.push(NumericalTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            edge: arcs[1][i],
            sense: Sense::Same,
            uv: [
                [Rational::from_i64(i as i64), q("3")],
                [Rational::from_i64(i as i64 + 1), q("3")],
            ],
        });
    }
    row.push(NumericalTrim {
        tag: catallax_topology::numerical::CopyTag::Independent,
        edge: seam,
        sense: Sense::Reversed,
        uv: [[q("4"), q("3")], [q("4"), q("0")]],
    });
    for i in (0..4).rev() {
        row.push(NumericalTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            edge: arcs[0][i],
            sense: Sense::Reversed,
            uv: [
                [Rational::from_i64(i as i64 + 1), q("0")],
                [Rational::from_i64(i as i64), q("0")],
            ],
        });
    }
    for (use_index, trim) in row.iter().enumerate() {
        bindings.push(NumericalNativeTrim {
            tag: catallax_topology::numerical::CopyTag::Independent,
            face: "bore".into(),
            boundary: 0,
            use_index,
            edge: definition.edges[trim.edge].key.clone(),
            curve: linear_uv(&definition.edges[trim.edge], trim, dec),
        });
    }
    definition.faces.push(NumericalFace {
        key: "bore".into(),
        surface: Surface::BSpline(SplineSupport::RationalAtlas(atlas)),
        sense: Sense::Reversed,
        loops: vec![row],
    });
}

fn source(bore: bool, gap: bool, dec: &mut Decider) -> (NumericalBody, String) {
    let value = json!({"format":"catallax-recipe/5","name":"fillet admission controls","steps":[
        {"id":"stock","op":"block","size":["3","3","3"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}
    ],"result":"stock"});
    let recipe = Recipe::parse(&value.to_string()).unwrap();
    let built = catallax_features::evaluate(&recipe, dec).unwrap();
    let tolerance = if gap { q("1/1000") } else { q("1/10000000") };
    let (mut definition, mut bindings) =
        catallax_features::numerical_embedding::definition_from_exact_model(
            &built.model,
            built.result,
            recipe.name.clone(),
            &tolerance,
            dec,
        )
        .unwrap();
    if bore {
        add_periodic_bore(&mut definition, &mut bindings, dec);
    }
    let selected = definition
        .edges
        .iter()
        .position(|edge| {
            let [a, b] = edge.vertices.map(|i| &definition.vertices[i].position);
            a[0] == q("0") && b[0] == q("0") && a[1] == q("0") && b[1] == q("0") && a[2] != b[2]
        })
        .unwrap();
    let key = definition.edges[selected].key.clone();
    for face in &mut definition.faces {
        let side = face.loops.iter().flatten().any(|t| t.edge == selected);
        if side {
            let Surface::Plane(plane) = face.surface.clone() else {
                panic!("box sides are planar")
            };
            let points: Vec<_> = face
                .loops
                .iter()
                .flatten()
                .flat_map(|t| t.uv.iter())
                .collect();
            let bounds: [[Rational; 2]; 2] = std::array::from_fn(|axis| {
                [
                    points.iter().map(|p| p[axis].clone()).min().unwrap(),
                    points.iter().map(|p| p[axis].clone()).max().unwrap(),
                ]
            });
            let vertical = plane.normal.cross(&plane.ref_dir);
            let plane = &plane;
            let vertical = &vertical;
            let controls = bounds[0]
                .iter()
                .flat_map(|u| {
                    bounds[1].iter().map(move |v| {
                        plane
                            .origin
                            .plus(&plane.ref_dir.scaled(&Real::from(u.clone())))
                            .plus(&vertical.scaled(&Real::from(v.clone())))
                    })
                })
                .collect();
            let knots = bounds
                .each_ref()
                .map(|b| vec![b[0].clone(), b[0].clone(), b[1].clone(), b[1].clone()]);
            let carrier =
                BSplineSurface::new([1, 1], knots, [2, 2], controls, vec![Real::one(); 4], dec)
                    .unwrap();
            let [atlas]: [_; 1] = RationalAtlas::partition(carrier, [false; 2], dec)
                .unwrap()
                .try_into()
                .unwrap();
            face.surface = Surface::BSpline(SplineSupport::RationalAtlas(atlas));
        }
        if side {
            bindings.retain(|binding| binding.face != face.key);
            for (boundary, row) in face.loops.iter().enumerate() {
                for (use_index, trim) in row.iter().enumerate() {
                    let edge = &definition.edges[trim.edge];
                    bindings.push(NumericalNativeTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        face: face.key.clone(),
                        boundary,
                        use_index,
                        edge: edge.key.clone(),
                        curve: linear_uv(edge, trim, dec),
                    });
                }
            }
        }
    }
    if gap {
        let border = definition
            .edges
            .iter_mut()
            .find(|edge| {
                let [a, b] = edge.vertices.map(|i| &definition.vertices[i].position);
                a[1] == q("0") && b[1] == q("0") && a[2] == q("0") && b[2] == q("0") && a[0] != b[0]
            })
            .unwrap();
        let Curve::Line(line) = &mut border.curve else {
            panic!("box border is linear")
        };
        line.origin.z = &line.origin.z + &Real::from(q("1/10000"));
    }
    let body = numerical::build_with_native_trims(
        definition,
        bindings,
        NumericalOptions {
            tolerance_mm: tolerance,
        },
        dec,
    )
    .unwrap();
    (body, key)
}

#[test]
fn a_stricter_fillet_request_refuses_a_loose_source_border() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (exact, key) = source(false, false, &mut dec);
        let loose = exact
            .rebuild(
                exact.definition().clone(),
                NumericalOptions {
                    tolerance_mm: q("1/1000"),
                },
                &mut dec,
            )
            .unwrap();
        let tight =
            fillet::build(loose, &key, "round", q("1/10"), q("1/10000000"), &mut dec).unwrap();
        assert!(
            tight
                .body()
                .exact_topology_faults()
                .iter()
                .all(|fault| fault.code == "topology.missing-incidence")
        );
        assert_eq!(tight.body().options().tolerance_mm, q("1/10000000"));
        let (source, key) = source(false, true, &mut dec);
        let refused =
            fillet::build(source, &key, "round", q("1/10"), q("1/10000000"), &mut dec).unwrap_err();
        assert_eq!(refused.code, "fillet.numerical-body");
    });
    assert!(!work.exhausted, "{work:?}");
    println!("strict tolerance control: {work:?}");
}

#[test]
fn rounding_an_outer_edge_preserves_both_untouched_periodic_seam_uses() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (source, key) = source(true, false, &mut dec);
        let seam = source
            .native_trims()
            .iter()
            .find(|a| {
                source.native_trims().iter().any(|b| {
                    a.face == b.face
                        && a.edge == b.edge
                        && (a.boundary, a.use_index) != (b.boundary, b.use_index)
                })
            })
            .unwrap();
        let face = seam.face.clone();
        let edge = seam.edge.clone();
        let sample = |body: &NumericalBody, dec: &mut Decider| {
            body.native_trims()
                .iter()
                .filter(|t| t.face == face && t.edge == edge)
                .map(|t| {
                    let range = t.curve.domain();
                    let at = Real::from((range[0] + range[1]).mul_pow2(-1));
                    let p = t.curve.point_at(&at, dec).unwrap();
                    (
                        t.boundary,
                        t.use_index,
                        p.x.as_rational().unwrap(),
                        p.y.as_rational().unwrap(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let before = sample(&source, &mut dec);
        assert_eq!(before.len(), 2);
        assert_ne!((&before[0].2, &before[0].3), (&before[1].2, &before[1].3));
        let rounded =
            fillet::build(source, &key, "round", q("1/10"), q("1/10000000"), &mut dec).unwrap();
        assert_eq!(sample(rounded.body(), &mut dec), before);
        assert!(
            rounded
                .body()
                .exact_topology_faults()
                .iter()
                .all(|fault| fault.code == "topology.missing-incidence")
        );
    });
    assert!(!work.exhausted, "{work:?}");
    println!("retained periodic seam control: {work:?}");
}

#[test]
fn cap_bridge_meets_the_adjacent_native_occurrences() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (body, key) = source(false, false, &mut dec);
        let mut definition = body.definition().clone();
        let mut native = body.native_trims().to_vec();
        let shift = q("1/100000000");
        for face in &mut definition.faces {
            let Surface::Plane(plane) = &face.surface else {
                continue;
            };
            if plane.normal.z.as_rational() == Some(q("0")) {
                continue;
            }
            for (boundary, row) in face.loops.iter_mut().enumerate() {
                for (use_index, trim) in row.iter_mut().enumerate() {
                    for uv in &mut trim.uv {
                        uv[0] = &uv[0] + &shift;
                    }
                    native.push(NumericalNativeTrim {
                        tag: catallax_topology::numerical::CopyTag::Independent,
                        face: face.key.clone(),
                        boundary,
                        use_index,
                        edge: definition.edges[trim.edge].key.clone(),
                        curve: linear_uv(&definition.edges[trim.edge], trim, &mut dec),
                    });
                }
            }
        }
        let body = body
            .rebuild_with_native_trims(definition, native, body.options().clone(), &mut dec)
            .unwrap();
        let rounded =
            fillet::build(body, &key, "round", q("1/10"), q("1/10000000"), &mut dec).unwrap();
        let body = rounded.body();
        let mut joins = 0;
        for face in body
            .definition()
            .faces
            .iter()
            .filter(|face| matches!(face.surface, Surface::Plane(_)))
        {
            for (boundary, row) in face.loops.iter().enumerate() {
                for (index, bridge) in row.iter().enumerate().filter(|(_, trim)| {
                    body.definition().edges[trim.edge]
                        .key
                        .starts_with("round/cap/")
                }) {
                    for (endpoint, neighbour, end) in [
                        (0, (index + row.len() - 1) % row.len(), 1),
                        (1, (index + 1) % row.len(), 0),
                    ] {
                        let trim = &row[neighbour];
                        let edge = &body.definition().edges[trim.edge];
                        let binding = body
                            .native_trims()
                            .iter()
                            .find(|t| {
                                t.face == face.key
                                    && t.boundary == boundary
                                    && t.use_index == neighbour
                            })
                            .unwrap();
                        let at = &edge.parameters[if trim.sense == Sense::Same {
                            end
                        } else {
                            1 - end
                        }];
                        let point = binding
                            .curve
                            .point_at(&Real::from(at.clone()), &mut dec)
                            .unwrap();
                        let residual = (&point.x.as_rational().unwrap() - &bridge.uv[endpoint][0])
                            .abs()
                            + (&point.y.as_rational().unwrap() - &bridge.uv[endpoint][1]).abs();
                        assert!(
                            residual <= q("1e-30"),
                            "native cap join residual: {residual}"
                        );
                        joins += 1;
                    }
                }
            }
        }
        assert_eq!(joins, 4);
    });
    assert!(!work.exhausted, "{work:?}");
    println!("cap native residual control: {work:?}");
}
