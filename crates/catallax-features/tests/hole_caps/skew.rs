use super::*;
use catallax_boolean::numerical::{
    NumericalCylinderSource, NumericalFacePartition, clip_cone_cylinder_face,
};
use catallax_features::{Built, evaluate, hole::append_planar_cap};
use catallax_numbers::transfer;
use catallax_recipe::Recipe;
use catallax_shapes::{BSplineCurve, Curve, Surface};
use catallax_topology::{
    Loop, Sense,
    numerical::{NumericalDefinition, NumericalEdge, NumericalNativeTrim, NumericalVertex},
};
use std::collections::{BTreeMap, BTreeSet};

pub(super) fn source() -> Recipe {
    let side = |segment: &str| json!({"step":"ring","role":{"side":{"segment":segment}}});
    Recipe::parse(&json!({"format":"catallax-recipe/3","name":"E07 collar skew prefix","steps":[
        {"id":"section","op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"loops":[[
            {"kind":"line","id":"bottom","to":["12","0"]},
            {"kind":"line","id":"right","to":["12","12"]},
            {"kind":"line","id":"top","to":["6","12"]},
            {"kind":"line","id":"left","to":["6","0"]}]]},
        {"id":"ring","op":"revolve","profile":"section","axis":{"point":["0","0","0"],"direction":["0","1","0"]},"angle":"360"},
        {"id":"collar","op":"chamfer","body":"ring","distance":"1","edges":[
            {"between":[side("right"),side("bottom")]},{"between":[side("right"),side("top")]}]}
    ],"result":"collar"}).to_string()).unwrap()
}

fn at(curve: &BSplineCurve, t: &Rational, dec: &mut Decider) -> [Rational; 3] {
    curve
        .evaluate(&Ball::from_rational(t, 128), dec)
        .unwrap()
        .map(|v| v.mid())
}

fn cap_vertex(
    definition: &mut NumericalDefinition,
    known: &mut BTreeMap<Rational, usize>,
    curve: &BSplineCurve,
    parameter: &Rational,
    dec: &mut Decider,
) -> usize {
    let key = if parameter == &q("4") {
        q("0")
    } else {
        parameter.clone()
    };
    if let Some(index) = known.get(&key) {
        return *index;
    }
    let index = definition.vertices.len();
    definition.vertices.push(NumericalVertex {
        key: format!("top/tool-seam-{key}"),
        position: at(curve, parameter, dec),
    });
    known.insert(key, index);
    index
}

fn polygon(
    definition: &NumericalDefinition,
    native: &[NumericalNativeTrim],
    face_index: usize,
    boundary: usize,
    dec: &mut Decider,
) -> Vec<[f64; 2]> {
    let face = &definition.faces[face_index];
    let mut points = Vec::new();
    for (use_index, trim) in face.loops[boundary].iter().enumerate() {
        let Some(binding) = native
            .iter()
            .find(|b| b.face == face.key && b.boundary == boundary && b.use_index == use_index)
        else {
            points.push(trim.uv[0].each_ref().map(Rational::to_f64_nearest));
            continue;
        };
        let edge = &definition.edges[trim.edge];
        assert_eq!(edge.key, binding.edge);
        let mut parameters = edge.parameters.clone();
        if trim.sense == Sense::Reversed {
            parameters.reverse();
        }
        for i in 0..64 {
            let t = &parameters[0]
                + &((&parameters[1] - &parameters[0]) * Rational::from_ratio(i, 64).unwrap());
            let [u, v, _] = at(&binding.curve, &t, dec);
            points.push([u.to_f64_nearest(), v.to_f64_nearest()]);
        }
    }
    points
}

fn inside(point: [f64; 2], polygon: &[[f64; 2]]) -> bool {
    let mut parity = false;
    for (a, b) in polygon
        .iter()
        .zip(polygon.iter().cycle().skip(1))
        .take(polygon.len())
    {
        if (a[1] > point[1]) != (b[1] > point[1])
            && point[0] < a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        {
            parity = !parity;
        }
    }
    parity
}

fn cap_material(
    definition: &NumericalDefinition,
    bindings: &[NumericalNativeTrim],
    dec: &mut Decider,
) {
    let face = definition.faces.last().unwrap();
    assert_eq!(face.key, "top");
    assert_eq!(face.loops.len(), 2);
    let Surface::Plane(plane) = &face.surface else {
        panic!("original plane")
    };
    let x = plane
        .ref_dir
        .scaled(&(&Real::one() / &plane.ref_dir.norm_squared().sqrt()));
    let normal = plane
        .normal
        .scaled(&(&Real::one() / &plane.normal.norm_squared().sqrt()));
    let y = normal.cross(&x);
    let outer = polygon(definition, bindings, definition.faces.len() - 1, 0, dec);
    let bore = polygon(definition, bindings, definition.faces.len() - 1, 1, dec);
    let origin = [&plane.origin.x, &plane.origin.y, &plane.origin.z].map(Real::to_f64);
    let axes = [&x, &y].map(|v| [&v.x, &v.y, &v.z].map(Real::to_f64));
    for binding in bindings.iter().filter(|binding| binding.face == face.key) {
        let edge = &definition.edges[face.loops[binding.boundary][binding.use_index].edge];
        let Curve::BSpline(curve) = &edge.curve else {
            panic!("retained cap spline")
        };
        for i in 0..=8 {
            let t = &edge.parameters[0]
                + &((&edge.parameters[1] - &edge.parameters[0])
                    * Rational::from_ratio(i, 8).unwrap());
            let xyz = at(curve, &t, dec).map(|value| value.to_f64_nearest());
            let uv = at(&binding.curve, &t, dec).map(|value| value.to_f64_nearest());
            let residual = (0..3)
                .map(|j| (xyz[j] - origin[j] - uv[0] * axes[0][j] - uv[1] * axes[1][j]).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(residual <= 1e-7, "{} at {t}: {residual}", edge.key);
        }
    }
    for (point, expected) in [
        (["8", "12", "0"], true),
        (["0", "12", "0"], false),
        (["21/2", "12", "0"], false),
        (["12", "12", "0"], false),
    ] {
        let [a, b, c] = point.map(|s| Real::from(q(s)));
        let relative = Point3::new(a, b, c).minus(&plane.origin);
        let uv = [relative.dot(&x).to_f64(), relative.dot(&y).to_f64()];
        assert_eq!(
            inside(uv, &outer) && !inside(uv, &bore),
            expected,
            "{point:?}"
        );
    }
}

fn sliver_material(
    definition: &NumericalDefinition,
    bindings: &[NumericalNativeTrim],
    dec: &mut Decider,
) -> Vec<String> {
    definition
        .faces
        .iter()
        .enumerate()
        .filter(|(_, face)| face.key.starts_with("upper/"))
        .filter_map(|(index, face)| {
            let outer = polygon(definition, bindings, index, 0, dec);
            let holes: Vec<_> = (1..face.loops.len())
                .map(|boundary| polygon(definition, bindings, index, boundary, dec))
                .collect();
            if !inside([0.001, 0.99], &outer)
                || holes.iter().any(|hole| inside([0.001, 0.99], hole))
            {
                return None;
            }
            let Surface::BSpline(support) = &face.surface else {
                panic!("retained cone source")
            };
            let xyz = support
                .carrier()
                .evaluate(
                    &Ball::from_rational(&q("1/1000"), 128),
                    &Ball::from_rational(&q("99/100"), 128),
                    dec,
                )
                .unwrap();
            let h = q("7/9997");
            let denominator = q("1") + &h * &h;
            let expected = [
                q("1199/100") * (q("1") - &h * &h) / denominator.clone(),
                q("1101/100"),
                q("1199/100") * q("2") * h / denominator,
            ];
            for (actual, expected) in xyz.iter().zip(expected) {
                assert!((actual.mid() - expected).abs() <= q("1e-7"));
            }
            Some(face.key.clone())
        })
        .collect()
}

pub(super) struct Component {
    pub definition: NumericalDefinition,
    pub bindings: Vec<NumericalNativeTrim>,
    pub cone: NumericalFacePartition,
    pub built: Built,
    pub span: CylinderSpan,
}

pub(super) fn component(diameter: &Rational, dec: &mut Decider) -> Component {
    let built = evaluate(&source(), dec).unwrap();
    let faces = built.model.faces_of_body(built.result);
    let upper = *faces.iter().find(|id| matches!(built.model.geom().surface(built.model.face(**id).surface),Surface::Cone(c) if c.apex.y.to_f64()>12.)).unwrap();
    let top = *faces.iter().find(|id| matches!(built.model.geom().surface(built.model.face(**id).surface),Surface::Plane(p) if p.origin.y.as_rational()==Some(q("12")))).unwrap();
    let Surface::Plane(plane) = built.model.geom().surface(built.model.face(top).surface) else {
        unreachable!()
    };
    let placement: HolePlacement = serde_json::from_value(
        json!({"id":"radial-0","point":["21/2","12","0"],"direction":["1","-1","1"]}),
    )
    .unwrap();
    let span = HoleTool::new(
        &placement,
        diameter,
        &HoleDepth::Through,
        &HoleTip::Flat,
        None,
        &q("26"),
        dec,
    )
    .unwrap()
    .cylinder_span()
    .unwrap();
    assert_eq!(span.interval_mm, [q("-99/2"), q("99/2")]);
    let options = NumericalOptions {
        tolerance_mm: q("1e-7"),
        vertex_tolerance_mm: None,
    };
    let part = clip_cone_cylinder_face(
        &built.model,
        upper,
        "upper",
        &span.cylinder,
        span.interval_mm.clone(),
        options.clone(),
        dec,
    )
    .unwrap()
    .unwrap();
    let original_cap_edges: BTreeSet<_> = built
        .model
        .face(top)
        .loops
        .iter()
        .flat_map(|lp| {
            let Loop::Coedges(uses) = built.model.loop_(*lp) else {
                panic!("circular cap")
            };
            uses.iter().map(|id| built.model.coedge(*id).edge)
        })
        .collect();
    let ring = part
        .rings
        .iter()
        .position(|ring| original_cap_edges.contains(&ring.edge))
        .unwrap();
    let mut outer = Vec::new();
    let mut ring_vertices = BTreeSet::new();
    for link in part.ring_links.iter().filter(|link| link.ring == ring) {
        let uses: Vec<_> = part
            .definition
            .faces
            .iter()
            .flat_map(|f| f.loops.iter().flatten())
            .filter(|t| t.edge == link.edge)
            .collect();
        assert_eq!(uses.len(), 1);
        outer.push((
            link.edge,
            if uses[0].sense == Sense::Same {
                Sense::Reversed
            } else {
                Sense::Same
            },
        ));
        ring_vertices.extend(part.definition.edges[link.edge].vertices);
    }
    let mut known = BTreeMap::new();
    for event in &part.events {
        if let Some(vertex) = event.vertex.filter(|v| ring_vertices.contains(v)) {
            for hit in &event.hits {
                known.insert(hit.uv[1][0].clone(), vertex);
            }
        }
    }
    assert_eq!(known.len(), 2);
    let mut definition = part.definition.clone();
    let mut bindings = part.native_trims.clone();
    let retained_keys: Vec<_> = definition.faces.iter().map(|f| f.key.clone()).collect();
    for quarter in part.cylinder_source.quarters() {
        let Surface::BSpline(support) = quarter.surface() else {
            panic!("native quarter")
        };
        let cap = cap_curve(
            &format!("top/tool-quarter-{}", quarter.id()),
            &span,
            support.carrier(),
            plane,
            &options,
            dec,
        )
        .unwrap()
        .unwrap();
        let [a, b] = cap.curve().domain();
        let mut nodes: BTreeMap<_, _> = known
            .iter()
            .filter(|(u, _)| *u >= a && *u <= b)
            .map(|(u, v)| (u.clone(), Some(*v)))
            .collect();
        nodes.entry(a.clone()).or_insert(None);
        nodes.entry(b.clone()).or_insert(None);
        let nodes: Vec<_> = nodes.into_iter().collect();
        for (segment, pair) in nodes.windows(2).enumerate() {
            let mid = (&pair[0].0 + &pair[1].0) / q("2");
            let p = at(cap.curve(), &mid, dec);
            let radial = &p[0] * &p[0] + &p[2] * &p[2];
            if radial <= q("36") || radial >= q("121") {
                continue;
            }
            let vertices = std::array::from_fn(|i| {
                pair[i].1.unwrap_or_else(|| {
                    cap_vertex(&mut definition, &mut known, cap.curve(), &pair[i].0, dec)
                })
            });
            let edge = definition.edges.len();
            definition.edges.push(NumericalEdge {
                key: format!("top/tool-quarter-{}/segment-{segment}", quarter.id()),
                curve: Curve::BSpline(cap.curve().clone()),
                vertices,
                parameters: [pair[0].0.clone(), pair[1].0.clone()],
            });
            outer.push((edge, Sense::Same));
        }
    }
    let bore = faces
        .iter()
        .find_map(
            |id| match built.model.geom().surface(built.model.face(*id).surface) {
                Surface::Cylinder(c) if c.radius_sq.as_rational() == Some(q("36")) => Some(c),
                _ => None,
            },
        )
        .unwrap();
    let axis = bore
        .axis
        .scaled(&(&Real::one() / &bore.axis.norm_squared().sqrt()));
    let height = |y: i64| {
        Point3::new(Real::zero(), Real::from(y), Real::zero())
            .minus(&bore.origin)
            .dot(&axis)
    };
    let mut heights = [
        dec.enclose_budgeted_at(&height(0), 128).unwrap().mid(),
        dec.enclose_budgeted_at(&height(12), 128).unwrap().mid(),
    ];
    heights.sort();
    let row = usize::from(dec.enclose_budgeted_at(&height(12), 128).unwrap().mid() == heights[1]);
    let bore = NumericalCylinderSource::new(bore, heights, options, dec).unwrap();
    let source = bore.source();
    let circle = BSplineCurve::new(
        source.degrees()[0],
        source.knots()[0].to_vec(),
        (0..9)
            .map(|i| source.controls()[2 * i + row].clone())
            .collect(),
        (0..9)
            .map(|i| source.weights()[2 * i + row].clone())
            .collect(),
        dec,
    )
    .unwrap();
    let first = definition.vertices.len();
    for i in 0..4 {
        definition.vertices.push(NumericalVertex {
            key: format!("top/bore-vertex-{i}"),
            position: at(&circle, &Rational::from_i64(i), dec),
        });
    }
    let material_normal =
        plane
            .normal
            .scaled(&Real::from(if built.model.face(top).sense == Sense::Same {
                1
            } else {
                -1
            }));
    let sense = if axis.dot(&material_normal).to_f64() > 0. {
        Sense::Reversed
    } else {
        Sense::Same
    };
    let mut inner = Vec::new();
    for i in 0..4 {
        let index = definition.edges.len();
        definition.edges.push(NumericalEdge {
            key: format!("top/bore-edge-{i}"),
            curve: Curve::BSpline(circle.clone()),
            vertices: [first + i, first + (i + 1) % 4],
            parameters: [
                Rational::from_i64(i as i64),
                Rational::from_i64(i as i64 + 1),
            ],
        });
        inner.push((index, sense));
    }
    append_planar_cap(
        &mut definition,
        &mut bindings,
        "top",
        plane,
        built.model.face(top).sense,
        &[outer, inner],
        dec,
    )
    .unwrap();
    assert_eq!(
        definition.faces[..retained_keys.len()]
            .iter()
            .map(|f| f.key.clone())
            .collect::<Vec<_>>(),
        retained_keys
    );
    assert!(retained_keys.contains(&"upper/quarter-0/region-0".to_owned()));
    Component {
        definition,
        bindings,
        cone: part,
        built,
        span,
    }
}

fn fixture(
    diameter: &Rational,
    dec: &mut Decider,
) -> (NumericalDefinition, Vec<NumericalNativeTrim>) {
    let Component {
        definition,
        bindings,
        ..
    } = component(diameter, dec);
    (definition, bindings)
}

#[test]
fn actual_skew_cap_keeps_the_bore_and_both_material_sides_after_edit_and_transfer() {
    let mut previous_keys = None;
    for diameter in [q("2"), q("3/2")] {
        let ((definition, bindings), work) = with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let result = fixture(&diameter, &mut dec);
            cap_material(&result.0, &result.1, &mut dec);
            let sliver = sliver_material(&result.0, &result.1, &mut dec);
            assert_eq!(sliver, ["upper/quarter-0/region-0"]);
            let mut missing = result.0.clone();
            missing.faces.retain(|face| face.key != sliver[0]);
            assert!(sliver_material(&missing, &result.1, &mut dec).is_empty());
            let bytes = transfer::encode(&result, &mut dec).unwrap();
            let mut refused = result.clone();
            let face = refused.0.faces.last().unwrap();
            let Surface::Plane(plane) = &face.surface else {
                panic!("original cap plane")
            };
            let plane = plane.clone();
            let sense = face.sense;
            let disconnected = face
                .loops
                .iter()
                .flatten()
                .map(|trim| (trim.edge, trim.sense))
                .collect();
            let error = append_planar_cap(
                &mut refused.0,
                &mut refused.1,
                "invalid-disconnected-cap",
                &plane,
                sense,
                &[disconnected],
                &mut dec,
            )
            .unwrap_err();
            assert_eq!(error.code, "hole.cap-boundary");
            assert_eq!(transfer::encode(&refused, &mut dec).unwrap(), bytes);
            let restored: (NumericalDefinition, Vec<NumericalNativeTrim>) =
                transfer::decode(&bytes, &mut dec).unwrap();
            cap_material(&restored.0, &restored.1, &mut dec);
            assert_eq!(sliver_material(&restored.0, &restored.1, &mut dec), sliver);
            result
        });
        assert!(!work.exhausted);
        assert_eq!(definition.faces.last().unwrap().loops[0].len(), 8);
        assert_eq!(definition.faces.last().unwrap().loops[1].len(), 4);
        assert_eq!(bindings.iter().filter(|b| b.face == "top").count(), 12);
        let keys = (
            definition
                .vertices
                .iter()
                .map(|v| v.key.clone())
                .collect::<Vec<_>>(),
            definition
                .edges
                .iter()
                .map(|e| e.key.clone())
                .collect::<Vec<_>>(),
            definition
                .faces
                .iter()
                .map(|f| f.key.clone())
                .collect::<Vec<_>>(),
        );
        if let Some(previous) = previous_keys {
            assert_eq!(keys, previous);
        }
        previous_keys = Some(keys);
        println!("actual skew cap diameter={diameter} work={}", work.used);
    }
}
