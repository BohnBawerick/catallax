use super::*;
use numerical::NumericalNativeTrim;

fn curve(controls: Vec<Point3>, end: Rational) -> BSplineCurve {
    let knots = std::iter::repeat_n(q("-1"), 4)
        .chain(std::iter::repeat_n(end, 4))
        .collect();
    BSplineCurve::new(
        3,
        knots,
        controls,
        vec![Real::one(); 4],
        &mut Decider::default(),
    )
    .unwrap()
}

fn input() -> (NumericalDefinition, Vec<NumericalNativeTrim>) {
    let mut definition = tube();
    let coefficients = [q("2"), q("-7/2"), q("9/2"), q("-1")];
    let edge = &mut definition.edges[0];
    let ends = edge
        .vertices
        .map(|i| definition.vertices[i].position.clone());
    let controls = coefficients
        .iter()
        .map(|c| {
            let p: [Rational; 3] =
                std::array::from_fn(|i| &ends[0][i] + &(c * &(&ends[1][i] - &ends[0][i])));
            Point3::new(
                p[0].clone().into(),
                p[1].clone().into(),
                p[2].clone().into(),
            )
        })
        .collect();
    edge.curve = Curve::BSpline(curve(controls, q("2")));
    edge.parameters = [q("0"), q("1")];
    let key = edge.key.clone();
    let mut bindings = Vec::new();
    for face in &definition.faces {
        for (boundary, row) in face.loops.iter().enumerate() {
            for (use_index, used) in row.iter().enumerate() {
                if used.edge != 0 {
                    continue;
                }
                let mut ends = used.uv.clone();
                if used.sense == Sense::Reversed {
                    ends.reverse();
                }
                let controls = coefficients
                    .iter()
                    .map(|c| {
                        let p: [Rational; 2] = std::array::from_fn(|i| {
                            &ends[0][i] + &(c * &(&ends[1][i] - &ends[0][i]))
                        });
                        Point3::new(
                            p[0].clone().into(),
                            p[1].clone().into(),
                            Real::from(2).sqrt(),
                        )
                    })
                    .collect();
                bindings.push(NumericalNativeTrim {
                    face: face.key.clone(),
                    boundary,
                    use_index,
                    edge: key.clone(),
                    curve: curve(controls, q("2")),
                });
            }
        }
    }
    assert_eq!(bindings.len(), 2);
    (definition, bindings)
}

fn controls() -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q("1e-7"),
    }
}

fn linear_native(
    ends: &[[Rational; 2]; 2],
    spike: bool,
    dec: &mut Decider,
) -> BSplineCurve {
    let ts = [q("0"), q("1/32"), q("1/16"), q("3/32"), q("1")];
    let direction: [Rational; 2] = std::array::from_fn(|i| &ends[1][i] - &ends[0][i]);
    let offset = q("1e-4");
    let points = ts.iter().enumerate().map(|(i, t)| {
        let mut p: [Rational; 2] = std::array::from_fn(|axis| &ends[0][axis] + &(t * &direction[axis]));
        if spike && i == 2 {
            p[0] = &p[0] - &(&direction[1] * &offset);
            p[1] = &p[1] + &(&direction[0] * &offset);
        }
        let [x, y] = p;
        Point3::new(x.into(), y.into(), Real::zero())
    }).collect();
    let knots = std::iter::once(q("0")).chain(ts).chain(std::iter::once(q("1"))).collect();
    BSplineCurve::new(1, knots, points, vec![Real::one(); 5], dec).unwrap()
}

#[test]
fn exact_exact_distance_pair_keeps_its_knot_grid_in_a_mixed_body() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let mut definition = tube();
        let numerical_face = definition.faces.iter().position(|face| {
            face.loops.len() == 1 && face.loops[0].len() == 4
                && face.loops[0].iter().all(|trim| trim.edge != 0)
        }).unwrap();
        let points: [Point3; 4] = std::array::from_fn(|i| {
            let trim = &definition.faces[numerical_face].loops[0][i];
            let edge = &definition.edges[trim.edge];
            let vertex = edge.vertices[usize::from(trim.sense == Sense::Reversed)];
            let [x, y, z] = definition.vertices[vertex].position.clone();
            Point3::new(x.into(), y.into(), z.into())
        });
        let [p0, p1, p2, p3] = points;
        let mut dec = Decider::default();
        let carrier = BSplineSurface::new(
            [1, 1], std::array::from_fn(|_| vec![q("0"), q("0"), q("1"), q("1")]), [2, 2],
            vec![p0, p3, p1, p2], vec![Real::one(); 4], &mut dec,
        ).unwrap();
        let atlases = catallax_shapes::RationalAtlas::partition(carrier, [false, false], &mut dec).unwrap();
        let [atlas]: [_; 1] = atlases.try_into().unwrap();
        definition.faces[numerical_face].surface = Surface::BSpline(SplineSupport::RationalAtlas(atlas));
        let corners = [[q("0"), q("0")], [q("1"), q("0")], [q("1"), q("1")], [q("0"), q("1")]];
        for (i, trim) in definition.faces[numerical_face].loops[0].iter_mut().enumerate() {
            trim.uv = [corners[i].clone(), corners[(i + 1) % 4].clone()];
        }
        let key = definition.faces[numerical_face].key.clone();
        let mut bindings = Vec::new();
        let mut target = None;
        for (face_index, face) in definition.faces.iter().enumerate() {
            for (boundary, row) in face.loops.iter().enumerate() {
                for (use_index, trim) in row.iter().enumerate() {
                    if trim.edge != 0 && face_index != numerical_face { continue; }
                    let mut ends = trim.uv.clone();
                    if trim.sense == Sense::Reversed { ends.reverse(); }
                    if face_index == 0 && trim.edge == 0 {
                        target = Some((bindings.len(), ends.clone()));
                    }
                    bindings.push(NumericalNativeTrim {
                        face: face.key.clone(), boundary, use_index,
                        edge: definition.edges[trim.edge].key.clone(),
                        curve: linear_native(&ends, false, &mut dec),
                    });
                }
            }
        }
        let body = numerical::build_with_sampled_carriers(
            definition, bindings, std::collections::BTreeSet::from([key]), controls(), &mut dec,
        ).unwrap();
        let encoded = catallax_numbers::transfer::encode(&body, &mut dec).unwrap();
        let body: numerical::NumericalBody = catallax_numbers::transfer::decode(&encoded, &mut dec).unwrap();
        assert_eq!(catallax_numbers::transfer::encode(&body, &mut dec).unwrap(), encoded);
        let (slot, ends) = target.unwrap();
        let mut changed = body.native_trims().to_vec();
        changed[slot].curve = linear_native(&ends, true, &mut dec);
        let refusal = body.rebuild_with_native_trims(
            body.definition().clone(), changed, controls(), &mut dec,
        ).unwrap_err();
        assert_eq!(refusal.code, "numbers.numerical-tolerance-exceeded");
    });
    assert!(!work.exhausted, "{work:?}");
}

#[test]
fn native_curves_survive_rebuild_vertex_repair_and_support_edit() {
    let ((), work) = with_operation_budget(20_000_000, || {
        let (definition, bindings) = input();
        let original = definition.vertices[0].position.clone();
        let mut dec = Decider::default();
        let body =
            numerical::build_with_native_trims(definition, bindings.clone(), controls(), &mut dec)
                .unwrap();
        let rebuilt = body
            .rebuild(body.definition().clone(), controls(), &mut dec)
            .unwrap();
        let mut moved = original.clone();
        moved[2] = &moved[2] + &q("1e-8");
        let repaired = rebuilt
            .repair_vertex("v0", moved.clone(), &mut dec)
            .unwrap();
        assert_eq!(repaired.evidence().repairs().len(), 1);
        assert_eq!(repaired.definition().vertices[0].position, moved);
        assert_eq!(repaired.source_definition().vertices[0].position, original);
        let face = &repaired.definition().faces[0];
        let Surface::Plane(mut plane) = face.surface.clone() else {
            panic!("plane");
        };
        plane.origin.z = &plane.origin.z + &Real::from(q("1e-8"));
        let edited = repaired
            .replace_surface(&face.key, Surface::Plane(plane), &mut dec)
            .unwrap();
        assert_eq!(edited.history().len(), 3);
        assert_eq!(edited.counts().edges, 24);
        for stage in edited.history() {
            assert_eq!(stage.native_trims.len(), bindings.len());
            assert!(
                stage.native_trims[0]
                    .curve
                    .same_definition(&bindings[0].curve)
            );
        }
        for retained in [edited.native_trims(), edited.source_native_trims()] {
            assert_eq!(retained.len(), bindings.len());
            for (a, b) in retained.iter().zip(&bindings) {
                assert!(a.curve.same_definition(&b.curve));
                assert_eq!(a.curve.domain(), [&q("-1"), &q("2")]);
            }
        }
        assert!(edited.evidence().checks().iter().any(|c|matches!(c,NumericalCheck::Samples{name,..} if name=="trim.native-endpoint-position")));
        assert!(edited.evidence().checks().iter().any(|c|matches!(c,NumericalCheck::Passed{name,..} if name=="trim.sampled-simple-oriented-parameter-loop")));
        let bytes = catallax_numbers::transfer::encode(&edited, &mut dec).unwrap();
        let restored: numerical::NumericalBody =
            catallax_numbers::transfer::decode(&bytes, &mut dec).unwrap();
        assert_eq!(restored.native_trims().len(), bindings.len());
        assert_eq!(restored.source_native_trims().len(), bindings.len());
        assert_eq!(restored.history().len(), 3);
        for stage in restored.history() {
            assert_eq!(stage.native_trims.len(), bindings.len());
        }
        assert_eq!(
            catallax_numbers::transfer::encode(&restored, &mut dec).unwrap(),
            bytes
        );
    });
    assert!(!work.exhausted, "{work:?}");
}

#[test]
fn duplicate_stale_reversed_and_partial_native_bindings_refuse() {
    let (definition, bindings) = input();
    for control in [
        "duplicate",
        "face",
        "address",
        "edge",
        "partial",
        "reverse",
        "off-curve",
    ] {
        let mut changed = bindings.clone();
        match control {
            "duplicate" => changed.push(bindings[0].clone()),
            "face" => changed[0].face = "missing".into(),
            "address" => changed[0].use_index = usize::MAX,
            "edge" => changed[0].edge = "missing".into(),
            "partial" => {
                changed[0].curve = curve(bindings[0].curve.controls().cloned().collect(), q("1/2"))
            }
            "reverse" => {
                let mut points = bindings[0].curve.controls().cloned().collect::<Vec<_>>();
                points.reverse();
                changed[0].curve = curve(points, q("2"));
            }
            "off-curve" => {
                let mut points = bindings[0].curve.controls().cloned().collect::<Vec<_>>();
                points[1].y = &points[1].y + &Real::from(q("1/100"));
                changed[0].curve = curve(points, q("2"));
            }
            _ => unreachable!(),
        }
        let (result, work) = with_operation_budget(20_000_000, || {
            numerical::build_with_native_trims(
                definition.clone(),
                changed,
                controls(),
                &mut Decider::default(),
            )
        });
        assert!(!work.exhausted, "{control}: {work:?}");
        let error = result.unwrap_err();
        assert!(
            matches!(
                error.code,
                "numerical.invalid-trim" | "numbers.numerical-tolerance-exceeded"
            ),
            "{control}: {error:?}"
        );
    }
}

#[test]
fn legacy_positional_native_body_and_history_payloads_refuse() {
    let mut dec = Decider::default();
    let body = catallax_numbers::transfer::decode::<numerical::NumericalBody>(
        include_bytes!("../fixtures/numerical-native-body-v1.handover"),
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(body.code, "handover.version-mismatch");
    let history = catallax_numbers::transfer::decode::<numerical::NumericalHistory>(
        include_bytes!("../fixtures/numerical-native-history-v1.handover"),
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(history.code, "handover.version-mismatch");
}
