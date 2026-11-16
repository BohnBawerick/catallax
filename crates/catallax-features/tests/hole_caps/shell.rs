use super::*;
use catallax_boolean::numerical::NumericalCylinderSource;
use catallax_numbers::{current_operation_work, transfer};
use catallax_shapes::{BSplineCurve, Curve, Line, Surface};
use catallax_step::numerical::{self as step, NumericalStepOptions};
use catallax_topology::{
    Sense,
    numerical::{
        self as n, NumericalDefinition, NumericalEdge, NumericalFace, NumericalNativeTrim,
        NumericalTrim, NumericalVertex,
    },
};

fn at(curve: &BSplineCurve, parameter: &Rational, dec: &mut Decider) -> [Rational; 3] {
    let [low, high] = curve.domain();
    let endpoint = if parameter == low {
        curve.controls().next()
    } else if parameter == high {
        curve.controls().last()
    } else {
        None
    };
    if let Some(point) = endpoint {
        return [&point.x, &point.y, &point.z]
            .map(|value| dec.enclose_budgeted_at(value, 128).unwrap().mid());
    }
    curve
        .evaluate(&Ball::from_rational(parameter, 128), dec)
        .unwrap()
        .map(|value| value.mid())
}

fn linear_uv(a: [Rational; 3], b: [Rational; 3], dec: &mut Decider) -> BSplineCurve {
    BSplineCurve::new(
        1,
        [0, 0, 1, 1].map(Rational::from_i64).to_vec(),
        [a, b]
            .map(|p| {
                let [x, y, z] = p.map(Real::from);
                Point3::new(x, y, z)
            })
            .to_vec(),
        vec![Real::one(); 2],
        dec,
    )
    .unwrap()
}

fn face(
    definition: &mut NumericalDefinition,
    bindings: &mut Vec<NumericalNativeTrim>,
    key: String,
    surface: Surface,
    sense: Sense,
    uses: Vec<(usize, Sense, BSplineCurve)>,
    dec: &mut Decider,
) {
    let mut trims = Vec::new();
    for (use_index, (index, sense, curve)) in uses.into_iter().enumerate() {
        let edge = &definition.edges[index];
        let mut parameters = edge.parameters.clone();
        if sense == Sense::Reversed {
            parameters.reverse();
        }
        let uv = parameters.map(|t| {
            let [u, v, _] = at(&curve, &t, dec);
            [u, v]
        });
        trims.push(NumericalTrim {
            edge: index,
            sense,
            uv,
        });
        bindings.push(NumericalNativeTrim {
            face: key.clone(),
            boundary: 0,
            use_index,
            edge: edge.key.clone(),
            curve,
        });
    }
    definition.faces.push(NumericalFace {
        key,
        surface,
        sense,
        loops: vec![trims],
    });
}

fn fixture(
    diameter: &Rational,
    options: &NumericalOptions,
    dec: &mut Decider,
) -> (NumericalDefinition, Vec<NumericalNativeTrim>) {
    let authored = span_with_diameter(true, true, diameter, dec);
    let source = NumericalCylinderSource::new(
        &authored.cylinder,
        authored.interval_mm.clone(),
        options.clone(),
        dec,
    )
    .unwrap();
    assert_eq!(source.axis_parameters_mm(), &authored.interval_mm);
    let span = CylinderSpan {
        cylinder: source.original().clone(),
        interval_mm: source.axis_parameters_mm().clone(),
        bands: authored.bands,
    };
    println!(
        "shell descriptor interval=[{},{}] source={} charts={}",
        span.interval_mm[0],
        span.interval_mm[1],
        source.source_work(),
        source.chart_work()
    );
    assert_eq!(source.quarters().len(), 4);
    let sources: Vec<_> = source
        .quarters()
        .iter()
        .enumerate()
        .map(|(i, quarter)| {
            assert_eq!(quarter.id(), i);
            assert_eq!(quarter.key(), format!("cylinder/quarter-{i}"));
            assert_eq!(quarter.axis_parameters_mm(), source.axis_parameters_mm());
            let Surface::BSpline(support) = quarter.surface() else {
                panic!("expected retained source quarter")
            };
            support.carrier()
        })
        .collect();
    let caps: Vec<Vec<_>> = [true, false]
        .into_iter()
        .map(|top| {
            sources
                .iter()
                .enumerate()
                .map(|(i, source)| {
                    cap_curve(
                        &format!("cap-{top}-{i}"),
                        &span,
                        source,
                        &cap(true, top),
                        options,
                        dec,
                    )
                    .unwrap()
                    .unwrap()
                })
                .collect()
        })
        .collect();
    let mut definition = NumericalDefinition {
        name: "native-hole-cap-shell".into(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    for (side, arcs) in caps.iter().enumerate() {
        for (i, arc) in arcs.iter().enumerate() {
            definition.vertices.push(NumericalVertex {
                key: format!("cap-{side}-vertex-{i}"),
                position: at(arc.curve(), &Rational::from_i64(i as i64), dec),
            });
            definition.edges.push(NumericalEdge {
                key: format!("cap-{side}-edge-{i}"),
                curve: Curve::BSpline(arc.curve().clone()),
                vertices: [4 * side + i, 4 * side + (i + 1) % 4],
                parameters: [i as i64, i as i64 + 1].map(Rational::from_i64),
            });
        }
    }
    for i in 0..4 {
        let point = |index: usize| {
            let [x, y, z] = definition.vertices[index].position.clone().map(Real::from);
            Point3::new(x, y, z)
        };
        let origin = point(i);
        let dir = point(i + 4).minus(&origin);
        definition.edges.push(NumericalEdge {
            key: format!("seam-{i}"),
            curve: Curve::Line(Line { origin, dir }),
            vertices: [i, i + 4],
            parameters: [q("0"), q("1")],
        });
    }
    let mut bindings = Vec::new();
    for side in 0..caps.len() {
        let indices: Vec<usize> = if side == 0 {
            (0..4).rev().collect()
        } else {
            (0..4).collect()
        };
        let sense = if side == 0 {
            Sense::Reversed
        } else {
            Sense::Same
        };
        catallax_features::hole::append_planar_cap(
            &mut definition,
            &mut bindings,
            &format!("cap-{side}"),
            &cap(true, side == 0),
            if side == 0 {
                Sense::Same
            } else {
                Sense::Reversed
            },
            &[indices.into_iter().map(|i| (4 * side + i, sense)).collect()],
            dec,
        )
        .unwrap();
    }
    for (i, quarter) in source.quarters().iter().enumerate() {
        let end = Rational::from_i64(i as i64 + 1);
        let start = Rational::from_i64(i as i64);
        let seam = |t: &Rational, dec: &mut Decider| {
            linear_uv(
                at(caps[0][i].cylinder_uv(), t, dec),
                at(caps[1][i].cylinder_uv(), t, dec),
                dec,
            )
        };
        let right = seam(&end, dec);
        let left = seam(&start, dec);
        face(
            &mut definition,
            &mut bindings,
            format!("wall-{i}"),
            quarter.surface().clone(),
            Sense::Same,
            vec![
                (i, Sense::Same, caps[0][i].cylinder_uv().clone()),
                (8 + (i + 1) % 4, Sense::Same, right),
                (4 + i, Sense::Reversed, caps[1][i].cylinder_uv().clone()),
                (8 + i, Sense::Reversed, left),
            ],
            dec,
        );
    }
    (definition, bindings)
}

#[test]
fn skew_caps_sew_edit_transfer_and_export_as_a_numerical_shell() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let options = NumericalOptions {
            tolerance_mm: q("1e-7"),
        };
        let step_options = NumericalStepOptions::default();
        let mut previous: Option<n::NumericalBody> = None;
        for (label, diameter) in [("initial", "2"), ("edited", "3/2")] {
            let before = current_operation_work().used;
            let (definition, bindings) = fixture(&q(diameter), &options, &mut dec);
            let prepared = current_operation_work().used;
            let body = if let Some(previous) = &previous {
                previous.rebuild_with_native_trims(definition, bindings, options.clone(), &mut dec)
            } else {
                n::build_with_native_trims(definition, bindings, options.clone(), &mut dec)
            }
            .unwrap();
            assert_eq!(
                (
                    body.counts().faces,
                    body.counts().edges,
                    body.counts().vertices
                ),
                (6, 12, 8)
            );
            assert_eq!(body.native_trims().len(), 24);
            assert_eq!(body.history().len(), usize::from(label == "edited"));
            let admitted = current_operation_work().used;
            let live = step::build(&body, &step_options, &mut dec).unwrap();
            let exported = current_operation_work().used;
            let saved = transfer::encode(&body, &mut dec).unwrap();
            let restored: n::NumericalBody = transfer::decode(&saved, &mut dec).unwrap();
            let reopened = step::build(&restored, &step_options, &mut dec).unwrap();
            assert_eq!(live.bytes(), reopened.bytes());
            assert_eq!(restored.native_trims().len(), 24);
            assert_eq!(restored.options(), &options);
            assert_eq!(restored.history().len(), body.history().len());
            let names = |body: &n::NumericalBody| {
                body.definition()
                    .faces
                    .iter()
                    .map(|face| face.key.clone())
                    .collect::<Vec<_>>()
            };
            assert_eq!(names(&body), names(&restored));
            if let Some(previous) = &previous {
                assert_eq!(names(previous), names(&body));
            }
            if let Some(directory) = std::env::var_os("CATALLAX_CAP_PROOF") {
                let directory = std::path::PathBuf::from(directory);
                std::fs::create_dir_all(&directory).unwrap();
                let path = directory.join(format!("{label}.step"));
                std::fs::write(&path, live.bytes()).unwrap();
                live.check_written(&std::fs::read(&path).unwrap(), &mut dec)
                    .unwrap();
                std::fs::write(directory.join(format!("{label}.handover")), saved).unwrap();
            }
            println!(
                "shell {label}: prepare={} admit={} export={} cumulative={}",
                prepared - before,
                admitted - prepared,
                exported - admitted,
                current_operation_work().used
            );
            previous = Some(body);
        }
    });
    assert!(!work.exhausted);
}
