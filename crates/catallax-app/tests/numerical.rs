//! Public numerical admission and continuation, independent of exact checkpoints.

use catallax_app::{NumericalBody, NumericalOptions, build_numerical};
use catallax_numbers::{
    Ball, BudgetPolicy, Decider, NumericalCheck, NumericalRepair, with_operation_budget,
};
use catallax_shapes::{Curve, Point3, Surface};
use catallax_topology::{Model, Sense};

mod numerical_cases;
use numerical_cases::{fixture, from_model, q, roof};

#[path = "../../catallax-import/tests/common/mod.rs"]
mod import_cases;

fn build(delta: &str) -> NumericalBody {
    build_numerical(
        fixture(delta),
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap()
}
fn samples(body: &NumericalBody) -> Vec<&catallax_numbers::NumericalSample> {
    body.evidence()
        .checks()
        .iter()
        .filter_map(|c| match c {
            NumericalCheck::Samples {
                name,
                scope,
                samples,
            } if name == "trim.edge-surface-distance" && scope.starts_with("face \"roof\",") => {
                Some(samples.as_slice())
            }
            _ => None,
        })
        .flatten()
        .collect()
}

#[test]
fn general_trim_returns_checked_geometry_with_scoped_samples_not_complete_bounds() {
    for delta in ["0", "0.001"] {
        let body = build(delta);
        let count = body.counts();
        assert_eq!(
            (
                count.lumps,
                count.shells,
                count.faces,
                count.edges,
                count.vertices
            ),
            (1, 1, 6, 12, 8)
        );
        assert_eq!(body.evidence().tolerance_mm(), &q("0.01"));
        assert!(body.evidence().repairs().is_empty());
        let roof_samples = samples(&body);
        assert_eq!(roof_samples.len(), 36);
        // Independent formula for the retained graph. Its unit-square trim is
        // strictly inside the carrier, rather than one of its complete rows.
        for sample in &roof_samples {
            let u = &sample.parameters[1];
            let v = &sample.parameters[2];
            let expected = q("4")
                * u.clone()
                * (q("1") - u.clone())
                * v.clone()
                * (q("1") - v.clone())
                * q(delta);
            assert!(sample.distance_upper_mm >= expected);
            assert!(&sample.distance_upper_mm - &expected < q("1e-25"));
        }
        let expected = if delta == "0" { q("0") } else { q("1/4500") };
        let measured = roof_samples
            .iter()
            .map(|s| s.distance_upper_mm.clone())
            .max()
            .unwrap();
        assert!(measured >= expected && &measured - &expected < q("1e-25"));
        let faults = body.exact_topology_faults();
        assert!(!faults.is_empty());
        assert!(
            faults
                .iter()
                .all(|f| f.code == "topology.missing-incidence")
        );
        // Material witnesses from the independently known graph and footprint.
        let Surface::BSpline(support) = &body
            .definition()
            .faces
            .iter()
            .find(|f| f.key == "roof")
            .unwrap()
            .surface
        else {
            panic!("patch")
        };
        let p = support
            .carrier()
            .evaluate(
                &Ball::from_rational(&q("1/2"), 128),
                &Ball::from_rational(&q("1/2"), 128),
                &Decider::default(),
            )
            .unwrap();
        assert!(p[2].lo() > q("1/2") && p[2].hi() < q("2"));
        let inside_planes = |p: Point3| {
            body.definition().faces.iter().all(|face| {
                let Surface::Plane(plane) = &face.surface else {
                    return true;
                };
                let side = p
                    .minus(&plane.origin)
                    .dot(&plane.normal)
                    .as_rational()
                    .unwrap();
                if face.sense == Sense::Same {
                    side <= q("0")
                } else {
                    side >= q("0")
                }
            })
        };
        assert!(inside_planes(Point3::new(
            q("1/2").into(),
            q("1/2").into(),
            q("1/2").into()
        )));
        assert!(!inside_planes(Point3::new(
            q("-1/2").into(),
            q("1/2").into(),
            q("1/2").into()
        )));
    }
}

#[test]
fn surface_edit_rechecks_samples_keeps_keys_and_refuses_atomically() {
    let original = build("0.001");
    let before = original.evidence().clone();
    let changed = original
        .replace_surface("roof", roof("0.002"), &mut Decider::default())
        .unwrap();
    assert_eq!(
        changed
            .definition()
            .faces
            .iter()
            .map(|f| &f.key)
            .collect::<Vec<_>>(),
        original
            .definition()
            .faces
            .iter()
            .map(|f| &f.key)
            .collect::<Vec<_>>()
    );
    assert_eq!(
        changed.evidence().tolerance_mm(),
        original.evidence().tolerance_mm()
    );
    let maximum = samples(&changed)
        .into_iter()
        .map(|s| s.distance_upper_mm.clone())
        .max()
        .unwrap();
    assert!(maximum >= q("1/2250") && &maximum - &q("1/2250") < q("1e-25"));
    let failed = original
        .replace_surface("roof", roof("0.1"), &mut Decider::default())
        .unwrap_err();
    assert_eq!(failed.code, "numbers.numerical-tolerance-exceeded");
    assert_eq!(original.evidence(), &before);
    let rebuilt = build_numerical(
        changed.definition().clone(),
        changed.options().clone(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(rebuilt.evidence(), changed.evidence());
}

#[test]
fn actual_vertex_repair_is_recorded_and_retained_without_a_composed_bound() {
    let mut definition = fixture("0.001");
    let before = definition.vertices[0].position.clone();
    definition.vertices[0].position[0] = &before[0] + &q("0.000001");
    let key = definition.vertices[0].key.clone();
    let body = build_numerical(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(body.evidence().repairs().is_empty());
    let repaired = body
        .repair_vertex(&key, before.clone(), &mut Decider::default())
        .unwrap();
    let [
        NumericalRepair::VertexMoved {
            vertex,
            after,
            distance_upper_mm,
            ..
        },
    ] = repaired.evidence().repairs()
    else {
        panic!("one vertex repair")
    };
    assert_eq!(vertex, &key);
    assert_eq!(after, &before);
    assert!(*distance_upper_mm >= q("0.000001") && distance_upper_mm - &q("0.000001") < q("1e-25"));
    let edited = repaired
        .replace_surface("roof", roof("0.002"), &mut Decider::default())
        .unwrap();
    assert_eq!(edited.evidence().repairs(), repaired.evidence().repairs());
    let mut far = before;
    far[0] = &far[0] + &q("0.1");
    assert_eq!(
        body.repair_vertex(&key, far, &mut Decider::default())
            .unwrap_err()
            .code,
        "numerical.repair-tolerance-exceeded"
    );
    assert!(body.evidence().repairs().is_empty());
}

#[test]
fn structural_and_trim_failures_are_not_hidden_by_numerical_acceptance() {
    let mut broken = fixture("0");
    broken
        .faces
        .iter_mut()
        .find(|f| f.key == "roof")
        .unwrap()
        .loops[0]
        .pop();
    assert_eq!(
        build_numerical(broken, NumericalOptions::default(), &mut Decider::default())
            .unwrap_err()
            .code,
        "topology.loop-not-closed"
    );
    for kind in ["crossing", "reversed", "outside", "missing"] {
        let mut definition = fixture("0");
        let face = definition
            .faces
            .iter_mut()
            .find(|f| f.key == "roof")
            .unwrap();
        let mut points: Vec<_> = face.loops[0].iter().map(|t| t.uv[0].clone()).collect();
        match kind {
            "crossing" => points.swap(1, 2),
            "reversed" => points.reverse(),
            "outside" => {
                for p in &mut points {
                    p[0] = &p[0] + &q("2");
                }
            }
            _ => (),
        }
        for i in 0..4 {
            face.loops[0][i].uv = [points[i].clone(), points[(i + 1) % 4].clone()];
        }
        if kind == "missing" {
            face.loops[0][0].uv[1][0] = q("1/2");
        }
        assert_eq!(
            build_numerical(
                definition,
                NumericalOptions::default(),
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
            "numerical.invalid-trim",
            "{kind}"
        );
    }
}

#[test]
fn work_exhaustion_and_invalid_definitions_return_no_body() {
    let definition = fixture("0.001");
    let (result, work) = with_operation_budget(0, || {
        build_numerical(
            definition,
            NumericalOptions::default(),
            &mut Decider::default(),
        )
    });
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
    let mut definition = fixture("0");
    definition.edges[0].vertices[0] = usize::MAX;
    assert_eq!(
        build_numerical(
            definition,
            NumericalOptions::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "numerical.invalid-definition"
    );
    assert_eq!(
        build_numerical(
            fixture("0"),
            NumericalOptions {
                tolerance_mm: q("-1")
            },
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "numerical.invalid-definition"
    );
}

#[test]
fn concave_interior_trim_keeps_all_six_edges_on_the_curved_carrier() {
    let footprint = [
        ["0", "0"],
        ["1", "0"],
        ["1", "0.4"],
        ["0.4", "0.4"],
        ["0.4", "1"],
        ["0", "1"],
    ];
    let lower: Vec<_> = footprint.iter().map(|p| [p[0], p[1], "0"]).collect();
    let upper: Vec<_> = footprint.iter().map(|p| [p[0], p[1], "1"]).collect();
    let mut rings = vec![lower.iter().rev().copied().collect(), upper.clone()];
    for i in 0..footprint.len() {
        let j = (i + 1) % footprint.len();
        rings.push(vec![lower[i], lower[j], upper[j], upper[i]]);
    }
    let source = import_cases::polyhedra_step(&[rings]);
    let file = catallax_import::read_step(source.as_bytes()).unwrap();
    let mut model = Model::new();
    let bodies = catallax_import::build_file(&mut model, &file, None).unwrap();
    assert!(model.check(bodies[0]).is_empty());
    let definition = from_model(&model, bodies[0], "0.001");
    let body = build_numerical(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(
        (
            body.counts().faces,
            body.counts().edges,
            body.counts().vertices
        ),
        (8, 18, 12)
    );
    let roof = body
        .definition()
        .faces
        .iter()
        .find(|f| f.key == "roof")
        .unwrap();
    assert_eq!(roof.loops[0].len(), 6);
    assert_eq!(samples(&body).len(), 54);
    // The retained boundary includes the reentrant physical corner (0.4, 0.4).
    assert!(
        roof.loops[0]
            .iter()
            .any(|t| t.uv[0] == [q("7/15"), q("7/15")])
    );
}

#[test]
fn reversing_a_carrier_and_its_parameters_preserves_the_checked_body() {
    let mut definition = fixture("0.001");
    let Curve::Line(line) = &mut definition.edges[0].curve else {
        panic!("line")
    };
    line.dir = line.dir.scaled(&q("-1").into());
    for t in &mut definition.edges[0].parameters {
        *t = -t.clone();
    }
    let reversed = build_numerical(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(reversed.counts(), build("0.001").counts());
    assert!(reversed.evidence().checks().iter().any(|c| matches!(c, NumericalCheck::Samples { samples, .. } if samples.iter().any(|s| s.parameters[0] < q("0")))));
}

#[test]
fn absent_or_multiple_loops_refuse_without_discarding_a_boundary() {
    for multiple in [false, true] {
        let mut definition = fixture("0");
        let face = &mut definition.faces[0];
        if multiple {
            face.loops.push(face.loops[0].clone());
        } else {
            face.loops.clear();
        }
        assert_eq!(
            build_numerical(
                definition,
                NumericalOptions::default(),
                &mut Decider::default()
            )
            .unwrap_err()
            .code,
            "numerical.invalid-trim"
        );
    }
}

#[test]
fn lower_precision_cap_is_respected_without_reinterpreting_the_tolerance() {
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 64,
        ..BudgetPolicy::default()
    });
    assert_eq!(
        build_numerical(fixture("0"), NumericalOptions::default(), &mut dec)
            .unwrap_err()
            .code,
        "numbers.precision-limit"
    );
}

#[test]
fn unsupported_support_is_a_refusal_instead_of_unchecked_acceptance() {
    use catallax_shapes::{EllipticCylinder, Vector3};
    let mut definition = fixture("0");
    // Spheres are admitted carriers now; an elliptic cylinder is not.
    definition.faces[0].surface = Surface::EllipticCylinder(
        EllipticCylinder::new(
            Point3::new(q("0").into(), q("0").into(), q("0").into()),
            Vector3::new(q("0").into(), q("0").into(), q("1").into()),
            Vector3::new(q("1").into(), q("0").into(), q("0").into()),
            q("4").into(),
            q("1").into(),
            &mut Decider::default(),
        )
        .unwrap(),
    );
    assert_eq!(
        build_numerical(
            definition,
            NumericalOptions::default(),
            &mut Decider::default()
        )
        .unwrap_err()
        .code,
        "numerical.unsupported-carrier"
    );
}
