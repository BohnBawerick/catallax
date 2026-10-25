//! Public exporter contracts and independent-read-back fixture generation.
use catallax_iges::{IgesError, IgesOptions, Units, write_iges};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model, block, cone, cylinder, empty, sphere, torus};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn fixture(name: &str) -> (Model, BodyId) {
    let mut m = Model::new();
    let b = match name {
        "box" => block(&mut m, [q("2"), q("3"), q("5")], &Pose::identity()).unwrap(),
        "cylinder" => cylinder(&mut m, q("2"), q("5"), &Pose::identity()).unwrap(),
        "cone" => cone(&mut m, q("3"), q("1"), q("2"), &Pose::identity()).unwrap(),
        "rounded-box" => block(&mut m, [q("1/3"), q("2/7"), q("5/11")], &Pose::identity()).unwrap(),
        "posed-cylinder" => {
            let pose = Pose::new(
                Rotation::from_quaternion(q("2"), q("1"), q("0"), q("0")).unwrap(),
                RVector3::new(q("1/3"), q("-7"), q("1000")).into(),
                &mut Decider::default(),
            )
            .unwrap();
            cylinder(&mut m, q("2/3"), q("5/7"), &pose).unwrap()
        }
        "P2-collar" => {
            let recipe =
                catallax_recipe::Recipe::parse(include_str!("../../../parts/P2-collar.json"))
                    .unwrap();
            let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
            return (built.model, built.result);
        }
        "bored-plate" => {
            let recipe=catallax_recipe::Recipe::parse(r#"{
              "format":"catallax-recipe/3","name":"bored plate",
              "steps":[{"id":"profile","op":"sketch",
              "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
              "loops":[[{"kind":"line","id":"a","to":["10","0"]},
                        {"kind":"line","id":"b","to":["10","8"]},
                        {"kind":"line","id":"c","to":["0","8"]},
                        {"kind":"line","id":"d","to":["0","0"]}],
                       [{"kind":"circle","id":"bore","centre":["5","4"],"radius":"2"}]]},
                       {"id":"plate","op":"extrude","profile":"profile","distance":"2"}],"result":"plate"}"#).unwrap();
            let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
            return (built.model, built.result);
        }
        _ => panic!("unknown fixture"),
    };
    (m, b)
}
fn export(m: &Model, b: BodyId, o: &IgesOptions) -> (Vec<u8>, catallax_iges::IgesReport) {
    let mut out = vec![];
    let report = write_iges(m, &[b], o, &mut Decider::default(), &mut out).unwrap();
    (out, report)
}
#[test]
fn public_fixtures_are_deterministic_shared_breps() {
    for (name, faces) in [
        ("box", 6),
        ("bored-plate", 7),
        ("P2-collar", 6),
        ("cylinder", 3),
        ("cone", 3),
    ] {
        let (m, b) = fixture(name);
        let (a, r) = export(&m, b, &IgesOptions::default());
        assert_eq!(a, export(&m, b, &IgesOptions::default()).0);
        assert_eq!(r.faces, faces);
        assert_eq!(r.edges - r.added_edges, m.edges_of_body(b).len());
        assert!(r.max_deviation <= IgesOptions::default().tolerance);
        assert_eq!(r.scale, q("1"));
        let s = String::from_utf8(a).unwrap();
        assert!(s.lines().all(|l| l.len() == 80));
        assert!(s.contains("8HCatallax"));
        assert_eq!(
            s.lines()
                .filter(|l| &l[72..73] == "D" && l[..8].trim() == "186")
                .count(),
            2
        );
    }
}
#[test]
fn strict_bound_gate_is_inclusive_and_never_touches_sink_on_refusal() {
    let mut m = Model::new();
    let b = block(&mut m, [q("1/3"), q("2/7"), q("5/11")], &Pose::identity()).unwrap();
    let mut options = IgesOptions {
        scale_digits: 5,
        tolerance: q("1"),
        ..Default::default()
    };
    let (_, r) = export(&m, b, &options);
    assert!(r.max_deviation > q("0"));
    assert_eq!(
        r.max_deviation,
        (&r.scale - &q("1")).abs() * r.coordinate_bound.clone()
    );
    options.tolerance = r.max_deviation.clone();
    assert_eq!(export(&m, b, &options).1, r);
    options.tolerance = &r.max_deviation / &q("2");
    let mut out = b"existing artifact".to_vec();
    let before = out.clone();
    assert!(matches!(
        write_iges(&m, &[b], &options, &mut Decider::default(), &mut out),
        Err(IgesError::ToleranceExceeded { .. })
    ));
    assert_eq!(out, before);
}
#[test]
fn unsupported_kinds_and_options_are_typed_and_atomic() {
    let mut m = Model::new();
    let bodies = [
        sphere(&mut m, q("2"), &Pose::identity()).unwrap(),
        torus(&mut m, q("3"), q("1"), &Pose::identity()).unwrap(),
    ];
    for b in bodies {
        let mut out = vec![];
        assert!(matches!(
            write_iges(
                &m,
                &[b],
                &IgesOptions::default(),
                &mut Decider::default(),
                &mut out
            ),
            Err(IgesError::UnsupportedSurface(_))
        ));
        assert!(out.is_empty());
    }
    let (m, b) = fixture("box");
    let mut out = vec![];
    let o = IgesOptions {
        scale_digits: 0,
        ..Default::default()
    };
    assert!(matches!(
        write_iges(&m, &[b], &o, &mut Decider::default(), &mut out),
        Err(IgesError::InvalidOptions)
    ));
    assert!(out.is_empty());
}
#[test]
fn rational_poses_units_and_scale_are_certified() {
    let mut m = Model::new();
    let pose = Pose::new(
        Rotation::from_quaternion(q("2"), q("1"), q("0"), q("0")).unwrap(),
        RVector3::new(q("1/3"), q("-7"), q("1000")).into(),
        &mut Decider::default(),
    )
    .unwrap();
    for b in [
        block(&mut m, [q("2"), q("3"), q("5")], &pose).unwrap(),
        cylinder(&mut m, q("2/3"), q("5/7"), &pose).unwrap(),
    ] {
        {
            let units = Units::Millimetres;
            let o = IgesOptions {
                units,
                ..Default::default()
            };
            let (bytes, r) = export(&m, b, &o);
            assert!(r.max_deviation <= o.tolerance);
            assert_eq!(r.units, units);
            assert!(String::from_utf8(bytes).unwrap().contains(match units {
                Units::Millimetres => ",2,2HMM,",
            }));
        }
    }
}
#[test]
fn empty_and_multiple_solids() {
    let mut m = Model::new();
    let z = empty(&mut m);
    let (_, r) = export(&m, z, &IgesOptions::default());
    assert_eq!(r.solids, 0);
    let a = block(&mut m, [q("1"), q("2"), q("3")], &Pose::identity()).unwrap();
    let b = cylinder(&mut m, q("1"), q("2"), &Pose::identity()).unwrap();
    let mut out = vec![];
    let r = write_iges(
        &m,
        &[a, b],
        &IgesOptions::default(),
        &mut Decider::default(),
        &mut out,
    )
    .unwrap();
    assert_eq!(r.solids, 2);
}
#[test]
#[ignore = "writes independent OpenCASCADE read-back inputs; set CATALLAX_IGES_FIXTURES"]
fn write_readback_fixtures() {
    let out = std::path::PathBuf::from(std::env::var("CATALLAX_IGES_FIXTURES").unwrap());
    std::fs::create_dir_all(&out).unwrap();
    for name in [
        "box",
        "bored-plate",
        "P2-collar",
        "cylinder",
        "cone",
        "rounded-box",
        "posed-cylinder",
    ] {
        let (m, b) = fixture(name);
        {
            let (suffix, units) = ("mm", Units::Millimetres);
            let (bytes, r) = export(
                &m,
                b,
                &IgesOptions {
                    units,
                    ..Default::default()
                },
            );
            std::fs::write(out.join(format!("{name}-{suffix}.igs")), bytes).unwrap();
            let c = m.counts(b);
            std::fs::write(out.join(format!("{name}-{suffix}.json")),format!(
                "{{\"max_deviation\":\"{}\",\"trim_deviation\":\"{}\",\"scale\":\"{}\",\"coordinate_bound\":\"{}\",\"faces\":{},\"edges\":{},\"vertices\":{},\"added_edges\":{},\"added_vertices\":{},\"source_edges\":{},\"source_vertices\":{}}}\n",
                r.max_deviation,r.trim_deviation,r.scale,r.coordinate_bound,r.faces,r.edges,r.vertices,r.added_edges,r.added_vertices,c.edges,c.vertices)).unwrap();
            println!("{name}-{suffix}: {r:?}");
        }
    }
}

#[test]
fn emitted_vertices_match_the_proved_whole_body_similarity() {
    for digits in [3, 7, 17] {
        let mut m = Model::new();
        let b = block(&mut m, [q("1/3"), q("2/7"), q("5/11")], &Pose::identity()).unwrap();
        let (bytes, r) = export(
            &m,
            b,
            &IgesOptions {
                scale_digits: digits,
                tolerance: q("1"),
                ..Default::default()
            },
        );
        let file = String::from_utf8(bytes).unwrap();
        let params: String = file
            .lines()
            .filter(|l| &l[72..73] == "P")
            .map(|l| l[..64].trim_end())
            .collect();
        let vertices = params.split(';').find(|s| s.starts_with("502,")).unwrap();
        let fields: Vec<_> = vertices.split(',').collect();
        assert_eq!(fields[1], "8");
        let emitted: std::collections::BTreeSet<Vec<_>> = fields[2..]
            .chunks(3)
            .map(|p| p.iter().map(|s| q(s)).collect())
            .collect();
        for vertex in m.vertices_of_body(b) {
            let p = m.geom().position(m.vertex(vertex).point);
            let exact = [p.x, p.y, p.z].map(|x| x.as_rational().unwrap());
            let mapped: Vec<_> = exact.iter().map(|x| x * &r.scale).collect();
            assert!(emitted.contains(&mapped));
            let distance_sq = (0..3)
                .map(|i| {
                    let d = &mapped[i] - &exact[i];
                    &d * &d
                })
                .fold(q("0"), |a, b| a + b);
            assert!(distance_sq <= &r.max_deviation * &r.max_deviation);
        }
    }
}

#[test]
fn cone_trim_error_is_composed_and_strictly_enforced() {
    let (m, b) = fixture("P2-collar");
    let (_, r) = export(&m, b, &IgesOptions::default());
    assert!(r.trim_deviation > q("0"));
    assert_eq!(r.scale, q("1"));
    assert_eq!(r.max_deviation, r.trim_deviation);
    let mut out = vec![];
    let o = IgesOptions {
        tolerance: &r.max_deviation / &q("2"),
        ..Default::default()
    };
    assert!(matches!(
        write_iges(&m, &[b], &o, &mut Decider::default(), &mut out),
        Err(IgesError::ToleranceExceeded { .. })
    ));
    assert!(out.is_empty());
}

#[test]
fn unavailable_numeric_encodings_and_work_refuse_without_output() {
    let mut m = Model::new();
    let b = cone(&mut m, q("3"), q("1"), q("3"), &Pose::identity()).unwrap();
    let mut out = vec![];
    assert!(matches!(
        write_iges(
            &m,
            &[b],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::UnrepresentableDefinition(_))
    ));
    assert!(out.is_empty());
    let (m, b) = fixture("box");
    let mut dec = Decider::new(catallax_numbers::BudgetPolicy {
        max_work: 0,
        ..Default::default()
    });
    assert!(matches!(
        write_iges(&m, &[b], &IgesOptions::default(), &mut dec, &mut out),
        Err(IgesError::Decision(_))
    ));
    assert!(out.is_empty());
}

#[test]
fn retained_step_bytes_and_source_model_are_unchanged() {
    for name in ["box", "bored-plate", "P2-collar"] {
        let (m, b) = fixture(name);
        let mut before = vec![];
        let mut after = vec![];
        let options = catallax_step::StepOptions::default();
        let report =
            catallax_step::write_step(&m, &[b], &options, &mut Decider::default(), &mut before)
                .unwrap();
        let _ = export(&m, b, &IgesOptions::default());
        let check =
            catallax_step::write_step(&m, &[b], &options, &mut Decider::default(), &mut after)
                .unwrap();
        assert_eq!(before, after);
        assert_eq!(report.value, check.value);
        assert_eq!(report.marks, check.marks);
    }
}

#[test]
fn marked_invalid_and_duplicate_bodies_are_rejected_atomically() {
    use catallax_numbers::{BoundKind, Mark};
    let (mut m, b) = fixture("box");
    let mut body = m.body(b).clone();
    body.marks.push(Mark::new(
        "test.bound",
        BoundKind::Distance,
        q("1/100"),
        "model bound",
    ));
    let marked = m.add_body(body);
    let mut invalid = m.body(b).clone();
    invalid.lumps.clear();
    let invalid = m.add_body(invalid);
    let mut out = b"retained".to_vec();
    assert!(matches!(
        write_iges(
            &m,
            &[marked],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::UncertifiedModel)
    ));
    assert!(matches!(
        write_iges(
            &m,
            &[invalid],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::InvalidTopology(_))
    ));
    assert!(matches!(
        write_iges(
            &m,
            &[b, b],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::InvalidTopology(_))
    ));
    assert_eq!(out, b"retained");
}

#[test]
fn algebraic_positions_are_not_silently_rationalised() {
    use catallax_numbers::Real;
    use catallax_shapes::Vector3;
    let pose = Pose::new(
        Rotation::identity(),
        Vector3::new(Real::from_i64(2).sqrt(), Real::zero(), Real::zero()),
        &mut Decider::default(),
    )
    .unwrap();
    let mut m = Model::new();
    let b = block(&mut m, [q("1"), q("2"), q("3")], &pose).unwrap();
    let mut out = vec![];
    assert!(matches!(
        write_iges(
            &m,
            &[b],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::UnrepresentableDefinition(_))
    ));
    assert!(out.is_empty());
}
