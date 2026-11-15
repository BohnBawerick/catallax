//! N03 numerical export and a deterministic exact legacy control.
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_shapes::Pose;
use catallax_topology::Model;
#[test]
fn n03_exact_legacy_control() {
    let mut model = Model::new();
    let body = catallax_topology::block(
        &mut model,
        [Rational::one(), Rational::one(), Rational::one()],
        &Pose::identity(),
    )
    .unwrap();
    let (answer, work) = with_operation_budget(20_000_000, || {
        let mut bytes = Vec::new();
        catallax_step::write_step(
            &model,
            &[body],
            &catallax_step::StepOptions::default(),
            &mut Decider::default(),
            &mut bytes,
        )
        .unwrap();
        bytes
    });
    assert_eq!(
        answer.as_slice(),
        include_bytes!("../../../evidence/phase1/N03/legacy-block.step")
    );
    assert_eq!(work.used, 2400);
    let mut seamed = Vec::new();
    let options = catallax_step::StepOptions {
        seams: catallax_step::Seams::Required,
        ..Default::default()
    };
    catallax_step::write_step(
        &model,
        &[body],
        &options,
        &mut Decider::default(),
        &mut seamed,
    )
    .unwrap();
    assert_eq!(seamed, answer);
    if let Ok(directory) = std::env::var("CATALLAX_N03_CAPTURE") {
        std::fs::write(
            std::path::Path::new(&directory).join("legacy-block.step"),
            &answer,
        )
        .unwrap();
    }
    println!(
        "N03 exact legacy bytes={} work={} exhausted={}",
        answer.len(),
        work.used,
        work.exhausted
    );
}

#[path = "common/n03_tube.rs"]
mod n02;
use catallax_step::numerical::{self, NumericalStepOptions};
use catallax_topology::numerical::{self as source, NumericalOptions};

fn control(gap: bool) -> catallax_topology::numerical::NumericalBody {
    let mut d = n02::tube();
    let opts = NumericalOptions {
        tolerance_mm: n02::q("0.01"),
    };
    if gap {
        let keeper = d.edges[0].vertices[0];
        let removed = d.vertices.len();
        let mut vertex = d.vertices[keeper].clone();
        vertex.key = "split-endpoint".into();
        vertex.position[0] = vertex.position[0].clone() + n02::q("0.005");
        d.vertices.push(vertex);
        d.edges[0].vertices[0] = removed;
        source::sew(&d, &[[keeper, removed]], opts, &mut Decider::default()).unwrap()
    } else {
        source::build(d, opts, &mut Decider::default()).unwrap()
    }
}
#[test]
fn n03_numerical_two_loop_and_repaired_controls_are_reread() {
    let capture = std::env::var("CATALLAX_N03_CAPTURE").ok();
    let directory = capture
        .as_ref()
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::env::temp_dir().join(format!("catallax-n03-{}", std::process::id()))
        });
    std::fs::create_dir_all(&directory).unwrap();
    for gap in [false, true] {
        let name = if gap { "repaired" } else { "two-loop" };
        let body = control(gap);
        let options = NumericalStepOptions::default();
        let (artifact, work) = with_operation_budget(20_000_000, || {
            numerical::build(&body, &options, &mut Decider::default())
        });
        let artifact = artifact.unwrap();
        assert!(!work.exhausted);
        assert_eq!(artifact.report().counts, [10, 24, 16]);
        assert_eq!(&artifact.report().source_evidence, body.evidence());
        assert!(artifact.report().emitted_evidence.repairs().is_empty());
        assert_eq!(
            artifact.report().source_evidence.repairs().len(),
            usize::from(gap)
        );
        let file = directory.join(format!("{name}.step"));
        std::fs::write(&file, artifact.bytes()).unwrap();
        let reread = std::fs::read(&file).unwrap();
        artifact
            .check_written(&reread, &mut Decider::default())
            .unwrap();
        let text = std::str::from_utf8(&reread).unwrap();
        assert!(text.contains("catallax numerical approximate"));
        assert!(text.contains("source_tolerance_mm=1/100"));
        assert!(text.contains("complete_body_bound=unestablished"));
        assert!(!text.contains("UNCERTAINTY_MEASURE_WITH_UNIT"));
        let mut corrupt = reread.clone();
        corrupt[0] = b'X';
        assert_eq!(
            artifact
                .check_written(&corrupt, &mut Decider::default())
                .unwrap_err()
                .code,
            "step.numerical-file-changed"
        );
        if capture.is_some() {
            std::fs::write(
                directory.join(format!("{name}-evidence.txt")),
                format!("{:?}\n", artifact.report()),
            )
            .unwrap();
            let report = serde_json::json!({"result":"numerical","exact":false,"complete_body_bound":null,"file":format!("{name}.step"),"source_name":artifact.report().source_name,"source_tolerance_mm":body.options().tolerance_mm.to_string(),"source_repairs":body.evidence().repairs().len(),"source_checks":body.evidence().checks().len(),"emitted_repairs":artifact.report().emitted_evidence.repairs().len(),"emitted_checks":artifact.report().emitted_evidence.checks().len(),"max_scalar_rounding_error_mm":artifact.report().max_scalar_rounding_error_mm.to_string(),"scalar_rounding_limit_mm":options.max_scalar_rounding_error_mm.to_string(),"bytes_reread_checked":true,"reference_check":"separate","counts":artifact.report().counts,"work":{"used":work.used,"limit":20_000_000,"exhausted":work.exhausted}});
            std::fs::write(
                directory.join(format!("{name}-report.json")),
                serde_json::to_vec_pretty(&report).unwrap(),
            )
            .unwrap();
        }
        println!(
            "N03 {name}: numerical bytes={} work={} reread=checked",
            reread.len(),
            work.used
        );
    }
    if capture.is_none() {
        std::fs::remove_dir_all(directory).unwrap();
    }
}
#[test]
fn n03_tolerated_source_gap_is_not_healed_on_export() {
    let mut definition = n02::tube();
    definition.vertices[0].position[0] =
        definition.vertices[0].position[0].clone() + n02::q("0.005");
    let body = source::build(
        definition,
        NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let before = format!("{:?}", body.definition());
    let refusal = numerical::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refusal.code, "step.numerical-gap");
    assert_eq!(before, format!("{:?}", body.definition()));
}
#[test]
fn n03_rounding_cannot_collapse_the_emitted_body() {
    let mut d = n02::tube();
    let scale = n02::q("0.000000001");
    let real_scale = scale.clone().into();
    for vertex in &mut d.vertices {
        for p in &mut vertex.position {
            *p = (&*p * &scale) + Rational::one();
        }
    }
    for edge in &mut d.edges {
        let catallax_shapes::Curve::Line(line) = &mut edge.curve else {
            panic!("line")
        };
        for p in [&mut line.origin.x, &mut line.origin.y, &mut line.origin.z] {
            *p = (&*p * &real_scale) + catallax_numbers::Real::one();
        }
        for p in [&mut line.dir.x, &mut line.dir.y, &mut line.dir.z] {
            *p = &*p * &real_scale;
        }
    }
    for face in &mut d.faces {
        let catallax_shapes::Surface::Plane(plane) = &mut face.surface else {
            panic!("plane")
        };
        for p in [
            &mut plane.origin.x,
            &mut plane.origin.y,
            &mut plane.origin.z,
        ] {
            *p = (&*p * &real_scale) + catallax_numbers::Real::one();
        }
        for trim in face.loops.iter_mut().flatten() {
            for p in trim.uv.iter_mut().flatten() {
                *p = &*p * &scale;
            }
        }
    }
    let body = source::build(d, NumericalOptions::default(), &mut Decider::default()).unwrap();
    let options = NumericalStepOptions {
        significant_digits: 8,
        ..Default::default()
    };
    assert_eq!(
        numerical::build(&body, &options, &mut Decider::default())
            .unwrap_err()
            .code,
        "numerical.invalid-definition"
    );
}
#[test]
fn n03_options_and_request_allowance_refuse_without_promotion() {
    let body = control(false);
    let options = NumericalStepOptions {
        significant_digits: 0,
        ..Default::default()
    };
    assert_eq!(
        numerical::build(&body, &options, &mut Decider::default())
            .unwrap_err()
            .code,
        "step.numerical-options"
    );
    let (result, work) = with_operation_budget(1, || {
        numerical::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(result.is_err());
    assert!(work.exhausted);
    let (combined, work) = with_operation_budget(20_000_000, || {
        let body = control(true);
        numerical::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
    });
    assert!(combined.is_ok());
    assert!(!work.exhausted);
    println!(
        "N03 combined fixture/source-sew/export work={} limit=20000000 exhausted={}",
        work.used, work.exhausted
    );
}

#[test]
fn n03_flat_spline_carrier_keeps_its_native_step_geometry_and_trims() {
    use catallax_numbers::Real;
    use catallax_shapes::{BSplineSurface, Point3, RationalPatch, Surface};
    let mut d = n02::tube();
    let mut dec = Decider::default();
    let controls = [(0, 0), (0, 4), (4, 0), (4, 4)]
        .map(|(x, y)| Point3::new(x.into(), y.into(), 1.into()))
        .to_vec();
    let knots = vec![n02::q("0"), n02::q("0"), n02::q("1"), n02::q("1")];
    let surface = BSplineSurface::new(
        [1, 1],
        [knots.clone(), knots],
        [2, 2],
        controls,
        vec![Real::one(); 4],
        &mut dec,
    )
    .unwrap();
    d.faces[0].surface = Surface::BSpline(RationalPatch::new(surface, &mut dec).unwrap().into());
    for trim in d.faces[0].loops.iter_mut().flatten() {
        for uv in trim.uv.iter_mut().flatten() {
            *uv = uv.clone() / n02::q("4");
        }
    }
    let body = source::build(d, NumericalOptions::default(), &mut dec).unwrap();
    let artifact = numerical::build(&body, &NumericalStepOptions::default(), &mut dec).unwrap();
    let text = std::str::from_utf8(artifact.bytes()).unwrap();
    assert!(text.contains("RATIONAL_B_SPLINE_SURFACE("));
    assert!(text.contains("=PCURVE("));
    assert_eq!(text.matches("=PLANE(").count(), body.counts().faces - 1);
    assert_eq!(text.matches("=FACE_BOUND(").count(), 2);
    assert!(artifact.report().spline_rounding.is_some());
    assert!(matches!(
        artifact.report().source_definition.faces[0].surface,
        Surface::BSpline(_)
    ));
    artifact.check_written(artifact.bytes(), &mut dec).unwrap();
}

#[test]
fn n03_product_name_round_trips_header_replacement_text() {
    for name in [
        "catallax export",
        "part GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))",
    ] {
        let options = NumericalStepOptions {
            product_name: name.into(),
            ..Default::default()
        };
        let artifact =
            numerical::build(&control(false), &options, &mut Decider::default()).unwrap();
        artifact
            .check_written(artifact.bytes(), &mut Decider::default())
            .unwrap();
        let text = std::str::from_utf8(artifact.bytes()).unwrap();
        let product = text
            .lines()
            .find_map(|line| line.strip_prefix("#4=PRODUCT("))
            .unwrap();
        assert_eq!(product.split('\'').nth(1).unwrap(), name);
        assert_eq!(product.split('\'').nth(3).unwrap(), name);
        let representation = text
            .lines()
            .find_map(|line| line.strip_prefix("#16=ADVANCED_BREP_SHAPE_REPRESENTATION("))
            .unwrap();
        assert_eq!(representation.split('\'').nth(1).unwrap(), name);
    }
}
