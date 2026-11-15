//! Cross-instance body adoption through the public byte API.

use catallax_app::handover::{read_body, read_numerical, write_body, write_numerical};
use catallax_features::{Built, evaluate};
use catallax_measure::Accuracy;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, RVector3};
use serde_json::{Value, json};

mod numerical_cases;
#[path = "../../catallax-step/tests/common/n03_tube.rs"]
mod numerical_tube;

fn recipe(text: &str) -> Recipe {
    let mut document: Value = serde_json::from_str(text).unwrap();
    if let Some(problem) = document.get("problem") {
        document = problem.clone();
    }
    Recipe::parse(&document.to_string()).unwrap()
}

fn block() -> Recipe {
    recipe(&json!({
        "format":"catallax-recipe/3", "name":"exact handover block", "result":"block",
        "steps":[{"id":"block","op":"block","size":["100000000000000000000000000000000000003/100000000000000000000000000000000000000","2","3"],
        "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}]
    }).to_string())
}

fn build(recipe: &Recipe) -> Built {
    evaluate(recipe, &mut Decider::default()).unwrap()
}

fn step(built: &Built) -> (String, Vec<u8>) {
    let mut bytes = Vec::new();
    let result = catallax_step::write_step(
        &built.model,
        &[built.result],
        &Default::default(),
        &mut Decider::default(),
        &mut bytes,
    );
    (format!("{result:?}"), bytes)
}

fn queries(built: &Built) -> Vec<String> {
    let model = &built.model;
    let body = built.result;
    let accuracy = Accuracy::default();
    let mut result = vec![
        format!("{:?}", model.counts(body)),
        format!("{:?}", model.check(body)),
        format!(
            "{:?}",
            catallax_measure::volume_area_centroid(model, body, &accuracy, &mut Decider::default())
        ),
        format!(
            "{:?}",
            catallax_measure::inertia(model, body, &accuracy, &mut Decider::default())
        ),
        format!(
            "{:?}",
            catallax_measure::bounding_box(model, body, &accuracy, &mut Decider::default())
        ),
    ];
    for face in model.faces_of_body(body) {
        result.push(format!(
            "{face:?}:{:?}:{:?}",
            model.face(face),
            model.face_origins(face)
        ));
        for origin in model.face_origins(face) {
            result.push(format!("{:?}", model.faces_with_origin(body, origin)));
        }
        result.push(format!(
            "{:?}",
            catallax_measure::face_area(model, face, &accuracy, &mut Decider::default())
        ));
        result.push(format!(
            "{:?}",
            catallax_measure::face_centroid(model, face, &accuracy, &mut Decider::default())
        ));
        result.push(format!(
            "{:?}",
            catallax_measure::face_bounding_box(model, face, &accuracy, &mut Decider::default())
        ));
    }
    for edge in model.edges_of_body(body) {
        let name = model.edge_name(body, edge);
        result.push(format!("{edge:?}:{:?}:{name:?}", model.edge(edge)));
        if let Some(name) = name {
            result.push(format!("{:?}", model.edges_with_name(body, &name)));
        }
        result.push(format!(
            "{:?}",
            catallax_measure::edge_length(model, edge, &accuracy, &mut Decider::default())
        ));
        result.push(format!(
            "{:?}",
            catallax_measure::edge_midpoint(model, edge, &accuracy, &mut Decider::default())
        ));
    }
    for vertex in model.vertices_of_body(body) {
        result.push(format!(
            "{vertex:?}:{:?}:{:?}",
            model.vertex(vertex),
            model.vertex_name(body, vertex)
        ));
        let point = model.geom().position(model.vertex(vertex).point);
        for coordinate in [&point.x, &point.y, &point.z] {
            result.push(format!(
                "{:?}",
                Decider::default().enclose(coordinate, &Rational::parse("1e-12").unwrap())
            ));
        }
    }
    let section = catallax_app::queries::section_curves(
        model,
        body,
        &catallax_app::queries::SectionPlane {
            origin: RPoint3::new(0.into(), 0.into(), Rational::parse("1/2").unwrap()),
            normal: RVector3::new(0.into(), 0.into(), 1.into()),
        },
        &accuracy,
        &mut Decider::default(),
    );
    match section {
        Ok(section) => {
            result.push(format!(
                "section:{}:{:?}",
                section.value.curves.len(),
                section.marks
            ));
            for curve in section.value.curves {
                result.push(format!(
                    "{:?}",
                    catallax_intersect::curve_form(
                        section.value.model.geom(),
                        curve.curve,
                        &mut Decider::default()
                    )
                ));
            }
        }
        Err(refusal) => result.push(format!("{refusal:?}")),
    }
    result
}

fn roundtrip(name: &str, recipe: &Recipe) {
    let source = build(recipe);
    let bytes = write_body(&source, &mut Decider::default()).unwrap();
    let started = std::time::Instant::now();
    let imported = read_body(&bytes, &mut Decider::default()).unwrap();
    let elapsed = started.elapsed();
    eprintln!(
        "handover {name}: bytes={} import_us={}",
        bytes.len(),
        elapsed.as_micros()
    );
    assert_eq!(
        write_body(&imported, &mut Decider::default()).unwrap(),
        bytes,
        "all retained fields: {name}"
    );
    assert_eq!(source.result, imported.result);
    assert_eq!(source.bodies, imported.bodies);
    assert_eq!(queries(&source), queries(&imported), "queries: {name}");
    assert_eq!(step(&source), step(&imported), "STEP: {name}");
    drop(source);
    assert!(
        imported.model.check(imported.result).is_empty(),
        "independent ownership: {name}"
    );
}

#[test]
fn block_retains_page_queries_and_step_bytes() {
    roundtrip("block", &block());
}

macro_rules! fixture {
    ($name:ident, $source:literal) => {
        #[test]
        fn $name() {
            roundtrip(stringify!($name), &recipe(include_str!($source)));
        }
    };
}

fixture!(p1_bracket, "../../../parts/P1-bracket.json");
fixture!(collar, "../../../parts/P2-collar.json");
fixture!(filleted_boss, "../../../parts/P3-filleted-boss.json");
fixture!(shelled_tray, "../../../parts/P4-shelled-tray.json");
fixture!(
    soliddraft_gui_bracket,
    "../../catallax-cli/tests/fixtures/soliddraft-gui-bracket.json"
);
fixture!(
    soliddraft_flush_bracket,
    "../../catallax-cli/tests/fixtures/soliddraft-flush-bracket.json"
);

#[test]
fn spline_boolean_result() {
    let mut recipe = recipe(include_str!(
        "../../catallax-features/tests/fixtures/spline-cut-cap-shell.json"
    ));
    recipe.result = catallax_recipe::StepId::new("cut").unwrap();
    roundtrip("spline-cut", &recipe);
}

#[test]
fn numerical_admission_tolerance_checks_and_history_survive() {
    let source = catallax_app::build_numerical(
        numerical_cases::fixture("0.001"),
        Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let source = source
        .rebuild(
            source.definition().clone(),
            source.options().clone(),
            &mut Decider::default(),
        )
        .unwrap();
    let bytes = write_numerical(&source, &mut Decider::default()).unwrap();
    let imported = read_numerical(&bytes, &mut Decider::default()).unwrap();
    assert_eq!(
        write_numerical(&imported, &mut Decider::default()).unwrap(),
        bytes
    );
    assert_eq!(source.evidence(), imported.evidence());
    assert_eq!(source.options(), imported.options());
    assert_eq!(source.counts(), imported.counts());
    assert_eq!(source.history().len(), imported.history().len());
    assert_eq!(
        source.exact_topology_faults(),
        imported.exact_topology_faults()
    );
    let options = catallax_step::numerical::NumericalStepOptions::default();
    let exported = |body: &catallax_app::NumericalBody| {
        catallax_step::numerical::build(body, &options, &mut Decider::default())
            .map(|file| file.bytes().to_vec())
            .map_err(|refusal| (refusal.code, refusal.detail))
    };
    assert_eq!(exported(&source), exported(&imported));
    assert!(read_body(&bytes, &mut Decider::default()).is_err());

    let planar = catallax_app::build_numerical(
        numerical_tube::tube(),
        Default::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let restored = read_numerical(
        &write_numerical(&planar, &mut Decider::default()).unwrap(),
        &mut Decider::default(),
    )
    .unwrap();
    let expected = exported(&planar).unwrap();
    assert_eq!(exported(&restored).unwrap(), expected);
}

#[test]
fn invalid_version_truncation_and_corruption_return_refusals() {
    let bytes = write_body(&build(&block()), &mut Decider::default()).unwrap();
    let mut version = bytes.clone();
    version[8] ^= 1;
    assert_eq!(
        read_body(&version, &mut Decider::default())
            .unwrap_err()
            .code,
        "handover.version-mismatch"
    );
    for length in [
        0,
        1,
        7,
        8,
        11,
        12,
        31,
        32,
        43,
        bytes.len() / 2,
        bytes.len() - 1,
    ] {
        assert_eq!(
            read_body(&bytes[..length], &mut Decider::default())
                .unwrap_err()
                .code,
            "handover.invalid-data"
        );
    }
    for index in [0, 12, bytes.len() / 2, bytes.len() - 1] {
        let mut corrupt = bytes.clone();
        corrupt[index] ^= 1;
        assert_eq!(
            read_body(&corrupt, &mut Decider::default())
                .unwrap_err()
                .code,
            "handover.invalid-data"
        );
    }
}

#[test]
fn import_spends_the_existing_request_allowance() {
    let bytes = write_body(&build(&block()), &mut Decider::default()).unwrap();
    let (answer, work) = with_operation_budget(0, || read_body(&bytes, &mut Decider::default()));
    assert!(answer.is_err());
    assert!(work.exhausted);
    let (answer, work) =
        with_operation_budget(100_000_000, || read_body(&bytes, &mut Decider::default()));
    answer.unwrap();
    assert!(!work.exhausted);
    assert!(work.used >= bytes.len() as u64);
}

#[test]
fn source_import_evidence_and_solid_keys_are_self_contained() {
    const HASH: &str = "1a602addaa4b534ef9427427d5231ba689a811a11159a5f9e7ebc17b3365c32b";
    let recipe = recipe(
        &json!({
            "format":"catallax-recipe/3", "name":"retained import", "result":"block",
            "steps":[{"id":"block","op":"import_step", "file":{"name":"source.step","sha256":HASH},
            "pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"decisions":{}}]
        })
        .to_string(),
    );
    let files = catallax_features::ImportFiles::from([(
        HASH.to_string(),
        include_bytes!("fixtures/block-seam-4um.step").as_slice(),
    )]);
    let source =
        catallax_features::evaluate_with_files(&recipe, &files, &mut Decider::default()).unwrap();
    let bytes = write_body(&source, &mut Decider::default()).unwrap();
    let imported = read_body(&bytes, &mut Decider::default()).unwrap();
    let step_id = catallax_recipe::StepId::new("block").unwrap();
    assert!(source.import_record(&step_id).is_some());
    assert_eq!(
        source.import_record(&step_id),
        imported.import_record(&step_id)
    );
    assert_eq!(
        format!("{:?}", source.import_solids(&step_id)),
        format!("{:?}", imported.import_solids(&step_id))
    );
    assert_eq!(queries(&source), queries(&imported));
    assert_eq!(step(&source), step(&imported));
}

#[test]
fn out_of_range_handles_are_a_refusal() {
    let mut source = build(&block());
    source.result = source.result.moved(1_000_000);
    let bytes = write_body(&source, &mut Decider::default()).unwrap();
    assert_eq!(
        read_body(&bytes, &mut Decider::default()).unwrap_err().code,
        "handover.invalid-data"
    );
}

#[test]
fn worker_build_and_mesh_can_be_adopted_with_warm_geometry_caches() {
    for recipe in [
        block(),
        recipe(include_str!("../../../parts/P3-filleted-boss.json")),
    ] {
        let options = catallax_mesh::MeshOptions::default();
        let source = catallax_app::build(&recipe, &options, &mut Decider::default()).unwrap();
        let bytes = write_body(&source.built, &mut Decider::default()).unwrap();
        let imported = read_body(&bytes, &mut Decider::default()).unwrap();
        assert_eq!(
            write_body(&imported, &mut Decider::default()).unwrap(),
            bytes
        );
        let display = catallax_mesh::mesh_body(
            &imported.model,
            imported.result,
            &options,
            &mut Decider::default(),
        )
        .unwrap();
        assert_eq!(format!("{:?}", source.display), format!("{display:?}"));
    }
}

#[test]
fn tiny_remaining_allowance_stops_before_a_large_root_is_materialized() {
    use catallax_numbers::transfer::{Codec, Writer};
    let mut source = build(&block());
    for _ in 0..100_000 {
        source.model.add_body(catallax_topology::Body {
            kind: catallax_topology::BodyKind::Empty,
            lumps: Vec::new(),
            marks: Default::default(),
        });
    }
    assert!(source.model.check(source.result).is_empty());
    let bytes = write_body(&source, &mut Decider::default()).unwrap();
    assert!(bytes.len() > 2_000_000);
    let ((result, written), work) = with_operation_budget(80, || {
        let mut decider = Decider::default();
        decider.charge_work(64).unwrap();
        let mut writer = Writer::new(&mut decider);
        let result = source.write(&mut writer);
        (result, writer.bytes_written())
    });
    assert!(
        result.is_err(),
        "root serialization escaped the request meter"
    );
    assert!(work.exhausted);
    assert!(
        written <= 16,
        "materialized {written} bytes with only 16 work units left"
    );
}
