//! Cross-instance body adoption through the public byte API.

use catallax_app::handover::{read_body, read_numerical, write_body, write_numerical};
use catallax_features::{Built, evaluate};
use catallax_measure::Accuracy;
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use catallax_shapes::{RPoint3, RVector3};
use catallax_topology::numerical::tolerance::{self, EntityKind};
use serde_json::{Value, json};

/// Whether the body's model stores `expected` for the entity, or the entity
/// is an edge re-tightened to exact, which stores nothing. A vertex is a
/// rounded point and always stores its tolerance.
fn stores_or_is_exact(
    body: &catallax_app::NumericalBody,
    kind: EntityKind,
    key: &str,
    expected: &Rational,
) -> bool {
    match body.stored_tolerance(kind, key) {
        Some(stored) => stored.mm() == expected,
        None => kind == EntityKind::Edge,
    }
}

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
    assert_eq!(source.topology_faults(), Vec::new());
    assert_eq!(imported.topology_faults(), Vec::new());
    assert_eq!(source.tolerances(), imported.tolerances());
    for entity in &source.tolerances().entities {
        assert!(stores_or_is_exact(
            &imported,
            entity.kind,
            &entity.key,
            &entity.tolerance_mm
        ));
    }
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
    // Not a version this reader takes.
    version[8] ^= 64;
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
            wires: Vec::new(),
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

fn stored_record(text: &str) -> Vec<u8> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .flat_map(str::split_ascii_whitespace)
        .map(|byte| u8::from_str_radix(byte, 16).unwrap())
        .collect()
}

#[test]
fn version_one_bodies_keep_geometry_history_and_empty_kind() {
    let bytes = stored_record(include_str!("handover_v1/certified.hex"));
    let imported = read_body(&bytes, &mut Decider::default()).unwrap();
    let source = build(&block());
    assert_eq!(queries(&source), queries(&imported));
    assert_eq!(step(&source), step(&imported));
    assert_eq!(
        &write_body(&imported, &mut Decider::default()).unwrap()[8..12],
        &3_u32.to_le_bytes()
    );

    let empty: catallax_topology::Body = catallax_numbers::transfer::decode(
        &stored_record(include_str!("handover_v1/empty.hex")),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(empty.kind, catallax_topology::BodyKind::Empty);
    assert!(empty.lumps.is_empty() && empty.wires.is_empty() && empty.marks.is_empty());

    let bytes = stored_record(include_str!("handover_v1/numerical.hex"));
    let imported = read_numerical(&bytes, &mut Decider::default()).unwrap();
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
    assert_eq!(source.counts(), imported.counts());
    assert_eq!(source.options(), imported.options());
    // A version-one record has no tolerance on any entity and no tolerance
    // report. Read back, every vertex of the body's model holds the largest
    // tolerance the record's own checks state, and so does every edge that
    // is not re-tightened to exact. The one checker passes it.
    assert!(imported.tolerances().entities.is_empty());
    assert_eq!(imported.topology_faults(), Vec::new());
    let largest =
        tolerance::largest_recorded(imported.evidence()).max(tolerance::floor_minimum_mm());
    let definition = imported.definition();
    let keys = definition
        .vertices
        .iter()
        .map(|vertex| (EntityKind::Vertex, vertex.key.as_str()))
        .chain(
            definition
                .edges
                .iter()
                .map(|edge| (EntityKind::Edge, edge.key.as_str())),
        );
    for (kind, key) in keys {
        assert!(stores_or_is_exact(&imported, kind, key, &largest), "{key}");
    }
    // The stored record holds the checks of the build that wrote it, sample
    // for sample. A rebuild measures the same definition as it is measured
    // today, with the same number of check records.
    assert_eq!(
        source.evidence().checks().len(),
        imported.evidence().checks().len()
    );
    let remeasured = imported
        .rebuild(
            imported.definition().clone(),
            imported.options().clone(),
            &mut Decider::default(),
        )
        .unwrap();
    assert_eq!(source.evidence(), remeasured.evidence());
    assert_eq!(source.history().len(), imported.history().len());
    let export = |body: &catallax_app::NumericalBody| {
        catallax_step::numerical::build(body, &Default::default(), &mut Decider::default())
            .unwrap()
            .bytes()
            .to_vec()
    };
    // The stamp of a file records its body's tolerance (check once, roadmap
    // item 3.0A step 5). The stored body holds the records of the build that
    // wrote it, from before admission measured edge gaps, so its file
    // records another number than a body measured today. Measured again it
    // writes the same bytes, and under the stamp the two files are the same.
    assert_eq!(export(&source), export(&remeasured));
    let under_the_stamp = |bytes: Vec<u8>| {
        let text = String::from_utf8(bytes).unwrap();
        text.lines()
            .filter(|line| !line.starts_with("FILE_DESCRIPTION("))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(
        under_the_stamp(export(&source)),
        under_the_stamp(export(&imported))
    );
    assert_eq!(
        &write_numerical(&imported, &mut Decider::default()).unwrap()[8..12],
        &3_u32.to_le_bytes()
    );
}

/// Version two is the format before a vertex, an edge and a face stored a
/// tolerance. The two records were written by the build of main named in
/// their first line.
#[test]
fn version_two_bodies_are_exact_and_numerical_ones_get_their_tolerances_back() {
    let bytes = stored_record(include_str!("handover_v2/certified.hex"));
    assert_eq!(&bytes[8..12], &2_u32.to_le_bytes());
    let imported = read_body(&bytes, &mut Decider::default()).unwrap();
    let source = build(&block());
    // Every entity of the record is exact, and the body answers and exports
    // as a build of today does.
    let model = &imported.model;
    assert!(
        model
            .vertices_of_body(imported.result)
            .into_iter()
            .all(|vertex| model.vertex(vertex).tolerance.is_none())
    );
    assert!(
        model
            .edges_of_body(imported.result)
            .into_iter()
            .all(|edge| model.edge(edge).tolerance.is_none())
    );
    assert!(
        model
            .faces_of_body(imported.result)
            .into_iter()
            .all(|face| model.face(face).tolerance.is_none())
    );
    assert!(model.check(imported.result).is_empty());
    assert_eq!(queries(&source), queries(&imported));
    assert_eq!(step(&source), step(&imported));
    // Written again, it is a version-three record that reads back the same.
    let rewritten = write_body(&imported, &mut Decider::default()).unwrap();
    assert_eq!(&rewritten[8..12], &3_u32.to_le_bytes());
    let again = read_body(&rewritten, &mut Decider::default()).unwrap();
    assert_eq!(queries(&again), queries(&imported));

    // A numerical body of that format stored its tolerance report beside a
    // model with no tolerance. Read back, its model holds the report's
    // values, as the model of a build of today does.
    let bytes = stored_record(include_str!("handover_v2/numerical.hex"));
    assert_eq!(&bytes[8..12], &2_u32.to_le_bytes());
    let imported = read_numerical(&bytes, &mut Decider::default()).unwrap();
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
    assert_eq!(source.tolerances(), imported.tolerances());
    assert_eq!(source.evidence(), imported.evidence());
    assert!(!imported.tolerances().entities.is_empty());
    for entity in &imported.tolerances().entities {
        assert!(
            stores_or_is_exact(&imported, entity.kind, &entity.key, &entity.tolerance_mm),
            "{}",
            entity.key
        );
    }
    assert_eq!(imported.topology_faults(), Vec::new());
    assert_eq!(source.counts(), imported.counts());
    let export = |body: &catallax_app::NumericalBody| {
        catallax_step::numerical::build(body, &Default::default(), &mut Decider::default())
            .unwrap()
            .bytes()
            .to_vec()
    };
    assert_eq!(export(&source), export(&imported));
    // Written again, it is a version-three record, and that record's model
    // holds the same tolerances without being built again.
    let rewritten = write_numerical(&imported, &mut Decider::default()).unwrap();
    assert_eq!(&rewritten[8..12], &3_u32.to_le_bytes());
    let again = read_numerical(&rewritten, &mut Decider::default()).unwrap();
    for entity in &imported.tolerances().entities {
        assert_eq!(
            again.stored_tolerance(entity.kind, &entity.key),
            imported.stored_tolerance(entity.kind, &entity.key)
        );
    }
    assert_eq!(again.topology_faults(), Vec::new());
}
