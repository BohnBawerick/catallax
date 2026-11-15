use catallax_app::build_incremental_with_files;
use catallax_features::ImportFiles;
use catallax_import::sha256_hex;
use catallax_measure::{Accuracy, area, centroid, inertia, mass_properties, volume_area_centroid};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Recipe, Role};
use catallax_topology::BodyKind;
use serde_json::json;

const SOURCE: &[u8] = include_bytes!("../../catallax-import/tests/fixtures/surface-square.step");

fn build(source: &[u8]) -> catallax_app::IncrementalApplicationBody {
    let hash = sha256_hex(source);
    let recipe = Recipe::parse(
        &json!({
            "format": "catallax-recipe/3", "name": "surface import",
            "steps": [{"id": "sheet", "op": "import_step",
                "file": {"name": "surface.step", "sha256": hash},
                "pose": {"rotation": ["1","0","0","0"], "translation": ["1","2","3"]}}],
            "result": "sheet"
        })
        .to_string(),
    )
    .unwrap();
    let files = ImportFiles::from([(hash, source)]);
    build_incremental_with_files(
        &recipe,
        None,
        &files,
        &Default::default(),
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn surface_import_displays_and_picks_through_the_application_entry() {
    let body = build(SOURCE);
    let built = body.checkpoint.built();
    assert_eq!(built.model.body(built.result).kind, BodyKind::Surface);
    assert!(!body.display.value.triangles.is_empty());
    assert_eq!(body.names.faces.len(), 1);
    assert_eq!(body.names.edges.len(), 4);
    assert_eq!(body.names.vertices.len(), 4);
    for triangle in &body.display.value.triangles {
        let origins = &body.names.faces[&triangle.face];
        assert!(origins.iter().any(|origin| origin.step.as_str() == "sheet"
            && matches!(&origin.role, Role::Imported { key } if key == "#30")));
    }
    for edge in body.names.edges.keys() {
        let picked = body.edge_curve(*edge, &mut Decider::default()).unwrap();
        assert_eq!(picked.edge, *edge);
    }
    assert_eq!(
        area(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut Decider::default()
        )
        .unwrap()
        .value
        .lo,
        Rational::from_i64(6)
    );
    for code in [
        centroid(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap_err()
        .code,
        inertia(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap_err()
        .code,
        mass_properties(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap_err()
        .code,
        volume_area_centroid(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut Decider::default(),
        )
        .unwrap_err()
        .code,
    ] {
        assert_eq!(code, "measure.not-solid");
    }
}

#[test]
fn a_mixed_import_collects_boundaries_as_surfaces_without_claiming_solid_material() {
    let surface_data = std::str::from_utf8(SOURCE)
        .unwrap()
        .split("DATA;\n")
        .nth(1)
        .unwrap()
        .split("ENDSEC;")
        .next()
        .unwrap()
        .replace('#', "#9009000");
    let solid = include_str!("fixtures/block-seam-4um.step");
    let source = solid.replace(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{surface_data}ENDSEC;\nEND-ISO-10303-21;"),
    );
    assert_ne!(source, solid);
    let body = build(source.as_bytes());
    let built = body.checkpoint.built();
    assert_eq!(built.model.body(built.result).kind, BodyKind::Surface);
    assert_eq!(body.names.faces.len(), 7);
    assert!(!body.display.value.triangles.is_empty());
    assert!(built.model.check(built.result).is_empty());
}
