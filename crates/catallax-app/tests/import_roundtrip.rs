//! Self-authored STEP bytes through the same import and recipe entries as SolidDraft.
use catallax_app::build_incremental_with_files;
use catallax_features::ImportFiles;
use catallax_import::{ImportOptions, apply, propose, read_step, rebuild, sha256_hex};
use catallax_measure::{Accuracy, bounding_box, volume};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_recipe::{Recipe, Role, StepId};
use catallax_shapes::Pose;
use catallax_step::{StepOptions, write_step};
use catallax_topology::{BodyId, Model, block, cylinder};
use serde_json::json;

fn recipe(hash: &str) -> Recipe {
    Recipe::parse(
        &json!({
            "format": "catallax-recipe/3",
            "name": "analytic roundtrip",
            "steps": [{
                "id": "part", "op": "import_step",
                "file": {"name": "self.step", "sha256": hash},
                "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]},
                "decisions": {}
            }],
            "result": "part"
        })
        .to_string(),
    )
    .unwrap()
}

fn roundtrip(original: &Model, body: BodyId, faces: usize, diagonal: f64) {
    let mut source = Vec::new();
    let export = write_step(
        original,
        &[body],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut source,
    )
    .unwrap()
    .value;
    assert!(export.max_deviation <= StepOptions::default().curve_tolerance);
    let hash = sha256_hex(&source);
    let options = ImportOptions::default();
    let proposals = propose(&source, &options).unwrap();
    assert!(
        proposals.proposals.is_empty(),
        "a clean export needs no repair"
    );
    let mut model = Model::new();
    let step = StepId::new("part").unwrap();
    let imported = apply(&mut model, &source, "self.step", &options, &[], Some(&step)).unwrap();
    assert!(imported.left_open.is_empty());
    assert_eq!(imported.bodies.len(), 1);
    assert!(imported.record.changes.is_empty());
    assert!(
        imported.derivative.is_none(),
        "an exact import must not move geometry"
    );
    let record = serde_json::from_str(&serde_json::to_string(&imported.record).unwrap()).unwrap();
    let reopened = rebuild(&mut Model::new(), &source, &options, &record, Some(&step)).unwrap();
    assert_eq!(reopened.record, imported.record);

    let definition = recipe(&hash);
    let files = ImportFiles::from([(hash.clone(), source.as_slice())]);
    let built = build_incremental_with_files(
        &definition,
        None,
        &files,
        &MeshOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let exact = built.checkpoint.built();
    assert!(exact.model.check(exact.result).is_empty());
    assert_eq!(built.names.faces.len(), faces);
    let file = read_step(&source).unwrap();
    for origins in built.names.faces.values() {
        let [origin] = origins.as_slice() else {
            panic!("one file face origin")
        };
        assert_eq!(origin.step, step);
        let Role::Imported { key } = &origin.role else {
            panic!("imported origin")
        };
        assert!(file.solids[0].faces.iter().any(|face| &face.name == key));
    }
    let expected = volume(
        original,
        body,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    let measured = volume(
        &exact.model,
        exact.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    assert!(expected.lo <= measured.hi && measured.lo <= expected.hi);
    assert!(&measured.hi - &measured.lo < Rational::parse("0.000001").unwrap());
    let bounds = |model: &Model, body| {
        bounding_box(model, body, &Accuracy::default(), &mut Decider::default())
            .unwrap()
            .value
            .unwrap()
    };
    assert_eq!(bounds(original, body), bounds(&exact.model, exact.result));

    // Prove carrier and boundary identity across the two models. These fixtures
    // have only planar polygons and complete circular rings, with no fitted trim.
    let mut comparison = exact.model.clone();
    for face in original.faces_of_body(body) {
        let carrier = comparison
            .geom_mut()
            .add_surface(original.geom().surface(original.face(face).surface).clone());
        assert_eq!(
            exact
                .model
                .faces_of_body(exact.result)
                .iter()
                .filter(|&&f| {
                    matches!(
                        comparison.geom().same_surface(
                            carrier,
                            comparison.face(f).surface,
                            &mut Decider::default()
                        ),
                        Outcome::ProvenCoincident(_)
                    )
                })
                .count(),
            1
        );
    }
    for edge in original.edges_of_body(body) {
        let curve = comparison
            .geom_mut()
            .add_curve(original.geom().curve(original.edge(edge).curve).clone());
        assert_eq!(
            exact
                .model
                .edges_of_body(exact.result)
                .iter()
                .filter(|&&e| {
                    matches!(
                        catallax_intersect::same_curve(
                            comparison.geom(),
                            curve,
                            comparison.edge(e).curve,
                            &mut Decider::default()
                        ),
                        Outcome::ProvenCoincident(_)
                    )
                })
                .count(),
            1
        );
    }
    for vertex in original.vertices_of_body(body) {
        let point = comparison
            .geom_mut()
            .add_point(original.geom().point(original.vertex(vertex).point).clone());
        assert!(exact.model.vertices_of_body(exact.result).iter().any(|&v| {
            matches!(
                comparison.geom().same_point(
                    point,
                    comparison.vertex(v).point,
                    &mut Decider::default()
                ),
                Outcome::ProvenCoincident(_)
            )
        }));
    }

    // SolidDraft auto uses 0.002 of the diagonal and a 15 degree normal control.
    // Exercise the app's coarse, standard and fine controls for this body's size.
    for factor in [4.0, 1.0, 0.25] {
        let options = MeshOptions {
            chord_tolerance: Rational::parse(&(0.002 * diagonal * factor).to_string()).unwrap(),
            angular_deviation: Some(Rational::parse("0.2617993877991494").unwrap()),
            ..MeshOptions::default()
        };
        let mesh = build_incremental_with_files(
            &definition,
            Some(built.checkpoint.clone()),
            &files,
            &options,
            &mut Decider::default(),
        )
        .unwrap()
        .display
        .value;
        assert_eq!(mesh.faces.len(), faces);
        assert!(
            mesh.faces
                .iter()
                .all(|face| face.max_chord_error <= options.chord_tolerance)
        );
        assert!(mesh.unavailable_controls.is_empty());
        for face in &mesh.faces {
            assert!(
                face.max_angular_deviation.as_ref().unwrap()
                    <= options.angular_deviation.as_ref().unwrap()
            );
        }
    }

    let mut output = Vec::new();
    let report = write_step(
        &exact.model,
        &[exact.result],
        &StepOptions::default(),
        &mut Decider::default(),
        &mut output,
    )
    .unwrap()
    .value;
    assert!(report.max_deviation <= StepOptions::default().curve_tolerance);
    let written = read_step(&output).unwrap();
    assert_eq!(written.solids.len(), 1);
    assert_eq!(written.solids[0].faces.len(), faces);
    let checked = apply(
        &mut Model::new(),
        &output,
        "reexport.step",
        &options,
        &[],
        None,
    )
    .unwrap();
    assert_eq!(checked.bodies.len(), 1);
    assert!(checked.left_open.is_empty());
}

#[test]
fn self_exported_planar_solid_imports_and_displays_without_movement() {
    let mut model = Model::new();
    let body = block(
        &mut model,
        [2, 3, 5].map(Rational::from_i64),
        &Pose::identity(),
    )
    .unwrap();
    roundtrip(&model, body, 6, 38.0_f64.sqrt());
}

#[test]
fn self_exported_cylinder_imports_and_displays_without_movement() {
    let mut model = Model::new();
    let body = cylinder(&mut model, 3.into(), 7.into(), &Pose::identity()).unwrap();
    roundtrip(&model, body, 3, 11.0);
}
