//! Recipe loft cases: exact laws, physical topology and application consumers.

use catallax_app::ApplicationBody;
use catallax_features::BuildError;
use catallax_measure::{Accuracy, Enclosure};
use catallax_mesh::MeshOptions;
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::{Recipe, Role};
use serde_json::{Value, json};
use std::collections::BTreeSet;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn accuracy() -> Accuracy {
    Accuracy {
        relative: Rational::zero(),
        absolute: q("1e-6"),
    }
}
fn options() -> MeshOptions {
    MeshOptions {
        chord_tolerance: q("1e-2"),
        max_vertices: 200_000,
        ..Default::default()
    }
}

fn polygon(prefix: &str, points: &[(&str, &str)]) -> Value {
    json!(
        points
            .iter()
            .enumerate()
            .map(|(i, (x, y))| json!({"kind":"line","id":format!("{prefix}_{i}"),"to":[x,y]}))
            .collect::<Vec<_>>()
    )
}
fn square(prefix: &str, half: &str) -> Value {
    let minus = format!("-{half}");
    polygon(
        prefix,
        &[
            (half, &minus),
            (half, half),
            (&minus, half),
            (&minus, &minus),
        ],
    )
}
fn circle(id: &str, radius: &str) -> Value {
    json!([{"kind":"circle","id":id,"centre":["0","0"],"radius":radius}])
}
fn spline(prefix: &str, cubic: bool) -> Value {
    json!([
        {"kind":"bezier","id":format!("{prefix}_roof"),"controls":if cubic {json!([["0","4"],["4","4"]])} else {json!([["2","4"]])},"to":["4","0"]},
        {"kind":"line","id":format!("{prefix}_floor"),"to":["0","0"]}
    ])
}
fn definition(loops: Vec<Value>, heights: &[&str], shape: &str) -> Value {
    let mut steps = Vec::new();
    let mut sections = Vec::new();
    for (i, (loops, height)) in loops.into_iter().zip(heights).enumerate() {
        let start: Vec<_> = loops
            .as_array()
            .unwrap()
            .iter()
            .map(|lp| lp[0]["id"].clone())
            .collect();
        steps.push(json!({"id":format!("section{i}"),"op":"sketch","plane":{"rotation":["1","0","0","0"],"translation":["0","0",height]},"loops":loops}));
        sections.push(json!({"profile":format!("section{i}"),"start":start}));
    }
    steps.push(json!({"id":"loft","op":"loft","sections":sections,"shape":shape}));
    json!({"format":"catallax-recipe/5","name":"rational loft cases","steps":steps,"result":"loft"})
}
fn recipe(value: &Value) -> Recipe {
    Recipe::parse(&value.to_string()).unwrap()
}
fn build(value: &Value) -> ApplicationBody {
    catallax_app::build(&recipe(value), &options(), &mut Decider::default()).unwrap()
}
fn encloses(actual: &Enclosure, expected: &str) {
    assert!(
        actual.contains(&q(expected)),
        "{actual:?} excludes {expected}"
    );
    assert!(
        actual.width() <= accuracy().absolute,
        "unmet volume/centroid width: {actual:?}"
    );
}
fn encloses_pi(actual: &Enclosure, factor: &str) {
    let expected = Ball::pi(256) * Ball::from_rational(&q(factor), 256);
    assert!(
        actual.contains(&expected.lo()) && actual.contains(&expected.hi()),
        "{actual:?} excludes {factor} pi"
    );
    assert!(
        actual.width() <= accuracy().absolute,
        "unmet volume width: {actual:?}"
    );
}
fn picking_is_physical(app: &ApplicationBody) {
    let model = &app.built.model;
    assert!(model.check(app.built.result).is_empty());
    assert_eq!(
        app.display
            .value
            .faces
            .iter()
            .map(|f| f.face)
            .collect::<BTreeSet<_>>(),
        app.names.faces.keys().copied().collect()
    );
    assert_eq!(
        app.display
            .value
            .edges
            .iter()
            .map(|e| e.edge)
            .collect::<BTreeSet<_>>(),
        app.names.edges.keys().copied().collect()
    );
    for triangle in &app.display.value.triangles {
        assert!(app.names.faces.contains_key(&triangle.face));
    }
    for face in &app.display.value.faces {
        assert!(face.max_chord_error <= options().chord_tolerance);
    }
}

fn step_regions_are_physical(app: &ApplicationBody) {
    let mut bytes = Vec::new();
    let options = catallax_step::StepOptions::default();
    let report = catallax_step::write_step(
        &app.built.model,
        &[app.built.result],
        &options,
        &mut Decider::default(),
        &mut bytes,
    )
    .unwrap()
    .value;
    assert!(report.max_deviation <= options.curve_tolerance);
    let parsed = catallax_import::read_step(&bytes).unwrap();
    assert_eq!(
        parsed
            .solids
            .iter()
            .map(|solid| solid.faces.len())
            .sum::<usize>(),
        report.faces
    );
    // The importer keys faces by #entity id, not their ADVANCED_FACE label.
    let text = String::from_utf8(bytes).unwrap();
    let file_faces: BTreeSet<_> = text
        .split("ADVANCED_FACE")
        .skip(1)
        .map(|tail| {
            tail.trim_start()
                .strip_prefix('(')
                .unwrap()
                .trim_start()
                .strip_prefix('\'')
                .unwrap()
                .split('\'')
                .next()
                .unwrap()
                .to_owned()
        })
        .collect();
    assert_eq!(file_faces.len(), report.faces);
    assert_eq!(report.face_map.len(), report.faces);
    assert_eq!(
        report
            .face_map
            .iter()
            .map(|region| region.file_region.clone())
            .collect::<BTreeSet<_>>(),
        file_faces
    );
    assert_eq!(
        report
            .face_map
            .iter()
            .map(|region| region.face)
            .collect::<BTreeSet<_>>(),
        app.names.faces.keys().copied().collect()
    );
    for region in &report.face_map {
        if let Some(domain) = &region.file_domain {
            assert!(domain.iter().all(|[lo, hi]| lo < hi));
        }
    }
}

#[test]
fn two_section_ruled_and_smooth_frustums_have_exact_mass_and_physical_topology() {
    for shape in ["ruled", "smooth"] {
        let app = build(&definition(
            vec![json!([square("a", "1")]), json!([square("b", "2")])],
            &["0", "3"],
            shape,
        ));
        picking_is_physical(&app);
        assert_eq!(
            (
                app.names.faces.len(),
                app.names.edges.len(),
                app.names.vertices.len()
            ),
            (6, 12, 8)
        );
        let mass = app
            .mass_properties(&accuracy(), &mut Decider::default())
            .unwrap();
        encloses(&mass.value.volume, "28");
        let center = mass.value.centroid.unwrap();
        encloses(&center[0], "0");
        encloses(&center[1], "0");
        encloses(&center[2], "51/28");
    }
}

fn circles(radii: &[&str], heights: &[&str]) -> Value {
    let loops = radii
        .iter()
        .enumerate()
        .map(|(i, radius)| json!([circle(&format!("c{i}"), radius)]))
        .collect();
    definition(loops, heights, "smooth")
}
/// Concentric circle sections give a solid of revolution. The natural spline
/// is linear in its station data, so section `v` is the circle of radius
/// `r(v)` at height `z(v)`, both natural cubics through the authored values.
/// The volume is `pi` times the integral of `r^2 z'`, which is positive here.
fn circular_loft_has_exact_mass(value: &Value, volume: &str, height: &str) {
    // The application entry applies its unchanged 20M allowance and reports
    // `ResourceLimit` when it is exhausted.
    let app = build(value);
    picking_is_physical(&app);
    assert_eq!(
        (
            app.names.faces.len(),
            app.names.edges.len(),
            app.names.vertices.len()
        ),
        (3, 2, 0)
    );
    let mass = app
        .mass_properties(&accuracy(), &mut Decider::default())
        .unwrap();
    encloses_pi(&mass.value.volume, volume);
    let center = mass.value.centroid.unwrap();
    encloses(&center[0], "0");
    encloses(&center[1], "0");
    encloses(&center[2], height);
    step_regions_are_physical(&app);
}

#[test]
fn two_circles_with_different_radii_build_inside_the_default_allowance() {
    circular_loft_has_exact_mass(&circles(&["1", "2"], &["0", "3"]), "7", "51/28");
}

#[test]
fn natural_three_square_loft_builds_inside_the_default_allowance() {
    let app = build(&definition(
        vec![
            json!([square("a", "1")]),
            json!([square("b", "2")]),
            json!([square("c", "1")]),
        ],
        &["0", "3", "6"],
        "smooth",
    ));
    picking_is_physical(&app);
    assert_eq!(
        (
            app.names.faces.len(),
            app.names.edges.len(),
            app.names.vertices.len()
        ),
        (6, 12, 8)
    );
    let mass = app
        .mass_properties(&accuracy(), &mut Decider::default())
        .unwrap();
    encloses(&mass.value.volume, "2298/35");
    let center = mass.value.centroid.unwrap();
    encloses(&center[0], "0");
    encloses(&center[1], "0");
    encloses(&center[2], "3");
}

#[test]
fn natural_three_circle_loft_builds_inside_the_default_allowance() {
    // The saved public recipe: radii 1, 2, 1 at heights 0, 1, 2.
    let saved: Value = serde_json::from_str(include_str!(
        "../../catallax-features/tests/fixtures/loft-natural-circles.json"
    ))
    .unwrap();
    circular_loft_has_exact_mass(&saved, "383/70", "1");
}

#[test]
fn unequally_spaced_natural_circles_build_inside_the_default_allowance() {
    // Uniform integer station spacing is retained: the stations are unevenly
    // placed in space, and no chord-length interpolation is substituted.
    circular_loft_has_exact_mass(
        &circles(&["1", "2", "1"], &["0", "2", "5"]),
        "383/28",
        "43411/18384",
    );
}

#[test]
fn quadratic_and_cubic_profiles_keep_exact_volume_and_area_centroid() {
    for (cubic, volume, y) in [(false, "16", "4/5"), (true, "144/5", "9/7")] {
        let app = build(&definition(
            vec![json!([spline("a", cubic)]), json!([spline("b", cubic)])],
            &["0", "3"],
            "ruled",
        ));
        picking_is_physical(&app);
        step_regions_are_physical(&app);
        assert_eq!(
            (
                app.names.faces.len(),
                app.names.edges.len(),
                app.names.vertices.len()
            ),
            (4, 6, 4)
        );
        let mass = app
            .mass_properties(&accuracy(), &mut Decider::default())
            .unwrap();
        encloses(&mass.value.volume, volume);
        let center = mass.value.centroid.unwrap();
        encloses(&center[0], "2");
        encloses(&center[1], y);
        encloses(&center[2], "3/2");
    }
}

#[test]
fn later_segment_splits_add_no_names_edges_or_vertices() {
    let top = polygon(
        "b",
        &[
            ("0", "-1"),
            ("1", "-1"),
            ("1", "0"),
            ("1", "1"),
            ("0", "1"),
            ("-1", "1"),
            ("-1", "0"),
            ("-1", "-1"),
        ],
    );
    let app = build(&definition(
        vec![json!([square("a", "1")]), json!([top])],
        &["0", "2"],
        "ruled",
    ));
    picking_is_physical(&app);
    assert_eq!(
        (
            app.names.faces.len(),
            app.names.edges.len(),
            app.names.vertices.len()
        ),
        (6, 12, 8)
    );
    let side_ids: BTreeSet<_> = app
        .names
        .faces
        .values()
        .flatten()
        .filter_map(|o| match &o.role {
            Role::Side { segment, .. } => Some(segment.as_str()),
            _ => None,
        })
        .collect();
    assert_eq!(side_ids, BTreeSet::from(["a_0", "a_1", "a_2", "a_3"]));
    encloses(
        &catallax_measure::volume(
            &app.built.model,
            app.built.result,
            &accuracy(),
            &mut Decider::default(),
        )
        .unwrap()
        .value,
        "8",
    );
}

#[test]
fn first_section_collinear_ids_survive_on_maximal_logical_faces() {
    let bottom = polygon(
        "a",
        &[
            ("0", "-1"),
            ("1", "-1"),
            ("1", "0"),
            ("1", "1"),
            ("0", "1"),
            ("-1", "1"),
            ("-1", "0"),
            ("-1", "-1"),
        ],
    );
    let app = build(&definition(
        vec![json!([bottom]), json!([square("b", "1")])],
        &["0", "2"],
        "smooth",
    ));
    assert_eq!(
        (
            app.names.faces.len(),
            app.names.edges.len(),
            app.names.vertices.len()
        ),
        (6, 12, 8)
    );
    let side_sets: Vec<_> = app
        .names
        .faces
        .values()
        .filter(|origins| origins.iter().any(|o| matches!(o.role, Role::Side { .. })))
        .collect();
    assert_eq!(
        side_sets
            .iter()
            .map(|origins| origins.len())
            .collect::<Vec<_>>(),
        vec![2, 2, 2, 2]
    );
    let ids: BTreeSet<_> = side_sets
        .into_iter()
        .flatten()
        .filter_map(|o| match &o.role {
            Role::Side { segment, .. } => Some(segment.as_str()),
            _ => None,
        })
        .collect();
    assert_eq!(
        ids,
        BTreeSet::from(["a_0", "a_1", "a_2", "a_3", "a_4", "a_5", "a_6", "a_7"])
    );
}

#[test]
fn nonaffine_parallel_sections_use_a_whole_polar_embedding_proof() {
    let top = polygon("b", &[("1", "0"), ("0", "2"), ("-2", "0"), ("0", "-1")]);
    let app = build(&definition(
        vec![json!([square("a", "1")]), json!([top])],
        &["0", "2"],
        "ruled",
    ));
    picking_is_physical(&app);
    encloses(
        &catallax_measure::volume(
            &app.built.model,
            app.built.result,
            &accuracy(),
            &mut Decider::default(),
        )
        .unwrap()
        .value,
        "23/3",
    );
}

#[test]
fn nonparallel_affine_sections_keep_the_spatial_volume_map() {
    let mut value = definition(
        vec![json!([square("a", "1")]), json!([square("b", "1")])],
        &["0", "4"],
        "ruled",
    );
    value["steps"][1]["plane"]["rotation"] = json!(["1", "0", "1/10", "0"]);
    let app = build(&value);
    picking_is_physical(&app);
    encloses(
        &catallax_measure::volume(
            &app.built.model,
            app.built.result,
            &accuracy(),
            &mut Decider::default(),
        )
        .unwrap()
        .value,
        "1600/101",
    );
}

#[test]
fn nested_holes_and_islands_remain_separate_material_lumps() {
    let app = build(&definition(
        vec![
            json!([
                circle("a_outer", "4"),
                circle("a_hole", "2"),
                circle("a_island", "1")
            ]),
            json!([
                circle("b_outer", "4"),
                circle("b_hole", "2"),
                circle("b_island", "1")
            ]),
        ],
        &["0", "3"],
        "ruled",
    ));
    picking_is_physical(&app);
    assert_eq!(app.built.model.body(app.built.result).lumps.len(), 2);
    assert_eq!(
        (
            app.names.faces.len(),
            app.names.edges.len(),
            app.names.vertices.len()
        ),
        (7, 6, 0)
    );
    encloses_pi(
        &catallax_measure::volume(
            &app.built.model,
            app.built.result,
            &accuracy(),
            &mut Decider::default(),
        )
        .unwrap()
        .value,
        "39",
    );
}

#[test]
fn coplanar_sections_regularize_to_empty_without_sheet_handles() {
    let app = build(&definition(
        vec![json!([square("a", "1")]), json!([square("b", "2")])],
        &["0", "0"],
        "ruled",
    ));
    assert_eq!(
        app.built.model.body(app.built.result).kind,
        catallax_topology::BodyKind::Empty
    );
    assert_eq!(
        (
            app.names.faces.len(),
            app.names.edges.len(),
            app.names.vertices.len()
        ),
        (0, 0, 0)
    );
    encloses(
        &catallax_measure::volume(
            &app.built.model,
            app.built.result,
            &accuracy(),
            &mut Decider::default(),
        )
        .unwrap()
        .value,
        "0",
    );
}

#[test]
fn reversed_correspondence_is_a_typed_refusal_not_a_twisted_approximation() {
    let top = polygon("b", &[("-1", "1"), ("1", "1"), ("1", "-1"), ("-1", "-1")]);
    let value = definition(
        vec![json!([square("a", "1")]), json!([top])],
        &["0", "2"],
        "ruled",
    );
    let error =
        catallax_app::build(&recipe(&value), &options(), &mut Decider::default()).unwrap_err();
    assert!(
        matches!(error, catallax_app::Error::Build(BuildError::Geometry { refusal, .. }) if refusal.code == "loft.admissibility-unproven")
    );
}

#[test]
fn legacy_arc_length_lofts_are_not_reinterpreted_as_rational_charts() {
    let mut value = definition(
        vec![json!([square("a", "1")]), json!([square("b", "1")])],
        &["0", "2"],
        "ruled",
    );
    value["format"] = json!("catallax-recipe/3");
    let error =
        catallax_app::build(&recipe(&value), &options(), &mut Decider::default()).unwrap_err();
    assert!(
        matches!(error, catallax_app::Error::Build(BuildError::Geometry { refusal, .. }) if refusal.code == "loft.legacy-parameterization")
    );
}
