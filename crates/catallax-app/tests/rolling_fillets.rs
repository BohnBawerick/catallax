//! Fillets through the ordinary no-display application build: the exact
//! builders where they finish, the rolling ball where they refuse.
use catallax_app::{EvaluatedApplicationOutcome, build_incremental_with_files_without_display};
use catallax_mesh::{MeshOptions, PreviewOptions, preview_step};
use catallax_numbers::{Decider, OperationWork, Rational, with_operation_budget};
use catallax_recipe::Recipe;
use serde_json::{Value, json};

/// The owner's part: an 80 by 60 profile with R25 top corners, 76 long, and a
/// 34 by 22.67 pocket across the line where the flat top meets one round.
const OWNER: &str = include_str!("../../catallax-features/tests/fixtures/owner-curved-pocket.json");

const WHOLE_RIM: [(&str, &str); 6] = [
    ("top", "near"),
    ("top", "left"),
    ("top", "far"),
    ("right-round", "near"),
    ("right-round", "right"),
    ("right-round", "far"),
];

/// The owner's part with rim edges rounded, each named by the outer face and
/// the pocket wall it lies between.
fn rim(radius: &str, edges: &[(&str, &str)]) -> Recipe {
    let mut value: Value = serde_json::from_str(OWNER).unwrap();
    let edges: Vec<Value> = edges
        .iter()
        .map(|(outer, wall)| {
            json!({"between": [
                {"step": "blank", "role": {"side": {"segment": outer}}},
                {"step": "pocket-tool", "role": {"side": {"segment": wall}}}
            ]})
        })
        .collect();
    value["steps"].as_array_mut().unwrap().push(json!({
        "id": "round", "op": "fillet", "body": "pocket", "radius": radius, "edges": edges
    }));
    value["result"] = json!("round");
    Recipe::parse(&value.to_string()).unwrap()
}

/// [`rim`] with the fillet's `overflow` and `corner` options written out.
fn rim_with(radius: &str, edges: &[(&str, &str)], options: Value) -> Recipe {
    let mut value: Value = serde_json::from_str(&rim(radius, edges).to_text()).unwrap();
    let fillet = value["steps"].as_array_mut().unwrap().last_mut().unwrap();
    for (name, option) in options.as_object().unwrap() {
        fillet[name] = option.clone();
    }
    Recipe::parse(&value.to_string()).unwrap()
}

fn build(
    recipe: &Recipe,
) -> (
    Result<EvaluatedApplicationOutcome, catallax_app::Error>,
    OperationWork,
) {
    with_operation_budget(100_000_000, || {
        build_incremental_with_files_without_display(
            recipe,
            None,
            &Default::default(),
            &mut Decider::default(),
        )
    })
}

/// The display SolidDraft asks for: the default chord and 15 degrees between normals.
fn display() -> PreviewOptions {
    PreviewOptions::from(&MeshOptions {
        angular_deviation: Some(Rational::parse("0.2617993875").unwrap()),
        ..Default::default()
    })
}

#[test]
fn the_owner_pocket_rim_builds_exports_and_displays() {
    // The smallest ordinary radius and the largest that fits beside the wall.
    for radius in ["1", "3.99"] {
        let (outcome, work) = build(&rim(radius, &WHOLE_RIM));
        println!("pocket rim R{radius}: {work:?}");
        assert!(!work.exhausted, "R{radius}: {work:?}");
        let EvaluatedApplicationOutcome::Numerical(result) = outcome.unwrap() else {
            panic!("R{radius}: the exact builders refuse this rim");
        };
        // The exact attempt gave its charges back, and they stay on the record.
        assert!(work.refunded > 0, "R{radius}: {work:?}");
        let body = result.body();
        assert_eq!(body.kind(), catallax_topology::BodyKind::Solid);
        assert_eq!(body.definition().faces.len(), 19, "R{radius}");
        let blends = body
            .definition()
            .faces
            .iter()
            .filter(|face| {
                result.face_origins()[&face.key].iter().any(|origin| {
                    origin.step.as_str() == "round"
                        && matches!(origin.role, catallax_recipe::Role::Blend { .. })
                })
            })
            .count();
        assert_eq!(blends, 6, "R{radius}");
        // The app shows and exports the body through its numerical STEP file.
        let (step, export) = with_operation_budget(100_000_000, || {
            catallax_step::numerical::build(body, &Default::default(), &mut Decider::default())
        });
        let step = step.unwrap();
        assert!(!export.exhausted, "R{radius}: {export:?}");
        let (shown, preview) = with_operation_budget(100_000_000, || {
            preview_step(step.bytes(), &display(), &mut Decider::default())
        });
        let shown = shown.unwrap();
        println!(
            "pocket rim R{radius}: export {export:?}, display {preview:?}, {} triangles",
            shown.mesh.triangles.len()
        );
        assert!(!preview.exhausted, "R{radius}: {preview:?}");
        assert_eq!(shown.mesh.faces.len(), 19, "R{radius}");
        // A million vertices is the display's cap; a part this plain needs few.
        assert!(shown.mesh.positions.len() < 250_000, "R{radius}");
    }
}

#[test]
fn a_rim_radius_past_the_wall_names_the_face_it_reaches_and_the_limit() {
    let refusal = |recipe: &Recipe| {
        let (outcome, work) = build(recipe);
        // The limit is proved from the supports before any strip is built.
        assert!(work.used < 2_000_000, "{work:?}");
        let catallax_app::Error::Build(catallax_features::BuildError::Geometry { refusal, .. }) =
            outcome.unwrap_err()
        else {
            panic!("a radius limit is a geometric refusal");
        };
        refusal
    };
    // The whole rim: 8 mm of wall under the R25 side leaves room for R4.
    let whole = refusal(&rim("5", &WHOLE_RIM));
    assert_eq!(whole.code, "fillet.radius-too-large");
    assert!(whole.detail.contains("; reaches face "), "{}", whole.detail);
    assert!(
        whole.detail.ends_with("; radius limit 4 mm"),
        "{}",
        whole.detail
    );
    // Without the edge along that wall the two arcs end on it, up to R8.
    let five: Vec<_> = WHOLE_RIM
        .into_iter()
        .filter(|edge| *edge != ("right-round", "right"))
        .collect();
    let open = refusal(&rim("8", &five));
    assert_eq!(open.code, "fillet.radius-too-large");
    assert!(
        open.detail.ends_with("; radius limit 8 mm"),
        "{}",
        open.detail
    );
}

#[test]
fn a_refused_exact_attempt_is_counted_the_same_at_every_thread_count() {
    // Past the radius limit the exact builders refuse and so does the rolling
    // ball, early: a short request that still makes a refund.
    let recipe = rim("5", &WHOLE_RIM);
    let mut reference = None;
    for threads in [1, 2, 4] {
        catallax_numbers::set_threads(std::num::NonZeroUsize::new(threads));
        let (outcome, work) = build(&recipe);
        assert!(outcome.is_err(), "{threads} threads");
        assert!(!work.exhausted, "{threads} threads: {work:?}");
        assert!(work.refunded > 0, "{threads} threads: {work:?}");
        let counted = (work.used, work.refunded);
        assert_eq!(
            *reference.get_or_insert(counted),
            counted,
            "{threads} threads"
        );
    }
    catallax_numbers::set_threads(None);
}

#[test]
fn an_open_rim_end_is_cut_by_the_wall_it_runs_into() {
    let five: Vec<_> = WHOLE_RIM
        .into_iter()
        .filter(|edge| *edge != ("right-round", "right"))
        .collect();
    // An ordinary radius, and one a tenth of a millimetre under the 8 mm at
    // which the arcs no longer end on the wall. There the cut runs almost
    // along the strip and tapers it to a sliver.
    for radius in ["5", "7.9"] {
        let (outcome, work) = build(&rim(radius, &five));
        println!("five-edge rim R{radius}: {work:?}");
        assert!(!work.exhausted, "R{radius}: {work:?}");
        let EvaluatedApplicationOutcome::Numerical(result) = outcome.unwrap() else {
            panic!("the exact builders refuse an arc strip ending on a wall");
        };
        assert_eq!(result.body().definition().faces.len(), 18, "R{radius}");
        let (step, export) = with_operation_budget(100_000_000, || {
            catallax_step::numerical::build(
                result.body(),
                &Default::default(),
                &mut Decider::default(),
            )
        });
        let step = step.unwrap();
        assert!(!export.exhausted, "R{radius}: {export:?}");
        let (shown, preview) = with_operation_budget(100_000_000, || {
            preview_step(step.bytes(), &display(), &mut Decider::default())
        });
        println!("five-edge rim R{radius}: export {export:?}, display {preview:?}");
        assert_eq!(shown.unwrap().mesh.faces.len(), 18, "R{radius}");
        assert!(!preview.exhausted, "R{radius}: {preview:?}");
    }
}

#[test]
fn naming_the_default_overflow_and_corner_builds_the_same_rim() {
    let (plain, plain_work) = build(&rim("1", &WHOLE_RIM));
    let EvaluatedApplicationOutcome::Numerical(plain) = plain.unwrap() else {
        panic!("the exact builders refuse this rim");
    };
    for options in [
        json!({"overflow": "roll_over"}),
        json!({"corner": "cliff"}),
        json!({"overflow": "roll_over", "corner": "cliff"}),
    ] {
        let (named, work) = build(&rim_with("1", &WHOLE_RIM, options.clone()));
        let EvaluatedApplicationOutcome::Numerical(named) = named.unwrap() else {
            panic!("{options}: the exact builders refuse this rim");
        };
        assert_eq!(
            named.body().definition().faces.len(),
            plain.body().definition().faces.len(),
            "{options}"
        );
        assert_eq!(work.used, plain_work.used, "{options}");
        assert_eq!(work.refunded, plain_work.refunded, "{options}");
    }
    // A smooth corner is a construction the rolling ball does not have, and
    // it says so rather than building a cliff.
    let (smooth, _) = build(&rim_with("1", &WHOLE_RIM, json!({"corner": "smooth"})));
    let catallax_app::Error::Build(catallax_features::BuildError::Geometry { refusal, .. }) =
        smooth.unwrap_err()
    else {
        panic!("a corner the construction lacks is a geometric refusal");
    };
    assert!(
        refusal.code == "features.numerical-fillet-closure"
            || refusal.detail.contains("features.numerical-fillet-closure"),
        "{refusal:?}"
    );
}

fn boxed(steps: Vec<Value>, result: &str) -> Recipe {
    let mut all = vec![
        json!({"id": "profile", "op": "sketch",
        "plane": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]},
        "loops": [[
            {"kind": "line", "id": "a", "to": ["30", "-20"]},
            {"kind": "line", "id": "b", "to": ["30", "20"]},
            {"kind": "line", "id": "c", "to": ["-30", "20"]},
            {"kind": "line", "id": "d", "to": ["-30", "-20"]}
        ]]}),
        json!({"id": "box", "op": "extrude", "profile": "profile", "distance": "20"}),
    ];
    all.extend(steps);
    Recipe::parse(
        &json!({"format": "catallax-recipe/5", "name": "box", "steps": all, "result": result})
            .to_string(),
    )
    .unwrap()
}

#[test]
fn a_fillet_the_exact_builders_finish_stays_exact_and_takes_later_features() {
    // All twelve edges of a box: cylinders and spherical corners, exactly.
    let all = boxed(
        vec![json!({"id": "round", "op": "fillet", "body": "box", "radius": "5", "edges": "all"})],
        "round",
    );
    let (outcome, work) = build(&all);
    assert!(!work.exhausted, "{work:?}");
    // Nothing was refused, so nothing was refunded.
    assert_eq!(work.refunded, 0, "{work:?}");
    let EvaluatedApplicationOutcome::Exact(body) = outcome.unwrap() else {
        panic!("a box fillet is exact");
    };
    let built = body.checkpoint.built();
    assert_eq!(built.model.faces_of_body(built.result).len(), 26);
    // An exact fillet keeps its continuation: a later cut reads it.
    let cut = boxed(
        vec![
            json!({"id": "round", "op": "fillet", "body": "box", "radius": "5", "edges": [
                {"between": [
                    {"step": "box", "role": {"side": {"segment": "a"}}},
                    {"step": "box", "role": {"side": {"segment": "b"}}}
                ]}
            ]}),
            json!({"id": "hole-profile", "op": "sketch",
                "plane": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]},
                "loops": [[{"kind": "circle", "id": "rim", "centre": ["0", "0"], "radius": "4"}]]}),
            json!({"id": "hole-tool", "op": "extrude", "profile": "hole-profile", "distance": "20"}),
            json!({"id": "bored", "op": "boolean", "kind": "subtract", "target": "round", "tool": "hole-tool"}),
        ],
        "bored",
    );
    let (outcome, work) = build(&cut);
    assert!(!work.exhausted, "{work:?}");
    assert!(matches!(
        outcome.unwrap(),
        EvaluatedApplicationOutcome::Exact(_)
    ));
}

#[test]
fn a_near_point_picks_between_two_arcs_of_one_circle() {
    // The sample part's boss foot: the flange cuts its circle into two arcs
    // between the same two faces, and the pick names one by a point on it.
    let recipe = Recipe::parse(include_str!("fixtures/soliddraft-demo-boss-arc.json")).unwrap();
    let (outcome, work) = build(&recipe);
    println!("boss foot arc R1: {work:?}");
    assert!(!work.exhausted, "{work:?}");
    let EvaluatedApplicationOutcome::Numerical(result) = outcome.unwrap() else {
        panic!("the exact builders refuse an arc strip ending on the flange");
    };
    // Eighteen faces of the sample part and the one strip.
    assert_eq!(result.body().definition().faces.len(), 19);
}

#[test]
fn a_near_point_picks_between_two_straight_edges_of_one_cylinder() {
    // The sample part's boss leaves the end of the flange along two straight
    // edges, at x = -sqrt(44) and x = +sqrt(44): the same cylinder against the
    // same plane, ending on the same faces. The kernel stores each as a branch
    // of the surface pair, not as a line, and the pick names one by a sampled
    // point on it.
    let fixture = include_str!("fixtures/soliddraft-demo-boss-side.json");
    for (sign, side) in [("-", -1.0), ("", 1.0)] {
        let mut value: Value = serde_json::from_str(fixture).unwrap();
        let fillet = value["steps"].as_array_mut().unwrap().last_mut().unwrap();
        let near = fillet["edges"][0]["near"][0].as_str().unwrap().to_owned();
        fillet["edges"][0]["near"][0] = json!(format!("{sign}{}", near.trim_start_matches('-')));
        let recipe = Recipe::parse(&value.to_string()).unwrap();
        let (outcome, work) = build(&recipe);
        println!("boss side at x = {sign}sqrt(44), R1: {work:?}");
        assert!(!work.exhausted, "{work:?}");
        let EvaluatedApplicationOutcome::Numerical(result) = outcome.unwrap() else {
            panic!("the exact builders refuse a strip ending on the flange and the boss top");
        };
        let body = result.body().definition();
        // Eighteen faces of the sample part and the one strip.
        assert_eq!(body.faces.len(), 19);
        // The part is symmetric in x, so its two sharp edges ended on four
        // vertices at (-sqrt(44), 10) and (+sqrt(44), 10). The strip replaces
        // the two of the edge that was picked and leaves the other two.
        let sharp = |at: f64| {
            body.vertices
                .iter()
                .filter(|vertex| {
                    let [x, y, _] = vertex.position.each_ref().map(Rational::to_f64_nearest);
                    (x - at * 44.0f64.sqrt()).abs() < 1e-6 && (y - 10.0).abs() < 1e-6
                })
                .count()
        };
        assert_eq!([sharp(side), sharp(-side)], [0, 2], "x = {sign}sqrt(44)");
    }
}
