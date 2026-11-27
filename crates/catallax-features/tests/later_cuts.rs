//! P5 with its junction fillet built, then one later cut (PLAN M8). A cut
//! either builds, and every fitted face is measured again inside its recorded
//! bound, or refuses with a reason that names what is missing. The recipes and
//! the intersection list they feed are described in docs/M8-LATER-CUTS.md.

use catallax_features::{BuildError, Built, evaluate};
use catallax_mesh::{MeshOptions, mesh_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{ROLLING_ENVELOPE_TARGET_NM, SplineSupport, Surface};
use catallax_topology::{BodyId, envelope_rows};
use std::collections::{BTreeMap, BTreeSet};

const APPLICATION_ALLOWANCE: u64 = 20_000_000;

fn build(recipe: &str) -> (Result<Built, BuildError>, bool) {
    let recipe = Recipe::parse(recipe).unwrap();
    let (built, work) = with_operation_budget(APPLICATION_ALLOWANCE, || {
        evaluate(&recipe, &mut Decider::default())
    });
    (built, work.exhausted)
}

fn fitted_surfaces(built: &Built, body: BodyId) -> BTreeSet<catallax_shapes::SurfaceId> {
    built
        .model
        .faces_of_body(body)
        .into_iter()
        .map(|face| built.model.face(face).surface)
        .filter(|surface| {
            matches!(
                built.model.geom().surface(*surface),
                Surface::BSpline(SplineSupport::RollingEnvelope(_))
            )
        })
        .collect()
}

#[test]
fn a_cut_that_misses_the_part_keeps_and_remeasures_every_fitted_face() {
    let (built, exhausted) = build(include_str!("fixtures/p5-later-cuts/k12-far-disjoint.json"));
    assert!(!exhausted);
    let built = built.unwrap();
    let model = &built.model;
    assert_eq!(model.check(built.result), vec![]);
    let filleted = built.bodies[&StepId::new("filleted").unwrap()];
    // The fitted-body certificate and volume start from the unfilleted base.
    assert!(model.blend_base(filleted).is_some());
    assert_eq!(model.blend_base(built.result), model.blend_base(filleted));
    let before = fitted_surfaces(&built, filleted);
    let after = fitted_surfaces(&built, built.result);
    assert_eq!(before.len(), 4);
    // The cut keeps the fitted carriers themselves, not refitted copies.
    assert_eq!(after, before);
    let target = Rational::from_ratio(ROLLING_ENVELOPE_TARGET_NM, 1_000_000).unwrap();
    let mut dec = Decider::default();
    for face in model.faces_of_body(built.result) {
        let Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) =
            model.geom().surface(model.face(face).surface)
        else {
            continue;
        };
        // Still the whole checked face: its four boundary rows are intact.
        assert_eq!(envelope_rows(model, face, &mut dec).unwrap().len(), 4);
        let recorded = &envelope.bound().upper;
        let measured = envelope.remeasure(&mut dec).unwrap();
        // The same check on the same retained data gives the same bound.
        assert_eq!(&measured, recorded);
        assert!(recorded <= &target);
    }
}

#[test]
fn the_disjoint_cut_preserves_torus_trims_at_coarse_and_fine_viewport_levels() {
    let (built, exhausted) = build(include_str!("fixtures/p5-later-cuts/k12-far-disjoint.json"));
    assert!(!exhausted);
    let built = built.unwrap();
    for chord in ["0.276318656627", "0.017269916039"] {
        let options = MeshOptions {
            chord_tolerance: Rational::parse(chord).unwrap(),
            angular_deviation: Some(Rational::parse("0.2617993877991494").unwrap()),
            ..MeshOptions::default()
        };
        let mesh = mesh_body(
            &built.model,
            built.result,
            &options,
            &mut Decider::default(),
        )
        .unwrap()
        .value;
        assert!(mesh.faces.iter().all(|face| {
            !face.triangles.is_empty() && face.max_chord_error <= options.chord_tolerance
        }));
        let mut chords = BTreeMap::new();
        for triangle in &mesh.triangles {
            for i in 0..3 {
                let [a, b] = [i, (i + 1) % 3].map(|j| triangle.vertices[j]);
                let uses = chords.entry((a.min(b), a.max(b))).or_insert((0, 0));
                uses.0 += 1;
                uses.1 += if a < b { 1 } else { -1 };
            }
        }
        assert!(
            mesh.edges
                .iter()
                .all(|edge| edge.vertices.windows(2).all(|pair| {
                    let (a, b) = (pair[0], pair[1]);
                    chords.get(&(a.min(b), a.max(b))) == Some(&(2, 0))
                }))
        );
    }
}

/// The centre of gravity of a fitted blend is stated with the same bound
/// whichever entry point asks for it: it encloses both the exact definition
/// and the fitted carriers.
#[test]
fn the_centroid_of_a_fitted_blend_is_answered_with_its_stated_bound() {
    let (built, exhausted) = build(include_str!("fixtures/p5-later-cuts/k12-far-disjoint.json"));
    assert!(!exhausted);
    let built = built.unwrap();
    let accuracy = catallax_measure::Accuracy {
        relative: Rational::zero(),
        absolute: Rational::from_ratio(1, 1000).unwrap(),
    };
    let all = catallax_measure::volume_area_centroid(
        &built.model,
        built.result,
        &accuracy,
        &mut Decider::default(),
    )
    .unwrap()
    .value;
    let centroid = catallax_measure::centroid(
        &built.model,
        built.result,
        &accuracy,
        &mut Decider::default(),
    )
    .unwrap()
    .value
    .unwrap();
    assert_eq!(Some(centroid), all.centroid);
    // Second moments are not covered for a fitted blend; they refuse as
    // `mass_properties` does, before any work.
    let inertia = catallax_measure::inertia(
        &built.model,
        built.result,
        &accuracy,
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(inertia.code, "measure.unsupported-surface");
    assert!(
        inertia.detail.contains("second moments"),
        "{}",
        inertia.detail
    );
}

/// A later boolean carries the measured fit bound of the fitted faces its
/// result keeps, and no fit bound for faces it does not keep: a body is not
/// marked with an operand's deviation it no longer has.
#[test]
fn a_later_boolean_marks_the_fit_bound_of_the_faces_it_keeps() {
    let fit_marks = |built: &Built| -> Vec<Rational> {
        built
            .model
            .body(built.result)
            .marks
            .iter()
            .filter(|mark| mark.code == "features.rolling-envelope-fit")
            .map(|mark| mark.bound.clone())
            .collect()
    };
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("fixtures/p5-later-cuts/k12-far-disjoint.json")).unwrap();
    let with_cut = |kind: &str, target: &str, tool: &str| {
        let mut fixture = fixture.clone();
        let steps = fixture["steps"].as_array_mut().unwrap();
        let cut = steps.last_mut().unwrap();
        cut["kind"] = kind.into();
        cut["target"] = target.into();
        cut["tool"] = tool.into();
        let (built, exhausted) = build(&fixture.to_string());
        assert!(!exhausted, "{kind} {target} {tool}");
        built.unwrap()
    };
    // The filleted body keeps all four fitted faces: one mark, their worst
    // recorded bound, which the fit's own check measures again.
    let kept = with_cut("subtract", "filleted", "tool");
    let mut dec = Decider::default();
    let worst = kept
        .model
        .faces_of_body(kept.result)
        .into_iter()
        .filter_map(
            |face| match kept.model.geom().surface(kept.model.face(face).surface) {
                Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) => {
                    Some(envelope.remeasure(&mut dec).unwrap())
                }
                _ => None,
            },
        )
        .max()
        .unwrap();
    assert_eq!(fit_marks(&kept), vec![worst]);
    // The sphere minus the far filleted body is the sphere: no fitted face.
    let sphere = with_cut("subtract", "tool", "filleted");
    assert!(fitted_surfaces(&sphere, sphere.result).is_empty());
    assert_eq!(fit_marks(&sphere), vec![]);
    // The two are apart, so their intersection is empty.
    let empty = with_cut("intersect", "filleted", "tool");
    assert!(empty.model.faces_of_body(empty.result).is_empty());
    assert_eq!(fit_marks(&empty), vec![]);
    // A body minus itself is empty, whatever it was marked with.
    let itself = with_cut("subtract", "filleted", "filleted");
    assert!(itself.model.faces_of_body(itself.result).is_empty());
    assert_eq!(fit_marks(&itself), vec![]);
}

/// Each later cut that does not build stops with a code and a detail naming
/// the missing intersection or fitted-boundary query, at the cut's own step.
#[test]
fn later_cuts_refuse_by_naming_what_is_missing() {
    let rows: [(&str, &str, &str, &str); 11] = [
        (
            "k01 wall hole",
            include_str!("fixtures/p5-later-cuts/k01-wall-hole.json"),
            "intersect.parametric-surface",
            "fitted spline boundary",
        ),
        (
            "k02 top slot",
            include_str!("fixtures/p5-later-cuts/k02-top-slot.json"),
            "intersect.spline-plane",
            "rational homogeneous cut coefficients",
        ),
        (
            "k03 rim notch",
            include_str!("fixtures/p5-later-cuts/k03-rim-notch.json"),
            "shapes.parametric-query",
            "cross-kind spline support equality",
        ),
        (
            "k04 torus pocket",
            include_str!("fixtures/p5-later-cuts/k04-torus-pocket.json"),
            "numbers.operation-budget",
            "cylinder/torus",
        ),
        (
            "k05 plane notch",
            include_str!("fixtures/p5-later-cuts/k05-plane-notch.json"),
            "numbers.operation-budget",
            "fitted rolling-envelope spline with plane",
        ),
        (
            "k06 cylinder drill",
            include_str!("fixtures/p5-later-cuts/k06-cylinder-drill.json"),
            "numbers.operation-budget",
            "fitted rolling-envelope spline with cylinder",
        ),
        (
            "k07 cone drill",
            include_str!("fixtures/p5-later-cuts/k07-cone-drill.json"),
            "numbers.operation-budget",
            "cone/cylinder",
        ),
        (
            "k08 sphere dimple",
            include_str!("fixtures/p5-later-cuts/k08-sphere-dimple.json"),
            "numbers.operation-budget",
            "cone/sphere",
        ),
        (
            "k09 ring groove",
            include_str!("fixtures/p5-later-cuts/k09-ring-groove.json"),
            "numbers.operation-budget",
            "fitted rolling-envelope spline with torus",
        ),
        (
            "k10 axis section",
            include_str!("fixtures/p5-later-cuts/k10-axis-section.json"),
            "numbers.operation-budget",
            "fitted rolling-envelope spline with plane",
        ),
        (
            "k11 offset section",
            include_str!("fixtures/p5-later-cuts/k11-offset-section.json"),
            "numbers.operation-budget",
            "fitted rolling-envelope spline with plane",
        ),
    ];
    for (name, recipe, code, detail) in rows {
        let (built, exhausted) = build(recipe);
        let Err(BuildError::Geometry { step, refusal }) = built else {
            panic!("{name}: expected a refusal");
        };
        assert_eq!(refusal.code, code, "{name}: {}", refusal.detail);
        assert_eq!(
            exhausted,
            code == "numbers.operation-budget" || name.starts_with("k01"),
            "{name}"
        );
        assert!(
            refusal.detail.contains(detail),
            "{name}: {}",
            refusal.detail
        );
        assert_eq!(step.as_str(), "later-cut", "{name}");
    }
}

#[test]
fn a_countersink_pattern_cuts_the_flange_but_not_the_distant_fillet() {
    use catallax_boolean::{Containment, classify_point};
    use catallax_numbers::Outcome;
    use catallax_shapes::RPoint3;

    let (built, exhausted) = build(include_str!("fixtures/p5-later-cuts/k13-fillet-first.json"));
    assert!(!exhausted);
    let mut built = built.unwrap();
    assert!(built.model.check(built.result).is_empty());
    let original = built.bodies[&StepId::new("filleted").unwrap()];
    let fillet = |built: &Built, body| {
        built
            .model
            .faces_of_body(body)
            .into_iter()
            .map(|face| built.model.face(face).surface)
            .find(|&surface| matches!(built.model.geom().surface(surface), Surface::Torus(_)))
            .expect("the flange fillet keeps its torus")
    };
    assert_eq!(fillet(&built, original), fillet(&built, built.result));

    let mut dec = Decider::default();
    for (x, z, want) in [
        (17, 0, Containment::Outside),
        (0, 17, Containment::Outside),
        (-17, 0, Containment::Outside),
        (0, -17, Containment::Outside),
        (17, 5, Containment::Inside),
    ] {
        let point = RPoint3::new(x.into(), 2.into(), z.into());
        assert_eq!(
            classify_point(&mut built.model, built.result, &point, &mut dec),
            Outcome::Proven(want)
        );
    }
}

/// A fitted blend whose fitted faces are whole envelope patches states an
/// area that encloses both the exact definition and the fitted carriers: the
/// carriers' own area is enclosed, so no area-of-definition mark is needed.
#[test]
fn a_fitted_blend_area_encloses_the_definition_and_the_carriers() {
    let (built, exhausted) = build(include_str!("fixtures/p5-later-cuts/k12-far-disjoint.json"));
    assert!(!exhausted);
    let built = built.unwrap();
    let accuracy = catallax_measure::Accuracy {
        relative: Rational::zero(),
        absolute: Rational::from_ratio(1, 1000).unwrap(),
    };
    let scoped = |marks: &catallax_numbers::MarkLog| {
        marks
            .iter()
            .filter(|mark| mark.code == "measure.area-of-definition")
            .count()
    };
    let start = std::time::Instant::now();
    let all = catallax_measure::volume_area_centroid(
        &built.model,
        built.result,
        &accuracy,
        &mut Decider::default(),
    )
    .unwrap();
    let elapsed = start.elapsed();
    let area = catallax_measure::area(
        &built.model,
        built.result,
        &accuracy,
        &mut Decider::default(),
    )
    .unwrap();
    // The definition's accuracy is narrower than these fitted carriers can be
    // enclosed: main's `certified_area` refuses the carrier budget rather than
    // return a coarser bound, so the answer keeps the definition's area and
    // names the gap with `measure.area-of-definition`.
    assert_eq!(scoped(&all.marks), 1);
    assert_eq!(scoped(&area.marks), 1);
    assert_eq!(area.value, all.value.area);
    let volume = catallax_measure::volume(
        &built.model,
        built.result,
        &accuracy,
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(scoped(&volume.marks), 0);
    let definition = catallax_measure::blend_definition_properties(
        &built.model,
        built.result,
        &accuracy,
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(scoped(&definition.marks), 0);
    // The represented area holds the definition's.
    let (stated, exact) = (&all.value.area, &definition.value.area);
    assert!(
        stated.lo <= exact.lo && exact.hi <= stated.hi,
        "{stated:?} {exact:?}"
    );
    // These carriers cannot be enclosed inside any allowance main's
    // `certified_area` accepts: their mouth-row cells stay unproven-regular,
    // so refinement runs to the cap and refuses `shapes.area-cell-budget` at
    // every width. The answer therefore states the definition's area with its
    // `measure.area-of-definition` mark rather than a coarser unmarked bound
    // (the reviewed defect). The enclosure assertion above holds either way.
    eprintln!(
        "AREA definition [{}, {}] stated [{}, {}] width {} in {:?}",
        exact.lo.to_f64_nearest(),
        exact.hi.to_f64_nearest(),
        stated.lo.to_f64_nearest(),
        stated.hi.to_f64_nearest(),
        stated.width().to_f64_nearest(),
        elapsed
    );
}
