//! P5 with its junction fillet built, then one later cut (PLAN M8). A cut
//! either builds, and every fitted face is measured again inside its recorded
//! bound, or refuses with a reason that names what is missing. The recipes and
//! the intersection list they feed are described in docs/M8-LATER-CUTS.md.

use catallax_features::{BuildError, Built, evaluate};
use catallax_numbers::{Decider, Rational, with_operation_budget};
use catallax_recipe::{Recipe, StepId};
use catallax_shapes::{ROLLING_ENVELOPE_TARGET_NM, SplineSupport, Surface};
use catallax_topology::{BodyId, envelope_rows};
use std::collections::BTreeSet;

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

/// Each later cut that does not build stops with a code and a detail naming
/// the missing intersection or the exhausted allowance, at the cut's own step.
#[test]
fn later_cuts_refuse_by_naming_what_is_missing() {
    let rows: [(&str, &str, &str, &str); 11] = [
        (
            "k02 top slot",
            include_str!("fixtures/p5-later-cuts/k02-top-slot.json"),
            "numbers.operation-budget",
            "allowance",
        ),
        (
            "k03 rim notch",
            include_str!("fixtures/p5-later-cuts/k03-rim-notch.json"),
            "boolean.parametric-query",
            "outgoing germ of a spline edge",
        ),
        (
            "k04 torus pocket",
            include_str!("fixtures/p5-later-cuts/k04-torus-pocket.json"),
            "intersect.traced-pair-required",
            "skew cylinder/torus",
        ),
        (
            "k05 plane notch",
            include_str!("fixtures/p5-later-cuts/k05-plane-notch.json"),
            "intersect.unsupported-pair",
            "fitted rolling-envelope spline with plane",
        ),
        (
            "k06 cylinder drill",
            include_str!("fixtures/p5-later-cuts/k06-cylinder-drill.json"),
            "intersect.unsupported-pair",
            "fitted rolling-envelope spline with cylinder",
        ),
        (
            "k07 cone drill",
            include_str!("fixtures/p5-later-cuts/k07-cone-drill.json"),
            "intersect.traced-pair-required",
            "non-coaxial cone/cylinder",
        ),
        (
            "k08 sphere dimple",
            include_str!("fixtures/p5-later-cuts/k08-sphere-dimple.json"),
            "intersect.traced-pair-required",
            "non-coaxial cone/sphere",
        ),
        (
            "k09 ring groove",
            include_str!("fixtures/p5-later-cuts/k09-ring-groove.json"),
            "intersect.unsupported-pair",
            "fitted rolling-envelope spline with torus",
        ),
        (
            "k10 axis section",
            include_str!("fixtures/p5-later-cuts/k10-axis-section.json"),
            "intersect.unsupported-pair",
            "fitted rolling-envelope spline with plane",
        ),
        (
            "k11 offset section",
            include_str!("fixtures/p5-later-cuts/k11-offset-section.json"),
            "intersect.unsupported-pair",
            "fitted rolling-envelope spline with plane",
        ),
        (
            "k13 fillet first",
            include_str!("fixtures/p5-later-cuts/k13-fillet-first.json"),
            "numbers.operation-budget",
            "allowance",
        ),
    ];
    for (name, recipe, code, detail) in rows {
        let (built, _) = build(recipe);
        let Err(BuildError::Geometry { step, refusal }) = built else {
            panic!("{name}: expected a refusal");
        };
        assert_eq!(refusal.code, code, "{name}: {}", refusal.detail);
        assert!(
            refusal.detail.contains(detail),
            "{name}: {}",
            refusal.detail
        );
        let expected = if name.starts_with("k13") {
            "drilled"
        } else {
            "later-cut"
        };
        assert_eq!(step.as_str(), expected, "{name}");
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
    assert_eq!(scoped(&all.marks), 0);
    assert_eq!(scoped(&area.marks), 0);
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
    // The represented area holds the definition's and is wider by what the
    // carriers and the rail strips can move it.
    let (stated, exact) = (&all.value.area, &definition.value.area);
    assert!(
        stated.lo <= exact.lo && exact.hi <= stated.hi,
        "{stated:?} {exact:?}"
    );
    assert!(stated.width() > exact.width());
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
