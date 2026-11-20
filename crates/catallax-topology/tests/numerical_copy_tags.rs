//! Copy tagging (decision Q8, section 4.4 of the Checkpoint 2 decisions):
//! what admission does with a tag, what resets one, and what the audit sees.
#[path = "common/analytic.rs"]
mod fixture;
use catallax_numbers::{Decider, NumericalCheck, NumericalSample, Rational};
use catallax_shapes::{Curve, Cylinder, Surface};
use catallax_topology::numerical::{self as n, CopyRule, CopySource, CopyTag};
use fixture::{q, quarter};

/// The quarter cylinder's curved face.
const WALL: usize = 4;
/// Its flat bottom.
const BOTTOM: usize = 0;
/// The record of a measured use, and the record of a derived one.
const MEASURED: &str = "trim.edge-surface-distance";
const BOUND: &str = "trim.edge-surface-derived-bound";

fn derived(bound: &str) -> CopyTag {
    CopyTag::Derived {
        from: CopySource::Edge3d,
        rule: CopyRule::ExactComposition,
        bound: q(bound),
    }
}

/// Tags every circle use of one face.
fn tag_circles(source: &mut n::NumericalDefinition, face: usize, tag: &CopyTag) {
    let circles: Vec<bool> = source
        .edges
        .iter()
        .map(|edge| matches!(edge.curve, Curve::Circle(_)))
        .collect();
    for trim in source.faces[face].loops.iter_mut().flatten() {
        if circles[trim.edge] {
            trim.tag = tag.clone();
        }
    }
}

/// The edge-surface records of one face's circle uses.
fn circle_records(body: &n::NumericalBody, face: usize, record: &str) -> Vec<NumericalSample> {
    let definition = body.definition();
    let face = &definition.faces[face];
    body.evidence()
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples {
                name,
                scope,
                samples,
            } if name == record && scope.starts_with(&format!("face {:?},", face.key)) => {
                Some((scope, samples))
            }
            _ => None,
        })
        .filter(|(scope, _)| {
            definition.edges.iter().any(|edge| {
                matches!(edge.curve, Curve::Circle(_))
                    && scope.contains(&format!("edge {:?},", edge.key))
            })
        })
        .flat_map(|(_, samples)| samples.clone())
        .collect()
}

fn build(source: n::NumericalDefinition) -> Result<n::NumericalBody, catallax_numbers::Refusal> {
    n::build(
        source,
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
}

fn radius(source: &mut n::NumericalDefinition, value: &str) {
    let Surface::Cylinder(cylinder) = &mut source.faces[WALL].surface else {
        unreachable!()
    };
    cylinder.radius_sq = (&q(value) * &q(value)).into();
}

#[test]
fn a_derived_copy_stores_its_bound_and_an_independent_one_is_measured() {
    // Untagged, the two circles on the wall are measured: each record has
    // the parameter of its largest gap.
    let measured = build(quarter()).unwrap();
    assert_eq!(circle_records(&measured, WALL, MEASURED).len(), 2);
    assert!(circle_records(&measured, WALL, BOUND).is_empty());
    assert!(
        n::audit_derived_copies(&measured, &mut Decider::default())
            .unwrap()
            .is_empty()
    );

    // Tagged, the record is the stated bound, with no margin and no place.
    let mut source = quarter();
    tag_circles(&mut source, WALL, &derived("1/1000"));
    let tagged = build(source).unwrap();
    assert!(circle_records(&tagged, WALL, MEASURED).is_empty());
    let records = circle_records(&tagged, WALL, BOUND);
    assert_eq!(records.len(), 2);
    for sample in &records {
        // The two ends of the use, and the bound itself.
        assert_eq!(sample.parameters.len(), 2);
        assert_eq!(sample.distance_upper_mm, q("1/1000"));
    }
    // The audit runs the full check on both and finds them inside the bound.
    let audit = n::audit_derived_copies(&tagged, &mut Decider::default()).unwrap();
    assert_eq!(audit.len(), 2);
    for row in &audit {
        assert_eq!(row.rule, CopyRule::ExactComposition);
        assert!(row.holds(), "{row:?}");
    }
}

#[test]
fn a_bound_over_the_limit_refuses_and_a_negative_bound_is_not_a_definition() {
    let mut source = quarter();
    tag_circles(&mut source, WALL, &derived("1/50"));
    assert_eq!(
        build(source).unwrap_err().code,
        "numbers.numerical-tolerance-exceeded"
    );
    let mut source = quarter();
    tag_circles(&mut source, WALL, &derived("-1/1000"));
    assert_eq!(
        build(source).unwrap_err().code,
        "numerical.invalid-definition"
    );
}

/// The risk the audit exists for: a wrong tag ships an unchecked gap. The
/// wall is 0.05 mm too wide, which admission refuses when it measures. With
/// its circles tagged as exact the body is admitted, and the audit is what
/// reports the copies as outside their bound.
#[test]
fn the_audit_catches_a_wrong_tag_that_admission_cannot() {
    let mut source = quarter();
    radius(&mut source, "1.05");
    assert_eq!(
        build(source.clone()).unwrap_err().code,
        "numbers.numerical-tolerance-exceeded"
    );
    tag_circles(&mut source, WALL, &derived("0"));
    // The lines on the wall still refuse when measured, so move only what
    // the tag covers out of admission's sight: check the audit on the uses.
    let audit = match build(source.clone()) {
        Ok(body) => n::audit_derived_copies(&body, &mut Decider::default()).unwrap(),
        Err(_) => {
            // Admission refused on an untagged use of the wall. Tag those as
            // well, wrongly, so the body is admitted.
            for trim in source.faces[WALL].loops.iter_mut().flatten() {
                trim.tag = derived("0");
            }
            let body = build(source).unwrap();
            n::audit_derived_copies(&body, &mut Decider::default()).unwrap()
        }
    };
    assert!(!audit.is_empty());
    let wrong: Vec<_> = audit.iter().filter(|row| !row.holds()).collect();
    assert!(!wrong.is_empty());
    for row in wrong {
        assert!(row.measured_mm > q("0.04"), "{row:?}");
    }
}

/// A circle on a plane is read against the plane whatever its tag says: that
/// reading is the only check that the 3D curve lies on the plane.
#[test]
fn a_plane_projected_use_is_measured_whatever_its_tag() {
    let mut source = quarter();
    tag_circles(&mut source, BOTTOM, &derived("0"));
    let body = build(source.clone()).unwrap();
    assert_eq!(circle_records(&body, BOTTOM, MEASURED).len(), 1);
    assert!(circle_records(&body, BOTTOM, BOUND).is_empty());
    assert!(
        n::audit_derived_copies(&body, &mut Decider::default())
            .unwrap()
            .is_empty()
    );
    // Lifted 0.05 mm off its plane, the tagged circle still refuses.
    for edge in &mut source.edges {
        let mut ends = edge.vertices;
        ends.sort();
        if let (Curve::Circle(circle), [1, 2]) = (&mut edge.curve, ends) {
            circle.center.z = q("1/20").into();
        }
    }
    // The same circle is also a use of the wall; tag that one as well so the
    // plane is the only face that can object.
    tag_circles(&mut source, WALL, &derived("0"));
    let refusal = build(source).unwrap_err();
    assert_eq!(refusal.code, "numbers.numerical-tolerance-exceeded");
}

/// Replacing a face's surface changes one copy of every use on it, so their
/// tags fall back to independent and the uses are measured again.
#[test]
fn replacing_a_surface_resets_the_tags_of_its_uses() {
    let mut source = quarter();
    tag_circles(&mut source, WALL, &derived("0"));
    let body = build(source).unwrap();
    assert_eq!(circle_records(&body, WALL, BOUND).len(), 2);
    let Surface::Cylinder(cylinder) = &body.definition().faces[WALL].surface else {
        unreachable!()
    };
    let wider = Surface::Cylinder(Cylinder {
        radius_sq: (&q("1.001") * &q("1.001")).into(),
        ..cylinder.clone()
    });
    let edited = body
        .replace_surface("f4", wider, &mut Decider::default())
        .unwrap();
    let wall = &edited.definition().faces[WALL];
    assert!(
        wall.loops
            .iter()
            .flatten()
            .all(|trim| trim.tag == CopyTag::Independent)
    );
    assert!(circle_records(&edited, WALL, BOUND).is_empty());
    let records = circle_records(&edited, WALL, MEASURED);
    assert_eq!(records.len(), 2);
    for sample in &records {
        // 1.1 times the 0.001 mm the wall moved, and its measurement error.
        let excess = &sample.distance_upper_mm - &q("0.0011");
        assert!(
            excess >= Rational::zero() && excess < q("1e-12"),
            "{excess}"
        );
    }
    // The body the edit started from keeps its tags.
    assert!(
        body.definition().faces[WALL]
            .loops
            .iter()
            .flatten()
            .any(|trim| trim.tag != CopyTag::Independent)
    );
}

/// An edit of one copy that leaves the tag as it was. The edited definition is
/// a clone of a tagged body's, so every use still states its derived tag; only
/// the 3D curve of the top circle changes. The rebuild cannot show that the
/// stated bound still holds for that use, so it is measured again, and the
/// untouched circle keeps its tag. Rebuilding the unedited definition keeps
/// both.
#[test]
fn a_rebuild_measures_a_use_whose_3d_copy_changed_under_an_unchanged_tag() {
    let mut source = quarter();
    tag_circles(&mut source, WALL, &derived("0"));
    let body = build(source).unwrap();
    assert_eq!(circle_records(&body, WALL, BOUND).len(), 2);

    let mut dec = Decider::default();
    let same = body
        .rebuild(body.definition().clone(), body.options().clone(), &mut dec)
        .unwrap();
    assert_eq!(circle_records(&same, WALL, BOUND).len(), 2);
    assert!(circle_records(&same, WALL, MEASURED).is_empty());

    // Widen the top circle (vertices 4 and 5) by 0.001 mm; the wall, the u/v
    // ends and every tag stay as they were.
    let mut edited = body.definition().clone();
    for edge in &mut edited.edges {
        let mut ends = edge.vertices;
        ends.sort();
        if let (Curve::Circle(circle), [4, 5]) = (&mut edge.curve, ends) {
            circle.radius_sq = (&q("1.001") * &q("1.001")).into();
        }
    }
    let rebuilt = body
        .rebuild(edited, body.options().clone(), &mut dec)
        .unwrap();
    assert_eq!(circle_records(&rebuilt, WALL, BOUND).len(), 1);
    let records = circle_records(&rebuilt, WALL, MEASURED);
    assert_eq!(records.len(), 1);
    // 1.1 times the 0.001 mm the circle moved, and its measurement error.
    let excess = &records[0].distance_upper_mm - &q("0.0011");
    assert!(
        excess >= Rational::zero() && excess < q("1e-12"),
        "{excess}"
    );
    let wall = &rebuilt.definition().faces[WALL];
    assert_eq!(
        wall.loops
            .iter()
            .flatten()
            .filter(|trim| trim.tag != CopyTag::Independent)
            .count(),
        1
    );
}
