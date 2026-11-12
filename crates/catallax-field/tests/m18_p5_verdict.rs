//! Exact field verdict for the authored P5 radius-2 fillet.
//!
//! The untrimmed torus plane rail crosses each countersink mouth. The numbers
//! below are derived from the public part and from the compiler's own admitted
//! radius-1/2 torus, and the exact rim-contact envelope is constructed at the
//! authored radius. No recipe answer, tolerance or allowance changes here.
use catallax_field::Plan;
use catallax_field::rolling_region::Region;
use catallax_numbers::{Decider, Outcome, Rational, StrictSign};
use catallax_recipe::{Recipe, Role};
use catallax_shapes::{Point3, RollingEnvelopeDefinition, Surface, Vector3};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn exact(value: &catallax_numbers::Real) -> Rational {
    value.as_rational().unwrap()
}
fn zero(value: &catallax_numbers::Real, dec: &mut Decider) {
    assert!(
        matches!(dec.sign(value), Outcome::ProvenCoincident(_)),
        "{value:?}"
    );
}
fn authored(radius: &str) -> Recipe {
    let mut value: serde_json::Value =
        serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    let step = value["steps"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|step| step["id"] == "filleted")
        .unwrap();
    step["radius"] = radius.into();
    Recipe::parse(&value.to_string()).unwrap()
}
fn step<'a>(value: &'a serde_json::Value, id: &str) -> &'a serde_json::Value {
    value["steps"]
        .as_array()
        .unwrap()
        .iter()
        .find(|step| step["id"] == id)
        .unwrap()
}
fn number(value: &serde_json::Value) -> Rational {
    Rational::parse(value.as_str().unwrap()).unwrap()
}

#[test]
fn the_plain_torus_rail_crosses_each_mouth_by_one_millimetre() {
    let value: serde_json::Value =
        serde_json::from_str(include_str!("../../../parts/P5-flanged-housing.json")).unwrap();
    let sink = step(&value, "sink");
    let mouth = number(&sink["base_radius"]);
    let offset = number(&sink["pose"]["translation"][0]);
    let mouth_plane = number(&sink["pose"]["translation"][1]);
    let to = step(&value, "section")["loops"][0]
        .as_array()
        .unwrap()
        .iter()
        .find(|segment| segment["id"] == "flange-top")
        .unwrap()["to"]
        .as_array()
        .unwrap()
        .clone();
    let body = number(&to[0]);
    assert_eq!(mouth_plane, number(&to[1]));
    assert_eq!(
        step(&value, "hole-pattern")["layout"]["circular"]["count"].as_u64(),
        Some(4)
    );
    let radius = number(&step(&value, "filleted")["radius"]);
    assert_eq!(
        [&body, &mouth, &offset, &radius],
        [&q("12"), &q("4"), &q("17"), &q("2")]
    );

    // The admitted radius-1/2 build fixes the rail from the compiler's torus.
    let mut dec = Decider::default();
    let plan = Plan::compile_checked(&authored("1/2"), &mut dec).unwrap();
    let torus = plan
        .faces()
        .iter()
        .find_map(|face| match &face.surface {
            Surface::Torus(torus) if matches!(face.origin.role, Role::Blend { .. }) => Some(torus),
            _ => None,
        })
        .expect("the admitted fillet is one torus");
    assert_eq!(exact(&(&torus.major - &torus.minor)), body);
    assert_eq!(exact(&torus.minor), q("1/2"));

    let rail = &body + &radius;
    let inner = &offset - &mouth;
    assert_eq!(rail, q("14"));
    assert_eq!(inner, q("13"));
    // The rail lies exactly 1 mm inside the mouth's nearest rim.
    assert_eq!(&rail - &inner, Rational::one());
    // Two transversal crossings per mouth: |rail - mouth| < offset < rail + mouth.
    assert!(&rail - &mouth < offset && offset < &rail + &mouth);

    let two = q("2");
    let x = (&offset * &offset + &rail * &rail - &mouth * &mouth) / (&offset * &two);
    assert_eq!(x, q("469/34"));
    let z_squared = &rail * &rail - &x * &x;
    assert_eq!(z_squared, q("6615/1156"));
    assert_eq!((&x - &offset) * (&x - &offset) + z_squared, &mouth * &mouth);
    // The crossing angle is the rolling region's exact transition angle.
    assert_eq!(&x / &rail, q("67/68"));
    assert_eq!(
        &x / &rail,
        (&rail * &rail + &offset * &offset - &mouth * &mouth) / (&rail * &offset * two)
    );
}

#[test]
fn the_rim_contact_continuation_is_exact_at_radius_two() {
    let mut dec = Decider::default();
    let region = Region::new(
        RollingEnvelopeDefinition::new(
            Point3::new(q("0").into(), q("5").into(), q("0").into()),
            Vector3::new(0.into(), 1.into(), 0.into()),
            Vector3::new(1.into(), 0.into(), 0.into()),
            q("12").into(),
            q("4").into(),
            q("17").into(),
            q("2").into(),
            [q("1").into(), q("1").into()],
            &mut dec,
        )
        .unwrap(),
        &mut dec,
    )
    .unwrap();

    // Deepest station: the ball centre sits 2 - sqrt(3) mm below the plane
    // rail's height and touches the mouth's nearest rim point (13, 5, 0).
    let middle = region.station(&q("0"), &mut dec).unwrap();
    let real = |value: &Rational| catallax_numbers::Real::from(value.clone());
    let five = real(&q("5"));
    let dip = &middle.centre.y - &five;
    zero(&(dip.square() - real(&q("3"))), &mut dec);
    assert_eq!(dec.sign(&dip), Outcome::Proven(StrictSign::Positive));
    let two = real(&q("2"));
    assert_eq!(
        dec.sign(&(&two - &dip)),
        Outcome::Proven(StrictSign::Positive)
    );
    zero(&(&middle.centre.x - &real(&q("14"))), &mut dec);
    zero(&(&middle.mouth.x - &real(&q("13"))), &mut dec);
    zero(&(&middle.mouth.y - &five), &mut dec);
    zero(&middle.mouth.z, &mut dec);

    // Both transitions touch the plane exactly where the rail crosses the
    // mouth: the two points (469/34, 5, +-21*sqrt(15)/34).
    let sqrt3 = real(&q("3")).sqrt();
    zero(&(&middle.centre.y - &(&five + &sqrt3)), &mut dec);
    for (u, sign) in [("-1", StrictSign::Negative), ("1", StrictSign::Positive)] {
        let corner = region.station(&q(u), &mut dec).unwrap();
        zero(&(&corner.mouth.x - &real(&q("469/34"))), &mut dec);
        zero(&(&corner.mouth.y - &five), &mut dec);
        zero(&(corner.mouth.z.square() - real(&q("6615/1156"))), &mut dec);
        assert_eq!(dec.sign(&corner.mouth.z), Outcome::Proven(sign));
    }
}

#[test]
fn only_the_strictly_clear_side_builds() {
    for (radius, builds) in [("1/2", true), ("3/4", true), ("1", false), ("2", false)] {
        let mut dec = Decider::default();
        match Plan::compile_checked(&authored(radius), &mut dec) {
            Ok(_) => assert!(builds, "radius {radius} unexpectedly built"),
            Err(refusal) => {
                assert!(!builds, "radius {radius} unexpectedly refused");
                assert_eq!(refusal.code, "fillet.unsupported-overflow", "{radius}");
            }
        }
    }
}
