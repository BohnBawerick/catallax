//! Field moments are checked against analytic solids, never a print mesh.
use catallax_field::{
    Plan,
    mass::{self, Accuracy, Options},
};
use catallax_measure::Accuracy as QuantityAccuracy;
use catallax_numbers::{Decider, Rational};
use catallax_recipe::Recipe;
fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
#[test]
fn translated_block_moments_are_exact_at_default_accuracy() {
    let recipe=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"moments","steps":[{"id":"box","op":"block","size":["1","2","3"],"pose":{"rotation":["1","0","0","0"],"translation":["5","6","7"]}}],"result":"box"}"#).unwrap();
    let result = mass::properties(
        &Plan::compile(&recipe).unwrap(),
        &Accuracy::default(),
        &Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(result.volume.lo, q("6"));
    assert!(result.volume.is_exact());
    let centroid = result.centroid.as_ref().unwrap();
    for (value, expected) in centroid.iter().zip(["5.5", "7", "8.5"]) {
        assert!(value.is_exact());
        assert_eq!(value.lo, q(expected));
    }
    let inertia = result.inertia.as_ref().unwrap();
    for (i, expected) in ["6.5", "5", "2.5"].into_iter().enumerate() {
        for (j, value) in inertia[i].iter().enumerate() {
            assert!(value.is_exact());
            assert_eq!(value.lo, if i == j { q(expected) } else { q("0") });
        }
    }
}
#[test]
fn curved_cavity_bounds_contain_independent_analytic_mass_and_inertia() {
    let recipe=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"sphere cavity","steps":[
        {"id":"box","op":"block","size":["2","2","2"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"cavity","op":"sphere","radius":"0.5","pose":{"rotation":["1","0","0","0"],"translation":["1","1","1"]}},
        {"id":"part","op":"boolean","kind":"subtract","target":"box","tool":"cavity"}],"result":"part"}"#).unwrap();
    let requested = QuantityAccuracy {
        relative: q("0"),
        absolute: q("0.25"),
    };
    let accuracy = Accuracy {
        volume: requested.clone(),
        centroid: requested.clone(),
        inertia: requested,
    };
    let result = mass::properties(
        &Plan::compile(&recipe).unwrap(),
        &accuracy,
        &Options::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let pi_lo = q("3.1415926535897932384626433832795028");
    let pi_hi = q("3.1415926535897932384626433832795029");
    let volume_lo = &q("8") - &pi_hi.checked_div(&q("6")).unwrap();
    let volume_hi = &q("8") - &pi_lo.checked_div(&q("6")).unwrap();
    assert!(result.volume.lo <= volume_lo && result.volume.hi >= volume_hi);
    assert!(result.volume.width() <= q("0.25"));
    for c in result.centroid.as_ref().unwrap() {
        assert!(c.contains(&q("1")));
        assert!(c.width() <= q("0.25"));
    }
    let inertia = result.inertia.as_ref().unwrap();
    let lo = &q("16/3") - &pi_hi.checked_div(&q("60")).unwrap();
    let hi = &q("16/3") - &pi_lo.checked_div(&q("60")).unwrap();
    for (i, row) in inertia.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            if i == j {
                assert!(value.lo <= lo && value.hi >= hi);
            } else {
                assert!(value.contains(&q("0")));
            }
            assert!(value.width() <= q("0.25"));
        }
    }
}
