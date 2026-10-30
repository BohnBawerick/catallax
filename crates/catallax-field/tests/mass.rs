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

#[test]
fn exchanged_box_bounds_cover_volume_centroid_and_inertia_perturbations() {
    use catallax_features::stand_in::{MassTarget, equivalent_blocks};
    let target = MassTarget {
        volume: q("6"),
        centroid: [q("0.5"), q("1"), q("1.5")],
        inertia: std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if i == j {
                    [q("6.5"), q("5"), q("2.5")][i].clone()
                } else {
                    q("0")
                }
            })
        }),
    };
    let proxy = equivalent_blocks(&mut catallax_topology::Model::new(), &target).unwrap();
    let bound = mass::surrogate_properties(&proxy, &q("0.001")).unwrap();
    // Enlarge by 0.00025 on each side and translate x by 0.00025.
    // Every corresponding corner moves less than the 0.001 deviation.
    let sides = [q("1.0005"), q("2.0005"), q("3.0005")];
    let volume = &(&sides[0] * &sides[1]) * &sides[2];
    assert!(bound.volume.contains(&volume));
    assert!(bound.centroid[0].contains(&q("0.50025")));
    for (value, expected) in bound.centroid.iter().zip(&target.centroid).skip(1) {
        assert!(value.contains(expected));
    }
    for (i, row) in bound.inertia.iter().enumerate() {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let inertia = (&volume * &(&(&sides[j] * &sides[j]) + &(&sides[k] * &sides[k])))
            .checked_div(&q("12"))
            .unwrap();
        assert!(row[i].contains(&inertia));
        for (j, value) in row.iter().enumerate() {
            if i != j {
                assert!(value.contains(&q("0")));
            }
        }
    }
}

#[test]
fn subtracting_lattice_cannot_use_its_region_as_an_outer_envelope() {
    let recipe=Recipe::parse(r#"{"format":"catallax-recipe/3","name":"negative lattice","steps":[
        {"id":"outer","op":"block","size":["8","8","8"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
        {"id":"region","op":"block","size":["4","4","4"],"pose":{"rotation":["1","0","0","0"],"translation":["2","2","2"]}},
        {"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"radius":{"constant":"0.25"}},
        {"id":"part","op":"boolean","kind":"subtract","target":"outer","tool":"beams"}],"result":"part"}"#).unwrap();
    let plan = Plan::compile(&recipe).unwrap();
    // Material between the beams survives. Subtracting the region instead
    // would erase it and therefore could not be an outer envelope.
    let point = catallax_shapes::RPoint3::new(q("3"), q("3"), q("3"));
    assert_eq!(
        plan.proven(&point, &mut catallax_field::Workspace::default())
            .location,
        catallax_field::Location::Inside
    );
    assert_eq!(
        mass::check_region_envelope(&plan).unwrap_err().code,
        "field.stand-in-envelope-coverage"
    );
}
#[test]
fn a_gyroid_has_no_mass_route_in_this_increment() {
    let recipe = Recipe::parse(r#"{"format":"catallax-recipe/3","name":"gyroid","steps":[{"id":"box","op":"block","size":["4","4","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"net","op":"gyroid","region":"box","period":"4","level":{"constant":"0"},"kind":"network","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"net"}"#).unwrap();
    let error = mass::properties(
        &Plan::compile(&recipe).unwrap(),
        &Accuracy::default(),
        &Options::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "field.print-coverage");
}
