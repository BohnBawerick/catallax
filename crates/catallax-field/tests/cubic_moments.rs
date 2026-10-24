//! Independent analytic moments of a cubic grid cell, including node overlaps.
use catallax_field::{
    Plan,
    mass::{self, Accuracy, Options},
};
use catallax_measure::Accuracy as QuantityAccuracy;
use catallax_numbers::{Ball, Decider, Rational};
use catallax_recipe::Recipe;
use serde_json::json;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

#[test]
fn cubic_cell_moments_include_pair_and_triple_intersections_after_translation() {
    for (size, radius, translation, subtract) in [
        ("2", "0.25", ["-3", "5", "11"], false),
        ("3", "0.3", ["0", "0", "0"], true),
    ] {
        let a = q(size);
        let r = q(radius);
        let center = translation.map(|v| &q(v) + &a.mul_pow2(-1));
        let mut recipe = json!({"format":"catallax-recipe/3","name":"analytic cubic cell","steps":[
            {"id":"box","op":"block","size":[size,size,size],"pose":{"rotation":["1","0","0","0"],"translation":translation}},
            {"id":"grid","op":"lattice","region":"box","cell":"cubic","size":size,
             "frame":{"rotation":["0","1","0","0"],"translation":center.each_ref().map(ToString::to_string)},
             "radius":{"constant":radius}}
        ],"result":"grid"});
        if subtract {
            recipe["steps"].as_array_mut().unwrap().push(json!({
                "id":"voids","op":"boolean","kind":"subtract","target":"box","tool":"grid"
            }));
            recipe["result"] = json!("voids");
        }
        let requested = QuantityAccuracy {
            relative: q("0"),
            absolute: q("1e-9"),
        };
        let result = mass::properties(
            &Plan::compile(&Recipe::parse(&recipe.to_string()).unwrap()).unwrap(),
            &Accuracy {
                volume: requested.clone(),
                centroid: requested.clone(),
                inertia: requested,
            },
            &Options::default(),
            &mut Decider::default(),
        )
        .unwrap();
        let b = |v: &Rational| Ball::from_rational(v, 128);
        let pi = Ball::pi(128);
        let root = b(&q("2")).sqrt().unwrap();
        // Three cylinders, less three bicylinder intersections, plus the
        // common tricylinder. These are moments of the curved field itself.
        let volume = &(&(&b(&q("3")) * &b(&a)) * &(&pi * &b(&(&r * &r))))
            - &(&(&b(&q("8")) * &root) * &b(&r.pow(3).unwrap()));
        let axial = (&pi * &b(&(&(&r * &r) * &a.pow(3).unwrap())))
            .checked_div(&b(&q("6")))
            .unwrap();
        let transverse = &pi * &b(&(&a * &r.pow(4).unwrap()));
        let overlap = &(&b(&q("184")) * &root) * &b(&r.pow(5).unwrap());
        let inertia = &(&axial + &transverse) - &overlap.checked_div(&b(&q("45"))).unwrap();
        let (volume, inertia) = if subtract {
            (
                &b(&a.pow(3).unwrap()) - &volume,
                &b(&a.pow(5).unwrap()).checked_div(&b(&q("6"))).unwrap() - &inertia,
            )
        } else {
            (volume, inertia)
        };
        assert!(result.volume.lo <= volume.lo() && result.volume.hi >= volume.hi());
        assert!(result.volume.width() <= q("1e-9"));
        for (value, expected) in result.centroid.as_ref().unwrap().iter().zip(&center) {
            assert!(value.contains(expected));
            assert!(value.width() <= q("1e-9"));
        }
        for (i, row) in result.inertia.as_ref().unwrap().iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                if i == j {
                    assert!(value.lo <= inertia.lo() && value.hi >= inertia.hi());
                } else {
                    assert!(value.contains(&q("0")));
                }
                assert!(value.width() <= q("1e-9"));
            }
        }
    }
}
