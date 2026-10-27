//! Independent section points for the immutable P3 recipe and its radius edit.

use catallax_boolean::{Containment, classify_point};
use catallax_features::evaluate;
use catallax_numbers::{BudgetPolicy, Decider, Outcome, Rational};
use catallax_recipe::{Operation, Recipe};
use catallax_shapes::RPoint3;

fn decider() -> Decider {
    Decider::new(BudgetPolicy {
        max_work: 1_000_000,
        ..Default::default()
    })
}

#[test]
fn p3_original_and_edited_material_matches_independent_sections() {
    #[derive(Clone, Copy, Debug)]
    enum Expected {
        Inside,
        Outside,
        Boundary,
    }
    use Expected::{Boundary, Inside, Outside};

    let mut recipe = Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap();
    for radius in ["2", "5/2"] {
        let Operation::Fillet { radius: value, .. } = &mut recipe.steps[3].op else {
            panic!("P3 junction must remain a fillet");
        };
        *value = serde_json::from_value(serde_json::json!(radius)).unwrap();
        let mut built = evaluate(&recipe, &mut decider()).unwrap();
        let mut probes = vec![
            // Both disputed grid points lie outside the finite stock and boss.
            (["40", "-10", "88/3"], Outside),
            (["80/3", "0", "22"], Outside),
            (["39999/1000", "-10", "88/3"], Outside),
            (["40001/1000", "-10", "88/3"], Outside),
            (["80/3", "-1/1000", "22"], Outside),
            (["80/3", "1/1000", "22"], Outside),
            (["80/3", "0", "21999/1000"], Outside),
            (["80/3", "0", "22001/1000"], Outside),
            // The four top rounds have radius 1. At z=49/5, the
            // inset is 1-sqrt(1-(4/5)^2)=2/5 on both horizontal axes.
            (["3/10", "15", "49/5"], Outside),
            (["2/5", "15", "49/5"], Boundary),
            (["1/2", "15", "49/5"], Inside),
            (["20", "3/10", "49/5"], Outside),
            (["20", "2/5", "49/5"], Boundary),
            (["20", "1/2", "49/5"], Inside),
            (["397/10", "15", "49/5"], Outside),
            (["198/5", "15", "49/5"], Boundary),
            (["79/2", "15", "49/5"], Inside),
            (["20", "297/10", "49/5"], Outside),
            (["20", "148/5", "49/5"], Boundary),
            (["20", "59/2", "49/5"], Inside),
            (["1/2", "1/2", "49/5"], Inside),
            (["3/10", "1/2", "49/5"], Outside),
            // Boss side and cap remain at rho=8 and z=22 in both variants.
            (["28", "15", "20"], Boundary),
            (["27999/1000", "15", "20"], Inside),
            (["28001/1000", "15", "20"], Outside),
            (["20", "15", "22"], Boundary),
            (["20", "15", "21999/1000"], Inside),
            (["20", "15", "22001/1000"], Outside),
            (["20", "15", "0"], Boundary),
            (["20", "15", "-1/1000"], Outside),
            (["20", "15", "1/1000"], Inside),
        ];
        if radius == "2" {
            // Junction centre is (rho,z)=(10,12). At z=54/5,
            // rho=10-sqrt(4-(6/5)^2)=42/5.
            probes.extend([
                (["283/10", "15", "54/5"], Inside),
                (["142/5", "15", "54/5"], Boundary),
                (["57/2", "15", "54/5"], Outside),
                (["284/10", "15", "11"], Outside),
            ]);
        } else {
            // Junction centre is (rho,z)=(21/2,25/2). At z=11,
            // rho=21/2-sqrt(25/4-9/4)=17/2.
            probes.extend([
                (["284/10", "15", "11"], Inside),
                (["57/2", "15", "11"], Boundary),
                (["286/10", "15", "11"], Outside),
            ]);
        }
        for (coordinates, expected) in probes {
            let [x, y, z] = coordinates.map(|v| Rational::parse(v).unwrap());
            let point = RPoint3::new(x, y, z);
            let actual = classify_point(&mut built.model, built.result, &point, &mut decider());
            assert!(
                matches!(
                    (&actual, expected),
                    (Outcome::Proven(Containment::Inside), Inside)
                        | (Outcome::Proven(Containment::Outside), Outside)
                        | (Outcome::ProvenCoincident(_), Boundary)
                ),
                "radius={radius}, point={coordinates:?}: expected {expected:?}, got {actual:?}"
            );
        }
    }
}
