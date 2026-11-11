//! Same answers as the print certificate's classifier wherever it has coverage.
use super::{BoxClassifier, BoxStorage, Classifier, bounds};
use crate::Plan;
use catallax_numbers::{Decider, Interval, Rational};
use catallax_recipe::Recipe;

fn part(name: &str) -> Plan {
    let path = format!("{}/../../parts/{name}.json", env!("CARGO_MANIFEST_DIR"));
    let recipe = Recipe::parse(&std::fs::read_to_string(path).unwrap()).unwrap();
    Plan::compile_checked(&recipe, &mut Decider::default()).unwrap()
}

#[test]
fn full_coverage_keeps_the_print_answers_where_print_has_coverage() {
    for name in ["P1-bracket", "P4-shelled-tray", "P6-lattice-bracket"] {
        let plan = part(name);
        let mut dec = Decider::default();
        let [x, y, z] = bounds(&plan, &mut dec).unwrap();
        let mut print = Classifier::for_print(&plan, &mut dec).unwrap();
        let full = BoxClassifier::new(&plan, &mut dec).unwrap();
        let mut storage = BoxStorage::default();
        let cells = 13;
        let at = |b: &[Rational; 2], k: i64| {
            &b[0] + &((&b[1] - &b[0]) * Rational::from_ratio(k, cells).unwrap())
        };
        let mut compared = 0;
        for i in 0..cells {
            for j in 0..cells {
                for k in 0..cells {
                    let lo = [at(&x, i), at(&y, j), at(&z, k)];
                    let hi = [at(&x, i + 1), at(&y, j + 1), at(&z, k + 1)];
                    let cover = std::array::from_fn(|a| {
                        Interval::new(
                            Interval::from_rational(&lo[a]).lo(),
                            Interval::from_rational(&hi[a]).hi(),
                        )
                        .unwrap()
                    });
                    assert_eq!(
                        full.classify(&lo, &hi, &mut storage, &mut dec).unwrap(),
                        print.classify(&cover, &mut dec).unwrap(),
                        "{name} box {lo:?} {hi:?}"
                    );
                    compared += 1;
                }
            }
        }
        assert_eq!(compared, cells.pow(3));
    }
}
