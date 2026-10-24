//! Independent solid decomposition checks for finite normal-sector integrals.

use catallax_features::evaluate;
use catallax_measure::{Accuracy, Enclosure, mass_properties};
use catallax_numbers::{Ball, BudgetPolicy, Decider, Rational};
use catallax_recipe::Recipe;
use catallax_shapes::Rotation;
use serde_json::json;

const BITS: u32 = 512;
fn b(n: i64) -> Ball {
    Ball::from_rational(&Rational::from_i64(n), BITS)
}
fn contains(actual: &Enclosure, expected: &Ball) {
    assert!(
        actual.lo <= expected.lo() && actual.hi >= expected.hi(),
        "{actual:?} does not contain {expected:?}"
    );
    let magnitude = actual.lo.abs().min(actual.hi.abs());
    let accuracy = Accuracy::default();
    assert!(actual.width() <= (accuracy.relative * magnitude).max(accuracy.absolute));
}

// Disjoint core, face slabs, quarter-cylinder edges and octant-ball corners.
// This volume decomposition does not use the surface divergence integrals.
fn rounded_box(size: [i64; 3], radius: i64) -> (Ball, Ball, [Ball; 3]) {
    let pi = Ball::pi(BITS);
    let [a, bb, c] = size.map(b);
    let r = b(radius);
    let r2 = &r * &r;
    let r3 = &r2 * &r;
    let r4 = &r3 * &r;
    let r5 = &r4 * &r;
    let core = &(&a * &bb) * &c;
    let surface = &b(2) * &(&(&a * &bb) + &(&(&bb * &c) + &(&c * &a)));
    let length = &(&a + &bb) + &c;
    let volume = &(&core + &(&surface * &r))
        + &(&(&pi * &r2) * &length)
        + (&(&b(4) * &pi) * &r3).checked_div(&b(3)).unwrap();
    let area = &surface + &(&(&b(2) * &pi) * &(&length * &r)) + (&(&b(4) * &pi) * &r2);
    let lengths = [a, bb, c];
    let second = std::array::from_fn::<_, 3, _>(|i| {
        let x = &lengths[i];
        let y = &lengths[(i + 1) % 3];
        let z = &lengths[(i + 2) % 3];
        let x2 = x * x;
        let x3 = &x2 * x;
        let core_moment = (&core * &x2).checked_div(&b(12)).unwrap();
        let normal_slabs = &(y * z)
            * &(&(&x2 * &r).checked_div(&b(2)).unwrap()
                + &(x * &r2)
                + (&b(2) * &r3).checked_div(&b(3)).unwrap());
        let other_slabs = (&(&r * &x3) * &(y + z)).checked_div(&b(6)).unwrap();
        let axial_edges = (&(&pi * &r2) * &x3).checked_div(&b(12)).unwrap();
        let radial_edges = &(y + z)
            * &(&(&(&pi * &r2) * &x2).checked_div(&b(4)).unwrap()
                + &(&(&b(4) * x) * &r3).checked_div(&b(3)).unwrap()
                + (&pi * &r4).checked_div(&b(4)).unwrap());
        let corners = &(&(&pi * &r3) * &x2).checked_div(&b(3)).unwrap()
            + &(&(&pi * &r4) * x).checked_div(&b(2)).unwrap()
            + (&(&b(4) * &pi) * &r5).checked_div(&b(15)).unwrap();
        &(&(&core_moment + &normal_slabs) + &other_slabs)
            + &(&axial_edges + &(&radial_edges + &corners))
    });
    let inertia = std::array::from_fn(|i| &second[(i + 1) % 3] + &second[(i + 2) % 3]);
    (volume, area, inertia)
}

#[test]
fn finite_offset_sectors_prove_translated_rotated_mass_and_tensor() {
    let size = [4, 7, 9];
    let quaternion = [2, 1, 3, 4];
    let translation = [-7, 11, 13];
    let radius = 2;
    let recipe = Recipe::parse(
        &json!({
            "format":"catallax-recipe/3", "name":"posed offset mass",
            "steps":[
                {"id":"stock","op":"block","size":size.map(|v|v.to_string()),
                 "pose":{"rotation":quaternion.map(|v|v.to_string()),
                         "translation":translation.map(|v|v.to_string())}},
                {"id":"rounded","op":"offset","body":"stock","distance":radius.to_string()}
            ], "result":"rounded"
        })
        .to_string(),
    )
    .unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let answer = mass_properties(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    let (volume, area, principal) = rounded_box(size, radius);
    let mass = answer.value;
    contains(&mass.volume, &volume);
    contains(&mass.area, &area);
    let [w, x, y, z] = quaternion.map(Rational::from_i64);
    let rotation = Rotation::from_quaternion(w, x, y, z).unwrap();
    let matrix = rotation.rational_matrix().unwrap();
    let center = mass.centroid.unwrap();
    let tensor = mass.inertia.unwrap();
    for i in 0..3 {
        let expected = (0..3).fold(Rational::from_i64(translation[i]), |sum, j| {
            sum + matrix[i][j] * &Rational::from_ratio(size[j], 2).unwrap()
        });
        assert!(center[i].contains(&expected));
        for j in 0..3 {
            let expected = (0..3).fold(b(0), |sum, k| {
                &sum + &(&principal[k] * &Ball::from_rational(&(matrix[i][k] * matrix[j][k]), BITS))
            });
            contains(&tensor[i][j], &expected);
        }
    }
    let mut exhausted = Decider::new(BudgetPolicy {
        max_work: 0,
        ..BudgetPolicy::default()
    });
    assert_eq!(
        mass_properties(
            &built.model,
            built.result,
            &Accuracy::default(),
            &mut exhausted
        )
        .unwrap_err()
        .code,
        "measure.budget-exhausted"
    );
}

#[test]
fn shelled_tray_algebraic_corner_trims_prove_mass_and_inertia() {
    let recipe = Recipe::parse(include_str!("../../../parts/P4-shelled-tray.json")).unwrap();
    let built = evaluate(&recipe, &mut Decider::default()).unwrap();
    let answer = mass_properties(
        &built.model,
        built.result,
        &Accuracy::default(),
        &mut Decider::default(),
    )
    .unwrap();
    assert!(answer.marks.is_empty(), "{:?}", answer.marks);
    let mass = answer.value;
    // Integrate horizontal sections of the clipped rectangle, then subtract
    // the rectangular cavity. These are not boundary flux formulas.
    let r2 = b(529);
    let x0 = b(304).sqrt().unwrap();
    let y0 = b(129).sqrt().unwrap();
    let area_primitive = |x: &Ball| {
        let y = (&r2 - &(x * x)).sqrt().unwrap();
        let angle = x.checked_div(&b(23)).unwrap().asin().unwrap();
        (&(x * &y) + &(&r2 * &angle)).checked_div(&b(2)).unwrap()
    };
    let outer_area = &b(4) * &(&(&b(15) * &x0) + &(area_primitive(&b(20)) - area_primitive(&x0)));
    let volume = &(&b(12) * &outer_area) - &b(9360);
    let first_z = &(&b(72) * &outer_area) - &b(65520);
    let raw_z2 = &(&b(576) * &outer_area) - &b(536640);
    let centered_z2 = &raw_z2 - &(&first_z * &first_z).checked_div(&volume).unwrap();
    let planar_second = |a: &Ball, bb: &Ball, cut: &Ball| {
        let primitive = |x: &Ball| {
            let y = (&r2 - &(x * x)).sqrt().unwrap();
            let angle = x.checked_div(&b(23)).unwrap().asin().unwrap();
            (&(&(x * &y) * &(&(&b(2) * &(x * x)) - &r2)) + &(&(&r2 * &r2) * &angle))
                .checked_div(&b(8))
                .unwrap()
        };
        &b(4)
            * &(&(&(bb * cut) * &(cut * cut)).checked_div(&b(3)).unwrap()
                + &(primitive(a) - primitive(cut)))
    };
    let raw_x2 = &(&b(12) * &planar_second(&b(20), &b(15), &x0)) - &b(1010880);
    let raw_y2 = &(&b(12) * &planar_second(&b(15), &b(20), &y0)) - &b(527280);
    let angle = b(20).checked_div(&b(23)).unwrap().asin().unwrap()
        - x0.checked_div(&b(23)).unwrap().asin().unwrap();
    let perimeter = &b(4) * &(&(&x0 + &y0) + &(&b(23) * &angle));
    let area = &(&b(2) * &outer_area) + &(&b(12) * &perimeter) + b(1240);
    contains(&mass.volume, &volume);
    contains(&mass.area, &area);
    let centroid = mass.centroid.unwrap();
    assert!(centroid[0].contains(&Rational::zero()));
    assert!(centroid[1].contains(&Rational::zero()));
    contains(&centroid[2], &first_z.checked_div(&volume).unwrap());
    let expected = [
        &raw_y2 + &centered_z2,
        &raw_x2 + &centered_z2,
        &raw_x2 + &raw_y2,
    ];
    let tensor = mass.inertia.unwrap();
    let zero = b(0);
    for (i, row) in tensor.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            contains(value, if i == j { &expected[i] } else { &zero });
        }
    }
}
