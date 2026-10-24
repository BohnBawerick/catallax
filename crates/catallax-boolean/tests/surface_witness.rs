//! Rational query witnesses remain close to rounded sphere poles.

use catallax_boolean::surface_witness;
use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
use catallax_shapes::{Point3, RPoint3, RVector3, Sphere, Surface, Vector3};

fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn sphere(center: &RPoint3, radius: &Rational) -> Surface {
    Surface::Sphere(
        Sphere::new(
            center.to_real(),
            radius.clone().into(),
            Vector3::new(Real::zero(), Real::zero(), Real::one()),
            Vector3::new(Real::one(), Real::zero(), Real::zero()),
            &mut Decider::default(),
        )
        .unwrap(),
    )
}

fn witness(surface: &Surface, target: &Point3, width: &Rational) -> RPoint3 {
    let point = surface_witness(surface, target, width, &mut Decider::default())
        .unwrap()
        .unwrap();
    let p = point.to_real();
    let residual = surface.implicit().unwrap().eval_real(&[p.x, p.y, p.z]);
    assert!(residual.rational_ref().unwrap().is_zero());
    point
}

#[test]
fn rounded_targets_on_both_sides_of_all_six_poles_keep_nearby_witnesses() {
    let width = Rational::one().mul_pow2(-64);
    let delta = Rational::one().mul_pow2(-48);
    for (center, radius) in [
        (RPoint3::origin(), r("1")),
        (RPoint3::new(r("3"), r("-2"), r("5")), r("2")),
    ] {
        let surface = sphere(&center, &radius);
        for axis in 0..3 {
            for sense in [-1, 1] {
                let mut xyz = [r("0"), r("0"), r("0")];
                xyz[axis] = Rational::from_i64(sense);
                let [x, y, z] = xyz;
                let direction = RVector3::new(x, y, z);
                for offset in [-delta.clone(), r("0"), delta.clone()] {
                    let target = center
                        .plus(&direction.scaled(&(&radius * &(Rational::one() + offset.clone()))));
                    let point = witness(&surface, &target.to_real(), &width);
                    let bound = &radius * &(&offset.abs() * &r("2") + &width * &r("8"));
                    assert!(point.minus(&target).norm_squared() <= &bound * &bound);
                    // Both axial poles must remain exact, including a rounded
                    // target outside the north pole's old singular chart.
                    if axis == 2 {
                        assert_eq!(point, center.plus(&direction.scaled(&radius)));
                    }
                }
            }
        }
    }
}

#[test]
fn nonpolar_rational_directions_keep_support_incidence_and_proximity() {
    let center = RPoint3::new(r("3"), r("-2"), r("5"));
    let radius = r("2");
    let surface = sphere(&center, &radius);
    let width = Rational::one().mul_pow2(-64);
    let delta = Rational::one().mul_pow2(-48);
    for [x, y, z] in [["3/5", "0", "4/5"], ["0", "-4/5", "-3/5"]] {
        let direction = RVector3::new(r(x), r(y), r(z));
        for offset in [-delta.clone(), r("0"), delta.clone()] {
            let target =
                center.plus(&direction.scaled(&(&radius * &(Rational::one() + offset.clone()))));
            let point = witness(&surface, &target.to_real(), &width);
            let bound = &radius * &(&offset.abs() * &r("2") + &width * &r("8"));
            assert!(point.minus(&target).norm_squared() <= &bound * &bound);
        }
    }
}

#[test]
fn translated_rounded_north_pole_matches_the_backward_cap_witness() {
    let surface = sphere(&RPoint3::new(r("0"), r("0"), r("5")), &r("2"));
    let target = RPoint3::new(r("0"), r("0"), r("3940649673949185/562949953421312"));
    assert_eq!(
        witness(&surface, &target.to_real(), &Rational::one().mul_pow2(-64)),
        RPoint3::new(r("0"), r("0"), r("7"))
    );
}

#[test]
fn sphere_query_witness_respects_whole_request_exhaustion() {
    let surface = sphere(&RPoint3::origin(), &r("1"));
    let target = RPoint3::new(r("0"), r("0"), r("2")).to_real();
    let (result, work) = with_operation_budget(0, || {
        surface_witness(
            &surface,
            &target,
            &Rational::one().mul_pow2(-64),
            &mut Decider::default(),
        )
    });
    assert!(work.exhausted);
    assert_eq!(result.unwrap_err().code, "numbers.operation-budget");
}
