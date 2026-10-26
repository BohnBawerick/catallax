//! Independent cubic power-basis oracle. It solves each coordinate derivative
//! algebraically; it never calls spline evaluation, jets, tracing or bounds.

use catallax_measure::{Accuracy, BoundingBox, Enclosure};
use catallax_numbers::{Decider, Outcome, Rational, Real, StrictSign};
use catallax_shapes::{BSplineCurve, OffsetSide, SplineSurface, Surface};

pub struct Extrema {
    pub min: [Enclosure; 3],
    pub max: [Enclosure; 3],
    pub endpoint_max: [Enclosure; 3],
}

fn sign(value: &Real, dec: &mut Decider) -> i8 {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => 1,
        Outcome::Proven(StrictSign::Negative) => -1,
        Outcome::ProvenCoincident(_) => 0,
        outcome => panic!("independent derivative-root predicate: {outcome:?}"),
    }
}

fn merge(slot: &mut Option<Enclosure>, value: Enclosure, minimum: bool) {
    match slot {
        None => *slot = Some(value),
        Some(old) if minimum => {
            old.lo = old.lo.clone().min(value.lo);
            old.hi = old.hi.clone().min(value.hi);
        }
        Some(old) => {
            old.lo = old.lo.clone().max(value.lo);
            old.hi = old.hi.clone().max(value.hi);
        }
    }
}

pub fn extrema(curve: &BSplineCurve) -> Extrema {
    assert_eq!(curve.degree(), 3);
    assert!(
        curve
            .weights()
            .all(|w| w.rational_ref() == Some(&Rational::one()))
    );
    let mut dec = Decider::default();
    let mut min = [None, None, None];
    let mut max = [None, None, None];
    let mut endpoint_max = [None, None, None];
    for (index, (span, _)) in curve.spans().enumerate() {
        // The checked family is composite cubic Bezier: each internal knot
        // occurs three times and each span has four Bernstein controls.
        assert_eq!(span, 3 + index * 3);
        let controls = std::array::from_fn::<_, 4, _>(|i| curve.control(index * 3 + i));
        for axis in 0..3 {
            let p = controls.map(|p| [&p.x, &p.y, &p.z][axis]);
            let a = &(p[3] - p[0]) + &(&(p[1] - p[2]) * &Real::from(3));
            let b = &(&(p[0] + p[2]) - &(p[1] * &Real::from(2))) * &Real::from(3);
            let c = &(p[1] - p[0]) * &Real::from(3);
            let qa = &a * &Real::from(3);
            let qb = &b * &Real::from(2);
            let mut roots = Vec::new();
            if sign(&qa, &mut dec) == 0 {
                if sign(&qb, &mut dec) != 0 {
                    roots.push(-(&c / &qb));
                }
            } else {
                let discriminant = &qb.square() - &(&(&qa * &c) * &Real::from(4));
                if sign(&discriminant, &mut dec) >= 0 {
                    let root = discriminant.sqrt();
                    let denominator = &qa * &Real::from(2);
                    let negative = -&qb;
                    roots.push(&(&negative - &root) / &denominator);
                    roots.push(&(&negative + &root) / &denominator);
                }
            }
            let mut values = vec![p[0].clone(), p[3].clone()];
            for t in roots {
                if sign(&t, &mut dec) > 0 && sign(&(&Real::one() - &t), &mut dec) > 0 {
                    let first = &a * &t + b.clone();
                    let second = &first * &t + c.clone();
                    values.push(&second * &t + p[0].clone());
                }
            }
            for (index, value) in values.iter().enumerate() {
                let ball = dec.enclose_budgeted_at(value, 512).unwrap();
                let value = Enclosure {
                    lo: ball.lo(),
                    hi: ball.hi(),
                };
                if index < 2 {
                    merge(&mut endpoint_max[axis], value.clone(), false);
                }
                merge(&mut min[axis], value.clone(), true);
                merge(&mut max[axis], value, false);
            }
        }
    }
    Extrema {
        min: min.map(Option::unwrap),
        max: max.map(Option::unwrap),
        endpoint_max: endpoint_max.map(Option::unwrap),
    }
}

pub fn check(actual: &BoundingBox, expected: &Extrema, accuracy: &Accuracy) {
    assert!(accuracy.relative.is_zero());
    for axis in 0..3 {
        assert!(
            actual.min[axis] <= expected.min[axis].lo,
            "minimum enclosure, axis {axis}"
        );
        assert!(
            actual.max[axis] >= expected.max[axis].hi,
            "maximum enclosure, axis {axis}"
        );
        assert!(
            &expected.min[axis].hi - &actual.min[axis] <= accuracy.absolute,
            "minimum width, axis {axis}"
        );
        assert!(
            &actual.max[axis] - &expected.max[axis].lo <= accuracy.absolute,
            "maximum width, axis {axis}"
        );
    }
}

/// Support function of the retained exact inward parallel ellipse at the lower
/// cap. This is a separate target control, never the represented-model oracle.
pub fn target_cap_max(spline: &SplineSurface) -> [Enclosure; 3] {
    assert_eq!(spline.definition().side(), OffsetSide::Inward);
    let Surface::EllipticCylinder(source) = spline.definition().original() else {
        unreachable!()
    };
    let major = source
        .major_dir
        .scaled(&(&Real::one() / &source.major_dir.norm_squared().sqrt()));
    let minor = spline.axis().cross(&major);
    let minor = minor.scaled(&(&Real::one() / &minor.norm_squared().sqrt()));
    let center = source.origin.plus(&spline.axis().scaled(spline.axial()[0]));
    std::array::from_fn(|i| {
        let a = [&major.x, &major.y, &major.z][i];
        let b = [&minor.x, &minor.y, &minor.z][i];
        let n = [&spline.axis().x, &spline.axis().y, &spline.axis().z][i];
        let radius = (&(&a.square() * &source.major_sq) + &(&b.square() * &source.minor_sq)).sqrt();
        let displacement = spline.definition().distance() * &(&Real::one() - &n.square()).sqrt();
        let value = [&center.x, &center.y, &center.z][i] + &(&radius - &displacement);
        let ball = Decider::default().enclose_budgeted_at(&value, 512).unwrap();
        Enclosure {
            lo: ball.lo(),
            hi: ball.hi(),
        }
    })
}
