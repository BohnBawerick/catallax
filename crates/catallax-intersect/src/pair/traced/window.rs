//! Finite cone heights when no cone generator is parallel to the cylinder axis.

use super::*;

pub(super) fn height_upper(
    slope: &Rational,
    radius: &Rational,
    radial: &RVector3,
    cosine: &Rational,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    // For p=A+h(d+s*u), u perpendicular to d, the cylinder-axis projection
    // has length at least h*|sin(theta)-s*|cos(theta)||. The cylinder radius
    // and the apex's radial distance bound that same length from above.
    let sine = Real::from(&r(1) - &(cosine * cosine)).sqrt();
    let gap = &sine - &Real::from(slope * &cosine.abs());
    let lower = if let Some(q) = gap.as_rational() {
        q.abs()
    } else {
        let bound = dec.enclose_budgeted_at(&gap, 128.min(dec.policy().max_precision_bits))?;
        if bound.lo().sign() == Sign::Positive {
            bound.lo()
        } else if bound.hi().sign() == Sign::Negative {
            -bound.hi()
        } else {
            r(0)
        }
    };
    if lower.sign() != Sign::Positive {
        return Err(Refusal::new(
            "intersect.trace-native-window",
            "a cone generator may parallel the cylinder axis; finite face domains are required",
        ));
    }
    dec.charge_work(1)?;
    let distance_upper = &(&radial.x.abs() + &radial.y.abs()) + &radial.z.abs();
    Ok((&distance_upper + radius) / lower)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parallel_and_oblique_windows_bound_height_without_perpendicular_axes() {
        let radial = RVector3 {
            x: Rational::from_ratio(21, 2).unwrap(),
            y: r(0),
            z: r(0),
        };
        for cosine in [r(1), r(-1)] {
            assert_eq!(
                height_upper(&r(1), &r(1), &radial, &cosine, &mut Decider::default()).unwrap(),
                Rational::from_ratio(23, 2).unwrap()
            );
        }
        for cosine in [
            Rational::from_ratio(3, 5).unwrap(),
            Rational::from_ratio(-4, 5).unwrap(),
        ] {
            let upper =
                height_upper(&r(1), &r(1), &radial, &cosine, &mut Decider::default()).unwrap();
            assert!(upper >= Rational::from_ratio(115, 2).unwrap() && upper < r(58));
        }
        let refusal = height_upper(
            &Rational::from_ratio(4, 3).unwrap(),
            &r(1),
            &radial,
            &Rational::from_ratio(3, 5).unwrap(),
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(refusal.code, "intersect.trace-native-window");
    }
}
