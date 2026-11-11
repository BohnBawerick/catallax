//! Limits on contact constructions do not certify a largest fitting feature.

use catallax_numbers::{Decider, Rational, Real, Refusal};

pub(crate) fn number(value: &Rational) -> String {
    value
        .to_decimal_exact()
        .unwrap_or_else(|| value.to_string())
}

/// Proof text and the lower endpoint, also used by the legacy radius hint.
pub(crate) fn bound(value: &Real, dec: &mut Decider) -> Option<(String, Rational)> {
    if let Some(exact) = value.as_rational() {
        return Some((format!("exactly {} mm", number(&exact)), exact));
    }
    let interval = dec
        .enclose(value, &Rational::from_ratio(1, 1_000_000)?)
        .ok()?;
    Some((
        format!(
            "in [{}, {}] mm (certified enclosure)",
            number(&interval.lo()),
            number(&interval.hi())
        ),
        interval.lo(),
    ))
}

pub(crate) fn unknown(mut refusal: Refusal) -> Refusal {
    if !refusal.detail.contains("; maximum unknown") {
        refusal.detail.push_str("; maximum unknown");
    }
    refusal
}

pub(crate) fn thickness(
    mut refusal: Refusal,
    limit: &Real,
    collision: &str,
    dec: &mut Decider,
) -> Refusal {
    if let Some((proof, _)) = bound(limit, dec) {
        refusal.detail.push_str(&format!(
            "; regular carrier thickness bound {proof}; collision: {collision}; equality needs continuation"
        ));
    }
    unknown(refusal)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn irrational_contact_bounds_are_enclosed_on_both_sides() {
        let mut dec = Decider::default();
        let value = Real::from_i64(2).sqrt();
        let (text, lo) = bound(&value, &mut dec).unwrap();
        assert!(text.contains("certified enclosure"));
        let (_, tail) = text.split_once(", ").unwrap();
        let hi = Rational::parse(tail.split(']').next().unwrap()).unwrap();
        let two = Rational::from_i64(2);
        assert!(&lo * &lo <= two);
        assert!(&hi * &hi >= two);
        assert!(&hi - &lo <= Rational::from_ratio(1, 1_000_000).unwrap());
    }
}
