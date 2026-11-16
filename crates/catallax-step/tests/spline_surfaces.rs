use catallax_numbers::{Decider, Rational};
use catallax_shapes::spline_tools::RationalSplineSurface;
use catallax_step::write_spline_surfaces_step;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}

fn rounding_witness(dec: &Decider) -> RationalSplineSurface {
    RationalSplineSurface::new(
        [1, 1],
        [
            vec![q("0"), q("0"), q("1"), q("1")],
            vec![q("0"), q("0"), q("1"), q("1")],
        ],
        [2, 2],
        vec![
            [q("0"), q("0"), q("0")],
            [q("0"), q("1"), q("0")],
            [q("1.04104"), q("0"), q("0")],
            [q("1.04104"), q("1"), q("0")],
        ],
        vec![q("1"); 4],
        dec,
    )
    .unwrap()
}

#[test]
fn serialized_uncertainty_covers_actual_control_rounding_without_exceeding_the_cap() {
    let mut dec = Decider::default();
    let source = rounding_witness(&dec);
    let cap = q("0.041049");
    let mut bytes = Vec::new();
    let report = write_spline_surfaces_step(&[source], &cap, 2, &mut dec, &mut bytes).unwrap();
    let file = std::str::from_utf8(&bytes).unwrap();
    let value = file
        .split("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(")
        .nth(1)
        .unwrap()
        .split(')')
        .next()
        .unwrap();
    let serialized = q(value);
    assert_eq!(report.max_deviation_mm, q("0.04104"));
    assert_eq!(serialized, report.uncertainty_mm);
    assert!(serialized >= q("0.04104"));
    assert!(serialized <= cap);
    assert!(!report.exact);
}

#[test]
fn excessive_control_rounding_refuses_without_writing_a_partial_file() {
    let mut dec = Decider::default();
    let source = rounding_witness(&dec);
    let mut bytes = b"original output".to_vec();
    let refusal =
        write_spline_surfaces_step(&[source], &q("0.04103"), 2, &mut dec, &mut bytes).unwrap_err();
    assert_eq!(refusal.code, "step.spline-surfaces");
    assert_eq!(bytes, b"original output");
}
