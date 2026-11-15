//! One shared cylinder seam carries two explicit periodic pcurves.
#[path = "../../catallax-topology/tests/common/periodic.rs"]
mod fixture;
use catallax_numbers::Decider;
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;

#[test]
fn a_complete_bore_emits_one_seam_curve_with_two_uses() {
    let body = n::build(
        fixture::bore(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let file = output::build(
        &body,
        &NumericalStepOptions::default(),
        &mut Decider::default(),
    )
    .unwrap();
    let text = std::str::from_utf8(file.bytes()).unwrap();
    assert_eq!(text.matches("=SEAM_CURVE(").count(), 1);
    assert_eq!(file.report().counts, [7, 15, 10]);
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&path).join("periodic-bore.step"),
            file.bytes(),
        )
        .unwrap();
    }
}
