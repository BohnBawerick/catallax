//! Numerical torus admission precedes STEP seam classification.
#[path = "../../catallax-topology/tests/common/toroids.rs"]
mod fixture;
use catallax_numbers::Decider;
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::{BodyKind, numerical as n};

#[test]
fn an_open_torus_sheet_keeps_surface_kind_and_toroidal_report_authority() {
    let body = n::build_surface(
        fixture::quarter(),
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
    assert_eq!(file.report().kind, BodyKind::Surface);
    assert_eq!(file.report().counts, [1, 4, 4]);
    assert!(
        file.report()
            .checks
            .iter()
            .any(|check| check == "file.explicit-toroidal-pcurves")
    );
    assert!(
        !file
            .report()
            .checks
            .iter()
            .any(|check| check == "file.explicit-cylindrical-pcurves"
                || check == "file.explicit-conical-pcurves")
    );
    let text = std::str::from_utf8(file.bytes()).unwrap();
    assert!(text.contains("SHELL_BASED_SURFACE_MODEL("));
    assert!(!text.contains("MANIFOLD_SOLID_BREP("));
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
    if let Ok(path) = std::env::var("CATALLAX_ANALYTIC_EVIDENCE") {
        std::fs::write(
            std::path::Path::new(&path).join("open-torus-sheet.step"),
            file.bytes(),
        )
        .unwrap();
    }
}

#[test]
fn every_wrapped_contact_requires_a_certificate_not_just_the_first_seam() {
    let error = n::build_surface(
        fixture::unpaired_wrapped_contact(),
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "numerical.torus-periodic-domain-unverified");
}

#[test]
fn a_subperiod_repeated_edge_refuses_before_a_complete_torus_seam_is_written() {
    let mut partial = fixture::band();
    let period = fixture::period();
    for trim in &mut partial.faces[0].loops[0] {
        for uv in &mut trim.uv {
            if uv[0] == period {
                uv[0] = &period / catallax_numbers::Rational::from_i64(4);
            }
        }
    }
    let error = n::build_surface(
        partial,
        n::NumericalOptions::default(),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(error.code, "numerical.torus-periodic-domain-unverified");
    let body = n::build_surface(
        fixture::band(),
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
    assert_eq!(
        std::str::from_utf8(file.bytes())
            .unwrap()
            .matches("=SEAM_CURVE(")
            .count(),
        1
    );
    file.check_written(file.bytes(), &mut Decider::default())
        .unwrap();
}
