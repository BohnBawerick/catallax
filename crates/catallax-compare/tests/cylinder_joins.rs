//! Bounded cylinder joins retain matching finite domains, normals and origins.
use catallax_compare::{Completion, Controller, Verdict};
use catallax_numbers::with_operation_budget;
use catallax_recipe::Recipe;
#[test]
fn cylinder_join_offsets_and_single_face_sweeps_agree() {
    for text in [
        include_str!("../../catallax-features/tests/fixtures/cylinder-joins/boss-plate.json"),
        include_str!("../../catallax-features/tests/fixtures/cylinder-joins/stepped-shaft.json"),
        include_str!(
            "../../catallax-features/tests/fixtures/cylinder-joins/boss-plate-thicken.json"
        ),
        include_str!(
            "../../catallax-features/tests/fixtures/cylinder-joins/stepped-shaft-thicken.json"
        ),
    ] {
        let r = Recipe::parse(text).unwrap();
        let controller = Controller::default();
        let revision = controller.begin();
        let (completion, work) =
            with_operation_budget(100_000_000, || revision.check(&r, &Default::default()));
        assert_eq!(completion, Completion::Published);
        assert!(!work.exhausted, "{} {work:?}", r.name);
        let report = controller.latest().unwrap();
        assert_eq!(
            report.verdict,
            Verdict::Agree,
            "{}: {:?} {:?}",
            r.name,
            report.findings,
            report.not_comparable
        );
        assert!(report.findings.is_empty());
        assert!(report.not_comparable.is_empty());
    }
}
