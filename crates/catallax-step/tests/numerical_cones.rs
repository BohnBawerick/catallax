//! A tapered bore exercises the reversed cone face and both periodic seam uses.
#[allow(dead_code)]
#[path = "../../catallax-topology/tests/common/cones.rs"]
mod fixture;

use catallax_numbers::Decider;
use catallax_step::numerical::{self as output, NumericalStepOptions};
use catallax_topology::numerical as n;

#[test]
fn inner_cone_seams_export_from_every_cyclic_loop_start() {
    for start in 0..4 {
        let mut d = fixture::positive_bore();
        d.faces.last_mut().unwrap().loops[0].rotate_left(start);
        let body = n::build(d, n::NumericalOptions::default(), &mut Decider::default()).unwrap();
        let file = output::build(
            &body,
            &NumericalStepOptions::default(),
            &mut Decider::default(),
        )
        .unwrap();
        assert!(
            std::str::from_utf8(file.bytes())
                .unwrap()
                .contains("CONICAL_SURFACE(")
        );
        assert!(
            file.report()
                .checks
                .iter()
                .any(|c| c == "file.cone-step-v-axial-height")
        );
        file.check_written(file.bytes(), &mut Decider::default())
            .unwrap();
        if let Ok(directory) = std::env::var("CATALLAX_CONE_EVIDENCE") {
            std::fs::create_dir_all(&directory).unwrap();
            std::fs::write(
                std::path::Path::new(&directory).join(format!("cone-bore-loop-{start}.step")),
                file.bytes(),
            )
            .unwrap();
        }
    }
}
