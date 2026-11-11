//! Original public controls; no private inputs or saved-answer inventory.
#[path = "../examples/reconstruction/fixture.rs"]
mod fixture;
use catallax_field::{
    Location, Workspace,
    reconstruction::{self, Sampling},
};
use catallax_numbers::{Decider, Outcome, Rational};
use catallax_shapes::RPoint3;
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn sampling() -> Sampling {
    Sampling {
        triangle_steps: 4,
        surface_steps: 8,
        angles: 128,
    }
}
#[test]
fn cylinder_and_cone_derive_clean_parts_and_keep_the_hole() {
    for top in [2, 3] {
        let bytes = fixture::fixture(top);
        let saved: &[u8] = if top == 2 {
            include_bytes!("../examples/reconstruction/box-hole.obj")
        } else {
            include_bytes!("../examples/reconstruction/cone-hole.obj")
        };
        assert_eq!(bytes, saved, "committed fixture must match its generator");
        let mut dec = Decider::default();
        let part = reconstruction::reconstruct(&bytes, &q(1), &mut dec).unwrap();
        assert_eq!(part.radii(), &[q(2), q(top)]);
        assert_eq!(part.centre(), &[q(0), q(0)]);
        assert_eq!(part.region_triangles().iter().sum::<usize>(), 192);
        let mut brep = part.brep().clone();
        assert!(brep.model.check(brep.result).is_empty());
        assert_eq!(brep.model.faces_of_body(brep.result).len(), 7);
        let volume =
            catallax_measure::volume(&brep.model, brep.result, &Default::default(), &mut dec)
                .unwrap()
                .value;
        let expected_volume = 600. - std::f64::consts::PI * if top == 2 { 24. } else { 38. };
        assert!((volume.lo.to_f64_nearest() - expected_volume).abs() < 1e-6);
        assert!((volume.hi.to_f64_nearest() - expected_volume).abs() < 1e-6);
        eprintln!(
            "top={top} B-rep volume mm^3=[{},{}] analytic={expected_volume}",
            volume.lo.to_f64_nearest(),
            volume.hi.to_f64_nearest()
        );
        for (p, expected) in [
            ([0, 0, 3], Location::Outside),
            ([4, 0, 3], Location::Inside),
            ([6, 0, 3], Location::Outside),
        ] {
            let p = RPoint3::new(q(p[0]), q(p[1]), q(p[2]));
            assert_eq!(
                part.field()
                    .proven_checked(&p, &mut Workspace::default(), &mut dec)
                    .unwrap()
                    .location,
                expected
            );
            let answer =
                catallax_boolean::classify_point(&mut brep.model, brep.result, &p, &mut dec);
            let wanted = if expected == Location::Inside {
                catallax_boolean::Containment::Inside
            } else {
                catallax_boolean::Containment::Outside
            };
            assert!(
                matches!(answer,Outcome::Proven(value) if value==wanted),
                "{answer:?}"
            );
        }
        let d = reconstruction::sampled_deviation(&bytes, &part, sampling(), &mut dec).unwrap();
        eprintln!("top={top} normal={d:?}");
        assert!(d.scan_to_part_mm > 0. && d.scan_to_part_mm < 0.04);
        assert!(d.part_to_scan_mm > 0. && d.part_to_scan_mm < 0.04);
        let adverse =
            reconstruction::filled_envelope_control(&bytes, &part, sampling(), &mut dec).unwrap();
        eprintln!("top={top} missing bore={adverse:?}");
        assert!((adverse.scan_to_part_mm - 3.).abs() < 1e-12);
        assert!(adverse.part_to_scan_mm > 1.9 && adverse.part_to_scan_mm < 3.);
        assert!((adverse.scan_to_part_mm - adverse.part_to_scan_mm).abs() > 0.05);
    }
}
#[test]
fn scale_and_translation_are_inferred_not_fixture_dimensions() {
    let original = String::from_utf8(fixture::fixture(2)).unwrap();
    let shifted = original
        .lines()
        .map(|line| {
            if let Some(row) = line.strip_prefix("v ") {
                let p: Vec<_> = row
                    .split_whitespace()
                    .map(|s| Rational::parse(s).unwrap())
                    .collect();
                format!(
                    "v {} {} {}",
                    catallax_recipe::Number::new(&p[0] + &q(7)),
                    catallax_recipe::Number::new(&p[1] - &q(4)),
                    catallax_recipe::Number::new(&p[2] + &q(2))
                )
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let part =
        reconstruction::reconstruct(shifted.as_bytes(), &q(2), &mut Decider::default()).unwrap();
    assert_eq!(part.centre(), &[q(14), q(-8)]);
    assert_eq!(part.radii(), &[q(4), q(4)]);
    assert_eq!(part.bounds()[2], [q(4), q(16)]);
}
#[test]
fn invalid_unmatched_and_source_identity_refuse() {
    let bytes = fixture::fixture(2);
    let text = String::from_utf8(bytes.clone()).unwrap();
    let mut lines: Vec<_> = text.lines().collect();
    lines.pop();
    assert!(
        reconstruction::reconstruct(lines.join("\n").as_bytes(), &q(1), &mut Decider::default())
            .is_err()
    );
    // Move one bore vertex radially; the mesh remains a source candidate, but
    // cannot satisfy the exact circle fit. Never discard the unmatched triangles.
    let changed = text.replacen("v 2 0 0", "v 2.01 0 0", 1);
    let error = reconstruction::reconstruct(changed.as_bytes(), &q(1), &mut Decider::default())
        .err()
        .unwrap();
    assert_eq!(error.code, "reconstruction.unsupported");
    let part = reconstruction::reconstruct(&bytes, &q(1), &mut Decider::default()).unwrap();
    assert!(
        reconstruction::sampled_deviation(
            changed.as_bytes(),
            &part,
            sampling(),
            &mut Decider::default()
        )
        .is_err()
    );
}
#[test]
fn independent_one_two_four_workers_keep_definition_and_sample_bits() {
    let bytes = fixture::fixture(2);
    let record = || {
        let mut dec = Decider::default();
        let part = reconstruction::reconstruct(&bytes, &q(1), &mut dec).unwrap();
        let d = reconstruction::sampled_deviation(&bytes, &part, sampling(), &mut dec).unwrap();
        (
            part.recipe().to_text(),
            d.scan_to_part_mm.to_bits(),
            d.part_to_scan_mm.to_bits(),
            d.scan_samples,
            d.part_samples,
            dec.stats().work,
        )
    };
    let expected = record();
    for workers in [1, 2, 4] {
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..workers).map(|_| scope.spawn(record)).collect();
            assert_eq!(handles.len(), workers);
            for handle in handles {
                assert_eq!(handle.join().unwrap(), expected);
            }
        });
    }
}

#[test]
fn source_polyhedron_and_fitted_solid_keep_distinct_meanings() {
    let bytes = fixture::fixture(2);
    let mut dec = Decider::default();
    let source = catallax_import::mesh::ClosedMesh::read(&bytes, &q(1), &mut dec).unwrap();
    let part = reconstruction::reconstruct(&bytes, &q(1), &mut dec).unwrap();
    let p = [
        Rational::parse("1.98").unwrap(),
        Rational::parse("0.28").unwrap(),
        q(3),
    ];
    assert_eq!(source.winding(&p, &mut dec).unwrap(), Some(1));
    assert_eq!(
        part.field()
            .proven_checked(
                &RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()),
                &mut Workspace::default(),
                &mut dec
            )
            .unwrap()
            .location,
        Location::Outside
    );
}
#[test]
fn extra_component_is_refused_instead_of_discarded() {
    let original = String::from_utf8(fixture::fixture(2)).unwrap();
    let count = original.lines().filter(|s| s.starts_with("v ")).count();
    let mut extra = String::new();
    for line in original.lines() {
        if let Some(row) = line.strip_prefix("v ") {
            let mut p: Vec<_> = row
                .split_whitespace()
                .map(|s| Rational::parse(s).unwrap())
                .collect();
            p[0] = &p[0] + &q(20);
            extra.push_str(&format!(
                "v {} {} {}\n",
                catallax_recipe::Number::new(p[0].clone()),
                catallax_recipe::Number::new(p[1].clone()),
                catallax_recipe::Number::new(p[2].clone())
            ));
        } else if let Some(row) = line.strip_prefix("f ") {
            let ids: Vec<_> = row
                .split_whitespace()
                .map(|s| s.parse::<usize>().unwrap() + count)
                .collect();
            extra.push_str(&format!("f {} {} {}\n", ids[0], ids[1], ids[2]));
        }
    }
    let error = reconstruction::reconstruct(
        (original + &extra).as_bytes(),
        &q(1),
        &mut Decider::default(),
    )
    .err()
    .unwrap();
    assert_eq!(error.code, "reconstruction.unsupported");
}
