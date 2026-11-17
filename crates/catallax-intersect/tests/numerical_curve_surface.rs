#![allow(missing_docs)]
#[path = "fixtures/p8_native_curve_patch.rs"]
mod fixture;

use catallax_intersect::{numerical_curve_surface::patch_meetings, trace};
use catallax_numbers::{Ball, Decider, Rational, with_operation_budget};
use catallax_shapes::{Curve, GeomStore, PointDef};

#[test]
fn captured_native_curve_accepts_crossing_and_refuses_tangency_and_overlap() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (source, cut, _) = fixture::patches(&mut dec);
        let section = source.sections(&cut, &mut dec).unwrap();
        assert_eq!(section.branches().len(), 1);
        let branch = &section.branches()[0];
        let [a, b] = branch.endpoints();
        let expected = [
            dec.enclose_budgeted_at(&b.x, 128).unwrap(),
            dec.enclose_budgeted_at(&b.y, 128).unwrap(),
            dec.enclose_budgeted_at(&b.z, 128).unwrap(),
        ];
        let mut geom = GeomStore::new();
        let curve = geom.add_curve(Curve::PatchSection(branch.clone()));
        let from = geom.add_point(PointDef::Explicit(a));
        let to = geom.add_point(PointDef::Explicit(b));
        let traced = trace(&geom, curve, Some(from), Some(to), &mut dec).unwrap();
        let tolerance = Rational::parse("0.01").unwrap();
        for (name, target) in fixture::variants(&source, branch, &mut dec) {
            let result = patch_meetings(&traced.value, &target, &tolerance, &mut dec);
            if name == "transverse-boundary" {
                let answer = result.unwrap_or_else(|error| panic!("{name}: {error:?}"));
                assert!(!answer.hits.is_empty());
                for hit in answer.hits {
                    let displacement = (0..3).fold(Rational::zero(), |sum, axis| {
                        let delta: Ball = &hit.position[axis] - &expected[axis];
                        sum + delta.lo().abs().max(delta.hi().abs())
                    });
                    assert!(displacement <= tolerance, "{name}: {displacement}");
                    assert!(!hit.normal_dot.contains_zero());
                }
                assert_eq!(answer.evidence.tolerance_mm(), &tolerance);
            } else {
                let error = result
                    .expect_err("unowned tangent/overlap contacts must not publish finite roots");
                assert!(
                    matches!(
                        error.code,
                        "intersect.numerical-nontransverse"
                            | "intersect.numerical-overlap"
                            | "intersect.numerical-search-incomplete"
                    ),
                    "{name}: {error:?}"
                );
                eprintln!("{name}: {error:?}");
            }
        }
    });
    assert!(!work.exhausted, "{work:?}");
}
