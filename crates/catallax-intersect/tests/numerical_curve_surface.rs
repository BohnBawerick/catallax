#![allow(missing_docs)]
#[path = "fixtures/p8_native_curve_patch.rs"]
mod fixture;

use catallax_intersect::{numerical_curve_surface::patch_meetings, trace};
use catallax_numbers::{Ball, Decider, Outcome, Rational, ZeroProof, with_operation_budget};
use catallax_shapes::{Curve, GeomStore, PointDef, SectionHits};

#[test]
fn captured_native_shared_boundary_preserves_exact_endpoint() {
    let ((), work) = with_operation_budget(100_000_000, || {
        let mut dec = Decider::default();
        let (source, cut, target) = fixture::patches(&mut dec);
        let section = source.sections(&cut, &mut dec).unwrap();
        let [branch] = section.branches() else {
            panic!("the captured source has one component");
        };
        let answer = branch.patch_hits(&target, &mut dec).unwrap().unwrap();
        let SectionHits::Points(points) = answer else {
            panic!("the adjacent native patch meets only at the shared endpoint");
        };
        assert_eq!(points.len(), 1);
        let expected = branch.endpoints()[1].clone();
        match dec.sign(&points[0].minus(&expected).norm_squared()) {
            Outcome::ProvenCoincident(proof) => {
                assert_ne!(proof.proof, ZeroProof::SeparationBound);
            }
            other => panic!("the original exact endpoint changed: {other:?}"),
        }
        let diagonal = fixture::diagonal_target(&mut dec);
        let Some(SectionHits::Points(points)) = branch.patch_hits(&diagonal, &mut dec).unwrap() else {
            panic!("the captured diagonal patch strictly misses the retained cylinder");
        };
        assert!(points.is_empty());
        let opposite = fixture::opposite_target(&mut dec);
        let bounds = branch.control_hull(&mut dec).unwrap();
        for endpoint in branch.endpoints() {
            for (axis, value) in [&endpoint.x, &endpoint.y, &endpoint.z].into_iter().enumerate() {
                let point = dec.enclose_budgeted_at(value, 128).unwrap().to_interval();
                assert!(
                    bounds[axis].lo() <= point.hi() && point.lo() <= bounds[axis].hi(),
                    "{bounds:?}: endpoint axis {axis} is outside the component enclosure"
                );
            }
        }
        let opposite_top = opposite
            .carrier()
            .controls()
            .iter()
            .map(|point| point.z.approx().hi())
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(bounds[2].lo() > opposite_top, "{bounds:?}");
        eprintln!("native component bounds: {bounds:?}");
        assert!(branch.patch_hits(&source, &mut dec).unwrap().is_none());
    });
    assert!(!work.exhausted, "{work:?}");
    eprintln!("exact shared-boundary work: {work:?}");
}

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
