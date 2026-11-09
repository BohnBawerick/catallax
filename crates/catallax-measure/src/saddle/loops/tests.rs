use super::*;
use catallax_numbers::{Real, with_operation_budget};
use catallax_shapes::{Point3, Vector3};

#[test]
fn crossing_densities_have_the_required_reflection_parities() {
    for cosine in ["0", "3/5", "-3/5"] {
        for branch in [-1, 1] {
            for narrow_positive in [false, true] {
                let ((), used) = with_operation_budget(1_000_000, || {
                    let value = |text: &str| Real::from(Rational::parse(text).unwrap());
                    let law = SaddleLaw {
                        foot: Point3::origin(),
                        across: Vector3::new(value("1"), value("0"), value("0")),
                        along: Vector3::new(value("0"), value("1"), value("0")),
                        axis: Vector3::new(value("0"), value("0"), value("1")),
                        wide: value("5"),
                        narrow: value("5"),
                        radius: value("1/2"),
                        gap: value("0"),
                        cosine: value(cosine),
                        sine: value(if cosine == "0" { "1" } else { "4/5" }),
                        sides: [1, 1],
                        branch,
                        ring: SaddleRing::Crossing { narrow_positive },
                        gap_sign: 1,
                    };
                    let k = Builder::default();
                    let reference = build(&k, &law, true, false, true);
                    let optimized = forms(&law, true, true);
                    let dec = Decider::default();
                    let policy = *dec.policy();
                    let mut work = Work::new(&policy, BITS, policy.max_work, &dec).unwrap();
                    let rule = GaussRule::new(&mut work).unwrap();
                    let a = Rational::from_ratio(1, 8).unwrap();
                    let b = Rational::from_ratio(3, 16).unwrap();
                    let positive = reference
                        .turning_gauss_panel(&a, &b, &rule, &mut work)
                        .unwrap();
                    let negative = reference
                        .turning_gauss_panel(&-&b, &-&a, &rule, &mut work)
                        .unwrap();
                    let pruned = optimized
                        .turning_gauss_panel(&a, &b, &rule, &mut work)
                        .unwrap();
                    for (index, (positive, negative)) in positive.iter().zip(&negative).enumerate()
                    {
                        let odd = [7, 10, 12].contains(&index);
                        let residual = if odd {
                            positive + negative
                        } else {
                            positive - negative
                        };
                        assert!(
                            residual.contains_zero(),
                            "{cosine} {branch} {narrow_positive} {index}"
                        );
                        assert!(
                            residual.hi() - residual.lo() < Rational::parse("1e-8").unwrap(),
                            "{cosine} {branch} {narrow_positive} {index}: [{}, {}]",
                            residual.lo(),
                            residual.hi()
                        );
                        if odd {
                            assert!(pruned[index].is_exact_zero());
                        } else {
                            assert!(
                                pruned[index].lo() <= positive.hi()
                                    && positive.lo() <= pruned[index].hi()
                            );
                        }
                    }
                });
                assert!(!used.exhausted);
            }
        }
    }
}
