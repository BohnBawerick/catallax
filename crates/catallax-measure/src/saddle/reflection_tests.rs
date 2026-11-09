use super::*;
use catallax_numbers::Real;
use catallax_shapes::{Point3, Vector3};

fn law(cosine: &str, branch: i64, positive: bool) -> SaddleLaw {
    let real = |s: &str| Real::from(Rational::parse(s).unwrap());
    SaddleLaw {
        foot: Point3::origin(),
        across: Vector3::new(real("1"), real("0"), real("0")),
        along: Vector3::new(real("0"), real("1"), real("0")),
        axis: Vector3::new(real("0"), real("0"), real("1")),
        wide: real("5"),
        narrow: real("5"),
        radius: real("1/2"),
        gap: real("0"),
        cosine: real(cosine),
        sine: real(if cosine == "0" { "1" } else { "4/5" }),
        sides: [1, 1],
        branch,
        ring: SaddleRing::Crossing {
            narrow_positive: positive,
        },
        gap_sign: 1,
    }
}

#[test]
fn reflected_crossing_keys_map_the_exact_centre_ellipse() {
    let position = |law: &SaddleLaw| {
        let SaddleRing::Crossing { narrow_positive } = law.ring else {
            unreachable!()
        };
        let radius = law.wide.as_rational().unwrap() + law.radius.as_rational().unwrap();
        let sin = Rational::from_ratio(3, 5).unwrap();
        let cos = Rational::from_ratio(4, 5).unwrap();
        let y = &radius * &Rational::from_i64(law.branch) * cos.clone();
        let narrow = &radius * &Rational::from_i64(if narrow_positive { 1 } else { -1 }) * cos;
        let z = (law.cosine.as_rational().unwrap() * y.clone() - narrow)
            / law.sine.as_rational().unwrap();
        [radius * sin, y, z]
    };
    for cosine in ["0", "3/5", "-3/5"] {
        let source = law(cosine, 1, true);
        for branch in [-1, 1] {
            for positive in [false, true] {
                let target = law(cosine, branch, positive);
                let reflection = crossing_reflection(&source, &target);
                if cosine != "0" && (branch == 1) != positive {
                    assert!(reflection.is_none());
                } else {
                    let signs = reflection.unwrap();
                    let reflected = position(&source);
                    let expected = position(&target);
                    for i in 0..3 {
                        assert_eq!(&reflected[i] * &Rational::from_i64(signs[i]), expected[i]);
                    }
                }
            }
        }
        let mut different = source.clone();
        different.radius = Real::one();
        assert!(crossing_reflection(&source, &different).is_none());
        different = source.clone();
        different.gap = Real::one();
        assert!(crossing_reflection(&source, &different).is_none());
        different = source.clone();
        different.ring = SaddleRing::Through;
        assert!(crossing_reflection(&source, &different).is_none());
    }
}
