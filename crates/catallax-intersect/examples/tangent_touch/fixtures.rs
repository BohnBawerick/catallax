use catallax_numbers::{Decider, MPoly, Rational, Real, Refusal};
use catallax_shapes::{BSplineSurface, Plane, Point3, RationalPatch, Vector3};

pub const CASES: &[&str] = &[
    "bowl",
    "saddle",
    "flat",
    "line",
    "flush",
    "off-centre",
    "several",
    "algebraic",
    "near-miss",
    "tiny-crossing",
    "transverse",
    "boundary",
    "hill",
    "rational-bowl",
    "touch-and-cut",
    "crease",
    "steep-off-centre",
];
pub const COPIES: &[&str] = &["original", "translated", "scaled", "rotated", "combined"];

pub fn q(n: i64, d: i64) -> Rational {
    Rational::from_ratio(n, d).unwrap()
}
fn constant(value: Rational) -> MPoly {
    MPoly::constant(2, value)
}
fn sq(value: &MPoly) -> MPoly {
    value * value
}

fn choose(n: usize, k: usize) -> Rational {
    (0..k).fold(Rational::one(), |v, i| {
        v * q((n - i) as i64, (i + 1) as i64)
    })
}

fn transform(value: &[Rational; 3], copy: &str, vector: bool) -> [Rational; 3] {
    let rotate = matches!(copy, "rotated" | "combined");
    let scale = if matches!(copy, "scaled" | "combined") {
        q(7, 3)
    } else {
        q(1, 1)
    };
    let mut value = if rotate {
        [
            &q(3, 5) * &value[0] + &q(4, 5) * &value[2],
            value[1].clone(),
            &q(-4, 5) * &value[0] + &q(3, 5) * &value[2],
        ]
    } else {
        value.clone()
    };
    for (i, v) in value.iter_mut().enumerate() {
        if !vector {
            *v = &*v * &scale;
        }
        if !vector && matches!(copy, "translated" | "combined") {
            *v = &*v + &[q(11, 3), q(-17, 5), q(29, 7)][i];
        }
    }
    value
}

fn transformed_point(value: &[Rational; 3], copy: &str) -> Point3 {
    let [x, y, z] = transform(value, copy, false).map(Real::from_rational);
    Point3::new(x, y, z)
}

pub fn build(case: &str, copy: &str, dec: &mut Decider) -> Result<(RationalPatch, Plane), Refusal> {
    let mut height = Rational::zero();
    let carrier = if case == "crease" {
        let mut controls = Vec::new();
        for (u, z) in [(-1, -1), (0, 0), (1, -1)] {
            for v in [-1, 1] {
                controls.push(transformed_point(&[q(u, 1), q(v, 1), q(z, 16)], copy));
            }
        }
        BSplineSurface::new(
            [1, 1],
            [
                [-1, -1, 0, 1, 1].map(Rational::from_i64).to_vec(),
                [-1, -1, 1, 1].map(Rational::from_i64).to_vec(),
            ],
            [3, 2],
            controls,
            vec![Real::one(); 6],
            dec,
        )?
    } else if case == "hill" {
        height = q(225, 256);
        let knots = [0, 0, 0, 0, 1, 2, 3, 3, 3, 3]
            .map(Rational::from_i64)
            .to_vec();
        let mut controls = Vec::new();
        for i in 0..6 {
            for j in 0..6 {
                controls.push(transformed_point(
                    &[
                        q(i, 1),
                        q(j, 1),
                        q(i64::from((2..=3).contains(&i) && (2..=3).contains(&j)), 1),
                    ],
                    copy,
                ));
            }
        }
        BSplineSurface::new(
            [3, 3],
            [knots.clone(), knots],
            [6, 6],
            controls,
            vec![Real::one(); 36],
            dec,
        )?
    } else {
        let u = MPoly::var(2, 0);
        let v = MPoly::var(2, 1);
        let bowl = &sq(&u) + &sq(&v);
        let f = match case {
            "bowl" | "rational-bowl" => bowl,
            "saddle" => &sq(&u) - &sq(&v),
            "flat" => &sq(&sq(&u)) + &sq(&sq(&v)),
            "line" => sq(&u),
            "flush" => MPoly::zero(2),
            "off-centre" | "steep-off-centre" => {
                &sq(&(&u - &constant(q(1, 3)))) + &sq(&(&v - &constant(q(2, 5))))
            }
            "several" => {
                &(&sq(&(&u + &constant(q(1, 2)))) + &sq(&v))
                    * &(&sq(&(&u - &constant(q(1, 2)))) + &sq(&v))
            }
            "algebraic" => &sq(&(&sq(&u) - &constant(q(1, 2)))) + &sq(&v),
            "near-miss" => &bowl + &constant(Rational::one().mul_pow2(-80)),
            "tiny-crossing" => &bowl - &constant(Rational::one().mul_pow2(-80)),
            "touch-and-cut" => &bowl * &(&u - &constant(q(1, 2))),
            "transverse" => &u + &v,
            "boundary" => &sq(&(&u - &constant(q(1, 1)))) + &sq(&v),
            _ => unreachable!(),
        };
        // Keep the graph inside the existing whole-patch chart certificate.
        // The steep copy remains in the round as a separate admission test.
        let f = if case == "steep-off-centre" {
            f
        } else {
            f.scale(&q(1, 16))
        };
        let w = if case == "rational-bowl" {
            &(&constant(q(2, 1)) + &u.scale(&q(1, 4))) + &v.scale(&q(1, 5))
        } else {
            constant(q(1, 1))
        };
        let mut homogeneous = [&u * &w, &v * &w, f, w];
        // Convert global [-1,1] powers to Bernstein controls on [0,1].
        for p in &mut homogeneous {
            for axis in 0..2 {
                *p = p.substitute(
                    axis,
                    &(&MPoly::var(2, axis).scale(&q(2, 1)) - &constant(q(1, 1))),
                );
            }
        }
        let degrees: [usize; 2] = std::array::from_fn(|axis| {
            homogeneous
                .iter()
                .flat_map(MPoly::terms)
                .map(|(m, _)| m.exponents()[axis] as usize)
                .max()
                .unwrap_or(1)
                .max(1)
        });
        let mut controls = Vec::new();
        let mut weights = Vec::new();
        for i in 0..=degrees[0] {
            for j in 0..=degrees[1] {
                let h: [Rational; 4] = homogeneous.each_ref().map(|p| {
                    p.terms()
                        .into_iter()
                        .filter_map(|(m, c)| {
                            let (a, b) = (m.exponents()[0] as usize, m.exponents()[1] as usize);
                            (a <= i && b <= j).then(|| {
                                c * choose(i, a) / choose(degrees[0], a) * choose(j, b)
                                    / choose(degrees[1], b)
                            })
                        })
                        .fold(Rational::zero(), |a, b| a + b)
                });
                let point = [0, 1, 2].map(|axis| &h[axis] / &h[3]);
                controls.push(transformed_point(&point, copy));
                weights.push(Real::from_rational(h[3].clone()));
            }
        }
        let knots = degrees.map(|d| {
            std::iter::repeat_n(q(-1, 1), d + 1)
                .chain(std::iter::repeat_n(q(1, 1), d + 1))
                .collect()
        });
        BSplineSurface::new(
            degrees,
            knots,
            degrees.map(|d| d + 1),
            controls,
            weights,
            dec,
        )?
    };
    let [nx, ny, nz] = transform(&[q(0, 1), q(0, 1), q(1, 1)], copy, true).map(Real::from_rational);
    let [rx, ry, rz] = transform(&[q(1, 1), q(0, 1), q(0, 1)], copy, true).map(Real::from_rational);
    let plane = Plane {
        origin: transformed_point(&[q(0, 1), q(0, 1), height], copy),
        normal: Vector3::new(nx, ny, nz),
        ref_dir: Vector3::new(rx, ry, rz),
    };
    Ok((RationalPatch::new(carrier, dec)?, plane))
}
