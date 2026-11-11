//! Quartic Bernstein substitution, keeping the torus equation factored.

use super::*;

pub(super) fn torus_net(
    torus: &Torus,
    controls: &[[Real; 4]],
    degrees: [usize; 2],
    dec: &mut Decider,
) -> Result<(Vec<Real>, [usize; 2]), Refusal> {
    let origin = [&torus.center.x, &torus.center.y, &torus.center.z];
    let axis = [&torus.axis.x, &torus.axis.y, &torus.axis.z];
    let mut local = Vec::with_capacity(controls.len());
    for h in controls {
        let mut d = h.clone();
        for i in 0..3 {
            d[i] = exact(&h[i], &(-exact(origin[i], &h[3], true, dec)?), false, dec)?;
        }
        local.push(d);
    }
    let major2 = exact(&torus.major, &torus.major, true, dec)?;
    let minor2 = exact(&torus.minor, &torus.minor, true, dec)?;
    let delta = exact(&major2, &(-minor2), false, dec)?;
    let mut a = (0..3).map(|i| ([i, i], Real::one())).collect::<Vec<_>>();
    a.push(([3, 3], delta));
    let (a, count) = quadric_net(&a, &local, degrees, dec)?;
    let norm = torus.axis.norm_squared();
    let mut radial = Vec::new();
    for i in 0..3 {
        for j in i..3 {
            let numerator = exact(axis[i], axis[j], true, dec)?;
            if let (Some(a), Some(b)) = (numerator.rational_ref(), norm.rational_ref()) {
                dec.reserve_exact_arithmetic(a, b, true)?;
            }
            let coefficient = -(&numerator / &norm);
            let coefficient = if i == j {
                exact(&Real::one(), &coefficient, false, dec)?
            } else {
                exact(&Real::from_i64(2), &coefficient, true, dec)?
            };
            radial.push(([i, j], coefficient));
        }
    }
    let (radial, _) = quadric_net(&radial, &local, degrees, dec)?;
    let (w2, _) = quadric_net(&[([3, 3], Real::one())], &local, degrees, dec)?;
    let (mut result, count4) = product(&a, &a, count, dec)?;
    let (radial_w2, _) = product(&radial, &w2, count, dec)?;
    let factor = exact(&Real::from_i64(4), &major2, true, dec)?;
    for (value, radial) in result.iter_mut().zip(radial_w2) {
        *value = exact(value, &(-exact(&factor, &radial, true, dec)?), false, dec)?;
    }
    Ok((result, count4))
}

fn binomial(n: usize, k: usize) -> Rational {
    let mut value = Rational::one();
    for i in 0..k {
        value = value * Rational::from_i64((n - i) as i64) / Rational::from_i64((i + 1) as i64);
    }
    value
}

/// Equal-degree tensor Bernstein product. The binomial ratios convert the
/// product basis to the Bernstein basis of twice the degree.
fn product(
    a: &[Real],
    b: &[Real],
    count: [usize; 2],
    dec: &mut Decider,
) -> Result<(Vec<Real>, [usize; 2]), Refusal> {
    let [m, n] = count.map(|c| c - 1);
    let result_count = [2 * m + 1, 2 * n + 1];
    let ratios = |degree: usize| {
        (0..=degree)
            .map(|i| {
                (0..=degree)
                    .map(|j| {
                        binomial(degree, i) * binomial(degree, j) / binomial(2 * degree, i + j)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let ru = ratios(m);
    let rv = ratios(n);
    let mut result = vec![Real::zero(); result_count[0] * result_count[1]];
    for i in 0..=m {
        dec.check_deadline()?;
        for j in 0..=n {
            for k in 0..=m {
                for l in 0..=n {
                    let weight = exact(
                        &ru[i][k].clone().into(),
                        &rv[j][l].clone().into(),
                        true,
                        dec,
                    )?;
                    let pair = exact(&a[i * count[1] + j], &b[k * count[1] + l], true, dec)?;
                    let term = exact(&weight, &pair, true, dec)?;
                    let index = (i + k) * result_count[1] + j + l;
                    result[index] = exact(&result[index], &term, false, dec)?;
                }
            }
        }
    }
    Ok((result, result_count))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn basis(axis: usize, degree: usize, index: usize) -> RealPoly {
        let x = RealPoly::var(2, axis);
        let complement = &RealPoly::constant(2, Real::one()) - &x;
        let mut p = RealPoly::constant(2, Real::from_rational(binomial(degree, index)));
        for _ in 0..index {
            p = &p * &x;
        }
        for _ in index..degree {
            p = &p * &complement;
        }
        p
    }

    #[test]
    fn quartic_net_and_span_equation_equal_the_homogeneous_implicit() {
        let mut dec = Decider::default();
        let controls = (0..3)
            .flat_map(|i| {
                (0..2).map(move |j| {
                    [
                        Real::from_i64(i - 1),
                        Real::from_i64(j + i),
                        Real::from_i64(i * j - 2),
                        Real::from_rational(Rational::from_ratio(3 + i + j, 3).unwrap()),
                    ]
                })
            })
            .collect::<Vec<_>>();
        let coordinates = std::array::from_fn(|axis| {
            (0..3).fold(RealPoly::zero(2), |sum, i| {
                sum + (0..2).fold(RealPoly::zero(2), |sum, j| {
                    sum + (basis(0, 2, i) * basis(1, 1, j)).scale(&controls[i * 2 + j][axis])
                })
            })
        });
        let span = Span::new(
            [0, 0],
            [
                [Rational::zero(), Rational::one()],
                [Rational::zero(), Rational::one()],
            ],
            coordinates,
        );
        let surface = SectionSurface::Torus(Torus {
            center: Point3::new(Real::from_i64(2), Real::from_i64(-3), Real::one()),
            axis: Vector3::new(Real::one(), Real::from_i64(2), Real::from_i64(-2)),
            major: Real::from_i64(3),
            minor: Real::from_rational(Rational::from_ratio(1, 2).unwrap()),
            ref_dir: Vector3::new(Real::from_i64(2), Real::from_i64(-1), Real::zero()),
        });
        let mut implicit = RealPoly::zero(2);
        for (monomial, coefficient) in surface.surface().implicit().unwrap().terms() {
            let mut term = RealPoly::constant(2, coefficient);
            let mut degree = 0;
            for (axis, &power) in monomial.exponents().iter().enumerate() {
                for _ in 0..power {
                    term = &term * &span.coordinates[axis];
                    degree += 1;
                }
            }
            for _ in degree..4 {
                term = &term * &span.coordinates[3];
            }
            implicit = implicit + term;
        }
        let (net, count) = surface.net(&controls, [2, 1], &mut dec).unwrap();
        assert_eq!(count, [9, 5]);
        let mut polynomial = RealPoly::zero(2);
        for i in 0..count[0] {
            for j in 0..count[1] {
                polynomial = polynomial
                    + (basis(0, count[0] - 1, i) * basis(1, count[1] - 1, j))
                        .scale(&net[i * count[1] + j]);
            }
        }
        assert!(
            (&surface.equation(&span) - &implicit)
                .is_zero(&mut dec)
                .unwrap()
        );
        assert!((&polynomial - &implicit).is_zero(&mut dec).unwrap());
    }
}
