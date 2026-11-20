//! Exact rational tools for building known-answer edges and their truths.

use super::*;
use catallax_shapes::{Point3, Vector3};

pub(super) fn q(text: &str) -> Rational {
    Rational::parse(text).expect("rational literal")
}
pub(super) fn qi(value: i64) -> Rational {
    Rational::from_i64(value)
}
pub(super) fn ratio(numerator: i64, denominator: i64) -> Rational {
    Rational::from_ratio(numerator, denominator).expect("nonzero denominator")
}
pub(super) fn pow2(exponent: i64) -> Rational {
    Rational::one().mul_pow2(exponent)
}
pub(super) fn point(p: [Rational; 3]) -> Point3 {
    let [x, y, z] = p;
    Point3::new(x.into(), y.into(), z.into())
}
pub(super) fn vector(v: [i64; 3]) -> Vector3 {
    let [x, y, z] = v.map(|c| Real::from(qi(c)));
    Vector3::new(x, y, z)
}

/// A clamped knot vector over `[lo, hi]` with the given interior knots.
pub(super) fn clamped(
    degree: usize,
    interior: &[Rational],
    lo: &Rational,
    hi: &Rational,
) -> Vec<Rational> {
    std::iter::repeat_n(lo.clone(), degree + 1)
        .chain(interior.iter().cloned())
        .chain(std::iter::repeat_n(hi.clone(), degree + 1))
        .collect()
}

/// Greville abscissae. Controls placed at these values reproduce the identity.
pub(super) fn greville(knots: &[Rational], degree: usize) -> Vec<Rational> {
    (0..knots.len() - degree - 1)
        .map(|i| {
            knots[i + 1..=i + degree]
                .iter()
                .fold(Rational::zero(), |sum, knot| sum + knot.clone())
                / qi(degree as i64)
        })
        .collect()
}

/// Cox-de Boor value of basis function `i`. The last knot belongs to the last
/// nonempty span, so a clamped spline is defined on its closed domain.
pub(super) fn basis(knots: &[Rational], degree: usize, i: usize, x: &Rational) -> Rational {
    if degree == 0 {
        let last = knots.last().expect("knots");
        let inside = &knots[i] <= x && x < &knots[i + 1];
        let closing = x == last && &knots[i + 1] == last && knots[i] < knots[i + 1];
        return if inside || closing {
            Rational::one()
        } else {
            Rational::zero()
        };
    }
    let mut value = Rational::zero();
    if knots[i + degree] > knots[i] {
        value = value
            + (x - &knots[i]) / (&knots[i + degree] - &knots[i]) * basis(knots, degree - 1, i, x);
    }
    if knots[i + degree + 1] > knots[i + 1] {
        value = value
            + (&knots[i + degree + 1] - x) / (&knots[i + degree + 1] - &knots[i + 1])
                * basis(knots, degree - 1, i + 1, x);
    }
    value
}

/// Exact value of the scalar spline with the given coefficients.
pub(super) fn spline_value(
    knots: &[Rational],
    degree: usize,
    coefficients: &[Rational],
    x: &Rational,
) -> Rational {
    coefficients
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .fold(Rational::zero(), |sum, (i, c)| {
            sum + c * &basis(knots, degree, i, x)
        })
}

/// Knots and coefficients of the derivative of a scalar spline.
pub(super) fn spline_derivative(
    knots: &[Rational],
    degree: usize,
    coefficients: &[Rational],
) -> (Vec<Rational>, Vec<Rational>) {
    let slopes = (0..coefficients.len() - 1)
        .map(|i| {
            let width = &knots[i + degree + 1] - &knots[i + 1];
            if width.is_zero() {
                Rational::zero()
            } else {
                qi(degree as i64) * (&coefficients[i + 1] - &coefficients[i]) / width
            }
        })
        .collect();
    (knots[1..knots.len() - 1].to_vec(), slopes)
}

/// The maximum of a function that rises and then falls on `[lo, hi]`, by
/// bisection on the sign of its exact derivative. `lipschitz` bounds the
/// derivative, which turns the final bracket into a rigorous enclosure.
pub(super) fn single_peak(
    value: &dyn Fn(&Rational) -> Rational,
    slope: &dyn Fn(&Rational) -> Rational,
    mut lo: Rational,
    mut hi: Rational,
    lipschitz: &Rational,
) -> Truth {
    for _ in 0..120 {
        let middle = (&lo + &hi).mul_pow2(-1);
        match slope(&middle).sign() {
            Sign::Positive => lo = middle,
            Sign::Negative => hi = middle,
            Sign::Zero => {
                return Truth::exact(value(&middle), middle);
            }
        }
    }
    let found = value(&lo).max(value(&hi));
    Truth {
        hi: &found + &(lipschitz * &(&hi - &lo)),
        lo: found,
        at: lo,
    }
}

/// The maximum of a nonnegative spline bump with one rise and one fall.
/// A B-spline basis function is such a bump (Curry and Schoenberg 1966: it is
/// a Polya frequency function), and so is a positive sum of neighbours whose
/// derivative coefficients change sign once (Schoenberg's variation
/// diminishing property).
pub(super) fn spline_peak(knots: &[Rational], degree: usize, coefficients: &[Rational]) -> Truth {
    let first = coefficients
        .iter()
        .position(|c| !c.is_zero())
        .expect("a bump");
    let last = coefficients
        .iter()
        .rposition(|c| !c.is_zero())
        .expect("a bump");
    let (slope_knots, slopes) = spline_derivative(knots, degree, coefficients);
    let signs: Vec<_> = slopes
        .iter()
        .map(Rational::sign)
        .filter(|sign| *sign != Sign::Zero)
        .collect();
    assert!(
        signs.windows(2).filter(|pair| pair[0] != pair[1]).count() == 1,
        "bump is not single-peaked"
    );
    let lipschitz = slopes
        .iter()
        .map(Rational::abs)
        .fold(Rational::zero(), Rational::max);
    single_peak(
        &|x| spline_value(knots, degree, coefficients, x),
        &|x| spline_value(&slope_knots, degree - 1, &slopes, x),
        knots[first].clone(),
        knots[last + degree + 1].clone(),
        &lipschitz,
    )
}

/// Boehm's knot insertion: the same curve on one more knot.
pub(super) fn insert_knot(
    degree: usize,
    knots: &mut Vec<Rational>,
    controls: &mut Vec<Vec<Rational>>,
    x: &Rational,
) {
    let span = knots
        .iter()
        .rposition(|knot| knot <= x)
        .expect("knot inside the domain")
        .min(controls.len() - 1);
    let mut next = Vec::with_capacity(controls.len() + 1);
    for i in 0..=controls.len() {
        if i + degree <= span {
            next.push(controls[i].clone());
        } else if i > span {
            next.push(controls[i - 1].clone());
        } else {
            let share = (x - &knots[i]) / (&knots[i + degree] - &knots[i]);
            let keep = Rational::one() - share.clone();
            next.push(
                controls[i - 1]
                    .iter()
                    .zip(&controls[i])
                    .map(|(a, b)| &keep * a + &share * b)
                    .collect(),
            );
        }
    }
    knots.insert(span + 1, x.clone());
    *controls = next;
}

pub(super) fn binomial(n: usize, k: usize) -> Rational {
    (0..k).fold(Rational::one(), |value, i| {
        value * qi((n - i) as i64) / qi(i as i64 + 1)
    })
}

/// The maximum of the Bernstein polynomial `B(k, n)`, attained at `k / n`.
pub(super) fn bernstein_peak(n: usize, k: usize) -> Rational {
    let at = ratio(k as i64, n as i64);
    binomial(n, k)
        * at.pow(k as i32).expect("power")
        * (Rational::one() - at.clone())
            .pow((n - k) as i32)
            .expect("power")
}

/// Bernstein coefficients of a product.
pub(super) fn bernstein_product(a: &[Rational], b: &[Rational]) -> Vec<Rational> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    (0..=m + n)
        .map(|k| {
            (k.saturating_sub(n)..=k.min(m)).fold(Rational::zero(), |sum, i| {
                sum + binomial(m, i) * binomial(n, k - i) / binomial(m + n, k)
                    * a[i].clone()
                    * b[k - i].clone()
            })
        })
        .collect()
}

/// Bernstein coefficients of the same polynomial at a higher degree.
pub(super) fn bernstein_elevated(a: &[Rational], degree: usize) -> Vec<Rational> {
    let mut current = a.to_vec();
    while current.len() <= degree {
        let n = current.len();
        current = (0..=n)
            .map(|i| {
                let share = ratio(i as i64, n as i64);
                let mut value = Rational::zero();
                if i > 0 {
                    value = value + &share * &current[i - 1];
                }
                if i < n {
                    value = value + (Rational::one() - share.clone()) * current[i].clone();
                }
                value
            })
            .collect();
    }
    current
}

/// Bernstein coefficients on `[0, 1]` from power-basis coefficients.
pub(super) fn bernstein_from_power(a: &[Rational]) -> Vec<Rational> {
    let n = a.len() - 1;
    (0..=n)
        .map(|i| {
            (0..=i).fold(Rational::zero(), |sum, k| {
                sum + binomial(i, k) / binomial(n, k) * a[k].clone()
            })
        })
        .collect()
}

pub(super) fn spline(
    degree: usize,
    knots: Vec<Rational>,
    controls: Vec<[Rational; 3]>,
    weights: Option<Vec<Rational>>,
) -> BSplineCurve {
    let weights = weights
        .map(|weights| weights.into_iter().map(Real::from).collect())
        .unwrap_or_else(|| vec![Real::one(); controls.len()]);
    BSplineCurve::new(
        degree,
        knots,
        controls.into_iter().map(point).collect(),
        weights,
        &mut Decider::default(),
    )
    .expect("valid fixture curve")
}

/// A u/v curve. The unused third coordinate is zero.
pub(super) fn pcurve(
    degree: usize,
    knots: Vec<Rational>,
    controls: Vec<[Rational; 2]>,
    weights: Option<Vec<Rational>>,
) -> BSplineCurve {
    spline(
        degree,
        knots,
        controls
            .into_iter()
            .map(|[u, v]| [u, v, Rational::zero()])
            .collect(),
        weights,
    )
}

/// A straight u/v curve over `[0, 1]`.
pub(super) fn straight(from: [Rational; 2], to: [Rational; 2]) -> BSplineCurve {
    pcurve(
        1,
        clamped(1, &[], &Rational::zero(), &Rational::one()),
        vec![from, to],
        None,
    )
}

/// A straight 3D spline over `[0, 1]`.
pub(super) fn segment(from: [Rational; 3], to: [Rational; 3]) -> Curve {
    Curve::BSpline(spline(
        1,
        clamped(1, &[], &Rational::zero(), &Rational::one()),
        vec![from, to],
        None,
    ))
}

/// A polynomial carrier in u-major control order.
pub(super) fn carrier(
    degrees: [usize; 2],
    knots: [Vec<Rational>; 2],
    controls: Vec<[Rational; 3]>,
) -> BSplineSurface {
    let count = [0, 1].map(|axis| knots[axis].len() - degrees[axis] - 1);
    let weights = vec![Real::one(); controls.len()];
    BSplineSurface::new(
        degrees,
        knots,
        count,
        controls.into_iter().map(point).collect(),
        weights,
        &mut Decider::default(),
    )
    .expect("valid fixture surface")
}

pub(super) fn plane(origin: [Rational; 3], normal: [i64; 3], ref_dir: [i64; 3]) -> Surface {
    Surface::Plane(catallax_shapes::Plane {
        origin: point(origin),
        normal: vector(normal),
        ref_dir: vector(ref_dir),
    })
}

pub(super) fn cylinder(
    origin: [Rational; 3],
    axis: [i64; 3],
    radius: &Rational,
    ref_dir: [i64; 3],
) -> Surface {
    Surface::Cylinder(catallax_shapes::Cylinder {
        origin: point(origin),
        axis: vector(axis),
        radius_sq: Real::from(radius * radius),
        ref_dir: vector(ref_dir),
    })
}

pub(super) fn circle(
    center: [Rational; 3],
    normal: [i64; 3],
    radius: &Rational,
    ref_dir: [i64; 3],
) -> Curve {
    Curve::Circle(catallax_shapes::Circle {
        center: point(center),
        normal: vector(normal),
        radius_sq: Real::from(radius * radius),
        ref_dir: vector(ref_dir),
    })
}

/// The face copy's exact u/v at an edge parameter.
pub(super) fn face_parameters(
    case: &EdgeCase,
    t: &Rational,
    dec: &mut Decider,
) -> Result<[Real; 2], Refusal> {
    Ok(match &case.face_copy {
        FaceCopy::Native(curve) => {
            let at = curve.point_at(&Real::from(t.clone()), dec)?;
            [at.x, at.y]
        }
        FaceCopy::Linear(ends) => {
            let [a, b] = &case.parameters;
            let share = (t - a) / (b - a);
            std::array::from_fn(|axis| {
                Real::from(&ends[0][axis] + &(&share * &(&ends[1][axis] - &ends[0][axis])))
            })
        }
    })
}

/// A 128-bit enclosure of the same-parameter gap at one rational parameter,
/// from the exact definitions. This is the oracle's only evaluator.
pub(super) fn gap_enclosure(
    case: &EdgeCase,
    t: &Rational,
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let p = curve_at(&case.curve, t, dec)?.expect("fixture curves have evaluators");
    let uv = face_parameters(case, t, dec)?;
    let q = match &case.surface {
        CaseSurface::Kernel(surface) => native::point(surface, &uv, dec)?,
        CaseSurface::Carrier(carrier) => {
            let u = dec.enclose_budgeted_at(&uv[0], BITS)?;
            let v = dec.enclose_budgeted_at(&uv[1], BITS)?;
            carrier.evaluate(&u, &v, dec)?
        }
    };
    let d = std::array::from_fn::<_, 3, _>(|i| &p[i] - &q[i]);
    let length = (d[0].square() + d[1].square() + d[2].square())
        .sqrt()
        .expect("sum of squares");
    Ok([length.lo().max(Rational::zero()), length.hi()])
}

#[test]
fn rational_tools_agree_with_their_definitions() {
    // Cardinal B-spline centre values 3/4 and 2/3 (Schoenberg 1946).
    for (degree, centre) in [(1, qi(1)), (2, ratio(3, 4)), (3, ratio(2, 3))] {
        let interior: Vec<_> = (1..=degree as i64 + 2).map(qi).collect();
        let knots = clamped(degree, &interior, &qi(0), &qi(degree as i64 + 3));
        let mut coefficients = vec![Rational::zero(); knots.len() - degree - 1];
        coefficients[degree + 1] = Rational::one();
        let peak = spline_peak(&knots, degree, &coefficients);
        assert!(peak.lo <= centre && centre <= peak.hi, "degree {degree}");
    }
    assert_eq!(bernstein_peak(2, 1), ratio(1, 2));
    assert_eq!(
        bernstein_product(&[qi(1), qi(0)], &[qi(0), qi(1)]),
        vec![qi(0), ratio(1, 2), qi(0)]
    );
    assert_eq!(
        bernstein_from_power(&[qi(0), qi(0), qi(1)]),
        vec![qi(0), qi(0), qi(1)]
    );
    assert_eq!(
        bernstein_elevated(&[qi(0), qi(1)], 3),
        vec![qi(0), ratio(1, 3), ratio(2, 3), qi(1)]
    );
    // Knot insertion keeps the function.
    let mut knots = clamped(3, &[ratio(1, 3), ratio(3, 5)], &qi(0), &qi(1));
    let mut controls: Vec<Vec<Rational>> = [3, -1, 4, 1, -5, 9].map(|c| vec![qi(c)]).into();
    let before = (knots.clone(), controls.clone());
    insert_knot(3, &mut knots, &mut controls, &ratio(1, 2));
    insert_knot(3, &mut knots, &mut controls, &ratio(1, 2));
    for i in 0..=20 {
        let x = ratio(i, 20);
        let old: Vec<_> = before.1.iter().map(|c| c[0].clone()).collect();
        let new: Vec<_> = controls.iter().map(|c| c[0].clone()).collect();
        assert_eq!(
            spline_value(&before.0, 3, &old, &x),
            spline_value(&knots, 3, &new, &x)
        );
    }
}
