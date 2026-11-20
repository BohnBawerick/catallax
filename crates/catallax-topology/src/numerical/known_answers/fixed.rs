//! The fourteen fixed known-answer edges of roadmap item 3.0A.
//!
//! In every case the two copies agree exactly except for one scalar
//! displacement along a fixed direction, so the same-parameter gap is that
//! scalar and its maximum has a closed form or an exact bisection enclosure.

use super::exact::*;
use super::*;
use catallax_shapes::SplineBoundary;

fn zero() -> Rational {
    Rational::zero()
}
fn one() -> Rational {
    Rational::one()
}
fn unit_knots(degree: usize, interior: &[Rational]) -> Vec<Rational> {
    clamped(degree, interior, &zero(), &one())
}
fn xy_plane() -> CaseSurface {
    CaseSurface::Kernel(plane([zero(), zero(), zero()], [0, 0, 1], [1, 0, 0]))
}

fn case(
    name: &str,
    source: &'static str,
    curve: Curve,
    surface: CaseSurface,
    face_copy: FaceCopy,
    truth: Truth,
) -> EdgeCase {
    EdgeCase {
        name: name.into(),
        source,
        geometry: EdgeGeometry {
            curve,
            parameters: [zero(), one()],
            surface,
            face_copy,
        },
        truth,
        domain_exit: None,
        zero: None,
        zero_ceiling: zero(),
    }
}

/// A strip `x = length u`, `y = width v`, `z = heights(u)` with a straight
/// v direction, and the straight edge along `v = 1/2` at height zero.
fn strip(
    name: &str,
    source: &'static str,
    degree: usize,
    knots: Vec<Rational>,
    heights: &[Rational],
    truth: Truth,
) -> EdgeCase {
    let (length, width) = (qi(100), qi(10));
    let controls = greville(&knots, degree)
        .iter()
        .zip(heights)
        .flat_map(|(u, z)| [zero(), width.clone()].map(|y| [&length * u, y, z.clone()]))
        .collect();
    let surface = carrier([degree, 1], [knots, unit_knots(1, &[])], controls);
    case(
        name,
        source,
        segment([zero(), ratio(5, 1), zero()], [length, ratio(5, 1), zero()]),
        CaseSurface::Carrier(surface),
        FaceCopy::Native(straight([zero(), ratio(1, 2)], [one(), ratio(1, 2)])),
        truth,
    )
}

/// Edge 1. A degree-1 curve with 130 interior knots carries a 1e-6 peak that
/// is 2^-30 wide and sits just past 5/16. OpenCASCADE fell back to 101 even
/// points above 101 knots (0028030) and lost this peak until PR #1640.
fn narrow_peak_between_eighths() -> EdgeCase {
    let peak = ratio(5, 16) + pow2(-30);
    let mut interior: Vec<_> = (1..128).map(|k| ratio(k, 128)).collect();
    interior.extend([&peak - &pow2(-31), peak.clone(), &peak + &pow2(-31)]);
    interior.sort();
    let knots = unit_knots(1, &interior);
    let height = q("1e-6");
    let controls = greville(&knots, 1)
        .into_iter()
        .map(|t| {
            let z = if t == peak { height.clone() } else { zero() };
            [qi(100) * t.clone(), zero(), z]
        })
        .collect();
    case(
        "01-narrow-peak-dense-knots",
        "OCCT 0028030 and PR #1640 DenseNonuniformKnotsPreserveNarrowMaximum",
        Curve::BSpline(spline(1, knots, controls, None)),
        xy_plane(),
        FaceCopy::Native(straight([zero(), zero()], [qi(100), zero()])),
        Truth::exact(height, peak),
    )
}

/// Edge 2. Both curve copies have one span. The peak is one basis function of
/// the surface, 2^-18 wide, so only the surface's knot lines reveal it.
fn peak_set_by_surface_knots() -> EdgeCase {
    let (centre, step) = (ratio(419, 1000), pow2(-20));
    let interior: Vec<_> = (-2..=2).map(|k| &centre + &(&step * &qi(k))).collect();
    let knots = unit_knots(3, &interior);
    let mut heights = vec![zero(); knots.len() - 4];
    heights[4] = q("3e-4");
    let peak = spline_peak(&knots, 3, &heights);
    strip(
        "02-surface-knot-peak",
        "OCCT PR #1640 SurfaceKnotsPreserveNarrowMaximum",
        3,
        knots,
        &heights,
        peak,
    )
}

/// Edge 3. One Bernstein function of a degree-25 surface peaks at 1/25.
fn narrow_peak_on_high_degree() -> EdgeCase {
    let mut heights = vec![zero(); 26];
    heights[1] = q("1e-3");
    strip(
        "03-degree-25-peak",
        "OCCT PR #1640 SurfaceDegreePreservesInteriorMaximum",
        25,
        unit_knots(25, &[]),
        &heights,
        Truth::exact(&heights[1] * &bernstein_peak(25, 1), ratio(1, 25)),
    )
}

/// Edge 4a. A rational cubic with end weight 2. Both copies share weights and
/// in-plane controls; one control of the 3D copy is lifted by `c`, so the gap
/// is `c` times that control's rational basis function,
/// `3 t (1 - t)^2 / (1 + (1 - t)^3)`.
fn rational_end_weight() -> EdgeCase {
    let lift = q("1e-4");
    let weights = vec![qi(2), one(), one(), one()];
    let flat = [[0, 0], [30, 20], [60, -10], [90, 0]].map(|p| p.map(qi));
    let controls = flat
        .iter()
        .enumerate()
        .map(|(i, [x, y])| {
            [
                x.clone(),
                y.clone(),
                if i == 1 { lift.clone() } else { zero() },
            ]
        })
        .collect();
    // With s = 1 - t the derivative's sign is that of 2 - 3 s - s^3.
    let value = |t: &Rational| {
        let s = one() - t.clone();
        &lift * &(qi(3) * t.clone() * s.clone() * s.clone() / (one() + &s * &s * s.clone()))
    };
    let slope = |t: &Rational| {
        let s = one() - t.clone();
        -(qi(2) - qi(3) * s.clone() - &s * &s * s.clone())
    };
    let truth = single_peak(&value, &slope, zero(), one(), &(qi(10) * lift.clone()));
    case(
        "04a-rational-end-weight",
        "OCCT PR #1640 CurveResidualPreservesEndSpanContinuation",
        Curve::BSpline(spline(
            3,
            unit_knots(3, &[]),
            controls,
            Some(weights.clone()),
        )),
        xy_plane(),
        FaceCopy::Native(pcurve(3, unit_knots(3, &[]), flat.into(), Some(weights))),
        truth,
    )
}

/// Edge 4b. A rational rail against a polynomial cubic face copy, as on a
/// numerical fillet. The rail is the cubic times the weight function
/// `(1, 4/5, 1)` over that weight, plus a lift on the weighted middle basis
/// function, whose value at 1/2 is 4/9.
fn rational_rail_against_cubic() -> EdgeCase {
    let lift = q("9e-5");
    let flat = [[0, 0], [20, 15], [50, 15], [70, 0]].map(|p| p.map(qi));
    let weight = [one(), ratio(4, 5), one()];
    let numerator = |axis: usize| {
        let cubic: Vec<_> = flat.iter().map(|p| p[axis].clone()).collect();
        bernstein_product(&cubic, &weight)
    };
    let (x, y) = (numerator(0), numerator(1));
    let z = bernstein_elevated(&[zero(), &lift * &weight[1], zero()], 5);
    let weights = bernstein_elevated(&weight, 5);
    let controls = (0..6)
        .map(|i| {
            [
                &x[i] / &weights[i],
                &y[i] / &weights[i],
                &z[i] / &weights[i],
            ]
        })
        .collect();
    let mut edge = case(
        "04b-rational-rail-against-cubic",
        "Catallax numerical fillet rails (numerical_carrier.rs) against cubic support curves",
        Curve::BSpline(spline(5, unit_knots(5, &[]), controls, Some(weights))),
        xy_plane(),
        FaceCopy::Native(pcurve(3, unit_knots(3, &[]), flat.into(), None)),
        Truth::exact(&lift * &ratio(4, 9), ratio(1, 2)),
    );
    edge.geometry.parameters = [zero(), ratio(7, 10)];
    edge
}

/// A cubic on the plane through `origin`, lifted along the normal by one
/// basis function of a non-uniform knot vector.
fn lifted_cubic(origin: &Rational, lift: &Rational) -> (Curve, CaseSurface, FaceCopy, Truth) {
    let interior = ["11/100", "23/100", "41/100", "62/100", "71/100", "88/100"].map(q);
    let knots = unit_knots(3, &interior);
    let count = knots.len() - 4;
    let mut heights = vec![zero(); count];
    heights[5] = lift.clone();
    let flat: Vec<[Rational; 2]> = greville(&knots, 3)
        .iter()
        .enumerate()
        .map(|(i, u)| [qi(50) * u.clone(), qi((i * i % 5) as i64)])
        .collect();
    let controls = flat
        .iter()
        .zip(&heights)
        .map(|([u, v], z)| [origin + u, origin + v, origin + z])
        .collect();
    (
        Curve::BSpline(spline(3, knots.clone(), controls, None)),
        CaseSurface::Kernel(plane(
            [origin.clone(), origin.clone(), origin.clone()],
            [0, 0, 1],
            [1, 0, 0],
        )),
        FaceCopy::Native(pcurve(3, knots.clone(), flat, None)),
        spline_peak(&knots, 3, &heights),
    )
}

/// Edge 5. The same small gap at coordinates of 1e6, 1e9 and 1e12 mm, where
/// one binary64 unit in the last place is up to 1.2e-4 mm.
fn far_origin(exponent: u32) -> EdgeCase {
    let origin = qi(10).pow(exponent as i32).expect("power");
    let (curve, surface, face_copy, truth) = lifted_cubic(&origin, &q("1.5e-5"));
    case(
        &format!("05-origin-1e{exponent}"),
        "OCCT PR #1640 (origin at 1e12) and the CompensatedReal residual",
        curve,
        surface,
        face_copy,
        truth,
    )
}

/// Edge 6. A parabola over `[0, 1/8]` peaks exactly midway between two equal
/// values. A lower, wider parabola over the rest of the edge tempts a checker
/// that refines only its best sample.
fn peak_midway_between_equal_samples() -> EdgeCase {
    let peak = q("2e-5");
    let knots = unit_knots(2, &[ratio(1, 8), ratio(1, 8)]);
    let heights = [
        zero(),
        qi(2) * peak.clone(),
        zero(),
        ratio(9, 5) * peak.clone(),
        zero(),
    ];
    let controls = greville(&knots, 2)
        .iter()
        .zip(&heights)
        .map(|(t, z)| [qi(100) * t.clone(), zero(), z.clone()])
        .collect();
    case(
        "06-peak-midway-between-equal-samples",
        "OCCT PR #1683 (tie between two equal samples)",
        Curve::BSpline(spline(2, knots, controls, None)),
        xy_plane(),
        FaceCopy::Native(straight([zero(), zero()], [qi(100), zero()])),
        Truth::exact(peak, ratio(1, 16)),
    )
}

/// Edge 7. A circle of radius 10 on its cylinder. The face copy keeps the
/// angle and lifts the height by one basis function. `shift` adds one stored
/// period to the face copy's angle.
fn seam(name: &str, start: Rational, end: Rational, shift: &Rational) -> EdgeCase {
    let (radius, height, lift) = (qi(10), qi(3), q("3e-5"));
    let span = &end - &start;
    let interior: Vec<_> = ["2/17", "5/17", "6/17", "9/17", "11/17", "15/17"]
        .iter()
        .map(|share| &start + &(&span * &q(share)))
        .collect();
    let knots = clamped(3, &interior, &start, &end);
    let mut lifts = vec![zero(); knots.len() - 4];
    lifts[4] = lift;
    let controls = greville(&knots, 3)
        .iter()
        .zip(&lifts)
        .map(|(u, dz)| [u + shift, &height + dz])
        .collect();
    let mut truth = spline_peak(&knots, 3, &lifts);
    // The stored period is a rational close to 2 pi, so a shifted copy also
    // differs along the circle by at most the radius times that difference.
    let pi = Ball::pi(BITS);
    let periods = (shift / &(pi.mid() * qi(2))).abs();
    truth.hi = &truth.hi + &(&radius * &periods * (pi.hi() - pi.lo()) * qi(2).clone());
    let mut edge = case(
        name,
        "OCCT 0025819 and 0026918 (periodic carriers); Catallax analytic.rs closed circles",
        circle(
            [zero(), zero(), height.clone()],
            [0, 0, 1],
            &radius,
            [1, 0, 0],
        ),
        CaseSurface::Kernel(cylinder(
            [zero(), zero(), zero()],
            [0, 0, 1],
            &radius,
            [1, 0, 0],
        )),
        FaceCopy::Native(pcurve(3, knots, controls, None)),
        truth,
    );
    edge.geometry.parameters = [start, end];
    edge
}

/// Edge 8. The face copy traces the same straight line with its parameter
/// stretched by 5.7e-6: `u = L (t + e (27/4) t^2 (1 - t))`. The geometric gap
/// is zero and the same-parameter gap is `e L` at `t = 2/3`.
fn parameter_stretch() -> EdgeCase {
    let (length, stretch) = (qi(10), q("5.7e-6"));
    let shares = [
        zero(),
        ratio(1, 3),
        ratio(2, 3) + ratio(9, 4) * stretch.clone(),
        one(),
    ];
    let controls = shares.iter().map(|s| [&length * s, zero()]).collect();
    case(
        "08-parameter-stretch",
        "solidcore #570 (EditParameterization stretched a trim by 1.0000057)",
        segment([zero(), zero(), zero()], [length.clone(), zero(), zero()]),
        xy_plane(),
        FaceCopy::Native(pcurve(3, unit_knots(3, &[]), controls, None)),
        Truth::exact(&length * &stretch, ratio(2, 3)),
    )
}

/// Edge 9. The face copy runs to `u = 21/20` on a patch that ends at `u = 1`.
/// The patch height is `h u^3`, and so is its polynomial continuation, which
/// lifts the far end to `h (21/20)^3` above the straight 3D edge.
fn face_copy_outside_domain() -> EdgeCase {
    let height = q("2e-3");
    let reach = ratio(21, 20);
    let heights = [zero(), zero(), zero(), height.clone()];
    let mut edge = strip(
        "09-face-copy-outside-domain",
        "solidcraft PR #18166 (trim at u = 71.94 lifted through extrapolation)",
        3,
        unit_knots(3, &[]),
        &heights,
        Truth::exact(&height * &reach.pow(3).expect("power"), one()),
    );
    edge.geometry.curve = segment(
        [zero(), ratio(5, 1), zero()],
        [qi(100) * reach.clone(), ratio(5, 1), zero()],
    );
    edge.geometry.face_copy = FaceCopy::Native(straight(
        [zero(), ratio(1, 2)],
        [reach.clone(), ratio(1, 2)],
    ));
    edge.domain_exit = Some(one() / reach);
    edge
}

/// Edge 10. A traced section between the plane `z = 0` and a plane tilted by
/// 2 atan(1/40), about 2.86 degrees. The trace wanders sideways by one basis
/// function; the face copy on the tilted plane is its exact projection, so
/// the gap is the wander times the sine of the angle, 80/1601.
fn near_tangent_section() -> EdgeCase {
    let wander = q("1e-2");
    let (sine, cosine) = (ratio(80, 1601), ratio(1599, 1601));
    let interior: Vec<_> = (1..32).map(|k| ratio(k * k, 1024)).collect();
    let knots = unit_knots(3, &interior);
    let mut sideways = vec![zero(); knots.len() - 4];
    sideways[12] = wander;
    let abscissae = greville(&knots, 3);
    let controls = abscissae
        .iter()
        .zip(&sideways)
        .map(|(t, y)| [qi(80) * t.clone(), y.clone(), zero()])
        .collect();
    let projected = abscissae
        .iter()
        .zip(&sideways)
        .map(|(t, y)| [qi(80) * t.clone(), y * &cosine])
        .collect();
    let mut truth = spline_peak(&knots, 3, &sideways);
    truth.lo = &truth.lo * &sine;
    truth.hi = &truth.hi * &sine;
    case(
        "10-near-tangent-section",
        "solidcraft #13630 and Catallax numerical Boolean sections (partition.rs)",
        Curve::BSpline(spline(3, knots.clone(), controls, None)),
        CaseSurface::Kernel(plane([zero(), zero(), zero()], [0, -80, 1599], [1, 0, 0])),
        FaceCopy::Native(pcurve(3, knots, projected, None)),
        truth,
    )
}

/// Edge 11. The gap is a degree-9 polynomial with a root at every eighth, so
/// a fit judged at eighths passes while the gap between them does not. Each
/// of the eight lobes has one extremum (Rolle), found by exact bisection.
fn passes_at_eighths() -> EdgeCase {
    let roots: Vec<_> = (0..=8).map(|k| ratio(k, 8)).collect();
    let polynomial = catallax_numbers::UPoly::from_roots(&roots);
    let derivative = polynomial.derivative();
    let lipschitz = derivative
        .coefficients()
        .iter()
        .map(Rational::abs)
        .fold(zero(), |sum, c| sum + c);
    let truth = roots
        .windows(2)
        .map(|lobe| {
            let sign = polynomial.eval(&(&lobe[0] + &lobe[1]).mul_pow2(-1)).sign();
            let flip = if sign == Sign::Negative {
                qi(-1)
            } else {
                one()
            };
            single_peak(
                &|t| &flip * &polynomial.eval(t),
                &|t| &flip * &derivative.eval(t),
                lobe[0].clone(),
                lobe[1].clone(),
                &lipschitz,
            )
        })
        .max_by(|a, b| a.lo.cmp(&b.lo))
        .expect("eight lobes");
    let heights = bernstein_from_power(polynomial.coefficients());
    let controls = heights
        .iter()
        .enumerate()
        .map(|(i, z)| [ratio(100 * i as i64, 9), zero(), z.clone()])
        .collect();
    case(
        "11-passes-at-eighths",
        "Catallax numerical_embedding.rs CHECK_SHARES (fit and admission share eighths)",
        Curve::BSpline(spline(9, unit_knots(9, &[]), controls, None)),
        xy_plane(),
        FaceCopy::Native(straight([zero(), zero()], [qi(100), zero()])),
        truth,
    )
}

/// A wavy bicubic patch with `spans` uniform u spans.
fn wavy_patch(spans: i64) -> (BSplineSurface, Vec<Rational>) {
    let interior: Vec<_> = (1..spans).map(|k| ratio(k, spans)).collect();
    let knots = unit_knots(3, &interior);
    let controls = greville(&knots, 3)
        .iter()
        .enumerate()
        .flat_map(|(i, u)| {
            (0..4).map(move |j| {
                [
                    qi(100) * u.clone(),
                    ratio(10 * j as i64, 3),
                    ratio(((i * 7 + j * 13) % 11) as i64 - 5, 10),
                ]
            })
        })
        .collect();
    (
        carrier([3, 3], [knots.clone(), unit_knots(3, &[])], controls),
        knots,
    )
}

/// The exact iso-curve controls of a patch at `v`, in u order.
fn row_at(surface: &BSplineSurface, v: &Rational) -> Vec<[Rational; 3]> {
    let [rows, columns] = surface.control_count();
    let [_, degree] = surface.degrees();
    let knots = surface.knots()[1];
    let shares: Vec<_> = (0..columns).map(|j| basis(knots, degree, j, v)).collect();
    (0..rows)
        .map(|i| {
            std::array::from_fn(|axis| {
                (0..columns).fold(zero(), |sum, j| {
                    let p = &surface.controls()[i * columns + j];
                    let coordinate = [&p.x, &p.y, &p.z][axis].as_rational().expect("rational");
                    sum + &shares[j] * &coordinate
                })
            })
        })
        .collect()
}

/// Edge 12. A 256-span edge on a 256-span bicubic patch, for work cost. The
/// 3D copy is the patch's own iso-curve with two neighbouring controls lifted
/// unequally, so the peak sits between knots.
fn many_spans() -> EdgeCase {
    let (surface, knots) = wavy_patch(256);
    let v = ratio(2, 5);
    let mut controls = row_at(&surface, &v);
    let mut lifts = vec![zero(); controls.len()];
    lifts[100] = q("4e-5");
    lifts[101] = q("2e-5");
    for (control, lift) in controls.iter_mut().zip(&lifts) {
        control[2] = &control[2] + lift;
    }
    let flat = greville(&knots, 3)
        .into_iter()
        .map(|u| [u, v.clone()])
        .collect();
    case(
        "12-256-spans",
        "Catallax R3.99 rim strip (324 stations) and OCCT issue #1372 (check cost)",
        Curve::BSpline(spline(3, knots.clone(), controls, None)),
        CaseSurface::Carrier(surface),
        FaceCopy::Native(pcurve(3, knots.clone(), flat, None)),
        spline_peak(&knots, 3, &lifts),
    )
}

/// Edge 13. Zero-gap controls. A checker must not raise a false alarm, and
/// may call the gap exactly zero only where the face copy follows by rule.
fn zero_gap_controls() -> Vec<EdgeCase> {
    let control = |mut edge: EdgeCase, zero: Zero| {
        edge.zero = Some(zero);
        edge.zero_ceiling = q("1e-10");
        edge
    };
    let line = Curve::Line(catallax_shapes::Line {
        origin: point([qi(1), qi(2), qi(3)]),
        dir: vector([3, 4, 0]),
    });
    let mut on_plane = case(
        "13a-zero-line-on-plane",
        "zero-gap control: line on plane",
        line,
        CaseSurface::Kernel(plane([zero(), zero(), qi(3)], [0, 0, 1], [1, 0, 0])),
        FaceCopy::Linear(Box::new([[qi(1), qi(2)], [qi(16), qi(22)]])),
        Truth::exact(zero(), zero()),
    );
    on_plane.geometry.parameters = [zero(), qi(5)];
    let mut on_cylinder = case(
        "13b-zero-circle-on-cylinder",
        "zero-gap control: circle on cylinder",
        circle([zero(), zero(), qi(3)], [0, 0, 1], &qi(10), [1, 0, 0]),
        CaseSurface::Kernel(cylinder(
            [zero(), zero(), zero()],
            [0, 0, 1],
            &qi(10),
            [1, 0, 0],
        )),
        FaceCopy::Linear(Box::new([[ratio(1, 4), qi(3)], [qi(2), qi(3)]])),
        Truth::exact(zero(), ratio(1, 4)),
    );
    on_cylinder.geometry.parameters = [ratio(1, 4), qi(2)];
    let (surface, knots) = wavy_patch(5);
    let rail = case(
        "13c-zero-boundary-rail",
        "zero-gap control: boundary row of the carrier, as a fillet rail",
        Curve::BSpline(BSplineCurve::boundary(&surface, SplineBoundary::VMin)),
        CaseSurface::Carrier(surface.clone()),
        FaceCopy::Native(straight([zero(), zero()], [one(), zero()])),
        Truth::exact(zero(), zero()),
    );
    let v = ratio(2, 5);
    let stored = case(
        "13d-zero-stored-iso-curve",
        "zero-gap control: an interior iso-curve stored as its own curve",
        Curve::BSpline(spline(3, knots, row_at(&surface, &v), None)),
        CaseSurface::Carrier(surface),
        FaceCopy::Native(straight([zero(), v.clone()], [one(), v])),
        Truth::exact(zero(), zero()),
    );
    vec![
        control(on_plane, Zero::Structural),
        control(on_cylinder, Zero::Structural),
        control(rail, Zero::Structural),
        control(stored, Zero::Measured),
    ]
}

/// Edge 14. An exact rational quadratic copy of a circular arc of radius 5
/// against a face copy that is linear in the cylinder's angle. With
/// `tau = tan(alpha / 2) = 1/3` the rational copy sits at angle
/// `2 atan(tau x)` for `x = 2 t - 1`, the face copy at `a x`. The difference
/// is stationary at `x^2 = (2 tau / a - 1) / tau^2`, and the gap is the chord
/// `2 R sin(d / 2)`.
fn circle_against_rational_copy() -> EdgeCase {
    let (radius, tau) = (qi(5), ratio(1, 3));
    let ball = |value: &Rational| Ball::from_rational(value, BITS);
    let half_angle = (ball(&tau).atan() * ball(&qi(2))).mid();
    let x = ball(&(qi(2) * tau.clone() / half_angle.clone() - one()))
        .sqrt()
        .expect("positive")
        * ball(&(one() / tau.clone()));
    let turned = (&x * &ball(&tau)).atan() * ball(&qi(2));
    let difference = &turned - &(&x * &ball(&half_angle));
    let chord = (difference * ball(&ratio(1, 2))).sin() * ball(&(qi(2) * radius.clone()));
    let at = (x.mid() + one()).mul_pow2(-1);
    // The stationary point is rounded to a ball; its width moves the chord by
    // a second-order amount that the relative slack below covers.
    let slack = chord.hi() * pow2(-100);
    case(
        "14-circle-against-rational-copy",
        "FreeCAD via solidcraft #8502 (circular edges written as rational splines)",
        Curve::BSpline(spline(
            2,
            unit_knots(2, &[]),
            vec![
                [qi(4), qi(-3), zero()],
                [ratio(25, 4), zero(), zero()],
                [qi(4), qi(3), zero()],
            ],
            Some(vec![one(), ratio(4, 5), one()]),
        )),
        CaseSurface::Kernel(cylinder(
            [zero(), zero(), zero()],
            [0, 0, 1],
            &radius,
            [1, 0, 0],
        )),
        FaceCopy::Linear(Box::new([
            [-half_angle.clone(), zero()],
            [half_angle, zero()],
        ])),
        Truth {
            lo: chord.lo() - slack.clone(),
            hi: chord.hi() + slack,
            at,
        },
    )
}

/// The names of the fixed cases, in order. An acceptance run must contain
/// exactly these: one or more cases for each of the fourteen numbered edges.
pub(super) const REQUIRED: [&str; 22] = [
    "01-narrow-peak-dense-knots",
    "02-surface-knot-peak",
    "03-degree-25-peak",
    "04a-rational-end-weight",
    "04b-rational-rail-against-cubic",
    "05-origin-1e6",
    "05-origin-1e9",
    "05-origin-1e12",
    "06-peak-midway-between-equal-samples",
    "07a-seam-crossing",
    "07b-one-period-shift",
    "07c-closed-circle",
    "08-parameter-stretch",
    "09-face-copy-outside-domain",
    "10-near-tangent-section",
    "11-passes-at-eighths",
    "12-256-spans",
    "13a-zero-line-on-plane",
    "13b-zero-circle-on-cylinder",
    "13c-zero-boundary-rail",
    "13d-zero-stored-iso-curve",
    "14-circle-against-rational-copy",
];

/// All fixed cases in receipt order.
pub(super) fn cases() -> Vec<EdgeCase> {
    let period = Ball::pi(BITS).mid() * qi(2);
    let mut cases = vec![
        narrow_peak_between_eighths(),
        peak_set_by_surface_knots(),
        narrow_peak_on_high_degree(),
        rational_end_weight(),
        rational_rail_against_cubic(),
        far_origin(6),
        far_origin(9),
        far_origin(12),
        peak_midway_between_equal_samples(),
        seam("07a-seam-crossing", ratio(-1, 2), ratio(1, 2), &zero()),
        seam("07b-one-period-shift", ratio(-1, 2), ratio(1, 2), &period),
        seam("07c-closed-circle", zero(), period.clone(), &zero()),
        parameter_stretch(),
        face_copy_outside_domain(),
        near_tangent_section(),
        passes_at_eighths(),
        many_spans(),
    ];
    cases.extend(zero_gap_controls());
    cases.push(circle_against_rational_copy());
    cases
}

#[test]
fn every_fixed_truth_is_attained_and_bounds_dense_samples() {
    let cases = cases();
    let numbers: std::collections::BTreeSet<_> = cases
        .iter()
        .map(|case| case.name[..2].to_string())
        .collect();
    assert_eq!(numbers.len(), 14, "one or more cases per numbered edge");
    let names: Vec<&str> = cases.iter().map(|case| case.name.as_str()).collect();
    assert_eq!(names, REQUIRED);
    for case in &cases {
        assert!(!case.source.is_empty());
        assert_eq!(case.truth.is_zero(), case.zero.is_some(), "{}", case.name);
        truth_is_attained_and_bounds(case, 192);
    }
}

#[test]
fn the_domain_exit_truth_is_the_polynomial_continuation() {
    let case = face_copy_outside_domain();
    let exit = case.domain_exit.as_ref().expect("exit");
    // Up to the exit the kernel's own evaluator agrees with h u^3.
    let mut dec = Decider::default();
    let [lo, hi] = gap_enclosure(&case.geometry, exit, &mut dec).expect("inside");
    let at_exit = q("2e-3");
    assert!(lo <= at_exit && at_exit <= hi);
    assert_eq!(case.truth.lo, q("2e-3") * ratio(9261, 8000));
}
