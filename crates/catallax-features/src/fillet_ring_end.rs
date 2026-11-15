//! The equal-radius torus cut where a mixed corner leaves its upright sharp.
//!
//! In units of the fillet radius, the end wall is one unit from the ring
//! axis and the torus has major radius two and minor radius one. In that
//! plane its equation is `(X²+Z²)²+8(X²-Z²)=0`. The branch from the planar
//! rail to the cylinder contact is the rational quartic
//! `X=4√2 u(1-u²)/D`, `Z=2√2(1-u⁴)/D`, `D=1+6u²+u⁴`,
//! for `√3-√2 <= u <= 1`. Substitution proves the complete curve lies on
//! both supports, including the singular carrier endpoint at X=Z=0.

use super::*;
use catallax_shapes::BSplineCurve;

/// Keep the exact pipe-plane section as the single boundary of both faces.
pub(super) fn pipe_geometry(
    model: &Model,
    run: &Run,
    wall: FaceId,
    old: VertexId,
    dec: &mut Decider,
) -> Result<Option<([Point3; 2], Curve)>, Refusal> {
    use catallax_shapes::{RationalAtlas, SectionSurface, spline_tools::RationalSplineSurface};

    let Surface::BSpline(SplineSupport::RollingPipe(pipe)) = &run.definition.surface else {
        return Ok(None);
    };
    let unsupported = || {
        Refusal::new(
            UNSUPPORTED_END,
            "a pipe end needs a regular vertical-wall section",
        )
    };
    let Surface::Plane(plane) = model.geom().surface(model.face(wall).surface) else {
        return Err(unsupported());
    };
    let law = pipe.law();
    if pipe.is_conical()
        || !run.definition.convex
        || sign(&plane.normal.cross(&law.across).norm_squared(), dec)? != Sign::Zero
        || sign(&law.base.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
    {
        return Err(unsupported());
    }
    let at = model.geom().position(model.vertex(old).point);
    let side = sign(&at.minus(&law.base).dot(&law.side), dec)?;
    let centre = match side {
        Sign::Positive => 1,
        Sign::Negative => 3,
        Sign::Zero => return Err(unsupported()),
    };
    let source = pipe.carrier();
    let mut rational = |value: &Real| dec.rational_value(value)?.ok_or_else(unsupported);
    let controls = source
        .controls()
        .iter()
        .map(|p| Ok([rational(&p.x)?, rational(&p.y)?, rational(&p.z)?]))
        .collect::<Result<Vec<_>, Refusal>>()?;
    let weights = source
        .weights()
        .iter()
        .map(&mut rational)
        .collect::<Result<Vec<_>, _>>()?;
    let source = RationalSplineSurface::new(
        source.degrees(),
        source.knots().map(<[_]>::to_vec),
        source.control_count(),
        controls,
        weights,
        dec,
    )?;
    let low = Rational::from_ratio(4 * centre - 1, 4).unwrap();
    let high = Rational::from_ratio(4 * centre + 1, 4).unwrap();
    let zero = Rational::zero();
    let one = Rational::one();
    let source = source.extracted([[&low, &high], [&zero, &one]], dec)?;
    // Local knots avoid cancellation in the retained section's power basis.
    let carrier = RationalSplineSurface::new(
        source.degrees(),
        [
            source.knots()[0]
                .iter()
                .map(|knot| (knot - &low) / (&high - &low))
                .collect(),
            source.knots()[1].to_vec(),
        ],
        source.control_count(),
        source.controls().to_vec(),
        source.weights().to_vec(),
        dec,
    )?
    .to_surface(dec)?;
    let atlases = RationalAtlas::partition_with_polynomial_bounds(carrier, [false; 2], dec)?;
    let [atlas] = atlases.as_slice() else {
        return Err(unsupported());
    };
    let section = atlas.sections(&SectionSurface::Plane(plane.clone()), dec)?;
    let [branch] = section.branches() else {
        return Err(unsupported());
    };
    let mut branch = branch.clone();
    let parameters = branch.parameter_endpoints();
    if sign(&parameters[0][1], dec)? != Sign::Zero {
        branch = branch.reversed();
    }
    let parameters = branch.parameter_endpoints();
    if branch.is_closed()
        || sign(&parameters[0][1], dec)? != Sign::Zero
        || sign(&(&parameters[1][1] - &Real::one()), dec)? != Sign::Zero
    {
        return Err(unsupported());
    }
    let radial_sq = &law.offset.square() - &(&law.radius * &law.tilt).square();
    if sign(&radial_sq, dec)? != Sign::Positive {
        return Err(unsupported());
    }
    let side = Real::from_i64(if side == Sign::Positive { 1 } else { -1 });
    // The plane rail has axial offset reach*r/along; the wall rail has none.
    // These are the same section ends without cancelled chart-root expressions.
    let mut feet = [
        law.base
            .plus(&law.axis.scaled(&(&(&law.reach * &law.radius) / &law.along)))
            .plus(&law.side.scaled(&(&side * &radial_sq.sqrt()))),
        law.base.plus(&law.side.scaled(&(&side * &law.body))),
    ];
    if !matches!(run.definition.supports[0], Surface::Plane(_)) {
        feet.swap(0, 1);
        branch = branch.reversed();
    }
    Ok(Some((feet, Curve::PatchSection(branch))))
}

pub(super) struct End {
    pub(super) feet: [Point3; 2],
    pub(super) curve: Curve,
    // The extended end direction, expressed on the input edge's circle.
    pub(super) hull_end: Point3,
}

fn choose(n: usize, k: usize) -> i64 {
    if k > n {
        return 0;
    }
    (0..k).fold(1, |a, i| a * (n - i) as i64 / (i + 1) as i64)
}

// Power coefficients in u, transformed to Bernstein coefficients on [lo,1].
fn bernstein(power: [Real; 5], lo: &Real) -> [Real; 5] {
    let width = &Real::one() - lo;
    let mut shifted: [Real; 5] = std::array::from_fn(|_| Real::zero());
    for (i, coefficient) in power.iter().enumerate() {
        for (j, value) in shifted.iter_mut().enumerate().take(i + 1) {
            let mut term = coefficient * &Real::from_i64(choose(i, j));
            for _ in 0..i - j {
                term = &term * lo;
            }
            for _ in 0..j {
                term = &term * &width;
            }
            *value = &*value + &term;
        }
    }
    std::array::from_fn(|i| {
        (0..=i).fold(Real::zero(), |sum, j| {
            &sum + &(&shifted[j]
                * &Real::from_rational(
                    Rational::from_ratio(choose(i, j), choose(4, j)).expect("nonzero binomial"),
                ))
        })
    })
}

/// A partial mixed end has a planar support normal to the ring axis, a
/// concave cylinder of the same radius, and an end wall tangent to that
/// cylinder. Other ring ends keep their existing construction/refusal.
pub(super) fn geometry(
    model: &Model,
    run: &Run,
    wall: FaceId,
    old: VertexId,
    dec: &mut Decider,
) -> Result<Option<End>, Refusal> {
    let Surface::Torus(torus) = &run.definition.surface else {
        return Ok(None);
    };
    if !run.definition.convex {
        return Ok(None);
    }
    let Some(planar) = run
        .definition
        .supports
        .iter()
        .position(|s| matches!(s, Surface::Plane(_)))
    else {
        return Ok(None);
    };
    let circular = 1 - planar;
    let (Surface::Cylinder(cylinder), Surface::Plane(end), Curve::Circle(input)) = (
        &run.definition.supports[circular],
        model.geom().surface(model.face(wall).surface),
        &run.definition.edge,
    ) else {
        return Ok(None);
    };
    let radius = &run.definition.radius;
    if sign(&(&torus.major - &(radius * &Real::from_i64(2))), dec)? != Sign::Zero
        || sign(&(&torus.minor - radius), dec)? != Sign::Zero
        || sign(&(&cylinder.radius_sq - &radius.square()), dec)? != Sign::Zero
        || sign(&end.normal.dot(&torus.axis), dec)? != Sign::Zero
        || !tangent_faces(
            model,
            run.faces[circular],
            wall,
            &model.geom().position(model.vertex(old).point),
            dec,
        )?
    {
        return Ok(None);
    }
    let e = model.edge(run.edge);
    let far = if e.start == Some(old) { e.end } else { e.start };
    let Some(far) = far else {
        return Ok(None);
    };
    let at = model.geom().position(model.vertex(old).point);
    let radial = at.minus(&input.center);
    if sign(&(&radial.norm_squared() - &radius.square()), dec)? != Sign::Zero {
        return Ok(None);
    }
    let q = radial.scaled(&(&Real::one() / radius));
    let across = input.center.minus(&torus.center);
    if sign(&(&across.norm_squared() - &radius.square()), dec)? != Sign::Zero
        || sign(&across.cross(&torus.axis).norm_squared(), dec)? != Sign::Zero
    {
        return Ok(None);
    }
    let v = across.scaled(&(&Real::one() / radius));
    let toward = model
        .geom()
        .position(model.vertex(far).point)
        .minus(&input.center);
    let away = toward
        .minus(&q.scaled(&toward.dot(&q)))
        .scaled(&Real::from_i64(-1));
    if sign(&away.norm_squared(), dec)? == Sign::Zero {
        return Ok(None);
    }
    let t = unit(&away);
    let two = Real::from_i64(2).sqrt();
    let three = Real::from_i64(3).sqrt();
    let lo = &three - &two;
    let zero = Real::zero;
    let weights = bernstein(
        [Real::one(), zero(), Real::from_i64(6), zero(), Real::one()],
        &lo,
    );
    let x = bernstein(
        [
            zero(),
            &Real::from_i64(4) * &two,
            zero(),
            -&(&Real::from_i64(4) * &two),
            zero(),
        ],
        &lo,
    );
    let z = bernstein(
        [
            &Real::from_i64(2) * &two,
            zero(),
            zero(),
            zero(),
            -&(&Real::from_i64(2) * &two),
        ],
        &lo,
    );
    let node = torus.center.plus(&q.scaled(radius));
    let foot = input
        .center
        .plus(&q.scaled(radius))
        .plus(&t.scaled(&(radius * &three)));
    let mut controls: Vec<_> = (0..5)
        .map(|i| {
            let numerator = node
                .minus(&Point3::origin())
                .scaled(&weights[i])
                .plus(&v.scaled(&(radius * &x[i])))
                .plus(&t.scaled(&(radius * &z[i])));
            [numerator.x, numerator.y, numerator.z, weights[i].clone()]
        })
        .collect();
    // Substituting u=√3-√2 and u=1 gives these simpler exact endpoint
    // coordinates. Keep the shared vertices and curve ends in that form.
    for (index, point) in [(0, &foot), (4, &node)] {
        controls[index] = [
            &point.x * &weights[index],
            &point.y * &weights[index],
            &point.z * &weights[index],
            weights[index].clone(),
        ];
    }
    let mut knots = vec![Rational::zero(); 5];
    knots.extend(vec![Rational::one(); 5]);
    let curve = BSplineCurve::from_homogeneous(4, knots, controls, dec)?;
    let mut feet = [node.clone(), node];
    feet[planar] = foot.clone();
    let hull_end = input.center.plus(
        &foot
            .minus(&input.center)
            .scaled(&(&Real::one() / &Real::from_i64(2))),
    );
    Ok(Some(End {
        feet,
        curve: Curve::BSpline(curve),
        hull_end,
    }))
}
