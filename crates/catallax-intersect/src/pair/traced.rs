//! Native analytic windows feeding the existing complete patch-pair engine.

use super::trace_candidates::Candidates;
use super::*;
use catallax_shapes::{
    BSplineSurface, Isometry, RPoint3, RVector3, RationalAtlas, SectionSurface, SplineSupport,
    Torus,
};

fn refused(detail: impl Into<String>) -> Refusal {
    Refusal::new("intersect.trace-native-cover", detail)
}

fn rational(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.clone());
    }
    dec.rational_value(value)?
        .ok_or_else(|| refused("the native window needs exact rational coefficient admission"))
}

fn point(value: &Point3, dec: &mut Decider) -> Result<RPoint3, Refusal> {
    Ok(RPoint3 {
        x: rational(&value.x, dec)?,
        y: rational(&value.y, dec)?,
        z: rational(&value.z, dec)?,
    })
}

fn axis(value: &Vector3, dec: &mut Decider) -> Result<RVector3, Refusal> {
    let length = value.norm_squared().sqrt();
    if sign(&length, dec)? != Sign::Positive {
        return Err(refused("the native window needs a nonzero axis"));
    }
    let direction = value.scaled(&(&Real::one() / &length));
    let direction = RVector3 {
        x: rational(&direction.x, dec)?,
        y: rational(&direction.y, dec)?,
        z: rational(&direction.z, dec)?,
    };
    if direction.norm_squared() != r(1) {
        return Err(refused(
            "the admitted rational frame axis is not exactly unit length",
        ));
    }
    Ok(direction)
}

fn frame(direction: &RVector3, center: &RPoint3) -> Result<Isometry, Refusal> {
    let origin = RPoint3::origin();
    let normal = RVector3 {
        x: -&direction.x,
        y: -&direction.y,
        z: &r(1) - &direction.z,
    };
    let rotation = if normal.is_zero() {
        Isometry::identity()
    } else {
        let flip_x = Isometry::reflection(
            &origin,
            &RVector3 {
                x: r(1),
                y: r(0),
                z: r(0),
            },
        )?;
        flip_x.then(&Isometry::reflection(&origin, &normal)?)
    };
    Ok(rotation.then(&Isometry::translation(center.minus(&origin))))
}

// Each quarter has half-angle h = 7u/(10-3u), whose inverse
// u = 10h/(7+3h) maps [0,1] onto [0,1]. Its positive denominator and
// X²+Y²=W² prove the whole quadrant, not just the control poles.
fn circle() -> [(Real, Real, Real); 9] {
    let poles = [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
        (1, 0),
    ];
    let middle = Real::from(Rational::from_ratio(7, 10).expect("positive denominator"));
    let next = Real::from(Rational::from_ratio(49, 50).expect("positive denominator"));
    let mut scale = Real::one();
    std::array::from_fn(|index| {
        let weight = if index % 2 == 0 {
            scale.clone()
        } else {
            &scale * &middle
        };
        let (x, y) = poles[index];
        if index % 2 == 1 {
            scale = &scale * &next;
        }
        (s(x), s(y), weight)
    })
}

fn knots() -> Vec<Rational> {
    [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
        .map(Rational::from_i64)
        .to_vec()
}

fn cylinder_source(
    radius: &Rational,
    bound: &Rational,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let radius = Real::from(radius.clone());
    let ends = [Real::from(-bound), Real::from(bound.clone())];
    let mut controls = Vec::with_capacity(18);
    let mut weights = Vec::with_capacity(18);
    for (x, y, w) in circle() {
        for z in &ends {
            controls.push(Point3::new(&radius * &x, &radius * &y, z.clone()));
            weights.push(w.clone());
        }
    }
    BSplineSurface::new(
        [2, 1],
        [knots(), [0, 0, 1, 1].map(Rational::from_i64).to_vec()],
        [9, 2],
        controls,
        weights,
        dec,
    )
}

fn torus_source(
    major: &Rational,
    minor: &Rational,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let (major, minor) = (Real::from(major.clone()), Real::from(minor.clone()));
    let circle = circle();
    let mut controls = Vec::with_capacity(81);
    let mut weights = Vec::with_capacity(81);
    for (x, y, w) in &circle {
        for (cv, sv, wv) in &circle {
            let radius = &major + &(&minor * cv);
            controls.push(Point3::new(&radius * x, &radius * y, &minor * sv));
            weights.push(w * wv);
        }
    }
    BSplineSurface::new([2, 2], [knots(), knots()], [9, 9], controls, weights, dec)
}

fn admitted(
    source: BSplineSurface,
    periodic: [bool; 2],
    carrier: &SectionSurface,
    map: &Isometry,
    dec: &mut Decider,
) -> Result<RationalAtlas, Refusal> {
    let mut atlases = RationalAtlas::partition_with_polynomial_bounds(source, periodic, dec)?;
    if atlases.len() != 1 {
        return Err(refused(
            "the fixed native tensor definition needs one complete connected atlas",
        ));
    }
    let atlas = atlases.pop().expect("one complete atlas was checked");
    atlas.certify_analytic_carrier(carrier, dec)?;
    match map.surface(&Surface::BSpline(SplineSupport::RationalAtlas(atlas))) {
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => Ok(atlas),
        _ => Err(refused(
            "exact rigid transport did not preserve the admitted atlas",
        )),
    }
}

fn same_oriented_carrier(a: &Surface, b: &Surface, dec: &mut Decider) -> Result<(), Refusal> {
    let a = a
        .implicit()
        .ok_or_else(|| refused("the mapped carrier has no native equation"))?;
    let b = b
        .implicit()
        .ok_or_else(|| refused("the requested carrier has no native equation"))?;
    let right = b.terms();
    for (monomial, coefficient) in a.terms() {
        if rational(&coefficient, dec)?.is_zero() {
            continue;
        }
        let Some((_, other)) = right.iter().find(|(term, _)| *term == monomial) else {
            return Err(refused(
                "the mapped native carrier has a different exact polynomial support",
            ));
        };
        let factor = &rational(other, dec)? / &rational(&coefficient, dec)?;
        if factor.sign() != Sign::Positive {
            return Err(refused(
                "the native implicit equations have no positive exact scale",
            ));
        }
        let difference = &a.scale(&Real::from(factor)) - &b;
        for (_, coefficient) in difference.terms() {
            if !rational(&coefficient, dec)?.is_zero() {
                return Err(refused(
                    "exact oriented equivalence to the requested native carrier did not prove",
                ));
            }
        }
        return Ok(());
    }
    Err(refused(
        "the mapped native carrier has a zero implicit polynomial",
    ))
}

fn zero_point() -> Point3 {
    Point3::new(s(0), s(0), s(0))
}

fn z_axis() -> Vector3 {
    Vector3::new(s(0), s(0), s(1))
}

/// Ring tori are compact. Relative to the cylinder-axis projection of the
/// torus centre, every meeting has |z| <= R+r. The strict extra unit is a
/// construction margin, not an intersection tolerance. Four exact circle
/// quarters on each periodic axis cover the whole native support in this window.
pub(super) fn cylinder_torus(
    cylinder: &Cylinder,
    torus: &Torus,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let mut stage = "bounded native window";
    let result = (|| -> Result<Plan, Refusal> {
        let radius = rational(&cylinder.radius(), dec)?;
        let major = rational(&torus.major, dec)?;
        let minor = rational(&torus.minor, dec)?;
        if radius.sign() != Sign::Positive || minor.sign() != Sign::Positive || major <= minor {
            return Err(refused(
                "the complete native window needs a positive cylinder and a strict ring torus",
            ));
        }
        let cylinder_axis = axis(&cylinder.axis, dec)?;
        let torus_axis = axis(&torus.axis, dec)?;
        let origin = point(&cylinder.origin, dec)?;
        let center = point(&torus.center, dec)?;
        let offset = center.minus(&origin).dot(&cylinder_axis);
        let projected = origin
            .to_real()
            .plus(&cylinder_axis.to_real().scaled(&Real::from(offset)));
        let projected = projected
            .as_rational()
            .ok_or_else(|| refused("the exact cylinder-axis projection was not admitted"))?;
        if !center.minus(&projected).dot(&cylinder_axis).is_zero() {
            return Err(refused(
                "the cylinder window centre has no exact orthogonal projection certificate",
            ));
        }
        let bound = &(&major + &minor) + &r(1);
        let cylinder_map = frame(&cylinder_axis, &projected)?;
        let torus_map = frame(&torus_axis, &center)?;
        let canonical_cylinder = SectionSurface::Cylinder(Cylinder {
            origin: zero_point(),
            axis: z_axis(),
            radius_sq: Real::from(&radius * &radius),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        });
        let canonical_torus = SectionSurface::Torus(Torus {
            center: zero_point(),
            axis: z_axis(),
            major: Real::from(major.clone()),
            minor: Real::from(minor.clone()),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        });
        stage = "exact oriented native-carrier equivalence";
        same_oriented_carrier(
            &cylinder_map.surface(&canonical_cylinder.surface()),
            &Surface::Cylinder(cylinder.clone()),
            dec,
        )?;
        same_oriented_carrier(
            &torus_map.surface(&canonical_torus.surface()),
            &Surface::Torus(torus.clone()),
            dec,
        )?;
        let cylinder_source = cylinder_source(&radius, &bound, dec)?;
        let torus_source = torus_source(&major, &minor, dec)?;
        let swap = sign(&sigma, dec)? == Sign::Negative;
        stage = "unproved branch proposals";
        let candidates = Candidates::generate(
            if swap {
                [
                    (&torus_source, &torus_map),
                    (&cylinder_source, &cylinder_map),
                ]
            } else {
                [
                    (&cylinder_source, &cylinder_map),
                    (&torus_source, &torus_map),
                ]
            },
            dec,
        )?;
        stage = "complete canonical cylinder admission";
        let cylinder = admitted(
            cylinder_source,
            [true, false],
            &canonical_cylinder,
            &cylinder_map,
            dec,
        )?;
        stage = "complete canonical torus admission";
        let torus = admitted(
            torus_source,
            [true, true],
            &canonical_torus,
            &torus_map,
            dec,
        )?;
        stage = "complete chart-pair preparation and certificates";
        let meeting = if swap {
            torus.patch_meeting_with_cached_spans(&cylinder, dec)?
        } else {
            cylinder.patch_meeting_with_cached_spans(&torus, dec)?
        };
        stage = "retained exact branch geometry and closed topology";
        let mut components = Vec::with_capacity(meeting.branches().len());
        for (index, branch) in meeting.branches().iter().enumerate() {
            if !branch.is_closed() {
                return Err(refused(
                    "a native component reaches the strict complete-window boundary",
                ));
            }
            components.push(meeting.curve_branch(index, dec)?);
        }
        stage = "exact proposal selection of the complete component cover";
        candidates.prove_selection(&components, dec)?;
        if components.is_empty() {
            return Ok(Plan::disjoint());
        }
        let mut plan = Plan::empty();
        plan.curves = components
            .into_iter()
            .map(|curve| Spec {
                form: CurveForm::PatchPair(curve),
                witness: None,
                contact: Contact::Crossing,
                singular: Vec::new(),
            })
            .collect();
        Ok(plan)
    })();
    result.map_err(|error| {
        Refusal::new(
            error.code,
            format!("traced cylinder/torus {stage}: {}", error.detail),
        )
    })
}
