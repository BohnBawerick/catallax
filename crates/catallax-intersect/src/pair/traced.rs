//! Native analytic windows feeding the existing complete patch-pair engine.

use super::trace_candidates::Candidates;
use super::*;
use catallax_shapes::{
    BSplineSurface, Cone, Isometry, RPoint3, RVector3, RationalAtlas, SectionSurface,
    SplineSupport, Torus,
};

pub(super) mod fitted;
pub(super) mod sphere;

pub(super) fn plan(a: &Surface, b: &Surface, dec: &mut Decider) -> Result<Plan, Refusal> {
    for (analytic, spline, swap) in [(a, b, false), (b, a, true)] {
        if let Surface::BSpline(SplineSupport::RollingEnvelope(envelope)) = spline
            && let Some(cutter) = SectionSurface::of(analytic)
        {
            return fitted::section(&cutter, envelope, swap, dec);
        }
    }
    match (a, b) {
        (Surface::Cylinder(c), Surface::Torus(t)) => cylinder_torus(c, t, s(1), dec),
        (Surface::Torus(t), Surface::Cylinder(c)) => cylinder_torus(c, t, s(-1), dec),
        (Surface::Cone(c), Surface::Cylinder(y)) => cone_cylinder(c, y, s(1), dec),
        (Surface::Cylinder(y), Surface::Cone(c)) => cone_cylinder(c, y, s(-1), dec),
        (Surface::Cone(c), Surface::Sphere(sphere)) => sphere::cone_sphere(c, sphere, s(1), dec),
        (Surface::Sphere(sphere), Surface::Cone(c)) => sphere::cone_sphere(c, sphere, s(-1), dec),
        _ => Err(Refusal::new(
            "intersect.trace-native-window",
            "this pair has no certified finite source window for exact branch tracing",
        )),
    }
}

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
        complete_pair(
            [
                NativeWindow {
                    source: cylinder_source,
                    periodic: [true, false],
                    carrier: canonical_cylinder,
                    map: cylinder_map,
                    admission_stage: "complete canonical cylinder admission",
                },
                NativeWindow {
                    source: torus_source,
                    periodic: [true, true],
                    carrier: canonical_torus,
                    map: torus_map,
                    admission_stage: "complete canonical torus admission",
                },
            ],
            &sigma,
            &mut stage,
            dec,
        )
    })();
    result.map_err(|error| {
        Refusal::new(
            error.code,
            format!("traced cylinder/torus {stage}: {}", error.detail),
        )
    })
}

struct NativeWindow {
    source: BSplineSurface,
    periodic: [bool; 2],
    carrier: SectionSurface,
    map: Isometry,
    admission_stage: &'static str,
}

fn complete_pair(
    [first, second]: [NativeWindow; 2],
    sigma: &Real,
    stage: &mut &'static str,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let swap = match sign(sigma, dec)? {
        Sign::Negative => true,
        Sign::Positive => false,
        Sign::Zero => return Err(refused("the native pair has no strict orientation")),
    };
    *stage = "unproved branch proposals";
    let candidates = Candidates::generate(
        if swap {
            [(&second.source, &second.map), (&first.source, &first.map)]
        } else {
            [(&first.source, &first.map), (&second.source, &second.map)]
        },
        dec,
    )?;
    *stage = first.admission_stage;
    let first = admitted(
        first.source,
        first.periodic,
        &first.carrier,
        &first.map,
        dec,
    )?;
    *stage = second.admission_stage;
    let second = admitted(
        second.source,
        second.periodic,
        &second.carrier,
        &second.map,
        dec,
    )?;
    *stage = "complete chart-pair preparation and certificates";
    let meeting = if swap {
        second.patch_meeting_with_cached_spans(&first, dec)?
    } else {
        first.patch_meeting_with_cached_spans(&second, dec)?
    };
    *stage = "retained exact branch geometry and closed topology";
    let mut components = Vec::with_capacity(meeting.branches().len());
    for (index, branch) in meeting.branches().iter().enumerate() {
        if !branch.is_closed() {
            return Err(refused(
                "a native component reaches the strict complete-window boundary",
            ));
        }
        components.push(meeting.curve_branch(index, dec)?);
    }
    *stage = "exact proposal selection of the complete component cover";
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
}

fn cone_source(
    slope: &Rational,
    interval: [&Rational; 2],
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let heights = interval.map(|height| Real::from(height.clone()));
    let radii = interval.map(|height| Real::from(slope * height));
    let mut controls = Vec::with_capacity(18);
    let mut weights = Vec::with_capacity(18);
    for (x, y, weight) in circle() {
        for end in 0..2 {
            controls.push(Point3::new(
                &radii[end] * &x,
                &radii[end] * &y,
                heights[end].clone(),
            ));
            weights.push(weight.clone());
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

/// For perpendicular unit axes d (cone) and n (cylinder), a point on both
/// has |(A-O).d + v| <= cylinder_radius. Thus every positive cone height
/// is bounded. A strict positive lower height excludes the cone's apex.
pub(super) fn cone_cylinder(
    cone: &Cone,
    cylinder: &Cylinder,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let mut stage = "bounded native window";
    let result = (|| -> Result<Plan, Refusal> {
        let slope = rational(&cone.slope, dec)?;
        let radius = rational(&cylinder.radius(), dec)?;
        if slope.sign() != Sign::Positive || radius.sign() != Sign::Positive {
            return Err(refused(
                "the native window needs positive cone slope and cylinder radius",
            ));
        }
        let cone_axis = axis(&cone.axis, dec)?;
        let cylinder_axis = axis(&cylinder.axis, dec)?;
        if !cone_axis.dot(&cylinder_axis).is_zero() {
            return Err(Refusal::new(
                "intersect.trace-native-window",
                "nonperpendicular cone/cylinder axes need a certified bounded native window",
            ));
        }
        let apex = point(&cone.apex, dec)?;
        let origin = point(&cylinder.origin, dec)?;
        let offset = apex.minus(&origin).dot(&cylinder_axis);
        let projected = origin
            .to_real()
            .plus(&cylinder_axis.to_real().scaled(&Real::from(offset)));
        let projected = projected
            .as_rational()
            .ok_or_else(|| refused("the exact cylinder-axis projection was not admitted"))?;
        let radial = apex.minus(&projected);
        if !radial.dot(&cylinder_axis).is_zero() {
            return Err(refused(
                "the cylinder window centre has no exact orthogonal projection certificate",
            ));
        }
        let height_offset = radial.dot(&cone_axis);
        let upper = &(-&height_offset) + &radius;
        if upper.sign() == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        if upper.is_zero() {
            return Err(Refusal::new(
                "intersect.trace-native-apex",
                "the only possible retained cone height reaches its singular apex",
            ));
        }
        let lower = &(-&height_offset) - &radius;
        let lower = if lower.sign() == Sign::Positive {
            &lower / &r(2)
        } else {
            let gap = (&radial.norm_squared() - &(&radius * &radius)).abs();
            if gap.is_zero() {
                return Err(Refusal::new(
                    "intersect.trace-native-apex",
                    "the cone apex lies on the cylinder and needs a singular native chart proof",
                ));
            }
            // |radial(A)-R| = |radial(A)^2-R^2|/(radial(A)+R).
            // The L1 bound and sqrt(1+s²) <= 1+s give a strict positive
            // cone height without admitting an approximate square root.
            let radial_upper = &(&radial.x.abs() + &radial.y.abs()) + &radial.z.abs();
            &gap / &(&r(2) * &(&(&(&radial_upper + &radius) + &r(1)) * &(&r(1) + &slope)))
        };
        if lower >= upper {
            return Ok(Plan::disjoint());
        }
        let cone_upper = &upper + &r(1);
        let cylinder_bound = &(&slope * &upper) + &r(1);
        let cone_map = frame(&cone_axis, &apex)?;
        let cylinder_map = frame(&cylinder_axis, &projected)?;
        let canonical_cone = SectionSurface::Cone(Cone {
            apex: zero_point(),
            axis: z_axis(),
            slope: Real::from(slope.clone()),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        });
        let canonical_cylinder = SectionSurface::Cylinder(Cylinder {
            origin: zero_point(),
            axis: z_axis(),
            radius_sq: Real::from(&radius * &radius),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        });
        stage = "exact oriented native-carrier equivalence";
        same_oriented_carrier(
            &cone_map.surface(&canonical_cone.surface()),
            &Surface::Cone(cone.clone()),
            dec,
        )?;
        same_oriented_carrier(
            &cylinder_map.surface(&canonical_cylinder.surface()),
            &Surface::Cylinder(cylinder.clone()),
            dec,
        )?;
        let cone_source = cone_source(&slope, [&lower, &cone_upper], dec)?;
        let cylinder_source = cylinder_source(&radius, &cylinder_bound, dec)?;
        complete_pair(
            [
                NativeWindow {
                    source: cone_source,
                    periodic: [true, false],
                    carrier: canonical_cone,
                    map: cone_map,
                    admission_stage: "complete canonical cone admission",
                },
                NativeWindow {
                    source: cylinder_source,
                    periodic: [true, false],
                    carrier: canonical_cylinder,
                    map: cylinder_map,
                    admission_stage: "complete canonical cylinder admission",
                },
            ],
            &sigma,
            &mut stage,
            dec,
        )
    })();
    result.map_err(|error| {
        Refusal::new(
            error.code,
            format!("traced cone/cylinder {stage}: {}", error.detail),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_cone_apex_on_the_cylinder_is_not_dropped_from_the_native_cover() {
        let cone = Cone {
            apex: Point3::new(s(12), s(0), s(0)),
            axis: Vector3::new(s(-1), s(0), s(0)),
            slope: Real::from(Rational::from_ratio(1, 12).unwrap()),
            ref_dir: Vector3::new(s(0), s(1), s(0)),
        };
        let cylinder = Cylinder {
            origin: zero_point(),
            axis: z_axis(),
            radius_sq: s(144),
            ref_dir: Vector3::new(s(1), s(0), s(0)),
        };
        let refusal =
            cone_cylinder(&cone, &cylinder, Real::one(), &mut Decider::default()).unwrap_err();
        assert_eq!(refusal.code, "intersect.trace-native-apex");
    }
}
