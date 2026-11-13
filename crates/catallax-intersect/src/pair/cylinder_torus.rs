//! Parallel, distinct axes: intersect the two exact radial intervals first.

use super::*;
use catallax_shapes::Torus;

pub(super) fn parallel(
    cylinder: &Cylinder,
    torus: &Torus,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    if !vector_zero(&cylinder.axis.cross(&torus.axis), dec)? {
        return Err(Refusal::new(
            "intersect.traced-pair-required",
            "skew cylinder/torus intersection requires certified branch tracing",
        ));
    }
    let origin = cylinder.origin.plus(&torus.axis.scaled(
        &(torus.center.minus(&cylinder.origin).dot(&torus.axis) / torus.axis.norm_squared()),
    ));
    let offset = origin.minus(&torus.center);
    let distance_sq = offset.norm_squared();
    let distance = root_real(&distance_sq);
    let radius = cylinder.radius();
    let cylinder_lower = if compare(&distance, &radius, dec)? == Sign::Negative {
        &radius - &distance
    } else {
        &distance - &radius
    };
    let cylinder_upper = &distance + &radius;
    let lower_radius = &torus.major - &torus.minor;
    let upper_radius = &torus.major + &torus.minor;
    let lower = if compare(&cylinder_lower, &lower_radius, dec)? == Sign::Positive {
        cylinder_lower.clone()
    } else {
        lower_radius.clone()
    };
    let upper = if compare(&cylinder_upper, &upper_radius, dec)? == Sign::Negative {
        cylinder_upper.clone()
    } else {
        upper_radius.clone()
    };
    let overlap = compare(&upper, &lower, dec)?;
    if overlap == Sign::Negative {
        return Ok(Plan::disjoint());
    }

    let e1 = offset.scaled(&(&s(1) / &distance));
    let e3 = unit(&torus.axis);
    let e2 = e3.cross(&e1);
    let radial_shift = &distance_sq + &cylinder.radius_sq;
    let radial_slope = s(2) * distance;
    let lower_sq = lower.square();
    let lo = &(&lower_sq - &radial_shift) / &radial_slope;
    if overlap == Sign::Zero {
        // Both radial endpoints coincide: y=z=0, one isolated contact.
        return Ok(Plan::point(origin.plus(&e1.scaled(&lo))));
    }
    let span = &(upper.square() - lower_sq) / &radial_slope;
    let y_ends = [
        compare(&lower, &cylinder_lower, dec)? == Sign::Zero,
        compare(&upper, &cylinder_upper, dec)? == Sign::Zero,
    ];
    let z_ends = [
        compare(&lower, &lower_radius, dec)? == Sign::Zero,
        compare(&upper, &upper_radius, dec)? == Sign::Zero,
    ];
    // Use the existing nodal-loop convention: an axial divider has priority
    // when both intervals coincide. Each loop then retains both shared nodes
    // and its distinct arriving/leaving tangents instead of crossing branches.
    let mode = if y_ends == [true, true] {
        CylinderTorusMode::Axial
    } else if z_ends == [true, true] {
        CylinderTorusMode::Transverse
    } else {
        CylinderTorusMode::Single
    };
    let divider = match mode {
        CylinderTorusMode::Single => None,
        CylinderTorusMode::Axial => Some(Divider::conic(origin.clone(), e3.clone())),
        CylinderTorusMode::Transverse => Some(Divider::conic(origin.clone(), e2.clone())),
    };
    let mut points = Vec::new();
    for end in 0..2 {
        if y_ends[end] && z_ends[end] {
            let x = if end == 0 { lo.clone() } else { &lo + &span };
            points.push(origin.plus(&e1.scaled(&x)));
        }
    }
    let x = &lo + &(&span / &s(2));
    let y = root_real(&(&cylinder.radius_sq - &x.square()));
    let radial = root_real(&(&radial_shift + &(&radial_slope * &x)));
    let z = root_real(&((&upper_radius - &radial) * (&radial - &lower_radius)));
    let base = origin.plus(&e1.scaled(&x));
    let positive = base.plus(&e2.scaled(&y)).plus(&e3.scaled(&z));
    let mut curves = Vec::with_capacity(if mode == CylinderTorusMode::Single {
        1
    } else {
        2
    });
    curves.push(Spec {
        form: CurveForm::Branch,
        witness: Some(positive),
        contact: Contact::Crossing,
        singular: (0..points.len()).collect(),
    });
    if mode != CylinderTorusMode::Single {
        let negative = match mode {
            CylinderTorusMode::Axial => base.plus(&e2.scaled(&y)).plus(&e3.scaled(&(-z))),
            CylinderTorusMode::Transverse => base.plus(&e2.scaled(&(-y))).plus(&e3.scaled(&z)),
            CylinderTorusMode::Single => unreachable!(),
        };
        curves.push(Spec {
            form: CurveForm::Branch,
            witness: Some(negative),
            contact: Contact::Crossing,
            singular: (0..points.len()).collect(),
        });
    }
    Ok(Plan {
        terminal: None,
        points,
        curves,
        divider,
        frame: Some(Frame::CylinderTorus {
            origin,
            e1,
            e2,
            e3,
            radius,
            radial_shift,
            radial_slope,
            lo,
            span,
            lower_radius,
            upper_radius,
            y_ends,
            z_ends,
            mode,
            sense: if sign(&sigma, dec)? == Sign::Positive {
                1
            } else {
                -1
            },
        }),
    })
}
