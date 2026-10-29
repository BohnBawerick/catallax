//! Exact analytic conics with a native rational arithmetic path.
//! Regular witnesses retain their supports; no rounded coordinate becomes a definition.

use super::{Divider, Plan, Spec, compare, is_zero, s};
use crate::{
    Contact, CurveForm,
    math::{root_real as root, *},
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{Cone, Cylinder, Plane, Point3, Sphere, Torus, Vector3};

fn traced(detail: &str) -> Refusal {
    Refusal::new("intersect.traced-pair-required", detail)
}
fn algebraic_direction(detail: &str) -> Refusal {
    Refusal::new("intersect.algebraic-direction-required", detail)
}
fn coaxial(
    origin: &Point3,
    axis: &Vector3,
    other: &Point3,
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    Ok(is_zero(&axis.cross(direction).norm_squared(), dec)?
        && is_zero(&other.minus(origin).cross(axis).norm_squared(), dec)?)
}
fn sq(value: &Real) -> Real {
    value.square()
}
fn positive_lower(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    if let Some(q) = value.rational_ref() {
        return Ok(q.clone());
    }
    let mut precision = 64;
    loop {
        let lower = dec.enclose_budgeted_at(value, precision)?.lo();
        if lower.sign() == Sign::Positive {
            return Ok(lower);
        }
        precision = precision.checked_mul(2).ok_or_else(|| {
            Refusal::new(
                "numbers.refinement-budget",
                "divider separation exceeds the precision index",
            )
        })?;
    }
}
fn circle(
    center: Point3,
    normal: Vector3,
    radius_sq: Real,
    ref_dir: Vector3,
    contact: Contact,
) -> Spec {
    let witness = if [&center.x, &center.y, &center.z]
        .iter()
        .all(|x| x.rational_ref().is_some())
        && radius_sq.rational_ref().is_some()
    {
        None
    } else {
        Some(
            center.plus(
                &ref_dir
                    .clone()
                    .scaled(&root_real(&(&radius_sq / &real(&ref_dir.norm_squared())))),
            ),
        )
    };
    Spec {
        form: CurveForm::Circle {
            center,
            normal,
            radius_sq,
            ref_dir,
        },
        witness,
        contact,
        singular: vec![],
    }
}
fn one_circle(
    center: Point3,
    normal: Vector3,
    radius_sq: Real,
    ref_dir: Vector3,
    contact: Contact,
) -> Plan {
    Plan {
        curves: vec![circle(center, normal, radius_sq, ref_dir, contact)],
        ..Plan::empty()
    }
}
fn form(form: CurveForm, witness: Point3) -> Plan {
    Plan {
        curves: vec![Spec {
            form,
            witness: Some(witness),
            contact: Contact::Crossing,
            singular: vec![],
        }],
        ..Plan::empty()
    }
}

pub(super) fn plane_cone(
    p: &Plane,
    c: &Cone,
    sigma: Real,
    rational: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let axis2 = c.axis.norm_squared();
    let normal2 = p.normal.norm_squared();
    let dot = c.axis.dot(&p.normal);
    let offset = p.origin.minus(&c.apex).dot(&p.normal);
    let slope2 = sq(&c.slope);
    if is_zero(&c.axis.cross(&p.normal).norm_squared(), dec)? {
        let t = &offset / &dot;
        return Ok(match sign(&t, dec)? {
            Sign::Negative => Plan::disjoint(),
            Sign::Zero => Plan::point(c.apex.clone()),
            Sign::Positive => one_circle(
                c.apex.plus(&c.axis.scaled(&t)).clone(),
                p.normal.scaled(&sigma),
                real(&(slope2 * axis2 * sq(&t))),
                c.ref_dir.clone(),
                Contact::Crossing,
            ),
        });
    }
    let projected = c.axis.minus(&p.normal.scaled(&(&dot / &normal2)));
    let projected2 = projected.norm_squared();
    let transverse = p.normal.cross(&c.axis);
    let k = &s(1) + &slope2;
    let quadratic = &axis2 - &(&k * &projected2);
    if is_zero(&offset, dec)? {
        if sign(&quadratic, dec)? == Sign::Positive {
            return Ok(Plan::point(c.apex.clone()));
        }
        let factor_sq = -&quadratic / (&axis2 * &normal2);
        let factor = if rational {
            factor_sq
                .rational_ref()
                .and_then(rational_sqrt)
                .map(Real::from_rational)
                .ok_or_else(|| {
                    algebraic_direction("cone generators require algebraic direction fields")
                })?
        } else {
            root(&factor_sq)
        };
        let mut plan = Plan::empty();
        plan.points.push(c.apex.clone());
        let tangent = is_zero(&quadratic, dec)?;
        let count = if tangent { 1 } else { 2 };
        for side in [s(1), s(-1)].into_iter().take(count) {
            let direction = projected.plus(&transverse.scaled(&(&factor * &side)));
            let gradient = direction
                .scaled(&axis2)
                .minus(&c.axis.scaled(&(&k * &direction.dot(&c.axis))));
            let orientation = &p.normal.cross(&gradient).dot(&direction) * &sigma;
            let orientation_sign = sign(&orientation, dec)?;
            let natural = if orientation_sign == Sign::Negative
                || (orientation_sign == Sign::Zero && sign(&sigma, dec)? == Sign::Negative)
            {
                direction.scaled(&s(-1))
            } else {
                direction.clone()
            };
            plan.curves.push(Spec {
                form: CurveForm::Line {
                    origin: c.apex.clone(),
                    dir: natural,
                },
                witness: Some(c.apex.plus(&direction).clone()),
                contact: if tangent {
                    Contact::Tangent
                } else {
                    Contact::Crossing
                },
                singular: vec![0],
            });
        }
        if !tangent {
            plan.divider = Some(Divider::conic(c.apex.clone(), transverse));
        }
        return Ok(plan);
    }
    let beta = &(&dot * &offset) / &normal2;
    let foot = c.apex.plus(&p.normal.scaled(&(&offset / &normal2)));
    if is_zero(&quadratic, dec)? {
        if sign(&beta, dec)? == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        let plane_term = &(&axis2 * &sq(&offset)) / &normal2;
        let constant = &plane_term - &(&k * &sq(&beta));
        let linear = &k * &beta;
        let shift = &constant / &(s(2) * (&linear * &projected2));
        let vertex = foot.plus(&projected.scaled(&shift)).clone();
        let focal = real(&(&linear / &(s(2) * axis2))) * root(&projected2);
        return Ok(form(
            CurveForm::Parabola {
                vertex: vertex.clone(),
                normal: p.normal.scaled(&(-sigma)),
                axis_dir: projected,
                focal,
            },
            vertex,
        ));
    }
    let center = foot.plus(&projected.scaled(&(&(&k * &beta) / &quadratic)));
    let scale = (&slope2 * &sq(&axis2)) * sq(&offset) / (&normal2 * &quadratic);
    let major2 = &scale / &quadratic;
    if sign(&quadratic, dec)? == Sign::Positive {
        if sign(&beta, dec)? == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        return Ok(Plan::form(
            CurveForm::Ellipse {
                center: center.clone(),
                normal: p.normal.scaled(&sigma),
                major_dir: projected,
                major_sq: real(&major2),
                minor_sq: real(&(&scale / &axis2)),
            },
            Contact::Crossing,
        ));
    }
    let witness = center
        .clone()
        .plus(&projected.clone().scaled(&root(&(&major2 / &projected2))));
    Ok(form(
        CurveForm::Hyperbola {
            center: center.clone(),
            normal: p.normal.scaled(&(-sigma)),
            axis_dir: projected,
            transverse_sq: real(&major2),
            conjugate_sq: real(&(-scale / axis2)),
            branch: StrictSign::Positive,
        },
        witness,
    ))
}

pub(super) fn cone_cylinder(
    c: &Cone,
    y: &Cylinder,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    if !coaxial(&c.apex, &c.axis, &y.origin, &y.axis, dec)? {
        return Err(traced(
            "non-coaxial cone/cylinder intersection requires certified branch tracing",
        ));
    }
    let t = root(&(&y.radius_sq / &(sq(&c.slope) * c.axis.norm_squared())));
    Ok(one_circle(
        c.apex.clone().plus(&c.axis.clone().scaled(&t)),
        c.axis.scaled(&(-sigma)),
        real(&y.radius_sq),
        c.ref_dir.clone(),
        Contact::Crossing,
    ))
}

// Sign of a + side sqrt(discriminant), without approximating the root.
fn rational_root_sign(a: &Rational, discriminant: &Rational, side: i64) -> Sign {
    if discriminant.is_zero() {
        return a.sign();
    }
    if side > 0 {
        if a.sign() != Sign::Negative {
            Sign::Positive
        } else {
            (discriminant - &(a * a)).sign()
        }
    } else if a.sign() != Sign::Positive {
        Sign::Negative
    } else {
        (&(a * a) - discriminant).sign()
    }
}
fn root_sign(a: &Real, discriminant: &Real, side: i64, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let (Some(a), Some(discriminant)) = (a.rational_ref(), discriminant.rational_ref()) {
        return Ok(rational_root_sign(a, discriminant, side));
    }
    if is_zero(discriminant, dec)? {
        return Ok(sign(a, dec)?);
    }
    let first = sign(a, dec)?;
    Ok(if side > 0 {
        if first != Sign::Negative {
            Sign::Positive
        } else {
            sign(&(discriminant - &sq(a)), dec)?
        }
    } else if first != Sign::Positive {
        Sign::Negative
    } else {
        sign(&(&sq(a) - discriminant), dec)?
    })
}

pub(super) fn cone_sphere(
    c: &Cone,
    sphere: &Sphere,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    if !is_zero(
        &sphere.center.minus(&c.apex).cross(&c.axis).norm_squared(),
        dec,
    )? {
        return Err(traced(
            "non-coaxial cone/sphere intersection requires certified branch tracing",
        ));
    }
    let axis2 = c.axis.norm_squared();
    let h = &sphere.center.minus(&c.apex).dot(&c.axis) / &axis2;
    let k = s(1) + sq(&c.slope);
    let discriminant = &(&k * &sq(&sphere.radius)) / &axis2 - sq(&c.slope) * sq(&h);
    if sign(&discriminant, dec)? == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset = root(&discriminant);
    for side in [1, -1]
        .into_iter()
        .take(if is_zero(&discriminant, dec)? { 1 } else { 2 })
    {
        match root_sign(&h, &discriminant, side, dec)? {
            Sign::Negative => continue,
            Sign::Zero => {
                plan.points.push(c.apex.clone());
                continue;
            }
            Sign::Positive => {}
        }
        let t = (real(&h) + &Real::from_i64(side) * &offset) / real(&k);
        let center = c.apex.clone().plus(&c.axis.clone().scaled(&t));
        plan.curves.push(circle(
            center,
            c.axis.scaled(&(&s(-side) * &sigma)),
            real(&(&sq(&c.slope) * &axis2)) * t.square(),
            c.ref_dir.clone(),
            if is_zero(&discriminant, dec)? {
                Contact::Tangent
            } else {
                Contact::Crossing
            },
        ));
    }
    if plan.curves.len() == 2 {
        plan.divider = Some(Divider::conic(
            c.apex.plus(&c.axis.scaled(&(h / k))),
            c.axis.clone(),
        ));
    }
    if plan.curves.is_empty() && plan.points.is_empty() {
        Ok(Plan::disjoint())
    } else {
        Ok(plan)
    }
}

pub(super) fn cones(a: &Cone, b: &Cone, dec: &mut Decider) -> Result<Plan, Refusal> {
    if !coaxial(&a.apex, &a.axis, &b.apex, &b.axis, dec)? {
        return Err(traced(
            "non-coaxial cone/cone intersection requires certified branch tracing",
        ));
    }
    let sense = if sign(&a.axis.dot(&b.axis), dec)? == Sign::Positive {
        s(1)
    } else {
        s(-1)
    };
    let h = b.apex.minus(&a.apex).dot(&a.axis) / a.axis.norm_squared();
    let signed_slope = &sense * &b.slope;
    let denominator = &a.slope - &signed_slope;
    if is_zero(&denominator, dec)? {
        return Ok(if is_zero(&h, dec)? {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        });
    }
    let t = &(-(&signed_slope * &h)) / &denominator;
    if sign(&t, dec)? == Sign::Negative || sign(&(&sense * &(&t - &h)), dec)? == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    if is_zero(&t, dec)? {
        return Ok(Plan::point(a.apex.clone()));
    }
    let normal = a
        .axis
        .scaled(&if sign(&denominator, dec)? == Sign::Positive {
            s(-1)
        } else {
            s(1)
        });
    Ok(one_circle(
        a.apex.plus(&a.axis.scaled(&t)).clone(),
        normal,
        real(&(sq(&a.slope) * a.axis.norm_squared() * sq(&t))),
        a.ref_dir.clone(),
        Contact::Crossing,
    ))
}

pub(super) fn plane_torus(
    p: &Plane,
    t: &Torus,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let offset = p.origin.minus(&t.center).dot(&p.normal);
    let dot = p.normal.dot(&t.axis);
    if is_zero(&p.normal.cross(&t.axis).norm_squared(), dec)? {
        let height = &offset / &dot;
        let center = t.center.plus(&t.axis.scaled(&height));
        let discriminant = sq(&t.minor) - t.axis.norm_squared() * sq(&height);
        if sign(&discriminant, dec)? == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        if is_zero(&discriminant, dec)? {
            return Ok(one_circle(
                center.clone(),
                p.normal.scaled(&sigma),
                real(&sq(&t.major)),
                t.ref_dir.clone(),
                Contact::Tangent,
            ));
        }
        let mut plan = Plan::empty();
        let offset = root(&discriminant);
        for side in [s(1), s(-1)] {
            let radius = real(&t.major) + &real(&side) * &offset;
            plan.curves.push(circle(
                center.clone(),
                p.normal.scaled(&(&side * &sigma)),
                radius.square(),
                t.ref_dir.clone(),
                Contact::Crossing,
            ));
        }
        plan.divider = Some(Divider::radial(
            center,
            t.axis.clone(),
            sq(&t.major) + discriminant,
            positive_lower(&(&t.major - &t.minor), dec)?,
        ));
        return Ok(plan);
    }
    if is_zero(&dot, dec)? && is_zero(&offset, dec)? {
        let radial = p.normal.cross(&t.axis);
        let shift = radial
            .clone()
            .scaled(&(real(&t.major) / root(&radial.norm_squared())));
        let mut plan = Plan::empty();
        for side in [1, -1] {
            let center = t.center.clone().plus(&shift.scaled(&Real::from_i64(side)));
            plan.curves.push(circle(
                center,
                p.normal.scaled(&sigma),
                real(&sq(&t.minor)),
                t.axis.clone(),
                Contact::Crossing,
            ));
        }
        plan.divider = Some(Divider::conic(t.center.clone(), radial));
        return Ok(plan);
    }
    // The support function of a ring torus is R |n_perp| + r |n|.
    // Compare the sum of square roots over the rational field.
    let radial_normal = p
        .normal
        .minus(&t.axis.scaled(&(&dot / &t.axis.norm_squared())));
    let radial2 = radial_normal.norm_squared();
    let normal2 = p.normal.norm_squared();
    let excess = sq(&offset) - (&sq(&t.major) * &radial2) - (&sq(&t.minor) * &normal2);
    let cross2 = s(4) * sq(&t.major) * sq(&t.minor) * (&radial2 * &normal2);
    if sign(&excess, dec)? == Sign::Positive {
        match sign(&(sq(&excess) - cross2), dec)? {
            Sign::Positive => return Ok(Plan::disjoint()),
            Sign::Zero => {
                let orientation = if sign(&offset, dec)? == Sign::Positive {
                    Real::one()
                } else {
                    Real::from_i64(-1)
                };
                let displacement = radial_normal
                    .clone()
                    .scaled(&root(&(sq(&t.major) / radial2)))
                    .plus(&p.normal.clone().scaled(&root(&(sq(&t.minor) / normal2))));
                return Ok(Plan::point(
                    t.center.clone().plus(&displacement.scaled(&orientation)),
                ));
            }
            Sign::Negative => {}
        }
    }
    Err(traced(
        "oblique or offset axial torus section requires certified branch tracing",
    ))
}

pub(super) fn cylinder_torus(
    c: &Cylinder,
    t: &Torus,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    if !coaxial(&t.center, &t.axis, &c.origin, &c.axis, dec)? {
        return super::cylinder_torus::parallel(c, t, sigma, dec);
    }
    let discriminant = sq(&t.minor) - sq(&(&c.radius() - &t.major));
    if sign(&discriminant, dec)? == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset = root(&(&discriminant / &t.axis.norm_squared()));
    for side in [1, -1]
        .into_iter()
        .take(if is_zero(&discriminant, dec)? { 1 } else { 2 })
    {
        let center = t
            .center
            .clone()
            .plus(&t.axis.clone().scaled(&(&Real::from_i64(side) * &offset)));
        plan.curves.push(circle(
            center,
            t.axis.scaled(&(&s(-side) * &sigma)),
            real(&c.radius_sq),
            c.ref_dir.clone(),
            if is_zero(&discriminant, dec)? {
                Contact::Tangent
            } else {
                Contact::Crossing
            },
        ));
    }
    if plan.curves.len() == 2 {
        plan.divider = Some(Divider::conic(t.center.clone(), t.axis.clone()));
    }
    Ok(plan)
}

pub(super) fn sphere_torus(
    sphere: &Sphere,
    t: &Torus,
    sigma: Real,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    let displacement = sphere.center.minus(&t.center);
    let axis2 = t.axis.norm_squared();
    if !is_zero(&displacement.cross(&t.axis).norm_squared(), dec)? {
        let axial = displacement.dot(&t.axis);
        let radial2 = displacement.norm_squared() - &sq(&axial) / &axis2;
        let distance_term = displacement.norm_squared() + sq(&t.major);
        let root2 = s(4) * (&sq(&t.major) * &radial2);
        // Distance to the core circle is 1-Lipschitz. These strict tests
        // certify separated surfaces, including a sphere inside the tube.
        if (compare(&t.minor, &sphere.radius, dec)? == Sign::Positive
            && root_sign(
                &(&distance_term - &sq(&(&t.minor - &sphere.radius))),
                &root2,
                -1,
                dec,
            )? == Sign::Negative)
            || root_sign(
                &(&distance_term - &sq(&(&t.minor + &sphere.radius))),
                &root2,
                -1,
                dec,
            )? == Sign::Positive
            || (compare(&sphere.radius, &(&t.major + &t.minor), dec)? == Sign::Positive
                && compare(
                    &displacement.norm_squared(),
                    &sq(&(&(&sphere.radius - &t.major) - &t.minor)),
                    dec,
                )? == Sign::Negative)
        {
            return Ok(Plan::disjoint());
        }
        if is_zero(&axial, dec)?
            && compare(&radial2, &sq(&t.major), dec)? == Sign::Zero
            && compare(&sphere.radius, &t.minor, dec)? == Sign::Zero
        {
            return Ok(one_circle(
                sphere.center.clone(),
                t.axis.cross(&displacement),
                real(&sq(&sphere.radius)),
                t.axis.clone(),
                Contact::Tangent,
            ));
        }
        return Err(traced(
            "non-coaxial sphere/torus intersection requires certified branch tracing",
        ));
    }
    let h = &displacement.dot(&t.axis) / &axis2;
    let height2 = &axis2 * &sq(&h);
    let denominator = &sq(&t.major) + &height2;
    let offset = (sq(&sphere.radius) - sq(&t.minor) - height2 - sq(&t.major)) / s(2);
    let base_radius = &t.major + &(&(&t.major * &offset) / &denominator);
    let base_height = -(&(&h * &offset) / &denominator);
    let discriminant = &(sq(&t.minor) - &sq(&offset) / &denominator) / &denominator;
    if sign(&discriminant, dec)? == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset_root = root(&(&discriminant / &axis2));
    for side in [1, -1]
        .into_iter()
        .take(if is_zero(&discriminant, dec)? { 1 } else { 2 })
    {
        let root = &Real::from_i64(side) * &offset_root;
        let center = t.center.clone().plus(
            &t.axis
                .clone()
                .scaled(&(real(&base_height) + &real(&t.major) * &root)),
        );
        let radius = real(&base_radius) + real(&(&h * &axis2)) * root;
        plan.curves.push(circle(
            center,
            t.axis.scaled(&(&s(-side) * &sigma)),
            radius.square(),
            t.ref_dir.clone(),
            if is_zero(&discriminant, dec)? {
                Contact::Tangent
            } else {
                Contact::Crossing
            },
        ));
    }
    if plan.curves.len() == 2 {
        plan.divider = Some(Divider::conic(
            t.center.plus(&t.axis.scaled(&base_height)),
            t.axis.clone(),
        ));
    }
    Ok(plan)
}

pub(super) fn cone_torus(
    c: &Cone,
    t: &Torus,
    sigma: Real,
    rational: bool,
    dec: &mut Decider,
) -> Result<Plan, Refusal> {
    if !coaxial(&c.apex, &c.axis, &t.center, &t.axis, dec)? {
        return super::cone_torus::parallel(c, t, sigma, dec);
    }
    let axis_sq = c.axis.norm_squared();
    let length = if rational {
        axis_sq
            .rational_ref()
            .and_then(rational_sqrt)
            .map(Real::from_rational)
            .ok_or_else(|| {
                algebraic_direction(
                    "cone/torus meridian requires the shared algebraic coefficient foundation",
                )
            })?
    } else {
        root(&axis_sq)
    };
    let h = &t.center.minus(&c.apex).dot(&c.axis) / &length;
    let k = s(1) + sq(&c.slope);
    let middle = &(&c.slope * &t.major) + &h;
    let discriminant = sq(&middle) - &k * &(sq(&t.major) + sq(&h) - sq(&t.minor));
    if sign(&discriminant, dec)? == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset = root(&discriminant);
    for side in [1, -1]
        .into_iter()
        .take(if is_zero(&discriminant, dec)? { 1 } else { 2 })
    {
        match root_sign(&middle, &discriminant, side, dec)? {
            Sign::Negative => continue,
            Sign::Zero => {
                plan.points.push(c.apex.clone());
                continue;
            }
            Sign::Positive => {}
        }
        let height = (real(&middle) + &Real::from_i64(side) * &offset) / real(&k);
        let center = c
            .apex
            .clone()
            .plus(&c.axis.clone().scaled(&(&height / &real(&length))));
        plan.curves.push(circle(
            center,
            c.axis.scaled(&(&s(-side) * &sigma)),
            real(&sq(&c.slope)) * height.square(),
            c.ref_dir.clone(),
            if is_zero(&discriminant, dec)? {
                Contact::Tangent
            } else {
                Contact::Crossing
            },
        ));
    }
    if plan.curves.len() == 2 {
        plan.divider = Some(Divider::conic(
            c.apex.plus(&c.axis.scaled(&(middle / (k * length)))),
            c.axis.clone(),
        ));
    }
    if plan.curves.is_empty() && plan.points.is_empty() {
        Ok(Plan::disjoint())
    } else {
        Ok(plan)
    }
}

pub(super) fn tori(a: &Torus, b: &Torus, dec: &mut Decider) -> Result<Plan, Refusal> {
    if !coaxial(&a.center, &a.axis, &b.center, &b.axis, dec)? {
        return Err(traced(
            "non-coaxial torus/torus intersection requires certified branch tracing",
        ));
    }
    let axis2 = a.axis.norm_squared();
    let h = &b.center.minus(&a.center).dot(&a.axis) / &axis2;
    let radial = &b.major - &a.major;
    let distance2 = sq(&radial) + &axis2 * &sq(&h);
    if is_zero(&distance2, dec)? {
        return Ok(if compare(&a.minor, &b.minor, dec)? == Sign::Zero {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        });
    }
    if compare(&distance2, &sq(&(&a.minor + &b.minor)), dec)? == Sign::Positive
        || compare(&distance2, &sq(&(&a.minor - &b.minor)), dec)? == Sign::Negative
    {
        return Ok(Plan::disjoint());
    }
    let fraction = (&distance2 + &sq(&a.minor) - sq(&b.minor)) / (&s(2) * &distance2);
    let base_radius = &a.major + &(&radial * &fraction);
    let base_height = &h * &fraction;
    let discriminant = &(sq(&a.minor) - &sq(&fraction) * &distance2) / &distance2;
    let mut plan = Plan::empty();
    let offset_root = root(&(&discriminant / &axis2));
    let first = if sign(&radial, dec)? == Sign::Negative
        || (is_zero(&radial, dec)? && sign(&h, dec)? == Sign::Positive)
    {
        -1
    } else {
        1
    };
    for side in [first, -first]
        .into_iter()
        .take(if is_zero(&discriminant, dec)? { 1 } else { 2 })
    {
        let root = &Real::from_i64(side) * &offset_root;
        let center = a.center.clone().plus(
            &a.axis
                .clone()
                .scaled(&(real(&base_height) + &real(&radial) * &root)),
        );
        let radius = real(&base_radius) - real(&(&h * &axis2)) * root;
        plan.curves.push(circle(
            center,
            a.axis.scaled(&s(-side)),
            radius.square(),
            a.ref_dir.clone(),
            if is_zero(&discriminant, dec)? {
                Contact::Tangent
            } else {
                Contact::Crossing
            },
        ));
    }
    if plan.curves.len() == 2 {
        let base = a.center.plus(&a.axis.scaled(&base_height));
        plan.divider = Some(if is_zero(&radial, dec)? {
            let first = &a.major - &a.minor;
            let second = &b.major - &b.minor;
            let lower = if compare(&first, &second, dec)? == Sign::Positive {
                second
            } else {
                first
            };
            Divider::radial(
                base,
                a.axis.clone(),
                sq(&base_radius) + axis2 * sq(&h) * discriminant,
                positive_lower(&lower, dec)?,
            )
        } else {
            Divider::conic(base, a.axis.clone())
        });
    }
    Ok(plan)
}
