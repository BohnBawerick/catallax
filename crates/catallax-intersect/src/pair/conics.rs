//! Rational-carrier conics. Each non-rational result retains its two supports
//! and a regular witness; no rounded coordinate becomes a definition.

use super::{Divider, Plan, Spec};
use crate::{Contact, CurveForm, math::*};
use catallax_numbers::{Rational, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{Cone, Cylinder, Plane, Point3, RPoint3, RVector3, Sphere, Torus};

fn traced(detail: &str) -> Refusal {
    Refusal::new("intersect.traced-pair-required", detail)
}
fn algebraic_direction(detail: &str) -> Refusal {
    Refusal::new("intersect.algebraic-direction-required", detail)
}
fn coaxial(origin: &RPoint3, axis: &RVector3, other: &RPoint3, direction: &RVector3) -> bool {
    axis.cross(direction).is_zero() && other.minus(origin).cross(axis).is_zero()
}
fn circle(
    center: Point3,
    normal: RVector3,
    radius_sq: Real,
    ref_dir: RVector3,
    contact: Contact,
) -> Spec {
    let witness = if point_rational(&center).is_some() && radius_sq.as_rational().is_some() {
        None
    } else {
        Some(
            center.plus(
                &ref_dir
                    .to_real()
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
    normal: RVector3,
    radius_sq: Real,
    ref_dir: RVector3,
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

pub(super) fn plane_cone(p: &Plane, c: &Cone, sigma: Rational) -> Result<Plan, Refusal> {
    let axis2 = c.axis.norm_squared();
    let normal2 = p.normal.norm_squared();
    let dot = c.axis.dot(&p.normal);
    let offset = p.origin.minus(&c.apex).dot(&p.normal);
    let slope2 = sq(&c.slope);
    if c.axis.cross(&p.normal).is_zero() {
        let t = &offset / &dot;
        return Ok(match t.sign() {
            Sign::Negative => Plan::disjoint(),
            Sign::Zero => Plan::point(c.apex.to_real()),
            Sign::Positive => one_circle(
                c.apex.plus(&c.axis.scaled(&t)).to_real(),
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
    let k = &r(1) + &slope2;
    let quadratic = &axis2 - &(&k * &projected2);
    if offset.is_zero() {
        if quadratic.sign() == Sign::Positive {
            return Ok(Plan::point(c.apex.to_real()));
        }
        let factor = rational_sqrt(&(-&quadratic / (&axis2 * &normal2))).ok_or_else(|| {
            algebraic_direction("cone generators require algebraic direction fields")
        })?;
        let mut plan = Plan::empty();
        plan.points.push(c.apex.to_real());
        let tangent = quadratic.is_zero();
        let count = if tangent { 1 } else { 2 };
        for side in [r(1), r(-1)].into_iter().take(count) {
            let direction = projected.plus(&transverse.scaled(&(&factor * &side)));
            let gradient = direction
                .scaled(&axis2)
                .minus(&c.axis.scaled(&(&k * &direction.dot(&c.axis))));
            let orientation = &p.normal.cross(&gradient).dot(&direction) * &sigma;
            let natural = if orientation.sign() == Sign::Negative
                || (orientation.is_zero() && sigma.sign() == Sign::Negative)
            {
                direction.scaled(&r(-1))
            } else {
                direction.clone()
            };
            plan.curves.push(Spec {
                form: CurveForm::Line {
                    origin: c.apex.to_real(),
                    dir: natural,
                },
                witness: Some(c.apex.plus(&direction).to_real()),
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
    if quadratic.is_zero() {
        if beta.sign() == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        let plane_term = &(&axis2 * &sq(&offset)) / &normal2;
        let constant = &plane_term - &(&k * &sq(&beta));
        let linear = &k * &beta;
        let shift = &constant / &(r(2) * (&linear * &projected2));
        let vertex = foot.plus(&projected.scaled(&shift)).to_real();
        let focal = real(&(&linear / &(r(2) * axis2))) * root(&projected2);
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
    if quadratic.sign() == Sign::Positive {
        if beta.sign() == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        return Ok(Plan::form(
            CurveForm::Ellipse {
                center: center.to_real(),
                normal: p.normal.scaled(&sigma),
                major_dir: projected,
                major_sq: real(&major2),
                minor_sq: real(&(&scale / &axis2)),
            },
            Contact::Crossing,
        ));
    }
    let witness = center
        .to_real()
        .plus(&projected.to_real().scaled(&root(&(&major2 / &projected2))));
    Ok(form(
        CurveForm::Hyperbola {
            center: center.to_real(),
            normal: p.normal.scaled(&(-sigma)),
            axis_dir: projected,
            transverse_sq: real(&major2),
            conjugate_sq: real(&(-scale / axis2)),
            branch: StrictSign::Positive,
        },
        witness,
    ))
}

pub(super) fn cone_cylinder(c: &Cone, y: &Cylinder, sigma: Rational) -> Result<Plan, Refusal> {
    if !coaxial(&c.apex, &c.axis, &y.origin, &y.axis) {
        return Err(traced(
            "non-coaxial cone/cylinder intersection requires certified branch tracing",
        ));
    }
    let t = root(&(sq(&y.radius) / (sq(&c.slope) * c.axis.norm_squared())));
    Ok(one_circle(
        c.apex.to_real().plus(&c.axis.to_real().scaled(&t)),
        c.axis.scaled(&(-sigma)),
        real(&sq(&y.radius)),
        c.ref_dir.clone(),
        Contact::Crossing,
    ))
}

// Sign of a + side sqrt(discriminant), without approximating the root.
fn root_sign(a: &Rational, discriminant: &Rational, side: i64) -> Sign {
    if discriminant.is_zero() {
        return a.sign();
    }
    if side > 0 {
        if a.sign() != Sign::Negative {
            Sign::Positive
        } else {
            (discriminant - &sq(a)).sign()
        }
    } else if a.sign() != Sign::Positive {
        Sign::Negative
    } else {
        (&sq(a) - discriminant).sign()
    }
}

pub(super) fn cone_sphere(c: &Cone, s: &Sphere, sigma: Rational) -> Result<Plan, Refusal> {
    if !s.center.minus(&c.apex).cross(&c.axis).is_zero() {
        return Err(traced(
            "non-coaxial cone/sphere intersection requires certified branch tracing",
        ));
    }
    let axis2 = c.axis.norm_squared();
    let h = &s.center.minus(&c.apex).dot(&c.axis) / &axis2;
    let k = r(1) + sq(&c.slope);
    let discriminant = &(&k * &sq(&s.radius)) / &axis2 - sq(&c.slope) * sq(&h);
    if discriminant.sign() == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset = root(&discriminant);
    for side in [1, -1]
        .into_iter()
        .take(if discriminant.is_zero() { 1 } else { 2 })
    {
        match root_sign(&h, &discriminant, side) {
            Sign::Negative => continue,
            Sign::Zero => {
                plan.points.push(c.apex.to_real());
                continue;
            }
            Sign::Positive => {}
        }
        let t = (real(&h) + &Real::from_i64(side) * &offset) / real(&k);
        let center = c.apex.to_real().plus(&c.axis.to_real().scaled(&t));
        plan.curves.push(circle(
            center,
            c.axis.scaled(&(&r(-side) * &sigma)),
            real(&(&sq(&c.slope) * &axis2)) * t.square(),
            c.ref_dir.clone(),
            if discriminant.is_zero() {
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

pub(super) fn cones(a: &Cone, b: &Cone) -> Result<Plan, Refusal> {
    if !coaxial(&a.apex, &a.axis, &b.apex, &b.axis) {
        return Err(traced(
            "non-coaxial cone/cone intersection requires certified branch tracing",
        ));
    }
    let sense = if a.axis.dot(&b.axis).sign() == Sign::Positive {
        r(1)
    } else {
        r(-1)
    };
    let h = b.apex.minus(&a.apex).dot(&a.axis) / a.axis.norm_squared();
    let signed_slope = &sense * &b.slope;
    let denominator = &a.slope - &signed_slope;
    if denominator.is_zero() {
        return Ok(if h.is_zero() {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        });
    }
    let t = &(-(&signed_slope * &h)) / &denominator;
    if t.sign() == Sign::Negative || (&sense * &(&t - &h)).sign() == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    if t.is_zero() {
        return Ok(Plan::point(a.apex.to_real()));
    }
    let normal = a.axis.scaled(&if denominator.sign() == Sign::Positive {
        r(-1)
    } else {
        r(1)
    });
    Ok(one_circle(
        a.apex.plus(&a.axis.scaled(&t)).to_real(),
        normal,
        real(&(sq(&a.slope) * a.axis.norm_squared() * sq(&t))),
        a.ref_dir.clone(),
        Contact::Crossing,
    ))
}

pub(super) fn plane_torus(p: &Plane, t: &Torus, sigma: Rational) -> Result<Plan, Refusal> {
    let offset = p.origin.minus(&t.center).dot(&p.normal);
    let dot = p.normal.dot(&t.axis);
    if p.normal.cross(&t.axis).is_zero() {
        let height = &offset / &dot;
        let center = t.center.plus(&t.axis.scaled(&height));
        let discriminant = sq(&t.minor) - t.axis.norm_squared() * sq(&height);
        if discriminant.sign() == Sign::Negative {
            return Ok(Plan::disjoint());
        }
        if discriminant.is_zero() {
            return Ok(one_circle(
                center.to_real(),
                p.normal.scaled(&sigma),
                real(&sq(&t.major)),
                t.ref_dir.clone(),
                Contact::Tangent,
            ));
        }
        let mut plan = Plan::empty();
        let offset = root(&discriminant);
        for side in [r(1), r(-1)] {
            let radius = real(&t.major) + &real(&side) * &offset;
            plan.curves.push(circle(
                center.to_real(),
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
            &t.major - &t.minor,
        ));
        return Ok(plan);
    }
    if dot.is_zero() && offset.is_zero() {
        let radial = p.normal.cross(&t.axis);
        let shift = radial
            .to_real()
            .scaled(&(real(&t.major) / root(&radial.norm_squared())));
        let mut plan = Plan::empty();
        for side in [1, -1] {
            let center = t
                .center
                .to_real()
                .plus(&shift.scaled(&Real::from_i64(side)));
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
    let cross2 = r(4) * sq(&t.major) * sq(&t.minor) * (&radial2 * &normal2);
    if excess.sign() == Sign::Positive {
        match (sq(&excess) - cross2).sign() {
            Sign::Positive => return Ok(Plan::disjoint()),
            Sign::Zero => {
                let orientation = if offset.sign() == Sign::Positive {
                    Real::one()
                } else {
                    Real::from_i64(-1)
                };
                let displacement = radial_normal
                    .to_real()
                    .scaled(&root(&(sq(&t.major) / radial2)))
                    .plus(&p.normal.to_real().scaled(&root(&(sq(&t.minor) / normal2))));
                return Ok(Plan::point(
                    t.center.to_real().plus(&displacement.scaled(&orientation)),
                ));
            }
            Sign::Negative => {}
        }
    }
    Err(traced(
        "oblique or offset axial torus section requires certified branch tracing",
    ))
}

pub(super) fn cylinder_torus(c: &Cylinder, t: &Torus, sigma: Rational) -> Result<Plan, Refusal> {
    if !coaxial(&t.center, &t.axis, &c.origin, &c.axis) {
        return Err(traced(
            "non-coaxial cylinder/torus intersection requires certified branch tracing",
        ));
    }
    let discriminant = sq(&t.minor) - sq(&(&c.radius - &t.major));
    if discriminant.sign() == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset = root(&(&discriminant / &t.axis.norm_squared()));
    for side in [1, -1]
        .into_iter()
        .take(if discriminant.is_zero() { 1 } else { 2 })
    {
        let center = t
            .center
            .to_real()
            .plus(&t.axis.to_real().scaled(&(&Real::from_i64(side) * &offset)));
        plan.curves.push(circle(
            center,
            t.axis.scaled(&(&r(-side) * &sigma)),
            real(&sq(&c.radius)),
            c.ref_dir.clone(),
            if discriminant.is_zero() {
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

pub(super) fn sphere_torus(s: &Sphere, t: &Torus, sigma: Rational) -> Result<Plan, Refusal> {
    let displacement = s.center.minus(&t.center);
    let axis2 = t.axis.norm_squared();
    if !displacement.cross(&t.axis).is_zero() {
        let axial = displacement.dot(&t.axis);
        let radial2 = displacement.norm_squared() - &sq(&axial) / &axis2;
        let distance_term = displacement.norm_squared() + sq(&t.major);
        let root2 = r(4) * (&sq(&t.major) * &radial2);
        // Distance to the core circle is 1-Lipschitz. These strict tests
        // certify separated surfaces, including a sphere inside the tube.
        if (t.minor > s.radius
            && root_sign(&(&distance_term - &sq(&(&t.minor - &s.radius))), &root2, -1)
                == Sign::Negative)
            || root_sign(&(&distance_term - &sq(&(&t.minor + &s.radius))), &root2, -1)
                == Sign::Positive
            || (s.radius > &t.major + &t.minor
                && displacement.norm_squared() < sq(&(&(&s.radius - &t.major) - &t.minor)))
        {
            return Ok(Plan::disjoint());
        }
        if axial.is_zero() && radial2 == sq(&t.major) && s.radius == t.minor {
            return Ok(one_circle(
                s.center.to_real(),
                t.axis.cross(&displacement),
                real(&sq(&s.radius)),
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
    let offset = (sq(&s.radius) - sq(&t.minor) - height2 - sq(&t.major)) / r(2);
    let base_radius = &t.major + &(&(&t.major * &offset) / &denominator);
    let base_height = -(&(&h * &offset) / &denominator);
    let discriminant = &(sq(&t.minor) - &sq(&offset) / &denominator) / &denominator;
    if discriminant.sign() == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset_root = root(&(&discriminant / &axis2));
    for side in [1, -1]
        .into_iter()
        .take(if discriminant.is_zero() { 1 } else { 2 })
    {
        let root = &Real::from_i64(side) * &offset_root;
        let center = t.center.to_real().plus(
            &t.axis
                .to_real()
                .scaled(&(real(&base_height) + &real(&t.major) * &root)),
        );
        let radius = real(&base_radius) + real(&(&h * &axis2)) * root;
        plan.curves.push(circle(
            center,
            t.axis.scaled(&(&r(-side) * &sigma)),
            radius.square(),
            t.ref_dir.clone(),
            if discriminant.is_zero() {
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

pub(super) fn cone_torus(c: &Cone, t: &Torus, sigma: Rational) -> Result<Plan, Refusal> {
    if !coaxial(&c.apex, &c.axis, &t.center, &t.axis) {
        return Err(traced(
            "non-coaxial cone/torus intersection requires certified branch tracing",
        ));
    }
    let length = rational_sqrt(&c.axis.norm_squared()).ok_or_else(|| {
        algebraic_direction(
            "cone/torus meridian requires the shared algebraic coefficient foundation",
        )
    })?;
    let h = &t.center.minus(&c.apex).dot(&c.axis) / &length;
    let k = r(1) + sq(&c.slope);
    let middle = &(&c.slope * &t.major) + &h;
    let discriminant = sq(&middle) - &k * &(sq(&t.major) + sq(&h) - sq(&t.minor));
    if discriminant.sign() == Sign::Negative {
        return Ok(Plan::disjoint());
    }
    let mut plan = Plan::empty();
    let offset = root(&discriminant);
    for side in [1, -1]
        .into_iter()
        .take(if discriminant.is_zero() { 1 } else { 2 })
    {
        match root_sign(&middle, &discriminant, side) {
            Sign::Negative => continue,
            Sign::Zero => {
                plan.points.push(c.apex.to_real());
                continue;
            }
            Sign::Positive => {}
        }
        let height = (real(&middle) + &Real::from_i64(side) * &offset) / real(&k);
        let center = c
            .apex
            .to_real()
            .plus(&c.axis.to_real().scaled(&(&height / &real(&length))));
        plan.curves.push(circle(
            center,
            c.axis.scaled(&(&r(-side) * &sigma)),
            real(&sq(&c.slope)) * height.square(),
            c.ref_dir.clone(),
            if discriminant.is_zero() {
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

pub(super) fn tori(a: &Torus, b: &Torus) -> Result<Plan, Refusal> {
    if !coaxial(&a.center, &a.axis, &b.center, &b.axis) {
        return Err(traced(
            "non-coaxial torus/torus intersection requires certified branch tracing",
        ));
    }
    let axis2 = a.axis.norm_squared();
    let h = &b.center.minus(&a.center).dot(&a.axis) / &axis2;
    let radial = &b.major - &a.major;
    let distance2 = sq(&radial) + &axis2 * &sq(&h);
    if distance2.is_zero() {
        return Ok(if a.minor == b.minor {
            Plan::coincident(true)
        } else {
            Plan::disjoint()
        });
    }
    if distance2 > sq(&(&a.minor + &b.minor)) || distance2 < sq(&(&a.minor - &b.minor)) {
        return Ok(Plan::disjoint());
    }
    let fraction = (&distance2 + &sq(&a.minor) - sq(&b.minor)) / (&r(2) * &distance2);
    let base_radius = &a.major + &(&radial * &fraction);
    let base_height = &h * &fraction;
    let discriminant = &(sq(&a.minor) - &sq(&fraction) * &distance2) / &distance2;
    let mut plan = Plan::empty();
    let offset_root = root(&(&discriminant / &axis2));
    let first =
        if radial.sign() == Sign::Negative || (radial.is_zero() && h.sign() == Sign::Positive) {
            -1
        } else {
            1
        };
    for side in [first, -first]
        .into_iter()
        .take(if discriminant.is_zero() { 1 } else { 2 })
    {
        let root = &Real::from_i64(side) * &offset_root;
        let center = a.center.to_real().plus(
            &a.axis
                .to_real()
                .scaled(&(real(&base_height) + &real(&radial) * &root)),
        );
        let radius = real(&base_radius) - real(&(&h * &axis2)) * root;
        plan.curves.push(circle(
            center,
            a.axis.scaled(&r(-side)),
            radius.square(),
            a.ref_dir.clone(),
            if discriminant.is_zero() {
                Contact::Tangent
            } else {
                Contact::Crossing
            },
        ));
    }
    if plan.curves.len() == 2 {
        let base = a.center.plus(&a.axis.scaled(&base_height));
        plan.divider = Some(if radial.is_zero() {
            Divider::radial(
                base,
                a.axis.clone(),
                sq(&base_radius) + axis2 * sq(&h) * discriminant,
                (&a.major - &a.minor).min(&b.major - &b.minor),
            )
        } else {
            Divider::conic(base, a.axis.clone())
        });
    }
    Ok(plan)
}
