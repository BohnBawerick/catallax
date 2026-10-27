//! Closed-circle integrals of the same surface forms used by quadrature.

use crate::integrate::{Builder, Expr, Program};
use crate::planar::real_direction as direction;
use catallax_numbers::{Decider, Outcome, Refusal, Sign};
use catallax_shapes::{Circle, Point3, RPoint3, Surface, Vector3};

fn dot<'a>(a: [Expr<'a>; 3], b: [Expr<'a>; 3]) -> Expr<'a> {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn unit<'a>(b: &'a Builder, v: &Vector3) -> [Expr<'a>; 3] {
    let length = b.exact(&v.norm_squared()).sqrt();
    direction(v).map(|c| b.exact(c) / length)
}
fn offset<'a>(b: &'a Builder, point: &Point3, anchor: &Point3) -> [Expr<'a>; 3] {
    direction(&point.minus(anchor)).map(|c| b.exact(c))
}
fn parallel(a: &Vector3, b: &Vector3) -> Option<bool> {
    Some(a.cross(b).norm_squared().as_rational()?.is_zero())
}
fn orientation<'a>(b: &'a Builder, normal: &Vector3, axis: &Vector3) -> Option<Expr<'a>> {
    Some(b.k(
        if normal.dot(axis).as_rational()?.sign() == Sign::Positive {
            1
        } else {
            -1
        },
    ))
}

pub(crate) fn compile(
    surface: &Surface,
    circle: &Circle,
    anchor: &RPoint3,
    indices: &[usize],
    property: crate::surface::Property,
    dec: &mut Decider,
) -> Result<Option<Program>, Refusal> {
    if let Surface::Cone(cone) = surface {
        if parallel(&cone.axis, &circle.normal) != Some(true) {
            return Ok(None);
        }
        let coaxial = circle
            .center
            .minus(&cone.apex)
            .cross(&cone.axis)
            .norm_squared();
        let proven = if let Some(value) = coaxial.as_rational() {
            value.is_zero()
        } else {
            match dec.sign(&coaxial) {
                Outcome::ProvenCoincident(_) => true,
                Outcome::Proven(_) | Outcome::NotExact(_) => false,
                Outcome::Refused(refusal) => return Err(refusal),
            }
        };
        if !proven {
            return Ok(None);
        }
    }
    Ok(compile_admitted(surface, circle, anchor, indices, property))
}

// Conical circle coaxiality and parallel normals are proved before any integral family.
fn compile_admitted(
    surface: &Surface,
    circle: &Circle,
    anchor: &RPoint3,
    indices: &[usize],
    property: crate::surface::Property,
) -> Option<Program> {
    if let crate::surface::Property::SecondMoment(i, j) = property {
        return second_moment(surface, circle, anchor, indices, i, j);
    }
    if property == crate::surface::Property::FaceCentroid {
        return face_centroid(surface, circle, anchor, indices);
    }
    let anchor = &anchor.to_real();
    let b = Builder::default();
    let pi = b.input(6);
    let radius2 = b.exact(&circle.radius_sq);
    let expressions = match surface {
        Surface::Plane(plane) => {
            if !parallel(&plane.normal, &circle.normal)? {
                return None;
            }
            let normal = unit(&b, &plane.normal);
            let center = offset(&b, &circle.center, anchor);
            let area = orientation(&b, &circle.normal, &plane.normal)? * pi * radius2;
            let volume = dot(center, normal) * area / b.k(3);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                normal[i]
                    * area
                    * (center[i].square() + radius2 * (b.k(1) - normal[i].square()) / b.k(4))
                    / b.k(2)
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::Cylinder(cylinder) => {
            if !parallel(&cylinder.axis, &circle.normal)?
                || !circle
                    .center
                    .minus(&cylinder.origin)
                    .cross(&cylinder.axis)
                    .norm_squared()
                    .as_rational()?
                    .is_zero()
                || !(&circle.radius_sq - &cylinder.radius_sq)
                    .as_rational()?
                    .is_zero()
            {
                return None;
            }
            let axis = unit(&b, &cylinder.axis);
            let z = dot(offset(&b, &circle.center, &cylinder.origin), axis);
            let origin = offset(&b, &cylinder.origin, anchor);
            let radius = b.exact(&cylinder.radius());
            let signed_pi = orientation(&b, &circle.normal, &cylinder.axis)? * pi;
            let area = -b.k(2) * signed_pi * radius * z;
            let volume = -b.k(2) * signed_pi * radius.square() * z / b.k(3);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                -signed_pi
                    * radius.square()
                    * (b.k(1) - axis[i].square())
                    * (origin[i] * z + axis[i] * z.square() / b.k(2))
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::Cone(cone) => {
            let axis = unit(&b, &cone.axis);
            let z = dot(offset(&b, &circle.center, &cone.apex), axis);
            let origin = offset(&b, &cone.apex, anchor);
            let slope = b.exact(&cone.slope);
            let signed_pi = orientation(&b, &circle.normal, &cone.axis)? * pi;
            let area = -signed_pi * slope * (b.k(1) + slope.square()).sqrt() * z.square();
            let volume = signed_pi * slope.square() * z.square() * dot(origin, axis) / b.k(3);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                let average_square = (b.k(1) - axis[i].square()) / b.k(2);
                let average = -slope * axis[i] * origin[i].square() / b.k(2)
                    + b.k(2) * origin[i] * z * slope * (average_square - axis[i].square()) / b.k(3)
                    + z.square()
                        * slope
                        * axis[i]
                        * ((b.k(2) - slope.square()) * average_square - axis[i].square())
                        / b.k(4);
                -slope * z.square() * signed_pi * average
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::Sphere(sphere) => {
            let axis = unit(&b, &circle.normal);
            let d = dot(offset(&b, &circle.center, &sphere.center), axis);
            let radius = b.exact(&sphere.radius);
            let origin = offset(&b, &sphere.center, anchor);
            let area = b.k(2) * pi * radius * (radius - d);
            let first = axis.map(|c| pi * radius * radius2 * c);
            let volume = (radius * area + dot(origin, first) / radius) / b.k(3);
            let second_integral = (radius.square() * radius - d.square() * d) / b.k(3);
            let third_integral = (radius.square().square() - d.square().square()) / b.k(4);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                let perpendicular = b.k(1) - axis[i].square();
                let second = b.k(2)
                    * pi
                    * radius
                    * (axis[i].square() * second_integral
                        + perpendicular * (radius.square() * (radius - d) - second_integral)
                            / b.k(2));
                let third = b.k(2)
                    * pi
                    * radius
                    * (axis[i].square() * axis[i] * third_integral
                        + b.k(3)
                            * axis[i]
                            * perpendicular
                            * (radius.square() * (radius.square() - d.square()) / b.k(2)
                                - third_integral)
                            / b.k(2));
                (third + b.k(2) * origin[i] * second + origin[i].square() * first[i])
                    / (b.k(2) * radius)
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        _ => return None,
    };
    Some(b.finish(&indices.iter().map(|i| expressions[*i]).collect::<Vec<_>>()))
}

fn second_moment(
    surface: &Surface,
    circle: &Circle,
    anchor: &RPoint3,
    indices: &[usize],
    i: usize,
    j: usize,
) -> Option<Program> {
    let b = Builder::default();
    let pi = b.input(6);
    let axis_vector = match surface {
        Surface::Plane(p) => &p.normal,
        Surface::Cylinder(c) => &c.axis,
        Surface::Cone(c) => &c.axis,
        Surface::Sphere(_) => &circle.normal,
        _ => return None,
    };
    if !matches!(surface, Surface::Cone(_)) && !parallel(axis_vector, &circle.normal)? {
        return None;
    }
    let axis = unit(&b, axis_vector);
    let delta = |a, c| b.k(if a == c { 1 } else { 0 });
    let perpendicular = |a: usize, c: usize| delta(a, c) - axis[a] * axis[c];
    let anchor = anchor.to_real();
    let radius2 = b.exact(&circle.radius_sq);
    let moment = match surface {
        Surface::Plane(p) => {
            let center = offset(&b, &circle.center, &anchor);
            let area = orientation(&b, &circle.normal, &p.normal)? * pi * radius2;
            dot(center, axis)
                * area
                * (center[i] * center[j] + radius2 * perpendicular(i, j) / b.k(4))
                / b.k(5)
        }
        Surface::Cylinder(c) => {
            if !circle
                .center
                .minus(&c.origin)
                .cross(&c.axis)
                .norm_squared()
                .as_rational()?
                .is_zero()
                || !(&circle.radius_sq - &c.radius_sq).as_rational()?.is_zero()
            {
                return None;
            }
            let z = dot(offset(&b, &circle.center, &c.origin), axis);
            let o = offset(&b, &c.origin, &anchor);
            let r = b.exact(&c.radius());
            let axial = dot(o, axis);
            let po: [Expr<'_>; 3] = std::array::from_fn(|a| o[a] - axial * axis[a]);
            let second = r * o[i] * o[j]
                + r * r.square() * perpendicular(i, j) / b.k(2)
                + r * (o[i] * po[j] + o[j] * po[i]) / b.k(2);
            let first = |a: usize| r * o[a] + r * po[a] / b.k(2);
            -b.k(2)
                * orientation(&b, &circle.normal, &c.axis)?
                * pi
                * r
                * (second * z
                    + (axis[j] * first(i) + axis[i] * first(j)) * z.square() / b.k(2)
                    + axis[i] * axis[j] * r * z.square() * z / b.k(3))
                / b.k(5)
        }
        Surface::Cone(c) => {
            let z = dot(offset(&b, &circle.center, &c.apex), axis);
            let o = offset(&b, &c.apex, &anchor);
            let slope = b.exact(&c.slope);
            let axial = dot(o, axis);
            let po: [Expr<'_>; 3] = std::array::from_fn(|a| o[a] - axial * axis[a]);
            let zero = -slope * axial;
            let first = |a: usize| slope * po[a] / b.k(2) - slope * axial * axis[a];
            let second = slope * (po[i] * axis[j] + po[j] * axis[i]) / b.k(2)
                - slope
                    * axial
                    * (slope.square() * perpendicular(i, j) / b.k(2) + axis[i] * axis[j]);
            -b.k(2)
                * orientation(&b, &circle.normal, &c.axis)?
                * pi
                * slope
                * z.square()
                * (zero * o[i] * o[j] / b.k(2)
                    + (o[i] * first(j) + o[j] * first(i)) * z / b.k(3)
                    + second * z.square() / b.k(4))
                / b.k(5)
        }
        Surface::Sphere(s) => {
            let d = dot(offset(&b, &circle.center, &s.center), axis);
            let r = b.exact(&s.radius);
            let o = offset(&b, &s.center, &anchor);
            let area = b.k(2) * pi * r * (r - d);
            let t2 = (r.square() * r - d.square() * d) / b.k(3);
            let t3 = (r.square().square() - d.square().square()) / b.k(4);
            let first = |a: usize| pi * r * radius2 * axis[a];
            let second = |a: usize, c: usize| {
                b.k(2)
                    * pi
                    * r
                    * (axis[a] * axis[c] * t2
                        + perpendicular(a, c) * (r.square() * (r - d) - t2) / b.k(2))
            };
            let third = |a: usize, c: usize, e: usize| {
                b.k(2)
                    * pi
                    * r
                    * (axis[a] * axis[c] * axis[e] * t3
                        + (axis[a] * perpendicular(c, e)
                            + axis[c] * perpendicular(a, e)
                            + axis[e] * perpendicular(a, c))
                            * (r.square() * (r.square() - d.square()) / b.k(2) - t3)
                            / b.k(2))
            };
            let mut value = r.square()
                * (o[i] * o[j] * area + o[i] * first(j) + o[j] * first(i) + second(i, j));
            for a in 0..3 {
                value = value
                    + o[a]
                        * (o[i] * o[j] * first(a)
                            + o[i] * second(a, j)
                            + o[j] * second(a, i)
                            + third(a, i, j));
            }
            return Some(
                b.finish(
                    &indices
                        .iter()
                        .map(|index| {
                            if *index == 0 {
                                area
                            } else {
                                value / (b.k(5) * r)
                            }
                        })
                        .collect::<Vec<_>>(),
                ),
            );
        }
        _ => return None,
    };
    Some(b.finish(&[moment]))
}

fn face_centroid(
    surface: &Surface,
    circle: &Circle,
    anchor: &RPoint3,
    indices: &[usize],
) -> Option<Program> {
    let b = Builder::default();
    let pi = b.input(6);
    let anchor = anchor.to_real();
    let radius2 = b.exact(&circle.radius_sq);
    let (area, moments) = match surface {
        Surface::Plane(p) => {
            if !parallel(&p.normal, &circle.normal)? {
                return None;
            }
            let area = orientation(&b, &circle.normal, &p.normal)? * pi * radius2;
            let center = offset(&b, &circle.center, &anchor);
            (area, center.map(|x| area * x))
        }
        Surface::Cylinder(c) => {
            if !parallel(&c.axis, &circle.normal)?
                || !circle
                    .center
                    .minus(&c.origin)
                    .cross(&c.axis)
                    .norm_squared()
                    .as_rational()?
                    .is_zero()
                || !(&circle.radius_sq - &c.radius_sq).as_rational()?.is_zero()
            {
                return None;
            }
            let axis = unit(&b, &c.axis);
            let z = dot(offset(&b, &circle.center, &c.origin), axis);
            let origin = offset(&b, &c.origin, &anchor);
            let density =
                -b.k(2) * orientation(&b, &circle.normal, &c.axis)? * pi * b.exact(&c.radius());
            (
                density * z,
                std::array::from_fn(|i| density * (origin[i] * z + axis[i] * z.square() / b.k(2))),
            )
        }
        Surface::Cone(c) => {
            let axis = unit(&b, &c.axis);
            let z = dot(offset(&b, &circle.center, &c.apex), axis);
            let origin = offset(&b, &c.apex, &anchor);
            let slope = b.exact(&c.slope);
            let area = -orientation(&b, &circle.normal, &c.axis)?
                * pi
                * slope
                * (b.k(1) + slope.square()).sqrt()
                * z.square();
            (
                area,
                std::array::from_fn(|i| area * (origin[i] + b.k(2) * axis[i] * z / b.k(3))),
            )
        }
        Surface::Sphere(s) => {
            let axis = unit(&b, &circle.normal);
            let d = dot(offset(&b, &circle.center, &s.center), axis);
            let radius = b.exact(&s.radius);
            let area = b.k(2) * pi * radius * (radius - d);
            let origin = offset(&b, &s.center, &anchor);
            (
                area,
                std::array::from_fn(|i| area * origin[i] + pi * radius * radius2 * axis[i]),
            )
        }
        _ => return None,
    };
    let values = [area, area, moments[0], moments[1], moments[2]];
    Some(b.finish(&indices.iter().map(|i| values[*i]).collect::<Vec<_>>()))
}
