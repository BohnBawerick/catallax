//! Closed-circle integrals of the same surface forms used by quadrature.

use crate::integrate::{Builder, Expr, Program};
use crate::planar::real_direction as direction;
use catallax_numbers::Sign;
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
) -> Option<Program> {
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
                || !(&circle.radius_sq - &cylinder.radius.square())
                    .as_rational()?
                    .is_zero()
            {
                return None;
            }
            let axis = unit(&b, &cylinder.axis);
            let z = dot(offset(&b, &circle.center, &cylinder.origin), axis);
            let origin = offset(&b, &cylinder.origin, anchor);
            let radius = b.exact(&cylinder.radius);
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
            if !parallel(&cone.axis, &circle.normal)?
                || !circle
                    .center
                    .minus(&cone.apex)
                    .cross(&cone.axis)
                    .norm_squared()
                    .as_rational()?
                    .is_zero()
            {
                return None;
            }
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
