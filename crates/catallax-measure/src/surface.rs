//! Boundary forms for analytic surfaces. No angle coordinate enters a decision.

use catallax_numbers::Real;
use catallax_shapes::{Point3, RPoint3, Surface, Vector3};

use crate::integrate::{Builder, Expr, Program};
use crate::planar::{real_components as components, real_direction as direction};

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Property {
    Area,
    Volume,
    Centroid,
    /// Raw second moment in the body's translated frame.
    SecondMoment(usize, usize),
}

pub(crate) struct Compiled {
    pub(crate) boundary: Program,
    pub(crate) full: Option<Program>,
    pub(crate) indices: Vec<usize>,
    pub(crate) closed: bool,
    pub(crate) sphere: bool,
}

type Vector<'a> = [Expr<'a>; 3];

fn dot<'a>(a: Vector<'a>, b: Vector<'a>) -> Expr<'a> {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross<'a>(a: Vector<'a>, b: Vector<'a>) -> Vector<'a> {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn constants<'a>(b: &'a Builder, v: [&Real; 3]) -> Vector<'a> {
    v.map(|r| b.exact(r))
}
fn unit<'a>(b: &'a Builder, v: &Vector3) -> Vector<'a> {
    let norm = b.exact(&v.norm_squared()).sqrt();
    constants(b, direction(v)).map(|x| x / norm)
}
fn offset<'a>(b: &'a Builder, p: &Point3, anchor: &RPoint3) -> Vector<'a> {
    let p = p.minus(&anchor.to_real());
    constants(b, direction(&p))
}
fn displacement<'a>(b: &'a Builder, p: &Point3) -> Vector<'a> {
    let o = constants(b, components(p));
    std::array::from_fn(|i| b.input(i) - o[i])
}

pub(crate) fn compile(
    surface: &Surface,
    anchor: &RPoint3,
    property: Property,
    closed: bool,
    gauge: &Vector3,
) -> Option<Compiled> {
    let b = Builder::default();
    let derivative = std::array::from_fn(|i| b.input(i + 3));
    let zero = b.k(0);
    let mut full = None;
    let mut sphere = false;
    if let Property::SecondMoment(i, j) = property {
        return second_moment(surface, anchor, closed, gauge, i, j, &b);
    }
    let expressions = match surface {
        Surface::Plane(p) => {
            let q = displacement(&b, &p.origin);
            let first = constants(&b, direction(&p.ref_dir));
            let second_vector = p.normal.cross(&p.ref_dir);
            let second = constants(&b, direction(&second_vector));
            let reference_sq = b.exact(&p.ref_dir.norm_squared());
            let u = dot(q, first) / reference_sq;
            let v = dot(q, second) / b.exact(&second_vector.norm_squared());
            let du = dot(derivative, first) / reference_sq;
            let normal = constants(&b, direction(&p.normal));
            let origin = offset(&b, &p.origin, anchor);
            let area = -b.exact(&p.normal.norm_squared()).sqrt() * reference_sq * v * du;
            let volume = -reference_sq * dot(normal, origin) * v * du / b.k(3);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                let base = origin[i] + u * first[i];
                let primitive = base.square() * v
                    + base * second[i] * v.square()
                    + second[i].square() * v.square() * v / b.k(3);
                -reference_sq * normal[i] * primitive * du / b.k(2)
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::Cylinder(c) => {
            let q = displacement(&b, &c.origin);
            let axis = unit(&b, &c.axis);
            let z = dot(q, axis);
            let radial = std::array::from_fn::<_, 3, _>(|i| q[i] - z * axis[i]);
            let radius = b.exact(&c.radius);
            let er = radial.map(|x| x / radius);
            let du = dot(axis, cross(q, derivative)) / radius.square();
            let origin = offset(&b, &c.origin, anchor);
            let area = -radius * z * du;
            let volume = -radius * (radius + dot(origin, er)) * z * du / b.k(3);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                let base = origin[i] + radius * er[i];
                let primitive = base.square() * z
                    + base * axis[i] * z.square()
                    + axis[i].square() * z.square() * z / b.k(3);
                -radius * er[i] * primitive * du / b.k(2)
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::EllipticCylinder(c) => {
            let q = displacement(&b, &c.origin);
            let axis = unit(&b, &c.axis);
            let first = unit(&b, &c.major_dir);
            let second = cross(axis, first);
            let major = b.exact(&c.major_sq).sqrt();
            let minor = b.exact(&c.minor_sq).sqrt();
            let x = dot(q, first);
            let y = dot(q, second);
            let z = dot(q, axis);
            let cosine = x / major;
            let sine = y / minor;
            let du = (x * dot(derivative, second) - y * dot(derivative, first)) / (major * minor);
            let density =
                (major.square() * sine.square() + minor.square() * cosine.square()).sqrt();
            let normal_density = std::array::from_fn::<_, 3, _>(|i| {
                minor * cosine * first[i] + major * sine * second[i]
            });
            let origin = offset(&b, &c.origin, anchor);
            let radial = std::array::from_fn::<_, 3, _>(|i| x * first[i] + y * second[i]);
            let area = -z * density * du;
            let volume = -z * (major * minor + dot(normal_density, origin)) * du / b.k(3);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                let base = origin[i] + radial[i];
                let primitive = base.square() * z
                    + base * axis[i] * z.square()
                    + axis[i].square() * z.square() * z / b.k(3);
                -normal_density[i] * primitive * du / b.k(2)
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::Cone(c) => {
            let q = displacement(&b, &c.apex);
            let axis = unit(&b, &c.axis);
            let z = dot(q, axis);
            let slope = b.exact(&c.slope);
            let radial = std::array::from_fn::<_, 3, _>(|i| q[i] - z * axis[i]);
            let er = radial.map(|x| x / (slope * z));
            let angular = dot(axis, cross(q, derivative));
            let z2du = angular / slope.square();
            let origin = offset(&b, &c.apex, anchor);
            let normal = std::array::from_fn::<_, 3, _>(|i| er[i] - slope * axis[i]);
            let length = (b.k(1) + slope.square()).sqrt();
            let area = -slope * length * z2du / b.k(2);
            let volume = -slope * dot(origin, normal) * z2du / b.k(6);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                let along = slope * er[i] + axis[i];
                let primitive = origin[i].square() / b.k(2)
                    + b.k(2) * origin[i] * along * z / b.k(3)
                    + along.square() * z.square() / b.k(4);
                -slope * normal[i] * primitive * z2du / b.k(2)
            });
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::Sphere(s) => {
            sphere = true;
            let q = displacement(&b, &s.center);
            let radius = b.exact(&s.radius);
            let pole = unit(&b, gauge);
            let area = radius * dot(pole, cross(q, derivative)) / (radius + dot(pole, q));
            let origin = offset(&b, &s.center, anchor);
            // Integrate spherical harmonics of degrees 1, 2 and 3 by their
            // Laplace-Beltrami eigenvalues. Only the scalar area needs a gauge.
            let boundary = cross(derivative, q);
            let volume = radius * area / b.k(3) - dot(origin, boundary) / b.k(6);
            let moments = std::array::from_fn::<_, 3, _>(|i| {
                origin[i] * radius * area / b.k(3)
                    - (radius.square() / b.k(2) + origin[i].square()) * boundary[i] / b.k(4)
                    - origin[i] * q[i] * boundary[i] / b.k(3)
                    - q[i].square() * boundary[i] / b.k(8)
            });
            let whole_area = b.k(4) * b.input(6) * radius.square();
            let whole_volume = whole_area * radius / b.k(3);
            full = Some([
                whole_area,
                whole_volume,
                whole_volume * origin[0],
                whole_volume * origin[1],
                whole_volume * origin[2],
            ]);
            [area, volume, moments[0], moments[1], moments[2]]
        }
        Surface::Torus(t) if closed => {
            let major = b.exact(&t.major);
            let minor = b.exact(&t.minor);
            let pi2 = b.input(6).square();
            let area = b.k(4) * pi2 * major * minor;
            let volume = b.k(2) * pi2 * major * minor.square();
            let origin = offset(&b, &t.center, anchor);
            full = Some([
                area,
                volume,
                volume * origin[0],
                volume * origin[1],
                volume * origin[2],
            ]);
            [zero; 5]
        }
        _ => return None,
    };
    let mut indices = match property {
        Property::Area => vec![0],
        Property::Volume => vec![1],
        Property::Centroid => vec![1, 2, 3, 4],
        Property::SecondMoment(_, _) => unreachable!(),
    };
    if sphere && !closed && property != Property::Area {
        indices.insert(0, 0);
    }
    let boundary = b.finish(&indices.iter().map(|i| expressions[*i]).collect::<Vec<_>>());
    let full = full
        .map(|expressions| b.finish(&indices.iter().map(|i| expressions[*i]).collect::<Vec<_>>()));
    Some(Compiled {
        boundary,
        full,
        indices,
        closed,
        sphere,
    })
}

/// Radial divergence field `F = p p_i p_j / 5`. Its divergence is `p_i p_j`.
/// All forms use the same translated frame as volume and first moments.
fn second_moment<'a>(
    surface: &Surface,
    anchor: &RPoint3,
    closed: bool,
    gauge: &Vector3,
    i: usize,
    j: usize,
    b: &'a Builder,
) -> Option<Compiled> {
    let dq = std::array::from_fn(|k| b.input(k + 3));
    let primitive = |base: Vector<'a>, along: Vector<'a>, t: Expr<'a>| {
        base[i] * base[j] * t
            + (base[i] * along[j] + base[j] * along[i]) * t.square() / b.k(2)
            + along[i] * along[j] * t.square() * t / b.k(3)
    };
    let mut full = None;
    let mut sphere = false;
    let mut area = b.k(0);
    let moment = match surface {
        Surface::Plane(p) => {
            let q = displacement(b, &p.origin);
            let first = constants(b, direction(&p.ref_dir));
            let second_vector = p.normal.cross(&p.ref_dir);
            let second = constants(b, direction(&second_vector));
            let reference_sq = b.exact(&p.ref_dir.norm_squared());
            let u = dot(q, first) / reference_sq;
            let v = dot(q, second) / b.exact(&second_vector.norm_squared());
            let du = dot(dq, first) / reference_sq;
            let origin = offset(b, &p.origin, anchor);
            let base = std::array::from_fn(|k| origin[k] + u * first[k]);
            -reference_sq
                * dot(constants(b, direction(&p.normal)), origin)
                * primitive(base, second, v)
                * du
                / b.k(5)
        }
        Surface::Cylinder(c) => {
            let q = displacement(b, &c.origin);
            let axis = unit(b, &c.axis);
            let z = dot(q, axis);
            let radius = b.exact(&c.radius);
            let er = std::array::from_fn(|k| (q[k] - z * axis[k]) / radius);
            let du = dot(axis, cross(q, dq)) / radius.square();
            let origin = offset(b, &c.origin, anchor);
            let base = std::array::from_fn(|k| origin[k] + radius * er[k]);
            -radius * (radius + dot(origin, er)) * primitive(base, axis, z) * du / b.k(5)
        }
        Surface::EllipticCylinder(c) => {
            let q = displacement(b, &c.origin);
            let axis = unit(b, &c.axis);
            let first = unit(b, &c.major_dir);
            let second = cross(axis, first);
            let major = b.exact(&c.major_sq).sqrt();
            let minor = b.exact(&c.minor_sq).sqrt();
            let x = dot(q, first);
            let y = dot(q, second);
            let z = dot(q, axis);
            let normal = std::array::from_fn(|k| {
                minor * x / major * first[k] + major * y / minor * second[k]
            });
            let origin = offset(b, &c.origin, anchor);
            let base = std::array::from_fn(|k| origin[k] + x * first[k] + y * second[k]);
            let du = (x * dot(dq, second) - y * dot(dq, first)) / (major * minor);
            -(major * minor + dot(normal, origin)) * primitive(base, axis, z) * du / b.k(5)
        }
        Surface::Cone(c) => {
            let q = displacement(b, &c.apex);
            let axis = unit(b, &c.axis);
            let z = dot(q, axis);
            let slope = b.exact(&c.slope);
            let er: Vector<'_> = std::array::from_fn(|k| (q[k] - z * axis[k]) / (slope * z));
            let origin = offset(b, &c.apex, anchor);
            let along: Vector<'_> = std::array::from_fn(|k| slope * er[k] + axis[k]);
            let normal = std::array::from_fn(|k| er[k] - slope * axis[k]);
            let polynomial = origin[i] * origin[j] / b.k(2)
                + (origin[i] * along[j] + origin[j] * along[i]) * z / b.k(3)
                + along[i] * along[j] * z.square() / b.k(4);
            -slope * dot(origin, normal) * polynomial * dot(axis, cross(q, dq))
                / slope.square()
                / b.k(5)
        }
        Surface::Sphere(s) => {
            sphere = true;
            let q = displacement(b, &s.center);
            let r = b.exact(&s.radius);
            let origin = offset(b, &s.center, anchor);
            let pole = unit(b, gauge);
            area = r * dot(pole, cross(q, dq)) / (r + dot(pole, q));
            let boundary = cross(dq, q);
            let delta = |a, c| b.k(if a == c { 1 } else { 0 });
            let first = |a: usize| -r * boundary[a] / b.k(2);
            let second = |a: usize, c: usize| {
                delta(a, c) * r.square() * area / b.k(3)
                    - r * (q[a] * boundary[c] + q[c] * boundary[a]) / b.k(6)
            };
            let third = |a: usize, c: usize, d: usize| {
                -r * (q[a] * q[c] * boundary[d]
                    + q[a] * q[d] * boundary[c]
                    + q[c] * q[d] * boundary[a])
                    / b.k(12)
                    - r.square()
                        * r
                        * (delta(a, c) * boundary[d]
                            + delta(a, d) * boundary[c]
                            + delta(c, d) * boundary[a])
                        / b.k(12)
            };
            let mut integral = r.square()
                * (origin[i] * origin[j] * area
                    + origin[i] * first(j)
                    + origin[j] * first(i)
                    + second(i, j));
            for a in 0..3 {
                integral = integral
                    + origin[a]
                        * (origin[i] * origin[j] * first(a)
                            + origin[i] * second(a, j)
                            + origin[j] * second(a, i)
                            + third(a, i, j));
            }
            let whole_area = b.k(4) * b.input(6) * r.square();
            let volume = whole_area * r / b.k(3);
            full = Some([
                whole_area,
                volume * (origin[i] * origin[j] + delta(i, j) * r.square() / b.k(5)),
            ]);
            integral / (b.k(5) * r)
        }
        Surface::Torus(t) if closed => {
            let major = b.exact(&t.major);
            let minor = b.exact(&t.minor);
            let axis = unit(b, &t.axis);
            let origin = offset(b, &t.center, anchor);
            let volume = b.k(2) * b.input(6).square() * major * minor.square();
            let transverse = major.square() / b.k(2) + b.k(3) * minor.square() / b.k(8);
            let axial = minor.square() / b.k(4);
            let delta = b.k(if i == j { 1 } else { 0 });
            full = Some([
                b.k(0),
                volume
                    * (origin[i] * origin[j]
                        + transverse * (delta - axis[i] * axis[j])
                        + axial * axis[i] * axis[j]),
            ]);
            b.k(0)
        }
        _ => return None,
    };
    let indices = if sphere && !closed {
        vec![0, 2]
    } else {
        vec![2]
    };
    let expressions = [area, moment];
    let select = |values: [Expr<'_>; 2]| {
        b.finish(
            &indices
                .iter()
                .map(|index| values[usize::from(*index == 2)])
                .collect::<Vec<_>>(),
        )
    };
    Some(Compiled {
        boundary: select(expressions),
        full: full.map(select),
        indices,
        closed,
        sphere,
    })
}
