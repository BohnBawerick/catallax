//! Boundary forms for analytic surfaces. No angle coordinate enters a decision.

use catallax_numbers::Rational;
use catallax_shapes::{RPoint3, RVector3, Surface};

use crate::integrate::{Builder, Expr, Program};
use crate::planar::{components, direction};

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Property {
    Area,
    Volume,
    Centroid,
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
fn constants<'a>(b: &'a Builder, v: [&Rational; 3]) -> Vector<'a> {
    v.map(|r| b.r(r.clone()))
}
fn unit<'a>(b: &'a Builder, v: &RVector3) -> Vector<'a> {
    let norm = b.r(v.norm_squared()).sqrt();
    constants(b, direction(v)).map(|x| x / norm)
}
fn offset<'a>(b: &'a Builder, p: &RPoint3, anchor: &RPoint3) -> Vector<'a> {
    let p = p.minus(anchor);
    constants(b, direction(&p))
}
fn displacement<'a>(b: &'a Builder, p: &RPoint3) -> Vector<'a> {
    let o = constants(b, components(p));
    std::array::from_fn(|i| b.input(i) - o[i])
}

pub(crate) fn compile(
    surface: &Surface,
    anchor: &RPoint3,
    property: Property,
    closed: bool,
    gauge: &RVector3,
) -> Option<Compiled> {
    let b = Builder::default();
    let derivative = std::array::from_fn(|i| b.input(i + 3));
    let zero = b.k(0);
    let mut full = None;
    let mut sphere = false;
    let expressions = match surface {
        Surface::Plane(p) => {
            let q = displacement(&b, &p.origin);
            let first = constants(&b, direction(&p.ref_dir));
            let second_vector = p.normal.cross(&p.ref_dir);
            let second = constants(&b, direction(&second_vector));
            let reference_sq = b.r(p.ref_dir.norm_squared());
            let u = dot(q, first) / reference_sq;
            let v = dot(q, second) / b.r(second_vector.norm_squared());
            let du = dot(derivative, first) / reference_sq;
            let normal = constants(&b, direction(&p.normal));
            let origin = offset(&b, &p.origin, anchor);
            let area = -b.r(p.normal.norm_squared()).sqrt() * reference_sq * v * du;
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
            let radius = b.r(c.radius.clone());
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
        Surface::Cone(c) => {
            let q = displacement(&b, &c.apex);
            let axis = unit(&b, &c.axis);
            let z = dot(q, axis);
            let slope = b.r(c.slope.clone());
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
            let radius = b.r(s.radius.clone());
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
            let major = b.r(t.major.clone());
            let minor = b.r(t.minor.clone());
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
