//! Positive construction identities only. A failed match says nothing about geometry.

use crate::{Point3, Surface, Vector3};

pub(crate) fn same_point(a: &Point3, b: &Point3) -> bool {
    a.x.same_definition(&b.x) && a.y.same_definition(&b.y) && a.z.same_definition(&b.z)
}

fn same_vector(a: &Vector3, b: &Vector3) -> bool {
    a.x.same_definition(&b.x) && a.y.same_definition(&b.y) && a.z.same_definition(&b.z)
}

// Negation is a definition operation, not a sign decision or a tolerance test.
pub(crate) fn parallel(a: &Vector3, b: &Vector3) -> bool {
    same_vector(a, b)
        || (a.x.same_definition(&(-&b.x))
            && a.y.same_definition(&(-&b.y))
            && a.z.same_definition(&(-&b.z)))
}

pub(crate) fn on_plane(origin: &Point3, normal: &Vector3, p: &Point3) -> bool {
    [
        (&normal.x, &origin.x, &p.x),
        (&normal.y, &origin.y, &p.y),
        (&normal.z, &origin.z, &p.z),
    ]
    .into_iter()
    .all(|(n, a, b)| n.rational_ref().is_some_and(|n| n.is_zero()) || a.same_definition(b))
}

pub(crate) fn on_axis(origin: &Point3, axis: &Vector3, p: &Point3) -> bool {
    if same_point(origin, p) {
        return true;
    }
    let coordinates = [
        (&axis.x, &origin.x, &p.x),
        (&axis.y, &origin.y, &p.y),
        (&axis.z, &origin.z, &p.z),
    ];
    // A coordinate axis leaves exactly one coordinate free. General axes
    // require the shared origin or an explicit builder fact instead.
    coordinates
        .iter()
        .filter(|(d, _, _)| !d.rational_ref().is_some_and(|d| d.is_zero()))
        .count()
        == 1
        && coordinates
            .into_iter()
            .all(|(d, a, b)| !d.rational_ref().is_some_and(|d| d.is_zero()) || a.same_definition(b))
}

pub(crate) fn axis(s: &Surface) -> Option<(&Point3, &Vector3)> {
    match s {
        Surface::Cylinder(s) => Some((&s.origin, &s.axis)),
        Surface::Cone(s) => Some((&s.apex, &s.axis)),
        Surface::Torus(s) => Some((&s.center, &s.axis)),
        Surface::EllipticCylinder(s) => Some((&s.origin, &s.axis)),
        _ => None,
    }
}

pub(crate) fn center(s: &Surface) -> Option<&Point3> {
    match s {
        Surface::Sphere(s) => Some(&s.center),
        Surface::Torus(s) => Some(&s.center),
        _ => None,
    }
}

pub(crate) fn direction(s: &Surface) -> Option<&Vector3> {
    match s {
        Surface::Plane(s) => Some(&s.normal),
        _ => axis(s).map(|(_, direction)| direction),
    }
}

pub(crate) fn coincident(a: &Surface, b: &Surface) -> bool {
    match (a, b) {
        (Surface::Plane(a), Surface::Plane(b)) => {
            parallel(&a.normal, &b.normal) && on_plane(&a.origin, &a.normal, &b.origin)
        }
        (Surface::Cylinder(a), Surface::Cylinder(b)) => {
            on_axis(&a.origin, &a.axis, &b.origin)
                && parallel(&a.axis, &b.axis)
                && a.radius.same_definition(&b.radius)
        }
        (Surface::Cone(a), Surface::Cone(b)) => {
            same_point(&a.apex, &b.apex)
                && same_vector(&a.axis, &b.axis)
                && a.slope.same_definition(&b.slope)
        }
        (Surface::Sphere(a), Surface::Sphere(b)) => {
            same_point(&a.center, &b.center) && a.radius.same_definition(&b.radius)
        }
        (Surface::Torus(a), Surface::Torus(b)) => {
            same_point(&a.center, &b.center)
                && parallel(&a.axis, &b.axis)
                && a.major.same_definition(&b.major)
                && a.minor.same_definition(&b.minor)
        }
        _ => a.same_definition(b),
    }
}
