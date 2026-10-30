//! Positive construction identities only. A failed match says nothing about geometry.

use crate::{Point3, Surface, Vector3};

pub(crate) fn same_point(a: &Point3, b: &Point3) -> bool {
    a.x.same_definition(&b.x) && a.y.same_definition(&b.y) && a.z.same_definition(&b.z)
}

fn same_vector(a: &Vector3, b: &Vector3) -> bool {
    a.x.same_definition(&b.x) && a.y.same_definition(&b.y) && a.z.same_definition(&b.z)
}

// Opposite direction definitions prove parallelism without allocating negations.
pub(crate) fn parallel(a: &Vector3, b: &Vector3) -> bool {
    same_vector(a, b)
        || (a.x.same_negated_definition(&b.x)
            && a.y.same_negated_definition(&b.y)
            && a.z.same_negated_definition(&b.z))
}

pub(crate) fn center(s: &Surface) -> Option<&Point3> {
    match s {
        Surface::Sphere(s) => Some(&s.center),
        Surface::Torus(s) => Some(&s.center),
        Surface::TorusSheet(s) => Some(&s.center),
        _ => None,
    }
}
