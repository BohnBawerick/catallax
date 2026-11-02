//! Cached hardware-interval exclusions; no decisions, refinement or geometry writes.

use catallax_numbers::Interval;
use catallax_shapes::{Curve, CurveId, GeomStore, Point3, Surface, Vector3};

type Box3 = [Interval; 3];
fn zero() -> Interval {
    Interval::point(0.0).expect("finite")
}
fn point(p: &Point3) -> Box3 {
    [p.x.approx(), p.y.approx(), p.z.approx()]
}
fn vector(v: &Vector3) -> Box3 {
    [v.x.approx(), v.y.approx(), v.z.approx()]
}
fn dot(a: Box3, b: Box3) -> Interval {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn squared(a: Box3) -> Interval {
    a[0].square() + a[1].square() + a[2].square()
}
fn cross(a: Box3, b: Box3) -> Box3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn around(center: &Point3, radius: Interval) -> Option<Box3> {
    let radius = Interval::new(-radius.hi(), radius.hi())?;
    Some(point(center).map(|coordinate| coordinate + radius))
}
fn surface_bounds(surface: &Surface) -> Option<Box3> {
    match surface {
        Surface::Sphere(s) if s.radius.approx().lo() > 0.0 => around(&s.center, s.radius.approx()),
        Surface::Cylinder(c)
            if c.radius_sq.approx().lo() > 0.0 && squared(vector(&c.axis)).lo() > 0.0 =>
        {
            let radius = c.radius_sq.approx().sqrt()?;
            let extent = around(&c.origin, radius)?;
            let axis = vector(&c.axis);
            Some(std::array::from_fn(|i| {
                if axis[i] == zero() {
                    extent[i]
                } else {
                    Interval::entire()
                }
            }))
        }
        Surface::Plane(p) => {
            let normal = vector(&p.normal);
            let axis = (0..3).find(|&i| {
                !normal[i].contains_zero() && (0..3).all(|j| i == j || normal[j] == zero())
            })?;
            Some(std::array::from_fn(|i| {
                if i == axis {
                    point(&p.origin)[i]
                } else {
                    Interval::entire()
                }
            }))
        }
        _ => None,
    }
}

/// Encloses the complete carrier, never just a sampled or retained arc.
pub(crate) fn enclosure(geom: &GeomStore, curve: CurveId) -> Option<Box3> {
    match geom.curve(curve) {
        Curve::Line(line) if squared(vector(&line.dir)).lo() > 0.0 => {
            let direction = vector(&line.dir);
            let origin = point(&line.origin);
            Some(std::array::from_fn(|i| {
                if direction[i] == zero() {
                    origin[i]
                } else {
                    Interval::entire()
                }
            }))
        }
        Curve::Circle(circle) if circle.radius_sq.approx().lo() > 0.0 => {
            around(&circle.center, circle.radius_sq.approx().sqrt()?)
        }
        Curve::Ellipse(ellipse)
            if ellipse.major_sq.approx().lo() > 0.0 && ellipse.minor_sq.approx().lo() > 0.0 =>
        {
            // Both semi-axes fit this sphere, independently of orientation.
            around(
                &ellipse.center,
                ellipse
                    .major_sq
                    .approx()
                    .hull(&ellipse.minor_sq.approx())
                    .sqrt()?,
            )
        }
        Curve::Intersection(curve) => {
            let a = surface_bounds(geom.surface(curve.surfaces[0]));
            let b = surface_bounds(geom.surface(curve.surfaces[1]));
            match (a, b) {
                (Some(a), Some(b)) => {
                    let mut result = [Interval::entire(); 3];
                    for i in 0..3 {
                        result[i] =
                            Interval::new(a[i].lo().max(b[i].lo()), a[i].hi().min(b[i].hi()))?;
                    }
                    Some(result)
                }
                (Some(bounds), None) | (None, Some(bounds)) => Some(bounds),
                (None, None) => None,
            }
        }
        _ => None,
    }
}

/// Strict exclusion only. Zero, tangency, overflow and unsupported supports fall through.
pub(crate) fn excludes(bounds: &Box3, surface: &Surface) -> bool {
    let value = match surface {
        Surface::Plane(p) => {
            let origin = point(&p.origin);
            dot(
                std::array::from_fn(|i| bounds[i] - origin[i]),
                vector(&p.normal),
            )
        }
        Surface::Sphere(s) if s.radius.approx().lo() > 0.0 => {
            let center = point(&s.center);
            squared(std::array::from_fn(|i| bounds[i] - center[i])) - s.radius.approx().square()
        }
        Surface::Cylinder(c) if c.radius_sq.approx().lo() > 0.0 => {
            let origin = point(&c.origin);
            let axis = vector(&c.axis);
            squared(cross(std::array::from_fn(|i| bounds[i] - origin[i]), axis))
                - c.radius_sq.approx() * squared(axis)
        }
        _ => return false,
    };
    !value.contains_zero()
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Decider, Rational, Real, with_operation_budget};
    use catallax_shapes::{
        Circle, Cylinder, Ellipse, IntersectionCurve, Plane, PointDef, RPoint3, RVector3, Sphere,
    };
    use catallax_topology::Model;

    fn p(x: i64, y: i64, z: i64) -> Point3 {
        RPoint3::new(
            Rational::from_i64(x),
            Rational::from_i64(y),
            Rational::from_i64(z),
        )
        .into()
    }
    fn v(x: i64, y: i64, z: i64) -> Vector3 {
        RVector3::new(
            Rational::from_i64(x),
            Rational::from_i64(y),
            Rational::from_i64(z),
        )
        .into()
    }
    fn plane(x: Real) -> Surface {
        Surface::Plane(Plane {
            origin: Point3::origin().plus(&v(1, 0, 0).scaled(&x)),
            normal: v(1, 0, 0),
            ref_dir: v(0, 1, 0),
        })
    }

    #[test]
    fn cached_exclusions_preserve_grazing_contacts_without_work() {
        let mut model = Model::new();
        let mut dec = Decider::default();
        let circle = model.geom_mut().add_curve(Curve::Circle(Circle {
            center: p(0, 0, 0),
            normal: v(0, 0, 1),
            radius_sq: Real::one(),
            ref_dir: v(1, 0, 0),
        }));
        let near = &Real::one() - &Real::from(Rational::parse("1e-100").unwrap());
        let tangent_cylinder = Surface::Cylinder(
            Cylinder::new(p(2, 0, 0), v(0, 0, 1), Real::one(), v(1, 0, 0), &mut dec).unwrap(),
        );
        let tangent_sphere = Surface::Sphere(
            Sphere::new(p(2, 0, 0), Real::one(), v(0, 0, 1), v(1, 0, 0), &mut dec).unwrap(),
        );
        let far_sphere = Surface::Sphere(
            Sphere::new(p(5, 0, 0), Real::one(), v(0, 0, 1), v(1, 0, 0), &mut dec).unwrap(),
        );
        let ellipse = model.geom_mut().add_curve(Curve::Ellipse(
            Ellipse::new(
                p(0, 0, 0),
                v(0, 0, 1),
                v(1, 0, 0),
                Real::from_i64(4),
                Real::one(),
                &mut dec,
            )
            .unwrap(),
        ));
        let first =
            Cylinder::new(p(8, 0, 4), v(0, 3, 0), Real::one(), v(1, 0, 0), &mut dec).unwrap();
        let second =
            Cylinder::new(p(8, 3, 4), v(-5, 0, 0), Real::one(), v(0, 1, 0), &mut dec).unwrap();
        let a = model.geom_mut().add_surface(Surface::Cylinder(first));
        let b = model.geom_mut().add_surface(Surface::Cylinder(second));
        let witness = RPoint3::new(
            Rational::parse("37/5").unwrap(),
            Rational::parse("12/5").unwrap(),
            Rational::parse("16/5").unwrap(),
        );
        let witness = model
            .geom_mut()
            .add_point(PointDef::Explicit(witness.into()));
        let branch = model
            .geom_mut()
            .add_curve(Curve::Intersection(IntersectionCurve {
                surfaces: [a, b],
                witness,
            }));
        let cases = [
            (circle, plane(Real::from_i64(2)), true),
            (circle, plane(Real::one()), false),
            (circle, plane(near), false),
            (circle, tangent_cylinder, false),
            (circle, tangent_sphere, false),
            (circle, far_sphere, true),
            (ellipse, plane(Real::from_i64(3)), true),
            (ellipse, plane(Real::from_i64(2)), false),
            (branch, plane(Real::zero()), true),
            (branch, plane(Real::from_i64(7)), false),
        ];
        let ((), work) = with_operation_budget(0, || {
            for (curve, surface, expected) in &cases {
                assert_eq!(
                    excludes(&enclosure(model.geom(), *curve).unwrap(), surface),
                    *expected
                );
            }
        });
        assert!(!work.exhausted, "{work:?}");
        assert_eq!(work.used, 0);
    }
}
