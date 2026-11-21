//! Re-tightening: an edge between two exact surfaces is exact again.
//!
//! Checkpoint 2 decisions, section 4.5: "Re-tighten tolerant edges between
//! exact surfaces after each operation, as Parasolid's `PK_EDGE_optimise` and
//! ACIS tweak do." Admission ends every numerical operation, and it gives
//! every edge the tolerance it measured between the edge's 3D curve and its
//! copy on each face. Where the 3D curve lies exactly on both surfaces the
//! edge needs none: it is an edge at default precision in Jackson's sense,
//! one 3D curve and no tube ("Boundary representation modelling with local
//! tolerances", 1995). Such an edge is stored exact, with the two `Incident`
//! facts an exact builder would have recorded.
//!
//! Parasolid and ACIS get there by intersecting the two surfaces again and
//! replacing the edge's curve. Nothing is replaced here: the curve is kept
//! and tested. A line or a circle lies on a plane, a cylinder, a cone, a
//! sphere or a torus exactly when a few polynomial identities hold between
//! their defining numbers, and where those numbers are rationals the
//! identities are decided by rational arithmetic. They are the identities of
//! the certified STEP import (`check_on_surface` in
//! `crates/catallax-import/src/build.rs`), which follow from the surfaces'
//! implicit equations. An edge that fails one, or whose numbers are not all
//! rational, keeps its measured tolerance.
//!
//! The vertices of such an edge keep theirs: a vertex of a numerical body is
//! a rounded point, and the rule vertex >= edge allows it.

use super::*;
use catallax_numbers::{Fact, FactKind, FactOrigin};
use catallax_shapes::{Circle, GeomKey, Line};

type V = [Rational; 3];

fn rational(value: &catallax_numbers::Real) -> Option<Rational> {
    value.rational_ref().cloned()
}

fn point(p: &Point3) -> Option<V> {
    Some([rational(&p.x)?, rational(&p.y)?, rational(&p.z)?])
}

fn vector(v: &catallax_shapes::Vector3) -> Option<V> {
    Some([rational(&v.x)?, rational(&v.y)?, rational(&v.z)?])
}

fn sub(a: &V, b: &V) -> V {
    [&a[0] - &b[0], &a[1] - &b[1], &a[2] - &b[2]]
}

fn dot(a: &V, b: &V) -> Rational {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}

fn cross(a: &V, b: &V) -> V {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

fn zero(v: &V) -> bool {
    v.iter().all(Rational::is_zero)
}

fn parallel(a: &V, b: &V) -> bool {
    zero(&cross(a, b))
}

/// Whether `curve` lies on `surface`, decided by rational identities.
/// `false` also where the pair has no identity here or a defining number is
/// not a rational: the edge then keeps its tolerance.
pub(super) fn lies_on(curve: &Curve, surface: &Surface) -> bool {
    match curve {
        Curve::Line(line) => line_on(line, surface),
        Curve::Circle(circle) => circle_on(circle, surface),
        _ => None,
    }
    .unwrap_or(false)
}

fn line_on(line: &Line, surface: &Surface) -> Option<bool> {
    let (p, d) = (point(&line.origin)?, vector(&line.dir)?);
    Some(match surface {
        // A point of the line on the plane, and its direction in the plane.
        Surface::Plane(plane) => {
            let (o, n) = (point(&plane.origin)?, vector(&plane.normal)?);
            dot(&sub(&p, &o), &n).is_zero() && dot(&d, &n).is_zero()
        }
        // A ruling: along the axis, at the cylinder's radius from it.
        Surface::Cylinder(cylinder) => {
            let (o, a) = (point(&cylinder.origin)?, vector(&cylinder.axis)?);
            let arm = cross(&sub(&p, &o), &a);
            parallel(&d, &a) && dot(&arm, &arm) == &rational(&cylinder.radius_sq)? * &dot(&a, &a)
        }
        _ => return None,
    })
}

fn circle_on(circle: &Circle, surface: &Surface) -> Option<bool> {
    let (c, m) = (point(&circle.center)?, vector(&circle.normal)?);
    let rho = rational(&circle.radius_sq)?;
    Some(match surface {
        // Centre on the plane, and the circle's plane parallel to it.
        Surface::Plane(plane) => {
            let (o, n) = (point(&plane.origin)?, vector(&plane.normal)?);
            dot(&sub(&c, &o), &n).is_zero() && parallel(&m, &n)
        }
        // A cross-section: centred on the axis, across it, at its radius.
        Surface::Cylinder(cylinder) => {
            let (o, a) = (point(&cylinder.origin)?, vector(&cylinder.axis)?);
            parallel(&m, &a) && parallel(&sub(&c, &o), &a) && rho == rational(&cylinder.radius_sq)?
        }
        // A cross-section of the nappe: its radius is the slope times the
        // distance from the apex along the unit axis.
        Surface::Cone(cone) => {
            let (apex, a) = (point(&cone.apex)?, vector(&cone.axis)?);
            let slope = rational(&cone.slope)?;
            let from_apex = sub(&c, &apex);
            let along = dot(&from_apex, &a);
            parallel(&m, &a)
                && parallel(&from_apex, &a)
                && along.sign() == catallax_numbers::Sign::Positive
                && &rho * &dot(&a, &a) == &(&slope * &slope) * &(&along * &along)
        }
        // The circle's axis passes through the centre of the sphere, and its
        // radius closes the right triangle.
        Surface::Sphere(sphere) => {
            let s = point(&sphere.center)?;
            let radius = rational(&sphere.radius)?;
            let offset = sub(&c, &s);
            parallel(&offset, &m) && &dot(&offset, &offset) + &rho == &radius * &radius
        }
        Surface::Torus(torus) => {
            let (t, a) = (point(&torus.center)?, vector(&torus.axis)?);
            let (major, minor) = (rational(&torus.major)?, rational(&torus.minor)?);
            let offset = sub(&c, &t);
            let along = dot(&offset, &a);
            if parallel(&m, &a) {
                // A parallel: centred on the axis at height h, with radius r
                // where (r - major)^2 + h^2 = minor^2. With rho = r^2 that is
                // k = 2 major r for k = rho + major^2 + h^2 - minor^2, so
                // k >= 0 and k^2 = 4 major^2 rho.
                let height_sq = &(&along * &along) / &dot(&a, &a);
                let k = &(&(&rho + &(&major * &major)) + &height_sq) - &(&minor * &minor);
                parallel(&offset, &a)
                    && k.sign() != catallax_numbers::Sign::Negative
                    && &k * &k == &(&Rational::from_i64(4) * &(&major * &major)) * &rho
            } else {
                // A meridian: centred on the major circle, in the plane
                // through the axis, at the minor radius.
                along.is_zero()
                    && dot(&offset, &offset) == &major * &major
                    && dot(&m, &a).is_zero()
                    && dot(&m, &offset).is_zero()
                    && rho == &minor * &minor
            }
        }
        _ => return None,
    })
}

/// The edges of `definition` whose 3D curve lies exactly on the surface of
/// every face that uses it, by index. An edge no face uses is not one.
pub(super) fn exact_edges(definition: &NumericalDefinition) -> Vec<bool> {
    let mut used = vec![false; definition.edges.len()];
    let mut exact = vec![true; definition.edges.len()];
    for face in &definition.faces {
        for trim in face.loops.iter().flatten() {
            used[trim.edge] = true;
            if exact[trim.edge] {
                exact[trim.edge] = lies_on(&definition.edges[trim.edge].curve, &face.surface);
            }
        }
    }
    exact
        .into_iter()
        .zip(used)
        .map(|(exact, used)| exact && used)
        .collect()
}

/// Records that `curve` lies on `surface`: the fact an exact builder records
/// for an edge of a face, here with the identities above as its proof.
pub(super) fn record(model: &mut Model, curve: CurveId, surface: SurfaceId) {
    let (a, b) = (GeomKey::Curve(curve), GeomKey::Surface(surface));
    if !model.geom().facts().holds(FactKind::Incident, a, b) {
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a,
            b,
            at: None,
            origin: FactOrigin::Construction("numerical.retighten"),
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Real;
    use catallax_shapes::{Cone, Cylinder, Plane, Sphere, Torus, Vector3};

    fn q(text: &str) -> Real {
        Real::from_rational(Rational::parse(text).unwrap())
    }
    fn p(x: &str, y: &str, z: &str) -> Point3 {
        Point3::new(q(x), q(y), q(z))
    }
    fn v(x: &str, y: &str, z: &str) -> Vector3 {
        Vector3 {
            x: q(x),
            y: q(y),
            z: q(z),
        }
    }
    fn line(origin: Point3, dir: Vector3) -> Curve {
        Curve::Line(Line { origin, dir })
    }
    fn circle(center: Point3, normal: Vector3, radius_sq: &str) -> Curve {
        Curve::Circle(Circle {
            center,
            normal,
            radius_sq: q(radius_sq),
            ref_dir: v("1", "0", "0"),
        })
    }
    fn z() -> Vector3 {
        v("0", "0", "1")
    }
    fn x() -> Vector3 {
        v("1", "0", "0")
    }

    #[test]
    fn a_line_or_a_circle_on_a_plane_or_a_cylinder() {
        let plane = Surface::Plane(Plane {
            origin: p("0", "0", "5"),
            normal: z(),
            ref_dir: x(),
        });
        assert!(lies_on(&line(p("1", "2", "5"), v("3", "4", "0")), &plane));
        assert!(!lies_on(
            &line(p("1", "2", "5"), v("3", "4", "1/1000000")),
            &plane
        ));
        assert!(!lies_on(
            &line(p("1", "2", "5000001/1000000"), v("3", "4", "0")),
            &plane
        ));
        assert!(lies_on(
            &circle(p("7", "8", "5"), v("0", "0", "-2"), "9"),
            &plane
        ));
        assert!(!lies_on(
            &circle(p("7", "8", "5"), v("0", "1/1000", "1"), "9"),
            &plane
        ));
        // An axis that is not a unit vector, as a pose leaves it.
        let cylinder = Surface::Cylinder(Cylinder {
            origin: p("1", "1", "0"),
            axis: v("0", "0", "3"),
            radius_sq: q("25"),
            ref_dir: x(),
        });
        assert!(lies_on(&line(p("4", "5", "-7"), z()), &cylinder));
        assert!(!lies_on(
            &line(p("4", "5", "-7"), v("0", "1/1000", "1")),
            &cylinder
        ));
        assert!(!lies_on(
            &line(p("4", "5000001/1000000", "-7"), z()),
            &cylinder
        ));
        assert!(lies_on(&circle(p("1", "1", "12"), z(), "25"), &cylinder));
        assert!(!lies_on(
            &circle(p("1", "1", "12"), z(), "25000001/1000000"),
            &cylinder
        ));
        assert!(!lies_on(
            &circle(p("1", "1000001/1000000", "12"), z(), "25"),
            &cylinder
        ));
    }

    #[test]
    fn a_circle_on_a_cone_a_sphere_or_a_torus() {
        // Radius 3/2 for each unit along an axis of length 2.
        let cone = Surface::Cone(Cone {
            apex: p("0", "0", "1"),
            axis: v("0", "0", "2"),
            slope: q("3/2"),
            ref_dir: x(),
        });
        assert!(lies_on(&circle(p("0", "0", "3"), z(), "9"), &cone));
        assert!(!lies_on(
            &circle(p("0", "0", "3"), z(), "9000001/1000000"),
            &cone
        ));
        // The other nappe is not the surface.
        assert!(!lies_on(&circle(p("0", "0", "-1"), z(), "9"), &cone));
        let sphere = Surface::Sphere(Sphere {
            center: p("1", "2", "3"),
            radius: q("5"),
            axis: z(),
            ref_dir: x(),
        });
        assert!(lies_on(&circle(p("1", "2", "6"), z(), "16"), &sphere));
        assert!(lies_on(
            &circle(p("1", "2", "3"), v("1", "1", "0"), "25"),
            &sphere
        ));
        assert!(!lies_on(&circle(p("1", "2", "6"), z(), "17"), &sphere));
        assert!(!lies_on(
            &circle(p("1", "2001/1000", "6"), z(), "16"),
            &sphere
        ));
        let torus = Surface::Torus(Torus {
            center: p("0", "0", "0"),
            axis: v("0", "0", "2"),
            major: q("5"),
            minor: q("2"),
            ref_dir: x(),
        });
        // The outer and the inner equator, and the top parallel.
        assert!(lies_on(&circle(p("0", "0", "0"), z(), "49"), &torus));
        assert!(lies_on(&circle(p("0", "0", "0"), z(), "9"), &torus));
        assert!(lies_on(&circle(p("0", "0", "2"), z(), "25"), &torus));
        assert!(!lies_on(&circle(p("0", "0", "2"), z(), "26"), &torus));
        assert!(!lies_on(&circle(p("0", "0", "1"), z(), "49"), &torus));
        // A meridian, and one off the major circle.
        assert!(lies_on(
            &circle(p("3", "4", "0"), v("-4", "3", "0"), "4"),
            &torus
        ));
        assert!(!lies_on(
            &circle(p("3", "4", "0"), v("-4", "3", "0"), "5"),
            &torus
        ));
        assert!(!lies_on(
            &circle(p("3", "4", "1/10"), v("-4", "3", "0"), "4"),
            &torus
        ));
    }

    #[test]
    fn a_number_that_is_not_rational_or_a_pair_with_no_identity_proves_nothing() {
        let plane = Surface::Plane(Plane {
            origin: p("0", "0", "0"),
            normal: z(),
            ref_dir: x(),
        });
        // The line is on the plane, but its origin is a root.
        let root = Point3::new(q("2").sqrt(), q("0"), q("0"));
        assert!(!lies_on(&line(root, x()), &plane));
        // A ruling of a cone is not among the identities.
        let cone = Surface::Cone(Cone {
            apex: p("0", "0", "0"),
            axis: z(),
            slope: q("1"),
            ref_dir: x(),
        });
        assert!(!lies_on(&line(p("0", "0", "0"), v("1", "0", "1")), &cone));
    }
}
