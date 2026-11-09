//! Exact sector rejection before intersecting a cliff's boundary carriers.

use catallax_numbers::{Decider, Refusal, Sign};
use catallax_shapes::{Curve, Plane, Point3, Vector3};
use catallax_topology::{EdgeId, Model};

use crate::{fillet::FilletDefinition, profile::sign};

pub(super) struct Sector([Plane; 2]);

fn rational_point(point: &Point3) -> bool {
    [&point.x, &point.y, &point.z]
        .iter()
        .all(|coordinate| coordinate.rational_ref().is_some())
}
fn rational_vector(vector: &Vector3) -> bool {
    [&vector.x, &vector.y, &vector.z]
        .iter()
        .all(|coordinate| coordinate.rational_ref().is_some())
}

impl Sector {
    pub(super) fn new(definition: &FilletDefinition) -> Result<Self, Refusal> {
        definition.sector_halfspaces().map(Self)
    }

    fn rational(&self) -> bool {
        self.0
            .iter()
            .all(|plane| rational_point(&plane.origin) && rational_vector(&plane.normal))
    }

    fn rejects_segment(&self, a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        if !self.rational() || !rational_point(a) || !rational_point(b) {
            return Ok(false);
        }
        for plane in &self.0 {
            if sign(&a.minus(&plane.origin).dot(&plane.normal), dec)? == Sign::Negative
                && sign(&b.minus(&plane.origin).dot(&plane.normal), dec)? == Sign::Negative
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn rejects_conic(&self, curve: &Curve, dec: &mut Decider) -> Result<bool, Refusal> {
        // Optional culling must not start a new algebraic decision before the
        // complete intersection. Nonrational carriers keep that original path.
        let rational = match curve {
            Curve::Circle(circle) => {
                rational_point(&circle.center)
                    && rational_vector(&circle.normal)
                    && circle.radius_sq.rational_ref().is_some()
            }
            Curve::Ellipse(ellipse) => {
                rational_point(&ellipse.center)
                    && rational_vector(&ellipse.normal)
                    && rational_vector(&ellipse.major_dir)
                    && ellipse.major_sq.rational_ref().is_some()
                    && ellipse.minor_sq.rational_ref().is_some()
            }
            _ => false,
        };
        if !self.rational() || !rational {
            return Ok(false);
        }
        let center = match curve {
            Curve::Circle(circle) => &circle.center,
            Curve::Ellipse(ellipse) => &ellipse.center,
            _ => return Ok(false),
        };
        for plane in &self.0 {
            let offset = center.minus(&plane.origin).dot(&plane.normal);
            if sign(&offset, dec)? != Sign::Negative {
                continue;
            }
            // The support of the whole conic in this normal direction. A
            // trimmed arc is a subset, so a strict separation rejects it too.
            let extent_squared = match curve {
                Curve::Circle(circle) => {
                    let height = plane.normal.dot(&circle.normal);
                    &circle.radius_sq
                        * &(&plane.normal.norm_squared()
                            - &(&height.square() / &circle.normal.norm_squared()))
                }
                Curve::Ellipse(ellipse) => {
                    let minor = ellipse.normal.cross(&ellipse.major_dir);
                    let along = plane.normal.dot(&ellipse.major_dir);
                    let across = plane.normal.dot(&minor);
                    &(&ellipse.major_sq * &(&along.square() / &ellipse.major_dir.norm_squared()))
                        + &(&ellipse.minor_sq * &(&across.square() / &minor.norm_squared()))
                }
                _ => unreachable!("the center belongs to a conic"),
            };
            if sign(&(&offset.square() - &extent_squared), dec)? == Sign::Positive {
                return Ok(true);
            }
        }
        Ok(false)
    }

    pub(super) fn rejects(
        &self,
        model: &Model,
        edge: EdgeId,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let edge = model.edge(edge);
        let curve = model.geom().curve(edge.curve);
        if let (Curve::Line(_), Some(a), Some(b)) = (curve, edge.start, edge.end) {
            let a = model.geom().position(model.vertex(a).point);
            let b = model.geom().position(model.vertex(b).point);
            return self.rejects_segment(&a, &b, dec);
        }
        self.rejects_conic(curve, dec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::Real;
    use catallax_shapes::{Circle, Ellipse};

    fn point(x: i64, y: i64, z: i64) -> Point3 {
        Point3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
    }
    fn vector(x: i64, y: i64, z: i64) -> Vector3 {
        Vector3::new(Real::from_i64(x), Real::from_i64(y), Real::from_i64(z))
    }
    fn sector() -> Sector {
        Sector([vector(1, 0, 0), vector(0, 1, 0)].map(|normal| Plane {
            origin: Point3::origin(),
            normal,
            ref_dir: vector(0, 0, 1),
        }))
    }

    #[test]
    fn halfspaces_match_the_original_sector_predicate() {
        use catallax_recipe::{Origin, Role, StepId};
        use catallax_shapes::{Cylinder, Line, Surface};
        let axis = vector(0, 0, 3);
        let definition = FilletDefinition {
            origin: Origin {
                step: StepId::new("sector").unwrap(),
                role: Role::Start,
            },
            supports: [vector(-2, 0, 0), vector(0, -3, 0)].map(|normal| {
                Surface::Plane(Plane {
                    origin: Point3::origin(),
                    normal,
                    ref_dir: axis.clone(),
                })
            }),
            same_sense: [true, true],
            edge: Curve::Line(Line {
                origin: Point3::origin(),
                dir: axis.clone(),
            }),
            radius: Real::from_i64(2),
            requested_radius: None,
            convex: true,
            surface: Surface::Cylinder(Cylinder {
                origin: point(2, 2, 0),
                axis: axis.clone(),
                radius_sq: Real::from_i64(4),
                ref_dir: vector(-1, 0, 0),
            }),
            contacts: [point(0, 2, 0), point(2, 0, 0)].map(|origin| {
                Curve::Line(Line {
                    origin,
                    dir: axis.clone(),
                })
            }),
        };
        let mut dec = Decider::default();
        for convex in [true, false] {
            let definition = FilletDefinition {
                convex,
                ..definition.clone()
            };
            let sector = Sector::new(&definition).unwrap();
            for x in [-1, 0, 1, 2, 3] {
                for y in [-1, 0, 1, 2, 3] {
                    for z in [-7, 0, 11] {
                        let point = point(x, y, z);
                        let expected = definition.in_sector(&point, &mut dec).unwrap();
                        let mut actual = true;
                        for plane in &sector.0 {
                            actual &=
                                sign(&point.minus(&plane.origin).dot(&plane.normal), &mut dec)
                                    .unwrap()
                                    != Sign::Negative;
                        }
                        assert_eq!(actual, expected, "{convex} {x} {y} {z}");
                    }
                }
            }
        }
    }

    #[test]
    fn segments_need_strict_separation_on_one_common_halfspace() {
        let sector = sector();
        let mut dec = Decider::default();
        assert!(
            sector
                .rejects_segment(&point(-2, 1, 0), &point(-1, 2, 3), &mut dec)
                .unwrap()
        );
        assert!(
            !sector
                .rejects_segment(&point(0, 1, 0), &point(-1, 2, 3), &mut dec)
                .unwrap()
        );
        assert!(
            !sector
                .rejects_segment(&point(-1, 2, 0), &point(2, -1, 3), &mut dec)
                .unwrap()
        );
    }

    #[test]
    fn complete_conic_support_bounds_keep_tangencies_and_nonunit_frames() {
        let sector = sector();
        let mut dec = Decider::default();
        for (normal, radius_sq) in [(vector(0, 0, 3), 1), (vector(1, 0, 1), 2)] {
            let circle = |x| {
                Curve::Circle(Circle {
                    center: point(x, 3, 0),
                    normal: normal.clone(),
                    radius_sq: Real::from_i64(radius_sq),
                    ref_dir: vector(0, 1, 0),
                })
            };
            assert!(sector.rejects_conic(&circle(-3), &mut dec).unwrap());
            assert!(!sector.rejects_conic(&circle(-1), &mut dec).unwrap());
            assert!(!sector.rejects_conic(&circle(3), &mut dec).unwrap());
        }
        let ellipse = |x| {
            Curve::Ellipse(Ellipse {
                center: point(x, 3, 0),
                normal: vector(0, 0, 3),
                major_dir: vector(0, 2, 0),
                major_sq: Real::from_i64(9),
                minor_sq: Real::one(),
            })
        };
        assert!(sector.rejects_conic(&ellipse(-3), &mut dec).unwrap());
        assert!(!sector.rejects_conic(&ellipse(-1), &mut dec).unwrap());
    }
}
