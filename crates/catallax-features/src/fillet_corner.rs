//! Exact spherical corners where three convex planar strips meet.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Circle, Curve, Point3, Sphere, Surface, Vector3};

use crate::{fillet::FilletDefinition, profile::sign};

/// The exact generating definition of one trihedral rolling-ball corner.
///
/// Three independent planes meet at a convex vertex and all three incident
/// straight edges carry the same radius. The ball touching all three planes
/// sweeps the sphere sector in the positive cone of their outward normals.
#[derive(Clone, Debug)]
pub struct CornerDefinition {
    /// Persistent `corner.vertex` origin, independent of the radius.
    pub origin: Origin,
    /// The three planar supports, each shared by two of the incident strips.
    pub supports: [Surface; 3],
    /// Unit outward material normals of the three supports.
    pub normals: [Vector3; 3],
    /// The input vertex where the three supports meet.
    pub vertex: Point3,
    /// Positive rolling radius in millimetres.
    pub radius: Real,
    /// The exact sphere about the three inward offsets.
    pub surface: Surface,
    /// The contact point on each support, in support order.
    pub feet: [Point3; 3],
    /// For each strip, the indices of its two supports in strip contact order.
    pub strips: [[usize; 2]; 3],
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-corner", detail)
}
fn unit(v: &Vector3) -> Vector3 {
    let squared = v.norm_squared();
    let norm = squared
        .as_rational()
        .and_then(|r| crate::extrude::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| squared.sqrt());
    v.scaled(&(&Real::one() / &norm))
}

impl CornerDefinition {
    /// Proves the three strips form an equal-radius convex corner of independent
    /// planes whose ball centre lies behind the vertex along every incident
    /// edge, and reduces it to its exact sphere. Every other junction is refused;
    /// this does not trim or publish a solid.
    pub fn new(
        origin: Origin,
        strips: [&FilletDefinition; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut supports: Vec<(Surface, bool)> = Vec::with_capacity(3);
        let mut indices = [[0; 2]; 3];
        for (strip, index) in strips.iter().zip(&mut indices) {
            if !strip.convex {
                return Err(bad(
                    "a concave strip at a three-strip vertex needs a general corner patch",
                ));
            }
            if sign(&(&strip.radius - &strips[0].radius), dec)? != Sign::Zero {
                return Err(bad("unequal radii at a vertex need a general corner patch"));
            }
            for (k, (support, same)) in strip.supports.iter().zip(strip.same_sense).enumerate() {
                if !matches!(support, Surface::Plane(_)) {
                    return Err(bad("an exact spherical corner needs three planar supports"));
                }
                let found = supports
                    .iter()
                    .position(|(s, _)| s.same_definition(support));
                index[k] = match found {
                    Some(i) if supports[i].1 != same => {
                        return Err(bad("a shared corner support changes orientation"));
                    }
                    Some(i) => i,
                    None => {
                        supports.push((support.clone(), same));
                        supports.len() - 1
                    }
                };
            }
        }
        if supports.len() != 3
            || indices.iter().any(|[a, b]| a == b)
            || (0..3).any(|i| indices.iter().filter(|s| s.contains(&i)).count() != 2)
        {
            return Err(bad(
                "three strips at a vertex must pair up three distinct supports",
            ));
        }
        let normals: Vec<Vector3> = supports
            .iter()
            .map(|(s, same)| {
                let Surface::Plane(p) = s else { unreachable!() };
                unit(&p.normal).scaled(&Real::from_i64(if *same { 1 } else { -1 }))
            })
            .collect();
        if sign(&normals[0].dot(&normals[1].cross(&normals[2])), dec)? == Sign::Zero {
            return Err(bad("an exact spherical corner needs independent supports"));
        }
        let radius = strips[0].radius.clone();
        // The ball centre lies at depth `radius` behind every support. On the
        // first strip's axis it is fixed by the third support; the other axes
        // must pass through the same point.
        let Surface::Cylinder(first) = &strips[0].surface else {
            return Err(bad("a planar strip must reduce to a cylinder"));
        };
        let third = (0..3).find(|i| !indices[0].contains(i)).unwrap();
        let Surface::Plane(p) = &supports[third].0 else {
            unreachable!()
        };
        let n = &normals[third];
        let divisor = first.axis.dot(n);
        if sign(&divisor, dec)? == Sign::Zero {
            return Err(bad("the third support is parallel to a strip axis"));
        }
        let depth = -(&radius + &first.origin.minus(&p.origin).dot(n));
        let centre = first.origin.plus(&first.axis.scaled(&(depth / divisor)));
        for strip in &strips[1..] {
            let Surface::Cylinder(c) = &strip.surface else {
                return Err(bad("a planar strip must reduce to a cylinder"));
            };
            if sign(&centre.minus(&c.origin).cross(&c.axis).norm_squared(), dec)? != Sign::Zero {
                return Err(bad("the three strip axes do not meet at one ball centre"));
            }
        }
        for (i, (support, _)) in supports.iter().enumerate() {
            let Surface::Plane(p) = support else {
                unreachable!()
            };
            let depth = &centre.minus(&p.origin).dot(&normals[i]) + &radius;
            if sign(&depth, dec)? != Sign::Zero {
                return Err(bad("the ball centre is not at radius behind every support"));
            }
        }
        // The vertex is on the first edge, at depth zero behind the third
        // support. A centre past it along any edge would extend that strip
        // beyond the vertex, which neither evaluator's strip region covers.
        let Curve::Line(edge) = &strips[0].edge else {
            return Err(bad("a planar strip must reduce to a cylinder"));
        };
        let reach = -&edge.origin.minus(&p.origin).dot(n) / edge.dir.dot(n);
        let vertex = edge.origin.plus(&edge.dir.scaled(&reach));
        for (strip, [a, b]) in strips.iter().zip(&indices) {
            let Curve::Line(edge) = &strip.edge else {
                return Err(bad("a planar strip must reduce to a cylinder"));
            };
            let far = (0..3).find(|i| i != a && i != b).unwrap();
            let along = centre.minus(&vertex).dot(&edge.dir) * -edge.dir.dot(&normals[far]);
            if sign(&along, dec)? != Sign::Positive {
                return Err(bad(
                    "a ball centre beyond the vertex needs a general corner patch",
                ));
            }
        }
        let feet: [Point3; 3] = std::array::from_fn(|i| centre.plus(&normals[i].scaled(&radius)));
        // Perpendicular contact normals already form the sphere's frame.
        let tilt = normals[0].dot(&normals[2]);
        let reference = if sign(&tilt, dec)? == Sign::Zero {
            normals[0].clone()
        } else {
            normals[0].minus(&normals[2].scaled(&tilt))
        };
        let surface = Surface::Sphere(Sphere::new(
            centre,
            radius.clone(),
            normals[2].clone(),
            reference,
            dec,
        )?);
        Ok(Self {
            origin,
            supports: [0, 1, 2].map(|i| supports[i].0.clone()),
            normals: [0, 1, 2].map(|i| normals[i].clone()),
            vertex,
            radius,
            surface,
            feet,
            strips: indices,
        })
    }

    /// The ball centre.
    pub fn centre(&self) -> &Point3 {
        let Surface::Sphere(s) = &self.surface else {
            unreachable!("corner carrier is a sphere")
        };
        &s.center
    }

    /// The great circle shared with one strip, whose shorter arc runs from the
    /// strip's first contact foot to its second along the sector boundary.
    pub fn boundary(&self, strip: usize) -> Curve {
        let [a, b] = self.strips[strip];
        Curve::Circle(Circle {
            center: self.centre().clone(),
            normal: self.normals[a].cross(&self.normals[b]),
            radius_sq: self.radius.square(),
            ref_dir: self.normals[a].clone(),
        })
    }

    /// For each normal, the vector whose dot product with a displacement gives
    /// that displacement's coefficient on the normal. A displacement is in the
    /// positive cone when every coefficient is nonnegative.
    pub fn dual_basis(&self) -> [Vector3; 3] {
        let [a, b, c] = &self.normals;
        let volume = Real::one() / a.dot(&b.cross(c));
        [b.cross(c), c.cross(a), a.cross(b)].map(|v| v.scaled(&volume))
    }

    /// Whether a point's direction from the centre lies in the closed sector
    /// spanned by the three outward normals.
    pub fn in_sector(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let d = p.minus(self.centre());
        for dual in self.dual_basis() {
            if sign(&d.dot(&dual), dec)? == Sign::Negative {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
