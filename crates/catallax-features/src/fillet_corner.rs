//! Exact rolling-ball sectors and offset normal fans of planar vertices.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Circle, Curve, Point3, Sphere, Surface, Vector3};

use crate::{
    fillet::{FilletDefinition, UNSUPPORTED_CORNER},
    profile::sign,
};

/// The corner patch of a ball touching three supports at once: the part of
/// its sphere whose directions from the centre lie in the positive cone of
/// the three unit directions to its contacts (recipe section 4.5.1). It
/// closes the junction of three strips that share that ball whatever their
/// supports, given the ball and its contact directions. Each pair of
/// contacts bounds it by the shorter great-circle arc between them. A convex
/// sector removes the material it bounds and a concave one adds it.
#[derive(Clone, Debug)]
pub struct SphereSector {
    /// The ball centre.
    pub centre: Point3,
    /// Positive ball radius in millimetres.
    pub radius: Real,
    /// Unit directions from the centre to the three contacts.
    pub directions: [Vector3; 3],
    /// Material is removed by a convex sector and added by a concave one.
    pub convex: bool,
    /// The exact sphere.
    pub surface: Surface,
    /// The three contacts, in direction order.
    pub feet: [Point3; 3],
}

impl SphereSector {
    /// The sector of the ball about `centre` spanned by three independent
    /// contact directions, which must be unit: they are not normalised here.
    /// This does not trim or publish a face.
    pub fn new(
        centre: Point3,
        radius: Real,
        directions: [Vector3; 3],
        convex: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if sign(&radius, dec)? != Sign::Positive {
            return Err(bad("a sphere sector needs a positive radius"));
        }
        if sign(
            &directions[0].dot(&directions[1].cross(&directions[2])),
            dec,
        )? == Sign::Zero
        {
            return Err(bad("an exact spherical corner needs independent supports"));
        }
        let feet = std::array::from_fn(|i| centre.plus(&directions[i].scaled(&radius)));
        // Perpendicular contact directions already form the sphere's frame.
        let tilt = directions[0].dot(&directions[2]);
        let reference = if sign(&tilt, dec)? == Sign::Zero {
            directions[0].clone()
        } else {
            directions[0].minus(&directions[2].scaled(&tilt))
        };
        let surface = Surface::Sphere(Sphere::new(
            centre.clone(),
            radius.clone(),
            directions[2].clone(),
            reference,
            dec,
        )?);
        Ok(Self {
            centre,
            radius,
            directions,
            convex,
            surface,
            feet,
        })
    }

    /// The great circle through contacts `a` and `b`, oriented so that its
    /// shorter arc runs from the first to the second along the sector boundary.
    pub fn boundary(&self, a: usize, b: usize) -> Curve {
        Curve::Circle(Circle {
            center: self.centre.clone(),
            normal: self.directions[a].cross(&self.directions[b]),
            radius_sq: self.radius.square(),
            ref_dir: self.directions[a].clone(),
        })
    }

    /// For each direction, the vector whose dot product with a displacement
    /// gives that displacement's coefficient on the direction. A displacement
    /// is in the positive cone when every coefficient is nonnegative.
    pub fn dual_basis(&self) -> [Vector3; 3] {
        let [a, b, c] = &self.directions;
        let volume = Real::one() / a.dot(&b.cross(c));
        [b.cross(c), c.cross(a), a.cross(b)].map(|v| v.scaled(&volume))
    }

    /// Whether a point's direction from the centre lies in the closed sector.
    pub fn in_sector(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let d = p.minus(&self.centre);
        for dual in self.dual_basis() {
            if sign(&d.dot(&dual), dec)? == Sign::Negative {
                return Ok(false);
            }
        }
        Ok(true)
    }
}

/// The exact generating definition of one trihedral rolling-ball corner.
///
/// Three independent planes meet at a vertex, all three incident straight
/// edges are convex or all three concave, and they carry the same radius. The
/// ball touching all three planes sweeps its [`SphereSector`]. A convex corner
/// removes that sector's material and a concave one adds it: the concave
/// corner is the convex corner of the complement.
#[derive(Clone, Debug)]
pub struct CornerDefinition {
    /// Persistent `corner.vertex` origin, independent of the radius.
    pub origin: Origin,
    /// The three planar supports, each shared by two of the incident strips.
    pub supports: [Surface; 3],
    /// The input vertex where the three supports meet.
    pub vertex: Point3,
    /// The sector, its contact directions in support order: the outward
    /// material normals at a convex corner and the inward ones at a concave
    /// corner.
    pub sector: SphereSector,
    /// For each strip, the indices of its two supports in strip contact order.
    pub strips: [[usize; 2]; 3],
}

fn bad(detail: &str) -> Refusal {
    Refusal::new(UNSUPPORTED_CORNER, detail)
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
    /// Proves the three strips form an equal-radius corner of independent
    /// planes, all convex or all concave, whose ball centre lies behind the
    /// vertex along every incident edge, and reduces it to its exact sphere.
    /// Every other junction is refused; this does not trim or publish a solid.
    pub fn new(
        origin: Origin,
        strips: [&FilletDefinition; 3],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut supports: Vec<(Surface, bool)> = Vec::with_capacity(3);
        let mut indices = [[0; 2]; 3];
        let convex = strips[0].convex;
        for (strip, index) in strips.iter().zip(&mut indices) {
            if strip.convex != convex {
                return Err(bad(
                    "convex and concave strips at one vertex need a general corner patch",
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
        // A concave corner's ball sits outside the material, against the
        // inward normals.
        let normals: Vec<Vector3> = supports
            .iter()
            .map(|(s, same)| {
                let Surface::Plane(p) = s else { unreachable!() };
                unit(&p.normal).scaled(&Real::from_i64(if *same == convex { 1 } else { -1 }))
            })
            .collect();
        if sign(&normals[0].dot(&normals[1].cross(&normals[2])), dec)? == Sign::Zero {
            return Err(bad("an exact spherical corner needs independent supports"));
        }
        let radius = strips[0].radius.clone();
        // The ball centre lies `radius` short of every support along its
        // contact direction. On the first strip's axis it is fixed by the third
        // support; the other axes must pass through the same point.
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
                return Err(bad("the ball centre is not at radius from every support"));
            }
        }
        // The vertex is on the first edge, on the third support. A centre past
        // it along any edge would extend that strip beyond the vertex, which
        // neither evaluator's strip region covers.
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
        let directions = [0, 1, 2].map(|i| normals[i].clone());
        Ok(Self {
            origin,
            supports: [0, 1, 2].map(|i| supports[i].0.clone()),
            vertex,
            sector: SphereSector::new(centre, radius, directions, convex, dec)?,
            strips: indices,
        })
    }

    /// The ball centre.
    pub fn centre(&self) -> &Point3 {
        &self.sector.centre
    }

    /// The great circle shared with one strip, whose shorter arc runs from the
    /// strip's first contact foot to its second along the sector boundary.
    pub fn boundary(&self, strip: usize) -> Curve {
        let [a, b] = self.strips[strip];
        self.sector.boundary(a, b)
    }

    /// Whether a point's direction from the centre lies in the closed sector.
    pub fn in_sector(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        self.sector.in_sector(p, dec)
    }
}

/// One vertex of the inward-offset four-plane contact complex.
#[derive(Clone, Debug)]
pub struct FanVertex {
    /// Centre of the ball tangent to the three active supports.
    pub centre: Point3,
    /// The three active supports, in the fan's support order.
    pub contacts: [usize; 3],
}

/// The normal fan is either one four-contact ball or two trihedral balls
/// connected by the constant-radius sheet of their shared support pair.
#[derive(Clone, Debug)]
pub enum FanTopology {
    /// All four radius-offset planes meet at one feasible centre.
    Single {
        /// Centre of the common four-contact ball.
        centre: Point3,
    },
    /// Two feasible triple-contact vertices share one internal offset edge.
    Split {
        /// Endpoints of the internal edge, with their active support triples.
        vertices: [FanVertex; 2],
        /// Shared support indices; this pair is not an original selected edge.
        bridge: [usize; 2],
        /// The offset vertex terminating each original selected pair.
        endpoints: [usize; 4],
    },
}

/// Exact offset-contact subdivision of a fully selected planar four-face fan.
/// It describes contact equations and normal cells, not a constructed B-rep.
#[derive(Clone, Debug)]
pub struct PlanarFan {
    /// Original common vertex of the four source support planes.
    pub vertex: Point3,
    /// Unit directions from a contact ball centre to the four source supports.
    pub directions: [Vector3; 4],
    /// Source support indices bounding each selected run, in run order.
    pub pairs: [[usize; 2]; 4],
    /// Common positive contact radius.
    pub radius: Real,
    /// Certified subdivision and attachment of the original selected runs.
    pub topology: FanTopology,
}

impl PlanarFan {
    /// Certify the exposed four-pair cycle and its inward-offset subdivision.
    /// The source planes pass through `vertex` with unit contact `directions`;
    /// `pairs` indexes those directions. Unproved or singular fans refuse.
    pub fn new(
        vertex: Point3,
        directions: [Vector3; 4],
        pairs: [[usize; 2]; 4],
        radius: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if sign(&radius, dec)? != Sign::Positive {
            return Err(bad("a contact fan needs a positive radius"));
        }
        for normal in &directions {
            if sign(&(normal.norm_squared() - Real::one()), dec)? != Sign::Zero {
                return Err(bad("a contact fan needs unit material-side normals"));
            }
        }
        let mut degree = [0; 4];
        for (i, &[a, b]) in pairs.iter().enumerate() {
            if a >= 4
                || b >= 4
                || a == b
                || pairs[..i]
                    .iter()
                    .any(|pair| pair.contains(&a) && pair.contains(&b))
            {
                return Err(bad("four selected pairs must form a simple contact cycle"));
            }
            degree[a] += 1;
            degree[b] += 1;
            let boundary = directions[a].cross(&directions[b]);
            let mut side = Sign::Zero;
            for (j, normal) in directions.iter().enumerate() {
                if j == a || j == b {
                    continue;
                }
                let next = sign(&boundary.dot(normal), dec)?;
                if next == Sign::Zero || side != Sign::Zero && side != next {
                    return Err(bad("a selected pair is not exposed in the contact cone"));
                }
                side = next;
            }
        }
        if degree != [2; 4] {
            return Err(bad("four selected pairs must close the contact cone"));
        }
        let mut vertices = [None, None];
        let mut count = 0;
        for (omitted, contacts) in [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
            .into_iter()
            .enumerate()
        {
            let [a, b, c] = contacts.map(|i| &directions[i]);
            let bisector = b.minus(a).cross(&c.minus(a));
            let divisor = bisector.dot(a);
            if sign(&divisor, dec)? == Sign::Zero {
                continue;
            }
            let centre = vertex.plus(&bisector.scaled(&(-&radius / divisor)));
            let depth = &centre.minus(&vertex).dot(&directions[omitted]) + &radius;
            match sign(&depth, dec)? {
                Sign::Positive => continue,
                Sign::Zero => {
                    return Ok(Self {
                        vertex,
                        directions,
                        pairs,
                        radius,
                        topology: FanTopology::Single { centre },
                    });
                }
                Sign::Negative => {}
            }
            if count == vertices.len() {
                return Err(bad("the offset fan has more than two simple vertices"));
            }
            vertices[count] = Some(FanVertex { centre, contacts });
            count += 1;
        }
        let [Some(first), Some(second)] = vertices else {
            return Err(bad(
                "the four-plane offset fan has no certified subdivision",
            ));
        };
        let vertices = [first, second];
        let mut common = [None, None];
        let mut count = 0;
        for support in 0..4 {
            if vertices.iter().all(|v| v.contacts.contains(&support)) {
                if count == common.len() {
                    return Err(bad("two offset vertices must share exactly two supports"));
                }
                common[count] = Some(support);
                count += 1;
            }
        }
        let [Some(a), Some(b)] = common else {
            return Err(bad("two offset vertices must share a contact edge"));
        };
        let bridge = [a, b];
        if pairs
            .iter()
            .any(|pair| pair.contains(&a) && pair.contains(&b))
        {
            return Err(bad(
                "the internal contact edge must not be a selected input edge",
            ));
        }
        let boundary = directions[a].cross(&directions[b]);
        let mut sides = [Sign::Zero; 2];
        for (i, v) in vertices.iter().enumerate() {
            let other = *v.contacts.iter().find(|j| !bridge.contains(j)).unwrap();
            sides[i] = sign(&boundary.dot(&directions[other]), dec)?;
        }
        if !matches!(
            sides,
            [Sign::Positive, Sign::Negative] | [Sign::Negative, Sign::Positive]
        ) {
            return Err(bad(
                "the two spherical normal cells do not partition the fan",
            ));
        }
        let mut endpoints = [0; 4];
        for (i, [a, b]) in pairs.iter().enumerate() {
            let present = vertices
                .each_ref()
                .map(|v| v.contacts.contains(a) && v.contacts.contains(b));
            endpoints[i] = match present {
                [true, false] => 0,
                [false, true] => 1,
                _ => return Err(bad("a selected pair has no unique offset endpoint")),
            };
        }
        Ok(Self {
            vertex,
            directions,
            pairs,
            radius,
            topology: FanTopology::Split {
                vertices,
                bridge,
                endpoints,
            },
        })
    }
}
