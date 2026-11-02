//! Elementary rolling-ball strips, retained with their source supports and edge.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{
    Circle, Curve, Cylinder, Isometry, Line, Point3, Surface, Torus, TorusSheet, Vector3,
};

use crate::profile::sign;

#[path = "fillet_parallel.rs"]
mod parallel;

/// The exact generating definition of one elementary fillet strip.
#[derive(Clone, Debug)]
pub struct FilletDefinition {
    /// Persistent blend origin, independent of the radius.
    pub origin: Origin,
    /// The two input supports, in input-face order.
    pub supports: [Surface; 2],
    /// Whether each input face follows its support's natural normal.
    pub same_sense: [bool; 2],
    /// The selected input edge carrier, not the ball-centre path.
    pub edge: Curve,
    /// Positive rolling radius in millimetres.
    pub radius: Real,
    /// Material is removed for a convex edge and added for a concave edge.
    pub convex: bool,
    /// Exact elementary carrier of the shorter contact arc.
    pub surface: Surface,
    /// Contact curves, each paired with the corresponding input support.
    pub contacts: [Curve; 2],
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-contact", detail)
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
fn normal(v: &Vector3, same: bool) -> Vector3 {
    unit(v).scaled(&Real::from_i64(if same { 1 } else { -1 }))
}

impl FilletDefinition {
    /// Reduces plane-plane, perpendicular plane-cylinder and parallel-cylinder
    /// contacts to an exact cylinder or torus. This does not publish a solid.
    pub fn new(
        origin: Origin,
        supports: [Surface; 2],
        same_sense: [bool; 2],
        edge: Curve,
        radius: Real,
        convex: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if sign(&radius, dec)? != Sign::Positive {
            return Err(Refusal::new(
                "fillet.non-positive-radius",
                "fillet radius must be positive",
            ));
        }
        let k = Real::from_i64(if convex { 1 } else { -1 });
        let (surface, contacts) = match (&supports[0], &supports[1], &edge) {
            (Surface::Plane(a), Surface::Plane(b), Curve::Line(e)) => {
                if sign(&e.dir.norm_squared(), dec)? == Sign::Zero {
                    return Err(Refusal::new(
                        "fillet.invalid-edge",
                        "line direction is zero",
                    ));
                }
                for p in [a, b] {
                    if sign(&p.normal.norm_squared(), dec)? == Sign::Zero
                        || sign(&p.normal.dot(&e.dir), dec)? != Sign::Zero
                        || sign(&e.origin.minus(&p.origin).dot(&p.normal), dec)? != Sign::Zero
                    {
                        return Err(Refusal::new(
                            "fillet.invalid-edge",
                            "line does not lie in both nondegenerate supports",
                        ));
                    }
                }
                let na = normal(&a.normal, same_sense[0]);
                let nb = normal(&b.normal, same_sense[1]);
                let dot = na.dot(&nb);
                if sign(&(Real::one() - dot.square()), dec)? != Sign::Positive {
                    return Err(bad("parallel planes have no regular rolling-ball edge"));
                }
                let shift = na
                    .plus(&nb)
                    .scaled(&(-(&k * &radius) / (Real::one() + dot)));
                let centre = e.origin.plus(&shift);
                let feet = [&na, &nb].map(|n| centre.plus(&n.scaled(&(&k * &radius))));
                (
                    Surface::Cylinder(Cylinder {
                        origin: centre,
                        axis: e.dir.clone(),
                        radius_sq: radius.square(),
                        ref_dir: na,
                    }),
                    feet.map(|origin| {
                        Curve::Line(Line {
                            origin,
                            dir: e.dir.clone(),
                        })
                    }),
                )
            }
            (Surface::Plane(p), Surface::Cylinder(c), Curve::Circle(e)) => {
                plane_cylinder(p, c, e, same_sense, &radius, &k, dec)?
            }
            (Surface::Cylinder(c), Surface::Plane(p), Curve::Circle(e)) => {
                let (surface, [a, b]) =
                    plane_cylinder(p, c, e, [same_sense[1], same_sense[0]], &radius, &k, dec)?;
                (surface, [b, a])
            }
            (Surface::Cylinder(a), Surface::Cylinder(b), Curve::Line(e)) => {
                parallel::cylinders(a, b, e, same_sense, &radius, &k, &origin.step, dec)?
            }
            _ => {
                return Err(bad(
                    "fillets need a line between planes or parallel cylinders, or a circle between a plane and perpendicular cylinder",
                ));
            }
        };
        Ok(Self {
            origin,
            supports,
            same_sense,
            edge,
            radius,
            convex,
            surface,
            contacts,
        })
    }
    /// Moves the retained definition and its exact reduction through an isometry.
    pub fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        let edge = map
            .curve(&self.edge)
            .ok_or_else(|| bad("edge has no exact isometry"))?;
        let a = map
            .curve(&self.contacts[0])
            .ok_or_else(|| bad("contact has no exact isometry"))?;
        let b = map
            .curve(&self.contacts[1])
            .ok_or_else(|| bad("contact has no exact isometry"))?;
        self.edge = edge;
        self.contacts = [a, b];
        self.supports = self.supports.each_ref().map(|s| map.surface(s));
        self.surface = map.surface(&self.surface);
        Ok(())
    }

    /// Proves the complete blend carrier disjoint from another support by
    /// analytic bounds. `false` leaves intersection or trimmed-domain work open.
    pub fn avoids_support(&self, surface: &Surface, dec: &mut Decider) -> Result<bool, Refusal> {
        separated(&self.surface, surface, dec)
    }

    /// Whether a cylinder point lies in the shorter arc's normal sector.
    /// End trimming is a separate topological constraint.
    pub fn in_sector(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let (Surface::Cylinder(c), [Curve::Line(a), Curve::Line(b)]) =
            (&self.surface, &self.contacts)
        else {
            return Err(bad(
                "line-strip sector query requires straight contact rails",
            ));
        };
        let (na, nb, sector_sign) = match &self.supports {
            [Surface::Plane(a), Surface::Plane(b)] => (
                normal(&a.normal, self.same_sense[0]),
                normal(&b.normal, self.same_sense[1]),
                Real::from_i64(if self.convex { 1 } else { -1 }),
            ),
            [Surface::Cylinder(_), Surface::Cylinder(_)] => (
                a.origin
                    .minus(&c.origin)
                    .scaled(&(&Real::one() / &self.radius)),
                b.origin
                    .minus(&c.origin)
                    .scaled(&(&Real::one() / &self.radius)),
                Real::one(),
            ),
            _ => return Err(bad("line-strip sector has no supported contact pair")),
        };
        let d = p.minus(&c.origin);
        let radial = d
            .minus(&c.axis.scaled(&(d.dot(&c.axis) / c.axis.norm_squared())))
            .scaled(&sector_sign);
        let dot = na.dot(&nb);
        let x = radial.dot(&na);
        let y = radial.dot(&nb);
        Ok(sign(&(&x - &(&dot * &y)), dec)? != Sign::Negative
            && sign(&(&y - &(&dot * &x)), dec)? != Sign::Negative)
    }

    /// Whether this strip and another consume overlapping arcs of one shared
    /// cylinder wall, given each strip's support index for that wall. Rails
    /// that pass each other can leave both feet in the open face.
    pub(crate) fn wall_arc_meets(
        &self,
        support: usize,
        other: &Self,
        other_support: usize,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        parallel::wall_arc_meets([(self, support), (other, other_support)], dec)
    }

    /// Proves two parallel plane-plane strips' blend regions disjoint, although
    /// their complete carriers may meet. Each region lies in the prism over its
    /// cross-section kite: ball centre, both contacts and the edge. Kites are
    /// convex, so one of their side normals separates disjoint pairs. A curved
    /// support can bulge out of the kite, so those strips get no such proof.
    pub(crate) fn sectors_apart(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        let kite = |strip: &Self| {
            let (
                Surface::Cylinder(c),
                [Curve::Line(a), Curve::Line(b)],
                Curve::Line(e),
                [Surface::Plane(_), Surface::Plane(_)],
            ) = (
                &strip.surface,
                &strip.contacts,
                &strip.edge,
                &strip.supports,
            )
            else {
                return None;
            };
            Some((
                c.axis.clone(),
                [&c.origin, &a.origin, &e.origin, &b.origin].map(Point3::clone),
            ))
        };
        let (Some((axis, first)), Some((other_axis, second))) = (kite(self), kite(other)) else {
            return Ok(false);
        };
        if sign(&axis.cross(&other_axis).norm_squared(), dec)? != Sign::Zero {
            return Ok(false);
        }
        for points in [&first, &second] {
            for i in 0..4 {
                // Every side normal is perpendicular to the common axis, so
                // projections ignore the kites' axial stations.
                let side = points[(i + 1) % 4].minus(&points[i]);
                let normal = axis.cross(&side);
                let range = |kite: &[Point3; 4], dec: &mut Decider| -> Result<_, Refusal> {
                    let values = kite.each_ref().map(|p| p.minus(&points[i]).dot(&normal));
                    let mut low = values[0].clone();
                    let mut high = values[0].clone();
                    for value in &values[1..] {
                        if sign(&(value - &low), dec)? == Sign::Negative {
                            low = value.clone();
                        }
                        if sign(&(value - &high), dec)? == Sign::Positive {
                            high = value.clone();
                        }
                    }
                    Ok((low, high))
                };
                let (al, ah) = range(&first, dec)?;
                let (bl, bh) = range(&second, dec)?;
                if sign(&(&bl - &ah), dec)? == Sign::Positive
                    || sign(&(&al - &bh), dec)? == Sign::Positive
                {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

fn plane_cylinder(
    plane: &catallax_shapes::Plane,
    cylinder: &Cylinder,
    edge: &Circle,
    same: [bool; 2],
    radius: &Real,
    k: &Real,
    dec: &mut Decider,
) -> Result<(Surface, [Curve; 2]), Refusal> {
    if sign(&plane.normal.norm_squared(), dec)? == Sign::Zero
        || sign(&cylinder.axis.norm_squared(), dec)? == Sign::Zero
        || sign(&cylinder.radius_sq, dec)? != Sign::Positive
        || sign(&edge.normal.norm_squared(), dec)? == Sign::Zero
        || sign(&edge.ref_dir.norm_squared(), dec)? == Sign::Zero
        || sign(&edge.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero
        || sign(&edge.ref_dir.dot(&edge.normal), dec)? != Sign::Zero
        || sign(
            &edge
                .center
                .minus(&cylinder.origin)
                .cross(&cylinder.axis)
                .norm_squared(),
            dec,
        )? != Sign::Zero
        || sign(&edge.center.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
        || sign(&(&edge.radius_sq - &cylinder.radius_sq), dec)? != Sign::Zero
    {
        return Err(Refusal::new(
            "fillet.invalid-edge",
            "circle does not lie in both nondegenerate supports",
        ));
    }
    if sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero {
        return Err(bad("oblique plane-cylinder contacts are not toroidal"));
    }
    let n = normal(&plane.normal, same[0]);
    let centre = edge.center.plus(&n.scaled(&(-(k * radius))));
    let shift = &(k * &Real::from_i64(if same[1] { 1 } else { -1 })) * radius;
    let major = &cylinder.radius() - &shift;
    // A wall contact on the axis side needs the ring torus. One away from
    // the axis keeps the tube angles of the outer quarter, where the distance
    // from the axis is at least `major`; for `0 < major <= radius` that
    // quarter lies on the regular outer sheet of a horn or spindle torus.
    let ring = sign(&(&major - radius), dec)? == Sign::Positive;
    if !ring && sign(&major, dec)? != Sign::Positive {
        return Err(bad(
            "this contact requires a horn or spindle torus and its sheet continuation",
        ));
    }
    let circle = |center: Point3, radius_sq: Real| {
        Curve::Circle(Circle {
            center,
            normal: edge.normal.clone(),
            radius_sq,
            ref_dir: edge.ref_dir.clone(),
        })
    };
    let surface = if ring {
        Surface::Torus(Torus {
            center: centre.clone(),
            axis: cylinder.axis.clone(),
            major: major.clone(),
            minor: radius.clone(),
            ref_dir: edge.ref_dir.clone(),
        })
    } else {
        Surface::TorusSheet(TorusSheet {
            center: centre.clone(),
            axis: cylinder.axis.clone(),
            major: major.clone(),
            minor: radius.clone(),
            ref_dir: edge.ref_dir.clone(),
        })
    };
    Ok((
        surface,
        [
            circle(edge.center.clone(), major.square()),
            circle(centre, cylinder.radius_sq.clone()),
        ],
    ))
}
fn bounds(
    surface: &Surface,
    axis: usize,
    dec: &mut Decider,
) -> Result<Option<(Real, Real)>, Refusal> {
    let component =
        |v: &catallax_shapes::Vector3| [v.x.clone(), v.y.clone(), v.z.clone()][axis].clone();
    let coordinate = |p: &Point3| [p.x.clone(), p.y.clone(), p.z.clone()][axis].clone();
    let (centre, extent) = match surface {
        Surface::Plane(p) => {
            let n = [&p.normal.x, &p.normal.y, &p.normal.z];
            for (i, value) in n.into_iter().enumerate() {
                if i != axis && sign(value, dec)? != Sign::Zero {
                    return Ok(None);
                }
            }
            (coordinate(&p.origin), Real::zero())
        }
        Surface::Cylinder(c) => {
            if sign(&component(&c.axis), dec)? != Sign::Zero {
                return Ok(None);
            }
            (coordinate(&c.origin), c.radius())
        }
        Surface::Sphere(s) => (coordinate(&s.center), s.radius.clone()),
        Surface::Torus(Torus {
            center,
            axis,
            major,
            minor,
            ..
        })
        | Surface::TorusSheet(TorusSheet {
            center,
            axis,
            major,
            minor,
            ..
        }) => {
            let axial = component(axis);
            let radial = (Real::one() - axial.square() / axis.norm_squared()).sqrt();
            (coordinate(center), &(major * &radial) + minor)
        }
        _ => return Ok(None),
    };
    Ok(Some((&centre - &extent, &centre + &extent)))
}

pub(crate) fn separated(a: &Surface, b: &Surface, dec: &mut Decider) -> Result<bool, Refusal> {
    for axis in 0..3 {
        if let (Some((al, ah)), Some((bl, bh))) = (bounds(a, axis, dec)?, bounds(b, axis, dec)?)
            && (sign(&(&bl - &ah), dec)? == Sign::Positive
                || sign(&(&al - &bh), dec)? == Sign::Positive)
        {
            return Ok(true);
        }
    }
    let pair = match (a, b) {
        (Surface::Cylinder(c), Surface::Torus(t)) => Some((c, t, b)),
        (Surface::Torus(t), Surface::Cylinder(c)) => Some((c, t, a)),
        _ => None,
    };
    if let Some((c, t, torus)) = pair {
        // Project the centre circle's box perpendicular to the cylinder axis.
        // Its distance is a lower bound; the two tube radii are added exactly.
        let mut lower_squared = Real::zero();
        for axis in 0..3 {
            let axial = [&c.axis.x, &c.axis.y, &c.axis.z][axis];
            if sign(axial, dec)? != Sign::Zero {
                continue;
            }
            let coordinate = [&c.origin.x, &c.origin.y, &c.origin.z][axis];
            let (lo, hi) = bounds(torus, axis, dec)?.expect("torus is bounded");
            let lo = &lo + &t.minor;
            let hi = &hi - &t.minor;
            let before = &lo - coordinate;
            let after = coordinate - &hi;
            let gap = if sign(&before, dec)? == Sign::Positive {
                before
            } else if sign(&after, dec)? == Sign::Positive {
                after
            } else {
                Real::zero()
            };
            lower_squared = lower_squared + gap.square();
        }
        if sign(&(lower_squared - (&c.radius() + &t.minor).square()), dec)? == Sign::Positive {
            return Ok(true);
        }
    }
    Ok(false)
}
