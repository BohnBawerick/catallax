//! Elementary rolling-ball strips, retained with their source supports and edge.

use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Circle, Curve, Cylinder, Isometry, Line, Point3, Surface, Torus, Vector3};

use crate::profile::sign;

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
    /// Reduces a regular plane-plane or perpendicular plane-cylinder contact to
    /// its exact cylinder or torus. This does not trim or publish a solid.
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
            _ => {
                return Err(bad(
                    "elementary fillets need a line between planes or a circle between a plane and perpendicular cylinder",
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
        let (Surface::Cylinder(c), [Surface::Plane(a), Surface::Plane(b)]) =
            (&self.surface, &self.supports)
        else {
            return Err(bad("line-strip sector query requires two planes"));
        };
        let na = normal(&a.normal, self.same_sense[0]);
        let nb = normal(&b.normal, self.same_sense[1]);
        let d = p.minus(&c.origin);
        let radial = d
            .minus(&c.axis.scaled(&(d.dot(&c.axis) / c.axis.norm_squared())))
            .scaled(&Real::from_i64(if self.convex { 1 } else { -1 }));
        let dot = na.dot(&nb);
        let x = radial.dot(&na);
        let y = radial.dot(&nb);
        Ok(sign(&(&x - &(&dot * &y)), dec)? != Sign::Negative
            && sign(&(&y - &(&dot * &x)), dec)? != Sign::Negative)
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
    if sign(&(&major - radius), dec)? != Sign::Positive {
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
    Ok((
        Surface::Torus(Torus {
            center: centre.clone(),
            axis: cylinder.axis.clone(),
            major: major.clone(),
            minor: radius.clone(),
            ref_dir: edge.ref_dir.clone(),
        }),
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
        Surface::Torus(t) => {
            let axial = component(&t.axis);
            let radial = (Real::one() - axial.square() / t.axis.norm_squared()).sqrt();
            (coordinate(&t.center), &(&t.major * &radial) + &t.minor)
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
