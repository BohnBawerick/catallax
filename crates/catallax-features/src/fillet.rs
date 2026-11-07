//! Elementary rolling-ball strips, retained with their source supports and edge.

use catallax_numbers::{Decider, OperationSpan, Rational, Real, Refusal, Sign};
use catallax_recipe::{Origin, StepId};
use catallax_shapes::{
    BallSide, Circle, Curve, Cylinder, Ellipse, Isometry, Line, Plane, Point3, RollingPipe,
    RollingSaddle, RollingSaddleDefinition, RollingTorus, Sphere, SplineSupport, Surface, Torus,
    TorusSheet, Vector3,
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
    /// The selected input edge carrier, not the ball-centre path. A ring between
    /// two curved supports keeps its intersection curve, whose ids name the ring
    /// in the store it came from: moved copies carry their geometry in `supports`
    /// and `surface`.
    pub edge: Curve,
    /// Positive rolling radius in millimetres.
    pub radius: Real,
    /// A larger declared radius whose contact branch ended at `radius`.
    /// `None` means this sheet is at the declared radius. `Some` retains the
    /// authored request while the exact homotopy keeps its last exposed sheet.
    pub requested_radius: Option<Real>,
    /// Material is removed for a convex edge and added for a concave edge.
    pub convex: bool,
    /// Exact elementary carrier of the shorter contact arc.
    pub surface: Surface,
    /// Contact curves, each paired with the corresponding input support. On
    /// a sphere strip, where the ball is as wide as its cylindrical support,
    /// the plane contact is a circle of zero radius: the one point where the
    /// cylinder's axis meets the plane.
    pub contacts: [Curve; 2],
}

/// The refusal code for a failed strip at or beyond a proved regular contact
/// limit, or whose contact with one face lies wholly outside it. Roll-over
/// geometry can still exist; this is not proof that the recipe is impossible.
pub const RADIUS_TOO_LARGE: &str = "fillet.radius-too-large";

/// The refusal code for a strip end the builder does not trim yet: an arc
/// strip on an end plane off its axis or on a curved end face, an end face
/// the rails do not cross, or a vertex of more than three faces.
pub const UNSUPPORTED_END: &str = "fillet.unsupported-end";

/// The refusal code for a junction of selected strips the builder does not
/// close yet: convex and concave strips at one vertex, curved strips meeting
/// at a sharp vertex, three strips without an exact spherical corner, or
/// incomplete, clipped or otherwise unproved equal-cylinder pole fans.
pub const UNSUPPORTED_CORNER: &str = "fillet.unsupported-corner";

// Marks the bound in a radius refusal's detail.
const LIMIT_PREFIX: &str = "; radius limit ";
const LIMIT_SUFFIX: &str = " mm";

/// The bound a [`RADIUS_TOO_LARGE`] refusal carries, in millimetres, when it
/// has one. The refusal proves every radius above an exact bound too large for
/// the same edge; this is that bound rounded down to three decimal places. A
/// radius below it is not proved to fit. `None` for other refusals and when
/// the refusal has no bound.
pub fn radius_limit(refusal: &Refusal) -> Option<Rational> {
    if refusal.code != RADIUS_TOO_LARGE {
        return None;
    }
    let (_, tail) = refusal.detail.rsplit_once(LIMIT_PREFIX)?;
    Rational::parse(tail.strip_suffix(LIMIT_SUFFIX)?).ok()
}

/// A [`RADIUS_TOO_LARGE`] refusal. A limit is shown rounded down to three
/// decimal places, and dropped when that leaves nothing positive.
pub(crate) fn too_large(reason: &str, limit: Option<&Real>, dec: &mut Decider) -> Refusal {
    let shown = limit.and_then(|limit| {
        let low = match limit.as_rational() {
            Some(exact) => exact,
            None => dec
                .enclose(limit, &Rational::from_ratio(1, 1_000_000)?)
                .ok()?
                .lo(),
        };
        let thousandths = Rational::from_integer((&low * &Rational::from_i64(1000)).floor());
        let shown = &thousandths / &Rational::from_i64(1000);
        (shown.sign() == Sign::Positive)
            .then(|| shown.to_decimal_exact())
            .flatten()
    });
    Refusal::new(
        RADIUS_TOO_LARGE,
        match shown {
            Some(shown) => format!("{reason}{LIMIT_PREFIX}{shown}{LIMIT_SUFFIX}"),
            None => reason.to_string(),
        },
    )
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-contact", detail)
}
/// The square root, exact for a rational square.
fn root(value: &Real) -> Real {
    value
        .as_rational()
        .and_then(|r| crate::extrude::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| value.sqrt())
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &root(&v.norm_squared())))
}
fn normal(v: &Vector3, same: bool) -> Vector3 {
    unit(v).scaled(&Real::from_i64(if same { 1 } else { -1 }))
}

impl FilletDefinition {
    /// Reduces plane-plane, parallel-cylinder and cylinder-plane generator
    /// contacts to an exact cylinder, perpendicular plane-cylinder rings and
    /// other circles where two surfaces of revolution about one axis meet to
    /// an exact torus, a plane-cylinder ring to an exact sphere when the ball
    /// is as wide as the cylinder, and oblique plane-cylinder ellipses to an
    /// exact rational band. This does not publish a solid.
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
            (Surface::Plane(p), Surface::Cylinder(c), Curve::Ellipse(e)) => {
                oblique(p, c, e, same_sense, &radius, convex, &origin.step, dec)?
            }
            (Surface::Cylinder(c), Surface::Plane(p), Curve::Ellipse(e)) => {
                let (surface, [a, b]) = oblique(
                    p,
                    c,
                    e,
                    [same_sense[1], same_sense[0]],
                    &radius,
                    convex,
                    &origin.step,
                    dec,
                )?;
                (surface, [b, a])
            }
            (Surface::Plane(p), Surface::Cylinder(c), Curve::Line(e)) => {
                parallel::plane_cylinder(p, c, e, same_sense, &radius, &k, &origin.step, dec)?
            }
            (Surface::Cylinder(c), Surface::Plane(p), Curve::Line(e)) => {
                let (surface, [a, b]) = parallel::plane_cylinder(
                    p,
                    c,
                    e,
                    [same_sense[1], same_sense[0]],
                    &radius,
                    &k,
                    &origin.step,
                    dec,
                )?;
                (surface, [b, a])
            }
            (Surface::Cylinder(a), Surface::Cylinder(b), Curve::Line(e)) => {
                parallel::cylinders(a, b, e, same_sense, &radius, &k, &origin.step, dec)?
            }
            (a, b, Curve::Circle(e)) => {
                coaxial([a, b], e, same_sense, &radius, convex, &origin.step, dec)?
            }
            _ => {
                return Err(bad(
                    "fillets need a line between planes, parallel cylinders or a cylinder and a parallel plane, a plane section of a cylinder, or a circle where two surfaces of revolution about one axis meet",
                ));
            }
        };
        Ok(Self {
            origin,
            supports,
            same_sense,
            edge,
            radius,
            requested_radius: None,
            convex,
            surface,
            contacts,
        })
    }
    /// The strip through `witness` on a cylinder contact: a regular closed
    /// intersection ring, or a pole-to-pole ellipse arc of equal cylinders
    /// whose offsets remain inside the centre ellipse's regular tube.
    /// Contacts are spline rails in the carrier's chart direction.
    #[allow(clippy::too_many_arguments)]
    pub fn ring(
        origin: Origin,
        supports: [Surface; 2],
        same_sense: [bool; 2],
        edge: Curve,
        witness: &Point3,
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
        let (Surface::Cylinder(first), Surface::Cylinder(second)) = (&supports[0], &supports[1])
        else {
            return Err(bad("a cylinder strip needs two cylinders"));
        };
        if !matches!(edge, Curve::Intersection(_) | Curve::Ellipse(_)) {
            return Err(bad(
                "a cylinder strip needs an intersection or elliptic arc",
            ));
        }
        let (surface, contacts) = saddle(
            first,
            second,
            same_sense,
            witness,
            &radius,
            convex,
            &origin.step,
            dec,
        )?;
        if let Curve::Ellipse(ellipse) = &edge {
            let Surface::BSpline(SplineSupport::RollingSaddle(band)) = &surface else {
                return Err(bad("an elliptic cylinder strip has no rolling band"));
            };
            if !band.is_edge(ellipse, dec)? {
                return Err(bad(
                    "the input ellipse is not this crossing band's source edge",
                ));
            }
        }
        Ok(Self {
            origin,
            supports,
            same_sense,
            edge,
            radius,
            requested_radius: None,
            convex,
            surface,
            contacts,
        })
    }

    /// The index of a contact that collapses to one point: the plane contact
    /// of a sphere strip.
    pub fn point_contact(&self) -> Option<usize> {
        if !matches!(self.surface, Surface::Sphere(_)) {
            return None;
        }
        self.supports
            .iter()
            .position(|support| matches!(support, Surface::Plane(_)))
    }

    /// Moves the retained definition and its exact reduction through an isometry.
    pub fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        let edge = match &self.edge {
            Curve::Intersection(_) => self.edge.clone(),
            edge => map
                .curve(edge)
                .ok_or_else(|| bad("edge has no exact isometry"))?,
        };
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
            _ => (
                a.origin
                    .minus(&c.origin)
                    .scaled(&(&Real::one() / &self.radius)),
                b.origin
                    .minus(&c.origin)
                    .scaled(&(&Real::one() / &self.radius)),
                Real::one(),
            ),
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

    /// The points where the edge and both contact rails of a straight strip
    /// between planes cross two cutting planes, each given by a point and a
    /// normal, the first three on the first. Each cross-section of the blend
    /// region lies in the triangle of its edge point and contacts, whose
    /// support lines are the ball's tangents there, so the region between the
    /// cuts lies in the convex hull of these points. `None` for other strips
    /// and for a cut parallel to the edge.
    pub fn cut_hull(
        &self,
        cuts: [(&Point3, &Vector3); 2],
        dec: &mut Decider,
    ) -> Result<Option<[Point3; 6]>, Refusal> {
        let (
            [Surface::Plane(_), Surface::Plane(_)],
            Curve::Line(edge),
            [Curve::Line(a), Curve::Line(b)],
        ) = (&self.supports, &self.edge, &self.contacts)
        else {
            return Ok(None);
        };
        let mut points = Vec::with_capacity(6);
        for (origin, normal) in cuts {
            let divisor = edge.dir.dot(normal);
            if sign(&divisor, dec)? == Sign::Zero {
                return Ok(None);
            }
            // The rails are parallel to the edge.
            for line in [edge, a, b] {
                let along = &origin.minus(&line.origin).dot(normal) / &divisor;
                points.push(line.origin.plus(&edge.dir.scaled(&along)));
            }
        }
        Ok(points.try_into().ok())
    }

    /// Points whose convex hull holds the region of a ring strip over one arc
    /// of its edge, from `from` to `to` turning positively about the edge
    /// circle's normal, both on the edge. Each meridian section of the region
    /// lies in the triangle of its edge point and contacts, as for
    /// [`Self::cut_hull`]. The edge and both contacts are coaxial circles, and
    /// each of their arcs lies in the triangle of its ends and the meeting
    /// point of their tangents. `None` for other strips and for an arc of half
    /// a turn or more.
    pub fn arc_hull(
        &self,
        [from, to]: [&Point3; 2],
        dec: &mut Decider,
    ) -> Result<Option<Vec<Point3>>, Refusal> {
        let (Curve::Circle(edge), [Curve::Circle(a), Curve::Circle(b)]) =
            (&self.edge, &self.contacts)
        else {
            return Ok(None);
        };
        let [u, w] = [from, to].map(|p| p.minus(&edge.center));
        if sign(&u.cross(&w).dot(&edge.normal), dec)? != Sign::Positive {
            return Ok(None);
        }
        let size = root(&edge.radius_sq);
        let [u, w] = [u, w].map(|v| v.scaled(&(&Real::one() / &size)));
        let lift = &Real::one() + &u.dot(&w);
        let mut points = Vec::with_capacity(9);
        for circle in [edge, a, b] {
            let radius = root(&circle.radius_sq);
            points.push(circle.center.plus(&u.scaled(&radius)));
            points.push(circle.center.plus(&w.scaled(&radius)));
            points.push(circle.center.plus(&u.plus(&w).scaled(&(&radius / &lift))));
        }
        Ok(Some(points))
    }

    /// Whether this strip and another consume overlapping bands of one shared
    /// support, given each strip's support index for it. Rails that pass each
    /// other can leave both feet in the open face. Plane-plane pairs are left
    /// to their carrier and kite proofs.
    pub(crate) fn bands_meet(
        &self,
        support: usize,
        other: &Self,
        other_support: usize,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        parallel::bands_meet([(self, support), (other, other_support)], dec)
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

/// The exact band of a closed contact between a plane and an oblique cylinder.
/// Rails keep the input edge's direction, so they replace it in each face loop.
#[allow(clippy::too_many_arguments)]
fn oblique(
    plane: &Plane,
    cylinder: &Cylinder,
    edge: &Ellipse,
    same: [bool; 2],
    radius: &Real,
    convex: bool,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Surface, [Curve; 2]), Refusal> {
    // The law builds its carrier's sections along with its contacts.
    let contact_span = OperationSpan::begin("brep", step, "fillet.contact", dec);
    let pipe = RollingPipe::new(plane, cylinder, same, radius.clone(), convex, dec).map_err(
        |refusal| {
            if refusal.code == "shapes.rolling-pipe" {
                bad(&refusal.detail)
            } else {
                refusal
            }
        },
    )?;
    if !pipe.is_edge(edge, dec)? {
        return Err(Refusal::new(
            "fillet.invalid-edge",
            "ellipse is not the section of this plane and cylinder",
        ));
    }
    contact_span.finish(dec, None, None);
    let section_span = OperationSpan::begin("brep", step, "fillet.sections", dec);
    let along_edge = |mut rail: Ellipse, dec: &mut Decider| -> Result<Curve, Refusal> {
        if sign(&rail.normal.dot(&edge.normal), dec)? == Sign::Negative {
            rail.normal = rail.normal.scaled(&Real::from_i64(-1));
        }
        Ok(Curve::Ellipse(rail))
    };
    let rails = [
        along_edge(pipe.plane_rail(dec)?, dec)?,
        along_edge(pipe.wall_rail(dec)?, dec)?,
    ];
    section_span.finish(dec, None, None);
    Ok((Surface::BSpline(SplineSupport::RollingPipe(pipe)), rails))
}

// The ball's side of a cylinder: a convex edge puts it in the material, and a
// face in the sense of the cylinder's normal has its material inside.
pub(crate) fn ball_side(same_sense: bool, convex: bool) -> BallSide {
    if same_sense == convex {
        BallSide::Inside
    } else {
        BallSide::Outside
    }
}

/// The exact torus of a circle where two surfaces of revolution about its axis
/// meet, other than a plane and a cylinder. Contacts keep the edge's direction.
fn coaxial(
    supports: [&Surface; 2],
    edge: &Circle,
    same: [bool; 2],
    radius: &Real,
    convex: bool,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Surface, [Curve; 2]), Refusal> {
    // The law builds its carrier's sections along with its contacts.
    let contact_span = OperationSpan::begin("brep", step, "fillet.contact", dec);
    let sides = same.map(|same| ball_side(same, convex));
    let torus = match RollingTorus::new(supports, sides, edge, radius, dec) {
        Ok(torus) => torus,
        Err(refusal) if refusal.code == "shapes.rolling-torus-edge" => {
            return Err(Refusal::new("fillet.invalid-edge", refusal.detail));
        }
        Err(refusal) if refusal.code == "shapes.rolling-torus" => {
            let Some(limit) = RollingTorus::radius_limit(supports, sides, edge, dec)? else {
                return Err(bad(&refusal.detail));
            };
            return Err(match sign(&(radius - &limit), dec)? {
                Sign::Positive => too_large(
                    "no ball of this radius rolls along the edge against both faces",
                    Some(&limit),
                    dec,
                ),
                Sign::Zero => too_large(
                    "the ball reaches this edge's regular contact limit",
                    Some(&limit),
                    dec,
                ),
                Sign::Negative => bad(&refusal.detail),
            });
        }
        Err(refusal) => return Err(refusal),
    };
    contact_span.finish(dec, None, None);
    Ok((torus.surface, torus.contacts.map(Curve::Circle)))
}

#[allow(clippy::too_many_arguments)]
fn saddle(
    first: &Cylinder,
    second: &Cylinder,
    same: [bool; 2],
    witness: &Point3,
    radius: &Real,
    convex: bool,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Surface, [Curve; 2]), Refusal> {
    let unsupported = |refusal: Refusal| {
        if refusal.code == "shapes.rolling-saddle" {
            bad(&refusal.detail)
        } else {
            refusal
        }
    };
    let first_wide = sign(&(&first.radius_sq - &second.radius_sq), dec)? != Sign::Negative;
    let ([wide, narrow], sides) = if first_wide {
        ([first, second], [same[0], same[1]])
    } else {
        ([second, first], [same[1], same[0]])
    };
    // The law builds its carrier's sections along with its contacts.
    let contact_span = OperationSpan::begin("brep", step, "fillet.contact", dec);
    let definition = RollingSaddleDefinition::crossing(
        wide,
        narrow,
        radius,
        sides.map(|same| ball_side(same, convex)),
        witness,
        dec,
    )
    .map_err(unsupported)?;
    contact_span.finish(dec, None, None);
    let fitting_span = (!matches!(
        definition.law().ring,
        catallax_shapes::SaddleRing::Crossing { .. }
    ))
    .then(|| OperationSpan::begin("brep", step, "fillet.fitting", dec));
    let saddle = RollingSaddle::fit(definition, dec).map_err(unsupported)?;
    if let Some(fitting_span) = fitting_span {
        fitting_span.finish(dec, None, None);
    }
    let [rail_wide, rail_narrow] = saddle.rails().map(Curve::BSpline);
    let rails = if first_wide {
        [rail_wide, rail_narrow]
    } else {
        [rail_narrow, rail_wide]
    };
    Ok((
        Surface::BSpline(SplineSupport::RollingSaddle(saddle)),
        rails,
    ))
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
    // quarter lies on the regular outer sheet of a horn or spindle torus. A
    // ball as wide as the cylinder sits on its axis and does not move: its
    // plane contact collapses to the point where the axis meets the plane.
    let ring = sign(&(&major - radius), dec)? == Sign::Positive;
    let sheet = if ring {
        Sign::Positive
    } else {
        sign(&major, dec)?
    };
    let collapsed = sheet == Sign::Zero;
    if sheet == Sign::Negative {
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
    // Either axis direction gives the same torus. A written torus turns its
    // tube angle from the outer equator toward its axis, and a strip whose
    // wall contact lies on that equator, with the ball inside the wall, is
    // charted only on the axis side of it: point the axis from the ball
    // centre to the plane contact there.
    let axis = if sign(&shift, dec)? == Sign::Positive
        && sign(&cylinder.axis.dot(&n.scaled(k)), dec)? == Sign::Negative
    {
        cylinder.axis.scaled(&Real::from_i64(-1))
    } else {
        cylinder.axis.clone()
    };
    let surface = if collapsed {
        Surface::Sphere(Sphere::new(
            centre.clone(),
            radius.clone(),
            n.scaled(k),
            edge.ref_dir.clone(),
            dec,
        )?)
    } else if ring {
        Surface::Torus(Torus {
            center: centre.clone(),
            axis,
            major: major.clone(),
            minor: radius.clone(),
            ref_dir: edge.ref_dir.clone(),
        })
    } else {
        Surface::TorusSheet(TorusSheet {
            center: centre.clone(),
            axis,
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
        Surface::BSpline(SplineSupport::RollingSaddle(saddle)) => {
            let [low, high] = saddle.extent(axis, dec)?;
            return Ok(Some((Real::from(low), Real::from(high))));
        }
        // Positive weights keep the band in the hull of its controls.
        Surface::BSpline(SplineSupport::RollingPipe(pipe)) => {
            let mut range: Option<(Rational, Rational)> = None;
            for (control, _) in pipe.enclosed_controls(dec)? {
                let (low, high) = (control[axis].lo(), control[axis].hi());
                range = Some(match range {
                    Some((least, most)) => (least.min(low), most.max(high)),
                    None => (low, high),
                });
            }
            return Ok(range.map(|(low, high)| (Real::from(low), Real::from(high))));
        }
        _ => return Ok(None),
    };
    Ok(Some((&centre - &extent, &centre + &extent)))
}

/// Proves a surface clear of the convex hull of points: a plane with every
/// point strictly on one side, or a surface whose bounds miss the points' box
/// along one axis. `false` leaves the question open.
pub fn hull_clear(
    points: &[Point3],
    surface: &Surface,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if let Surface::Plane(plane) = surface {
        let mut sides = Vec::with_capacity(points.len());
        for p in points {
            sides.push(sign(&p.minus(&plane.origin).dot(&plane.normal), dec)?);
        }
        return Ok(sides.iter().all(|s| *s == Sign::Positive)
            || sides.iter().all(|s| *s == Sign::Negative));
    }
    let Some(first) = points.first() else {
        return Ok(true);
    };
    for axis in 0..3 {
        let Some((low, high)) = bounds(surface, axis, dec)? else {
            continue;
        };
        let coordinate = |p: &Point3| [&p.x, &p.y, &p.z][axis].clone();
        let mut least = coordinate(first);
        let mut most = least.clone();
        for p in &points[1..] {
            let value = coordinate(p);
            if sign(&(&value - &least), dec)? == Sign::Negative {
                least = value.clone();
            }
            if sign(&(&value - &most), dec)? == Sign::Positive {
                most = value;
            }
        }
        if sign(&(&low - &most), dec)? == Sign::Positive
            || sign(&(&least - &high), dec)? == Sign::Positive
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Proves two hull point sets apart by a plane through three points of
/// either set: its own points on one side, the other's strictly on the other
/// side. A [`FilletDefinition::cut_hull`] tries its five faces; other sets,
/// such as [`FilletDefinition::arc_hull`] gives, try every triple, which
/// includes every face of their hull.
pub(crate) fn hulls_apart(
    first: &[Point3],
    second: &[Point3],
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    // Both cuts, both supports and the chord plane through the two rails.
    const FACES: [[usize; 3]; 5] = [[0, 1, 2], [3, 4, 5], [0, 1, 3], [0, 2, 3], [1, 2, 4]];
    for (own, other) in [(first, second), (second, first)] {
        let n = own.len();
        let faces: Vec<[usize; 3]> = if n == 6 {
            FACES.to_vec()
        } else {
            (0..n)
                .flat_map(|i| (i + 1..n).flat_map(move |j| (j + 1..n).map(move |k| [i, j, k])))
                .collect()
        };
        'faces: for [i, j, k] in faces {
            let normal = own[j].minus(&own[i]).cross(&own[k].minus(&own[i]));
            let mut side = Sign::Zero;
            for p in own {
                match (sign(&p.minus(&own[i]).dot(&normal), dec)?, side) {
                    (Sign::Zero, _) => {}
                    (s, Sign::Zero) => side = s,
                    (s, side) if s != side => continue 'faces,
                    _ => {}
                }
            }
            if side == Sign::Zero {
                continue;
            }
            let mut beyond = true;
            for p in other {
                let s = sign(&p.minus(&own[i]).dot(&normal), dec)?;
                beyond &= s != Sign::Zero && s != side;
            }
            if beyond {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

/// Proves two hull point sets that share the point `touch` meet nowhere
/// else: a plane through it keeps the first set on one closed side and the
/// second, but for points at `touch`, strictly on the other, so the second
/// hull meets that plane only there. The normals of planes through `touch`
/// with each set on its own closed side form a cone. Its extreme rays are
/// normal to two directions from `touch`, and their sum is strict wherever
/// any normal in the cone is.
pub(crate) fn hulls_touch(
    first: &[Point3],
    second: &[Point3],
    touch: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let mut directions = Vec::with_capacity(first.len() + second.len());
    for (points, own) in [(first, true), (second, false)] {
        for p in points {
            let d = p.minus(touch);
            if sign(&d.norm_squared(), dec)? != Sign::Zero {
                directions.push((d, own));
            }
        }
    }
    // Whether each set stays on its own closed side, the second strictly
    // when asked.
    let holds = |n: &Vector3, strict: bool, dec: &mut Decider| -> Result<bool, Refusal> {
        for (d, own) in &directions {
            match (sign(&n.dot(d), dec)?, own) {
                (Sign::Negative, true) | (Sign::Positive, false) => return Ok(false),
                (Sign::Zero, false) if strict => return Ok(false),
                _ => {}
            }
        }
        Ok(true)
    };
    let mut sum: Option<Vector3> = None;
    for (i, (a, _)) in directions.iter().enumerate() {
        for (b, _) in &directions[i + 1..] {
            let ray = a.cross(b);
            if sign(&ray.norm_squared(), dec)? == Sign::Zero {
                continue;
            }
            for ray in [ray.clone(), ray.scaled(&Real::from_i64(-1))] {
                if holds(&ray, false, dec)? {
                    sum = Some(match sum {
                        Some(sum) => sum.plus(&ray),
                        None => ray,
                    });
                }
            }
        }
    }
    match sum {
        Some(n) => holds(&n, true, dec),
        None => Ok(false),
    }
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
