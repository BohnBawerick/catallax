//! The exact torus a ball sweeps rolling round a circle where two surfaces of
//! revolution about one axis meet.
//!
//! A plane normal to the axis, and a cylinder, cone, sphere or torus about it,
//! cut a half-plane through the axis in a profile that is a line or a circle.
//! Write a point of that half-plane `(rho, z)`: its distance from the axis and
//! its height along it. A ball of radius `r` touches a support where its centre
//! lies on the profile's offset by `r` to the ball's side, so ball centres lie
//! where the two offsets cross. Turning a crossing about the axis sweeps the
//! circle of centres, and the ball an exact torus whose contacts are parallel
//! circles of the two supports. Nothing is fitted.
//!
//! The edge circle fixes the axis and the profile point `e` where the crossing
//! starts at radius zero. An offset line moves, and an offset circle grows or
//! shrinks, in proportion to the radius. Two offset lines cross once, at a
//! point affine in the radius. A line and a circle, or two circles, cross
//! twice while a few conditions affine in the radius hold, and each crossing
//! keeps its side of the line through the circles' centres, or its side along
//! the offset line. The crossing that starts at `e` is therefore the one on
//! `e`'s side at radius `r`: the radius homotopy never changes branch. Where a
//! condition fails the offsets stop crossing, and that branch has no ball of
//! the radius.
//!
//! A torus's tube angle starts on its outer equator and turns toward its axis.
//! The axis is chosen so that a band ending on that equator lies on the axis
//! side of it. A band across the equator crosses the tube's periodic seam.

use catallax_numbers::{Decider, Real, Refusal, Sign};

use crate::{BallSide, Circle, Point3, Surface, Torus, TorusSheet, Vector3, math};

fn unsupported(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.rolling-torus", detail)
}

// A point or direction of the half-plane through the axis, `[rho, z]`.
type Pair = [Real; 2];

fn dot(a: &Pair, b: &Pair) -> Real {
    &(&a[0] * &b[0]) + &(&a[1] * &b[1])
}
fn cross(a: &Pair, b: &Pair) -> Real {
    &(&a[0] * &b[1]) - &(&a[1] * &b[0])
}
fn minus(a: &Pair, b: &Pair) -> Pair {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
fn plus(a: &Pair, b: &Pair) -> Pair {
    [&a[0] + &b[0], &a[1] + &b[1]]
}
fn scaled(a: &Pair, by: &Real) -> Pair {
    [&a[0] * by, &a[1] * by]
}
// A quarter turn from `rho` toward `z`.
fn turned(a: &Pair) -> Pair {
    [-&a[1], a[0].clone()]
}
fn side(side: BallSide) -> Real {
    Real::from_i64(match side {
        BallSide::Outside => 1,
        BallSide::Inside => -1,
    })
}
fn positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(math::sign(value, dec)? == Sign::Positive)
}
fn zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(math::sign(value, dec)? == Sign::Zero)
}

// A support's profile in the half-plane through the axis.
enum Profile {
    // `normal . q = offset` for a unit `normal`, the support's normal.
    Line { normal: Pair, offset: Real },
    // `|q - centre| = radius`, the support's normal pointing away from the centre.
    Circle { centre: Pair, radius: Real },
}
impl Profile {
    fn normal_at(&self, q: &Pair) -> Pair {
        match self {
            Self::Line { normal, .. } => normal.clone(),
            Self::Circle { centre, radius } => scaled(&minus(q, centre), &(&Real::one() / radius)),
        }
    }
    // Zero exactly on the profile.
    fn residual(&self, q: &Pair) -> Real {
        match self {
            Self::Line { normal, offset } => &dot(normal, q) - offset,
            Self::Circle { centre, radius } => {
                let d = minus(q, centre);
                &dot(&d, &d) - &radius.square()
            }
        }
    }
    // The profile moved by `shift` along its normal.
    fn offset(&self, shift: &Real) -> Self {
        match self {
            Self::Line { normal, offset } => Self::Line {
                normal: normal.clone(),
                offset: offset + shift,
            },
            Self::Circle { centre, radius } => Self::Circle {
                centre: centre.clone(),
                radius: radius + shift,
            },
        }
    }
    // The contact of a ball centred on the offset by `shift`.
    fn foot(&self, centre: &Pair, shift: &Real) -> Pair {
        match self {
            Self::Line { normal, .. } => minus(centre, &scaled(normal, shift)),
            Self::Circle {
                centre: middle,
                radius,
            } => plus(
                middle,
                &scaled(&minus(centre, middle), &(radius / &(radius + shift))),
            ),
        }
    }
}

// The axis of the edge circle, with a unit direction.
struct Frame {
    origin: Point3,
    axis: Vector3,
}
impl Frame {
    fn height(&self, p: &Point3) -> Real {
        p.minus(&self.origin).dot(&self.axis)
    }
    fn parallel(&self, v: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
        zero(&v.cross(&self.axis).norm_squared(), dec)
    }
    fn on_axis(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        self.parallel(&p.minus(&self.origin), dec)
    }
    // The support's profile, when it turns about this axis.
    fn profile(&self, surface: &Surface, dec: &mut Decider) -> Result<Profile, Refusal> {
        let apart = || unsupported("the supports do not turn about the edge circle's axis");
        match surface {
            Surface::Plane(plane) => {
                if !self.parallel(&plane.normal, dec)? {
                    return Err(apart());
                }
                let up = if positive(&plane.normal.dot(&self.axis), dec)? {
                    Real::one()
                } else {
                    Real::from_i64(-1)
                };
                Ok(Profile::Line {
                    offset: &up * &self.height(&plane.origin),
                    normal: [Real::zero(), up],
                })
            }
            Surface::Cylinder(cylinder) => {
                if !self.parallel(&cylinder.axis, dec)? || !self.on_axis(&cylinder.origin, dec)? {
                    return Err(apart());
                }
                Ok(Profile::Line {
                    normal: [Real::one(), Real::zero()],
                    offset: math::exact_sqrt(&cylinder.radius_sq),
                })
            }
            Surface::Cone(cone) => {
                if !self.parallel(&cone.axis, dec)? || !self.on_axis(&cone.apex, dec)? {
                    return Err(apart());
                }
                // The generator from the apex runs `(slope, 1)` along the
                // cone's own axis, and the normal turns away from that axis.
                let along = if positive(&cone.axis.dot(&self.axis), dec)? {
                    cone.slope.clone()
                } else {
                    -&cone.slope
                };
                let length = math::exact_sqrt(&(Real::one() + cone.slope.square()));
                let normal = [&Real::one() / &length, -&(&along / &length)];
                Ok(Profile::Line {
                    offset: &normal[1] * &self.height(&cone.apex),
                    normal,
                })
            }
            Surface::Sphere(sphere) => {
                if !self.on_axis(&sphere.center, dec)? {
                    return Err(apart());
                }
                Ok(Profile::Circle {
                    centre: [Real::zero(), self.height(&sphere.center)],
                    radius: sphere.radius.clone(),
                })
            }
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
                if !self.parallel(axis, dec)? || !self.on_axis(center, dec)? {
                    return Err(apart());
                }
                Ok(Profile::Circle {
                    centre: [major.clone(), self.height(center)],
                    radius: minor.clone(),
                })
            }
            _ => Err(unsupported("a support is not a surface of revolution")),
        }
    }
}

// The axis, the edge's profile point and both profiles.
fn setup(
    supports: [&Surface; 2],
    edge: &Circle,
    dec: &mut Decider,
) -> Result<(Frame, Pair, [Profile; 2]), Refusal> {
    if !positive(&edge.radius_sq, dec)? || zero(&edge.normal.norm_squared(), dec)? {
        return Err(Refusal::new(
            "shapes.rolling-torus-edge",
            "the edge circle is degenerate",
        ));
    }
    let frame = Frame {
        origin: edge.center.clone(),
        axis: edge
            .normal
            .scaled(&(&Real::one() / &math::exact_sqrt(&edge.normal.norm_squared()))),
    };
    let e = [math::exact_sqrt(&edge.radius_sq), Real::zero()];
    let profiles = [
        frame.profile(supports[0], dec)?,
        frame.profile(supports[1], dec)?,
    ];
    for profile in &profiles {
        if !zero(&profile.residual(&e), dec)? {
            return Err(Refusal::new(
                "shapes.rolling-torus-edge",
                "the edge circle does not lie on both supports",
            ));
        }
    }
    Ok((frame, e, profiles))
}

/// `start + slope r > 0`, a condition on the radius `r`.
struct Affine {
    start: Real,
    slope: Real,
}

// The conditions under which the offsets by `signs` times the radius keep
// crossing twice. Two lines always cross once and need none.
fn conditions(profiles: &[Profile; 2], signs: &[Real; 2]) -> Vec<Affine> {
    let affine = |start: Real, slope: Real| Affine { start, slope };
    match profiles {
        [Profile::Line { .. }, Profile::Line { .. }] => Vec::new(),
        [
            Profile::Line { normal, offset },
            Profile::Circle { centre, radius },
        ]
        | [
            Profile::Circle { centre, radius },
            Profile::Line { normal, offset },
        ] => {
            let (line, circle) = if matches!(profiles[0], Profile::Line { .. }) {
                (&signs[0], &signs[1])
            } else {
                (&signs[1], &signs[0])
            };
            // The line lies `across` from the circle's centre along its normal.
            let across = offset - &dot(normal, centre);
            vec![
                affine(radius - &across, circle - line),
                affine(radius + &across, circle + line),
            ]
        }
        [
            Profile::Circle {
                centre: first,
                radius: a,
            },
            Profile::Circle {
                centre: second,
                radius: b,
            },
        ] => {
            let d = minus(second, first);
            let apart = math::exact_sqrt(&dot(&d, &d));
            let [s, t] = signs;
            vec![
                affine(&(a + b) - &apart, s + t),
                affine(&(&apart - a) + b, t - s),
                affine(&(&apart + a) - b, s - t),
            ]
        }
    }
}

// Where two offsets cross on the branch through `e`. Offsetting keeps each
// line's normal and each circle's centre, so the side of the line through the
// centres, or along the offset line, on which the unmoved profiles cross at
// `e` picks the same crossing at every radius.
fn crossing(offsets: &[Profile; 2], e: &Pair, dec: &mut Decider) -> Result<Pair, Refusal> {
    let tangent = || unsupported("the supports are tangent along the edge");
    let pick = |base: &Pair, across: &Pair, dec: &mut Decider| -> Result<Real, Refusal> {
        match math::sign(&dot(&minus(e, base), across), dec)? {
            Sign::Positive => Ok(Real::one()),
            Sign::Negative => Ok(Real::from_i64(-1)),
            Sign::Zero => Err(tangent()),
        }
    };
    match offsets {
        [
            Profile::Line {
                normal: a,
                offset: p,
            },
            Profile::Line {
                normal: b,
                offset: q,
            },
        ] => {
            let det = cross(a, b);
            if zero(&det, dec)? {
                return Err(tangent());
            }
            Ok([
                &(&(p * &b[1]) - &(q * &a[1])) / &det,
                &(&(&a[0] * q) - &(&b[0] * p)) / &det,
            ])
        }
        [
            Profile::Line { normal, offset },
            Profile::Circle { centre, radius },
        ]
        | [
            Profile::Circle { centre, radius },
            Profile::Line { normal, offset },
        ] => {
            let across = offset - &dot(normal, centre);
            let along = turned(normal);
            let way = pick(centre, &along, dec)?;
            let half = math::exact_sqrt(&(radius.square() - across.square()));
            Ok(plus(
                &plus(centre, &scaled(normal, &across)),
                &scaled(&along, &(&way * &half)),
            ))
        }
        [
            Profile::Circle {
                centre: first,
                radius: a,
            },
            Profile::Circle {
                centre: second,
                radius: b,
            },
        ] => {
            let d = minus(second, first);
            let squared = dot(&d, &d);
            let along = turned(&d);
            let way = pick(first, &along, dec)?;
            let part =
                &(&(&squared + &a.square()) - &b.square()) / &(&Real::from_i64(2) * &squared);
            let height = math::exact_sqrt(&(&(&a.square() / &squared) - &part.square()));
            Ok(plus(
                &plus(first, &scaled(&d, &part)),
                &scaled(&along, &(&way * &height)),
            ))
        }
    }
}

// Where two offsets that touch meet: the foot of the perpendicular from a
// circle's centre to a tangent line, or the point of two tangent circles on
// the line through their centres. The square root that separates two
// crossings is zero here and is left out. `None` for two lines, which cross.
fn touch(offsets: &[Profile; 2]) -> Option<Pair> {
    match offsets {
        [Profile::Line { .. }, Profile::Line { .. }] => None,
        [
            Profile::Line { normal, offset },
            Profile::Circle { centre, .. },
        ]
        | [
            Profile::Circle { centre, .. },
            Profile::Line { normal, offset },
        ] => {
            let across = offset - &dot(normal, centre);
            Some(plus(centre, &scaled(normal, &across)))
        }
        [
            Profile::Circle {
                centre: first,
                radius: a,
            },
            Profile::Circle {
                centre: second,
                radius: b,
            },
        ] => {
            let d = minus(second, first);
            let squared = dot(&d, &d);
            let part =
                &(&(&squared + &a.square()) - &b.square()) / &(&Real::from_i64(2) * &squared);
            Some(plus(first, &scaled(&d, &part)))
        }
    }
}

// Whether a direction lies on the shorter arc from `p` to `q`, ends included:
// a combination of the two with no negative weight.
fn on_arc(direction: &Pair, p: &Pair, q: &Pair, dec: &mut Decider) -> Result<bool, Refusal> {
    let turn = math::sign(&cross(p, q), dec)?;
    let from = math::sign(&cross(p, direction), dec)?;
    let to = math::sign(&cross(direction, q), dec)?;
    Ok([from, to].iter().all(|s| *s == turn || *s == Sign::Zero))
}

/// The torus a ball sweeps rolling round a circle where two surfaces of
/// revolution about the circle's axis meet, with its contact circles.
#[derive(Clone, Debug)]
pub struct RollingTorus {
    /// A ring torus, or the outer sheet of a horn or spindle torus when the
    /// circle of ball centres is no wider than the ball and the band between
    /// the contacts stays off the axis.
    pub surface: Surface,
    /// The contact circle on each support, with the edge circle's normal and
    /// seam direction, so each runs the edge's way round.
    pub contacts: [Circle; 2],
}

/// The half tube of a ball at the radius where two coaxial supports' offsets
/// touch: a band of the ring torus between two parallels that are opposite
/// ends of the tube's diameters.
#[derive(Clone, Debug)]
pub struct InscribedTorus {
    /// The whole ring torus. The band is the half through `middle`.
    pub torus: Torus,
    /// The contact circle on each support, with the edge circle's normal and
    /// seam direction.
    pub contacts: [Circle; 2],
    /// The parallel halfway along the band, on the edge's side of the
    /// contacts. The other half belongs to the crossing branch opposite.
    pub middle: Circle,
}

impl RollingTorus {
    /// The torus of a ball of `radius` on the given side of each support:
    /// `Outside` is the side its normal points to. Supports are planes normal
    /// to the edge circle's axis, or cylinders, cones, spheres, tori or torus
    /// sheets about it, and the edge lies on both. Refuses with
    /// `shapes.rolling-torus-edge` when the edge circle does not lie on both
    /// supports, and with `shapes.rolling-torus` for supports that do not turn
    /// about its axis, tangent supports, offsets that stop crossing, a ball
    /// centre or contact on the axis, or a band that would cross the axis.
    pub fn new(
        supports: [&Surface; 2],
        sides: [BallSide; 2],
        edge: &Circle,
        radius: &Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if !positive(radius, dec)? {
            return Err(unsupported("the radius is not positive"));
        }
        let (frame, e, profiles) = setup(supports, edge, dec)?;
        let normals = profiles.each_ref().map(|profile| profile.normal_at(&e));
        if zero(&cross(&normals[0], &normals[1]), dec)? {
            return Err(unsupported("the supports are tangent along the edge"));
        }
        let signs = sides.map(side);
        for condition in conditions(&profiles, &signs) {
            if !positive(&condition.start, dec)? {
                return Err(unsupported("the supports do not cross at the edge"));
            }
            if !positive(&(&condition.start + &(&condition.slope * radius)), dec)? {
                return Err(unsupported(
                    "the offsets of the supports stop crossing before this radius",
                ));
            }
        }
        let shifts = signs.each_ref().map(|s| s * radius);
        let offsets = [
            profiles[0].offset(&shifts[0]),
            profiles[1].offset(&shifts[1]),
        ];
        let centre = crossing(&offsets, &e, dec)?;
        if !positive(&centre[0], dec)? {
            return Err(unsupported("the ball centre reaches the axis"));
        }
        let feet = [
            profiles[0].foot(&centre, &shifts[0]),
            profiles[1].foot(&centre, &shifts[1]),
        ];
        for foot in &feet {
            if !positive(&foot[0], dec)? {
                return Err(unsupported("a contact reaches the axis"));
            }
        }
        let [p, q] = feet.each_ref().map(|foot| minus(foot, &centre));
        if zero(&cross(&p, &q), dec)? {
            return Err(unsupported(
                "the ball has no shorter arc between its contacts",
            ));
        }
        let inward = [Real::from_i64(-1), Real::zero()];
        let ring = positive(&(&centre[0] - radius), dec)?;
        if !ring && on_arc(&inward, &p, &q, dec)? {
            return Err(unsupported("the band would cross the axis"));
        }
        // A band ending on the outer equator keeps to the axis side of it.
        let flip = if zero(&p[1], dec)? && positive(&p[0], dec)? {
            !positive(&q[1], dec)?
        } else if zero(&q[1], dec)? && positive(&q[0], dec)? {
            !positive(&p[1], dec)?
        } else {
            false
        };
        let on_axis = |height: &Real| frame.origin.plus(&frame.axis.scaled(height));
        let circle = |at: &Pair| Circle {
            center: on_axis(&at[1]),
            normal: edge.normal.clone(),
            radius_sq: at[0].square(),
            ref_dir: edge.ref_dir.clone(),
        };
        let axis = if flip {
            edge.normal.scaled(&Real::from_i64(-1))
        } else {
            edge.normal.clone()
        };
        let surface = if ring {
            Surface::Torus(Torus {
                center: on_axis(&centre[1]),
                axis,
                major: centre[0].clone(),
                minor: radius.clone(),
                ref_dir: edge.ref_dir.clone(),
            })
        } else {
            Surface::TorusSheet(TorusSheet {
                center: on_axis(&centre[1]),
                axis,
                major: centre[0].clone(),
                minor: radius.clone(),
                ref_dir: edge.ref_dir.clone(),
            })
        };
        Ok(Self {
            surface,
            contacts: feet.each_ref().map(circle),
        })
    }

    /// The half tube a ball of exactly the limit `radius` sweeps where the two
    /// offsets stop crossing and touch instead. `None` unless one crossing
    /// condition is exactly zero at `radius` and every other one holds
    /// strictly, so the centres form one circle shared by both crossing
    /// branches. There the contacts are opposite ends of a diameter, no arc
    /// is shorter, and the radius homotopy keeps the half on the edge's side:
    /// its one-sided limit. `None` too when the tube would reach the axis or
    /// the edge's side is undecided. Refuses as [`Self::new`] does for an
    /// edge off its supports or supports that do not turn about its axis.
    pub fn inscribed(
        supports: [&Surface; 2],
        sides: [BallSide; 2],
        edge: &Circle,
        radius: &Real,
        dec: &mut Decider,
    ) -> Result<Option<InscribedTorus>, Refusal> {
        if !positive(radius, dec)? {
            return Ok(None);
        }
        let (frame, e, profiles) = setup(supports, edge, dec)?;
        let normals = profiles.each_ref().map(|profile| profile.normal_at(&e));
        if zero(&cross(&normals[0], &normals[1]), dec)? {
            return Ok(None);
        }
        let signs = sides.map(side);
        let mut touching = 0;
        for condition in conditions(&profiles, &signs) {
            if !positive(&condition.start, dec)? {
                return Ok(None);
            }
            let at = &condition.start + &(&condition.slope * radius);
            if zero(&at, dec)? {
                touching += 1;
            } else if !positive(&at, dec)? {
                return Ok(None);
            }
        }
        if touching != 1 {
            return Ok(None);
        }
        let shifts = signs.each_ref().map(|s| s * radius);
        let offsets = [
            profiles[0].offset(&shifts[0]),
            profiles[1].offset(&shifts[1]),
        ];
        let Some(centre) = touch(&offsets) else {
            return Ok(None);
        };
        // A torus whose tube stays off the axis: every centre, contact and
        // band point is farther out than zero.
        if !positive(&(&centre[0] - radius), dec)? {
            return Ok(None);
        }
        let feet = [
            profiles[0].foot(&centre, &shifts[0]),
            profiles[1].foot(&centre, &shifts[1]),
        ];
        let [p, q] = feet.each_ref().map(|foot| minus(foot, &centre));
        if !zero(&cross(&p, &q), dec)? || positive(&dot(&p, &q), dec)? {
            return Ok(None);
        }
        for (foot, profile) in feet.iter().zip(&profiles) {
            if !zero(&profile.residual(foot), dec)? {
                return Ok(None);
            }
        }
        // The half from `p` that turns toward the edge point.
        let turn = match math::sign(&cross(&p, &minus(&e, &centre)), dec)? {
            Sign::Positive => Real::one(),
            Sign::Negative => Real::from_i64(-1),
            Sign::Zero => return Ok(None),
        };
        let middle = plus(&centre, &scaled(&turned(&p), &turn));
        let on_axis = |height: &Real| frame.origin.plus(&frame.axis.scaled(height));
        let circle = |at: &Pair| Circle {
            center: on_axis(&at[1]),
            normal: edge.normal.clone(),
            radius_sq: at[0].square(),
            ref_dir: edge.ref_dir.clone(),
        };
        Ok(Some(InscribedTorus {
            torus: Torus {
                center: on_axis(&centre[1]),
                axis: edge.normal.clone(),
                major: centre[0].clone(),
                minor: radius.clone(),
                ref_dir: edge.ref_dir.clone(),
            },
            contacts: feet.each_ref().map(circle),
            middle: circle(&middle),
        }))
    }

    /// The radius past which no ball rests on both supports along the edge,
    /// on the given sides, when the offsets prove one. A cone's offset also
    /// rounds its apex, so its profile proves no bound.
    pub fn radius_limit(
        supports: [&Surface; 2],
        sides: [BallSide; 2],
        edge: &Circle,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        if supports.iter().any(|s| matches!(s, Surface::Cone(_))) {
            return Ok(None);
        }
        let (_, _, profiles) = setup(supports, edge, dec)?;
        let mut limit: Option<Real> = None;
        for condition in conditions(&profiles, &sides.map(side)) {
            if !positive(&condition.start, dec)? {
                return Ok(None);
            }
            if positive(&condition.slope, dec)? || zero(&condition.slope, dec)? {
                continue;
            }
            let at = &condition.start / &(-&condition.slope);
            limit = Some(match limit {
                Some(limit) if positive(&(&limit - &at), dec)? => at,
                Some(limit) => limit,
                None => at,
            });
        }
        Ok(limit)
    }
}
