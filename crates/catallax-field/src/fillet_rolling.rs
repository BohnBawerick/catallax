//! Fillets the elementary route does not section, compiled as rolling-ball
//! path pieces the field derives itself (4.2a steps 2 and 3).
//!
//! A plane with a cylinder in any pose, a cone or a sphere becomes a
//! [`Revolved`] piece; two cylinders with parallel axes become [`Generator`]
//! pieces, one per sharp run; two cylinders, cones or spheres about one axis
//! become a [`Coaxial`] piece; a cylinder with a second cylinder or a sphere
//! off its axis becomes a [`Saddle`] piece per ring or loop. The elementary route keeps every pair and edge it
//! sections today, so its answers are unchanged; these pieces only take what
//! it skips or refuses. The ball side of each support comes from the same
//! exact material quadrant the elementary route uses. A piece's region must be
//! clear of every other support, proven on its box; a support that reaches it
//! needs trimmed continuation and is refused as before.

use catallax_features::fillet::FilletDefinition;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{EdgeName, Origin, Role, StepId};
use catallax_shapes::{Curve, CurveId, GeomStore, Line, Point3, Surface, Vector3};

use super::{bad, coordinates, square_root, unit};
use crate::rolling_region::coaxial::{Coaxial, Profile};
use crate::rolling_region::elementary::Flat;
use crate::rolling_region::generator::{Generator, Wall};
use crate::rolling_region::piece::{Kind, Piece};
use crate::rolling_region::revolved::{Revolution, Revolved};
use crate::rolling_region::saddle::{Branch, Partner, Saddle, Shape};
use crate::{Plan, geometry};

/// One sharp edge run compiled as a rolling-ball path piece.
pub(super) struct Run {
    pub(super) piece: Piece,
    faces: [usize; 2],
    // Faces proven not to overflow the run: end faces of a straight run,
    // which bound it, a rim's mouths, narrower rounds that cut a run's
    // extension past a cliff, and the faces an earlier round rounded where
    // a sector over it has a fan of one material.
    ends: Vec<usize>,
    // Half-spaces `n . (x - p) >= 0` that hold the region, each bounded by
    // a junction with a tangent neighbour piece. A face on the far side of
    // one cannot enter the region.
    junctions: Vec<(Point3, Vector3)>,
    // A solid cylinder that holds the region: a ring whose ball lies inside
    // its cylinder support (4.2a step 5).
    inside: Option<Solid>,
}
/// The closed solid cylinder of squared radius `radius_sq` about the axis
/// through `point` along `axis`.
#[derive(Clone, Debug)]
pub(super) struct Solid {
    pub(super) point: Point3,
    pub(super) axis: Vector3,
    pub(super) radius_sq: Real,
}
impl Run {
    pub(super) fn new(piece: Piece, faces: [usize; 2], ends: Vec<usize>) -> Self {
        Self {
            piece,
            faces,
            ends,
            junctions: Vec::new(),
            inside: None,
        }
    }
    /// The run with the half-spaces of its junctions.
    pub(super) fn with_junctions(mut self, junctions: Vec<(Point3, Vector3)>) -> Self {
        self.junctions = junctions;
        self
    }
    /// The run whose region lies in the closed solid cylinder `solid`.
    pub(super) fn inside(mut self, solid: Solid) -> Self {
        self.inside = Some(solid);
        self
    }
    /// The run with the faces that cut its extension past a cliff (RECIPE
    /// 4.5.1, 4.2a step 5): narrower rounds it does not roll onto, each past
    /// the cliff. By the cliff rule they cut the region there, so the
    /// clearance leaves them to that cut, as it leaves end faces.
    pub(super) fn cut_by(mut self, faces: Vec<usize>) -> Self {
        self.ends.extend(faces);
        self
    }
    /// The run with faces that cross its fan where the fan is proven all
    /// material (convex) or all air (concave), so they bound nothing there
    /// (4.2a step 5: an earlier round's two faces under a sector whose ball
    /// lies outside that round's cylinder).
    pub(super) fn buried(mut self, faces: Vec<usize>) -> Self {
        self.ends.extend(faces);
        self
    }
}

/// Whether the rolling route takes this support pair: a plane with a cylinder,
/// cone or sphere, or two cylinders with parallel axes (step 2); and, for
/// named edges, two cylinders on crossing axes, a cylinder with a sphere, or
/// two cylinders, cones or spheres about one axis (step 3). An all-edges
/// selection still skips the step 3 pairs, as it did before them, so its
/// answers do not change.
pub(super) fn admits(
    supports: &[Surface; 2],
    named: bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    Ok(match (&supports[0], &supports[1]) {
        (Surface::Plane(_), Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_))
        | (Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_), Surface::Plane(_)) => true,
        (Surface::Cylinder(a), Surface::Cylinder(b))
            if geometry::sign(&a.axis.cross(&b.axis).norm_squared(), dec)? == Sign::Zero =>
        {
            true
        }
        _ if !named => false,
        // A plane with a torus (step 6), named edges only, as step 3's pairs.
        (Surface::Plane(_), Surface::Torus(_)) | (Surface::Torus(_), Surface::Plane(_)) => true,
        (Surface::Cylinder(_), Surface::Cylinder(_) | Surface::Sphere(_))
        | (Surface::Sphere(_), Surface::Cylinder(_)) => true,
        (a, b) => coaxial(a, b, dec)?,
    })
}

/// The axis of a cylinder or cone, as a point and a direction.
fn axis_of(surface: &Surface) -> Option<(&Point3, &Vector3)> {
    match surface {
        Surface::Cylinder(c) => Some((&c.origin, &c.axis)),
        Surface::Cone(c) => Some((&c.apex, &c.axis)),
        _ => None,
    }
}

/// Whether two curved surfaces of revolution share an axis: a cylinder or
/// cone has its axis, a sphere any line through its centre.
fn coaxial(a: &Surface, b: &Surface, dec: &mut Decider) -> Result<bool, Refusal> {
    let on = |p: &Point3, (q, d): (&Point3, &Vector3), dec: &mut Decider| {
        Ok::<bool, Refusal>(geometry::sign(&p.minus(q).cross(d).norm_squared(), dec)? == Sign::Zero)
    };
    Ok(match (a, b) {
        (Surface::Sphere(s), other) | (other, Surface::Sphere(s)) => match other {
            Surface::Sphere(t) => {
                geometry::sign(&s.center.minus(&t.center).norm_squared(), dec)? != Sign::Zero
            }
            _ => match axis_of(other) {
                Some(line) => on(&s.center, line, dec)?,
                None => false,
            },
        },
        _ => match (axis_of(a), axis_of(b)) {
            (Some(first), Some((q, d))) => {
                geometry::sign(&first.1.cross(d).norm_squared(), dec)? == Sign::Zero
                    && on(q, first, dec)?
            }
            _ => false,
        },
    })
}

fn reach_refusal(detail: &str) -> Refusal {
    Refusal::new("field.rolling-path-reach", detail)
}

pub(super) fn oriented(value: &Real, detail: &str, dec: &mut Decider) -> Result<bool, Refusal> {
    match geometry::sign(value, dec)? {
        Sign::Positive => Ok(true),
        Sign::Negative => Ok(false),
        Sign::Zero => Err(bad(detail)),
    }
}

/// The component of `v` perpendicular to `axis`.
pub(super) fn radial(v: &Vector3, axis: &Vector3) -> Vector3 {
    v.minus(&axis.scaled(&(&v.dot(axis) / &axis.norm_squared())))
}

/// A cylinder, cone, sphere or torus as a support of revolution, with the ball
/// direction's component along its outward normal at `sample`: positive
/// when the ball lies outside it. `None` for any other surface.
pub(super) fn revolution(
    other: &Surface,
    sample: &Point3,
    ball: &Vector3,
) -> Option<(Revolution, Real)> {
    let at = |origin: &Point3| sample.minus(origin);
    Some(match other {
        Surface::Cylinder(c) => (
            Revolution::Cylinder {
                point: c.origin.clone(),
                axis: c.axis.clone(),
                radius: square_root(&c.radius_sq),
            },
            ball.dot(&radial(&at(&c.origin), &c.axis)),
        ),
        Surface::Cone(c) => {
            let w = radial(&at(&c.apex), &c.axis);
            let along = unit(&c.axis).scaled(&(&c.slope * &w.norm_squared().sqrt()));
            (
                Revolution::Cone {
                    apex: c.apex.clone(),
                    axis: c.axis.clone(),
                    slope: c.slope.clone(),
                },
                ball.dot(&w.minus(&along)),
            )
        }
        Surface::Sphere(s) => (
            Revolution::Sphere {
                centre: s.center.clone(),
                radius: s.radius.clone(),
            },
            ball.dot(&at(&s.center)),
        ),
        // The outward normal points from the core circle's point in the
        // sample's meridian (4.2a step 6).
        Surface::Torus(t) => {
            let w = radial(&at(&t.center), &t.axis);
            let core = t.center.plus(&unit(&w).scaled(&t.major));
            (
                Revolution::Torus {
                    centre: t.center.clone(),
                    axis: t.axis.clone(),
                    major: t.major.clone(),
                    minor: t.minor.clone(),
                    edge: sample.clone(),
                },
                ball.dot(&sample.minus(&core)),
            )
        }
        _ => return None,
    })
}

/// A point on a closed edge, the end of its first axis, when it has no
/// rational point.
pub(super) fn algebraic_sample(edge: &Curve) -> Option<Point3> {
    match edge {
        Curve::Circle(c) => Some(
            c.center
                .plus(&unit(&c.ref_dir).scaled(&square_root(&c.radius_sq))),
        ),
        Curve::Ellipse(e) => Some(
            e.center
                .plus(&unit(&e.major_dir).scaled(&square_root(&e.major_sq))),
        ),
        _ => None,
    }
}

/// A rational point on a closed edge: an end of either axis.
pub(super) fn loop_sample(edge: &Curve) -> Option<Point3> {
    let (centre, axes) = match edge {
        Curve::Circle(c) => {
            let r = square_root(&c.radius_sq);
            let second = c.normal.cross(&c.ref_dir);
            (
                &c.center,
                [unit(&c.ref_dir).scaled(&r), unit(&second).scaled(&r)],
            )
        }
        Curve::Ellipse(e) => {
            let second = e.normal.cross(&e.major_dir);
            (
                &e.center,
                [
                    unit(&e.major_dir).scaled(&square_root(&e.major_sq)),
                    unit(&second).scaled(&square_root(&e.minor_sq)),
                ],
            )
        }
        _ => return None,
    };
    for axis in &axes {
        for sense in [1, -1] {
            let point = centre.plus(&axis.scaled(&Real::from_i64(sense)));
            if coordinates(&point).is_ok() {
                return Some(point);
            }
        }
    }
    None
}

impl Plan {
    pub(super) fn blend_origin(&self, faces: [usize; 2], step: &StepId) -> Result<Origin, Refusal> {
        let edge = EdgeName::new(
            self.faces[faces[0]].origin.clone(),
            self.faces[faces[1]].origin.clone(),
        )
        .map_err(|e| bad(&e.to_string()))?;
        Ok(Origin {
            step: step.clone(),
            role: Role::Blend { edge },
        })
    }

    /// The ball's direction at an edge point: the quadrant of `material`,
    /// rebuilt from the same gradients. For a convex edge the senses mark the
    /// supports whose gradient points away from the ball, for a concave one
    /// those whose gradient points toward it.
    pub(super) fn ball_direction(
        &self,
        faces: [usize; 2],
        p: &Point3,
        convex: bool,
        senses: [bool; 2],
    ) -> Result<Vector3, Refusal> {
        // The same exact gradients as `material`, at a point that may be
        // algebraic.
        let [a, b] = faces.map(|i| self.gradient_at(i, p));
        let dot = a.dot(&b);
        let da = a.minus(&b.scaled(&(&dot / &b.norm_squared())));
        let db = b.minus(&a.scaled(&(&dot / &a.norm_squared())));
        let [sa, sb] = senses.map(|sense| Real::from_i64(if convex == sense { -1 } else { 1 }));
        Ok(da.scaled(&sa).plus(&db.scaled(&sb)))
    }

    pub(super) fn flat(
        surface: &Surface,
        ball: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<Flat>, Refusal> {
        let Surface::Plane(plane) = surface else {
            return Ok(None);
        };
        let toward = oriented(
            &ball.dot(&plane.normal),
            "the ball lies in a support plane",
            dec,
        )?;
        Ok(Some(Flat {
            point: plane.origin.clone(),
            normal: if toward {
                plane.normal.clone()
            } else {
                plane.normal.scaled(&Real::from_i64(-1))
            },
        }))
    }

    /// Pieces for one support-intersection component the elementary route does
    /// not section. `None` when no rolling piece exists for the component.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn rolling_component(
        &self,
        root: usize,
        faces: [usize; 2],
        supports: &[Surface; 2],
        edge: &Curve,
        curve: CurveId,
        geom: &mut GeomStore,
        source_faces: &[usize],
        bound: &Rational,
        radius: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Run>>, Refusal> {
        match edge {
            Curve::Circle(_) | Curve::Ellipse(_)
                if !matches!(supports, [Surface::Cylinder(_), Surface::Cylinder(_)]) =>
            {
                let sample = loop_sample(edge)
                    .or_else(|| algebraic_sample(edge))
                    .ok_or_else(|| bad("a closed edge needs a point for its material"))?;
                self.revolved_runs(
                    root,
                    faces,
                    supports,
                    Some(edge),
                    &sample,
                    radius,
                    step,
                    dec,
                )
            }
            // A plane normal to a torus's axis meets it in circles about the
            // axis; an irrational one is kept as an intersection curve, whose
            // witness is a point of the edge (4.2a step 6).
            Curve::Intersection(intersection)
                if matches!(
                    supports,
                    [Surface::Plane(_), Surface::Torus(_)] | [Surface::Torus(_), Surface::Plane(_)]
                ) =>
            {
                let sample = geom.position(intersection.witness).clone();
                self.revolved_runs(root, faces, supports, None, &sample, radius, step, dec)
            }
            Curve::Ellipse(ellipse)
                if matches!(supports, [Surface::Cylinder(_), Surface::Cylinder(_)]) =>
            {
                self.crossing_run(root, faces, supports, ellipse, radius, step, dec)
            }
            Curve::Intersection(intersection)
                if matches!(
                    supports,
                    [
                        Surface::Cylinder(_),
                        Surface::Cylinder(_) | Surface::Sphere(_)
                    ] | [Surface::Sphere(_), Surface::Cylinder(_)]
                ) =>
            {
                let sample = geom.position(intersection.witness).clone();
                let Some((convex, senses)) = self.material_point(root, &sample, faces, dec)? else {
                    return Ok(Some(vec![]));
                };
                let ball = self.ball_direction(faces, &sample, convex, senses)?;
                Ok(Some(vec![self.saddle_run(
                    faces, supports, &sample, &ball, convex, step, radius, dec,
                )?]))
            }
            Curve::Line(line)
                if matches!(supports, [Surface::Cylinder(_), Surface::Cylinder(_)]) =>
            {
                let samples = self.line_samples(line, curve, geom, source_faces, bound, dec)?;
                let mut found: Vec<([Rational; 2], bool, [bool; 2], Point3)> = Vec::new();
                for (sample, interval) in samples {
                    let interval = interval.expect("a line sample has its stations");
                    let Some((convex, senses)) = self.material(root, &sample, faces, dec)? else {
                        continue;
                    };
                    if let Some(last) = found.last_mut()
                        && last.0[1] == interval[0]
                        && last.1 == convex
                        && last.2 == senses
                    {
                        last.0[1] = interval[1].clone();
                        continue;
                    }
                    found.push((interval, convex, senses, sample));
                }
                let mut runs = Vec::with_capacity(found.len());
                for (interval, convex, senses, sample) in found {
                    runs.push(self.generator_run(
                        root,
                        faces,
                        supports,
                        line,
                        &interval,
                        (convex, senses),
                        &sample,
                        source_faces,
                        radius,
                        step,
                        dec,
                    )?);
                }
                Ok(Some(runs))
            }
            _ => Ok(None),
        }
    }

    /// The ring piece of a closed edge where a plane meets a cylinder, cone,
    /// sphere or torus, or the coaxial piece of two curved supports about one
    /// axis, from an edge point `sample`. `edge` is the edge's circle or
    /// ellipse when the intersection gave one; a coaxial piece needs its
    /// circle. An empty list when the edge is no sharp input boundary.
    #[allow(clippy::too_many_arguments)]
    fn revolved_runs(
        &self,
        root: usize,
        faces: [usize; 2],
        supports: &[Surface; 2],
        edge: Option<&Curve>,
        sample: &Point3,
        radius: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Run>>, Refusal> {
        let Some((convex, senses)) = self.material_point(root, sample, faces, dec)? else {
            return Ok(Some(vec![]));
        };
        let ball = self.ball_direction(faces, sample, convex, senses)?;
        let [first, second] = supports;
        let (plane, other) = match (
            Self::flat(first, &ball, dec)?,
            Self::flat(second, &ball, dec)?,
        ) {
            (Some(plane), None) => (plane, second),
            (None, Some(plane)) => (plane, first),
            (None, None) => {
                let Some(Curve::Circle(circle)) = edge else {
                    return Ok(None);
                };
                let run = self.coaxial_run(
                    faces, supports, circle, sample, &ball, convex, step, radius, dec,
                )?;
                return Ok(Some(vec![run]));
            }
            _ => return Ok(None),
        };
        let Some((support, outside)) = revolution(other, sample, &ball) else {
            return Ok(None);
        };
        let outside = oriented(&outside, "the ball lies in a support's tangent plane", dec)?;
        let path = Revolved::new(plane, support, outside, radius.clone(), dec)?;
        let piece = Piece::new(
            Kind::Revolved(path),
            self.blend_origin(faces, step)?,
            convex,
            [first, second],
            dec,
        )?;
        Ok(Some(vec![Run::new(piece, faces, vec![])]))
    }

    /// A ring or loop where a cylinder meets a second cylinder or a sphere,
    /// as a [`Saddle`] piece. A ring is parametrized about the cylinder it
    /// winds round; a loop about the first cylinder it is one loop for.
    #[allow(clippy::too_many_arguments)]
    fn saddle_run(
        &self,
        faces: [usize; 2],
        supports: &[Surface; 2],
        sample: &Point3,
        ball: &Vector3,
        convex: bool,
        step: &StepId,
        radius: &Real,
        dec: &mut Decider,
    ) -> Result<Run, Refusal> {
        let side = |surface: &Surface, dec: &mut Decider| -> Result<bool, Refusal> {
            let toward = match surface {
                Surface::Cylinder(c) => ball.dot(&radial(&sample.minus(&c.origin), &c.axis)),
                Surface::Sphere(s) => ball.dot(&sample.minus(&s.center)),
                _ => return Err(bad("a saddle piece needs cylinders or a sphere")),
            };
            oriented(&toward, "the ball lies in a support's tangent plane", dec)
        };
        let partner = |surface: &Surface| match surface {
            Surface::Cylinder(c) => Partner::Cylinder {
                point: c.origin.clone(),
                axis: c.axis.clone(),
                radius: square_root(&c.radius_sq),
            },
            Surface::Sphere(s) => Partner::Sphere {
                centre: s.center.clone(),
                radius: s.radius.clone(),
            },
            _ => unreachable!("a saddle partner is a cylinder or a sphere"),
        };
        // Candidate parametrizations: each cylinder, with the other support
        // as its partner.
        let mut candidates = Vec::new();
        for (k, support) in supports.iter().enumerate() {
            if let Surface::Cylinder(c) = support {
                candidates.push((k, c, partner(&supports[1 - k])));
            }
        }
        let mut chosen = None;
        for wanted in [Shape::Rings, Shape::Loop] {
            for (k, c, other) in &candidates {
                let body = square_root(&c.radius_sq);
                if Saddle::edge_shape((&c.origin, &c.axis, &body), other, dec)? == wanted {
                    chosen = Some((*k, (*c).clone(), other.clone(), wanted));
                    break;
                }
            }
            if chosen.is_some() {
                break;
            }
        }
        let (k, cylinder, other, shape) = chosen.ok_or_else(|| {
            reach_refusal("the edge is neither a ring round a cylinder nor one loop")
        })?;
        let body = square_root(&cylinder.radius_sq);
        let branch = match shape {
            Shape::Rings => Branch::Ring(Saddle::edge_sign(
                (&cylinder.origin, &cylinder.axis, &body),
                &other,
                sample,
                dec,
            )?),
            _ => Branch::Loop,
        };
        let outside = [side(&supports[k], dec)?, side(&supports[1 - k], dec)?];
        let path = Saddle::new(
            (cylinder.origin.clone(), cylinder.axis.clone(), body),
            other,
            outside,
            branch,
            radius.clone(),
            dec,
        )?;
        let piece = Piece::new(
            Kind::Saddle(path),
            self.blend_origin(faces, step)?,
            convex,
            [&supports[0], &supports[1]],
            dec,
        )?;
        Ok(Run::new(piece, faces, vec![]))
    }

    /// Two equal cylinders whose axes meet cross in two ellipses through two
    /// poles; one [`Saddle`] crossing piece holds all four bands. It is built
    /// at the ellipse in the plane normal to the sum of the unit axes; the
    /// other ellipse adds nothing. The material is probed on all four bands,
    /// halfway between the poles, and must agree.
    #[allow(clippy::too_many_arguments)]
    fn crossing_run(
        &self,
        root: usize,
        faces: [usize; 2],
        supports: &[Surface; 2],
        ellipse: &catallax_shapes::Ellipse,
        radius: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Run>>, Refusal> {
        let [Surface::Cylinder(a), Surface::Cylinder(b)] = supports else {
            unreachable!("a crossing of two cylinders")
        };
        let (ua, ub) = (unit(&a.axis), unit(&b.axis));
        if geometry::sign(&ua.plus(&ub).cross(&ellipse.normal).norm_squared(), dec)? != Sign::Zero {
            return Ok(Some(vec![]));
        }
        let normal = ua.cross(&ub);
        let t = &b.origin.minus(&a.origin).cross(&ub).dot(&normal) / &normal.norm_squared();
        let meeting = a.origin.plus(&ua.scaled(&t));
        let across = unit(&normal.cross(&ua));
        let body = square_root(&a.radius_sq);
        let mut found: Option<((bool, [bool; 2]), Point3)> = None;
        let mut missing = 0;
        for nu in [ua.plus(&ub), ua.minus(&ub)] {
            for sense in [1, -1] {
                let m = across.scaled(&Real::from_i64(sense));
                let rise = -&(&(&body * &m.dot(&nu)) / &ua.dot(&nu));
                let sample = meeting.plus(&m.scaled(&body)).plus(&ua.scaled(&rise));
                match self.material_point(root, &sample, faces, dec)? {
                    None => missing += 1,
                    Some(material) => match &found {
                        None => found = Some((material, sample)),
                        Some((first, _)) if *first == material => {}
                        Some(_) => {
                            return Err(bad(
                                "the four bands of a crossing differ in their material side",
                            ));
                        }
                    },
                }
            }
        }
        let Some(((convex, senses), sample)) = found else {
            return Ok(Some(vec![]));
        };
        if missing > 0 {
            return Err(bad(
                "a crossing clipped to fewer than four bands needs ends at its poles",
            ));
        }
        let ball = self.ball_direction(faces, &sample, convex, senses)?;
        let mut outside = Vec::with_capacity(2);
        for c in [a, b] {
            outside.push(oriented(
                &ball.dot(&radial(&sample.minus(&c.origin), &c.axis)),
                "the ball lies in a support's tangent plane",
                dec,
            )?);
        }
        let path = Saddle::new(
            (meeting.clone(), a.axis.clone(), body),
            Partner::Cylinder {
                point: meeting,
                axis: b.axis.clone(),
                radius: square_root(&b.radius_sq),
            },
            [outside[0], outside[1]],
            Branch::Crossing,
            radius.clone(),
            dec,
        )?;
        let piece = Piece::new(
            Kind::Saddle(path),
            self.blend_origin(faces, step)?,
            convex,
            [&supports[0], &supports[1]],
            dec,
        )?;
        Ok(Some(vec![Run::new(piece, faces, vec![])]))
    }

    /// The circle where two curved supports about one axis meet, as a
    /// [`Coaxial`] piece. The profiles are written in the circle's frame: its
    /// centre is the origin of stations and its normal the axis.
    #[allow(clippy::too_many_arguments)]
    fn coaxial_run(
        &self,
        faces: [usize; 2],
        supports: &[Surface; 2],
        circle: &catallax_shapes::Circle,
        sample: &Point3,
        ball: &Vector3,
        convex: bool,
        step: &StepId,
        radius: &Real,
        dec: &mut Decider,
    ) -> Result<Run, Refusal> {
        let axis = unit(&circle.normal);
        let point = &circle.center;
        let station = |p: &Point3| p.minus(point).dot(&axis);
        let mut profiles = Vec::with_capacity(2);
        let mut sides = Vec::with_capacity(2);
        for support in supports {
            let (profile, outside) = match support {
                Surface::Cylinder(c) => (
                    Profile::Cylinder {
                        radius: square_root(&c.radius_sq),
                    },
                    ball.dot(&radial(&sample.minus(&c.origin), &c.axis)),
                ),
                Surface::Cone(c) => {
                    let w = radial(&sample.minus(&c.apex), &c.axis);
                    let along = unit(&c.axis).scaled(&(&c.slope * &w.norm_squared().sqrt()));
                    let sense = if oriented(&c.axis.dot(&axis), "a cone axis across the edge", dec)?
                    {
                        Real::one()
                    } else {
                        Real::from_i64(-1)
                    };
                    (
                        Profile::Cone {
                            apex: station(&c.apex),
                            slope: &sense * &c.slope,
                        },
                        ball.dot(&w.minus(&along)),
                    )
                }
                Surface::Sphere(s) => (
                    Profile::Sphere {
                        centre: station(&s.center),
                        radius: s.radius.clone(),
                    },
                    ball.dot(&sample.minus(&s.center)),
                ),
                _ => return Err(bad("a coaxial piece needs cylinders, cones or spheres")),
            };
            profiles.push(profile);
            sides.push(oriented(
                &outside,
                "the ball lies in a support's tangent plane",
                dec,
            )?);
        }
        let profiles: [Profile; 2] = profiles.try_into().map_err(|_| bad("two supports"))?;
        let path = Coaxial::new(
            point.clone(),
            axis,
            profiles,
            [sides[0], sides[1]],
            [square_root(&circle.radius_sq), Real::zero()],
            radius.clone(),
            dec,
        )?;
        let piece = Piece::new(
            Kind::Coaxial(path),
            self.blend_origin(faces, step)?,
            convex,
            [&supports[0], &supports[1]],
            dec,
        )?;
        Ok(Run::new(piece, faces, vec![]))
    }

    /// The two supports of a straight piece, each with the ball's side.
    pub(super) fn generator_walls(
        supports: &[Surface; 2],
        ball: &Vector3,
        sample: &Point3,
        dec: &mut Decider,
    ) -> Result<[Wall; 2], Refusal> {
        let mut walls = Vec::with_capacity(2);
        for support in supports {
            walls.push(match support {
                Surface::Plane(_) => Wall::Plane(Self::flat(support, ball, dec)?.expect("a plane")),
                Surface::Cylinder(c) => Wall::Cylinder {
                    point: c.origin.clone(),
                    axis: c.axis.clone(),
                    radius: square_root(&c.radius_sq),
                    outside: oriented(
                        &ball.dot(&radial(&sample.minus(&c.origin), &c.axis)),
                        "the ball lies in a support's tangent plane",
                        dec,
                    )?,
                },
                _ => return Err(bad("a generator piece needs planes or cylinders")),
            });
        }
        walls.try_into().map_err(|_| bad("two supports"))
    }

    /// A straight run between supports parallel to its line, over the
    /// stations `interval` of `line`, bounded by an end plane at each end.
    /// An end plane normal to the line ends the piece there. A plane face
    /// that leans across the line instead ends it at a normal plane past the
    /// lean; see [`Plan::oblique_stations`].
    #[allow(clippy::too_many_arguments)]
    pub(super) fn generator_run(
        &self,
        root: usize,
        faces: [usize; 2],
        supports: &[Surface; 2],
        line: &Line,
        interval: &[Rational; 2],
        (convex, senses): (bool, [bool; 2]),
        sample: &Point3,
        source_faces: &[usize],
        radius: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Run, Refusal> {
        let ball = self.ball_direction(faces, sample, convex, senses)?;
        let walls = Self::generator_walls(supports, &ball, sample, dec)?;
        let mut ends = Vec::new();
        let mut found = [false; 2];
        // End planes per end, each flagged when it leans across the line.
        let mut planes: [Vec<(usize, bool)>; 2] = [Vec::new(), Vec::new()];
        for &i in source_faces {
            if faces.contains(&i) {
                continue;
            }
            let Surface::Plane(p) = &self.faces[i].surface else {
                continue;
            };
            let across = line.dir.dot(&p.normal);
            if geometry::sign(&across, dec)? == Sign::Zero {
                continue;
            }
            let station = &p.origin.minus(&line.origin).dot(&p.normal) / &across;
            for (end, t) in interval.iter().enumerate() {
                if geometry::sign(&(&station - &Real::from(t.clone())), dec)? != Sign::Zero {
                    continue;
                }
                let leans =
                    geometry::sign(&line.dir.cross(&p.normal).norm_squared(), dec)? != Sign::Zero;
                planes[end].push((i, leans));
                found[end] = true;
                ends.push(i);
            }
        }
        let oblique = planes.iter().flatten().any(|(_, leans)| *leans);
        if !oblique && found != [true, true] {
            return Err(bad("finite line strips need proven transversal end planes"));
        }
        let stations = if oblique {
            if found != [true, true] {
                return Err(bad("a generator piece needs end planes normal to its line"));
            }
            self.oblique_stations(root, &walls, line, interval, &planes, sample, radius, dec)?
        } else {
            interval.clone().map(Real::from)
        };
        let path = Generator::new(
            walls,
            line.origin.clone(),
            line.dir.clone(),
            stations,
            radius.clone(),
            dec,
        )?;
        let piece = Piece::new(
            Kind::Generator(path),
            self.blend_origin(faces, step)?,
            convex,
            [&supports[0], &supports[1]],
            dec,
        )?;
        Ok(Run::new(piece, faces, ends))
    }

    /// The stations a straight piece runs to when an end face leans across
    /// its line (4.2a step 4, the corpus row `vertex-4-one`), as the
    /// elementary strip's domain widens for one.
    ///
    /// Both supports must be planes. Each section of the region then lies in
    /// the triangle of its edge point and two feet: the region is the
    /// quadrilateral of the edge point, the feet and the centre, less the
    /// ball, and the ball holds the triangle of the centre and the feet. Past
    /// the station where one end plane has passed all three points, the
    /// section lies beyond that plane, so each end takes the plane that passes
    /// them first. The body must lie on the run's side of every leaning end
    /// plane: it is an intersection of certified convex planar leaves whose
    /// planes include the end plane. Then no material lies past an end plane,
    /// and the region removes the strip as the end faces cut it, as the B-rep
    /// cuts its exact cylinder by the unselected end walls. An irrational
    /// station is rounded outward to a dyadic rational, which adds only air.
    #[allow(clippy::too_many_arguments)]
    fn oblique_stations(
        &self,
        root: usize,
        walls: &[Wall; 2],
        line: &Line,
        interval: &[Rational; 2],
        planes: &[Vec<(usize, bool)>; 2],
        sample: &Point3,
        radius: &Real,
        dec: &mut Decider,
    ) -> Result<[Real; 2], Refusal> {
        if !matches!(walls, [Wall::Plane(_), Wall::Plane(_)]) {
            return Err(bad(
                "a straight piece with a leaning end plane needs two plane supports",
            ));
        }
        let body = match self.convex_planar_input(root, dec) {
            Ok(body) => body,
            Err(refusal) if refusal.code.starts_with("numbers.") => return Err(refusal),
            Err(_) => return Err(bad("a leaning end plane needs a convex planar body")),
        };
        let plane_of = |i: usize| match &self.faces[i].surface {
            Surface::Plane(p) => p,
            _ => unreachable!("end planes are planes"),
        };
        for &(i, leans) in planes.iter().flatten() {
            if !leans {
                continue;
            }
            let p = plane_of(i);
            let mut bounds = false;
            for q in &body {
                if geometry::sign(&p.normal.cross(&q.normal).norm_squared(), dec)? == Sign::Zero
                    && geometry::sign(&p.origin.minus(&q.origin).dot(&q.normal), dec)? == Sign::Zero
                    && geometry::sign(&sample.minus(&q.origin).dot(&q.normal), dec)? != Sign::Zero
                {
                    bounds = true;
                    break;
                }
            }
            if !bounds {
                return Err(bad("a leaning end plane must bound the whole body"));
            }
        }
        let section = Generator::new(
            walls.clone(),
            line.origin.clone(),
            line.dir.clone(),
            interval.clone().map(Real::from),
            radius.clone(),
            dec,
        )?
        .section_at(&Real::zero());
        let mut stations = interval.clone().map(Real::from);
        for (end, list) in planes.iter().enumerate() {
            // Moving out of the run: down from the first end, up past the
            // second.
            let out = if end == 0 {
                Sign::Negative
            } else {
                Sign::Positive
            };
            let mut nearest: Option<Real> = None;
            for &(i, _) in list {
                let p = plane_of(i);
                let across = line.dir.dot(&p.normal);
                // The station where the plane has passed all three points.
                let mut passed: Option<Real> = None;
                for x in &section {
                    let t = &p.origin.minus(x).dot(&p.normal) / &across;
                    let further = match &passed {
                        None => true,
                        Some(s) => geometry::sign(&(&t - s), dec)? == out,
                    };
                    if further {
                        passed = Some(t);
                    }
                }
                let passed = passed.expect("three section points");
                let nearer = match &nearest {
                    None => true,
                    Some(s) => geometry::sign(&(s - &passed), dec)? == out,
                };
                if nearer {
                    nearest = Some(passed);
                }
            }
            let nearest = nearest.expect("an end plane at each end");
            stations[end] = match nearest.rational_ref() {
                Some(exact) => Real::from(exact.clone()),
                None => {
                    let ball = geometry::enclose(&nearest, 64, dec)?;
                    Real::from(if end == 0 { ball.lo() } else { ball.hi() })
                }
            };
        }
        Ok(stations)
    }

    /// A straight strip the elementary route refuses, as a generator run from
    /// its supports, ball sides, edge and stations. Its section, contacts and
    /// carrier are not read.
    pub(super) fn converted_run(
        &self,
        root: usize,
        definition: &FilletDefinition,
        interval: &Option<[Rational; 2]>,
        source_faces: &[usize],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Run, Refusal> {
        let Curve::Line(line) = &definition.edge else {
            return Err(bad("only straight runs become generator pieces"));
        };
        let interval = interval
            .as_ref()
            .ok_or_else(|| bad("a line strip needs a finite input interval"))?;
        let mut faces = [0; 2];
        for (slot, support) in faces.iter_mut().zip(&definition.supports) {
            *slot = self
                .faces
                .iter()
                .position(|f| f.surface.same_definition(support))
                .ok_or_else(|| bad("a strip support is not a source face"))?;
        }
        let sample = line.origin.plus(&line.dir.scaled(&Real::from(
            (&interval[0] + &interval[1]) / Rational::from_i64(2),
        )));
        self.generator_run(
            root,
            faces,
            &definition.supports,
            line,
            interval,
            (definition.convex, definition.same_sense),
            &sample,
            source_faces,
            &definition.radius,
            step,
            dec,
        )
    }

    /// Every other support is clear of the piece's box, so of its region.
    /// Supports equal to the piece's own or to its end planes are skipped.
    pub(super) fn rolling_clear(
        &self,
        run: &Run,
        source_faces: &[usize],
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let overflow = || {
            Refusal::new(
                "fillet.unsupported-overflow",
                "the field contact region reaches another support and needs trimmed continuation",
            )
        };
        let bounds = run.piece.bounds().ok_or_else(overflow)?;
        let region = bounds
            .each_ref()
            .map(|[lo, hi]| Ball::from_rational(lo, 128).hull(&Ball::from_rational(hi, 128)));
        let own: Vec<usize> = run.faces.iter().chain(&run.ends).copied().collect();
        // Boxes of the panels of a closed path, made when the whole box is
        // too coarse: a box round a ring holds its axis.
        let mut boxes: Option<Vec<[Ball; 3]>> = None;
        for &i in source_faces {
            let face = &self.faces[i];
            if own
                .iter()
                .any(|&j| j == i || face.support.same_definition(&self.faces[j].support))
            {
                continue;
            }
            dec.charge_work(1)?;
            if !face.polynomial.eval_ball(&region, dec)?.contains_zero() {
                continue;
            }
            if boxes.is_none() {
                boxes = Some(run.piece.boxes(dec)?);
            }
            for panel in boxes.as_deref().unwrap_or_default() {
                dec.charge_work(1)?;
                if face.polynomial.eval_ball(panel, dec)?.contains_zero() {
                    if self.beyond_junction(run, i, dec)?
                        || self.outside_solid(run, i, dec)?
                        || self.beside_prism(run, i, dec)?
                        || self.beside_box(run, i, dec)?
                        || self.on_own_plane(run, i, dec)?
                    {
                        break;
                    }
                    return Err(overflow());
                }
            }
        }
        Ok(())
    }

    /// Whether face `i` is a plane that keeps the run's solid cylinder on
    /// one side: parallel to its axis and at least its radius from it. The
    /// region lies in that closed solid (a ring's wall `rho^2 <= R^2` when
    /// its ball lies inside its cylinder support), so the plane meets it at
    /// most on its boundary, along the line where a tangent plane touches
    /// the cylinder: the face of an earlier round's tangent neighbour, which
    /// passes through the box of the ring that rolls over that round (4.2a
    /// step 5). A plane that supports a convex set meets it only on its
    /// boundary.
    fn outside_solid(&self, run: &Run, i: usize, dec: &mut Decider) -> Result<bool, Refusal> {
        let (Some(solid), Surface::Plane(plane)) = (&run.inside, &self.faces[i].surface) else {
            return Ok(false);
        };
        dec.charge_work(1)?;
        if geometry::sign(&plane.normal.dot(&solid.axis), dec)? != Sign::Zero {
            return Ok(false);
        }
        let offset = plane.normal.dot(&solid.point.minus(&plane.origin));
        let reach = &solid.radius_sq * &plane.normal.norm_squared();
        Ok(geometry::sign(&(&offset.square() - &reach), dec)? != Sign::Negative)
    }

    /// Whether face `i` is a plane with a straight piece's prism on one
    /// closed side: every vertex of [`Generator::plane_prism`] is on that
    /// side. A plane that supports a convex set holding the region meets the
    /// region at most on its boundary. The test is exact, so it holds where
    /// the piece's axis-aligned box is too loose: a wide strip in a rotated
    /// pose, whose box reaches the part's far faces (4.2a step 5).
    fn beside_prism(&self, run: &Run, i: usize, dec: &mut Decider) -> Result<bool, Refusal> {
        let (Kind::Generator(path), Surface::Plane(plane)) =
            (run.piece.kind(), &self.faces[i].surface)
        else {
            return Ok(false);
        };
        let Some(points) = path.plane_prism() else {
            return Ok(false);
        };
        let mut side = Sign::Zero;
        for p in &points {
            dec.charge_work(1)?;
            match geometry::sign(&plane.normal.dot(&p.minus(&plane.origin)), dec)? {
                Sign::Zero => {}
                found if side == Sign::Zero => side = found,
                found if found != side => return Ok(false),
                _ => {}
            }
        }
        Ok(true)
    }

    /// Whether face `i` is a plane that holds one of the run's own plane
    /// supports, a face of another leaf on the same plane, such as a boss's
    /// bottom face on the plate it stands on (4.2a step 5). The region lies
    /// on the ball's closed side of each support it rolls on, so such a
    /// plane meets it at most on its boundary.
    fn on_own_plane(&self, run: &Run, i: usize, dec: &mut Decider) -> Result<bool, Refusal> {
        let Surface::Plane(plane) = &self.faces[i].surface else {
            return Ok(false);
        };
        for &own in &run.faces {
            let Surface::Plane(support) = &self.faces[own].surface else {
                continue;
            };
            dec.charge_work(1)?;
            if geometry::sign(&plane.normal.cross(&support.normal).norm_squared(), dec)?
                == Sign::Zero
                && geometry::sign(&plane.normal.dot(&support.origin.minus(&plane.origin)), dec)?
                    == Sign::Zero
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Whether face `i`'s leaf puts the whole face beyond one side of the
    /// run's box, where the region cannot be: the face meets the region at
    /// most in that side's plane. A curved face evaluated over a long box
    /// can seem to reach it by interval dependency alone, such as a far
    /// corner round over a strip that runs past a cliff along the whole edge
    /// (4.2a step 5). Proven from the face's leaf, as at a junction.
    fn beside_box(&self, run: &Run, i: usize, dec: &mut Decider) -> Result<bool, Refusal> {
        let (geometry::Leaf::Feature(feature), Some(bounds)) =
            (&self.leaves[self.faces[i].leaf], run.piece.bounds())
        else {
            return Ok(false);
        };
        for (axis, [lo, hi]) in bounds.iter().enumerate() {
            for (sense, level) in [(-1, -lo), (1, hi.clone())] {
                let mut normal = [Real::zero(), Real::zero(), Real::zero()];
                normal[axis] = Real::from_i64(sense);
                let [x, y, z] = normal;
                let normal = Vector3::new(x, y, z);
                if feature.face_above(self.definition_role(i), &normal, &Real::from(level), dec)? {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// Whether face `i` lies on the far side of one of the run's junctions,
    /// `n . (x - p) <= 0`, while the region lies on the near side: then the
    /// face meets the region at most in the junction plane, which no face of
    /// the region lies on (4.2a step 4). Proven from the face's leaf.
    fn beyond_junction(&self, run: &Run, i: usize, dec: &mut Decider) -> Result<bool, Refusal> {
        let geometry::Leaf::Feature(feature) = &self.leaves[self.faces[i].leaf] else {
            return Ok(false);
        };
        for (point, normal) in &run.junctions {
            let away = normal.scaled(&Real::from_i64(-1));
            let level = away.dot(&point.minus(&Point3::origin()));
            if feature.face_above(self.definition_role(i), &away, &level, dec)? {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Compiles one run's region and applies it to `result`.
    pub(super) fn apply_run(
        &mut self,
        result: usize,
        run: Run,
        source_faces: &[usize],
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        self.rolling_clear(&run, source_faces, dec)?;
        let convex = run.piece.convex();
        // The face joins the registry under the piece's name when it lies on
        // an exact surface (a ring torus, a sphere or a cylinder).
        let faces = run
            .piece
            .face()
            .map(|surface| geometry::face(surface, run.piece.origin().clone(), dec))
            .transpose()?
            .into_iter()
            .collect();
        let node = self.leaf(geometry::Leaf::Rolling(Box::new(run.piece)), faces, dec)?;
        self.boolean(
            if convex {
                catallax_recipe::BooleanKind::Subtract
            } else {
                catallax_recipe::BooleanKind::Union
            },
            result,
            node,
            dec,
        )
    }
}
