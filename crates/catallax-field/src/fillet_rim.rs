//! Rim contact (4.2a step 4): a concave ring between a plane and a
//! perpendicular cylinder whose plane contact meets circular mouths.
//!
//! The elementary ring refuses when another support reaches its region. This
//! route runs only then. It proves, from the field's own supports and
//! leaves, that every other source face either misses the ring's fan or is
//! the side of a hole whose mouth lies in the plane, and builds one [`Rim`]
//! piece whose ball drops onto each mouth's edge where the plane contact
//! would fall inside the mouth. Nothing here reads a B-rep fillet module, a
//! fillet definition's section, contacts or carrier, or the B-rep's
//! `RollingEnvelopeDefinition`.
//!
//! The fan (`R <= rho <= C`, `0 <= y <= r` in the ring's frame) holds the
//! region. A face misses it when its support misses every panel of the fan,
//! tested in the ring's frame by the support's exact second-order form about
//! each panel's centre, or when the face's leaf does: a face lies in its
//! leaf. A mouth is a Round primitive (a cylinder or cone frustum) whose axis
//! is parallel to the ring's, lying wholly on the material side of the plane
//! with one cap in it: that cap's circle is the edge, and its side face can
//! meet the fan, which lies on the ball side, only along that circle. A
//! probe just below the plane, inside the mouth and the fan, shows whether
//! the cap is a hole in the plane's face (air below) or buried in material.
//! Everything else that could bound the plane's face inside the fan would
//! reach the fan and be refused, so the face within the fan is the fan's
//! annulus less the hole mouths.

use catallax_features::fillet::FilletDefinition;
use catallax_numbers::{Ball, Decider, Rational, Real, RealPoly, Refusal, Sign};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Point3, Surface, Vector3};

use super::rolling::{Run, algebraic_sample, loop_sample, oriented, radial};
use super::{square_root, unit};
use crate::metric::Kind as Shape;
use crate::rolling_region::piece::{Kind, Piece};
use crate::rolling_region::revolved::{basis_of, hull_points};
use crate::rolling_region::rim::{Mouth, Rim, world};
use crate::{Instruction, Location, Plan, geometry};

const BITS: u32 = 96;

fn overflow(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-overflow", detail)
}

/// A Round primitive's cap circle in the ring's plane.
struct Cap {
    centre: Point3,
    radius: Real,
    depth: Real,
    flare: Real,
}

/// The ring's frame for the reach tests.
struct Frame {
    origin: [Ball; 3],
    // First basis direction, axis, second basis direction.
    axes: [[Ball; 3]; 3],
    // The fan's panels in the frame.
    fan: Vec<[Ball; 3]>,
}

fn mid(b: &Ball) -> Ball {
    let m = b.mid();
    Ball::from_rational(&m, BITS)
}

/// Whether a support's zero set may meet one of the fan's panels. A support
/// of degree at most two is its exact second-order form about each panel's
/// centre `m`: `p(m + J h) = p(m) + (J^T grad p(m)).h + h^T (J^T H J) h / 2`
/// with `J` the frame's axes, so a small panel is not widened by the
/// support's terms in far coordinates. Higher degrees are evaluated on the
/// panel's world box.
fn reaches(poly: &RealPoly, frame: &Frame, dec: &mut Decider) -> Result<bool, Refusal> {
    let quadric = poly.total_degree(dec)?.is_none_or(|d| d <= 2);
    let gradient: [RealPoly; 3] = std::array::from_fn(|i| poly.partial(i));
    let zero = [0, 1, 2].map(|_| Ball::from_rational(&Rational::zero(), BITS));
    let hessian: Option<[[Ball; 3]; 3]> = if quadric {
        let mut h: [[Ball; 3]; 3] = std::array::from_fn(|_| zero.clone());
        for (i, row) in h.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = gradient[i].partial(j).eval_ball(&zero, dec)?;
            }
        }
        Some(h)
    } else {
        None
    };
    for panel in &frame.fan {
        dec.charge_work(1)?;
        let Some(h) = &hessian else {
            let boxed = world(&frame.origin, &frame.axes, panel);
            if poly.eval_ball(&boxed, dec)?.contains_zero() {
                return Ok(true);
            }
            continue;
        };
        let centre = panel.each_ref().map(mid);
        let step: [Ball; 3] = std::array::from_fn(|k| &panel[k] - &centre[k]);
        let at = world(&frame.origin, &frame.axes, &centre);
        let mut value = poly.eval_ball(&at, dec)?;
        let g: Vec<Ball> = gradient
            .iter()
            .map(|p| p.eval_ball(&at, dec))
            .collect::<Result<_, _>>()?;
        for (k, axis) in frame.axes.iter().enumerate() {
            let slope = (0..3).fold(zero[0].clone(), |sum, i| &sum + &(&g[i] * &axis[i]));
            value = &value + &(&slope * &step[k]);
        }
        let half = Ball::from_rational(&(Rational::one() / Rational::from_i64(2)), BITS);
        for k in 0..3 {
            for l in k..3 {
                let mut a = zero[0].clone();
                for (i, row) in h.iter().enumerate() {
                    for (j, cell) in row.iter().enumerate() {
                        a = &a + &(&(&frame.axes[k][i] * cell) * &frame.axes[l][j]);
                    }
                }
                let term = if k == l {
                    &(&a * &half) * &step[k].square()
                } else {
                    &(&a * &step[k]) * &step[l]
                };
                value = &value + &term;
            }
        }
        if value.contains_zero() {
            return Ok(true);
        }
    }
    Ok(false)
}

impl Plan {
    /// The location just past `p` along `d`, from the leaves' exact ray
    /// probes at a point that may be algebraic.
    fn probe(
        &self,
        root: usize,
        p: &Point3,
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let mut locations = Vec::with_capacity(root + 1);
        for instruction in &self.instructions[..=root] {
            let location = match *instruction {
                Instruction::Empty => Location::Outside,
                Instruction::Leaf(i) => crate::distance::ray(&self.leaves[i], p, d, dec)?,
                Instruction::Boolean { kind, a, b } => {
                    crate::combine(kind, locations[a], locations[b])
                }
            };
            locations.push(location);
        }
        Ok(locations[root])
    }

    /// The rim route for a concave plane-cylinder ring whose elementary
    /// region another support reaches. `None` when the ring is not of this
    /// kind, so the caller keeps its refusal.
    pub(super) fn rim_run(
        &self,
        root: usize,
        definition: &FilletDefinition,
        source_faces: &[usize],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<Run>, Refusal> {
        let Curve::Circle(_) = &definition.edge else {
            return Ok(None);
        };
        if definition.convex {
            return Ok(None);
        }
        let supports = &definition.supports;
        let (plane_at, cylinder_at) = match supports {
            [Surface::Plane(_), Surface::Cylinder(_)] => (0, 1),
            [Surface::Cylinder(_), Surface::Plane(_)] => (1, 0),
            _ => return Ok(None),
        };
        let Surface::Cylinder(cylinder) = &supports[cylinder_at] else {
            unreachable!("a cylinder support")
        };
        let mut faces = [0; 2];
        for (slot, support) in faces.iter_mut().zip(supports) {
            *slot = self
                .faces
                .iter()
                .position(|f| f.surface.same_definition(support))
                .ok_or_else(|| super::bad("a ring support is not a source face"))?;
        }
        let sample = loop_sample(&definition.edge)
            .or_else(|| algebraic_sample(&definition.edge))
            .ok_or_else(|| super::bad("a ring needs a point for its material"))?;
        let ball_side = self.ball_direction(faces, &sample, false, definition.same_sense)?;
        let Some(plane) = Self::flat(&supports[plane_at], &ball_side, dec)? else {
            return Ok(None);
        };
        let outside = oriented(
            &ball_side.dot(&radial(&sample.minus(&cylinder.origin), &cylinder.axis)),
            "the ball lies in the cylinder's tangent plane",
            dec,
        )?;
        if !outside
            || geometry::sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)?
                != Sign::Zero
        {
            return Ok(None);
        }
        // The ring's frame: the axis meets the plane at the origin, and the
        // unit axis points to the ball.
        let axis = unit(&plane.normal);
        let level = plane.normal.dot(&plane.point.minus(&Point3::origin()));
        let along = &(&level - &plane.normal.dot(&cylinder.origin.minus(&Point3::origin())))
            / &plane.normal.dot(&cylinder.axis);
        let origin = cylinder.origin.plus(&cylinder.axis.scaled(&along));
        let basis = basis_of(&axis, dec)?;
        let body = square_root(&cylinder.radius_sq);
        let radius = definition.radius.clone();
        let major = &body + &radius;
        let plain = Rim::new(
            origin.clone(),
            axis.clone(),
            basis.clone(),
            body.clone(),
            radius.clone(),
            Vec::new(),
            dec,
        )?;
        let (frame_origin, axes) = plain.frame(dec)?;
        let frame = Frame {
            origin: frame_origin,
            axes,
            fan: plain.fan(dec)?,
        };
        let boxes = plain.boxes(dec)?;
        let height = |p: &Point3| axis.dot(&p.minus(&origin));
        let across = |p: &Point3| {
            let d = p.minus(&origin);
            [d.dot(&basis[0]), d.dot(&basis[1])]
        };
        // Every other source face is shown to miss the fan or is a mouth's
        // side; the run then needs no further clearance.
        let mut ends = Vec::new();
        let mut caps = Vec::new();
        for &i in source_faces {
            let face = &self.faces[i];
            if faces.contains(&i)
                || faces
                    .iter()
                    .any(|&j| face.support.same_definition(&self.faces[j].support))
            {
                continue;
            }
            ends.push(i);
            // Another face in the ring's plane bounds nothing inside it.
            if let Surface::Plane(p) = &face.surface
                && geometry::sign(&p.normal.cross(&axis).norm_squared(), dec)? == Sign::Zero
                && geometry::sign(&height(&p.origin), dec)? == Sign::Zero
            {
                continue;
            }
            if !reaches(&face.polynomial, &frame, dec)? {
                continue;
            }
            let geometry::Leaf::Primitive(primitive) = &self.leaves[face.leaf] else {
                return Err(overflow(
                    "a support that reaches the ring is not a primitive's face",
                ));
            };
            if primitive.kind != Shape::Round {
                return Err(overflow(
                    "a support that reaches the ring is not a cylinder's or cone's face",
                ));
            }
            // Every face of the leaf lies in the leaf, so a leaf that misses
            // the fan holds no face that reaches it.
            let placement = primitive.inverse_placement().inverse();
            let [r0, r1, h] = primitive.parameters.clone();
            let wide = r0.clone().max(r1.clone());
            let base = placement.apply_point(&Point3::origin());
            let top = placement.apply_point(&Point3::new(
                Real::zero(),
                Real::zero(),
                Real::from(h.clone()),
            ));
            let direction = top.minus(&base);
            let parallel =
                geometry::sign(&direction.cross(&axis).norm_squared(), dec)? == Sign::Zero;
            if parallel {
                let levels = [height(&base), height(&top)];
                let [x, z] = across(&base);
                let o = square_root(&(&x.square() + &z.square()));
                let wide = Real::from(wide.clone());
                let mut below = true;
                let mut over = true;
                for level in &levels {
                    below &= geometry::sign(level, dec)? == Sign::Negative;
                    over &= geometry::sign(&(level - &radius), dec)? == Sign::Positive;
                }
                if below
                    || over
                    || geometry::sign(&(&(&o - &wide) - &major), dec)? == Sign::Positive
                    || geometry::sign(&(&body - &(&o + &wide)), dec)? == Sign::Positive
                {
                    continue;
                }
            } else {
                let mut corners = Vec::with_capacity(8);
                for corner in 0..8 {
                    let pick = |bit: usize, lo: &Rational, hi: &Rational| {
                        if (corner >> bit) & 1 == 0 {
                            lo.clone()
                        } else {
                            hi.clone()
                        }
                    };
                    let local = Point3::new(
                        pick(0, &-&wide, &wide).into(),
                        pick(1, &-&wide, &wide).into(),
                        pick(2, &Rational::zero(), &h).into(),
                    );
                    let p = placement.apply_point(&local);
                    corners.push([
                        geometry::enclose(&p.x, BITS, dec)?,
                        geometry::enclose(&p.y, BITS, dec)?,
                        geometry::enclose(&p.z, BITS, dec)?,
                    ]);
                }
                let extent = hull_points(&corners);
                if boxes.iter().all(|b| {
                    (0..3).any(|k| extent[k].hi() < b[k].lo() || extent[k].lo() > b[k].hi())
                }) {
                    continue;
                }
                return Err(overflow(
                    "a cylinder or cone across the ring's axis reaches its region",
                ));
            }
            if !matches!(face.surface, Surface::Cylinder(_) | Surface::Cone(_)) {
                return Err(overflow(
                    "a cap of a cylinder or cone reaches the ring's region",
                ));
            }
            let levels = [height(&base), height(&top)];
            let signs = [
                geometry::sign(&levels[0], dec)?,
                geometry::sign(&levels[1], dec)?,
            ];
            if signs.contains(&Sign::Positive) {
                return Err(overflow(
                    "a cylinder or cone reaches the ring's region above its plane",
                ));
            }
            let at = signs
                .iter()
                .position(|s| *s == Sign::Zero)
                .expect("a leaf that is not below the plane has a cap in it");
            let (centre, radius_at, other) = if at == 0 {
                (base, r0, r1)
            } else {
                (top, r1, r0)
            };
            caps.push(Cap {
                centre,
                radius: Real::from(radius_at.clone()),
                depth: Real::from(h.clone()),
                flare: Real::from(other - radius_at),
            });
        }
        // Mouths in the plane must keep apart, so each hole's mouth is one
        // whole disc of the plane.
        for (k, a) in caps.iter().enumerate() {
            for b in &caps[..k] {
                dec.charge_work(1)?;
                let apart =
                    &a.centre.minus(&b.centre).norm_squared() - &(&a.radius + &b.radius).square();
                if geometry::sign(&apart, dec)? != Sign::Positive {
                    return Err(overflow("two mouths in the ring's plane meet"));
                }
            }
        }
        let mut mouths = Vec::new();
        for cap in &caps {
            let offset = cap.centre.minus(&origin);
            let o = square_root(&offset.norm_squared());
            // A mouth clear of the centre circle, or touching it at one point,
            // leaves the plane contact on the plane; one under the cylinder
            // misses the fan.
            if geometry::sign(&(&(&o - &cap.radius) - &major), dec)? != Sign::Negative
                || geometry::sign(&(&(&o + &cap.radius) - &body), dec)? != Sign::Positive
            {
                continue;
            }
            if geometry::sign(&(&o - &major), dec)? != Sign::Positive
                || geometry::sign(&(&(&o - &cap.radius) - &body), dec)? != Sign::Positive
            {
                return Err(Refusal::new(
                    "field.rolling-path-reach",
                    "a mouth the centre circle crosses must lie beyond it and outside the cylinder",
                ));
            }
            // Below the plane at a point inside the mouth and the fan: air
            // means a hole, material a buried cap.
            let at = &(&(&o - &cap.radius) + &major) / &(&Real::from_i64(2) * &o);
            let inside = origin.plus(&offset.scaled(&at));
            let down = axis.scaled(&Real::from_i64(-1));
            match self.probe(root, &inside, &down, dec)? {
                Location::Outside => {}
                Location::Inside => continue,
                _ => return Err(overflow("the material under a mouth is unresolved")),
            }
            if self.probe(root, &inside, &axis, dec)? != Location::Outside {
                return Err(overflow("the ring's side of a mouth is not open"));
            }
            let [x, z] = across(&cap.centre);
            mouths.push(Mouth {
                centre: [x, z],
                radius: cap.radius.clone(),
                depth: cap.depth.clone(),
                flare: cap.flare.clone(),
            });
        }
        let rim = Rim::new(origin, axis, basis, body, radius, mouths, dec)?;
        let piece = Piece::new(
            Kind::Rim(rim),
            self.blend_origin(faces, step)?,
            false,
            [&supports[0], &supports[1]],
            dec,
        )?;
        Ok(Some(Run::new(piece, faces, ends)))
    }
}
