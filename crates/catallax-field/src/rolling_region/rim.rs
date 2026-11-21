//! A concave ring whose plane contact crosses circular mouths (4.2a step 4).
//!
//! The ball rolls on a plane and outside a perpendicular cylinder, so its
//! centre circle has radius `C = R + r` about the cylinder axis, `r` above
//! the plane (the ring of step 1). Where the plane holds a circular mouth, a
//! hole of radius `m` whose centre lies `o` from the axis, and the plane
//! contact would fall inside it, the ball leaves the plane for the mouth's
//! edge. Its centre stays on the cylinder's offset, so it keeps the radius
//! `C` and drops to the height `H = sqrt(r^2 - g^2)`, where `g = m - D` and
//! `D` is the horizontal distance from the centre to the mouth's axis: the
//! ball touches the edge circle at the point of it nearest the centre. This
//! is the P5 rim contact of M18 (`docs/M18-ROLLING-REGION.md`), there taken
//! from the B-rep's `RollingEnvelopeDefinition` for one mouth; here the field
//! derives it from its own supports for every mouth of the ring at once, so
//! one path covers the whole ring and no seam between a mouth and the ring
//! needs its own proof. With no mouth crossed it is the ring, with the
//! ring's closed forms.
//!
//! OpenCASCADE rolls a ball on a surface and an edge with
//! `BRepBlend_SurfRstConstRad` (its constant-radius function with the edge's
//! curve on the other face as the second support), walked by
//! `BRepBlend_SurfRstLineBuilder` (read at V7_8_0). SMLib names the case a
//! cliff rollover (`SM_FG_CLIFF_ROLLOVER` in `SmFilletGeom.h`), found where
//! the rail crosses an edge (`SM_FV_RAIL_X_EDGEUSE`), and steps it like its
//! other blends. Both take the edge as the second support, as here; the field
//! keeps the closed form of that contact, the tube of radius `r` about the
//! edge circle, and replaces walking and stepping with exact walls and a
//! whole-interval cover, because a traced spine gives no guaranteed range.
//!
//! The region is M18's: in the frame of the ring (origin where the axis
//! meets the plane, `y` along the axis toward the ball), points with
//! `0 <= y <= r` and `R <= rho <= C` that lie below the lower envelope of
//! every ball. A point above a ball's centre, within `r` of its axis line, is
//! above that envelope.

use super::path::{Cover, Frame, Path};
use super::revolved::{angle, hull_points, panels};
use super::{BITS, gap, interval};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Isometry, Point3, Vector3};

fn reach(detail: &str) -> Refusal {
    Refusal::new("field.rolling-path-reach", detail)
}
fn positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Positive)
}

/// A circular hole in the ring's plane that the plane contact crosses.
#[derive(Clone, Debug)]
pub struct Mouth {
    /// Centre in the ring's frame, along its two basis directions.
    pub centre: [Real; 2],
    /// Radius of the edge circle in mm.
    pub radius: Real,
    /// Depth of the hole's side face below the plane, in mm.
    pub depth: Real,
    /// The side face's radius at that depth less its radius at the plane:
    /// zero for a straight bore, negative for a countersink that narrows.
    pub flare: Real,
}

/// The ring path with its mouth contacts.
#[derive(Clone, Debug)]
pub struct Rim {
    origin: Point3,
    axis: Vector3,
    basis: [Vector3; 2],
    body: Real,
    radius: Real,
    major: Real,
    mouths: Vec<Mouth>,
    band: [[Rational; 2]; 3],
    boxes: Vec<[Ball; 3]>,
}

impl Rim {
    /// The ring of radius `body + radius` above the plane through `origin`
    /// normal to the unit `axis` (which points to the ball), with an
    /// orthonormal `basis` of the plane, and its mouths.
    ///
    /// Each mouth must lie outside the cylinder (`o - m > R`), on the far
    /// side of the centre circle from the axis (`o > C`), and be crossed by
    /// it (`o - m < C`), with a deepest gap `C - (o - m)` below `r`. The
    /// ball must stay on the edge: at the deepest station its centre lies
    /// strictly inside the edge's normal cone, `depth H > -flare g`, and it
    /// stays above the side face's far end, `r - H < depth`. The mouths'
    /// closed sectors of the ring must be pairwise apart.
    pub fn new(
        origin: Point3,
        axis: Vector3,
        basis: [Vector3; 2],
        body: Real,
        radius: Real,
        mouths: Vec<Mouth>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        if !positive(&radius, dec)? || !positive(&body, dec)? {
            return Err(reach("the rolling and cylinder radii must be positive"));
        }
        let major = &body + &radius;
        let mut sectors: Vec<(Real, Real, [Real; 2])> = Vec::with_capacity(mouths.len());
        for mouth in &mouths {
            dec.charge_work(1)?;
            let [x, z] = &mouth.centre;
            let o = (&x.square() + &z.square()).sqrt();
            let m = &mouth.radius;
            if !positive(m, dec)? || !positive(&mouth.depth, dec)? {
                return Err(reach("a mouth needs a positive radius and depth"));
            }
            if !positive(&(&(&o - m) - &body), dec)? {
                return Err(reach("a mouth must lie wholly outside the cylinder"));
            }
            if !positive(&(&o - &major), dec)? {
                return Err(reach(
                    "the centre circle must cross a mouth on its near side",
                ));
            }
            let deepest = &(&major + m) - &o;
            if !positive(&deepest, dec)? {
                return Err(reach("the centre circle does not cross the mouth"));
            }
            if !positive(&(&radius - &deepest), dec)? {
                return Err(reach("the ball would fall through the mouth"));
            }
            let height = (&radius.square() - &deepest.square()).sqrt();
            // In the mouth's meridian the side face's outward normal is
            // (-depth, -flare) in (radial, up); the centre lies at
            // (-g, H) from the edge. It stays inside the normal cone when
            // its direction is not past that normal.
            if !positive(
                &(&(&mouth.depth * &height) + &(&mouth.flare * &deepest)),
                dec,
            )? {
                return Err(reach(
                    "the ball would leave the mouth's edge for its side face",
                ));
            }
            if !positive(&(&mouth.depth - &(&radius - &height)), dec)? {
                return Err(reach("the ball would reach below the mouth's side face"));
            }
            // Half-angle of the mouth's sector: the centre circle meets the
            // edge circle where cos = (C^2 + o^2 - m^2) / (2 C o).
            let cosine = &(&(&major.square() + &o.square()) - &m.square())
                / &(&(&Real::from_i64(2) * &major) * &o);
            let sine = (&Real::one() - &cosine.square()).sqrt();
            let toward = [x / &o, z / &o];
            for (other_cos, other_sin, other) in &sectors {
                dec.charge_work(1)?;
                // Two closed sectors are apart when the angle between their
                // directions exceeds the sum of their half-angles, each in
                // (0, pi/2), as `Sectors::new` proves it.
                let limit = &(&cosine * other_cos) - &(&sine * other_sin);
                let dot = &(&toward[0] * &other[0]) + &(&toward[1] * &other[1]);
                if !positive(&(&limit - &dot), dec)? {
                    return Err(reach("two mouths' sectors of the ring touch or overlap"));
                }
            }
            sectors.push((cosine, sine, toward));
        }
        let mut rim = Self {
            origin,
            axis,
            basis,
            body,
            radius,
            major,
            mouths,
            band: std::array::from_fn(|_| [Rational::zero(), Rational::zero()]),
            boxes: Vec::new(),
        };
        rim.boxes = rim.panel_boxes(dec)?;
        let band = hull_points(&rim.boxes);
        rim.band = band.map(|b| [b.lo(), b.hi()]);
        Ok(rim)
    }

    /// Ball radius in mm.
    pub fn radius(&self) -> &Real {
        &self.radius
    }

    /// The crossed mouths.
    pub fn mouths(&self) -> &[Mouth] {
        &self.mouths
    }

    /// The fan over each of 64 panels in the ring's frame: `R <= rho <= C`,
    /// `0 <= y <= r` and the panel's angles.
    pub fn fan(&self, dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        let pi = Ball::pi_checked(BITS, dec)?;
        let zero = interval(&Rational::zero(), &Rational::zero());
        let one = interval(&Rational::one(), &Rational::one());
        let plane = [
            [one.clone(), zero.clone(), zero.clone()],
            [zero.clone(), zero.clone(), one],
        ];
        let rho = enclose(&self.body, BITS, dec)?.hull(&enclose(&self.major, BITS, dec)?);
        let height = zero.hull(&enclose(&self.radius, BITS, dec)?);
        Ok(panels()
            .map(|u| {
                let (cos, sin, _) = angle(&pi, &plane, &u);
                [&rho * &cos, height.clone(), &rho * &sin]
            })
            .collect())
    }

    /// The ring's frame as world balls: origin, then the first basis
    /// direction, the axis and the second basis direction.
    pub fn frame(&self, dec: &mut Decider) -> Result<([Ball; 3], [[Ball; 3]; 3]), Refusal> {
        Ok((
            point_balls(&self.origin, dec)?,
            [
                vector_balls(&self.basis[0], dec)?,
                vector_balls(&self.axis, dec)?,
                vector_balls(&self.basis[1], dec)?,
            ],
        ))
    }

    /// The fan's panels as world boxes.
    fn panel_boxes(&self, dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        let (origin, frame) = self.frame(dec)?;
        Ok(self
            .fan(dec)?
            .into_iter()
            .map(|local| world(&origin, &frame, &local))
            .collect())
    }

    /// Applies an exact rigid motion; the mouths are in the ring's frame and
    /// move with it.
    pub fn transform(&mut self, pose: &Isometry) {
        self.origin = pose.apply_point(&self.origin);
        self.axis = pose.apply_vector(&self.axis);
        for e in &mut self.basis {
            *e = pose.apply_vector(e);
        }
        self.band = super::moved_box(pose, &self.band);
        self.boxes = self
            .boxes
            .iter()
            .map(|b| {
                super::moved_box(pose, &b.each_ref().map(|x| [x.lo(), x.hi()]))
                    .map(|[lo, hi]| interval(&lo, &hi))
            })
            .collect();
    }

    /// An axis-aligned box holding the region.
    pub fn hull(&self, _dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        Ok(self.band.each_ref().map(|[lo, hi]| interval(lo, hi)))
    }

    /// Boxes that together hold the region, one per panel.
    pub fn boxes(&self, _dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        Ok(self.boxes.clone())
    }

    fn rough(&self) -> Rough {
        let f = |x: &Real| x.to_f64();
        let fv = |v: &Vector3| [f(&v.x), f(&v.y), f(&v.z)];
        Rough {
            origin: [f(&self.origin.x), f(&self.origin.y), f(&self.origin.z)],
            frame: [fv(&self.basis[0]), fv(&self.axis), fv(&self.basis[1])],
            body: f(&self.body),
            radius: f(&self.radius),
            major: f(&self.major),
            mouths: self
                .mouths
                .iter()
                .map(|m| ([f(&m.centre[0]), f(&m.centre[1])], f(&m.radius)))
                .collect(),
        }
    }

    /// Tentative binary64 value for display: negative inside the region.
    /// The square root is the caller's.
    pub(crate) fn fast(&self, p: [f64; 3], sqrt: fn(f64) -> f64) -> f64 {
        let rough = self.rough();
        let local = rough.local(p);
        let r = rough.radius;
        let mut nearest = f64::INFINITY;
        let mut above = false;
        const SAMPLES: usize = 512;
        for i in 0..SAMPLES {
            let u = -1.0 + 2.0 * (i as f64 + 0.5) / SAMPLES as f64;
            let c = rough.centre(u, sqrt);
            let flat = (local[0] - c[0]).powi(2) + (local[2] - c[2]).powi(2);
            if flat < r * r && local[1] > c[1] {
                above = true;
            }
            nearest = nearest.min(sqrt(flat + (local[1] - c[1]).powi(2)));
        }
        let rho = sqrt(local[0] * local[0] + local[2] * local[2]);
        let walls = [local[1], r - local[1], rho - rough.body, rough.major - rho];
        let fan = walls.iter().all(|w| *w >= 0.0);
        if fan && !above && nearest > r {
            walls.iter().map(|w| -w).fold(r - nearest, f64::max)
        } else {
            (nearest - r).abs().max(f64::MIN_POSITIVE)
        }
    }

    fn local_vector(&self, point: &Point3) -> [Real; 3] {
        let d = point.minus(&self.origin);
        [
            d.dot(&self.basis[0]),
            d.dot(&self.axis),
            d.dot(&self.basis[1]),
        ]
    }
}

/// A box of the ring's frame as a world box.
pub fn world(origin: &[Ball; 3], frame: &[[Ball; 3]; 3], local: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| {
        (0..3).fold(origin[i].clone(), |sum, j| {
            &sum + &(&local[j] * &frame[j][i])
        })
    })
}

fn point_balls(p: &Point3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        enclose(&p.x, BITS, dec)?,
        enclose(&p.y, BITS, dec)?,
        enclose(&p.z, BITS, dec)?,
    ])
}
fn vector_balls(v: &Vector3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        enclose(&v.x, BITS, dec)?,
        enclose(&v.y, BITS, dec)?,
        enclose(&v.z, BITS, dec)?,
    ])
}

/// A [`Rim`] in binary64.
struct Rough {
    origin: [f64; 3],
    frame: [[f64; 3]; 3],
    body: f64,
    radius: f64,
    major: f64,
    mouths: Vec<([f64; 2], f64)>,
}
impl Rough {
    fn local(&self, p: [f64; 3]) -> [f64; 3] {
        let d = [
            p[0] - self.origin[0],
            p[1] - self.origin[1],
            p[2] - self.origin[2],
        ];
        self.frame.map(|e| d[0] * e[0] + d[1] * e[1] + d[2] * e[2])
    }
    /// The centre at `u` in the ring's frame.
    fn centre(&self, u: f64, sqrt: fn(f64) -> f64) -> [f64; 3] {
        let t = std::f64::consts::PI * u;
        let (x, z) = (self.major * t.cos(), self.major * t.sin());
        let mut g: f64 = 0.0;
        for ([mx, mz], m) in &self.mouths {
            let d = sqrt((x - mx).powi(2) + (z - mz).powi(2));
            g = g.max(m - d);
        }
        let h = sqrt((self.radius * self.radius - g * g).max(0.0));
        [x, h, z]
    }
    /// The parameter of the centre nearest `p` (in the frame): the best of
    /// 1024 samples, polished by a golden-section search.
    fn nearest(&self, p: [f64; 3]) -> f64 {
        const SAMPLES: usize = 1024;
        let gap = |u: f64| {
            let c = self.centre(u, f64::sqrt);
            (p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2) + (p[2] - c[2]).powi(2)
        };
        let step = 2.0 / SAMPLES as f64;
        let mut best = (f64::INFINITY, 0.0);
        for i in 0..SAMPLES {
            let u = -1.0 + step * (i as f64 + 0.5);
            let g = gap(u);
            if g < best.0 {
                best = (g, u);
            }
        }
        let ratio = (5f64.sqrt() - 1.0) / 2.0;
        let (mut lo, mut hi) = ((best.1 - step).max(-1.0), (best.1 + step).min(1.0));
        for _ in 0..80 {
            let a = hi - ratio * (hi - lo);
            let b = lo + ratio * (hi - lo);
            if gap(a) < gap(b) {
                hi = b;
            } else {
                lo = a;
            }
        }
        (lo + hi) / 2.0
    }
}

/// Interval cover of a [`Rim`]; `u` in `[-1, 1]` is the angle `pi u` from
/// the first basis direction, in the ring's frame.
#[derive(Clone, Debug)]
pub struct RimCover {
    pi: Ball,
    body: Ball,
    radius: Ball,
    major: Ball,
    // Centre, radius and least distance `o - C` from the centre circle of
    // each mouth.
    mouths: Vec<([Ball; 2], Ball, Ball)>,
    // The deepest gap of all mouths: no station's gap exceeds it.
    cap: Ball,
}

/// One panel of the chart.
struct Chart {
    cos: Ball,
    sin: Ball,
    height: Ball,
    slope: Ball,
    foot: [Ball; 2],
}

impl RimCover {
    fn chart(&self, u: &Ball) -> Result<Chart, Refusal> {
        let zero = interval(&Rational::zero(), &Rational::zero());
        let one = interval(&Rational::one(), &Rational::one());
        let plane = [
            [one.clone(), zero.clone(), zero.clone()],
            [zero.clone(), zero.clone(), one.clone()],
        ];
        let (cos, sin, _) = angle(&self.pi, &plane, u);
        let x = &self.major * &cos;
        let z = &self.major * &sin;
        // d/dangle of the horizontal centre.
        let dx = -&(&self.major * &sin);
        let dz = &self.major * &cos;
        let mut g = zero.clone();
        let mut lean = zero.clone();
        let mut foot = [x.clone(), z.clone()];
        for ([mx, mz], m, near) in &self.mouths {
            let wx = &x - mx;
            let wz = &z - mz;
            let distance = (&wx.square() + &wz.square())
                .sqrt()
                .ok_or_else(|| gap("undefined mouth distance"))?;
            // No station of the circle comes nearer the mouth's axis than
            // `o - C`, which is positive.
            let distance = interval(&distance.lo().max(near.lo()), &distance.hi().max(near.lo()));
            let raw = m - &distance;
            if raw.hi() <= Rational::zero() {
                continue;
            }
            // The gap is m - D where the centre is over the mouth and zero
            // elsewhere; at most one mouth is crossed at any station, so the
            // sum over mouths holds the station's gap.
            let part = interval(&raw.lo().max(Rational::zero()), &raw.hi());
            g = &g + &part;
            // The edge contact is the edge point nearest the centre:
            // M + m (c - M) / D = c + (c - M)(m / D - 1).
            let ratio = m
                .checked_div(&distance)
                .ok_or_else(|| gap("undefined mouth ratio"))?;
            let over = &ratio - &one;
            let over = interval(
                &over.lo().max(Rational::zero()),
                &over.hi().max(Rational::zero()),
            );
            foot[0] = &foot[0] + &(&wx * &over);
            foot[1] = &foot[1] + &(&wz * &over);
            // g dD/dangle / D, with dD/dangle = (c - M).c' / D.
            let along = &(&wx * &dx) + &(&wz * &dz);
            let rate = along
                .checked_div(&distance)
                .ok_or_else(|| gap("undefined mouth rate"))?;
            lean = &lean + &(&part * &rate);
        }
        // No station's gap exceeds the deepest one.
        let g = interval(&g.lo(), &g.hi().min(self.cap.hi()));
        let height = (&self.radius.square() - &g.square())
            .sqrt()
            .ok_or_else(|| gap("undefined rim height"))?;
        // dH/dangle = -g dg/dangle / H = g (dD/dangle) / H.
        let slope = lean
            .checked_div(&height)
            .ok_or_else(|| gap("undefined rim slope"))?;
        Ok(Chart {
            cos,
            sin,
            height,
            slope,
            foot,
        })
    }
}

impl Cover for RimCover {
    fn radius(&self) -> &Ball {
        &self.radius
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        let c = self.chart(u)?;
        Ok([&self.major * &c.cos, c.height, &self.major * &c.sin])
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        let c = self.chart(u)?;
        let zero = interval(&Rational::zero(), &Rational::zero());
        let [fx, fz] = c.foot;
        Ok(Frame {
            centre: [&self.major * &c.cos, c.height.clone(), &self.major * &c.sin],
            tangent: [-&(&self.major * &c.sin), c.slope, &self.major * &c.cos],
            feet: [
                [&self.body * &c.cos, c.height, &self.body * &c.sin],
                [fx, zero, fz],
            ],
            side: None,
        })
    }
    // M18's lower envelope: a point within `r` of a ball's axis line and
    // above its centre is above the envelope.
    fn excludes(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = self.centre(u)?;
        let flat = (&p[0] - &c[0]).square() + (&p[2] - &c[2]).square();
        let all = &flat + &(&p[1] - &c[1]).square();
        let r2 = self.radius.square().lo();
        Ok(flat.hi() < r2 && (p[1].lo() > c[1].hi() || all.hi() < r2))
    }
    fn below(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = self.centre(u)?;
        let flat = (&p[0] - &c[0]).square() + (&p[2] - &c[2]).square();
        let all = &flat + &(&p[1] - &c[1]).square();
        let r2 = self.radius.square().hi();
        Ok(flat.lo() > r2 || (p[1].hi() < c[1].lo() && all.lo() > r2))
    }
}

impl Path for Rim {
    type Cover = RimCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        self.local_vector(point)
    }
    fn walls(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let [x, y, z] = local;
        let rho = &x.square() + &z.square();
        Ok(vec![
            y.clone(),
            &self.radius - y,
            &rho - &self.body.square(),
            &self.major.square() - &rho,
        ])
    }
    fn walls_are_faces(&self) -> bool {
        true
    }
    fn uniform(&self) -> bool {
        self.mouths.is_empty()
    }
    fn separation(&self, local: &[Real; 3]) -> Option<Real> {
        if !self.mouths.is_empty() {
            return None;
        }
        // The ring-torus polynomial, with the sign of (rho - C)^2 + h^2 - r^2
        // since C > r.
        let [x, y, z] = local;
        let rho = &x.square() + &z.square();
        let height = y - &self.radius;
        let major = self.major.square();
        let s = &(&(&rho + &height.square()) + &major) - &self.radius.square();
        Some(&s.square() - &(&(&Real::from_i64(4) * &major) * &rho))
    }
    fn cover(&self, dec: &mut Decider) -> Result<RimCover, Refusal> {
        let mut mouths = Vec::with_capacity(self.mouths.len());
        let mut cap = interval(&Rational::zero(), &Rational::zero());
        for mouth in &self.mouths {
            let [x, z] = &mouth.centre;
            let o = (&x.square() + &z.square()).sqrt();
            let deepest = &(&self.major + &mouth.radius) - &o;
            let deepest = enclose(&deepest, BITS, dec)?;
            cap = interval(&cap.lo(), &cap.hi().max(deepest.hi()));
            mouths.push((
                [enclose(x, BITS, dec)?, enclose(z, BITS, dec)?],
                enclose(&mouth.radius, BITS, dec)?,
                enclose(&(&o - &self.major), BITS, dec)?,
            ));
        }
        Ok(RimCover {
            pi: Ball::pi_checked(BITS, dec)?,
            body: enclose(&self.body, BITS, dec)?,
            radius: enclose(&self.radius, BITS, dec)?,
            major: enclose(&self.major, BITS, dec)?,
            mouths,
            cap,
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        if !self.mouths.is_empty() {
            return Ok(None);
        }
        // The ring: the meridian through the point; on the axis any one.
        let [x, _, z] = local;
        let rho = &x.square() + &z.square();
        let (cos, sin) = if sign(&rho, dec)? == Sign::Zero {
            (Real::one(), Real::zero())
        } else {
            let length = rho.sqrt();
            (x / &length, z / &length)
        };
        let ball = |v: &Real, dec: &mut Decider| enclose(v, BITS, dec);
        let radius = &self.radius;
        Ok(Some(Frame {
            centre: [
                ball(&(&self.major * &cos), dec)?,
                ball(radius, dec)?,
                ball(&(&self.major * &sin), dec)?,
            ],
            tangent: [
                ball(&-&sin, dec)?,
                ball(&Real::zero(), dec)?,
                ball(&cos, dec)?,
            ],
            feet: [
                [
                    ball(&(&self.body * &cos), dec)?,
                    ball(radius, dec)?,
                    ball(&(&self.body * &sin), dec)?,
                ],
                [
                    ball(&(&self.major * &cos), dec)?,
                    ball(&Real::zero(), dec)?,
                    ball(&(&self.major * &sin), dec)?,
                ],
            ],
            side: None,
        }))
    }
    fn guess(&self, local: &[Real; 3]) -> Option<Rational> {
        if self.mouths.is_empty() {
            return None;
        }
        let p = local.each_ref().map(|x| x.to_f64());
        Rational::from_f64_exact(self.rough().nearest(p))
    }
}
