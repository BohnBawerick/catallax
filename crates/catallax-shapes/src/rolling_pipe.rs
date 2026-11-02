//! The exact surface a ball sweeps rolling between a plane and an oblique
//! circular cylinder.
//!
//! Take the unit axis `a`, the outward plane normal `n = alpha a + beta w` and
//! `s = a x w`, and write `m(t) = cos t w + sin t s`. Ball centres lie on the
//! plane's offset and on a cylinder of radius `rho`, so they form an ellipse:
//!
//! ```text
//! c(t) = base + cos t (rho w + lift a) + sin t rho s
//! ```
//!
//! Its tangent is normal to both `n` and `m(t)`, so the contact rays `n` and
//! `m(t)` span the normal plane and each contact curve is an ellipse too.
//! `m(t)` is a rational unit vector in `tan(t/2)`, and a circle through a known
//! rational point has a rational chart by projection from its antipode. The
//! pipe is therefore exactly rational, of degree (4, 2), and nothing is fitted.
//! Its construction deviation from the rolling-ball definition is zero.
//!
//! A chart row `v = 0` is the plane contact and `v = 1` is the wall contact.
//! `u` runs over four spans, one quarter turn of `t` each, from `t = 0`.

use std::sync::Arc;

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};

use crate::{
    BSplineSurface, Cylinder, Ellipse, EnvelopeGrid, Isometry, Plane, Point3, Pose, Vector3, math,
};

mod grid;
mod line;
mod query;

const ENCLOSURE_BITS: u32 = 256;

fn invalid(detail: impl Into<String>) -> Refusal {
    Refusal::new("shapes.rolling-pipe", detail)
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &math::exact_sqrt(&v.norm_squared())))
}
fn binomial(n: usize, k: usize) -> i64 {
    (0..k).fold(1, |value, i| value * (n - i) as i64 / (i as i64 + 1))
}

// A polynomial on the unit square in the tensor Bernstein basis, `u`-major.
// Products keep every coefficient exact, and positive coefficients prove a
// positive polynomial over the whole square.
#[derive(Clone)]
struct Bernstein {
    degrees: [usize; 2],
    values: Vec<Real>,
}
impl Bernstein {
    fn new(degrees: [usize; 2], values: Vec<Real>) -> Self {
        debug_assert_eq!(values.len(), (degrees[0] + 1) * (degrees[1] + 1));
        Self { degrees, values }
    }
    fn in_u(values: [Real; 3]) -> Self {
        Self::new([2, 0], values.into())
    }
    fn in_v(values: &[i64]) -> Self {
        Self::new(
            [0, values.len() - 1],
            values.iter().map(|v| Real::from_i64(*v)).collect(),
        )
    }
    fn at(&self, i: usize, j: usize) -> &Real {
        &self.values[i * (self.degrees[1] + 1) + j]
    }
    fn times(&self, other: &Self) -> Self {
        let degrees = [0, 1].map(|axis| self.degrees[axis] + other.degrees[axis]);
        let mut values = vec![Real::zero(); (degrees[0] + 1) * (degrees[1] + 1)];
        for i in 0..=self.degrees[0] {
            for j in 0..=self.degrees[1] {
                for p in 0..=other.degrees[0] {
                    for q in 0..=other.degrees[1] {
                        let scale = Rational::from_i64(
                            binomial(self.degrees[0], i)
                                * binomial(self.degrees[1], j)
                                * binomial(other.degrees[0], p)
                                * binomial(other.degrees[1], q),
                        ) / Rational::from_i64(
                            binomial(degrees[0], i + p) * binomial(degrees[1], j + q),
                        );
                        let at = (i + p) * (degrees[1] + 1) + j + q;
                        values[at] =
                            &values[at] + &(&(self.at(i, j) * other.at(p, q)) * &Real::from(scale));
                    }
                }
            }
        }
        Self::new(degrees, values)
    }
    fn plus(&self, other: &Self) -> Self {
        assert_eq!(self.degrees, other.degrees);
        Self::new(
            self.degrees,
            self.values
                .iter()
                .zip(&other.values)
                .map(|(a, b)| a + b)
                .collect(),
        )
    }
    fn scaled(&self, factor: &Real) -> Self {
        Self::new(
            self.degrees,
            self.values.iter().map(|v| v * factor).collect(),
        )
    }
}

/// The exact frame and numbers of a pipe, for consumers that integrate over it.
/// `m(t) = cos t across + sin t side`, and ball centres are
/// `base + cos t (offset across + lift axis) + sin t offset side`.
#[derive(Debug)]
pub struct PipeLaw {
    /// Point on the cylinder axis in the plane's offset.
    pub base: Point3,
    /// Unit cylinder axis.
    pub axis: Vector3,
    /// Unit vector square to the axis, in the plane of the axis and normal.
    pub across: Vector3,
    /// `axis x across`.
    pub side: Vector3,
    /// Unit plane normal, `along axis + tilt across`.
    pub normal: Vector3,
    /// Cylinder radius.
    pub body: Real,
    /// Radius of the centre path's cylinder.
    pub offset: Real,
    /// Ball radius.
    pub radius: Real,
    /// Axial height of the centre path at `t = 0`.
    pub lift: Real,
    /// Sign of the plane contact ray: the ray is `reach normal`.
    pub reach: Real,
    /// Sign of the wall contact ray: the ray is `wall m(t)`.
    pub wall: Real,
    /// Component of the plane normal along the axis.
    pub along: Real,
    /// Component of the plane normal across the axis.
    pub tilt: Real,
}

#[derive(Debug)]
struct Data {
    law: PipeLaw,
    carrier: BSplineSurface,
}

/// The exact blend of a closed regular contact branch between a plane and a
/// cylinder that meet obliquely.
#[derive(Clone, Debug)]
pub struct RollingPipe {
    data: Arc<Data>,
}

impl RollingPipe {
    /// Reduces the contact to its three ellipses and exact carrier. `same[i]`
    /// says whether the input face follows support `i`'s natural normal, and
    /// `convex` puts the ball on the material side. Refuses a plane that is not
    /// strictly oblique and a radius at which the swept band could fold.
    pub fn new(
        plane: &Plane,
        cylinder: &Cylinder,
        same: [bool; 2],
        radius: Real,
        convex: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&plane.origin, dec)?;
        math::validate_vector(&plane.normal, dec)?;
        math::validate_point(&cylinder.origin, dec)?;
        math::validate_vector(&cylinder.axis, dec)?;
        math::validate(&radius, dec)?;
        math::positive(&radius, "shapes.non-positive-radius", dec)?;
        math::positive(&cylinder.radius_sq, "shapes.non-positive-radius", dec)?;
        math::directions(&plane.normal, &cylinder.axis, dec)?;
        let reach = Real::from_i64(if convex { 1 } else { -1 });
        let toward = |same: bool| Real::from_i64(if same { 1 } else { -1 });
        let axis = unit(&cylinder.axis);
        let normal = unit(&plane.normal).scaled(&toward(same[0]));
        let along = normal.dot(&axis);
        if math::sign(&along, dec)? == Sign::Zero {
            return Err(invalid(
                "a plane parallel to the axis has straight contacts",
            ));
        }
        let tilt_sq = &Real::one() - &along.square();
        if math::sign(&tilt_sq, dec)? != Sign::Positive {
            return Err(invalid("a plane normal to the axis has circular contacts"));
        }
        let tilt = math::exact_sqrt(&tilt_sq);
        let across = normal
            .minus(&axis.scaled(&along))
            .scaled(&(&Real::one() / &tilt));
        let side = axis.cross(&across);
        let body = cylinder.radius();
        let wall = &reach * &toward(same[1]);
        let offset = &body - &(&wall * &radius);
        if math::sign(&offset, dec)? != Sign::Positive {
            return Err(invalid(
                "the cylinder offset reaches its axis and needs contact continuation",
            ));
        }
        // A ball on the far side of the wall rolls inside the centre ellipse,
        // whose smallest curvature radius is `offset |along|`.
        if math::sign(&wall, dec)? == Sign::Negative
            && math::sign(
                &(&(&offset.square() * &along.square()) - &radius.square()),
                dec,
            )? != Sign::Positive
        {
            return Err(invalid(
                "the rolling ball folds against the tightest curve of its centre path",
            ));
        }
        let height = cylinder.origin.minus(&plane.origin).dot(&normal);
        let base = cylinder
            .origin
            .plus(&axis.scaled(&(&-(&(&reach * &radius) + &height) / &along)));
        let lift = -(&(&tilt * &offset) / &along);
        let law = PipeLaw {
            base,
            axis,
            across,
            side,
            normal,
            body,
            offset,
            radius,
            lift,
            reach,
            wall,
            along,
            tilt,
        };
        let carrier = law.carrier(dec)?;
        Ok(Self {
            data: Arc::new(Data { law, carrier }),
        })
    }

    /// Positive expression-node identity for deterministic cache reuse.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
    /// The exact law the carrier was built from.
    pub fn law(&self) -> &PipeLaw {
        &self.data.law
    }
    /// The exact rational carrier over `[0, 4] x [0, 1]`.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }
    /// Positive rolling radius in millimetres.
    pub fn radius(&self) -> &Real {
        &self.data.law.radius
    }
    /// The ball-centre ellipse.
    pub fn spine(&self, dec: &mut Decider) -> Result<Ellipse, Refusal> {
        let d = &self.data.law;
        d.ellipse(d.base.clone(), &d.offset, d.normal.clone(), dec)
    }
    /// The contact ellipse on the plane.
    pub fn plane_rail(&self, dec: &mut Decider) -> Result<Ellipse, Refusal> {
        let d = &self.data.law;
        let center = d.base.plus(&d.normal.scaled(&(&d.reach * &d.radius)));
        d.ellipse(center, &d.offset, d.normal.clone(), dec)
    }
    /// The contact ellipse on the cylinder.
    pub fn wall_rail(&self, dec: &mut Decider) -> Result<Ellipse, Refusal> {
        let d = &self.data.law;
        let normal = d.major(&d.body).cross(&d.side);
        d.ellipse(d.base.clone(), &d.body, normal, dec)
    }
    /// The ball centre at chart `u`, read from the law and not the carrier.
    pub fn centre(&self, u: &Rational) -> Point3 {
        let d = &self.data.law;
        let [cos, sin] = turn(u);
        d.base
            .plus(&d.major(&d.offset).scaled(&cos))
            .plus(&d.side.scaled(&(&d.offset * &sin)))
    }
    /// The chart point, read from the law and not the carrier. Both rays are
    /// unit, and the section runs from the plane contact to the wall contact.
    pub fn point(&self, u: &Rational, v: &Rational) -> Point3 {
        let d = &self.data.law;
        let [cos, sin] = turn(u);
        let first = d.normal.scaled(&d.reach);
        let last = d
            .across
            .scaled(&cos)
            .plus(&d.side.scaled(&sin))
            .scaled(&d.wall);
        let v = Real::from(v.clone());
        let ray = first
            .scaled(&(&Real::from_i64(2) - &v))
            .plus(&last.scaled(&v));
        let chord = first.dot(&ray);
        let unit = ray
            .scaled(&(&(&Real::from_i64(2) * &chord) / &ray.norm_squared()))
            .minus(&first);
        self.centre(u).plus(&unit.scaled(&d.radius))
    }
    /// Whether `edge` is the exact input edge of this contact, the plane section
    /// of the cylinder on the sheet this pipe rolls along.
    pub fn is_edge(&self, edge: &Ellipse, dec: &mut Decider) -> Result<bool, Refusal> {
        let d = &self.data.law;
        // The edge sits one rolling height from the plane rail along the axis.
        let shift = &(&d.reach * &d.radius) / &d.along;
        let center = d.base.plus(&d.axis.scaled(&shift));
        let semi = d
            .across
            .scaled(&d.body)
            .plus(&d.axis.scaled(&-(&(&d.tilt * &d.body) / &d.along)));
        let zero = |value: Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&value, dec)? == Sign::Zero)
        };
        Ok(zero(edge.center.minus(&center).norm_squared(), dec)?
            && zero(edge.normal.cross(&d.normal).norm_squared(), dec)?
            && zero(&edge.minor_sq - &d.body.square(), dec)?
            && zero(&edge.major_sq - &semi.norm_squared(), dec)?
            && zero(edge.major_dir.cross(&semi).norm_squared(), dec)?)
    }
    /// Whether `path` is exactly the ball-centre ellipse, as a point set.
    pub fn centres_on(&self, path: &Ellipse, dec: &mut Decider) -> Result<bool, Refusal> {
        let spine = self.spine(dec)?;
        let zero = |value: Real, dec: &mut Decider| -> Result<bool, Refusal> {
            Ok(math::sign(&value, dec)? == Sign::Zero)
        };
        Ok(zero(path.center.minus(&spine.center).norm_squared(), dec)?
            && zero(path.normal.cross(&spine.normal).norm_squared(), dec)?
            && zero(&path.major_sq - &spine.major_sq, dec)?
            && zero(&path.minor_sq - &spine.minor_sq, dec)?
            && zero(path.major_dir.cross(&spine.major_dir).norm_squared(), dec)?)
    }
    /// Whether every control of the carrier lies strictly on one side of a
    /// plane, which keeps the whole band off that plane.
    pub fn clear_of(&self, plane: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        let mut side = None;
        for control in self.carrier().controls() {
            let value = math::sign(&control.minus(&plane.origin).dot(&plane.normal), dec)?;
            if value == Sign::Zero || side.replace(value).is_some_and(|last| last != value) {
                return Ok(false);
            }
        }
        Ok(true)
    }
    /// A uniform grid with `chord` at most `tolerance`. Column zero and the
    /// last column are the seam, line zero is the plane rail and the last line
    /// is the wall rail. `u` increases counter-clockwise about `chart_normal`.
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<EnvelopeGrid, Refusal> {
        grid::fit(self, tolerance, max_points, dec)
    }
    /// Every control and weight in world coordinates, enclosed in `u`-major order.
    pub fn enclosed_controls(&self, dec: &Decider) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        let carrier = &self.data.carrier;
        carrier
            .controls()
            .iter()
            .zip(carrier.weights())
            .map(|(control, weight)| {
                Ok((
                    [
                        dec.enclose_at(&control.x, ENCLOSURE_BITS)?,
                        dec.enclose_at(&control.y, ENCLOSURE_BITS)?,
                        dec.enclose_at(&control.z, ENCLOSURE_BITS)?,
                    ],
                    dec.enclose_at(weight, ENCLOSURE_BITS)?,
                ))
            })
            .collect()
    }
    /// The direction `u` turns about, which is the sense of the chart's rails.
    pub fn chart_normal(&self) -> Vector3 {
        let d = &self.data.law;
        d.major(&d.body).cross(&d.side)
    }
    /// Exact rigid transport of the law and carrier.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }
    pub(crate) fn mapped(&self, pose: &Isometry) -> Self {
        let d = &self.data.law;
        Self {
            data: Arc::new(Data {
                law: PipeLaw {
                    base: pose.apply_point(&d.base),
                    axis: pose.apply_vector(&d.axis),
                    across: pose.apply_vector(&d.across),
                    side: pose.apply_vector(&d.side),
                    normal: pose.apply_vector(&d.normal),
                    body: d.body.clone(),
                    offset: d.offset.clone(),
                    radius: d.radius.clone(),
                    lift: d.lift.clone(),
                    reach: d.reach.clone(),
                    wall: d.wall.clone(),
                    along: d.along.clone(),
                    tilt: d.tilt.clone(),
                },
                carrier: self.data.carrier.mapped(pose),
            }),
        }
    }
}

// The cosine and sine at chart `u`: span `floor(u)` is a quarter turn, with the
// half-angle tangent `u - floor(u)` measured from its start.
fn turn(u: &Rational) -> [Real; 2] {
    let floor = u.floor();
    let tangent = u - &Rational::from(floor.clone());
    let span = floor.to_i64().map(|value| value.rem_euclid(4));
    let wide = Rational::one() + &tangent * &tangent;
    let cos = (Rational::one() - &tangent * &tangent) / wide.clone();
    let sin = (&tangent * &Rational::from_i64(2)) / wide;
    let [cos, sin] = [cos, sin].map(Real::from);
    match span {
        Some(0) => [cos, sin],
        Some(1) => [-&sin, cos],
        Some(2) => [-&cos, -&sin],
        _ => [sin, -&cos],
    }
}

impl PipeLaw {
    // The semi-axis over `w`, for an ellipse at radial scale `scale`. The lift
    // stays that of the ball centre, because every contact keeps its height.
    fn major(&self, scale: &Real) -> Vector3 {
        self.across
            .scaled(scale)
            .plus(&self.axis.scaled(&self.lift))
    }
    fn ellipse(
        &self,
        center: Point3,
        scale: &Real,
        normal: Vector3,
        dec: &mut Decider,
    ) -> Result<Ellipse, Refusal> {
        let minor_sq = scale.square();
        Ellipse::new(
            center,
            normal,
            self.major(scale),
            &minor_sq + &self.lift.square(),
            minor_sq,
            dec,
        )
    }

    // One Bezier patch per quarter turn. With W = 1 + tau^2, a unit ray
    // M / W toward the wall, `g = n.M`, and the chord direction
    // D = (2 - v) W first + v M, the second meeting of the line from `-first`
    // along D with the unit circle has denominator W Q, where
    // Q = ((2 - v)^2 + v^2) W + 2 v (2 - v) (first.M).
    // |first.M| <= tilt W keeps Q's Bernstein coefficients positive, and the
    // tilt is under one because the plane is strictly oblique.
    fn carrier(&self, dec: &mut Decider) -> Result<BSplineSurface, Refusal> {
        let at = |values: [i64; 3]| values.map(Real::from_i64);
        let wide = Bernstein::in_u(at([1, 1, 2]));
        let flat = [1, 1, 0];
        let rise = [0, 1, 2];
        let start = Bernstein::in_v(&[2, 1]);
        let stop = Bernstein::in_v(&[0, 1]);
        let mut controls = Vec::with_capacity(17 * 3);
        let mut weights = Vec::with_capacity(17 * 3);
        let mut scale = Real::one();
        for (span, (c, s)) in [(1, 0), (0, 1), (-1, 0), (0, -1)].into_iter().enumerate() {
            let cos = Bernstein::in_u(std::array::from_fn(|i| {
                Real::from_i64(c * flat[i] - s * rise[i])
            }));
            let sin = Bernstein::in_u(std::array::from_fn(|i| {
                Real::from_i64(s * flat[i] + c * rise[i])
            }));
            // `first` is `reach n` and the wall ray `wall m`, so `first.M` is
            // `reach wall tilt cos`.
            let meet = cos.scaled(&(&(&self.reach * &self.wall) * &self.tilt));
            let ends = Bernstein::in_v(&[4, 2, 2]).times(&wide);
            let cross = Bernstein::in_v(&[0, 2, 2]).times(&meet);
            let q = ends.plus(&cross);
            let h = start.times(&wide).plus(&stop.times(&meet));
            let denominator = wide.times(&q);
            let double = &self.radius * &Real::from_i64(2);
            // Chord components along `first`, `w` and `s`.
            let chord_first = start.times(&wide);
            let chord = [&cos, &sin].map(|part| stop.times(part).scaled(&self.wall));
            let along = [&cos, &sin].map(|part| part.times(&q).scaled(&self.offset));
            let major = h.times(&chord[0]).scaled(&double).plus(&along[0]);
            let minor = h.times(&chord[1]).scaled(&double).plus(&along[1]);
            let height = cos.times(&q).scaled(&self.lift);
            let toward = h
                .times(&chord_first)
                .scaled(&double)
                .plus(&denominator.scaled(&-&self.radius));
            let first = self.normal.scaled(&self.reach);
            for i in usize::from(span != 0)..=4 {
                for j in 0..=2 {
                    let weight = denominator.at(i, j);
                    let per = |value: &Bernstein| value.at(i, j) / weight;
                    controls.push(
                        self.base
                            .plus(&self.across.scaled(&per(&major)))
                            .plus(&self.side.scaled(&per(&minor)))
                            .plus(&self.axis.scaled(&per(&height)))
                            .plus(&first.scaled(&per(&toward))),
                    );
                    weights.push(weight * &scale);
                }
            }
            // A span ends with W = 2 and the next starts with W = 1, and both
            // factors of the denominator double across the joint.
            scale = &scale * &Real::from_i64(4);
        }
        let knot = |value: i64, count: usize| vec![Rational::from_i64(value); count];
        let u = [knot(0, 5), knot(1, 4), knot(2, 4), knot(3, 4), knot(4, 5)].concat();
        let v = [knot(0, 3), knot(1, 3)].concat();
        BSplineSurface::new([4, 2], [u, v], [17, 3], controls, weights, dec)
    }
}
