//! Whole-face definition checks. Display triangles are not the measuring tool.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Circle, CurveId, Line, Point3, PointId, Surface, Vector3};
use catallax_topology::{FaceId, Loop, Model};

const BITS: u32 = 128;

/// The exact path of a rolling ball's centre. An arc retains its selected branch
/// and endpoints, including a degree-four intersection carrier.
#[derive(Clone, Debug)]
pub enum CentrePath {
    /// The complete line, not a finite segment.
    Line(Line),
    /// The complete circle.
    Circle(Circle),
    /// A finite, certified arc in the measured model's geometry store.
    Arc {
        /// Exact carrier and branch selector.
        curve: CurveId,
        /// Arc start, or absent with `to` for a whole closed curve.
        from: Option<PointId>,
        /// Arc end.
        to: Option<PointId>,
    },
}

/// A face's own definition, independent of any model-wide min/max residual.
/// Absence of a record means the face is exact by construction.
#[derive(Clone, Debug)]
pub enum FaceDefinition {
    /// Tube at `radius` from the selected ball-centre path. Distances to the two
    /// rolling supports define the *centre path*, not points on the fillet face.
    Fillet {
        /// Intersection of the outward-offset rolling supports.
        path: CentrePath,
        /// Positive rolling radius in millimetres.
        radius: Real,
    },
    /// Unsigned distance from the original analytic supporting surface.
    Offset {
        /// Original face's exact supporting definition, not a Boolean field.
        original: Surface,
        /// Non-negative distance in millimetres.
        distance: Real,
    },
}

/// Deterministic refinement limits. Exhausting cells returns the attained
/// enclosure, never a false claim that the requested width was reached.
#[derive(Clone, Debug)]
pub struct MeasurementOptions {
    /// Desired width of the worst-deviation enclosure, in millimetres.
    pub width: Rational,
    /// Maximum number of live parameter cells.
    pub max_cells: usize,
    /// Subdivisions of each certified centre-path trace piece.
    pub path_segments: usize,
}

impl Default for MeasurementOptions {
    fn default() -> Self {
        Self {
            width: Rational::from_ratio(1, 10_000_000).unwrap(),
            max_cells: 4096,
            path_segments: 64,
        }
    }
}

/// Proven enclosure of `sup(p in face) |distance(p, definition) - radius|`.
/// Both endpoints are millimetres. The upper endpoint is the local tolerance.
#[derive(Clone, Debug)]
pub struct FaceMeasurement {
    /// Lower bound on the face's worst deviation, established at exact witnesses.
    pub lower: Rational,
    /// Upper bound over the entire face, including unsampled regions.
    pub upper: Rational,
    /// Whether the requested enclosure width was achieved.
    pub converged: bool,
    /// Number of live cells covering the whole supporting parameter rectangle.
    pub cells: usize,
    /// Number of exact, trim-checked witnesses used for the lower bound.
    pub witnesses: usize,
}

fn refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new("field.face-measurement", detail)
}
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn ball(x: &Rational) -> Ball {
    Ball::from_rational(x, BITS)
}
fn enclosure(lo: Rational, hi: Rational) -> Ball {
    ball(&lo).hull(&ball(&hi))
}
fn abs(x: &Ball) -> Ball {
    enclosure(
        if x.contains_zero() {
            q(0)
        } else {
            x.lo().abs().min(x.hi().abs())
        },
        x.lo().abs().max(x.hi().abs()),
    )
}
fn norm(p: &[Ball; 3]) -> Ball {
    (&p[0].square() + &p[1].square() + p[2].square())
        .sqrt()
        .unwrap()
}
fn point(p: &Point3, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_at(&p.x, BITS)?,
        dec.enclose_at(&p.y, BITS)?,
        dec.enclose_at(&p.z, BITS)?,
    ])
}
fn unit(v: &Vector3, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    let length = v.norm_squared().sqrt();
    Ok([
        dec.enclose_at(&(&v.x / &length), BITS)?,
        dec.enclose_at(&(&v.y / &length), BITS)?,
        dec.enclose_at(&(&v.z / &length), BITS)?,
    ])
}
fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &a[0] * &b[0] + &a[1] * &b[1] + &a[2] * &b[2]
}
fn relative(p: &[Ball; 3], origin: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &p[i] - &origin[i])
}
fn radial(p: &[Ball; 3], axis: &[Ball; 3]) -> Ball {
    // Cross product avoids subtracting correlated squared lengths.
    norm(&[
        &p[1] * &axis[2] - &p[2] * &axis[1],
        &p[2] * &axis[0] - &p[0] * &axis[2],
        &p[0] * &axis[1] - &p[1] * &axis[0],
    ])
}

struct Segment {
    cover: [Ball; 3],
    witness: [Ball; 3],
}
struct Frame {
    origin: [Ball; 3],
    axis: [Ball; 3],
}
impl Frame {
    fn new(origin: &Point3, axis: &Vector3, dec: &Decider) -> Result<Self, Refusal> {
        Ok(Self {
            origin: point(origin, dec)?,
            axis: unit(axis, dec)?,
        })
    }
}
enum RadialKind {
    Line,
    Circle,
    Cylinder,
    Cone,
    Torus,
}
enum Metric {
    Plane(Frame),
    Sphere {
        center: [Ball; 3],
        radius: Ball,
    },
    Radial {
        frame: Frame,
        kind: RadialKind,
        radius: Ball,
        minor: Ball,
    },
    Arc(Vec<Segment>),
    ExtrudedArc {
        frame: Frame,
        segments: Vec<Segment>,
    },
}
impl Metric {
    fn support(
        surface: &Surface,
        options: &MeasurementOptions,
        dec: &Decider,
    ) -> Result<Self, Refusal> {
        let zero = ball(&q(0));
        let (origin, axis, kind, radius, minor) = match surface {
            Surface::Plane(s) => return Ok(Self::Plane(Frame::new(&s.origin, &s.normal, dec)?)),
            Surface::Sphere(s) => {
                return Ok(Self::Sphere {
                    center: point(&s.center, dec)?,
                    radius: dec.enclose_at(&s.radius, BITS)?,
                });
            }
            Surface::Cylinder(s) => (
                &s.origin,
                &s.axis,
                RadialKind::Cylinder,
                dec.enclose_at(&s.radius, BITS)?,
                zero,
            ),
            Surface::Cone(s) => (
                &s.apex,
                &s.axis,
                RadialKind::Cone,
                dec.enclose_at(&s.slope, BITS)?,
                zero,
            ),
            Surface::Torus(s) => (
                &s.center,
                &s.axis,
                RadialKind::Torus,
                dec.enclose_at(&s.major, BITS)?,
                dec.enclose_at(&s.minor, BITS)?,
            ),
            Surface::EllipticCylinder(s) => {
                let tau = Ball::pi_checked(BITS, dec)? * ball(&q(2));
                let mut segments = Vec::with_capacity(options.path_segments);
                for i in 0..options.path_segments {
                    dec.check_deadline()?;
                    let lo = Rational::from_ratio(i as i64, options.path_segments as i64).unwrap();
                    let hi =
                        Rational::from_ratio((i + 1) as i64, options.path_segments as i64).unwrap();
                    let mid = (&lo + &hi).mul_pow2(-1);
                    segments.push(Segment {
                        cover: surface.evaluate(&(&enclosure(lo, hi) * &tau), &zero, dec)?,
                        witness: surface.evaluate(&(&ball(&mid) * &tau), &zero, dec)?,
                    });
                }
                return Ok(Self::ExtrudedArc {
                    frame: Frame::new(&s.origin, &s.axis, dec)?,
                    segments,
                });
            }
            _ => {
                return Err(refusal(
                    "original surface has no certified Euclidean distance formula",
                ));
            }
        };
        Ok(Self::Radial {
            frame: Frame::new(origin, axis, dec)?,
            kind,
            radius,
            minor,
        })
    }

    fn distance(&self, p: &[Ball; 3], dec: &Decider) -> Result<Ball, Refusal> {
        dec.check_deadline()?;
        Ok(match self {
            Self::Plane(frame) => abs(&dot(&relative(p, &frame.origin), &frame.axis)),
            Self::Sphere { center, radius } => abs(&(&norm(&relative(p, center)) - radius)),
            Self::Radial {
                frame,
                kind,
                radius,
                minor,
            } => {
                let d = relative(p, &frame.origin);
                let rho = radial(&d, &frame.axis);
                match kind {
                    RadialKind::Line => rho,
                    RadialKind::Cylinder => abs(&(&rho - radius)),
                    RadialKind::Circle | RadialKind::Torus => {
                        let distance = ((&rho - radius).square() + dot(&d, &frame.axis).square())
                            .sqrt()
                            .unwrap();
                        if matches!(kind, RadialKind::Torus) {
                            abs(&(&distance - minor))
                        } else {
                            distance
                        }
                    }
                    RadialKind::Cone => {
                        let z = dot(&d, &frame.axis);
                        let t = (&(&rho * radius) + &z)
                            .checked_div(&(radius.square() + ball(&q(1))))
                            .unwrap();
                        let t = enclosure(t.lo().max(q(0)), t.hi().max(q(0)));
                        ((rho - radius * &t).square() + (z - t).square())
                            .sqrt()
                            .unwrap()
                    }
                }
            }
            Self::Arc(segments) => segment_distance(p, segments, dec)?,
            Self::ExtrudedArc { frame, segments } => {
                let axial = dot(&relative(p, &frame.origin), &frame.axis);
                let projected = std::array::from_fn(|i| &p[i] - &(&frame.axis[i] * &axial));
                segment_distance(&projected, segments, dec)?
            }
        })
    }
}

fn segment_distance(p: &[Ball; 3], segments: &[Segment], dec: &Decider) -> Result<Ball, Refusal> {
    let mut lo: Option<Rational> = None;
    let mut hi: Option<Rational> = None;
    for segment in segments {
        dec.check_deadline()?;
        let lower = norm(&relative(p, &segment.cover)).lo();
        let upper = norm(&relative(p, &segment.witness)).hi();
        lo = Some(lo.map_or(lower.clone(), |x| x.min(lower)));
        hi = Some(hi.map_or(upper.clone(), |x| x.min(upper)));
    }
    Ok(enclosure(
        lo.ok_or_else(|| refusal("empty centre path"))?,
        hi.unwrap(),
    ))
}

fn segments(pieces: &[TracePiece], count: usize, dec: &Decider) -> Result<Vec<Segment>, Refusal> {
    let mut result = Vec::new();
    for piece in pieces {
        for i in 0..count {
            dec.check_deadline()?;
            let lo = Rational::from_ratio(i as i64, count as i64).unwrap();
            let hi = Rational::from_ratio((i + 1) as i64, count as i64).unwrap();
            let mid = (&lo + &hi).mul_pow2(-1);
            result.push(Segment {
                cover: piece.evaluate(&enclosure(lo, hi), dec)?,
                witness: piece.evaluate(&ball(&mid), dec)?,
            });
        }
    }
    Ok(result)
}

fn prepare(
    definition: &FaceDefinition,
    model: &Model,
    options: &MeasurementOptions,
    dec: &mut Decider,
) -> Result<(Metric, Ball), Refusal> {
    let (metric, distance) = match definition {
        FaceDefinition::Offset { original, distance } => {
            (Metric::support(original, options, dec)?, distance)
        }
        FaceDefinition::Fillet { path, radius } => {
            let metric = match path {
                CentrePath::Line(line) => Metric::Radial {
                    frame: Frame::new(&line.origin, &line.dir, dec)?,
                    kind: RadialKind::Line,
                    radius: ball(&q(0)),
                    minor: ball(&q(0)),
                },
                CentrePath::Circle(circle) => Metric::Radial {
                    frame: Frame::new(&circle.center, &circle.normal, dec)?,
                    kind: RadialKind::Circle,
                    radius: dec.enclose_at(&circle.radius_sq.sqrt(), BITS)?,
                    minor: ball(&q(0)),
                },
                CentrePath::Arc { curve, from, to } => {
                    let traced = trace(model.geom(), *curve, *from, *to, dec)?;
                    if !traced.marks.is_empty() {
                        return Err(refusal("centre path trace was not exact"));
                    }
                    Metric::Arc(segments(traced.value.pieces(), options.path_segments, dec)?)
                }
            };
            (metric, radius)
        }
    };
    let distance = dec.enclose_at(distance, BITS)?;
    if distance.lo().sign() == Sign::Negative
        || (matches!(definition, FaceDefinition::Fillet { .. }) && distance.lo() <= q(0))
    {
        return Err(refusal(
            "offset distance must be non-negative and fillet radius positive",
        ));
    }
    Ok((metric, distance))
}

// A constant identity proves the entire supporting surface, including every
// possible trim, and avoids numerical cancellation at an exact rolling surface.
fn exact_zero(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    match dec.sign(value) {
        catallax_numbers::Outcome::ProvenCoincident(proof) => {
            Ok(proof.proof != catallax_numbers::ZeroProof::SeparationBound)
        }
        catallax_numbers::Outcome::Refused(error) => Err(error),
        _ => Ok(false),
    }
}
fn same_axis(
    a: &Point3,
    u: &Vector3,
    b: &Point3,
    v: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    Ok(exact_zero(&u.cross(v).norm_squared(), dec)?
        && exact_zero(&a.minus(b).cross(u).norm_squared(), dec)?)
}
fn constant(
    surface: &Surface,
    definition: &FaceDefinition,
    dec: &mut Decider,
) -> Result<Option<Ball>, Refusal> {
    let deviation = match (surface, definition) {
        (
            Surface::Torus(s),
            FaceDefinition::Fillet {
                path: CentrePath::Circle(c),
                radius,
            },
        ) if exact_zero(&s.center.minus(&c.center).norm_squared(), dec)?
            && exact_zero(&s.axis.cross(&c.normal).norm_squared(), dec)?
            && exact_zero(&(&s.major.square() - &c.radius_sq), dec)? =>
        {
            Some(dec.enclose_at(&(&s.minor - radius), BITS)?)
        }
        (
            Surface::Cylinder(s),
            FaceDefinition::Fillet {
                path: CentrePath::Line(l),
                radius,
            },
        ) if same_axis(&s.origin, &s.axis, &l.origin, &l.dir, dec)? => {
            Some(dec.enclose_at(&(&s.radius - radius), BITS)?)
        }
        (_, FaceDefinition::Offset { original, distance }) => {
            let gap = match (surface, original) {
                (Surface::Plane(a), Surface::Plane(b))
                    if exact_zero(&a.normal.cross(&b.normal).norm_squared(), dec)? =>
                {
                    Some(a.origin.minus(&b.origin).dot(&b.normal) / b.normal.norm_squared().sqrt())
                }
                (Surface::Cylinder(a), Surface::Cylinder(b))
                    if same_axis(&a.origin, &a.axis, &b.origin, &b.axis, dec)? =>
                {
                    Some(&a.radius - &b.radius)
                }
                (Surface::Sphere(a), Surface::Sphere(b))
                    if exact_zero(&a.center.minus(&b.center).norm_squared(), dec)? =>
                {
                    Some(&a.radius - &b.radius)
                }
                (Surface::Torus(a), Surface::Torus(b))
                    if exact_zero(&a.center.minus(&b.center).norm_squared(), dec)?
                        && exact_zero(&a.axis.cross(&b.axis).norm_squared(), dec)?
                        && exact_zero(&(&a.major - &b.major), dec)? =>
                {
                    Some(&a.minor - &b.minor)
                }
                _ if original.same_definition(surface) => Some(Real::zero()),
                _ => None,
            };
            gap.map(|gap| -> Result<Ball, Refusal> {
                let gap = abs(&dec.enclose_at(&gap, BITS)?);
                Ok(&gap - &dec.enclose_at(distance, BITS)?)
            })
            .transpose()?
        }
        _ => None,
    };
    Ok(deviation.map(|x| abs(&x)))
}

#[derive(Clone)]
struct Cell {
    domain: [[Rational; 2]; 2],
    upper: Rational,
    sequence: usize,
    depth: usize,
}
impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        self.upper == other.upper && self.sequence == other.sequence
    }
}
impl Eq for Cell {}
impl PartialOrd for Cell {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Cell {
    fn cmp(&self, other: &Self) -> Ordering {
        self.upper
            .cmp(&other.upper)
            .then_with(|| other.sequence.cmp(&self.sequence))
    }
}

struct Cover {
    domain: [[Rational; 2]; 2],
    angular: [bool; 2],
}

fn cover(
    model: &Model,
    face: FaceId,
    surface: &Surface,
    dec: &mut Decider,
) -> Result<Cover, Refusal> {
    let periodic = [q(0), q(1)];
    match surface {
        Surface::Torus(_) => {
            return Ok(Cover {
                domain: [periodic.clone(), periodic],
                angular: [true, true],
            });
        }
        Surface::Sphere(_) => {
            return Ok(Cover {
                domain: [
                    periodic,
                    [
                        Rational::from_ratio(-1, 4).unwrap(),
                        Rational::from_ratio(1, 4).unwrap(),
                    ],
                ],
                angular: [true, true],
            });
        }
        _ => {}
    }
    let (origin, first, second, angular) = match surface {
        Surface::Plane(s) => (
            &s.origin,
            unit(&s.ref_dir, dec)?,
            unit(&s.normal.cross(&s.ref_dir), dec)?,
            [false, false],
        ),
        Surface::Cylinder(s) => (
            &s.origin,
            unit(&s.ref_dir, dec)?,
            unit(&s.axis, dec)?,
            [true, false],
        ),
        Surface::Cone(s) => (
            &s.apex,
            unit(&s.ref_dir, dec)?,
            unit(&s.axis, dec)?,
            [true, false],
        ),
        Surface::EllipticCylinder(s) => (
            &s.origin,
            unit(&s.major_dir, dec)?,
            unit(&s.axis, dec)?,
            [true, false],
        ),
        _ => return Err(refusal("unknown face support")),
    };
    let origin = point(origin, dec)?;
    let mut bounds: Option<[Ball; 2]> = None;
    let mut include = |p: [Ball; 3]| {
        let d = relative(&p, &origin);
        let values = [dot(&d, &first), dot(&d, &second)];
        bounds = Some(match bounds.take() {
            None => values,
            Some(b) => std::array::from_fn(|i| b[i].hull(&values[i])),
        });
    };
    for lp in &model.face(face).loops {
        match model.loop_(*lp) {
            Loop::Vertex(v) => include(point(&model.geom().position(model.vertex(*v).point), dec)?),
            Loop::Coedges(cs) => {
                for c in cs {
                    let edge = model.edge(model.coedge(*c).edge);
                    let traced = trace(
                        model.geom(),
                        edge.curve,
                        edge.start.map(|v| model.vertex(v).point),
                        edge.end.map(|v| model.vertex(v).point),
                        dec,
                    )?;
                    if !traced.marks.is_empty() {
                        return Err(refusal("face boundary trace was not exact"));
                    }
                    for piece in traced.value.pieces() {
                        // Boundary-coordinate extrema enclose every linear chart
                        // coordinate on a compact plane/cylinder/cone face. Angular
                        // coordinates instead cover their complete periodic domain.
                        for part in segments(std::slice::from_ref(piece), 8, dec)? {
                            include(part.cover);
                        }
                    }
                }
            }
        }
    }
    let bounds = bounds.ok_or_else(|| refusal("non-compact face has no finite chart cover"))?;
    let domain = std::array::from_fn(|i| {
        if angular[i] {
            periodic.clone()
        } else {
            [bounds[i].lo(), bounds[i].hi()]
        }
    });
    if domain.iter().any(|d| d[0] >= d[1]) {
        return Err(refusal("face chart cover is degenerate"));
    }
    Ok(Cover { domain, angular })
}

// Global bounds on the Euclidean norms of the chart partial derivatives.
// The frame is orthonormal. Integrating along the two parameter directions
// bounds every point of a cell, including points never evaluated by the checker.
fn speeds(surface: &Surface, cover: &Cover, dec: &Decider) -> Result<[Rational; 2], Refusal> {
    let upper = |x: &Real| dec.enclose_at(x, BITS).map(|b| b.hi());
    Ok(match surface {
        Surface::Plane(_) => [q(1), q(1)],
        Surface::Cylinder(s) => [upper(&s.radius)?, q(1)],
        Surface::Cone(s) => {
            let vmax = cover.domain[1][0].abs().max(cover.domain[1][1].abs());
            [
                upper(&s.slope)? * vmax,
                upper(&(s.slope.square() + Real::one()).sqrt())?,
            ]
        }
        Surface::Sphere(s) => {
            let radius = upper(&s.radius)?;
            [radius.clone(), radius]
        }
        Surface::Torus(s) => [upper(&(&s.major + &s.minor))?, upper(&s.minor)?],
        Surface::EllipticCylinder(s) => [upper(&s.major_sq.sqrt())?, q(1)],
        _ => return Err(refusal("surface has no proven chart derivative bounds")),
    })
}

struct Check<'a> {
    surface: &'a Surface,
    metric: &'a Metric,
    distance: &'a Ball,
    angular: [bool; 2],
    tau: Ball,
    speeds: [Rational; 2],
}
impl Check<'_> {
    fn evaluate(
        &self,
        domain: &[[Rational; 2]; 2],
        dec: &Decider,
    ) -> Result<([Ball; 3], Ball), Refusal> {
        let parameters: [Ball; 2] = std::array::from_fn(|i| {
            let b = ball(&(&domain[i][0] + &domain[i][1]).mul_pow2(-1));
            if self.angular[i] { &b * &self.tau } else { b }
        });
        let p = self.surface.evaluate(&parameters[0], &parameters[1], dec)?;
        let deviation = abs(&(&self.metric.distance(&p, dec)? - self.distance));
        let radius = (0..2).fold(q(0), |sum, i| {
            let width = &domain[i][1] - &domain[i][0];
            let scale = if self.angular[i] { self.tau.hi() } else { q(1) };
            sum + (&(&self.speeds[i] * &width) * &scale).mul_pow2(-1)
        });
        // Distance to a set and |distance-r| are 1-Lipschitz. This uniform
        // chart-radius proof, not a sample maximum, certifies the whole cell.
        Ok((
            p,
            enclosure(
                (&deviation.lo() - &radius).max(q(0)),
                deviation.hi() + radius,
            ),
        ))
    }
}

/// Measures the whole B-rep face against its own definition. The upper proof
/// covers a complete supporting chart rectangle derived from certified boundary
/// traces. It may conservatively include trimmed-away regions. Lower witnesses
/// must belong to the actual trimmed face. No display mesh or min/max field is
/// used. Refinement cannot invalidate an earlier upper bound.
pub fn measure_face(
    model: &mut Model,
    face: FaceId,
    definition: Option<&FaceDefinition>,
    options: &MeasurementOptions,
    dec: &mut Decider,
) -> Result<FaceMeasurement, Refusal> {
    dec.check_deadline()?;
    if options.width <= q(0)
        || options.max_cells == 0
        || options.path_segments == 0
        || options.path_segments > i64::MAX as usize
    {
        return Err(refusal(
            "positive width, cell limit and path segment count required",
        ));
    }
    let Some(definition) = definition else {
        return Ok(FaceMeasurement {
            lower: q(0),
            upper: q(0),
            converged: true,
            cells: 0,
            witnesses: 0,
        });
    };
    let surface = model.geom().surface(model.face(face).surface).clone();
    let (metric, distance) = prepare(definition, model, options, dec)?;
    if let Some(value) = constant(&surface, definition, dec)? {
        let lower = value.lo().max(q(0));
        let upper = value.hi();
        return Ok(FaceMeasurement {
            converged: &upper - &lower <= options.width,
            lower,
            upper,
            cells: 0,
            witnesses: 0,
        });
    }
    let cover = cover(model, face, &surface, dec)?;
    let check = Check {
        surface: &surface,
        metric: &metric,
        distance: &distance,
        angular: cover.angular,
        tau: Ball::pi_checked(BITS, dec)? * ball(&q(2)),
        speeds: speeds(&surface, &cover, dec)?,
    };
    let mut lower = q(0);
    let mut witnesses = 0;
    // This witness is an exact point of the trimmed face, not a mesh vertex.
    let seed = catallax_boolean::face_point(model, face, dec)?;
    if seed.marks.is_empty() {
        let value = abs(&(&metric.distance(&point(&seed.value, dec)?, dec)? - &distance));
        lower = value.lo().max(q(0));
        witnesses += 1;
    }
    let initial = check.evaluate(&cover.domain, dec)?.1.hi();
    let mut cells = BinaryHeap::from([Cell {
        domain: cover.domain,
        upper: initial,
        sequence: 0,
        depth: 0,
    }]);
    let mut sequence = 1;
    while cells.len() < options.max_cells && &cells.peek().unwrap().upper - &lower > options.width {
        dec.check_deadline()?;
        let cell = cells.pop().unwrap();
        // Alternate normalised chart dimensions. This guarantees shrinking in
        // both dimensions even when one is angular and the other linear.
        let dimension = cell.depth % 2;
        let mid = (&cell.domain[dimension][0] + &cell.domain[dimension][1]).mul_pow2(-1);
        for side in 0..2 {
            let mut domain = cell.domain.clone();
            domain[dimension][1 - side] = mid.clone();
            let (p, value) = check.evaluate(&domain, dec)?;
            // Samples only improve the lower bound. Exact rational projection
            // plus exact trim membership prevents a supporting-surface hit from
            // being treated as a face witness.
            let target = Point3::new(
                Real::from_rational(p[0].mid()),
                Real::from_rational(p[1].mid()),
                Real::from_rational(p[2].mid()),
            );
            if let Some(witness) =
                catallax_boolean::surface_witness(&surface, &target, &options.width, dec)?
                && catallax_boolean::face_contains(model, face, &witness, dec)?
            {
                let p = point(&witness.to_real(), dec)?;
                lower = lower.max(
                    abs(&(&metric.distance(&p, dec)? - &distance))
                        .lo()
                        .max(q(0)),
                );
                witnesses += 1;
            }
            cells.push(Cell {
                domain,
                upper: cell.upper.clone().min(value.hi()),
                sequence,
                depth: cell.depth + 1,
            });
            sequence += 1;
        }
    }
    let upper = cells.peek().unwrap().upper.clone();
    Ok(FaceMeasurement {
        converged: &upper - &lower <= options.width,
        lower,
        upper,
        cells: cells.len(),
        witnesses,
    })
}
