//! Deterministic interval branch-and-bound on compact trimmed analytic surfaces.
//!
//! Lower bounds cover every retained parameter cell. Upper bounds use only exact
//! points admitted by the caller's trim predicate. Exhaustion returns the attained
//! enclosure, never a claim that the requested width was achieved.

use crate::{Point3, Surface, Vector3, math};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use std::{cmp::Ordering, collections::BinaryHeap};

const BITS: u32 = 128;
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn ball(x: &Rational) -> Ball {
    Ball::from_rational(x, BITS)
}
fn range(a: &Rational, b: &Rational) -> Ball {
    ball(a).hull(&ball(b))
}
fn norm(x: &[Ball; 3]) -> Ball {
    (&x[0].square() + &x[1].square() + x[2].square())
        .sqrt()
        .unwrap()
}
fn failure(detail: &str) -> Refusal {
    Refusal::new("shapes.distance-domain", detail)
}
fn point_box(p: &Point3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_budgeted_at(&p.x, BITS)?,
        dec.enclose_budgeted_at(&p.y, BITS)?,
        dec.enclose_budgeted_at(&p.z, BITS)?,
    ])
}
fn vector_box(v: &Vector3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        dec.enclose_budgeted_at(&v.x, BITS)?,
        dec.enclose_budgeted_at(&v.y, BITS)?,
        dec.enclose_budgeted_at(&v.z, BITS)?,
    ])
}

/// Refinement controls. No clock or approximate comparison determines traversal.
#[derive(Clone, Debug)]
pub struct Options {
    /// Requested absolute enclosure width, in millimetres.
    pub width: Rational,
    /// Maximum number of cells visited for trim classification.
    pub max_cells: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            width: Rational::from_ratio(1, 1_000_000).unwrap(),
            max_cells: 4096,
        }
    }
}
/// Proof supplied by a trim adapter for one connected chart cell.
#[derive(Clone, Copy, Debug)]
pub struct Membership {
    /// The exact chart centre belongs to the closed trimmed face.
    pub witness: bool,
    /// The entire cell is disjoint from the trimmed face. This requires proof;
    /// failure of the centre's membership test alone never suffices.
    pub excluded: bool,
}
/// Attained enclosure of the minimum unsigned Euclidean distance.
#[derive(Clone, Debug)]
pub struct Bounds {
    /// Proven nonnegative lower bound.
    pub lower: Rational,
    /// Proven upper bound, absent until an admitted witness pair is found.
    pub upper: Option<Rational>,
    /// Exact admitted witnesses corresponding to the upper bound.
    pub witnesses: Option<[Point3; 2]>,
    /// Whether the requested width was reached.
    pub converged: bool,
    /// Number of visited cells, including rejected cells.
    pub cells: usize,
}
/// A compact analytic carrier and a caller-certified enclosing world box.
/// The box must contain the entire trimmed face, not just its boundary or samples.
#[derive(Clone, Debug)]
pub struct Domain {
    surface: Surface,
    bounds: [Ball; 3],
    origin: Point3,
    frame: [Vector3; 3],
    origin_box: [Ball; 3],
    frame_box: [[Ball; 3]; 3],
    dimensions: [Ball; 2],
    charts: Vec<Chart>,
}
#[derive(Clone, Debug)]
struct Chart {
    quadrant: [usize; 2],
    range: [[Rational; 2]; 2],
}
fn circle(t: &Real, quadrant: usize) -> [Real; 2] {
    let t2 = t.square();
    let denominator = &Real::one() + &t2;
    let c = &(Real::one() - t2) / &denominator;
    let s = (&Real::from_i64(2) * t) / denominator;
    match quadrant {
        0 => [c, s],
        1 => [-s, c],
        2 => [-c, -s],
        _ => [s, -c],
    }
}
fn circle_box(t: &Ball, quadrant: usize) -> [Ball; 2] {
    let denominator = ball(&q(1)) + t.square();
    let c = (ball(&q(1)) - t.square())
        .checked_div(&denominator)
        .unwrap();
    let s = (&ball(&q(2)) * t).checked_div(&denominator).unwrap();
    match quadrant {
        0 => [c, s],
        1 => [-s, c],
        2 => [-c, -s],
        _ => [s, -c],
    }
}
impl Domain {
    /// Prepares all periodic charts and finite linear parameter ranges. Frames
    /// and dimensions may be algebraic; chart witnesses retain those definitions.
    pub fn new(surface: &Surface, bounds: [Ball; 3], dec: &mut Decider) -> Result<Self, Refusal> {
        let (origin, axis, reference) = surface.definition_frame();
        let normalize = |v: &Vector3| v.scaled(&(Real::one() / v.norm_squared().sqrt()));
        let u = normalize(reference);
        let axis = normalize(axis);
        let v = axis.cross(&u);
        let frame = [u, v, axis];
        let origin_box = point_box(origin, dec)?;
        let frame_box = [
            vector_box(&frame[0], dec)?,
            vector_box(&frame[1], dec)?,
            vector_box(&frame[2], dec)?,
        ];
        let dimensions = match surface {
            Surface::Plane(_) => [Real::zero(), Real::zero()],
            Surface::Cylinder(s) => [s.radius.clone(), Real::zero()],
            Surface::Cone(s) => [s.slope.clone(), Real::zero()],
            Surface::Sphere(s) => [s.radius.clone(), Real::zero()],
            Surface::Torus(s) => [s.major.clone(), s.minor.clone()],
            Surface::EllipticCylinder(s) => [s.major_sq.sqrt(), s.minor_sq.sqrt()],
        };
        let dimensions = [
            dec.enclose_budgeted_at(&dimensions[0], BITS)?,
            dec.enclose_budgeted_at(&dimensions[1], BITS)?,
        ];
        let displacement: [Ball; 3] = std::array::from_fn(|i| &bounds[i] - &origin_box[i]);
        let linear = |axis: usize| -> [Rational; 2] {
            let value = (0..3).fold(ball(&q(0)), |sum, i| {
                sum + &displacement[i] * &frame_box[axis][i]
            });
            [value.lo(), value.hi()]
        };
        let u = if matches!(surface, Surface::Plane(_)) {
            linear(0)
        } else {
            [q(0), q(1)]
        };
        let mut v = match surface {
            Surface::Sphere(_) => [q(-1), q(1)],
            Surface::Torus(_) => [q(0), q(1)],
            Surface::Plane(_) => linear(1),
            _ => linear(2),
        };
        if matches!(surface, Surface::Cone(_)) {
            v[0] = v[0].clone().max(q(0));
        }
        if u[0] > u[1] || v[0] > v[1] {
            return Err(failure("empty parameter cover"));
        }
        let mut charts = Vec::new();
        for i in 0..if matches!(surface, Surface::Plane(_)) {
            1
        } else {
            4
        } {
            for j in 0..if matches!(surface, Surface::Torus(_)) {
                4
            } else {
                1
            } {
                charts.push(Chart {
                    quadrant: [i, j],
                    range: [u.clone(), v.clone()],
                });
            }
        }
        Ok(Self {
            surface: surface.clone(),
            bounds,
            origin: origin.clone(),
            frame,
            origin_box,
            frame_box,
            dimensions,
            charts,
        })
    }
    /// Carrier retained by this domain.
    pub fn surface(&self) -> &Surface {
        &self.surface
    }
    fn point(&self, chart: &Chart) -> Point3 {
        let [u, v] = chart
            .range
            .each_ref()
            .map(|r| Real::from((&r[0] + &r[1]).mul_pow2(-1)));
        let local = if matches!(self.surface, Surface::Plane(_)) {
            [u, v, Real::zero()]
        } else {
            let [c, s] = circle(&u, chart.quadrant[0]);
            match &self.surface {
                Surface::Cylinder(a) => [&c * &a.radius, &s * &a.radius, v],
                Surface::EllipticCylinder(a) => {
                    [&c * &a.major_sq.sqrt(), &s * &a.minor_sq.sqrt(), v]
                }
                Surface::Cone(a) => [&(&c * &a.slope) * &v, &(&s * &a.slope) * &v, v],
                Surface::Sphere(a) => {
                    let [cv, sv] = circle(&v, 0);
                    [
                        &(&c * &cv) * &a.radius,
                        &(&s * &cv) * &a.radius,
                        &sv * &a.radius,
                    ]
                }
                Surface::Torus(a) => {
                    let [cv, sv] = circle(&v, chart.quadrant[1]);
                    let radius = &a.major + &(&a.minor * &cv);
                    [&c * &radius, &s * &radius, &sv * &a.minor]
                }
                _ => unreachable!(),
            }
        };
        self.origin
            .plus(&self.frame[0].scaled(&local[0]))
            .plus(&self.frame[1].scaled(&local[1]))
            .plus(&self.frame[2].scaled(&local[2]))
    }
    fn cover(&self, chart: &Chart, dec: &mut Decider) -> Result<Option<[Ball; 3]>, Refusal> {
        let [u, v] = chart.range.each_ref().map(|r| range(&r[0], &r[1]));
        dec.check_deadline()?;
        let local = if matches!(self.surface, Surface::Plane(_)) {
            [u, v, ball(&q(0))]
        } else {
            let [c, s] = circle_box(&u, chart.quadrant[0]);
            match &self.surface {
                Surface::Cylinder(_) => [&c * &self.dimensions[0], &s * &self.dimensions[0], v],
                Surface::EllipticCylinder(_) => {
                    [&c * &self.dimensions[0], &s * &self.dimensions[1], v]
                }
                Surface::Cone(_) => [
                    &(&c * &self.dimensions[0]) * &v,
                    &(&s * &self.dimensions[0]) * &v,
                    v,
                ],
                Surface::Sphere(_) => {
                    let [cv, sv] = circle_box(&v, 0);
                    [
                        &(&c * &cv) * &self.dimensions[0],
                        &(&s * &cv) * &self.dimensions[0],
                        &sv * &self.dimensions[0],
                    ]
                }
                Surface::Torus(_) => {
                    let [cv, sv] = circle_box(&v, chart.quadrant[1]);
                    let radius = &self.dimensions[0] + &(&self.dimensions[1] * &cv);
                    [&c * &radius, &s * &radius, &sv * &self.dimensions[1]]
                }
                _ => unreachable!(),
            }
        };
        let mut p = self.origin_box.clone();
        for (j, coefficient) in local.iter().enumerate() {
            for (i, coordinate) in p.iter_mut().enumerate() {
                *coordinate = &*coordinate + &(&self.frame_box[j][i] * coefficient);
            }
        }
        if p.iter()
            .zip(&self.bounds)
            .any(|(coordinate, bound)| coordinate.hi() < bound.lo() || bound.hi() < coordinate.lo())
        {
            return Ok(None);
        }
        Ok(Some(p))
    }
    fn clipped(&self, mut cover: [Ball; 3]) -> [Ball; 3] {
        for (coordinate, bound) in cover.iter_mut().zip(&self.bounds) {
            *coordinate = range(
                &coordinate.lo().max(bound.lo()),
                &coordinate.hi().min(bound.hi()),
            );
        }
        cover
    }
}

/// An exact point on the carrier near a target. This is a witness candidate,
/// not a closest-point proof, especially on an elliptic cylinder.
pub fn project(surface: &Surface, target: &Point3, dec: &mut Decider) -> Result<Point3, Refusal> {
    let (origin, axis, reference) = surface.definition_frame();
    let unit = |v: &Vector3| v.scaled(&(Real::one() / v.norm_squared().sqrt()));
    let axis = unit(axis);
    let d = target.minus(origin);
    let h = d.dot(&axis);
    if matches!(surface, Surface::Plane(_)) {
        return Ok(target.plus(&axis.scaled(&(-h))));
    }
    let planar = d.minus(&axis.scaled(&h));
    let rho2 = planar.norm_squared();
    let radial = if math::sign(&rho2, dec)? == Sign::Zero {
        unit(reference)
    } else {
        planar.scaled(&(Real::one() / rho2.sqrt()))
    };
    Ok(match surface {
        Surface::Cylinder(a) => origin
            .plus(&axis.scaled(&h))
            .plus(&radial.scaled(&a.radius)),
        Surface::Cone(a) => {
            let t = (&h + &(&a.slope * &rho2.sqrt())) / (Real::one() + a.slope.square());
            let t = if math::sign(&t, dec)? == Sign::Negative {
                Real::zero()
            } else {
                t
            };
            origin
                .plus(&axis.scaled(&t))
                .plus(&radial.scaled(&(&a.slope * &t)))
        }
        Surface::Sphere(a) => {
            let n = d.norm_squared();
            let direction = if math::sign(&n, dec)? == Sign::Zero {
                unit(reference)
            } else {
                d.scaled(&(Real::one() / n.sqrt()))
            };
            origin.plus(&direction.scaled(&a.radius))
        }
        Surface::Torus(a) => {
            let x = &rho2.sqrt() - &a.major;
            let n = x.square() + h.square();
            let [c, s] = if math::sign(&n, dec)? == Sign::Zero {
                [Real::one(), Real::zero()]
            } else {
                let length = n.sqrt();
                [&x / &length, &h / &length]
            };
            origin
                .plus(&radial.scaled(&(&a.major + &(&a.minor * &c))))
                .plus(&axis.scaled(&(&a.minor * &s)))
        }
        Surface::EllipticCylinder(a) => {
            let u = unit(reference);
            let v = axis.cross(&u);
            let x = d.dot(&u) / a.major_sq.sqrt();
            let y = d.dot(&v) / a.minor_sq.sqrt();
            let n = x.square() + y.square();
            let [c, s] = if math::sign(&n, dec)? == Sign::Zero {
                [Real::one(), Real::zero()]
            } else {
                let length = n.sqrt();
                [&x / &length, &y / &length]
            };
            origin
                .plus(&axis.scaled(&h))
                .plus(&u.scaled(&(a.major_sq.sqrt() * c)))
                .plus(&v.scaled(&(a.minor_sq.sqrt() * s)))
        }
        _ => unreachable!(),
    })
}

/// Trim classification callback used by the distance engine. The box encloses
/// the whole connected chart cell. A caller may conservatively never exclude it.
pub trait Trim {
    /// Tests the exact centre and optionally proves the entire cell absent.
    fn classify(
        &mut self,
        centre: &Point3,
        cover: &[Ball; 3],
        dec: &mut Decider,
    ) -> Result<Membership, Refusal>;
}
impl<F> Trim for F
where
    F: FnMut(&Point3, &[Ball; 3], &mut Decider) -> Result<Membership, Refusal>,
{
    fn classify(
        &mut self,
        p: &Point3,
        b: &[Ball; 3],
        d: &mut Decider,
    ) -> Result<Membership, Refusal> {
        self(p, b, d)
    }
}
#[derive(Clone)]
struct Node {
    charts: [Chart; 2],
    lower: Rational,
    witness_upper: Option<Rational>,
    sequence: usize,
}
impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.lower == other.lower && self.sequence == other.sequence
    }
}
impl Eq for Node {}
impl PartialOrd for Node {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Node {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .lower
            .cmp(&self.lower)
            .then_with(|| other.sequence.cmp(&self.sequence))
    }
}
fn separation(a: &[Ball; 3], b: &[Ball; 3]) -> Rational {
    norm(&std::array::from_fn(|i| {
        let gap = (a[i].lo() - b[i].hi()).max(b[i].lo() - a[i].hi()).max(q(0));
        ball(&gap)
    }))
    .lo()
    .max(q(0))
}
fn witness_distance(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<Rational, Refusal> {
    let d = a.minus(b);
    Ok(dec
        .enclose_budgeted_at(&d.norm_squared(), BITS)?
        .sqrt()
        .unwrap()
        .hi())
}
fn support_lower(
    surface: &Surface,
    p: &Point3,
    projected: &Point3,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let value = if let Surface::EllipticCylinder(s) = surface {
        let d = p.minus(&s.origin);
        let u = s
            .major_dir
            .scaled(&(Real::one() / s.major_dir.norm_squared().sqrt()));
        let axis = s.axis.scaled(&(Real::one() / s.axis.norm_squared().sqrt()));
        let v = axis.cross(&u);
        let x = d.dot(&u).square() / s.major_sq.clone();
        let y = d.dot(&v).square() / s.minor_sq.clone();
        let residual = (x + y).sqrt() - Real::one();
        // This residual is 1/minor-Lipschitz in the transverse plane.
        let residual = if math::sign(&residual, dec)? == Sign::Negative {
            -residual
        } else {
            residual
        };
        s.minor_sq.sqrt() * residual
    } else {
        p.minus(projected).norm_squared().sqrt()
    };
    Ok(dec.enclose_budgeted_at(&value, BITS)?.lo().max(q(0)))
}
struct Search<'a> {
    domains: [&'a Domain; 2],
    trims: [&'a mut dyn Trim; 2],
    fixed: Option<&'a Point3>,
    result: Bounds,
    queue: BinaryHeap<Node>,
    sequence: usize,
}
impl Search<'_> {
    fn visit(&mut self, charts: [Chart; 2], dec: &mut Decider) -> Result<(), Refusal> {
        dec.check_deadline()?;
        self.result.cells += 1;
        let mut covers = [
            std::array::from_fn(|_| ball(&q(0))),
            std::array::from_fn(|_| ball(&q(0))),
        ];
        let mut points = [Point3::origin(), Point3::origin()];
        let mut admitted = [false; 2];
        for i in 0..2 {
            if i == 0
                && let Some(p) = self.fixed
            {
                points[i] = p.clone();
                covers[i] = point_box(p, dec)?;
                admitted[i] = true;
                continue;
            }
            let Some(cover) = self.domains[i].cover(&charts[i], dec)? else {
                return Ok(());
            };
            points[i] = self.domains[i].point(&charts[i]);
            let membership = self.trims[i].classify(&points[i], &cover, dec)?;
            if membership.excluded {
                return Ok(());
            }
            admitted[i] = membership.witness;
            covers[i] = self.domains[i].clipped(cover);
        }
        let lower = separation(&covers[0], &covers[1]);
        let witness_upper = if admitted.into_iter().all(|x| x) {
            Some(self.accept(points, dec)?)
        } else {
            None
        };
        if self.result.upper.as_ref().is_none_or(|u| lower < *u) {
            self.queue.push(Node {
                charts,
                lower,
                witness_upper,
                sequence: self.sequence,
            });
            self.sequence += 1;
        }
        Ok(())
    }
    fn accept(&mut self, points: [Point3; 2], dec: &mut Decider) -> Result<Rational, Refusal> {
        let upper = witness_distance(&points[0], &points[1], dec)?;
        if self.result.upper.as_ref().is_none_or(|u| upper < *u) {
            self.result.upper = Some(upper.clone());
            self.result.witnesses = Some(points);
        }
        Ok(upper)
    }
    fn run(mut self, options: &Options, dec: &mut Decider) -> Result<Bounds, Refusal> {
        let mut iterations = 0usize;
        while let Some(node) = self.queue.peek() {
            let lower = self
                .result
                .upper
                .as_ref()
                .map_or_else(|| node.lower.clone(), |u| node.lower.clone().min(u.clone()));
            self.result.lower = lower.max(self.result.lower);
            if self
                .result
                .upper
                .as_ref()
                .is_some_and(|u| u - &self.result.lower <= options.width)
            {
                break;
            }
            if self.result.cells + 2 > options.max_cells {
                break;
            }
            iterations += 1;
            // Boundary-straddling cells may keep the best lower bound forever
            // while their centres lie outside the trim. Periodically refine the
            // best admitted cell too, so feasible witnesses can approach that
            // same closed boundary. This changes work order, never a proof bound.
            let feasible = if iterations.is_multiple_of(8) {
                self.queue
                    .iter()
                    .filter(|n| n.witness_upper.is_some())
                    .min_by(|a, b| {
                        a.witness_upper
                            .cmp(&b.witness_upper)
                            .then_with(|| a.sequence.cmp(&b.sequence))
                    })
                    .cloned()
            } else {
                None
            };
            let node = if let Some(node) = feasible {
                self.queue.retain(|n| n.sequence != node.sequence);
                node
            } else {
                self.queue.pop().unwrap()
            };
            // Split the largest world-space cover reduction, with fixed tie order.
            let mut choice = None;
            for i in if self.fixed.is_some() { 1 } else { 0 }..2 {
                for axis in 0..2 {
                    let span = &node.charts[i].range[axis][1] - &node.charts[i].range[axis][0];
                    if span.is_zero() {
                        continue;
                    }
                    let mut halves = node.charts.clone();
                    let parent = self.domains[i]
                        .cover(&node.charts[i], dec)?
                        .map_or(q(0), |b| b.iter().fold(q(0), |s, x| s + x.hi() - x.lo()));
                    halves[i].range[axis][1] = (&node.charts[i].range[axis][0]
                        + &node.charts[i].range[axis][1])
                        .mul_pow2(-1);
                    let left = self.domains[i]
                        .cover(&halves[i], dec)?
                        .map_or(q(0), |b| b.iter().fold(q(0), |s, x| s + x.hi() - x.lo()));
                    halves[i].range[axis][1] = node.charts[i].range[axis][1].clone();
                    halves[i].range[axis][0] = (&node.charts[i].range[axis][0]
                        + &node.charts[i].range[axis][1])
                        .mul_pow2(-1);
                    let right = self.domains[i]
                        .cover(&halves[i], dec)?
                        .map_or(q(0), |b| b.iter().fold(q(0), |s, x| s + x.hi() - x.lo()));
                    let reduction = parent - left.max(right);
                    if choice
                        .as_ref()
                        .is_none_or(|(_, _, best): &(usize, usize, Rational)| reduction > *best)
                    {
                        choice = Some((i, axis, reduction));
                    }
                }
            }
            let Some((i, axis, _)) = choice else {
                self.queue.push(node);
                break;
            };
            let mut a = node.charts.clone();
            let mut b = node.charts;
            let mid = (&a[i].range[axis][0] + &a[i].range[axis][1]).mul_pow2(-1);
            a[i].range[axis][1] = mid.clone();
            b[i].range[axis][0] = mid;
            self.visit(a, dec)?;
            self.visit(b, dec)?;
        }
        if self.queue.is_empty() {
            self.result.lower = self.result.upper.clone().unwrap_or(self.result.lower);
        }
        self.result.converged = self
            .result
            .upper
            .as_ref()
            .is_some_and(|u| u - &self.result.lower <= options.width);
        Ok(self.result)
    }
}
fn validate(options: &Options) -> Result<(), Refusal> {
    if options.width <= q(0) || options.max_cells == 0 {
        Err(failure("positive width and cell count required"))
    } else {
        Ok(())
    }
}
fn empty_bounds() -> Bounds {
    Bounds {
        lower: q(0),
        upper: None,
        witnesses: None,
        converged: false,
        cells: 0,
    }
}
/// Bounds distance from a point to a closed trimmed face. Projected and nearby
/// exact carrier witnesses are tried before interval subdivision. Every witness
/// is trim-checked; a support hit alone can never establish an upper bound.
pub fn point_distance(
    domain: &Domain,
    point: &Point3,
    options: &Options,
    trim: &mut dyn Trim,
    dec: &mut Decider,
) -> Result<Bounds, Refusal> {
    validate(options)?;
    let mut ignored = |_: &Point3, _: &[Ball; 3], _: &mut Decider| {
        Ok(Membership {
            witness: true,
            excluded: false,
        })
    };
    let mut search = Search {
        domains: [domain, domain],
        trims: [&mut ignored, trim],
        fixed: Some(point),
        result: empty_bounds(),
        queue: BinaryHeap::new(),
        sequence: 0,
    };
    let step = options.width.mul_pow2(-3);
    for offset in 0..28 {
        let target = if offset == 0 {
            point.clone()
        } else {
            let n = offset - 1;
            point.plus(&Vector3::new(
                Real::from(&q(n % 3 - 1) * &step),
                Real::from(&q(n / 3 % 3 - 1) * &step),
                Real::from(&q(n / 9 - 1) * &step),
            ))
        };
        let p = project(&domain.surface, &target, dec)?;
        if offset == 0 {
            search.result.lower = support_lower(&domain.surface, point, &p, dec)?;
        }
        let cover = point_box(&p, dec)?;
        if search.trims[1].classify(&p, &cover, dec)?.witness {
            search.accept([point.clone(), p], dec)?;
        }
        if search
            .result
            .upper
            .as_ref()
            .is_some_and(|u| u - &search.result.lower <= options.width)
        {
            search.result.converged = true;
            return Ok(search.result);
        }
    }
    if options.max_cells < domain.charts.len() {
        return Ok(search.result);
    }
    for chart in &domain.charts {
        search.visit([chart.clone(), chart.clone()], dec)?;
    }
    search.run(options, dec)
}
/// Bounds the minimum distance between two compact trimmed faces. This is a
/// minimum-distance query, not a Hausdorff or whole-face deviation measurement.
pub fn face_distance(
    a: &Domain,
    b: &Domain,
    options: &Options,
    trim_a: &mut dyn Trim,
    trim_b: &mut dyn Trim,
    dec: &mut Decider,
) -> Result<Bounds, Refusal> {
    validate(options)?;
    let mut search = Search {
        domains: [a, b],
        trims: [trim_a, trim_b],
        fixed: None,
        result: empty_bounds(),
        queue: BinaryHeap::new(),
        sequence: 0,
    };
    if options.max_cells < a.charts.len() * b.charts.len() {
        return Ok(search.result);
    }
    for ca in &a.charts {
        for cb in &b.charts {
            search.visit([ca.clone(), cb.clone()], dec)?;
        }
    }
    search.run(options, dec)
}

#[cfg(test)]
mod tests;
