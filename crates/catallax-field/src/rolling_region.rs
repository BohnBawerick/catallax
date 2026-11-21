//! Exact rolling-ball path pieces, independently of fitted B-spline carriers.
//!
//! [`path`] holds the one general routine for a path piece: membership of its
//! region and a guaranteed distance range to its rolling-ball face. [`Region`]
//! is the P5 rim-contact piece; [`elementary`] holds the line and circle pieces
//! of plane-plane and plane-cylinder fillets, derived from their supports.
//! [`revolved`] holds a plane with a cylinder, cone or sphere in any pose and
//! [`generator`] the straight paths between supports parallel to one line
//! (4.2a step 2). [`coaxial`] holds two curved supports about one axis and
//! [`saddle`] a cylinder with a second cylinder or a sphere off its axis
//! (4.2a step 3). [`rim`] holds a ring whose plane contact crosses circular
//! mouths, the field's own P5 rim contact (4.2a step 4). [`piece`] compiles
//! any of them as a field leaf.
//!
//! Only [`piece`] admits recipe fillets; the other building blocks are not
//! admissions. See `docs/M18-ROLLING-REGION.md` for the rim region's explicit set and its
//! remaining incidence obligations, and `docs/M42A-FIELD-FILLETS.md` for the
//! general routine.

use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, OperationAllowance, Rational, Real, Refusal, Sign};
use catallax_shapes::{Point3, RollingEnvelopeDefinition};

pub mod coaxial;
pub mod elementary;
pub mod generator;
pub mod path;
pub mod piece;
pub mod revolved;
pub mod rim;
pub mod saddle;
mod sectors;
pub use path::{Answer, Range};
pub use sectors::{Sector, Sectors, Transition};

// Local refusal limits of the M18 rolling region (docs/M18-ROLLING-REGION.md):
// 96-bit outward rounding, at most 4096 evaluated panels and depth 64. They are
// limits of this building block, not changes to any allowance or accuracy.
const BITS: u32 = 96;
const PANELS: usize = 4096;
const DEPTH: usize = 64;
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn b(n: i64) -> Ball {
    Ball::from_rational(&q(n), BITS)
}
fn interval(lo: &Rational, hi: &Rational) -> Ball {
    Ball::from_rational(lo, BITS).hull(&Ball::from_rational(hi, BITS))
}
fn gap(detail: &str) -> Refusal {
    Refusal::new("field.rolling-region-unresolved", detail)
}
/// The box of a rational box's image under an exact rigid motion: the box
/// of its moved corners, each widened to its certified binary64 bounds. A
/// corner without finite bounds widens the box to one no region reaches past.
pub(crate) fn moved_box(
    pose: &catallax_shapes::Isometry,
    b: &[[Rational; 2]; 3],
) -> [[Rational; 2]; 3] {
    let far = Rational::from_i64(1 << 40);
    let mut moved: Option<[[Rational; 2]; 3]> = None;
    for corner in 0..8 {
        let [x, y, z] = std::array::from_fn(|i| b[i][(corner >> i) & 1].clone());
        let image = pose.apply_point(&catallax_shapes::RPoint3::new(x, y, z).to_real());
        let c = [image.x, image.y, image.z].map(|v| {
            if let Some(v) = v.as_rational() {
                return [v.clone(), v];
            }
            let e = v.approx();
            match (
                Rational::from_f64_exact(e.lo()),
                Rational::from_f64_exact(e.hi()),
            ) {
                (Some(lo), Some(hi)) => [lo, hi],
                _ => [-&far, far.clone()],
            }
        });
        moved = Some(match moved {
            None => c,
            Some(m) => std::array::from_fn(|i| {
                [
                    m[i][0].clone().min(c[i][0].clone()),
                    m[i][1].clone().max(c[i][1].clone()),
                ]
            }),
        });
    }
    moved.expect("eight corners")
}
fn div(a: &Ball, b: &Ball) -> Result<Ball, Refusal> {
    a.checked_div(b)
        .ok_or_else(|| gap("rolling chart divisor is unresolved"))
}

/// Membership in a path piece's explicitly defined region. Paths without a
/// closed-form separation return a typed refusal at boundaries and
/// insufficiently separated cases, never `Boundary`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Relation {
    /// Strictly within the defined region.
    Interior,
    /// Strictly outside the defined region.
    Exterior,
    /// Exactly on the region's boundary; only closed-form paths prove this.
    Boundary,
}

/// Exact contacts at one station of the finite centre curve.
#[derive(Clone, Debug)]
pub struct Station {
    /// Centre of the rolling ball.
    pub centre: Point3,
    /// Contact with the body cylinder.
    pub cylinder: Point3,
    /// Contact with the circular mouth edge in its plane.
    pub mouth: Point3,
}

/// Lower envelope of the finite family of balls, clipped by the mouth sector,
/// the body/offset cylinders and the plane/offset plane. All lengths are mm.
/// This set is distinct from a complete torus and from a fitted carrier.
#[derive(Clone, Debug)]
pub struct Region {
    definition: RollingEnvelopeDefinition,
    centre: Real,
    cosine: Real,
    sine: Real,
    half_tangent: Real,
}
impl Region {
    /// Uses the validated frame and rim normal-cone regime of the definition.
    pub fn new(definition: RollingEnvelopeDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        let centre = definition.body() + definition.radius();
        let cosine = (&(&centre.square() + &definition.offset().square())
            - &definition.mouth().square())
            / (&(&Real::from_i64(2) * definition.offset()) * &centre);
        let sine = (&Real::one() - &cosine.square()).sqrt();
        if sign(&cosine, dec)? != Sign::Positive || sign(&sine, dec)? != Sign::Positive {
            return Err(gap(
                "finite mouth sector must be strictly acute and nonempty",
            ));
        }
        let half_tangent = &sine / &(&Real::one() + &cosine);
        Ok(Self {
            definition,
            centre,
            cosine,
            sine,
            half_tangent,
        })
    }
    /// Exact source definition. No fitting tolerance is introduced.
    pub fn definition(&self) -> &RollingEnvelopeDefinition {
        &self.definition
    }
    fn world(&self, [x, y, z]: [Real; 3]) -> Point3 {
        let d = &self.definition;
        d.origin().plus(
            &d.toward()
                .scaled(&x)
                .plus(&d.axis().scaled(&y))
                .plus(&d.across().scaled(&z)),
        )
    }
    /// Normalized half-angle station in [-1,1], including both transitions.
    pub fn station(&self, u: &Rational, dec: &mut Decider) -> Result<Station, Refusal> {
        dec.charge_work(1)?;
        if u < &q(-1) || u > &q(1) {
            return Err(gap("station is outside the finite mouth sector"));
        }
        let d = &self.definition;
        let t = &self.half_tangent * &Real::from(u.clone());
        let denominator = &Real::one() + &t.square();
        let cosine = (&Real::one() - &t.square()) / denominator.clone();
        let sine = (&Real::from_i64(2) * &t) / denominator;
        let x = &self.centre * &cosine;
        let z = &self.centre * &sine;
        let dx = &x - d.offset();
        let distance = (&dx.square() + &z.square()).sqrt();
        let height = (d.radius().square() - (d.mouth() - &distance).square()).sqrt();
        let scale = d.mouth() / &distance;
        Ok(Station {
            centre: self.world([x, height.clone(), z.clone()]),
            cylinder: self.world([d.body() * &cosine, height, d.body() * &sine]),
            mouth: self.world([d.offset() + &(&scale * &dx), Real::zero(), scale * z]),
        })
    }
    /// Strict membership with a complete interval cover of the centre curve.
    /// A midpoint can prove exclusion, never inclusion. No distance or named
    /// face certificate is implied by this result.
    pub fn classify(&self, point: &Point3, dec: &mut Decider) -> Result<Relation, Refusal> {
        path::classify(self, point, dec)
    }
    /// Guaranteed range for the distance from `point` to the rim envelope face,
    /// refined until it is at most `width` wide or the cover's limits are met.
    pub fn distance(
        &self,
        point: &Point3,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<Range, Refusal> {
        path::distance(self, point, width, dec)
    }
}

impl path::Path for Region {
    type Cover = Constants;
    fn local(&self, point: &Point3) -> [Real; 3] {
        let d = &self.definition;
        let delta = point.minus(d.origin());
        [
            delta.dot(d.toward()),
            delta.dot(d.axis()),
            delta.dot(d.across()),
        ]
    }
    fn walls(&self, [x, y, z]: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let d = &self.definition;
        let radial = &x.square() + &z.square();
        Ok(vec![
            y.clone(),
            d.radius() - y,
            &radial - &d.body().square(),
            &self.centre.square() - &radial,
            x.clone(),
            &(x * &self.sine) - &(z * &self.cosine),
            &(x * &self.sine) + &(z * &self.cosine),
        ])
    }
    fn cover(&self, dec: &mut Decider) -> Result<Constants, Refusal> {
        Constants::new(self, dec)
    }
    fn frames(&self, dec: &mut Decider) -> Result<Constants, Refusal> {
        let mut constants = Constants::new(self, dec)?;
        constants.body = Some(enclose(self.definition.body(), BITS, dec)?);
        Ok(constants)
    }
}
/// Interval enclosures of the rim chart, prepared under the caller's decider.
#[derive(Clone, Debug)]
pub struct Constants {
    c: Ball,
    o: Ball,
    m: Ball,
    r: Ball,
    t: Ball,
    // Body radius; enclosed only for distance queries, so membership work is
    // unchanged.
    body: Option<Ball>,
}
impl Constants {
    fn new(region: &Region, dec: &mut Decider) -> Result<Self, Refusal> {
        let d = &region.definition;
        Ok(Self {
            c: enclose(&region.centre, BITS, dec)?,
            o: enclose(d.offset(), BITS, dec)?,
            m: enclose(d.mouth(), BITS, dec)?,
            r: enclose(d.radius(), BITS, dec)?,
            t: enclose(&region.half_tangent, BITS, dec)?,
            body: None,
        })
    }
    // Angle, mouth distance D, mouth gap m-D and height H over the panel.
    fn chart(&self, u: &Ball) -> Result<[Ball; 6], Refusal> {
        let t = &self.t * u;
        let denominator = &b(1) + &t.square();
        let cosine = div(&(&b(1) - &t.square()), &denominator)?;
        let sine = div(&(&b(2) * &t), &denominator)?;
        let x = &self.c * &cosine;
        let z = &self.c * &sine;
        let distance = ((&x - &self.o).square() + z.square())
            .sqrt()
            .ok_or_else(|| gap("undefined mouth distance"))?;
        // The validated finite chart has 0 <= m-D <= m-(o-C).
        // Intersecting these enclosing ranges only tightens a known invariant.
        let gap_ball = &self.m - &distance;
        let upper = (&self.m - &(&self.o - &self.c)).hi();
        let lo = gap_ball.lo().max(q(0));
        let hi = gap_ball.hi().min(upper);
        if lo > hi {
            return Err(gap("mouth-gap enclosure is inconsistent"));
        }
        let gap_ball = interval(&lo, &hi);
        let height = (&self.r.square() - &gap_ball.square())
            .sqrt()
            .ok_or_else(|| gap("undefined rolling height"))?;
        Ok([cosine, sine, distance, gap_ball, height, x])
    }
}
impl path::Cover for Constants {
    fn radius(&self) -> &Ball {
        &self.r
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        let [_, sine, _, _, height, x] = self.chart(u)?;
        Ok([x, height, &self.c * &sine])
    }
    /// With `theta` the centre angle, `dD/dtheta = C o sin / D` and
    /// `dH/dtheta = (m-D) dD/dtheta / H`. The tangent below is that derivative
    /// times the positive factor `D H`, so it needs no division.
    fn frame(&self, u: &Ball) -> Result<path::Frame, Refusal> {
        let body = self
            .body
            .as_ref()
            .ok_or_else(|| gap("rim frames need the body radius enclosure"))?;
        let [cosine, sine, distance, gap_ball, height, x] = self.chart(u)?;
        let z = &self.c * &sine;
        let scale = &distance * &height;
        let tangent = [
            -&(&z * &scale),
            &(&(&gap_ball * &self.c) * &self.o) * &sine,
            &x * &scale,
        ];
        let ratio = div(&self.m, &distance)?;
        let mouth = [&self.o + &(&ratio * &(&x - &self.o)), b(0), &ratio * &z];
        let cylinder = [body * &cosine, height.clone(), body * &sine];
        Ok(path::Frame {
            centre: [x, height, z],
            tangent,
            feet: [cylinder, mouth],
            side: None,
        })
    }
    fn excludes(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = path::Cover::centre(self, u)?;
        let xy = (&p[0] - &c[0]).square() + (&p[2] - &c[2]).square();
        let all = &xy + &(&p[1] - &c[1]).square();
        Ok(xy.hi() < self.r.square().lo()
            && (p[1].lo() > c[1].hi() || all.hi() < self.r.square().lo()))
    }
    fn below(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = path::Cover::centre(self, u)?;
        let xy = (&p[0] - &c[0]).square() + (&p[2] - &c[2]).square();
        let all = &xy + &(&p[1] - &c[1]).square();
        Ok(xy.lo() > self.r.square().hi()
            || (p[1].hi() < c[1].lo() && all.lo() > self.r.square().hi()))
    }
}

/// Complete per-point semantic output, including work and full refusal evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct Record<T = Relation> {
    /// Exact set membership (or another per-point answer) or typed refusal.
    pub result: Result<T, Refusal>,
    /// Actual numeric work in this point's transferred allowance.
    pub work: u64,
    /// Whether that allowance was exhausted.
    pub exhausted: bool,
}
/// Ordered results and number of nonempty shared-pool assignments.
#[derive(Clone, Debug, PartialEq)]
pub struct Batch<T = Relation> {
    /// One complete record per input, in input order.
    pub records: Vec<Record<T>>,
    /// Number of nonempty assignments admitted to the shared field pool.
    pub workers: usize,
}
impl Region {
    /// Production batch service. Defaults to all available cores, capped by
    /// input count. Splits one caller-owned allowance by input index, independent
    /// of worker count. No cache or unused allowance is transferred between jobs.
    /// Propagates the caller's numeric policy and absolute deadline to every
    /// worker, and checks cancellation before and after the complete batch.
    pub fn classify_batch(
        &self,
        points: &[Point3],
        workers: Option<usize>,
        allowance: OperationAllowance,
        dec: &mut Decider,
    ) -> Result<Batch, Refusal> {
        batch(
            points,
            workers,
            allowance,
            dec,
            "discarded membership after allowance exhaustion",
            |point, child| self.classify(point, child),
        )
    }
}

/// Batch membership plus distance range for any path piece, with the same
/// index-based allowance split and ordered collection as
/// [`Region::classify_batch`].
pub fn query_batch<P: path::Path + Sync>(
    path: &P,
    points: &[Point3],
    width: &Rational,
    workers: Option<usize>,
    allowance: OperationAllowance,
    dec: &mut Decider,
) -> Result<Batch<Answer>, Refusal> {
    batch(
        points,
        workers,
        allowance,
        dec,
        "discarded answer after allowance exhaustion",
        |point, child| path::query(path, point, width, child),
    )
}

/// Batch distance ranges for any path piece, with the same split and order.
/// A point on the face needs no membership answer, which an interval cover
/// cannot settle there.
pub fn distance_batch<P: path::Path + Sync>(
    path: &P,
    points: &[Point3],
    width: &Rational,
    workers: Option<usize>,
    allowance: OperationAllowance,
    dec: &mut Decider,
) -> Result<Batch<Range>, Refusal> {
    batch(
        points,
        workers,
        allowance,
        dec,
        "discarded range after allowance exhaustion",
        |point, child| path::distance(path, point, width, child),
    )
}

fn batch<T: Send, F>(
    points: &[Point3],
    workers: Option<usize>,
    mut allowance: OperationAllowance,
    dec: &mut Decider,
    discarded: &str,
    run: F,
) -> Result<Batch<T>, Refusal>
where
    F: Fn(&Point3, &mut Decider) -> Result<T, Refusal> + Sync,
{
    dec.check_deadline()?;
    let count = workers.unwrap_or_else(crate::workers::worker_count);
    if count == 0 {
        return Err(gap("worker count must be positive"));
    }
    if points.is_empty() {
        return Ok(Batch {
            records: vec![],
            workers: 0,
        });
    }
    let count = count.min(crate::workers::worker_count()).min(points.len());
    let policy = *dec.policy();
    let deadline = dec.deadline();
    let share = allowance.limit() / points.len() as u64;
    let remainder = allowance.limit() % points.len() as u64;
    let mut partitions: Vec<Vec<_>> = (0..count).map(|_| vec![]).collect();
    for (index, point) in points.iter().enumerate() {
        let token = allowance
            .split_off(share + u64::from((index as u64) < remainder))
            .expect("partitioned allowance");
        partitions[index % count].push((index, point, token));
    }
    let partitions = partitions
        .into_iter()
        .map(|jobs| std::sync::Mutex::new(Some(jobs)))
        .collect::<Vec<_>>();
    let mut records = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        crate::workers::with_workers(count, || {
            crate::workers::map(&partitions, |jobs| {
                jobs.lock()
                    .unwrap()
                    .take()
                    .expect("single-use rolling-region jobs")
                    .into_iter()
                    .map(|(index, point, token)| {
                        let (mut result, work) = token.run(|| {
                            let mut child = Decider::new(policy);
                            child.set_deadline(deadline);
                            run(point, &mut child)
                        });
                        if work.exhausted && result.is_ok() {
                            result = Err(gap(discarded));
                        }
                        (
                            index,
                            Record {
                                result,
                                work: work.used,
                                exhausted: work.exhausted,
                            },
                        )
                    })
                    .collect::<Vec<_>>()
            })
        })
    }))
    .map_err(|_| gap("rolling-region worker panicked"))??
    .into_iter()
    .flatten()
    .collect::<Vec<_>>();
    dec.check_deadline()?;
    records.sort_by_key(|row| row.0);
    Ok(Batch {
        records: records.into_iter().map(|row| row.1).collect(),
        workers: count,
    })
}
