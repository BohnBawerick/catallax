//! Exact lower rolling-ball regions, independently of fitted B-spline carriers.
//!
//! This opt-in geometric building block is NOT a recipe fillet admission. See
//! `docs/M18-ROLLING-REGION.md` for the explicit set and the remaining incidence
//! obligations before it can replace a P5 recipe region.

use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, OperationAllowance, Rational, Real, Refusal, Sign};
use catallax_shapes::{Point3, RollingEnvelopeDefinition};

mod sectors;
pub use sectors::{Sector, Sectors, Transition};

const BITS: u32 = 96;
const PANELS: usize = 4096;
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
fn div(a: &Ball, b: &Ball) -> Result<Ball, Refusal> {
    a.checked_div(b)
        .ok_or_else(|| gap("rolling chart divisor is unresolved"))
}

/// Strict membership in the explicitly defined lower-envelope region.
/// Boundary and insufficiently separated cases return a typed refusal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Relation {
    /// Strictly within the defined clipped lower region.
    Interior,
    /// Strictly outside the defined region.
    Exterior,
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
        dec.charge_work(1)?;
        let d = &self.definition;
        let delta = point.minus(d.origin());
        let local = [
            delta.dot(d.toward()),
            delta.dot(d.axis()),
            delta.dot(d.across()),
        ];
        let [x, y, z] = &local;
        let radial = &x.square() + &z.square();
        let constraints = [
            y.clone(),
            d.radius() - y,
            &radial - &d.body().square(),
            &self.centre.square() - &radial,
            x.clone(),
            &(x * &self.sine) - &(z * &self.cosine),
            &(x * &self.sine) + &(z * &self.cosine),
        ];
        let mut wall = false;
        for c in constraints {
            match sign(&c, dec)? {
                Sign::Negative => return Ok(Relation::Exterior),
                Sign::Zero => wall = true,
                Sign::Positive => {}
            }
        }
        let point = [
            enclose(x, BITS, dec)?,
            enclose(y, BITS, dec)?,
            enclose(z, BITS, dec)?,
        ];
        let constants = Constants::new(self, dec)?;
        let mut pending = vec![(q(-1), q(1), 0)];
        let mut panels = 0;
        while let Some((lo, hi, depth)) = pending.pop() {
            dec.charge_work(1)?;
            panels += 1;
            if panels > PANELS {
                return Err(gap("finite centre-curve cover did not separate the query"));
            }
            let mid = (&lo + &hi) / q(2);
            // A whole-ball witness is a sound exclusion from the lower region.
            if constants.excludes(&point, &interval(&mid, &mid))? {
                return Ok(Relation::Exterior);
            }
            if constants.below(&point, &interval(&lo, &hi))? {
                continue;
            }
            if depth == 64 {
                return Err(gap("finite centre-curve subdivision depth exhausted"));
            }
            pending.push((mid.clone(), hi, depth + 1));
            pending.push((lo, mid, depth + 1));
        }
        if wall {
            Err(gap(
                "membership touches a clipping wall; incidence is not certified",
            ))
        } else {
            Ok(Relation::Interior)
        }
    }
}
struct Constants {
    c: Ball,
    o: Ball,
    m: Ball,
    r: Ball,
    t: Ball,
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
        })
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        let t = &self.t * u;
        let denominator = &b(1) + &t.square();
        let x = &self.c * &div(&(&b(1) - &t.square()), &denominator)?;
        let z = &self.c * &div(&(&b(2) * &t), &denominator)?;
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
        Ok([x, height, z])
    }
    fn excludes(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = self.centre(u)?;
        let xy = (&p[0] - &c[0]).square() + (&p[2] - &c[2]).square();
        let all = &xy + &(&p[1] - &c[1]).square();
        Ok(xy.hi() < self.r.square().lo()
            && (p[1].lo() > c[1].hi() || all.hi() < self.r.square().lo()))
    }
    fn below(&self, p: &[Ball; 3], u: &Ball) -> Result<bool, Refusal> {
        let c = self.centre(u)?;
        let xy = (&p[0] - &c[0]).square() + (&p[2] - &c[2]).square();
        let all = &xy + &(&p[1] - &c[1]).square();
        Ok(xy.lo() > self.r.square().hi()
            || (p[1].hi() < c[1].lo() && all.lo() > self.r.square().hi()))
    }
}

/// Complete per-point semantic output, including work and full refusal evidence.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    /// Exact set membership or typed refusal.
    pub result: Result<Relation, Refusal>,
    /// Actual numeric work in this point's transferred allowance.
    pub work: u64,
    /// Whether that allowance was exhausted.
    pub exhausted: bool,
}
/// Ordered results and actual number of nonempty scoped workers.
#[derive(Clone, Debug, PartialEq)]
pub struct Batch {
    /// One complete record per input, in input order.
    pub records: Vec<Record>,
    /// Number of nonempty workers actually started.
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
        mut allowance: OperationAllowance,
        dec: &mut Decider,
    ) -> Result<Batch, Refusal> {
        dec.check_deadline()?;
        let count =
            workers.unwrap_or_else(|| std::thread::available_parallelism().map_or(1, usize::from));
        if count == 0 {
            return Err(gap("worker count must be positive"));
        }
        if points.is_empty() {
            return Ok(Batch {
                records: vec![],
                workers: 0,
            });
        }
        let count = count.min(points.len());
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
        let mut records = std::thread::scope(|scope| {
            let mut handles = vec![];
            let mut failure = None;
            for jobs in partitions {
                match std::thread::Builder::new().spawn_scoped(scope, move || {
                    jobs.into_iter()
                        .map(|(index, point, token)| {
                            let (mut result, work) = token.run(|| {
                                let mut child = Decider::new(policy);
                                child.set_deadline(deadline);
                                self.classify(point, &mut child)
                            });
                            if work.exhausted && result.is_ok() {
                                result =
                                    Err(gap("discarded membership after allowance exhaustion"));
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
                }) {
                    Ok(h) => handles.push(h),
                    Err(_) => {
                        failure = Some(gap("could not start rolling-region worker"));
                        break;
                    }
                }
            }
            let mut out = vec![];
            for handle in handles {
                match handle.join() {
                    Ok(rows) => out.extend(rows),
                    Err(_) => failure = Some(gap("rolling-region worker panicked")),
                }
            }
            failure.map_or(Ok(out), Err)
        })?;
        dec.check_deadline()?;
        records.sort_by_key(|row| row.0);
        Ok(Batch {
            records: records.into_iter().map(|row| row.1).collect(),
            workers: count,
        })
    }
}
