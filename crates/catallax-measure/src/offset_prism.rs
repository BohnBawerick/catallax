//! Normal-offset cross-section integrals, independent of the fitted carrier.
//!
//! Reflection symmetry reduces Green's area/second-moment integrals and arc
//! length to one principal quadrant. The existing validated Taylor integrator
//! supplies the whole-panel remainder; no meshed or fitted values enter here.

use catallax_intersect::{TracePiece, trace_elliptic_offset_quarter};
use catallax_numbers::{Answer, Answered, Ball, Decider, MarkLog, Outcome, Rational, Real, Refusal, StrictSign};
use catallax_shapes::{EllipticCylinder, EllipticOffsetDefinition, OffsetSide, Point3, Surface, Vector3};

use crate::{Accuracy, Enclosure, MassProperties, VolumeAreaCentroid};
use crate::integrate::{Builder, Program, Stop, Work};

/// Exact-target volume, area, centroid and centroidal unit-density inertia of
/// the product of an elliptic normal-offset profile and an axial interval.
/// Coordinates use the original cylinder's unit axis. An empty interval denotes
/// an empty solid. This does not measure the fitted spline carrier or validate a
/// body certificate; a caller binding this definition to a body must check that
/// certificate separately. Nonregular inward profiles refuse.
pub fn elliptic_offset_prism_mass_properties(
    definition: &EllipticOffsetDefinition,
    axial: [&Real; 2],
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<MassProperties> {
    measure(definition, axial, accuracy, true, dec)
}

/// The three original exact-target properties, without computing inertia.
/// See [`elliptic_offset_prism_mass_properties`] for the definition and scope.
pub fn elliptic_offset_prism_volume_area_centroid(
    definition: &EllipticOffsetDefinition,
    axial: [&Real; 2],
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<VolumeAreaCentroid> {
    let answer = measure(definition, axial, accuracy, false, dec)?;
    Ok(Answer::with_marks(VolumeAreaCentroid {
        volume: answer.value.volume, area: answer.value.area, centroid: answer.value.centroid,
    }, answer.marks))
}

fn program(inertia: bool) -> Program {
    let b = Builder::default();
    let x = b.input(0);
    let dx = b.input(3);
    let dy = b.input(4);
    let area = x * dy;
    let perimeter = (dx * dx + dy * dy).sqrt();
    if inertia {
        let y = b.input(1);
        b.finish(&[area, perimeter, x * x * x * dy / b.k(3), -y * y * y * dx / b.k(3)])
    } else {
        b.finish(&[area, perimeter])
    }
}

fn integrate(
    program: &Program, trace: &TracePiece, depth: u32, work: &mut Work,
) -> Result<Vec<Ball>, Stop> {
    let count = 1u64.checked_shl(depth).ok_or(Stop::Budget)?;
    let step = Rational::one().mul_pow2(-i64::from(depth));
    let mut sum: Option<Vec<Ball>> = None;
    for index in 0..count {
        let lo = &step * &Rational::from_integer(catallax_numbers::Integer::from_u64(index));
        let hi = &lo + &step;
        let values = program.panel(trace, &lo, &hi, work)?;
        if let Some(sum) = &mut sum {
            for (sum, value) in sum.iter_mut().zip(values) { *sum = work.add(sum, &value)?; }
        } else {
            sum = Some(values);
        }
    }
    let four = work.rational(&Rational::from_i64(4))?;
    let mut sum = sum.unwrap();
    for value in &mut sum { *value = work.mul(value, &four)?; }
    Ok(sum)
}

fn real3(values: &[Real; 3], work: &mut Work) -> Result<[Ball; 3], Stop> {
    Ok([work.real(&values[0])?, work.real(&values[1])?, work.real(&values[2])?])
}

fn inertia_entries(tensor: &crate::InertiaTensor) -> [&Enclosure; 6] {
    [&tensor[0][0], &tensor[0][1], &tensor[0][2], &tensor[1][1], &tensor[1][2], &tensor[2][2]]
}

struct Frame {
    axes: Option<[[Real; 3]; 3]>,
    center: [Real; 3],
    height: Real,
}

impl Frame {
    fn new(original: &EllipticCylinder, axial: [&Real; 2], inertia: bool) -> Self {
        let unit = |v: &Vector3| v.scaled(&(&Real::one() / &v.norm_squared().sqrt()));
        let axis = unit(&original.axis);
        let midpoint = (axial[0] + axial[1]) / Real::from_i64(2);
        let center = original.origin.plus(&axis.scaled(&midpoint));
        let axes = inertia.then(|| {
            let major = unit(&original.major_dir);
            let minor = axis.cross(&major);
            [major, minor, axis].map(|v| [v.x, v.y, v.z])
        });
        Self {
            axes,
            center: [center.x, center.y, center.z], height: axial[1] - axial[0],
        }
    }

    fn properties(&self, values: &[Ball], work: &mut Work) -> Result<MassProperties, Stop> {
        let height = work.real(&self.height)?;
        let volume = work.mul(&values[0], &height)?;
        let lateral = work.mul(&values[1], &height)?;
        let two = work.rational(&Rational::from_i64(2))?;
        let caps = work.mul(&values[0], &two)?;
        let area = work.add(&caps, &lateral)?;
        let center = real3(&self.center, work)?;
        let inertia = if values.len() == 4 {
            let mu = work.mul(&values[2], &height)?;
            let mv = work.mul(&values[3], &height)?;
            let height2 = work.mul(&height, &height)?;
            let axial = work.mul(&volume, &height2)?;
            let twelve = work.rational(&Rational::from_i64(12))?;
            let mw = work.div(&axial, &twelve)?;
            let principal = [work.add(&mv, &mw)?, work.add(&mu, &mw)?, work.add(&mu, &mv)?];
            let axes = self.axes.as_ref().unwrap();
            let basis = [real3(&axes[0], work)?, real3(&axes[1], work)?, real3(&axes[2], work)?];
            let zero = Enclosure::exact(Rational::zero());
            let mut tensor: [[Enclosure; 3]; 3] = std::array::from_fn(|_| std::array::from_fn(|_| zero.clone()));
            for i in 0..3 {
                for j in i..3 {
                    let mut sum = work.rational(&Rational::zero())?;
                    for (axis, moment) in basis.iter().zip(&principal) {
                        let coefficient = work.mul(&axis[i], &axis[j])?;
                        let value = work.mul(&coefficient, moment)?;
                        sum = work.add(&sum, &value)?;
                    }
                    tensor[i][j] = Enclosure::from_ball(&sum);
                    if i != j { tensor[j][i] = tensor[i][j].clone(); }
                }
            }
            Some(tensor)
        } else { None };
        Ok(MassProperties {
            volume: Enclosure::from_ball(&volume), area: Enclosure::from_ball(&area),
            centroid: Some(center.each_ref().map(Enclosure::from_ball)), inertia,
        })
    }
}

fn reached(result: &MassProperties, accuracy: &Accuracy) -> bool {
    accuracy.reached(&result.volume) && accuracy.reached(&result.area)
        && result.centroid.as_ref().unwrap().iter().all(|c| accuracy.reached(c))
        && result.inertia.as_ref().is_none_or(|tensor| inertia_entries(tensor).into_iter().all(|c| accuracy.reached(c)))
}

fn budget() -> Refusal {
    Refusal::new("measure.budget-exhausted", "no exact-target offset-prism enclosure reached within the measurement budget")
}

fn measure(
    definition: &EllipticOffsetDefinition, axial: [&Real; 2], accuracy: &Accuracy,
    inertia: bool, dec: &mut Decider,
) -> Answered<MassProperties> {
    accuracy.validate();
    dec.check_deadline()?;
    match dec.sign(&(axial[1] - axial[0])) {
        Outcome::Proven(StrictSign::Positive) => {},
        Outcome::Proven(_) | Outcome::ProvenCoincident(_) => return Ok(Answer::exact(MassProperties {
            volume: Enclosure::exact(Rational::zero()), area: Enclosure::exact(Rational::zero()),
            centroid: None, inertia: None,
        })),
        Outcome::NotExact(n) => return Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => return Err(r),
    }
    let Surface::EllipticCylinder(original) = definition.original() else {
        unreachable!("validated elliptic-offset definition");
    };
    let signed = match definition.side() {
        OffsetSide::Outward => definition.distance().clone(), OffsetSide::Inward => -definition.distance(),
    };
    // Integrate in the principal plane, avoiding translation-dependent interval
    // cancellation. The exact source frame transports centroid and inertia.
    let canonical = EllipticOffsetDefinition::new(EllipticCylinder {
        origin: Point3::origin(), axis: Vector3::new(Real::zero(), Real::zero(), Real::one()),
        major_dir: Vector3::new(Real::one(), Real::zero(), Real::zero()),
        major_sq: original.major_sq.clone(), minor_sq: original.minor_sq.clone(),
    }, signed, dec)?;
    let trace = trace_elliptic_offset_quarter(&canonical, &Real::zero(), 0, dec)?.value;
    let frame = Frame::new(original, axial, inertia);
    let program = program(inertia);
    let policy = *dec.policy();
    if policy.max_precision_bits == 0 || policy.max_work == 0 { return Err(budget()); }
    let mut precision = 128u32.min(policy.max_precision_bits);
    // These halves separate one of the two normal-denominator coordinates
    // from zero, before Taylor evaluation; refinement never changes geometry.
    let mut depth = 1;
    let mut remaining = policy.max_work;
    let mut best = None;
    let mut split_precision = false;
    loop {
        let mut work = match Work::new(&policy, precision, remaining, dec) {
            Ok(work) => work, Err(Stop::Refusal(r)) => return Err(r), Err(_) => break,
        };
        let attempt = integrate(&program, &trace, depth, &mut work)
            .and_then(|values| frame.properties(&values, &mut work));
        remaining = work.remaining();
        dec.check_deadline()?;
        match attempt {
            Ok(result) => {
                let done = reached(&result, accuracy);
                best = Some(result);
                if done { return Ok(Answer::exact(best.unwrap())); }
                if work.quadrature_dominates { depth += 1; }
                else if precision < policy.max_precision_bits {
                    precision = precision.saturating_mul(2).min(policy.max_precision_bits);
                } else { break; }
            },
            Err(Stop::Refusal(r)) => return Err(r),
            Err(Stop::Budget) => break,
            Err(Stop::Split | Stop::GaugeSplit) => {
                if !split_precision && precision < policy.max_precision_bits {
                    precision = precision.saturating_mul(2).min(policy.max_precision_bits);
                    split_precision = true;
                } else { depth += 1; split_precision = false; }
            },
        }
        if remaining == 0 { break; }
    }
    let result = best.ok_or_else(budget)?;
    let mut marks = MarkLog::new();
    accuracy.mark(&mut marks, &[&result.volume], "volume");
    accuracy.mark(&mut marks, &[&result.area], "area");
    accuracy.mark(&mut marks, &result.centroid.as_ref().unwrap().each_ref(), "centroid");
    if let Some(tensor) = &result.inertia {
        accuracy.mark(&mut marks, &inertia_entries(tensor), "inertia");
    }
    Ok(Answer::with_marks(result, marks))
}
