//! Arc-length queries over certified traces, independent of display tessellation.

use catallax_intersect::{TracePiece, trace};
use catallax_numbers::{Answer, Answered, Ball, Decider, Rational, Refusal};
use catallax_topology::{EdgeId, Model, Sense};

use crate::{
    Accuracy, Enclosure,
    integrate::{Builder, Program, Stop, Work},
};

struct Panel {
    piece: usize,
    lo: Rational,
    hi: Rational,
    length: Ball,
}

fn panel(
    pieces: &[TracePiece],
    program: &Program,
    piece: usize,
    lo: Rational,
    hi: Rational,
    work: &mut Work,
) -> Result<Panel, Stop> {
    let length = program.panel(&pieces[piece], &lo, &hi, work)?.remove(0);
    Ok(Panel {
        piece,
        lo,
        hi,
        length,
    })
}

#[allow(clippy::too_many_arguments)]
fn initial(
    pieces: &[TracePiece],
    program: &Program,
    piece: usize,
    lo: Rational,
    hi: Rational,
    depth: u32,
    work: &mut Work,
    panels: &mut Vec<Panel>,
) -> Result<(), Stop> {
    match panel(pieces, program, piece, lo.clone(), hi.clone(), work) {
        Ok(p) => {
            panels.push(p);
            Ok(())
        }
        Err(Stop::Split) if depth < 48 => {
            let mid = (&lo + &hi).mul_pow2(-1);
            initial(
                pieces,
                program,
                piece,
                lo,
                mid.clone(),
                depth + 1,
                work,
                panels,
            )?;
            initial(pieces, program, piece, mid, hi, depth + 1, work, panels)
        }
        Err(e) => Err(e),
    }
}

pub(crate) fn measure(
    model: &Model,
    edge_id: EdgeId,
    accuracy: &Accuracy,
    midpoint: bool,
    dec: &mut Decider,
) -> Answered<(Enclosure, Option<[Enclosure; 3]>)> {
    accuracy.validate();
    if !midpoint && let Some(answer) = analytic(model, edge_id, accuracy, dec)? {
        return Ok(Answer::with_marks((answer.value, None), answer.marks));
    }
    let edge = model.edge(edge_id);
    let (mut from, mut to) = (edge.start, edge.end);
    if edge.sense == Sense::Reversed {
        std::mem::swap(&mut from, &mut to);
    }
    let traced = trace(
        model.geom(),
        edge.curve,
        from.map(|v| model.vertex(v).point),
        to.map(|v| model.vertex(v).point),
        dec,
    )?;
    let pieces = traced.value.pieces();
    let b = Builder::default();
    let speed = (b.input(3).square() + b.input(4).square() + b.input(5).square()).sqrt();
    let program = b.finish(&[speed]);
    let policy = *dec.policy();
    if policy.max_precision_bits == 0 || policy.max_work == 0 {
        return Err(Refusal::new(
            "measure.budget-exhausted",
            "no edge enclosure within budget",
        ));
    }
    let mut work = Work::new(
        &policy,
        policy.max_precision_bits.min(256),
        policy.max_work,
        dec,
    )
    .map_err(refusal)?;
    let mut panels = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        let (lo, hi) = piece.range();
        initial(pieces, &program, i, lo, hi, 0, &mut work, &mut panels).map_err(refusal)?;
    }
    let zero = Ball::from_rational(&Rational::zero(), work.precision);
    let mut best = None;
    for _ in 0..256 {
        dec.check_deadline()?;
        let total = panels.iter().fold(zero.clone(), |s, p| &s + &p.length);
        let length = Enclosure::from_ball(&total);
        let target_lo = total.lo().mul_pow2(-1);
        let target_hi = total.hi().mul_pow2(-1);
        let mut cumulative = zero.clone();
        let mut candidates = Vec::new();
        let mut positions: Option<[Ball; 3]> = None;
        for (i, p) in panels.iter().enumerate() {
            let next = &cumulative + &p.length;
            if midpoint && cumulative.lo() <= target_hi && next.hi() >= target_lo {
                let parameter = Ball::from_rational(&p.lo, work.precision)
                    .hull(&Ball::from_rational(&p.hi, work.precision));
                let xyz = pieces[p.piece].evaluate(&parameter, dec)?;
                positions = Some(match positions {
                    Some(old) => std::array::from_fn(|j| old[j].hull(&xyz[j])),
                    None => xyz,
                });
                candidates.push(i);
            }
            cumulative = next;
        }
        let center = positions.map(|xyz| xyz.each_ref().map(Enclosure::from_ball));
        let length_reached = accuracy.reached(&length)
            && (!midpoint
                || accuracy.absolute <= Rational::zero()
                || length.width() <= accuracy.absolute.mul_pow2(-4));
        let center_reached = !midpoint
            || center
                .as_ref()
                .is_some_and(|xyz| xyz.iter().all(|v| accuracy.reached(v)));
        best = Some((length, center));
        if length_reached && center_reached {
            break;
        }
        let mut refined = Vec::new();
        let mut failure = None;
        for (i, p) in panels.into_iter().enumerate() {
            if !length_reached || candidates.binary_search(&i).is_ok() {
                let mid = (&p.lo + &p.hi).mul_pow2(-1);
                if let Err(e) = initial(
                    pieces,
                    &program,
                    p.piece,
                    p.lo,
                    mid.clone(),
                    0,
                    &mut work,
                    &mut refined,
                )
                .and_then(|()| {
                    initial(
                        pieces,
                        &program,
                        p.piece,
                        mid,
                        p.hi,
                        0,
                        &mut work,
                        &mut refined,
                    )
                }) {
                    failure = Some(e);
                    break;
                }
            } else {
                refined.push(p);
            }
        }
        match failure {
            Some(Stop::Refusal(r)) => return Err(r),
            Some(_) => break,
            None => panels = refined,
        }
    }
    let (length, center) = best.ok_or_else(|| refusal(Stop::Budget))?;
    let mut marks = traced.marks;
    if !midpoint {
        accuracy.mark(&mut marks, &[&length], "edge length");
    }
    if let Some(xyz) = &center {
        accuracy.mark(&mut marks, &xyz.each_ref(), "edge midpoint");
    }
    Ok(Answer::with_marks((length, center), marks))
}

fn refusal(stop: Stop) -> Refusal {
    match stop {
        Stop::Refusal(r) => r,
        _ => Refusal::new(
            "measure.budget-exhausted",
            "no edge enclosure within budget",
        ),
    }
}

fn analytic(
    model: &Model,
    edge_id: EdgeId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Result<Option<Answer<Enclosure>>, Refusal> {
    use catallax_shapes::Curve;
    let edge = model.edge(edge_id);
    let bits = dec.policy().max_precision_bits.min(256);
    if bits == 0 || dec.policy().max_work == 0 {
        return Ok(None);
    }
    let policy = *dec.policy();
    let mut work = Work::new(&policy, bits, policy.max_work, dec).map_err(refusal)?;
    let length = match model.geom().curve(edge.curve) {
        Curve::Line(_) => {
            let (Some(a), Some(b)) = (edge.start, edge.end) else {
                return Ok(None);
            };
            let a = model.geom().position(model.vertex(a).point);
            let b = model.geom().position(model.vertex(b).point);
            work.real(&b.minus(&a).norm_squared().sqrt()).map_err(refusal)?
        }
        Curve::Circle(c) if edge.start == edge.end => {
            let radius = work.real(&c.radius_sq.sqrt()).map_err(refusal)?;
            let pi = work.pi().map_err(refusal)?;
            let two = work.rational(&Rational::from_i64(2)).map_err(refusal)?;
            let product = work.mul(&pi, &radius).map_err(refusal)?;
            work.mul(&product, &two).map_err(refusal)?
        }
        _ => return Ok(None),
    };
    let length = Enclosure::from_ball(&length);
    let mut marks = catallax_numbers::MarkLog::new();
    accuracy.mark(&mut marks, &[&length], "edge length");
    Ok(Some(Answer::with_marks(length, marks)))
}

pub(crate) fn midpoint(
    model: &Model,
    edge_id: EdgeId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<[Enclosure; 3]> {
    use catallax_shapes::{Curve, Point3};
    accuracy.validate();
    let edge = model.edge(edge_id);
    let two = catallax_numbers::Real::from(2);
    let point = match model.geom().curve(edge.curve) {
        Curve::Line(_) => edge.start.zip(edge.end).map(|(a, b)| {
            let a = model.geom().position(model.vertex(a).point);
            let b = model.geom().position(model.vertex(b).point);
            Point3::new(
                &(&a.x + &b.x) / &two,
                &(&a.y + &b.y) / &two,
                &(&a.z + &b.z) / &two,
            )
        }),
        Curve::Circle(c) if edge.start == edge.end => {
            let start = if let Some(v) = edge.start {
                model.geom().position(model.vertex(v).point)
            } else {
                c.center.plus(
                    &c.ref_dir
                        .scaled(&(&c.radius_sq / &c.ref_dir.norm_squared()).sqrt()),
                )
            };
            Some(Point3::new(
                &(&c.center.x * &two) - &start.x,
                &(&c.center.y * &two) - &start.y,
                &(&c.center.z * &two) - &start.z,
            ))
        }
        _ => None,
    };
    if let Some(point) = point {
        let bits = dec.policy().max_precision_bits.min(256);
        let coordinates = [&point.x, &point.y, &point.z].map(|x| {
            if let Some(value) = x.as_rational() {
                Ok(Enclosure::exact(value))
            } else if bits == 0 || dec.policy().max_work == 0 {
                Err(refusal(Stop::Budget))
            } else {
                dec.enclose_at(x, bits).map(|b| Enclosure::from_ball(&b))
            }
        });
        let [x, y, z] = coordinates;
        let values = [x?, y?, z?];
        let mut marks = catallax_numbers::MarkLog::new();
        accuracy.mark(&mut marks, &values.each_ref(), "edge midpoint");
        return Ok(Answer::with_marks(values, marks));
    }
    let answer = measure(model, edge_id, accuracy, true, dec)?;
    Ok(Answer::with_marks(answer.value.1.unwrap(), answer.marks))
}
