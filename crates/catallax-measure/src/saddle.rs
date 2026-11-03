//! What a rolling saddle changes on the body under it, as integrals round the
//! ring. Rails trim two cylinders, which the engine cannot chart, so the blend is
//! read as its unfilleted base plus this change.
//!
//! The law is the one in `catallax_shapes::rolling_saddle`, in the frame of the
//! narrow cylinder: `u1` along `across`, `u2` along `C along - S axis` and `h` up
//! the narrow axis `S along + C axis`, all from its foot. The pocket between the
//! ball and the corner stands on the wide cylinder in columns along the narrow
//! axis, and it stops on the narrow cylinder, which those columns never cross. So
//! for any `g` of `(u1, u2)`,
//!
//! ```text
//! integral of g over the pocket = shadow integral of g (band h - floor)
//! ```
//!
//! where `floor` is the wide cylinder's height. `floor = (w - C u2) / S`, for
//! `w = branch sqrt(R^2 - (gap + u1)^2)`, carries `u2` only linearly, so a potential
//! in `u2` turns its shadow integral by Green's theorem into two integrals round
//! closed curves in the shadow, the wide rail's and the narrow circle. Both have
//! `u1` rational in `cos t` and `sin t` and a known `w`.
//!
//! The band is a canal surface about the centre path, so with `T` that path's unit
//! tangent, a point is `c + r e` for `e` in the normal plane from the wide contact
//! ray `A` to the narrow one `B`, and `dA = r (|c'| - r e . T') ds dt`. Rays are
//! square to `T`, so `A . T' = -A' . T`, and the law gives `A'` and `B'` directly.
//! The band is a graph over the shadow, so `|e . d| dA` is its shadow's area, and
//! `e . d` keeps the sign of `A . d`. Each integrand over `s` is a polynomial in
//! `cos s` and `sin s`, as on a rolling pipe, which leaves one integral in `t`.
//!
//! Area works the same way: each cylinder loses the strip between the input ring
//! and its rail, and on a cylinder's flat chart that strip's area is a closed
//! integral of height against angle.
//!
//! These enclose the exact definition. `Scope::Represented` in `junction` adds the
//! tube bound its module note derives, so its volume also holds for the carrier.

use catallax_numbers::{Ball, Decider, Integer, Rational, Refusal};
use catallax_shapes::{RollingSaddle, SaddleLaw, SaddleRing};
use std::collections::BTreeMap;

use crate::integrate::{Builder, Program, Stop, Work};
use crate::pipe::{Primitives, Section};

// A fitted band's own bound is far wider than these integrals, so 96 bits leave
// room. Eight panels split the section angle's enclosure, and sixteen do not, as
// on a rolling pipe. Each further depth doubles them.
const BITS: u32 = 96;
const PANELS: u64 = 16;
const MAX_DEPTH: u32 = 6;

/// What one band changes, before its sign: the pocket's volume and first moments
/// about the world origin, and the change of area.
pub(crate) struct Change {
    pub(crate) pocket: Ball,
    pub(crate) moment: Option<[Ball; 3]>,
    pub(crate) area: Ball,
    /// The band's area, its integrated curvature bound, and its rails' length.
    pub(crate) band: Ball,
    pub(crate) turning: Ball,
    pub(crate) length: Ball,
}

fn refused(detail: &'static str) -> Refusal {
    Refusal::new("measure.unsupported-surface", detail)
}

fn stop(stop: Stop) -> Refusal {
    match stop {
        Stop::Refusal(refusal) => refusal,
        Stop::Budget => Refusal::new(
            "numbers.work-budget",
            "a saddle's integrals ran out of work",
        ),
        Stop::Split | Stop::GaugeSplit => refused("a saddle's integrand has no enclosure"),
    }
}

// Outputs in order: signed volume, band area, the wide strip's two closed
// integrals, the narrow strip, the turning bound and the rail length bound, then,
// when asked for, the signed moments along `u1`, `u2` and `h`.
fn forms(law: &SaddleLaw, moments: bool) -> Program {
    let k = Builder::default();
    let [wide, narrow, ball, gap, c0, s0] = [
        &law.wide,
        &law.narrow,
        &law.radius,
        &law.gap,
        &law.cosine,
        &law.sine,
    ]
    .map(|value| k.exact(value));
    let [sw, sn] = law.sides.map(|side| k.k(side));
    let branch = k.k(law.branch);
    let (cos, sin) = (k.input(0), k.input(1));
    let two = k.k(2);
    let three = k.k(3);
    let big = wide + sw * ball;
    let small = narrow + sn * ball;
    let ratio = wide / big;

    // The centre path and its speed, in the narrow frame.
    let x = gap + small * cos;
    let root = (big.square() - x.square()).sqrt();
    let y = branch * root;
    let dx = -(small * sin);
    let dy = -(x * dx) / y;
    let height = (y - small * c0 * sin) / s0;
    let rise = (dy - small * c0 * cos) / s0;
    let speed = (small.square() + rise.square()).sqrt();
    let centre = [small * cos, small * sin, height];

    // The contact rays, the second axis of their plane, and the path's turning
    // along each.
    let first = [-(sw * x) / big, -(sw * c0 * y) / big, -(sw * s0 * y) / big];
    let toward = [-(sn * cos), -(sn * sin), k.k(0)];
    let cosine = first[0] * toward[0] + first[1] * toward[1] + first[2] * toward[2];
    let sine = (k.k(1) - cosine.square()).sqrt();
    let second: [_; 3] = std::array::from_fn(|i| (toward[i] - cosine * first[i]) / sine);
    let turn_first = sw * (small.square() * sin.square() + c0 * dy * small * cos + s0 * dy * rise)
        / (big * speed);
    let turn_toward = sn * small / speed;
    let turn_second = (turn_toward - cosine * turn_first) / sine;

    // The band.
    let zero = k.k(0);
    let density = Section::linear(
        ball * speed,
        -(ball.square() * turn_first),
        -(ball.square() * turn_second),
    );
    let position: [Section; 3] =
        std::array::from_fn(|i| Section::linear(centre[i], ball * first[i], ball * second[i]));
    let shadow = Section::linear(zero, first[2], second[2]).times(&density);
    let tilt = -(sw * branch);
    let angle = cosine.angle_to(sine);
    let mut primitives = Primitives {
        builder: &k,
        angle,
        cosine,
        sine,
        known: BTreeMap::new(),
    };
    let volume_band = tilt * primitives.integral(&position[2].times(&shadow));
    let moment_band = moments.then(|| {
        [
            tilt * primitives.integral(&position[0].times(&position[2]).times(&shadow)),
            tilt * primitives.integral(&position[1].times(&position[2]).times(&shadow)),
            tilt * primitives.integral(&position[2].times(&position[2]).times(&shadow)) / two,
        ]
    });
    let band_area = primitives.integral(&density);

    // The wide rail and the input ring, in the shadow. Green's theorem takes each
    // shadow integral to these two, with the rail outermost for a ball outside
    // the narrow cylinder.
    let rail_u1 = x * ratio - gap;
    let rail_w = y * ratio;
    let rail_z = (c0 * y - small * sin) / s0;
    let rail_u2 = c0 * rail_w - s0 * rail_z;
    let rail_du = dx * ratio;
    let ring_u1 = narrow * cos;
    let ring_u2 = narrow * sin;
    let ring_x = gap + ring_u1;
    let ring_root = (wide.square() - ring_x.square()).sqrt();
    let ring_w = branch * ring_root;
    let ring_du = -(narrow * sin);
    let potentials = |u1, u2, w| {
        let plain = (w * u2 - c0 * u2 * u2 / two) / s0;
        [
            plain,
            u1 * plain,
            (w * u2 * u2 / two - c0 * u2 * u2 * u2 / three) / s0,
            (w * w * u2 - w * c0 * u2 * u2 + c0 * c0 * u2 * u2 * u2 / three) / (two * s0 * s0),
        ]
    };
    let on_rail = potentials(rail_u1, rail_u2, rail_w);
    let on_ring = potentials(ring_u1, ring_u2, ring_w);
    let closed = |i: usize| sn * (on_rail[i] * rail_du - on_ring[i] * ring_du);

    // The strips each cylinder loses.
    let wide_rail = wide * rail_z * dx / root;
    let ring_z = (c0 * ring_w - narrow * sin) / s0;
    let wide_ring = wide * ring_z * ring_du / ring_root;
    let narrow_strip = narrow * (height - (ring_w - c0 * narrow * sin) / s0);

    // Smooth upper bounds for the slack, by `|v| <= (1 + |v|^2) / 2`.
    let half = |squared| (k.k(1) + squared) / two;
    let turning = ball * angle * half(turn_first.square() + turn_second.square());
    let rail_rise = (c0 * dy - small * cos) / s0;
    let length = half(ratio * ratio * (dx * dx + dy * dy) + rail_rise * rail_rise)
        + half(narrow.square() + rise.square());

    let mut outputs = vec![
        volume_band + closed(0),
        band_area,
        wide_rail,
        wide_ring,
        narrow_strip,
        turning,
        length,
    ];
    if let Some([along, round, up]) = moment_band {
        outputs.extend([along + closed(1), round + closed(2), up + closed(3)]);
    }
    k.finish(&outputs)
}

fn width(value: &Ball) -> Rational {
    value.hi() - value.lo()
}

fn magnitude(value: &Ball) -> Ball {
    if value.hi() < Rational::zero() {
        &Ball::from_rational(&Rational::zero(), BITS) - value
    } else {
        value.clone()
    }
}

fn panels(program: &Program, depth: u32, work: &mut Work) -> Result<Vec<Ball>, Stop> {
    let count = PANELS.checked_shl(depth).ok_or(Stop::Budget)?;
    let step = Rational::from_i64(2) / Rational::from_integer(Integer::from_u64(count));
    let mut sums: Option<Vec<Ball>> = None;
    for i in 0..count {
        let a = &step * &Rational::from_integer(Integer::from_u64(i));
        let b = &a + &step;
        let values = program.turning_panel(&a, &b, work)?;
        sums = Some(match sums {
            None => values,
            Some(sums) => sums
                .iter()
                .zip(&values)
                .map(|(sum, value)| work.add(sum, value))
                .collect::<Result<_, _>>()?,
        });
    }
    let pi = work.pi()?;
    sums.ok_or(Stop::Budget)?
        .iter()
        .map(|sum| work.mul(sum, &pi))
        .collect()
}

/// The band's change, each integral no wider than `target` where the budget
/// reaches it.
pub(crate) fn change(
    saddle: &RollingSaddle,
    target: &Rational,
    moments: bool,
    dec: &mut Decider,
) -> Result<Change, Refusal> {
    let law = saddle.definition().law();
    // The integrals run on `t` round the narrow cylinder, which a loop does not.
    if law.ring != SaddleRing::Through {
        return Err(refused("a loop's change has no integrals yet"));
    }
    let program = forms(&law, moments);
    let policy = *dec.policy();
    let bits = BITS.min(policy.max_precision_bits);
    let mut remaining = policy.max_work;
    let mut values = None;
    for depth in 0..=MAX_DEPTH {
        let mut work = Work::new(&policy, bits, remaining, dec).map_err(stop)?;
        let attempt = panels(&program, depth, &mut work);
        remaining = work.remaining();
        match attempt {
            Ok(found) => {
                let done = found.iter().all(|value| width(value) <= *target);
                values = Some(found);
                if done {
                    break;
                }
            }
            Err(Stop::Split | Stop::GaugeSplit) if depth < MAX_DEPTH => {}
            Err(other) => return Err(stop(other)),
        }
    }
    let values = values.ok_or_else(|| refused("a saddle's integrand has no enclosure"))?;
    let [
        volume,
        band,
        rail,
        ring,
        strip,
        turning,
        length,
        moment @ ..,
    ] = values.as_slice()
    else {
        return Err(refused("a saddle's program has a fixed output count"));
    };
    let sign = Ball::from_rational(&Rational::from_i64(law.sides[0] * law.branch), BITS);
    let pocket = volume * &sign;
    let vector = |v: &catallax_shapes::Vector3| -> Result<[Ball; 3], Refusal> {
        Ok([
            dec.enclose_at(&v.x, BITS)?,
            dec.enclose_at(&v.y, BITS)?,
            dec.enclose_at(&v.z, BITS)?,
        ])
    };
    let [across, side, axis] = [&law.across, &law.along, &law.axis].map(vector);
    let (across, side, axis) = (across?, side?, axis?);
    let (cosine, sine) = (
        dec.enclose_at(&law.cosine, BITS)?,
        dec.enclose_at(&law.sine, BITS)?,
    );
    let foot = [
        dec.enclose_at(&law.foot.x, BITS)?,
        dec.enclose_at(&law.foot.y, BITS)?,
        dec.enclose_at(&law.foot.z, BITS)?,
    ];
    let moment = moment.first_chunk::<3>().map(|[along, round, up]| {
        std::array::from_fn(|i| {
            let round_axis = &(&side[i] * &cosine) - &(&axis[i] * &sine);
            let narrow_axis = &(&side[i] * &sine) + &(&axis[i] * &cosine);
            let local = &(&(&across[i] * along) + &(&round_axis * round)) + &(&narrow_axis * up);
            &(&foot[i] * &pocket) + &(&local * &sign)
        })
    });
    let area = &(band - &magnitude(&(rail - ring))) - &magnitude(strip);
    Ok(Change {
        pocket,
        moment,
        area,
        band: band.clone(),
        turning: turning.clone(),
        length: length.clone(),
    })
}
