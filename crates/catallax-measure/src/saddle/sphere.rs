//! What a rolling-ball band between a sphere and a cylinder changes on the
//! body under it, as integrals round the ring.
//!
//! The law is the one in `catallax_shapes::rolling_sphere`, in its own frame:
//! the sphere's centre at the origin, `x` along `across`, `y` along `along` and
//! `z` up the cylinder axis, which runs through `(gap, 0, z)`. With offset
//! radii `P` of the sphere and `p` of the cylinder, the ball centre at angle `t`
//! is `c = (gap + p cos t, p sin t, branch sqrt(P^2 - x^2 - y^2))`.
//!
//! The pocket between the ball and the corner stands on the sphere in columns
//! along the cylinder axis, and it stops on the cylinder, which those columns
//! never cross. Its shadow lies between the ring's, which is the cylinder's
//! circle of radius `c` about `(gap, 0)`, and the sphere rail's, which is the
//! centre path's shadow scaled by `R / P`: a circle of radius `p R / P` about
//! `(gap R / P, 0)`. So for any `g` of `(x, y)`,
//!
//! ```text
//! integral of g over the pocket = shadow integral of g (band z - floor)
//! ```
//!
//! up to the pocket's sign, where `floor = branch s` and `s = sqrt(R^2 - x^2 - y^2)`
//! is the sphere's height. The band term is the canal-surface integral of
//! `crate::saddle`: the band is a graph over the shadow, and `e . z` keeps the
//! sign of the sphere contact ray's `z`, which is `-side branch`.
//!
//! The floor term goes round the two circles by Green's theorem, with closed
//! potentials. A radial `f(rho)` integrates as the flux of `F(rho) / rho^2 (x, y)`
//! with `F(rho) = integral of f(u) u du from 0`, so
//!
//! ```text
//! shadow integral of s     = closed integral of (R^2 + R s + s^2) / (3 (R + s)) (x dy - y dx)
//! shadow integral of R / s = closed integral of R / (R + s) (x dy - y dx)
//! shadow integral of x s   = closed integral of -s^3 / 3 dy
//! shadow integral of y s   = closed integral of  s^3 / 3 dx
//! shadow integral of s^2/2 = closed integral of (R^2 x - x^3 / 3 - x y^2) / 2 dy
//! ```
//!
//! The second is the area of the sphere over the shadow. The cylinder loses the
//! strip between the ring and its rail, which on its flat chart has area
//! `c times the integral of (rail z - ring z) dt`. The definition requires
//! `gap + c < R` and `gap + p < P`, so `s` stays positive on both circles and no
//! potential divides by zero.
//!
//! On the canal surface, with rays `A = -side c / P` toward the sphere and
//! `B = -side (cos t, sin t, 0)` toward the cylinder, `A . T' = side |c'| / P` and
//! `B . T' = side p / |c'|` for the centre path's unit tangent `T`.
//!
//! These enclose the exact definition. `Scope::Represented` in `junction` adds the
//! tube bound its module note derives, so its volume also holds for the carrier.

use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_shapes::{RollingSphere, RollingSphereDefinition};
use std::collections::BTreeMap;

use super::{BITS, Change, integrals, magnitude, refused, stop};
use crate::integrate::{Builder, Program, Work};
use crate::pipe::{Primitives, Section};

// Outputs in order: signed volume, band area, the sphere's area over the
// shadow, the cylinder strip, the turning bound and the rail length bound,
// then, when asked for, the signed moments along `x`, `y` and `z`.
fn forms(law: &RollingSphereDefinition, moments: bool) -> Program {
    let k = Builder::default();
    let [sphere, cylinder, ball, gap] = [
        law.sphere_radius(),
        law.cylinder_radius(),
        law.radius(),
        law.gap(),
    ]
    .map(|value| k.exact(value));
    let [ss, sc] = law.sides().map(|side| k.k(side));
    let branch = k.k(law.branch());
    let (cos, sin) = (k.input(0), k.input(1));
    let zero = k.k(0);
    let two = k.k(2);
    let three = k.k(3);
    let big = sphere + ss * ball;
    let small = cylinder + sc * ball;
    let ratio = sphere / big;

    // The centre path and its speed.
    let x = gap + small * cos;
    let y = small * sin;
    let root = (big.square() - x.square() - y.square()).sqrt();
    let z = branch * root;
    let dx = -(small * sin);
    let dy = small * cos;
    let dz = gap * small * sin / z;
    let speed = (small.square() + dz.square()).sqrt();
    let centre = [x, y, z];

    // The contact rays, the second axis of their plane, and the path's turning
    // along each.
    let first = [-(ss * x) / big, -(ss * y) / big, -(ss * z) / big];
    let toward = [-(sc * cos), -(sc * sin), zero];
    let cosine = first[0] * toward[0] + first[1] * toward[1];
    let sine = (k.k(1) - cosine.square()).sqrt();
    let second: [_; 3] = std::array::from_fn(|i| (toward[i] - cosine * first[i]) / sine);
    let turn_first = ss * speed / big;
    let turn_toward = sc * small / speed;
    let turn_second = (turn_toward - cosine * turn_first) / sine;

    // The band.
    let density = Section::linear(
        ball * speed,
        -(ball.square() * turn_first),
        -(ball.square() * turn_second),
    );
    let position: [Section; 3] =
        std::array::from_fn(|i| Section::linear(centre[i], ball * first[i], ball * second[i]));
    let shadow = Section::linear(zero, first[2], second[2]).times(&density);
    let tilt = -(ss * branch);
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

    // The floor and the sphere's area over the shadow, round the sphere rail's
    // shadow and the ring, which is the cylinder's circle. The rail is outermost
    // for a ball outside the cylinder.
    let potentials = |u, v, du, dv| {
        let s = (sphere.square() - u * u - v * v).sqrt();
        let flux = u * dv - v * du;
        let cube = s * s * s;
        [
            branch * (sphere.square() + sphere * s + s * s) / (three * (sphere + s)) * flux,
            branch * -(cube / three) * dv,
            branch * (cube / three) * du,
            (sphere.square() * u - u * u * u / three - u * v * v) / two * dv,
            sphere / (sphere + s) * flux,
        ]
    };
    let on_rail = potentials(ratio * x, ratio * y, ratio * dx, ratio * dy);
    let ring_u = gap + cylinder * cos;
    let ring_v = cylinder * sin;
    let on_ring = potentials(ring_u, ring_v, -(cylinder * sin), cylinder * cos);
    let closed = |i: usize| sc * (on_rail[i] - on_ring[i]);

    // The strip the cylinder loses, between the ring and the rail at the centre's
    // height.
    let ring_z = branch * (sphere.square() - ring_u.square() - ring_v.square()).sqrt();
    let cylinder_strip = cylinder * (z - ring_z);

    // Smooth upper bounds for the slack, by `|v| <= (1 + |v|^2) / 2`.
    let half = |squared| (k.k(1) + squared) / two;
    let turning = ball * angle * half(turn_first.square() + turn_second.square());
    let length = half(ratio * ratio * speed.square()) + half(cylinder.square() + dz.square());

    let mut outputs = vec![
        volume_band - closed(0),
        band_area,
        closed(4),
        cylinder_strip,
        turning,
        length,
    ];
    if let Some([along, round, up]) = moment_band {
        outputs.extend([along - closed(1), round - closed(2), up - closed(3)]);
    }
    k.finish(&outputs)
}

/// The band's change, each integral no wider than `target` where the budget
/// reaches it.
pub(crate) fn change(
    band: &RollingSphere,
    target: &Rational,
    moments: bool,
    dec: &mut Decider,
) -> Result<Change, Refusal> {
    let law = band.definition();
    let program = forms(law, moments);
    let range = [Rational::zero(), Rational::from_i64(2)];
    let policy = *dec.policy();
    let bits = BITS.min(policy.max_precision_bits);
    let values = {
        let mut work = Work::new(&policy, bits, policy.max_work, dec).map_err(stop)?;
        integrals(&program, target, range, false, &mut work).map_err(stop)?
    };
    let [
        volume,
        band_area,
        sphere_area,
        cylinder_strip,
        turning,
        length,
        moment @ ..,
    ] = values.as_slice()
    else {
        return Err(refused("a sphere band's program has a fixed output count"));
    };
    // Above the sphere the band stands over its floor; below it, under it.
    let [sphere_side, _] = law.sides();
    let sign = Ball::from_rational(&Rational::from_i64(sphere_side * law.branch()), BITS);
    let pocket = volume * &sign;
    if pocket.lo() <= Rational::zero() {
        return Err(refused("a sphere band's pocket has no certified volume"));
    }
    let vector = |v: &catallax_shapes::Vector3| -> Result<[Ball; 3], Refusal> {
        Ok([
            dec.enclose_at(&v.x, BITS)?,
            dec.enclose_at(&v.y, BITS)?,
            dec.enclose_at(&v.z, BITS)?,
        ])
    };
    let origin = law.centre().minus(&catallax_shapes::Point3::origin());
    let [origin, across, along, axis] =
        [&origin, law.across(), law.along(), law.axis()].map(vector);
    let (origin, across, along, axis) = (origin?, across?, along?, axis?);
    let moment = moment.first_chunk::<3>().map(|[mx, my, mz]| {
        std::array::from_fn(|i| {
            let local = &(&(&across[i] * mx) + &(&along[i] * my)) + &(&axis[i] * mz);
            &(&origin[i] * &pocket) + &(&local * &sign)
        })
    });
    let area = &(band_area - &magnitude(sphere_area)) - &magnitude(cylinder_strip);
    Ok(Change {
        pocket,
        moment,
        local_moment: None,
        area,
        band: band_area.clone(),
        turning: turning.clone(),
        length: length.clone(),
    })
}
