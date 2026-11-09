//! What a rolling saddle changes on a loop, by the divergence theorem.
//!
//! A loop crosses the wide cylinder's outline seen along the narrow axis, so the
//! columns the through ring stands on the wide cylinder graze it there. Here the
//! pocket between the ball and the corner is bounded by the band and by a strip
//! of each wall, between the input loop and the rail. With `O` the law's origin,
//! on the wide axis, and `n` the pocket's outward normal,
//!
//! ```text
//! volume = (1/3) integral of (x - O) . n,    moment_k = integral of (x - O)_k^2 n_k / 2
//! ```
//!
//! over those three pieces. On the band `x = c + r e` and `n = -e`, with the area
//! element of a canal surface, so each integrand over the section is a
//! polynomial in its angle, as on a through ring. On the wide wall
//! `(x - O) . n = -s_w R`. On the narrow one it is `-s_n (gap cos t + r)`. Each
//! strip integral is `integral of g dh dangle` on that cylinder's flat chart.
//! Green's theorem turns it into the rail's closed integral of `G dangle` less
//! the input loop's, with `G` a primitive of `g` in the height. Both loops run the
//! way `w` grows, and the area's own sign orients the difference.
//!
//! Both loops follow the law in `catallax_shapes::rolling_saddle`. With
//! `c = (P + g - p) / 2` and `d = (P - g + p) / 2`, `X = c + d cos 2w`,
//! `y = A sin w` and `s = B cos w`, where `A^2 = 4 d (P - d sin^2 w)` and
//! `B^2 = 4 d (p - d cos^2 w)`. The input loop takes the cylinders' own radii.

use super::*;
use crate::integrate::Expr;

/// Outputs before any moments: the band's flux and area, then the wide strip's
/// area, the narrow strip's area and flux, the turning bound and the rail length
/// bound, all before their signs.
pub(super) const OUTPUTS: usize = 7;

// A loop of the law, in the law's frame, with its derivatives in `w`.
struct Path<'a> {
    x: Expr<'a>,
    y: Expr<'a>,
    z: Expr<'a>,
    s: Expr<'a>,
    dx: Expr<'a>,
    dy: Expr<'a>,
    dz: Expr<'a>,
    ds: Expr<'a>,
}

// One half of an ellipse, with both source rails meeting at its poles.
fn crossing_path<'a>(
    radius: Expr<'a>,
    [a, b]: [Expr<'a>; 2],
    [c0, s0]: [Expr<'a>; 2],
    cos: Expr<'a>,
    sin: Expr<'a>,
) -> Path<'a> {
    let (x, y, s) = (radius * sin, a * radius * cos, b * radius * cos);
    let (dx, dy, ds) = (radius * cos, -(a * radius * sin), -(b * radius * sin));
    Path {
        x,
        y,
        s,
        z: (c0 * y - s) / s0,
        dx,
        dy,
        ds,
        dz: (c0 * dy - ds) / s0,
    }
}

#[allow(clippy::too_many_arguments)]
fn path<'a>(
    k: &'a Builder,
    big: Expr<'a>,
    small: Expr<'a>,
    reach: Expr<'a>,
    facing: Expr<'a>,
    [c0, s0]: [Expr<'a>; 2],
    cos: Expr<'a>,
    sin: Expr<'a>,
) -> Path<'a> {
    let two = k.k(2);
    let middle = (big + reach - small) / two;
    let half = (big - reach + small) / two;
    let x = facing * (middle + half * (cos * cos - sin * sin));
    let dx = -(facing * k.k(4) * half * sin * cos);
    let wide = half * (big - half * sin * sin);
    let narrow = half * (small - half * cos * cos);
    let (a, b) = (two * wide.sqrt(), two * narrow.sqrt());
    let da = -(two * half * half * sin * cos) / wide.sqrt();
    let db = two * half * half * sin * cos / narrow.sqrt();
    let y = sin * a;
    let dy = cos * a + sin * da;
    let s = cos * b;
    let ds = cos * db - sin * b;
    Path {
        x,
        y,
        z: (c0 * y - s) / s0,
        s,
        dx,
        dy,
        dz: (c0 * dy - ds) / s0,
        ds,
    }
}

// A loop's turn on the wide wall's flat chart, for a cylinder of radius `rho`.
fn wide_turn<'a>(p: &Path<'a>, rho: Expr<'a>) -> Expr<'a> {
    (p.x * p.dy - p.y * p.dx) / (rho * rho)
}

// A loop's turn on the narrow wall's flat chart.
fn narrow_turn<'a>(p: &Path<'a>, gap: Expr<'a>, rho: Expr<'a>) -> Expr<'a> {
    ((p.x - gap) * p.ds - p.s * p.dx) / (rho * rho)
}

fn cube(value: Expr<'_>) -> Expr<'_> {
    value * value * value
}

// The rail's term less the input loop's, each given its loop and its offsets.
fn difference<'a>(
    loops: [(&Path<'a>, Expr<'a>, Expr<'a>); 2],
    term: impl Fn(&Path<'a>, Expr<'a>, Expr<'a>) -> Expr<'a>,
) -> Expr<'a> {
    let [(a, a_wide, a_narrow), (b, b_wide, b_narrow)] = loops;
    term(a, a_wide, a_narrow) - term(b, b_wide, b_narrow)
}

// Outputs in order: the band's flux of `x - O` along `e` and its area, the wide
// strip's chart area, the narrow strip's chart area and flux, the turning bound
// and the rail length bound, then, when asked for, the band's three moment fluxes,
// the wide strip's two and the narrow strip's three.
pub(super) fn forms(law: &SaddleLaw, moments: bool, slack: bool) -> Program {
    let k = Builder::default();
    build(&k, law, moments, true, slack)
}

fn build<'a>(
    k: &'a Builder,
    law: &SaddleLaw,
    moments: bool,
    symmetry: bool,
    slack: bool,
) -> Program {
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
    let facing = k.k(law.gap_sign);
    let reach = facing * gap;
    let (cos, sin) = (k.input(0), k.input(1));
    let (zero, two, three) = (k.k(0), k.k(2), k.k(3));
    let big = wide + sw * ball;
    let small = narrow + sn * ball;
    let crossing = match law.ring {
        SaddleRing::Crossing { narrow_positive } => {
            Some([k.k(law.branch), k.k(if narrow_positive { 1 } else { -1 })])
        }
        _ => None,
    };
    let on_path = |big, small| match crossing {
        Some(signs) => crossing_path(big, signs, [c0, s0], cos, sin),
        None => path(k, big, small, reach, facing, [c0, s0], cos, sin),
    };
    let centre = on_path(big, small);
    let edge = on_path(wide, narrow);

    // The band. Its contact rays run from the centre to each axis.
    let (x, y, z, s) = (centre.x, centre.y, centre.z, centre.s);
    let (dx, dy, dz, ds) = (centre.dx, centre.dy, centre.dz, centre.ds);
    let speed = (dx * dx + dy * dy + dz * dz).sqrt();
    let u = x - gap;
    let first = [-(sw * x) / big, -(sw * y) / big, zero];
    let (cosine, sine, second, turn_first, turn_second) = if let Some([a, b]) = crossing {
        let kappa = k.k(1) - a * b * c0;
        let cosine = k.k(1) - kappa * cos.square();
        let root = (kappa * (two - kappa * cos.square())).sqrt();
        // Factor the closing cos(phi) before dividing. The section basis and
        // its curvature then have finite limits at both collapsed poles.
        let second = [
            -(sw * sin * kappa * cos) / root,
            -(sw * (b * c0 - a * cosine)) / root,
            sw * b * s0 / root,
        ];
        let turn_first = sw * big / speed;
        (
            cosine,
            cos * root,
            second,
            turn_first,
            turn_first * kappa * cos / root,
        )
    } else {
        let toward = [
            -(sn * u) / small,
            -(sn * s * c0) / small,
            sn * s * s0 / small,
        ];
        let cosine = first[0] * toward[0] + first[1] * toward[1] + first[2] * toward[2];
        let sine = (k.k(1) - cosine.square()).sqrt();
        let second = std::array::from_fn(|i| (toward[i] - cosine * first[i]) / sine);
        let turn_first = sw * (dx * dx + dy * dy) / (big * speed);
        let swing = u * ds - s * dx;
        let turn_toward = sn * swing * swing / (small * small * small * speed);
        (
            cosine,
            sine,
            second,
            turn_first,
            (turn_toward - cosine * turn_first) / sine,
        )
    };
    let density = Section::linear(
        ball * speed,
        -(ball.square() * turn_first),
        -(ball.square() * turn_second),
    );
    let at = [x, y, z];
    let reach_out = Section::linear(
        ball,
        at[0] * first[0] + at[1] * first[1],
        at[0] * second[0] + at[1] * second[1] + at[2] * second[2],
    );
    let angle = if crossing.is_some() {
        // The primitive's angle has an open (0, pi) branch. Shift the
        // half-angle into its interior so a zero-width pole remains regular.
        let one = k.k(1);
        let half_cosine = ((one + cosine) / two).sqrt();
        let half_sine = sine / (two * half_cosine);
        two * ((-half_sine).angle_to(half_cosine) - zero.angle_to(one))
    } else {
        cosine.angle_to(sine)
    };
    let mut primitives = Primitives {
        builder: k,
        angle,
        cosine,
        sine,
        known: BTreeMap::new(),
    };
    let band_flux = primitives.integral(&reach_out.times(&density));
    let band_area = primitives.integral(&density);

    // Each wall's angle, and its height, on the rail and on the input loop.
    let narrow_turn = |p: &Path<'a>, rho: Expr<'a>| narrow_turn(p, gap, rho);
    let up = |p: &Path<'a>| (p.y - c0 * p.s) / s0;
    let loops = [(&centre, big, small), (&edge, wide, narrow)];
    let wide_area = difference(loops, |p, rho, _| p.z * wide_turn(p, rho));
    let narrow_area = difference(loops, |p, _, rho| up(p) * narrow_turn(p, rho));
    let narrow_flux = difference(loops, |p, _, rho| {
        (gap * (p.x - gap) / rho + narrow) * narrow * up(p) * narrow_turn(p, rho)
    });

    // Smooth upper bounds for the slack, by `|v| <= (1 + |v|^2) / 2`.
    let half = |squared| (k.k(1) + squared) / two;
    let turning = ball * angle * half(turn_first.square() + turn_second.square());
    let rise = (dy - c0 * ds) / s0;
    let length = half((wide / big).square() * (dx * dx + dy * dy) + dz * dz)
        + half((narrow / small).square() * (dx * dx + ds * ds) + rise * rise);

    let mut outputs = vec![
        band_flux,
        band_area,
        wide_area,
        narrow_area,
        narrow_flux,
        turning,
        length,
    ];
    if moments {
        for i in 0..3 {
            let position = Section::linear(at[i], ball * first[i], ball * second[i]);
            let ray = Section::linear(zero, first[i], second[i]);
            let flux = position.times(&position).times(&ray).times(&density);
            outputs.push(-(primitives.integral(&flux) / two));
        }
        outputs.push(difference(loops, |p, rho, _| {
            p.z * wide.square() * cube(p.x / rho) / two * wide_turn(p, rho)
        }));
        outputs.push(difference(loops, |p, rho, _| {
            p.z * wide.square() * cube(p.y / rho) / two * wide_turn(p, rho)
        }));
        outputs.push(difference(loops, |p, _, rho| {
            let (cos_t, h) = ((p.x - gap) / rho, up(p));
            let across = gap + narrow * cos_t;
            narrow * across * across * h / two * cos_t * narrow_turn(p, rho)
        }));
        // Along and up the wide axis, a narrow point is `a + climb h`.
        for (lean, climb) in [(c0, s0), (-s0, c0)] {
            outputs.push(difference(loops, |p, _, rho| {
                let (sin_t, h) = (p.s / rho, up(p));
                let a = narrow * sin_t * lean;
                let primitive = a * a * h + a * climb * h * h + climb * climb * h * h * h / three;
                narrow * primitive / two * sin_t * lean * narrow_turn(p, rho)
            }));
        }
    }
    // Reflection x -> -x leaves a crossing band fixed. These three first
    // moment densities are odd; all other outputs are even. Their complete
    // integrals on [-1/2, 1/2] vanish exactly, including at the collapsed ends.
    if symmetry && crossing.is_some() && moments {
        for index in [7, 10, 12] {
            outputs[index] = zero;
        }
    }
    // An exact carrier has no tube displacement. Its curvature and boundary
    // length cannot contribute to that zero tube, so discard their DAGs.
    if !slack {
        outputs[5] = zero;
        outputs[6] = zero;
    }
    k.finish(&outputs)
}

#[cfg(test)]
#[path = "loops/tests.rs"]
mod tests;
