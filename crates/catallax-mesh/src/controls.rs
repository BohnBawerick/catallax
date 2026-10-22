//! Certificates refer to the rounded coordinates actually returned to callers.

use crate::{
    chart::{self, BITS, Chart, Uv},
    failure,
};
use catallax_numbers::{Ball, Decider, Rational, Refusal};
use catallax_shapes::{Surface, Vector3};

fn dot(a: &[Ball; 3], b: &[Ball; 3]) -> Ball {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}
fn unit(v: &Vector3, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    let norm = v.norm_squared().sqrt();
    let values = [&v.x, &v.y, &v.z].map(|x| dec.enclose_at(&(x / &norm), BITS));
    let [x, y, z] = values;
    Ok([x?, y?, z?])
}
fn normal(chart: &Chart, uv: Uv, dec: &Decider) -> Result<[Ball; 3], Refusal> {
    let (axis, reference) = match &chart.surface {
        Surface::Plane(p) => return unit(&p.normal, dec),
        Surface::Cylinder(c) => (&c.axis, &c.ref_dir),
        Surface::Cone(c) => (&c.axis, &c.ref_dir),
        Surface::Sphere(s) => (&s.axis, &s.ref_dir),
        Surface::Torus(t) => (&t.axis, &t.ref_dir),
        _ => return Err(failure("surface has no normal certificate")),
    };
    let axis_unit = unit(axis, dec)?;
    let first = unit(reference, dec)?;
    let second = unit(&axis.cross(reference), dec)?;
    let tau = &Ball::pi_checked(BITS, dec)? * &Ball::from_rational(&Rational::from_i64(2), BITS);
    let angle = |i| &Ball::from_rational(&chart::q(uv[i]), BITS) * &tau;
    let u = angle(0);
    let cosine = u.cos();
    let sine = u.sin();
    let radial: [Ball; 3] = std::array::from_fn(|i| &(&first[i] * &cosine) + &(&second[i] * &sine));
    match &chart.surface {
        Surface::Cylinder(_) => Ok(radial),
        Surface::Cone(c) => {
            let slope = dec.enclose_at(&c.slope, BITS)?;
            let norm = dec.enclose_at(
                &(catallax_numbers::Real::from(1) + c.slope.square()).sqrt(),
                BITS,
            )?;
            let values = std::array::from_fn(|i| {
                (&radial[i] - &(&axis_unit[i] * &slope))
                    .checked_div(&norm)
                    .unwrap()
            });
            Ok(values)
        }
        Surface::Sphere(_) | Surface::Torus(_) => {
            let v = angle(1);
            let cosine = v.cos();
            let sine = v.sin();
            Ok(std::array::from_fn(|i| {
                &(&radial[i] * &cosine) + &(&axis_unit[i] * &sine)
            }))
        }
        _ => unreachable!(),
    }
}

pub(crate) fn length_squared(a: [f64; 3], b: [f64; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| {
        let d = chart::q(a[i]) - chart::q(b[i]);
        sum + &d * &d
    })
}

pub(crate) fn cosine_limit(angle: &Rational, dec: &Decider) -> Result<Rational, Refusal> {
    if *angle >= Ball::pi_checked(BITS, dec)?.hi() {
        Ok(Rational::from_i64(-1))
    } else {
        Ok(Ball::from_rational(angle, BITS).cos().hi())
    }
}

pub(crate) fn quality(
    chart: &Chart,
    uvs: &[Uv; 3],
    positions: [[f64; 3]; 3],
    angular_limit: Option<&Rational>,
    dec: &Decider,
) -> Result<(Rational, bool), Refusal> {
    let [a, b, c] = positions.map(|p| p.map(chart::q));
    let ab: [Rational; 3] = std::array::from_fn(|i| &b[i] - &a[i]);
    let ac: [Rational; 3] = std::array::from_fn(|i| &c[i] - &a[i]);
    let cross = [
        &ab[1] * &ac[2] - &ab[2] * &ac[1],
        &ab[2] * &ac[0] - &ab[0] * &ac[2],
        &ab[0] * &ac[1] - &ab[1] * &ac[0],
    ];
    let norm2 = cross.iter().fold(Rational::zero(), |s, x| s + x * x);
    if norm2.is_zero() {
        return Err(failure("display triangle has zero area"));
    }
    let norm = Ball::from_rational(&norm2, BITS).sqrt().unwrap();
    if norm.lo() <= Rational::zero() {
        return Err(failure("triangle normal exceeds certificate precision"));
    }
    let longest = (0..3)
        .map(|i| length_squared(positions[i], positions[(i + 1) % 3]))
        .max()
        .unwrap();
    let aspect = &longest / &norm.lo();
    let angular_reached = if let Some(limit) = angular_limit {
        let center = std::array::from_fn(|i| {
            let lo = uvs.iter().map(|p| p[i]).fold(f64::INFINITY, f64::min);
            let hi = uvs.iter().map(|p| p[i]).fold(f64::NEG_INFINITY, f64::max);
            lo * 0.5 + hi * 0.5
        });
        let spans: [Rational; 2] = std::array::from_fn(|i| {
            uvs.iter()
                .map(|p| (chart::q(p[i]) - chart::q(center[i])).abs())
                .max()
                .unwrap()
        });
        let parameter_span = match chart.surface {
            Surface::Plane(_) => Rational::zero(),
            Surface::Cylinder(_) | Surface::Cone(_) => spans[0].clone(),
            Surface::Sphere(_) | Surface::Torus(_) => &spans[0] + &spans[1],
            _ => return Err(failure("surface has no normal variation certificate")),
        };
        // Each angular parameter rotates a unit normal. The spherical triangle
        // inequality bounds its variation by the sum of the rotation magnitudes.
        let variation =
            parameter_span * (Ball::pi_checked(BITS, dec)?.hi() * Rational::from_i64(2));
        let margin = limit - &variation;
        if margin <= Rational::zero() {
            false
        } else {
            let facet = cross
                .each_ref()
                .map(|x| Ball::from_rational(x, BITS).checked_div(&norm).unwrap());
            let cosine = dot(&facet, &normal(chart, center, dec)?)
                .lo()
                .max(Rational::from_i64(-1))
                .min(Rational::one());
            cosine >= cosine_limit(&margin, dec)?
        }
    } else {
        true
    };
    Ok((aspect, angular_reached))
}
