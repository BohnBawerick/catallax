//! Local classification over every open cell of a central plane arrangement.
use crate::{Location, geometry};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::Vector3;

pub(crate) fn classify(
    planes: &[Vector3],
    dec: &mut Decider,
    mut ray: impl FnMut(&Vector3, &mut Decider) -> Result<Location, Refusal>,
) -> Result<Option<Location>, Refusal> {
    use Location::*;
    let Some(n) = planes.first() else {
        return Ok(None);
    };
    if planes.len() == 1 {
        return Ok(
            match (ray(n, dec)?, ray(&n.scaled(&Real::from_i64(-1)), dec)?) {
                (Inside, Inside) => Some(Inside),
                (Outside, Outside) => Some(Outside),
                (Inside, Outside) | (Outside, Inside) => Some(Boundary),
                _ => None,
            },
        );
    }
    let mut seen = 0;
    // Every spherical chamber has a vertex where two nonparallel planes meet.
    // At each vertex, each incident plane supplies a transverse line. Perturb
    // to both sides of both ends of that line to visit its adjacent sectors.
    for (i, a) in planes.iter().enumerate() {
        for b in &planes[i + 1..] {
            let v = a.cross(b);
            if geometry::sign(&v.norm_squared(), dec)? == Sign::Zero {
                continue;
            }
            for n in planes {
                if geometry::sign(&n.dot(&v), dec)? != Sign::Zero {
                    continue;
                }
                let w = n.cross(&v);
                for sv in [-1, 1] {
                    let v = v.scaled(&Real::from_i64(sv));
                    for sw in [-1, 1] {
                        let w = w.scaled(&Real::from_i64(sw));
                        for sn in [-1, 1] {
                            let n = n.scaled(&Real::from_i64(sn));
                            let d = perturb(planes, [&v, &w, &n], dec)?;
                            seen |= match ray(&d, dec)? {
                                Inside => 1,
                                Outside => 2,
                                _ => return Ok(None),
                            };
                            if seen == 3 {
                                return Ok(Some(Boundary));
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(match seen {
        1 => Some(Inside),
        2 => Some(Outside),
        _ => None,
    })
}

fn perturb(
    planes: &[Vector3],
    basis: [&Vector3; 3],
    dec: &mut Decider,
) -> Result<Vector3, Refusal> {
    let mut epsilon = Rational::one();
    for plane in planes {
        let c = basis.map(|v| plane.dot(v));
        let mut first = None;
        for (i, value) in c.iter().enumerate() {
            let sign = geometry::sign(value, dec)?;
            if sign != Sign::Zero {
                first = Some((
                    i,
                    if sign == Sign::Negative {
                        -value
                    } else {
                        value.clone()
                    },
                ));
                break;
            }
        }
        let (i, magnitude) = first.expect("the perturbation basis spans three dimensions");
        let mut denominator = &magnitude + &Real::one();
        for tail in &c[i + 1..] {
            denominator = denominator
                + if let Some(q) = tail.rational_ref() {
                    Real::from(q.abs())
                } else {
                    tail * tail + Real::one() // |x| <= x^2 + 1.
                };
        }
        // With epsilon <= this bound and epsilon <= 1, the tail cannot change
        // the sign of the first nonzero coefficient of n.(v + e*w + e^2*u).
        let bound = magnitude / denominator;
        let mut precision = crate::metric::PRECISION;
        let lower = loop {
            let lower = geometry::enclose(&bound, precision, dec)?.lo();
            if lower.sign() == Sign::Positive {
                break lower;
            }
            precision = precision.checked_mul(2).ok_or_else(|| {
                Refusal::new(
                    "numbers.refinement-budget",
                    "planar cell separation exceeds the precision index",
                )
            })?;
        };
        epsilon = epsilon.min(lower);
    }
    let e = Real::from(epsilon);
    Ok(basis[0]
        .plus(&basis[1].scaled(&e))
        .plus(&basis[2].scaled(&(&e * &e))))
}
