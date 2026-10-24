//! Exact rational chart restrictions and Bernstein signs on open trim cells.
//! Chart denominators are products of 1+t², hence strictly positive. Clearing
//! them preserves signs, including exact cancellation on a curved carrier.

use super::{Chart, Kind, q};
use crate::display::Cell;
use catallax_numbers::{Decider, MPoly, Rational, Refusal, Sign, UPoly};
use std::collections::BTreeMap;

fn powers(value: &MPoly, degree: usize, dec: &mut Decider) -> Result<Vec<MPoly>, Refusal> {
    let mut result = vec![MPoly::constant(2, q(1))];
    for i in 0..degree {
        dec.check_deadline()?;
        result.push(&result[i] * value);
    }
    Ok(result)
}

fn circle(axis: usize, quadrant: usize) -> [MPoly; 3] {
    let t = MPoly::var(2, axis);
    let one = MPoly::constant(2, q(1));
    let square = &t * &t;
    let c = &one - &square;
    let s = t.scale(&q(2));
    let denominator = &one + &square;
    let negative = |p: &MPoly| p.scale(&q(-1));
    let [c, s] = match quadrant {
        0 => [c, s],
        1 => [negative(&s), c],
        2 => [negative(&c), negative(&s)],
        _ => [s, negative(&c)],
    };
    [c, s, denominator]
}

fn choose(n: usize, k: usize) -> Rational {
    let mut result = q(1);
    for j in 0..k.min(n - k) {
        result = result * q(i64::try_from(n - j).expect("analytic degree"))
            / q(i64::try_from(j + 1).expect("analytic degree"));
    }
    result
}

// Bernstein coefficients of x^power on [lo,hi], elevated to degree n.
fn basis(n: usize, range: &[Rational; 2]) -> Vec<Vec<Rational>> {
    let span = &range[1] - &range[0];
    let mut lo = vec![q(1)];
    let mut width = vec![q(1)];
    for i in 0..n {
        lo.push(&lo[i] * &range[0]);
        width.push(&width[i] * &span);
    }
    (0..=n)
        .map(|control| {
            (0..=n)
                .map(|power| {
                    (0..=power.min(control)).fold(q(0), |sum, k| {
                        sum + choose(power, k)
                            * lo[power - k].clone()
                            * width[k].clone()
                            * choose(control, k)
                            / choose(n, k)
                    })
                })
                .collect()
        })
        .collect()
}

pub(in crate::display) struct Pullback {
    terms: Vec<([usize; 2], Rational)>,
    degree: [usize; 2],
}
impl Pullback {
    fn new(polynomial: &MPoly) -> Self {
        let terms: Vec<_> = polynomial
            .terms()
            .into_iter()
            .map(|(m, c)| ([m.exponents()[0] as usize, m.exponents()[1] as usize], c))
            .collect();
        let degree =
            std::array::from_fn(|axis| terms.iter().map(|(e, _)| e[axis]).max().unwrap_or(0));
        Self { terms, degree }
    }

    pub(in crate::display) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub(in crate::display) fn sign(
        &self,
        cell: &Cell,
        dec: &mut Decider,
    ) -> Result<Option<Sign>, Refusal> {
        if self.is_zero() {
            return Ok(Some(Sign::Zero));
        }
        if cell.range.iter().any(|r| r[0] >= r[1]) {
            return Ok(None);
        }
        let u = basis(self.degree[0], &cell.range[0]);
        let v = basis(self.degree[1], &cell.range[1]);
        let mut positive = false;
        let mut negative = false;
        for u_control in &u {
            for v_control in &v {
                dec.check_deadline()?;
                let coefficient = self.terms.iter().fold(q(0), |sum, (e, c)| {
                    sum + c * &u_control[e[0]] * v_control[e[1]].clone()
                });
                positive |= coefficient > q(0);
                negative |= coefficient < q(0);
                if positive && negative {
                    return Ok(None);
                }
            }
        }
        // Every tensor-product Bernstein basis function is strictly positive
        // in the open rectangle. Nonnegative coefficients, at least one
        // positive, therefore prove a strict sign even at an aligned trim.
        Ok(match (positive, negative) {
            (true, false) => Some(Sign::Positive),
            (false, true) => Some(Sign::Negative),
            _ => None,
        })
    }

    pub(in crate::display) fn cuts(
        &self,
        cell: &Cell,
        dec: &mut Decider,
    ) -> Result<[Vec<Rational>; 2], Refusal> {
        let mut result = [Vec::new(), Vec::new()];
        for axis in 0..2 {
            let mut coefficients = BTreeMap::<usize, Vec<Rational>>::new();
            for (e, c) in &self.terms {
                let v = coefficients.entry(e[1 - axis]).or_default();
                v.resize(v.len().max(e[axis] + 1), q(0));
                v[e[axis]] = &v[e[axis]] + c;
            }
            let polynomials: Vec<_> = coefficients.into_values().map(UPoly::new).collect();
            // A whole parameter rail must be a root of every coefficient in
            // the other parameter. Use one linear/quadratic coefficient to
            // propose exact rational roots, then verify all coefficients.
            // Higher degree rails retain ordinary certified subdivision.
            let Some(p) = polynomials
                .iter()
                .find(|p| matches!(p.degree(), Some(1 | 2)))
            else {
                continue;
            };
            dec.check_deadline()?;
            let c = p.coefficients();
            let roots = if c.len() == 2 {
                vec![-&c[0] / c[1].clone()]
            } else {
                let discriminant = &c[1] * &c[1] - q(4) * c[2].clone() * c[0].clone();
                if let Some(root) = crate::profile::sqrt(&discriminant) {
                    vec![
                        (-&c[1] - root.clone()) / (&c[2] * &q(2)),
                        (-&c[1] + root) / (&c[2] * &q(2)),
                    ]
                } else {
                    Vec::new()
                }
            };
            for root in roots {
                dec.check_deadline()?;
                if root > cell.range[axis][0]
                    && root < cell.range[axis][1]
                    && polynomials.iter().all(|p| p.eval(&root).is_zero())
                {
                    result[axis].push(root);
                }
            }
        }
        Ok(result)
    }
}

impl Chart {
    fn homogeneous(&self) -> ([MPoly; 3], MPoly) {
        let one = MPoly::constant(2, q(1));
        let zero = MPoly::zero(2);
        let u = MPoly::var(2, 0);
        let v = MPoly::var(2, 1);
        let [c, s, du] = circle(0, self.quadrant[0]);
        let (local, denominator) = match self.kind {
            Kind::Plane => ([u, v, zero], one),
            Kind::Disk => ([&c * &v, &s * &v, zero], du),
            Kind::Cylinder => ([c.scale(&self.radius), s.scale(&self.radius), &v * &du], du),
            Kind::Cone => (
                [
                    (&c * &v).scale(&self.radius),
                    (&s * &v).scale(&self.radius),
                    &v * &du,
                ],
                du,
            ),
            Kind::Sphere => {
                let [cv, sv, dv] = circle(1, 0);
                (
                    [
                        (&c * &cv).scale(&self.radius),
                        (&s * &cv).scale(&self.radius),
                        (&sv * &du).scale(&self.radius),
                    ],
                    &du * &dv,
                )
            }
            Kind::Torus => {
                let [cv, sv, dv] = circle(1, self.quadrant[1]);
                let radial = &dv.scale(&self.radius) + &cv.scale(&self.minor);
                (
                    [&c * &radial, &s * &radial, (&sv * &du).scale(&self.minor)],
                    &du * &dv,
                )
            }
        };
        let world = std::array::from_fn(|i| {
            (0..3).fold(denominator.scale(&self.origin[i]), |sum, j| {
                &sum + &local[j].scale(&self.frame[j][i])
            })
        });
        (world, denominator)
    }

    pub(in crate::display) fn pullback(
        &self,
        polynomial: &MPoly,
        dec: &mut Decider,
    ) -> Result<Pullback, Refusal> {
        let degree = polynomial.total_degree().unwrap_or(0) as usize;
        let (coordinates, denominator) = self.homogeneous();
        let denominator = powers(&denominator, degree, dec)?;
        let mut coordinate_powers = Vec::new();
        for coordinate in &coordinates {
            coordinate_powers.push(powers(coordinate, degree, dec)?);
        }
        let mut result = MPoly::zero(2);
        for (m, coefficient) in polynomial.terms() {
            dec.check_deadline()?;
            let mut term = denominator[degree - m.degree() as usize].scale(&coefficient);
            for axis in 0..3 {
                term = &term * &coordinate_powers[axis][m.exponents()[axis] as usize];
            }
            result = &result + &term;
        }
        Ok(Pullback::new(&result))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::{Cone, Cylinder, Plane, RPoint3, RVector3, Sphere, Surface, Torus};

    fn power(value: &Rational, n: u32) -> Rational {
        (0..n).fold(q(1), |p, _| &p * value)
    }

    #[test]
    fn homogeneous_restrictions_match_independent_chart_points() {
        let origin = RPoint3::new(2.into(), (-3).into(), 4.into()).to_real();
        let axis = RVector3::new(0.into(), 3.into(), 4.into()).to_real();
        let reference = RVector3::new(1.into(), 0.into(), 0.into()).to_real();
        let surfaces = [
            Surface::Plane(Plane {
                origin: origin.clone(),
                normal: axis.clone(),
                ref_dir: reference.clone(),
            }),
            Surface::Cylinder(Cylinder {
                origin: origin.clone(),
                axis: axis.clone(),
                radius_sq: 25.into(),
                ref_dir: reference.clone(),
            }),
            Surface::Cone(Cone {
                apex: origin.clone(),
                axis: axis.clone(),
                slope: Rational::from_ratio(2, 3).unwrap().into(),
                ref_dir: reference.clone(),
            }),
            Surface::Sphere(Sphere {
                center: origin.clone(),
                axis: axis.clone(),
                radius: 5.into(),
                ref_dir: reference.clone(),
            }),
            Surface::Torus(Torus {
                center: origin,
                axis,
                major: 7.into(),
                minor: 2.into(),
                ref_dir: reference,
            }),
        ];
        let base = MPoly::from_terms(
            3,
            vec![
                (vec![2, 0, 0], q(1)),
                (vec![0, 1, 1], q(2)),
                (vec![0, 0, 1], q(3)),
                (vec![0, 0, 0], q(-7)),
            ],
        );
        let polynomial = &base * &base;
        let mut dec = Decider::default();
        for surface in surfaces {
            for chart in Chart::all(&surface, &q(20)).unwrap() {
                assert!(
                    chart
                        .pullback(surface.implicit().unwrap().rational().unwrap(), &mut dec)
                        .unwrap()
                        .is_zero()
                );
                let restriction = chart.pullback(&polynomial, &mut dec).unwrap();
                let (_, denominator) = chart.homogeneous();
                for uv in [
                    [
                        Rational::from_ratio(1, 3).unwrap(),
                        Rational::from_ratio(2, 5).unwrap(),
                    ],
                    [q(0), q(1)],
                ] {
                    let p = chart.point(&uv);
                    let actual = restriction.terms.iter().fold(q(0), |sum, (e, c)| {
                        sum + c * &power(&uv[0], e[0] as u32) * power(&uv[1], e[1] as u32)
                    });
                    let d = denominator.eval(&uv);
                    assert!(d > q(0));
                    assert_eq!(
                        actual,
                        polynomial.eval(&[p.x, p.y, p.z])
                            * power(&d, polynomial.total_degree().unwrap())
                    );
                }
            }
        }
    }

    #[test]
    fn bernstein_proves_open_signs_at_aligned_zero_edges() {
        let p = MPoly::from_terms(
            2,
            vec![(vec![1, 0], q(1)), (vec![2, 0], q(-1)), (vec![0, 1], q(1))],
        );
        let cell = Cell {
            range: [[q(0), q(1)], [q(0), q(1)]],
            depth: 0,
        };
        let mut dec = Decider::default();
        assert_eq!(
            Pullback::new(&p).sign(&cell, &mut dec).unwrap(),
            Some(Sign::Positive)
        );
        assert_eq!(
            Pullback::new(&p.scale(&q(-1)))
                .sign(&cell, &mut dec)
                .unwrap(),
            Some(Sign::Negative)
        );
        let crossing = MPoly::from_terms(
            2,
            vec![
                (vec![1, 0], q(1)),
                (vec![0, 0], Rational::from_ratio(-1, 2).unwrap()),
            ],
        );
        assert_eq!(
            Pullback::new(&crossing).sign(&cell, &mut dec).unwrap(),
            None
        );
        let half = Cell {
            range: [[Rational::from_ratio(1, 2).unwrap(), q(1)], [q(0), q(1)]],
            depth: 1,
        };
        assert_eq!(
            Pullback::new(&crossing).sign(&half, &mut dec).unwrap(),
            Some(Sign::Positive)
        );
    }

    #[test]
    fn rational_rails_are_verified_against_every_other_parameter_coefficient() {
        let a = UPoly::from_roots(&[q(2), q(3)]);
        let b = UPoly::from_roots(&[q(2), q(4)]);
        let polynomial = MPoly::from_terms(
            2,
            a.coefficients()
                .iter()
                .enumerate()
                .map(|(j, c)| (vec![0, j as u32], c.clone()))
                .chain(
                    b.coefficients()
                        .iter()
                        .enumerate()
                        .map(|(j, c)| (vec![1, j as u32], c.clone())),
                )
                .collect(),
        );
        let cell = Cell {
            range: [[q(0), q(1)], [q(0), q(5)]],
            depth: 0,
        };
        let cuts = Pullback::new(&polynomial)
            .cuts(&cell, &mut Decider::default())
            .unwrap();
        assert_eq!(cuts[1], vec![q(2)]);
        assert!(cuts[0].is_empty());
    }
}
