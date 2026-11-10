//! Definition-bound flux, area and moments of complete rational patch faces.

use crate::integrate::{Builder, Expr, Program, Stop, Work};
use crate::surface::Property;
use catallax_numbers::{Ball, Decider, Integer, Rational, Real, RealPoly, Refusal};
use catallax_shapes::RPoint3;
use catallax_topology::{FaceId, Model, RationalPatchFace, Sense};

pub(crate) struct Band {
    spans: Vec<Program>,
    best: Option<Vec<Ball>>,
}

fn power<'a>(mut value: Expr<'a>, mut exponent: u32, one: Expr<'a>) -> Expr<'a> {
    let mut result = one;
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = result * value;
        }
        exponent >>= 1;
        if exponent != 0 {
            value = value.square();
        }
    }
    result
}

fn polynomial<'a>(k: &'a Builder, p: &RealPoly, parameters: [Expr<'a>; 2]) -> Expr<'a> {
    let one = k.k(1);
    p.terms()
        .into_iter()
        .fold(k.k(0), |sum, (monomial, coefficient)| {
            sum + k.exact(&coefficient)
                * power(parameters[0], monomial.exponents()[0], one)
                * power(parameters[1], monomial.exponents()[1], one)
        })
}

fn casteljau(mut values: Vec<[RealPoly; 4]>, axis: usize) -> [RealPoly; 4] {
    let variable = RealPoly::var(2, axis);
    let complement = &RealPoly::constant(2, Real::one()) - &variable;
    for length in (1..values.len()).rev() {
        for i in 0..length {
            let (before, after) = values.split_at_mut(i + 1);
            for (value, next) in before[i].iter_mut().zip(&after[0]) {
                *value = &complement * &*value + &variable * next;
            }
        }
    }
    values.swap_remove(0)
}

fn unit_equations(controls: &[[Real; 4]], columns: usize, anchor: &[Real; 3]) -> [RealPoly; 4] {
    let rows = controls
        .chunks(columns)
        .map(|row| {
            let values = row
                .iter()
                .map(|h| {
                    std::array::from_fn(|axis| {
                        let value = if axis == 3 {
                            h[3].clone()
                        } else {
                            &h[axis] - &(&anchor[axis] * &h[3])
                        };
                        RealPoly::constant(2, value)
                    })
                })
                .collect();
            casteljau(values, 1)
        })
        .collect();
    casteljau(rows, 0)
}

fn cross_polynomial(a: [&RealPoly; 3], b: [&RealPoly; 3]) -> [RealPoly; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normal_numerator(h: &[RealPoly; 4]) -> [RealPoly; 3] {
    let du = h.each_ref().map(|p| p.partial(0));
    let dv = h.each_ref().map(|p| p.partial(1));
    let xyz = std::array::from_fn(|i| &h[i]);
    let u = std::array::from_fn(|i| &du[i]);
    let v = std::array::from_fn(|i| &dv[i]);
    let uv = cross_polynomial(u, v);
    let uh = cross_polynomial(u, xyz);
    let hv = cross_polynomial(xyz, v);
    // (W Hu - H Wu) cross (W Hv - H Wv) has a known factor W.
    // Cancel it symbolically before interval/Taylor evaluation; W is positive.
    std::array::from_fn(|i| &h[3] * &uv[i] - &dv[3] * &uh[i] - &du[3] * &hv[i])
}

fn oriented_flux<'a>(
    k: &'a Builder,
    position: [Expr<'a>; 3],
    normal: [Expr<'a>; 3],
    sense: Sense,
) -> Expr<'a> {
    (0..3).fold(k.k(0), |sum, i| sum + position[i] * normal[i])
        * k.k(if sense == Sense::Same { 1 } else { -1 })
}

fn forms(coordinates: &[RealPoly; 4], property: Property, sense: Sense) -> Program {
    let k = Builder::default();
    let zero = k.k(0);
    let parameters = [k.input(0), k.input(1)];
    let weight = polynomial(&k, &coordinates[3], parameters);
    let denominator = weight.square() * weight;
    let normal = normal_numerator(coordinates)
        .each_ref()
        .map(|p| polynomial(&k, p, parameters) / denominator);
    let position = || std::array::from_fn(|i| polynomial(&k, &coordinates[i], parameters) / weight);
    let outputs = match property {
        Property::Volume => vec![oriented_flux(&k, position(), normal, sense) / k.k(3)],
        Property::Area | Property::FaceCentroid => {
            let area = normal.iter().fold(zero, |sum, x| sum + x.square()).sqrt();
            if matches!(property, Property::Area) {
                vec![area]
            } else {
                let position = position();
                vec![
                    area,
                    position[0] * area,
                    position[1] * area,
                    position[2] * area,
                ]
            }
        }
        Property::Centroid => {
            let position = position();
            let orientation = k.k(if sense == Sense::Same { 1 } else { -1 });
            vec![
                oriented_flux(&k, position, normal, sense) / k.k(3),
                position[0].square() * normal[0] * orientation / k.k(2),
                position[1].square() * normal[1] * orientation / k.k(2),
                position[2].square() * normal[2] * orientation / k.k(2),
            ]
        }
        Property::SecondMoment(i, j) => {
            let position = position();
            vec![position[i] * position[j] * oriented_flux(&k, position, normal, sense) / k.k(5)]
        }
    };
    k.finish(&outputs)
}

pub(crate) fn prepare(
    model: &Model,
    face: FaceId,
    anchor: &RPoint3,
    property: Property,
    dec: &mut Decider,
) -> Result<Band, Refusal> {
    let checked = RationalPatchFace::of_face(model, face, dec)?.ok_or_else(|| {
        Refusal::new(
            "measure.rational-patch-domain",
            "a rational face has no checked complete domain",
        )
    })?;
    let columns = checked.support().carrier().degrees()[1] + 1;
    let anchor = [&anchor.x, &anchor.y, &anchor.z].map(|value| Real::from_rational(value.clone()));
    let spans = checked
        .support()
        .bezier_spans(dec)?
        .map(|(_, controls)| {
            forms(
                &unit_equations(controls, columns, &anchor),
                property,
                model.face(face).sense,
            )
        })
        .collect();
    Ok(Band { spans, best: None })
}

impl Band {
    pub(crate) fn integrate(
        &mut self,
        rest: &[Ball],
        depth: u32,
        rectangle: bool,
        work: &mut Work,
    ) -> Result<Vec<Ball>, Stop> {
        if !rectangle
            && let Some(best) = &self.best
            && best.iter().zip(rest).all(|(own, other)| {
                (own.hi() - own.lo()) * Rational::from_i64(8) <= other.hi() - other.lo()
            })
        {
            return Ok(best.clone());
        }
        let count = 1usize.checked_shl(depth).ok_or(Stop::Budget)?;
        if count
            .checked_mul(count)
            .and_then(|cells| cells.checked_mul(self.spans.len()))
            .is_none_or(|cells| cells as u64 > work.remaining())
        {
            return Err(Stop::Budget);
        }
        let rule = if rectangle {
            None
        } else {
            Some(work.gauss_rule()?)
        };
        let axis: Vec<_> = (0..=count)
            .map(|i| {
                Rational::from_integer(Integer::from_u64(i as u64)).mul_pow2(-i64::from(depth))
            })
            .collect();
        let mut sums: Option<Vec<Ball>> = None;
        for program in &self.spans {
            for u in axis.windows(2) {
                for v in axis.windows(2) {
                    let values = if let Some(rule) = &rule {
                        program.parameter_tensor_panel(
                            [&u[0], &u[1]],
                            [&v[0], &v[1]],
                            rule,
                            work,
                        )?
                    } else {
                        program.parameter_tensor_rectangle([&u[0], &u[1]], [&v[0], &v[1]], work)?
                    };
                    sums = Some(match sums {
                        None => values,
                        Some(previous) => previous
                            .iter()
                            .zip(values)
                            .map(|(sum, value)| work.add(sum, &value))
                            .collect::<Result<_, _>>()?,
                    });
                }
            }
        }
        let sums = sums.ok_or(Stop::Budget)?;
        if !rectangle {
            self.best = Some(sums.clone());
        }
        Ok(sums)
    }
}
