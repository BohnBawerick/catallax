//! Definition-bound flux, area and moments of complete rational patch faces.

use crate::integrate::{Builder, Expr, NativeArea, Program, Stop, Work};
use crate::surface::Property;
use catallax_numbers::{Ball, Decider, Integer, Rational, Real, RealPoly, Refusal};
use catallax_shapes::{RPoint3, RationalPatch};
use catallax_topology::{FaceId, Model, RationalPatchFace, Sense};
use std::cell::RefCell;
use std::collections::{BTreeMap, VecDeque};
use std::sync::Arc;

pub(crate) struct Band {
    spans: Vec<PatchProgram>,
    fitted: bool,
    best: Option<Arc<[Ball]>>,
    cache: Option<FaceKey>,
}

#[derive(Clone)]
struct FaceKey {
    patch: RationalPatch,
    anchor: RPoint3,
    property: Property,
    sense: Sense,
}

struct Contribution {
    face: FaceKey,
    precision: u32,
    depth: u32,
    cost: u64,
    checkpoints: Option<u64>,
    values: Arc<[Ball]>,
}

const CACHE_CAPACITY: usize = 64;

thread_local! {
    static CONTRIBUTIONS: RefCell<VecDeque<Contribution>> = const { RefCell::new(VecDeque::new()) };
}

impl FaceKey {
    fn lookup(&self, depth: u32, work: &mut Work) -> Result<Option<Arc<[Ball]>>, Stop> {
        let budgeted = catallax_numbers::operation_budget_active();
        let found = CONTRIBUTIONS.with_borrow(|cache| {
            cache
                .iter()
                .rev()
                .find(|entry| {
                    entry.precision == work.precision
                        && entry.depth == depth
                        && (!budgeted || entry.checkpoints.is_some())
                        && self.property == entry.face.property
                        && self.sense == entry.face.sense
                        && self.patch.shares_definition_nodes(&entry.face.patch)
                        && self.anchor == entry.face.anchor
                })
                .map(|entry| (entry.cost, entry.checkpoints, entry.values.clone()))
        });
        if let Some((cost, checkpoints, values)) = found
            && work.reuse(cost, checkpoints)?
        {
            return Ok(Some(values));
        }
        Ok(None)
    }

    fn remember(
        &self,
        depth: u32,
        precision: u32,
        cost: u64,
        checkpoints: Option<u64>,
        values: Arc<[Ball]>,
    ) {
        CONTRIBUTIONS.with_borrow_mut(|cache| {
            if cache.len() == CACHE_CAPACITY {
                cache.pop_front();
            }
            cache.push_back(Contribution {
                face: self.clone(),
                precision,
                depth,
                cost,
                checkpoints,
                values,
            });
        });
    }
}

enum PatchProgram {
    Tensor(Program, Option<NativeArea>),
    Marginal(Program),
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

fn casteljau(mut values: Vec<[RealPoly; 4]>, axis: usize, centered: bool) -> [RealPoly; 4] {
    let variable = RealPoly::var(2, axis);
    let variable = if centered {
        variable + RealPoly::constant(2, Real::from_rational(Rational::one().mul_pow2(-1)))
    } else {
        variable
    };
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

fn unit_equations(
    controls: &[[Real; 4]],
    columns: usize,
    anchor: &[Real; 3],
    centered: bool,
) -> [RealPoly; 4] {
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
            casteljau(values, 1, centered)
        })
        .collect();
    casteljau(rows, 0, centered)
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

fn horner<'a>(k: &'a Builder, p: &RealPoly, parameters: [Expr<'a>; 2]) -> Expr<'a> {
    let mut rows = BTreeMap::<u32, Vec<(u32, Real)>>::new();
    for (power, value) in p.terms() {
        let [u, v] = <[u32; 2]>::try_from(power.exponents()).expect("two parameters");
        rows.entry(u).or_default().push((v, value));
    }
    let fold = |terms: Vec<(u32, Expr<'a>)>, parameter| {
        let mut result = k.k(0);
        let mut previous = terms.last().map_or(0, |(power, _)| *power);
        for (degree, coefficient) in terms.into_iter().rev() {
            result = result * power(parameter, previous - degree, k.k(1)) + coefficient;
            previous = degree;
        }
        result * power(parameter, previous, k.k(1))
    };
    let rows = rows
        .into_iter()
        .map(|(u, mut row)| {
            row.sort_by_key(|(v, _)| *v);
            let terms = row.into_iter().map(|(v, c)| (v, k.exact(&c))).collect();
            (u, fold(terms, parameters[1]))
        })
        .collect();
    fold(rows, parameters[0])
}

fn integrate_travel(poly: RealPoly) -> RealPoly {
    // The centered travel coordinate runs from -1/2 to 1/2. Integrate
    // each monomial in that coordinate while holding the profile fixed.
    let terms = poly.terms().into_iter().filter_map(|(p, c)| {
        let [u, v] = <[u32; 2]>::try_from(p.exponents()).expect("two parameters");
        (v % 2 == 0).then(|| {
            let factor = Rational::from_ratio(1, i64::from(v) + 1)
                .unwrap()
                .mul_pow2(-i64::from(v));
            (vec![u, 0], c * Real::from_rational(factor))
        })
    });
    RealPoly::from_terms(2, terms.collect())
}

// All carried copies of a profile span have the same weight. Sum their
// signed polynomial moments before quadrature so small centroid coordinates
// retain the cancellation across faces.
#[derive(Default)]
pub(crate) struct Marginals(BTreeMap<Vec<(u32, Rational)>, Vec<RealPoly>>);

impl Marginals {
    fn add(
        &mut self,
        weight: &RealPoly,
        coordinates: &[RealPoly; 4],
        normal: &[RealPoly; 3],
        property: Property,
        sense: Sense,
    ) -> bool {
        let Some(weight) = weight.rational() else {
            return false;
        };
        let mut key = Vec::new();
        for (power, coefficient) in weight.terms() {
            if power.exponents()[1] != 0 {
                return false;
            }
            key.push((power.exponents()[0], coefficient));
        }
        let orientation = Real::from_i64(if sense == Sense::Same { 1 } else { -1 });
        let flux = (0..3).fold(RealPoly::constant(2, Real::zero()), |sum, i| {
            sum + &coordinates[i] * &normal[i]
        });
        let mut numerators = vec![integrate_travel(flux).scale(&orientation)];
        if property == Property::Centroid {
            for i in 0..3 {
                let square = &coordinates[i] * &coordinates[i];
                numerators.push(integrate_travel(&square * &normal[i]).scale(&orientation));
            }
        }
        match self.0.entry(key) {
            std::collections::btree_map::Entry::Vacant(entry) => {
                entry.insert(numerators);
            }
            std::collections::btree_map::Entry::Occupied(mut entry) => {
                for (sum, term) in entry.get_mut().iter_mut().zip(numerators) {
                    *sum = &*sum + &term;
                }
            }
        }
        true
    }

    pub(crate) fn finish(self) -> Option<Band> {
        let spans: Vec<_> = self
            .0
            .into_iter()
            .map(|(key, numerators)| {
                let weight = RealPoly::from_terms(
                    2,
                    key.into_iter()
                        .map(|(power, value)| (vec![power, 0], Real::from_rational(value)))
                        .collect(),
                );
                let k = Builder::default();
                let parameters = [k.input(0), k.input(1)];
                let inverse = k.k(1) / horner(&k, &weight, parameters);
                let inv4 = inverse.square().square();
                let outputs: Vec<_> = numerators
                    .iter()
                    .enumerate()
                    .map(|(i, numerator)| {
                        let form = horner(&k, numerator, parameters) * inv4;
                        if i == 0 {
                            form / k.k(3)
                        } else {
                            form * inverse / k.k(2)
                        }
                    })
                    .collect();
                PatchProgram::Marginal(k.finish(&outputs))
            })
            .collect();
        (!spans.is_empty()).then_some(Band {
            spans,
            fitted: true,
            best: None,
            cache: None,
        })
    }
}

fn fitted_forms(
    coordinates: &[RealPoly; 4],
    property: Property,
    sense: Sense,
    marginals: &mut Marginals,
    native_supported: bool,
) -> Option<PatchProgram> {
    let weight_poly = &coordinates[3];
    let normal = normal_numerator(coordinates);
    if matches!(property, Property::Volume | Property::Centroid)
        && marginals.add(weight_poly, coordinates, &normal, property, sense)
    {
        return None;
    }
    if property != Property::Area {
        return Some(PatchProgram::Tensor(
            forms(coordinates, property, sense),
            None,
        ));
    }
    let k = Builder::default();
    let parameters = [k.input(0), k.input(1)];
    let inverse = k.k(1) / horner(&k, weight_poly, parameters);
    let square = normal
        .iter()
        .fold(RealPoly::constant(2, Real::zero()), |sum, n| sum + n * n);
    let area = horner(&k, &square, parameters).sqrt() * inverse.square() * inverse;
    let native = if native_supported {
        let cube = &(weight_poly * weight_poly) * weight_poly;
        NativeArea::new(&square, &cube)
    } else {
        None
    };
    Some(PatchProgram::Tensor(k.finish(&[area]), native))
}

pub(crate) fn prepare(
    model: &Model,
    face: FaceId,
    anchor: &RPoint3,
    property: Property,
    fitted: bool,
    marginals: &mut Marginals,
    dec: &mut Decider,
) -> Result<Option<Band>, Refusal> {
    let checked = RationalPatchFace::of_face(model, face, dec)?.ok_or_else(|| {
        Refusal::new(
            "measure.rational-patch-domain",
            "a rational face has no checked complete domain",
        )
    })?;
    let degrees = checked.support().carrier().degrees();
    let columns = degrees[1] + 1;
    // Bounded polynomial degrees keep a binary64 panel's working set small.
    let native_supported = degrees.iter().all(|&degree| degree <= 3);
    let sense = model.face(face).sense;
    // Check the complete face domain even when a contribution is cached.
    let coordinates = [&anchor.x, &anchor.y, &anchor.z].map(|r| Real::from_rational(r.clone()));
    let spans: Vec<_> = checked
        .support()
        .bezier_spans(dec)?
        .filter_map(|(_, controls)| {
            let coordinates = unit_equations(controls, columns, &coordinates, fitted);
            if fitted {
                fitted_forms(&coordinates, property, sense, marginals, native_supported)
            } else {
                Some(PatchProgram::Tensor(
                    forms(&coordinates, property, sense),
                    None,
                ))
            }
        })
        .collect();
    let cache = (!fitted
        && spans.iter().all(
            |span| matches!(span, PatchProgram::Tensor(program, None) if program.is_rational()),
        ))
    .then(|| FaceKey {
        patch: checked.support().clone(),
        anchor: anchor.clone(),
        property,
        sense,
    });
    Ok((!spans.is_empty()).then_some(Band {
        spans,
        fitted,
        best: None,
        cache,
    }))
}

impl Band {
    pub(crate) fn has_native_area(&self) -> bool {
        self.spans
            .iter()
            .any(|program| matches!(program, PatchProgram::Tensor(_, Some(_))))
    }

    pub(crate) fn integrate(
        &mut self,
        rest: &[Ball],
        depth: u32,
        rectangle: bool,
        work: &mut Work,
    ) -> Result<Vec<Ball>, Stop> {
        if !self.fitted
            && !rectangle
            && let Some(best) = &self.best
            && best.iter().zip(rest).all(|(own, other)| {
                (own.hi() - own.lo()) * Rational::from_i64(8) <= other.hi() - other.lo()
            })
        {
            return Ok(best.to_vec());
        }
        let marginal = self
            .spans
            .iter()
            .all(|span| matches!(span, PatchProgram::Marginal(_)));
        let native = work.native_area && self.has_native_area();
        let depths = if !rectangle && native {
            work.tensor_depth
        } else if native {
            // A broad range of the rational profile is the first domain
            // obstacle. Subsequent Gauss splits use measured remainders.
            [depth, depth.saturating_sub(1)]
        } else if marginal {
            [depth, 0]
        } else {
            [depth; 2]
        };
        let counts = [
            1usize.checked_shl(depths[0]).ok_or(Stop::Budget)?,
            1usize.checked_shl(depths[1]).ok_or(Stop::Budget)?,
        ];
        if counts[0]
            .checked_mul(counts[1])
            .and_then(|cells| cells.checked_mul(self.spans.len()))
            .is_none_or(|cells| cells as u64 > work.remaining())
        {
            return Err(Stop::Budget);
        }
        let rule = if rectangle || marginal {
            None
        } else {
            Some(work.gauss_rule()?)
        };
        // Obtain the rule even on a hit: its proof and its charge belong to
        // this attempt and it may be needed by other faces in the same sum.
        let cache = self.cache.as_ref().filter(|_| {
            !rectangle && work.dec.deadline().is_none() && catallax_numbers::shared_caches()
        });
        if let Some(face) = cache
            && let Some(values) = face.lookup(depth, work)?
        {
            let result = values.to_vec();
            self.best = Some(values);
            return Ok(result);
        }
        let remaining = work.remaining();
        // The only charges to the surrounding allowance in these rational
        // programs are single-unit deadline checks. Keep those too, so a warm
        // contribution does not turn an exhausted request into an answer.
        let before = (cache.is_some() && catallax_numbers::operation_budget_active())
            .then(catallax_numbers::current_operation_work);
        let axis = |index| -> Vec<_> {
            (0..=counts[index])
                .map(|i| {
                    let value = Rational::from_integer(Integer::from_u64(i as u64))
                        .mul_pow2(-i64::from(depths[index]));
                    if self.fitted {
                        value - Rational::one().mul_pow2(-1)
                    } else {
                        value
                    }
                })
                .collect()
        };
        let profile = axis(0);
        let travel = (depths[0] != depths[1]).then(|| axis(1));
        let travel = travel.as_ref().unwrap_or(&profile);
        let mut sums: Option<Vec<Ball>> = None;
        for program in &self.spans {
            for u in profile.windows(2) {
                match program {
                    PatchProgram::Marginal(program) => {
                        let values = program.parameter_panel(&u[0], &u[1], rectangle, work)?;
                        accumulate(&mut sums, values, work)?;
                    }
                    PatchProgram::Tensor(program, area) => {
                        for v in travel.windows(2) {
                            let u = [&u[0], &u[1]];
                            let v = [&v[0], &v[1]];
                            let values = match (&rule, area) {
                                (Some(rule), Some(area)) if native => {
                                    area.integrate(u, v, rule, work)?
                                }
                                (Some(rule), _) => {
                                    program.parameter_tensor_panel(u, v, rule, self.fitted, work)?
                                }
                                (None, _) => program.parameter_tensor_rectangle(u, v, work)?,
                            };
                            accumulate(&mut sums, values, work)?;
                        }
                    }
                }
            }
        }
        let sums = sums.ok_or(Stop::Budget)?;
        if !self.fitted && !rectangle {
            let values = Arc::from(sums.as_slice());
            if let Some(face) = cache {
                let checkpoints = before.and_then(|before| {
                    let after = catallax_numbers::current_operation_work();
                    (!before.exhausted && !after.exhausted).then(|| after.used - before.used)
                });
                if before.is_none() || checkpoints.is_some() {
                    face.remember(
                        depth,
                        work.precision,
                        remaining - work.remaining(),
                        checkpoints,
                        Arc::clone(&values),
                    );
                }
            }
            self.best = Some(values);
        }
        Ok(sums)
    }
}

fn accumulate(
    sums: &mut Option<Vec<Ball>>,
    values: Vec<Ball>,
    work: &mut Work,
) -> Result<(), Stop> {
    *sums = Some(match sums.take() {
        None => values,
        Some(previous) => previous
            .iter()
            .zip(values)
            .map(|(sum, value)| work.add(sum, &value))
            .collect::<Result<_, _>>()?,
    });
    Ok(())
}
