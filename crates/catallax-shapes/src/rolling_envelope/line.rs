//! What the line solver in [`crate::rational_sheet`] needs from the envelope.

use super::*;
use crate::rational_sheet::{Query, Rect, line_hits, rational_dot};
use catallax_numbers::{MPoly, PolySystem, RealPoly, RootBox};

fn unresolved_line(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

impl RollingEnvelope {
    /// Every point where a whole line meets the closed finite carrier, in line
    /// order. Each is an exact transverse crossing and appears once, including
    /// at a knot and on a boundary row. A tangency, or a contact this cannot
    /// separate, refuses. It never reports a hit or a miss from an enclosure.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        line_hits(self, origin, direction, dec)
    }
}

impl LineSheet for RollingEnvelope {
    fn query_bits(&self, dec: &Decider) -> Result<u32, Refusal> {
        RollingEnvelope::query_bits(self, dec)
    }
    fn budgeted_frame(&self, bits: u32, dec: &mut Decider) -> Result<LocalFrame, Refusal> {
        self.data.definition.budgeted_frame(bits, dec)
    }
    fn local_line(&self, origin: &Point3, direction: &Vector3) -> [[Real; 3]; 2] {
        let def = &self.data.definition;
        [
            def.local(origin),
            [
                direction.dot(&def.toward),
                direction.dot(&def.axis),
                direction.dot(&def.across),
            ],
        ]
    }
    fn spans(&self, bits: u32) -> Vec<Patch> {
        Patch::at(
            &self.data.carrier,
            &self.data.net.controls,
            self.data.net.bits,
            DEGREE,
            bits,
        )
    }
    fn chart(&self) -> Rect {
        [
            [Rational::from_i64(-1), Rational::one()],
            [Rational::zero(), Rational::one()],
        ]
    }
    fn knots(&self) -> Vec<Rational> {
        self.data.nodes.iter().map(|node| node.s.clone()).collect()
    }
    fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        RollingEnvelope::on_boundary(self, point, dec)
    }
    fn row_contacts(
        &self,
        high: bool,
        line: &[[Real; 3]; 2],
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        let [o, d] = line;
        if high {
            if math::sign(&d[1], dec)? == Sign::Zero {
                if math::sign(&o[1], dec)? != Sign::Zero {
                    return Ok(Vec::new());
                }
                return Err(unresolved_line(
                    "a line in the mouth plane needs a circle contact solver",
                ));
            }
            return Ok(vec![-(&o[1] / &d[1])]);
        }
        let body = &self.data.definition.body;
        let quadratic = RealPoly::from_terms(
            1,
            vec![
                (vec![0], &(&o[0].square() + &o[2].square()) - &body.square()),
                (
                    vec![1],
                    &Real::from_i64(2) * &(&(&o[0] * &d[0]) + &(&o[2] * &d[2])),
                ),
                (vec![2], &d[0].square() + &d[2].square()),
            ],
        );
        if quadratic.is_zero(dec)? {
            return Err(unresolved_line(
                "a line on the rail cylinder needs a tangency solver",
            ));
        }
        quadratic.real_roots(dec)
    }
    fn knot_section(&self, index: usize) -> [[Real; 4]; 3] {
        let data = &self.data;
        let (pair, column) = if index + 1 < data.nodes.len() {
            (&data.nodes[index..index + 2], 0)
        } else {
            (&data.nodes[index - 1..index + 1], DEGREE)
        };
        let t = data.definition.transition();
        span_rows(&data.definition, pair)
            .map(|row| std::array::from_fn(|k| row[k].at(column, &t.t)))
    }
    // The exact equations of one span in (line parameter, span parameter, v,
    // t, end weight). The carrier's only irrationals are the transition tangent
    // t and the shared end weight, so they enter as variables fixed by their
    // squares, and each root coordinate stays a leaf of a rational system.
    fn span_system(
        &self,
        query: &Query<'_, Self>,
        index: usize,
        rect: &Rect,
        dec: &mut Decider,
    ) -> Result<Option<(Arc<PolySystem>, RootBox)>, Refusal> {
        let data = &self.data;
        let def = &data.definition;
        let t = def.transition();
        let rational = |value: &Real| value.rational_ref().cloned();
        let vector = |value: &[Real; 3]| -> Option<[Rational; 3]> {
            Some([
                rational(&value[0])?,
                rational(&value[1])?,
                rational(&value[2])?,
            ])
        };
        let (Some(origin), Some(direction), Some(first), Some(second)) = (
            vector(&query.exact[0]),
            vector(&query.exact[1]),
            vector(&query.normals[0]),
            vector(&query.normals[1]),
        ) else {
            return Ok(None);
        };
        let Some(tangent_square) = rational(&t.square) else {
            return Ok(None);
        };
        let pair = &data.nodes[index..index + 2];
        let ends = pair.iter().any(|node| node.s.abs() == Rational::one());
        let with_end = |value: i64| -> Vec<Node> {
            pair.iter()
                .map(|node| {
                    let mut node = node.clone();
                    if node.s.abs() == Rational::one() {
                        node.values[WEIGHT] = Real::from_i64(value);
                    }
                    node
                })
                .collect()
        };
        let nvars = if ends { 5 } else { 4 };
        let (Some(base), Some(end)) = (
            span_terms(def, &with_end(0), nvars),
            span_terms(def, &with_end(1), nvars),
        ) else {
            return Ok(None);
        };
        let homogeneous: [MPoly; 4] = std::array::from_fn(|k| {
            if ends {
                let weight = MPoly::var(nvars, 4);
                &base[k] + &(&(&end[k] - &base[k]) * &weight)
            } else {
                base[k].clone()
            }
        });
        let combine = |normal: &[Rational; 3]| {
            let along = (0..3).fold(MPoly::zero(nvars), |sum, k| {
                &sum + &homogeneous[k].scale(&normal[k])
            });
            let offset = rational_dot(normal, &origin);
            &along - &homogeneous[3].scale(&offset)
        };
        let lambda = MPoly::var(nvars, 0);
        let length = rational_dot(&direction, &direction);
        let square = |variable: usize, value: &Rational| {
            let x = MPoly::var(nvars, variable);
            &(&x * &x) - &MPoly::constant(nvars, value.clone())
        };
        let mut equations = vec![
            combine(&first),
            combine(&second),
            &(&lambda * &homogeneous[3]).scale(&length) - &combine(&direction),
            square(3, &tangent_square),
        ];
        let tangent = dec.enclose_budgeted_at(&t.t, query.bits)?;
        let slack = Rational::one().mul_pow2(-i64::from(query.bits / 2));
        let mut lo = vec![Rational::zero(); nvars];
        let mut hi = vec![Rational::zero(); nvars];
        lo[3] = tangent.lo() - slack.clone();
        hi[3] = tangent.hi() + slack.clone();
        if ends {
            let t_mouth = rational(&t.mouth);
            let Some(mouth) = t_mouth else {
                return Ok(None);
            };
            let value = (Rational::one() + tangent_square.clone() * mouth.clone() * mouth)
                / (q(2) * (Rational::one() + tangent_square.clone()));
            equations.push(square(4, &value));
            let weight = pair
                .iter()
                .find(|node| node.s.abs() == Rational::one())
                .expect("an end span has an end node");
            let weight = dec.enclose_budgeted_at(&weight.values[WEIGHT], query.bits)?;
            lo[4] = weight.lo() - slack.clone();
            hi[4] = weight.hi() + slack.clone();
        }
        let width = &data.nodes[index + 1].s - &data.nodes[index].s;
        lo[1] = (&rect[0][0] - &data.nodes[index].s) / width.clone();
        hi[1] = (&rect[0][1] - &data.nodes[index].s) / width;
        lo[2] = rect[1][0].clone();
        hi[2] = rect[1][1].clone();
        // Along the line, the rectangle's whole image bounds the parameter.
        let hull = query.spans[index].clip(rect, query.bits).hull()?;
        let relative: [Ball; 3] = std::array::from_fn(|i| &hull[i] - &query.enclosed[0][i]);
        let along = dot(&relative, &query.enclosed[1]);
        let scale = dot(&query.enclosed[1], &query.enclosed[1]);
        let Some(parameter) = along.checked_div(&scale) else {
            return Ok(None);
        };
        lo[0] = parameter.lo() - slack.clone();
        hi[0] = parameter.hi() + slack;
        let Some(system) = PolySystem::new(equations) else {
            return Ok(None);
        };
        Ok(RootBox::new(lo, hi).map(|candidate| (Arc::new(system), candidate)))
    }
}

// One span's homogeneous local coordinates as polynomials in the span
// parameter (variable 1), v (2) and t (3), or `None` off rational node data.
fn span_terms(def: &RollingEnvelopeDefinition, pair: &[Node], nvars: usize) -> Option<[MPoly; 4]> {
    let rows = span_rows(def, pair);
    let bernstein = |variable: usize, degree: usize, i: usize| -> MPoly {
        let x = MPoly::var(nvars, variable);
        let rest = &MPoly::constant(nvars, Rational::one()) - &x;
        let mut value = MPoly::constant(nvars, binomial(degree, i));
        for _ in 0..i {
            value = &value * &x;
        }
        for _ in i..degree {
            value = &value * &rest;
        }
        value
    };
    let along: Vec<MPoly> = (0..=DEGREE).map(|i| bernstein(1, DEGREE, i)).collect();
    let across: Vec<MPoly> = (0..3).map(|j| bernstein(2, 2, j)).collect();
    let tangent = MPoly::var(nvars, 3);
    let mut components = Vec::with_capacity(4);
    for k in 0..4 {
        let mut sum = MPoly::zero(nvars);
        for (j, row) in rows.iter().enumerate() {
            let basis = &across[j];
            let mut line = MPoly::zero(nvars);
            for (i, along) in along.iter().enumerate() {
                let even = row[k].even[i].rational_ref()?;
                let odd = row[k].odd[i].rational_ref()?;
                let coefficient = &MPoly::constant(nvars, even.clone()) + &tangent.scale(odd);
                line = &line + &(along * &coefficient);
            }
            sum = &sum + &(&line * basis);
        }
        components.push(sum);
    }
    components.try_into().ok()
}
