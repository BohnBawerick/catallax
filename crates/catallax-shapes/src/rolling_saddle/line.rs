//! What the line solver in [`crate::rational_sheet`] needs from a saddle.
//!
//! The cylinders may lean at an angle with no rational unit axis, such as 45
//! degrees, so the carrier's world controls are rational only over three square
//! roots: the lengths of the two axes and of their cross product. Each root
//! enters a span's equations as a variable fixed by its square, which keeps the
//! system rational and leaves each hit a leaf of it.

use super::*;
use crate::rational_sheet::{LineSheet, Query, Rect, line_hits};
use crate::rolling_envelope::binomial;
use catallax_numbers::{MPoly, PolySystem, RealPoly, RootBox};

fn unresolved_line(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

impl RollingSaddle {
    /// Every point where a whole line meets the closed carrier, in line order.
    /// Each is an exact transverse crossing and appears once. A tangency, or a
    /// contact this cannot separate, refuses. It never reports a hit or a miss
    /// from an enclosure alone.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        line_hits(self, origin, direction, dec)
    }
}

impl LineSheet for RollingSaddle {
    fn query_bits(&self, dec: &Decider) -> Result<u32, Refusal> {
        RollingSaddle::query_bits(self, dec)
    }
    fn budgeted_frame(&self, bits: u32, dec: &mut Decider) -> Result<LocalFrame, Refusal> {
        self.data.definition.budgeted_frame(bits, dec)
    }
    fn local_line(&self, origin: &Point3, direction: &Vector3) -> [[Real; 3]; 2] {
        let def = &self.data.definition;
        [
            def.local(origin),
            [
                direction.dot(&def.across),
                direction.dot(&def.along),
                direction.dot(&def.axis),
            ],
        ]
    }
    fn spans(&self, bits: u32) -> Vec<Patch> {
        self.patches(bits)
    }
    fn chart(&self) -> Rect {
        [
            [Rational::from_i64(-1), Rational::from_i64(3)],
            [Rational::zero(), Rational::one()],
        ]
    }
    fn knots(&self) -> Vec<Rational> {
        self.data.nodes.iter().map(|node| node.u.clone()).collect()
    }
    fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        RollingSaddle::on_boundary(self, point, dec)
    }
    // A row lies on its cylinder, so its contacts are the line's roots there.
    fn row_contacts(
        &self,
        high: bool,
        line: &[[Real; 3]; 2],
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        let def = &self.data.definition;
        let [o, d] = line;
        let (radius, start, axis) = if high {
            (
                &def.narrow,
                [&o[0] - &def.gap, o[1].clone(), o[2].clone()],
                [Real::zero(), def.sine.clone(), def.cosine.clone()],
            )
        } else {
            (
                &def.wide,
                o.clone(),
                [Real::zero(), Real::zero(), Real::one()],
            )
        };
        let dot = |a: &[Real; 3], c: &[Real; 3]| {
            &(&(&a[0] * &c[0]) + &(&a[1] * &c[1])) + &(&a[2] * &c[2])
        };
        let flat = |a: &[Real; 3]| -> [Real; 3] {
            let along = dot(a, &axis);
            std::array::from_fn(|i| &a[i] - &(&along * &axis[i]))
        };
        let (start, along) = (flat(&start), flat(d));
        let quadratic = RealPoly::from_terms(
            1,
            vec![
                (vec![0], &dot(&start, &start) - &radius.square()),
                (vec![1], &Real::from_i64(2) * &dot(&start, &along)),
                (vec![2], dot(&along, &along)),
            ],
        );
        if quadratic.is_zero(dec)? {
            return Err(unresolved_line(
                "a line on a rail's cylinder needs a tangency solver",
            ));
        }
        quadratic.real_roots(dec)
    }
    fn knot_section(&self, index: usize) -> [[Real; 4]; 3] {
        let data = &self.data;
        let def = &data.definition;
        let (pair, column) = if index + 1 < data.nodes.len() {
            (&data.nodes[index..index + 2], 0)
        } else {
            (&data.nodes[index - 1..index + 1], DEGREE)
        };
        let [wide, middle, narrow] = span_rows(def, pair);
        // The narrow frame's axes in the wide one, about its foot.
        let local = |row: &[Vec<Real>; 4]| -> [Real; 4] {
            let [x, e, z, w] = [0, 1, 2, 3].map(|k| row[k][column].clone());
            [
                &x + &(&def.gap * &w),
                &(&def.cosine * &e) + &(&def.sine * &z),
                &(&def.cosine * &z) - &(&def.sine * &e),
                w,
            ]
        };
        [
            [0, 1, 2, 3].map(|k| wide[k][column].clone()),
            local(&middle),
            local(&narrow),
        ]
    }
    fn span_system(
        &self,
        query: &Query<'_, Self>,
        index: usize,
        rect: &Rect,
        dec: &mut Decider,
    ) -> Result<Option<(Arc<PolySystem>, RootBox)>, Refusal> {
        let data = &self.data;
        let def = &data.definition;
        let rational = |value: &Real| value.as_rational();
        let triple = |x: &Real, y: &Real, z: &Real| -> Option<[Rational; 3]> {
            Some([rational(x)?, rational(y)?, rational(z)?])
        };
        let point = |p: &Point3| triple(&p.x, &p.y, &p.z);
        let vector = |v: &Vector3| triple(&v.x, &v.y, &v.z);
        let (Some(origin), Some(direction)) = (point(query.origin), vector(query.direction)) else {
            return Ok(None);
        };
        let (Some(low), Some(foot), Some(a), Some(c)) = (
            point(&def.origin),
            point(&def.foot),
            vector(&def.axes[0]),
            vector(&def.axes[1]),
        ) else {
            return Ok(None);
        };
        let dot = |x: &[Rational; 3], y: &[Rational; 3]| -> Rational {
            (0..3).fold(Rational::zero(), |sum, i| sum + &x[i] * &y[i])
        };
        let cross = |x: &[Rational; 3], y: &[Rational; 3]| -> [Rational; 3] {
            std::array::from_fn(|i| {
                let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                &x[j] * &y[k] - &x[k] * &y[j]
            })
        };
        let meet = cross(&c, &a);
        let inner = dot(&a, &c);
        let [square_a, square_c, square_meet] = [dot(&a, &a), dot(&c, &c), dot(&meet, &meet)];
        // A root that is rational stays a constant, and the rest get variables.
        let squares = [&square_a, &square_c, &square_meet];
        let mut nvars = 3;
        let mut slots = [None; 3];
        let roots: [Option<Rational>; 3] =
            squares.map(|square| math::exact_sqrt(&Real::from(square.clone())).as_rational());
        for (slot, root) in slots.iter_mut().zip(&roots) {
            if root.is_none() {
                *slot = Some(nvars);
                nvars += 1;
            }
        }
        let generator = |which: usize| match (&roots[which], slots[which]) {
            (Some(root), _) => MPoly::constant(nvars, root.clone()),
            (None, Some(slot)) => MPoly::var(nvars, slot),
            _ => unreachable!("every root is rational or has a slot"),
        };
        let [ga, gc, gm] = [0, 1, 2].map(generator);
        let (ga_gm, gc_gm) = (&ga * &gm, &gc * &gm);
        let scaled = |factor: &MPoly, values: [Rational; 3]| -> [MPoly; 3] {
            values.map(|value| factor.scale(&value))
        };
        // A reflection moves the frame, but the moved axes' cross product turns
        // the other way.
        let turn = Rational::from_i64(if self.data.reversing { -1 } else { 1 });
        let across = scaled(
            &gm,
            meet.clone().map(|v| v * turn.clone() / square_meet.clone()),
        );
        let along = scaled(
            &ga_gm,
            std::array::from_fn(|k| {
                (&c[k] - &(&inner * &a[k] / square_a.clone())) / square_meet.clone()
            }),
        );
        let axis = scaled(&ga, a.clone().map(|v| v / square_a.clone()));
        let side = scaled(
            &gc_gm,
            std::array::from_fn(|k| {
                (&inner * &c[k] - &square_c * &a[k]) / (&square_c * &square_meet)
            }),
        );
        let narrow = scaled(&gc, c.clone().map(|v| v / square_c.clone()));
        let pair = &data.nodes[index..index + 2];
        let rows = span_rows(def, pair);
        let t = MPoly::var(nvars, 1);
        let v = MPoly::var(nvars, 2);
        let one = MPoly::constant(nvars, Rational::one());
        let bernstein = |x: &MPoly, degree: usize, i: usize| -> MPoly {
            let rest = &one - x;
            let mut value = MPoly::constant(nvars, binomial(degree, i));
            for _ in 0..i {
                value = &value * x;
            }
            for _ in i..degree {
                value = &value * &rest;
            }
            value
        };
        let along_u: Vec<MPoly> = (0..=DEGREE).map(|i| bernstein(&t, DEGREE, i)).collect();
        let mut homogeneous: [MPoly; 4] = std::array::from_fn(|_| MPoly::zero(nvars));
        for (j, row) in rows.iter().enumerate() {
            let mut lines: [MPoly; 4] = std::array::from_fn(|_| MPoly::zero(nvars));
            for (line, values) in lines.iter_mut().zip(row) {
                for (basis, value) in along_u.iter().zip(values) {
                    let Some(value) = rational(value) else {
                        return Ok(None);
                    };
                    *line = &*line + &basis.scale(&value);
                }
            }
            let across_v = bernstein(&v, 2, j);
            let (base, second, third) = if j == 0 {
                (&low, &along, &axis)
            } else {
                (&foot, &side, &narrow)
            };
            for k in 0..3 {
                let world = &(&(&lines[3].scale(&base[k]) + &(&lines[0] * &across[k]))
                    + &(&lines[1] * &second[k]))
                    + &(&lines[2] * &third[k]);
                homogeneous[k] = &homogeneous[k] + &(&world * &across_v);
            }
            homogeneous[3] = &homogeneous[3] + &(&lines[3] * &across_v);
        }
        let least = (0..3)
            .min_by_key(|&i| direction[i].abs())
            .expect("three axes");
        let unit: [Rational; 3] = std::array::from_fn(|i| {
            if i == least {
                Rational::one()
            } else {
                Rational::zero()
            }
        });
        let first = cross(&direction, &unit);
        let second = cross(&direction, &first);
        let combine = |normal: &[Rational; 3]| {
            let along = (0..3).fold(MPoly::zero(nvars), |sum, k| {
                &sum + &homogeneous[k].scale(&normal[k])
            });
            &along - &homogeneous[3].scale(&dot(normal, &origin))
        };
        let lambda = MPoly::var(nvars, 0);
        let mut equations = vec![
            combine(&first),
            combine(&second),
            &(&lambda * &homogeneous[3]).scale(&dot(&direction, &direction)) - &combine(&direction),
        ];
        let slack = Rational::one().mul_pow2(-i64::from(query.bits / 2));
        let mut lo = vec![Rational::zero(); nvars];
        let mut hi = vec![Rational::zero(); nvars];
        for (which, slot) in slots.iter().enumerate() {
            let Some(slot) = *slot else {
                continue;
            };
            let x = MPoly::var(nvars, slot);
            equations.push(&(&x * &x) - &MPoly::constant(nvars, squares[which].clone()));
            let root =
                dec.enclose_budgeted_at(&Real::from(squares[which].clone()).sqrt(), query.bits)?;
            lo[slot] = root.lo() - slack.clone();
            hi[slot] = root.hi() + slack.clone();
        }
        let width = &pair[1].u - &pair[0].u;
        lo[1] = (&rect[0][0] - &pair[0].u) / width.clone();
        hi[1] = (&rect[0][1] - &pair[0].u) / width;
        lo[2] = rect[1][0].clone();
        hi[2] = rect[1][1].clone();
        // Along the line, the rectangle's whole image bounds the parameter.
        let hull = query.spans[index].clip(rect, query.bits).hull()?;
        let relative: [Ball; 3] = std::array::from_fn(|i| &hull[i] - &query.enclosed[0][i]);
        let inner_ball = |x: &[Ball; 3], y: &[Ball; 3]| {
            &(&(&x[0] * &y[0]) + &(&x[1] * &y[1])) + &(&x[2] * &y[2])
        };
        let along_line = inner_ball(&relative, &query.enclosed[1]);
        let scale = inner_ball(&query.enclosed[1], &query.enclosed[1]);
        let Some(parameter) = along_line.checked_div(&scale) else {
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
