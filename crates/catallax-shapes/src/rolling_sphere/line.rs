//! What the line solver in [`crate::rational_sheet`] needs from a sphere band.
//!
//! Every row is placed in one frame about the sphere's centre, so a span's
//! world controls are rational when that frame and the rows are. A band with
//! an irrational frame or row has no span equations and its crossings refuse.

use super::*;
use crate::rational_sheet::{LineSheet, Query, Rect, line_hits};
use crate::rolling_envelope::binomial;
use catallax_numbers::{MPoly, PolySystem, RealPoly, RootBox};

fn unresolved_line(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

impl RollingSphere {
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

impl LineSheet for RollingSphere {
    fn query_bits(&self, dec: &Decider) -> Result<u32, Refusal> {
        RollingSphere::query_bits(self, dec)
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
        RollingSphere::knots(self)
    }
    fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        RollingSphere::on_boundary(self, point, dec)
    }
    // The low row lies on the sphere and the high row on the cylinder, so
    // their contacts are the line's roots there.
    fn row_contacts(
        &self,
        high: bool,
        line: &[[Real; 3]; 2],
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        let def = &self.data.definition;
        let [o, d] = line;
        let (radius, start, along) = if high {
            (
                &def.cylinder,
                [&o[0] - &def.gap, o[1].clone(), Real::zero()],
                [d[0].clone(), d[1].clone(), Real::zero()],
            )
        } else {
            (&def.sphere, o.clone(), d.clone())
        };
        let dot = |a: &[Real; 3], c: &[Real; 3]| {
            &(&(&a[0] * &c[0]) + &(&a[1] * &c[1])) + &(&a[2] * &c[2])
        };
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
                "a line on a rail's support needs a tangency solver",
            ));
        }
        quadratic.real_roots(dec)
    }
    fn knot_section(&self, index: usize) -> [[Real; 4]; 3] {
        let nodes = &self.data.nodes;
        let (pair, column) = if index + 1 < nodes.len() {
            (&nodes[index..index + 2], 0)
        } else {
            (&nodes[index - 1..index + 1], DEGREE)
        };
        let rows = span_rows(&self.data.definition, pair);
        rows.map(|row| [0, 1, 2, 3].map(|k| row[k][column].clone()))
    }
    fn span_system(
        &self,
        query: &Query<'_, Self>,
        index: usize,
        rect: &Rect,
        _: &mut Decider,
    ) -> Result<Option<(Arc<PolySystem>, RootBox)>, Refusal> {
        let def = &self.data.definition;
        let rational = |value: &Real| value.as_rational();
        let triple = |x: &Real, y: &Real, z: &Real| -> Option<[Rational; 3]> {
            Some([rational(x)?, rational(y)?, rational(z)?])
        };
        let point = |p: &Point3| triple(&p.x, &p.y, &p.z);
        let vector = |v: &Vector3| triple(&v.x, &v.y, &v.z);
        let (Some(origin), Some(direction)) = (point(query.origin), vector(query.direction)) else {
            return Ok(None);
        };
        let (Some(centre), Some(across), Some(along), Some(axis)) = (
            point(&def.centre),
            vector(&def.across),
            vector(&def.along),
            vector(&def.axis),
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
        let nvars = 3;
        let one = MPoly::constant(nvars, Rational::one());
        let t = MPoly::var(nvars, 1);
        let v = MPoly::var(nvars, 2);
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
        let rows = span_rows(def, &self.data.nodes[index..index + 2]);
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
            for k in 0..3 {
                let world = &(&(&lines[3].scale(&centre[k]) + &lines[0].scale(&across[k]))
                    + &lines[1].scale(&along[k]))
                    + &lines[2].scale(&axis[k]);
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
        let equations = vec![
            combine(&first),
            combine(&second),
            &(&lambda * &homogeneous[3]).scale(&dot(&direction, &direction)) - &combine(&direction),
        ];
        let slack = Rational::one().mul_pow2(-i64::from(query.bits / 2));
        let mut lo = vec![Rational::zero(); nvars];
        let mut hi = vec![Rational::zero(); nvars];
        let knots = RollingSphere::knots(self);
        let start = &knots[index];
        let width = &knots[index + 1] - start;
        lo[1] = (&rect[0][0] - start) / width.clone();
        hi[1] = (&rect[0][1] - start) / width;
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
