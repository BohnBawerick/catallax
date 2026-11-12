//! The fillet face as a surface of a body: rigid transport, chart preimages,
//! exact rail tests, normals, sides and the complete line solver.
//!
//! Every copy keeps what the fit proved in the fit's frame, and a rigid frame
//! that places it. Queries map into the fit's frame, where the carrier's
//! enclosed controls and exact rows live. The carrier is a sheet over `tau`
//! and a quadratic second parameter, so the line solver of
//! [`crate::rational_sheet`] applies to it: its span equations come from the
//! exact rows, and each rail's contacts with a line come from its support.

use super::*;
use crate::rational_sheet::{LineSheet, LocalFrame, Patch, Query, Rect, cover, line_hits};
use crate::rolling_envelope::binomial;
use crate::{BSplineCurve, Pose, SplineBoundary, Vector3, math};
use catallax_numbers::{MPoly, PolySystem, RootBox};

/// Precision of the stored control enclosures.
const NET_BITS: u32 = 128;

fn unresolved(detail: &str) -> Refusal {
    Refusal::new("numbers.geometric-gap-unproven", detail)
}

/// The carrier's homogeneous controls, enclosed.
pub(super) fn net(carrier: &BSplineSurface) -> Vec<[Ball; 4]> {
    carrier
        .controls()
        .iter()
        .zip(carrier.weights())
        .map(|(control, weight)| {
            let w = weight.as_rational().expect("rational carrier weights");
            let coordinate = |value: &Real| {
                let value = value.as_rational().expect("rational carrier controls");
                Ball::from_rational(&(value * w.clone()), NET_BITS)
            };
            [
                coordinate(&control.x),
                coordinate(&control.y),
                coordinate(&control.z),
                Ball::from_rational(&w, NET_BITS),
            ]
        })
        .collect()
}

/// The plane `normal . x = offset` in the fit's frame, when the frame is
/// rational.
pub(super) fn local_plane(
    frame: &Isometry,
    normal: &[Rational; 3],
    offset: &Rational,
) -> Option<([Rational; 3], Rational)> {
    let matrix = frame.rational_matrix()?;
    let translation = [
        frame.translation.x.as_rational()?,
        frame.translation.y.as_rational()?,
        frame.translation.z.as_rational()?,
    ];
    // x = M y + t, so n . x = c becomes (M^T n) . y = c - n . t.
    let local: [Rational; 3] = std::array::from_fn(|j| {
        (0..3).fold(Rational::zero(), |sum, i| sum + matrix[i][j] * &normal[i])
    });
    let moved = (0..3).fold(Rational::zero(), |sum, i| {
        sum + &normal[i] * &translation[i]
    });
    Some((local, offset - &moved))
}

impl SplineFillet {
    /// Exact rigid transport of the face, its supports and its proofs.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }

    pub(crate) fn mapped(&self, map: &Isometry) -> Self {
        Self {
            data: Arc::new(Data {
                fit: self.data.fit.clone(),
                carrier: self.data.carrier.mapped(map),
                supports: self.data.supports.each_ref().map(|patch| patch.mapped(map)),
                frame: self.data.frame.then(map),
            }),
        }
    }

    /// Positive identity of the fitted law and carrier, for cache reuse.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data.fit, &other.data.fit)
    }

    fn query_bits(&self, dec: &Decider) -> Result<u32, Refusal> {
        match NET_BITS.min(dec.policy().max_precision_bits) {
            0 => Err(Refusal::new(
                "numbers.precision-limit",
                "no fillet query precision available",
            )),
            bits => Ok(bits),
        }
    }

    fn patches(&self, bits: u32) -> Vec<Patch> {
        let fit = &self.data.fit;
        Patch::at(
            &fit.carrier,
            &fit.net,
            NET_BITS,
            fit.carrier.degrees()[0],
            bits,
        )
    }

    /// The frame's origin and the fit's axes in it: `x = M y + t` has
    /// origin `t` and the columns of `M` as axes.
    fn axes(&self) -> (Point3, [Vector3; 3]) {
        let frame = &self.data.frame;
        let matrix = frame.matrix();
        let column = |j: usize| {
            Vector3::new(
                matrix[0][j].clone(),
                matrix[1][j].clone(),
                matrix[2][j].clone(),
            )
        };
        let origin = Point3::new(
            frame.translation.x.clone(),
            frame.translation.y.clone(),
            frame.translation.z.clone(),
        );
        (origin, [column(0), column(1), column(2)])
    }

    fn local_frame(&self, bits: u32, dec: &mut Decider) -> Result<LocalFrame, Refusal> {
        let (origin, [a, b, c]) = self.axes();
        LocalFrame::budgeted(&origin, [&a, &b, &c], bits, dec)
    }

    /// Enclosed world controls and weights, for a writer that rounds them.
    pub fn enclosed_controls(&self, dec: &mut Decider) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        let bits = self.query_bits(dec)?;
        let frame = self.local_frame(bits, dec)?;
        self.data
            .fit
            .net
            .iter()
            .map(|h| {
                let local = [0, 1, 2].map(|k| {
                    h[k].with_precision(bits)
                        .checked_div(&h[3].with_precision(bits))
                        .ok_or_else(|| unresolved("a fillet weight is not positive"))
                });
                let [x, y, z] = local;
                Ok((frame.world(&[x?, y?, z?]), h[3].clone()))
            })
            .collect()
    }

    /// Encloses every chart preimage of a spatial box.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        dec.check_deadline()?;
        let bits = self.query_bits(dec)?;
        let ceiling = dec.policy().max_precision_bits;
        let point = point
            .each_ref()
            .map(|value| value.with_precision(value.precision().min(ceiling)));
        let (origin, [a, b, c]) = self.axes();
        let frame = LocalFrame::new(&origin, [&a, &b, &c], bits, dec)?;
        let found = cover(self.patches(bits), &frame.local(&point), dec)?;
        Ok(found.map(|cell| {
            cell.map(|[lo, hi]| {
                Ball::from_rational(&lo, bits).hull(&Ball::from_rational(&hi, bits))
            })
        }))
    }

    /// The rails, the carrier's own boundary rows, so a neighbour can share
    /// them exactly.
    pub fn rails(&self) -> [BSplineCurve; 2] {
        [SplineBoundary::VMin, SplineBoundary::VMax]
            .map(|side| BSplineCurve::boundary(&self.data.carrier, side))
    }

    fn rail_of(&self, point: &Point3, dec: &mut Decider) -> Result<Option<usize>, Refusal> {
        for (side, rail) in self.rails().iter().enumerate() {
            if rail.parameter_of(point, dec)?.is_some() {
                return Ok(Some(side));
            }
        }
        Ok(None)
    }

    /// Whether a point is exactly on a rail. `Ok(false)` says only that.
    pub fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        math::validate_point(point, dec)?;
        Ok(self.rail_of(point, dec)?.is_some())
    }

    /// Exact normal at a rail point: its support's normal, pointing to the
    /// ball's centre. Interior points have no exact chart parameters and
    /// refuse.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        let Some(side) = self.rail_of(point, dec)? else {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "a fillet normal needs exact chart parameters, and only its rails have them",
            ));
        };
        let normal = self.data.supports[side].normal_at(point, dec)?;
        // The centre lies at P + side r N/|N| in the fit's frame; a reversing
        // frame turns the support's own normal round.
        let toward = self.data.fit.definition.sides[side] == Sign::Positive;
        Ok(if toward != self.data.frame.reversing() {
            normal
        } else {
            normal.scaled(&Real::from_i64(-1))
        })
    }

    /// Whether a point is on the carrier. A nonzero sign proves it off, and
    /// names no global side of the finite support. Exact zero is decided on
    /// the rails only, so an interior point within the enclosure's width
    /// refuses.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        if self.on_boundary(point, dec)? {
            return Ok(Sign::Zero);
        }
        let bits = self.query_bits(dec)?;
        let at = math::budgeted_point(point, bits, dec)?;
        if self.parameters_of_box(&at, dec)?.is_none() {
            return Ok(Sign::Positive);
        }
        Err(unresolved(
            "a point is within the carrier's enclosure width",
        ))
    }

    /// Every point where a whole line meets the closed carrier, in line
    /// order. Each is an exact transverse crossing and appears once. A
    /// tangency, or a contact this cannot separate, refuses.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        line_hits(self, origin, direction, dec)
    }
}

impl LineSheet for SplineFillet {
    fn query_bits(&self, dec: &Decider) -> Result<u32, Refusal> {
        SplineFillet::query_bits(self, dec)
    }
    fn budgeted_frame(&self, bits: u32, dec: &mut Decider) -> Result<LocalFrame, Refusal> {
        self.local_frame(bits, dec)
    }
    fn local_line(&self, origin: &Point3, direction: &Vector3) -> [[Real; 3]; 2] {
        let back = self.data.frame.inverse();
        let o = back.apply_point(origin);
        let d = back.apply_vector(direction);
        [[o.x, o.y, o.z], [d.x, d.y, d.z]]
    }
    fn spans(&self, bits: u32) -> Vec<Patch> {
        self.patches(bits)
    }
    fn chart(&self) -> Rect {
        let nodes = &self.data.fit.nodes;
        [
            [nodes[0].tau.clone(), nodes[nodes.len() - 1].tau.clone()],
            [Rational::zero(), Rational::one()],
        ]
    }
    fn knots(&self) -> Vec<Rational> {
        self.data
            .fit
            .nodes
            .iter()
            .map(|node| node.tau.clone())
            .collect()
    }
    fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        SplineFillet::on_boundary(self, point, dec)
    }
    // Each rail lies on its support, so its contacts with a line are among
    // the line's hits on that support, in the fit's frame.
    fn row_contacts(
        &self,
        high: bool,
        line: &[[Real; 3]; 2],
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        let support = &self.data.fit.definition.supports[usize::from(high)];
        let [o, d] = line;
        let origin = Point3::new(o[0].clone(), o[1].clone(), o[2].clone());
        let direction = Vector3::new(d[0].clone(), d[1].clone(), d[2].clone());
        let scale = direction.dot(&direction);
        // The direction is nonzero, so the division is defined.
        Ok(support
            .line_hits(&origin, &direction, dec)?
            .iter()
            .map(|hit| &hit.minus(&origin).dot(&direction) / &scale)
            .collect())
    }
    fn knot_section(&self, index: usize) -> [[Real; 4]; 3] {
        let rows = &self.data.fit.rows;
        let (span, column) = if index < rows.len() {
            (&rows[index], 0)
        } else {
            let last = &rows[rows.len() - 1];
            (last, last[0][0].len() - 1)
        };
        span.each_ref()
            .map(|row| [0, 1, 2, 3].map(|k| Real::from_rational(row[k][column].clone())))
    }
    fn span_system(
        &self,
        query: &Query<'_, Self>,
        index: usize,
        rect: &Rect,
        _: &mut Decider,
    ) -> Result<Option<(Arc<PolySystem>, RootBox)>, Refusal> {
        let rational = |value: &Real| value.as_rational();
        let triple = |v: &[Real; 3]| -> Option<[Rational; 3]> {
            Some([rational(&v[0])?, rational(&v[1])?, rational(&v[2])?])
        };
        let (Some(origin), Some(direction)) = (triple(&query.exact[0]), triple(&query.exact[1]))
        else {
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
        let rows = &self.data.fit.rows[index];
        let degree = rows[0][0].len() - 1;
        let along_u: Vec<MPoly> = (0..=degree).map(|i| bernstein(&t, degree, i)).collect();
        let mut homogeneous: [MPoly; 4] = std::array::from_fn(|_| MPoly::zero(nvars));
        for (j, row) in rows.iter().enumerate() {
            let across_v = bernstein(&v, 2, j);
            for (slot, values) in homogeneous.iter_mut().zip(row) {
                let line = along_u
                    .iter()
                    .zip(values)
                    .fold(MPoly::zero(nvars), |sum, (basis, value)| {
                        &sum + &basis.scale(value)
                    });
                *slot = &*slot + &(&line * &across_v);
            }
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
        let nodes = &self.data.fit.nodes;
        let start = &nodes[index].tau;
        let width = &nodes[index + 1].tau - start;
        lo[1] = (&rect[0][0] - start) / width.clone();
        hi[1] = (&rect[0][1] - start) / width;
        lo[2] = rect[1][0].clone();
        hi[2] = rect[1][1].clone();
        // Along the line, the rectangle's whole image bounds the parameter.
        let hull = query.spans[index].clip(rect, query.bits).hull()?;
        let relative: [Ball; 3] = std::array::from_fn(|i| &hull[i] - &query.enclosed[0][i]);
        let inner = |x: &[Ball; 3], y: &[Ball; 3]| {
            &(&(&x[0] * &y[0]) + &(&x[1] * &y[1])) + &(&x[2] * &y[2])
        };
        let along_line = inner(&relative, &query.enclosed[1]);
        let scale = inner(&query.enclosed[1], &query.enclosed[1]);
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
