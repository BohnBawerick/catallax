//! Complete transverse crossings of a line with a checked offset sweep.
//!
//! The carrier sweeps a closed convex composite cubic C straight along its axis
//! a. A line with direction d off that axis stays in the plane through it with
//! normal n = d x a, so it can only meet the carrier where g = n.(C - o) is
//! zero, and g'(u) = 0 is exactly where the line is tangent there. On a span
//! g is a cubic whose Bernstein coefficients come from the actual controls.
//!
//! Signs at the span ends are exact, so a crossing at a join is seen once. Each
//! span is then cut until every cell keeps g from zero or keeps g' from zero; a
//! monotone cell owns at most one zero, decided by its end signs. A cell that
//! never settles, a tangency included, refuses. A crossing clear of the joins is
//! a certified root of the span's own cubic and the line parameter it gives, and
//! every crossing still has to land inside the sweep's axial interval. Nothing
//! is accepted from an enclosure.

use super::*;
use catallax_numbers::RealPoly;

// Proposal bounds, not tolerances: running out refuses and proves nothing.
const MAX_DEPTH: u32 = 48;
const MAX_CELLS: usize = 4096;
const MAX_NARROWING: usize = 16;

fn unresolved(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

fn strict(value: &Ball) -> Option<Sign> {
    if value.lo() > Rational::zero() {
        Some(Sign::Positive)
    } else if value.hi() < Rational::zero() {
        Some(Sign::Negative)
    } else {
        None
    }
}

// One strict sign across a whole Bernstein row.
fn settled(row: &[Ball]) -> bool {
    let first = strict(&row[0]);
    first.is_some() && row[1..].iter().all(|value| strict(value) == first)
}

fn opposite(a: Sign, c: Sign) -> bool {
    matches!(
        (a, c),
        (Sign::Positive, Sign::Negative) | (Sign::Negative, Sign::Positive)
    )
}

// A part of one span, in the span's own parameter on [0,1].
struct Cell {
    from: Rational,
    to: Rational,
    net: [Ball; 4],
    ends: [Sign; 2],
    depth: u32,
}
impl Cell {
    fn monotone(&self) -> bool {
        let slopes: [Ball; 3] = std::array::from_fn(|i| &self.net[i + 1] - &self.net[i]);
        settled(&slopes)
    }
}

// Where a crossing was found, in its span's own parameter.
enum Root {
    Exact(Rational),
    Between(Rational, Rational),
}

struct Query<'a> {
    bits: u32,
    half: Ball,
    origin: &'a Point3,
    direction: &'a Vector3,
    // g on each span's controls, exactly and enclosed.
    exact: Vec<[Real; 4]>,
    nets: Vec<[Ball; 4]>,
    // The controls along the line's own shadow in the section plane.
    along: Vec<[Real; 4]>,
    shadow: Real,
    enclosed_shadow: Ball,
}

impl SplineSurface {
    /// Every point where a whole line meets the finite carrier, in line order.
    /// Each is one exact transverse crossing, counted once at a join or seam.
    /// A tangency, an overlap or a contact this cannot separate refuses.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        dec.check_deadline()?;
        if let Some(transport) = &self.data.transport {
            let hits = transport.source.line_hits(
                &transport.inverse.apply_point(origin),
                &transport.inverse.apply_vector(direction),
                dec,
            )?;
            return Ok(hits
                .iter()
                .map(|hit| transport.pose.apply_point(hit))
                .collect());
        }
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        if math::sign(&direction.norm_squared(), dec)? == Sign::Zero {
            return Err(Refusal::new("shapes.zero-vector", "direction is zero"));
        }
        let normal = direction.cross(self.axis());
        if math::sign(&normal.norm_squared(), dec)? == Sign::Zero {
            // Parallel to the sweep, the line lies wholly on one side of the section.
            if self.cross_section_side(origin, dec)? == Sign::Zero {
                return Err(unresolved(
                    "a line along the sweep overlaps a finite generator",
                ));
            }
            return Ok(Vec::new());
        }
        let query = Query::new(self, origin, direction, &normal, dec)?;
        let height = self.axial_coordinate(origin);
        let lean = direction.dot(self.axis());
        let level = math::sign(&lean, dec)? == Sign::Zero;
        let mut parameters = Vec::with_capacity(2);
        for (span, root) in query.roots(dec)? {
            let parameter = query.parameter(span, &root, dec)?;
            // A level line keeps one exact height, with no root in the comparison.
            let axial = if level {
                height.clone()
            } else {
                &height + &(&parameter * &lean)
            };
            if math::sign(&(&axial - &self.data.axial[0]), dec)? == Sign::Negative
                || math::sign(&(&axial - &self.data.axial[1]), dec)? == Sign::Positive
            {
                continue;
            }
            parameters.push(parameter);
        }
        if let [first, second] = parameters.as_mut_slice()
            && math::sign(&(&*first - &*second), dec)? == Sign::Positive
        {
            std::mem::swap(first, second);
        }
        Ok(parameters.iter().map(|at| query.world(at)).collect())
    }
}

impl<'a> Query<'a> {
    fn new(
        surface: &SplineSurface,
        origin: &'a Point3,
        direction: &'a Vector3,
        normal: &Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(precision_refusal());
        }
        let axis = surface.axis();
        let across = direction.minus(&axis.scaled(&direction.dot(axis)));
        let shadow = across.norm_squared();
        let enclosed_shadow = dec.enclose_budgeted_at(&shadow, bits)?;
        let curve = surface.lower_edge();
        let mut exact = Vec::with_capacity(curve.control_count() / 3);
        let mut nets = Vec::with_capacity(curve.control_count() / 3);
        let mut along = Vec::with_capacity(curve.control_count() / 3);
        for (span, _) in curve.spans() {
            dec.check_deadline()?;
            let relative: [Vector3; 4] =
                std::array::from_fn(|i| curve.control(span - 3 + i).minus(origin));
            let values = relative.each_ref().map(|control| control.dot(normal));
            let net = [
                dec.enclose_budgeted_at(&values[0], bits)?,
                dec.enclose_budgeted_at(&values[1], bits)?,
                dec.enclose_budgeted_at(&values[2], bits)?,
                dec.enclose_budgeted_at(&values[3], bits)?,
            ];
            exact.push(values);
            nets.push(net);
            along.push(relative.each_ref().map(|control| control.dot(&across)));
        }
        Ok(Self {
            bits,
            half: Ball::from_rational(&Rational::one().mul_pow2(-1), bits),
            origin,
            direction,
            exact,
            nets,
            along,
            shadow,
            enclosed_shadow,
        })
    }

    fn world(&self, parameter: &Real) -> Point3 {
        self.origin.plus(&self.direction.scaled(parameter))
    }

    // Every crossing of the closed curve, once. A convex loop has two or none.
    fn roots(&self, dec: &mut Decider) -> Result<Vec<(usize, Root)>, Refusal> {
        let count = self.exact.len();
        let mut joins = Vec::with_capacity(count);
        for (exact, net) in self.exact.iter().zip(&self.nets) {
            joins.push(match strict(&net[0]) {
                Some(sign) => sign,
                None => math::sign(&exact[0], dec)?,
            });
        }
        let mut roots = Vec::with_capacity(2);
        let mut cells = 0;
        for span in 0..count {
            let previous = (span + count - 1) % count;
            let ends = [joins[span], joins[(span + 1) % count]];
            if ends[0] == Sign::Zero {
                self.transverse_join(previous, span, dec)?;
                roots.push((span, Root::Exact(Rational::zero())));
            }
            self.span_roots(span, ends, &mut cells, &mut roots, dec)?;
        }
        if !matches!(roots.len(), 0 | 2) {
            return Err(unresolved(
                "a line contact on this convex section is not a separated crossing pair",
            ));
        }
        Ok(roots)
    }

    // The curve is C1 at a join, so a crossing there keeps one slope sign.
    fn transverse_join(
        &self,
        previous: usize,
        span: usize,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let arriving = &self.exact[previous];
        let leaving = &self.exact[span];
        let before = math::sign(&(&arriving[3] - &arriving[2]), dec)?;
        let after = math::sign(&(&leaving[1] - &leaving[0]), dec)?;
        if before == Sign::Zero || before != after {
            return Err(unresolved(
                "a line touching this spline at a join is not a transverse crossing",
            ));
        }
        Ok(())
    }

    fn span_roots(
        &self,
        span: usize,
        ends: [Sign; 2],
        cells: &mut usize,
        roots: &mut Vec<(usize, Root)>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let mut pending = vec![Cell {
            from: Rational::zero(),
            to: Rational::one(),
            net: self.nets[span].clone(),
            ends,
            depth: 0,
        }];
        while let Some(cell) = pending.pop() {
            dec.check_deadline()?;
            *cells += 1;
            if settled(&cell.net) {
                continue;
            }
            if cell.monotone() {
                if opposite(cell.ends[0], cell.ends[1]) {
                    let [from, to] = self.narrowed(cell, dec)?;
                    roots.push((span, Root::Between(from, to)));
                }
                continue;
            }
            if cell.depth >= MAX_DEPTH || *cells > MAX_CELLS {
                return Err(unresolved(
                    "a line contact on this finite spline is not a separated transverse crossing",
                ));
            }
            let [left, right] = split_bernstein(cell.net.clone(), &self.half);
            let middle = (&cell.from + &cell.to).mul_pow2(-1);
            let sign = match strict(&left[3]) {
                Some(sign) => sign,
                None => {
                    let at = Real::from_rational(middle.clone());
                    math::sign(&bernstein(&self.exact[span], &at), dec)?
                }
            };
            if sign == Sign::Zero {
                roots.push((span, Root::Exact(middle.clone())));
            }
            pending.push(Cell {
                from: middle.clone(),
                to: cell.to,
                net: right,
                ends: [sign, cell.ends[1]],
                depth: cell.depth + 1,
            });
            pending.push(Cell {
                from: cell.from,
                to: middle,
                net: left,
                ends: [cell.ends[0], sign],
                depth: cell.depth + 1,
            });
        }
        Ok(())
    }

    // Cuts a one-crossing cell down, far enough for its root certificate.
    fn narrowed(&self, mut cell: Cell, dec: &mut Decider) -> Result<[Rational; 2], Refusal> {
        for _ in 0..MAX_NARROWING {
            dec.check_deadline()?;
            let [left, right] = split_bernstein(cell.net.clone(), &self.half);
            let Some(sign) = strict(&left[3]) else {
                break;
            };
            let middle = (&cell.from + &cell.to).mul_pow2(-1);
            if opposite(cell.ends[0], sign) {
                cell.to = middle;
                cell.ends[1] = sign;
                cell.net = left;
            } else {
                cell.from = middle;
                cell.ends[0] = sign;
                cell.net = right;
            }
        }
        Ok([cell.from, cell.to])
    }

    // The line parameter of a crossing, from its span parameter.
    fn parameter(&self, span: usize, root: &Root, dec: &mut Decider) -> Result<Real, Refusal> {
        let along = &self.along[span];
        match root {
            // A join is a control itself, with no cubic to expand around it.
            Root::Exact(at) if at.is_zero() => Ok(&along[0] / &self.shadow),
            Root::Exact(at) => {
                Ok(&bernstein(along, &Real::from_rational(at.clone())) / &self.shadow)
            }
            Root::Between(from, to) => self.certified(span, [from, to], dec),
        }
    }

    // The span's cubic and the line parameter it fixes, as one rational system,
    // so the parameter is a leaf and a hit is exactly on the line that made it.
    fn certified(
        &self,
        span: usize,
        cell: [&Rational; 2],
        dec: &mut Decider,
    ) -> Result<Real, Refusal> {
        let cubic = |row: &[Real; 4]| {
            let [a, b, c, d] = row;
            let three = Real::from_i64(3);
            let powers = [
                a.clone(),
                &three * &(b - a),
                &three * &(&(a - b) + &(c - b)),
                &(d - a) + &(&three * &(b - c)),
            ];
            RealPoly::from_terms(2, (0..).zip(powers).map(|(i, c)| (vec![i, 0], c)).collect())
        };
        let equations = [
            cubic(&self.exact[span]),
            RealPoly::var(2, 1).scale(&self.shadow) - cubic(&self.along[span]),
        ];
        let row = self.along[span]
            .each_ref()
            .map(|value| dec.enclose_budgeted_at(value, self.bits));
        let [a, b, c, d] = row;
        let hull = a?.hull(&b?).hull(&c?).hull(&d?);
        let line = hull
            .checked_div(&self.enclosed_shadow)
            .ok_or_else(precision_refusal)?;
        let bounds = [[cell[0].clone(), cell[1].clone()], [line.lo(), line.hi()]];
        let solved = RealPoly::simple_root_in(&equations, &bounds, self.bits, dec)?;
        let mut solved = solved
            .ok_or_else(|| unresolved("this finite spline crossing has no rational equations"))?;
        Ok(solved.swap_remove(1))
    }
}

fn bernstein(row: &[Real; 4], at: &Real) -> Real {
    let rest = &Real::one() - at;
    let three = Real::from_i64(3);
    let basis = [
        &rest.square() * &rest,
        &three * &(&rest.square() * at),
        &three * &(&rest * &at.square()),
        &at.square() * at,
    ];
    (0..4).fold(Real::zero(), |sum, i| &sum + &(&row[i] * &basis[i]))
}
