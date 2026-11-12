//! Certified surface area of a whole rational B-spline patch.
//!
//! Each knot span is cut into equal cells. On a cell with midpoint `m`, the
//! area element `f = |S_u x S_v|` is `f(m)` plus a linear term, which
//! integrates to zero over the symmetric cell, plus a Lagrange remainder. With
//! cell sides `a` and `b` the remainder integrates to at most
//! `|cell| (sup|f_uu| a^2 / 24 + sup|f_uv| a b / 16 + sup|f_vv| b^2 / 24)`.
//! The second derivatives of `f` come from interval enclosures of the
//! surface's derivatives up to third order over the cell, so the enclosure
//! narrows with the square of the cell size. A cell where the area element is
//! not separated from zero is halved a fixed number of times; one still not
//! separated is enclosed between zero and its largest area element, because
//! `f` is a norm and its integral over the cell is at most the cell's measure
//! times that maximum. This covers a patch whose area element vanishes on a
//! boundary curve, such as a rolling envelope's mouth row.

use super::*;

// Derivative multi-indices up to order three, in this order.
const ORDERS: [(usize, usize); 10] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (2, 0),
    (1, 1),
    (0, 2),
    (3, 0),
    (2, 1),
    (1, 2),
    (0, 3),
];

// Halvings of a cell whose enclosures do not separate the area element from
// zero. A cell still unseparated at this depth is enclosed between zero and
// its largest area element, which bounds its integral without dividing by a
// value that may be zero (a rolling envelope's mouth row leaves the area
// element equal to zero on that boundary curve).
const MAX_SPLITS: u32 = 3;

fn index(order: (usize, usize)) -> usize {
    ORDERS
        .iter()
        .position(|&o| o == order)
        .expect("order up to three")
}

fn choose(n: usize, k: usize) -> i64 {
    (0..k).fold(1i64, |c, i| c * (n - i) as i64 / (i as i64 + 1))
}

type Vector = [Ball; 3];

fn cross(a: &Vector, b: &Vector) -> Vector {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

fn dot(a: &Vector, b: &Vector) -> Ball {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}

fn add(a: &Vector, b: &Vector) -> Vector {
    std::array::from_fn(|k| &a[k] + &b[k])
}

fn scale(a: &Vector, factor: &Ball) -> Vector {
    std::array::from_fn(|k| &a[k] * factor)
}

fn magnitude(value: &Ball) -> Rational {
    value.lo().abs().max(value.hi().abs())
}

// A homogeneous derivative net with its controls enclosed once.
struct BallNet {
    net: Net,
    controls: Vec<[Ball; 4]>,
}

impl BallNet {
    fn new(net: Net, precision: u32, dec: &Decider) -> Result<Self, Refusal> {
        let controls = net
            .controls
            .iter()
            .map(|point| {
                dec.check_deadline()?;
                Ok([
                    dec.enclose_at(&point[0], precision)?,
                    dec.enclose_at(&point[1], precision)?,
                    dec.enclose_at(&point[2], precision)?,
                    dec.enclose_at(&point[3], precision)?,
                ])
            })
            .collect::<Result<_, Refusal>>()?;
        Ok(Self { net, controls })
    }

    // Encloses the net over an exact parameter box inside one knot span of
    // the surface. Spans of a derivative net cover the same intervals.
    fn enclose(
        &self,
        cell: [[&Rational; 2]; 2],
        precision: u32,
        dec: &Decider,
    ) -> Result<[Ball; 4], Refusal> {
        let net = &self.net;
        let span = |axis: usize| {
            net.spans(axis)
                .find(|&i| {
                    net.knots[axis][i] <= *cell[axis][0] && *cell[axis][1] <= net.knots[axis][i + 1]
                })
                .ok_or_else(|| invalid("an area cell leaves its knot span"))
        };
        let span = [span(0)?, span(1)?];
        let mut columns = Vec::with_capacity(net.degrees[1] + 1);
        let mut row = Vec::with_capacity(net.degrees[0] + 1);
        for j in span[1] - net.degrees[1]..=span[1] {
            row.clear();
            for i in span[0] - net.degrees[0]..=span[0] {
                row.push(self.controls[i * net.count[1] + j].clone());
            }
            columns.push(de_boor(
                &mut row,
                &net.knots[0],
                span[0],
                cell[0],
                precision,
                dec,
            )?);
        }
        de_boor(
            &mut columns,
            &net.knots[1],
            span[1],
            cell[1],
            precision,
            dec,
        )
    }
}

// Every homogeneous derivative net up to order three; `None` where a degree
// runs out, so that derivative is zero.
struct Nets {
    nets: [Option<BallNet>; 10],
    weight_lower: Rational,
}

impl Nets {
    fn new(surface: &BSplineSurface, precision: u32, dec: &Decider) -> Result<Self, Refusal> {
        let [p, u, v, uu, uv, vv] = &surface.data.jets;
        let third = |net: &Option<Net>, axis: usize| net.as_ref().and_then(|n| n.derivative(axis));
        let nets = [
            p.clone(),
            u.clone(),
            v.clone(),
            uu.clone(),
            uv.clone(),
            vv.clone(),
            third(uu, 0),
            third(uu, 1),
            third(uv, 1),
            third(vv, 1),
        ];
        let mut enclosed: [Option<BallNet>; 10] = std::array::from_fn(|_| None);
        for (slot, net) in enclosed.iter_mut().zip(nets) {
            if let Some(net) = net {
                *slot = Some(BallNet::new(net, precision, dec)?);
            }
        }
        Ok(Self {
            nets: enclosed,
            weight_lower: surface.data.weight_lower.clone(),
        })
    }

    // Cartesian derivatives up to `order` over the cell, by Leibniz's rule on
    // `H = w S`.
    fn derivatives(
        &self,
        cell: [[&Rational; 2]; 2],
        order: usize,
        precision: u32,
        dec: &Decider,
    ) -> Result<Vec<Vector>, Refusal> {
        let count = ORDERS.iter().filter(|(a, b)| a + b <= order).count();
        let mut homogeneous = Vec::with_capacity(count);
        for net in &self.nets[..count] {
            homogeneous.push(match net {
                Some(net) => net.enclose(cell, precision, dec)?,
                None => std::array::from_fn(|_| zero(precision)),
            });
        }
        let weight = &homogeneous[0][3];
        let weight = bounds(
            &weight.lo().max(self.weight_lower.clone()),
            &weight.hi(),
            precision,
        );
        let mut cartesian: Vec<Vector> = Vec::with_capacity(count);
        for (position, &(a, b)) in ORDERS[..count].iter().enumerate() {
            let mut value: Vector = std::array::from_fn(|k| homogeneous[position][k].clone());
            for i in 0..=a {
                for j in 0..=b {
                    if i == 0 && j == 0 {
                        continue;
                    }
                    let factor = Ball::from_rational(
                        &Rational::from_i64(choose(a, i) * choose(b, j)),
                        precision,
                    );
                    let w = &homogeneous[index((i, j))][3] * &factor;
                    let lower = &cartesian[index((a - i, b - j))];
                    value = std::array::from_fn(|k| &value[k] - &(&w * &lower[k]));
                }
            }
            let divided = value
                .iter()
                .map(|component| {
                    component.checked_div(&weight).ok_or_else(|| {
                        invalid("positive rational chart denominator was not separated")
                    })
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            cartesian.push(divided.try_into().expect("three components"));
        }
        Ok(cartesian)
    }

    // Enclosure of the area over one cell, or `None` where the cell's
    // enclosures do not separate the area element from zero.
    fn cell(
        &self,
        cell: [[&Rational; 2]; 2],
        precision: u32,
        dec: &Decider,
        coarse: bool,
    ) -> Result<Option<Ball>, Refusal> {
        let width = [cell[0][1] - cell[0][0], cell[1][1] - cell[1][0]];
        let mid = [
            (cell[0][0] + cell[0][1]) / Rational::from_i64(2),
            (cell[1][0] + cell[1][1]) / Rational::from_i64(2),
        ];
        let point =
            self.derivatives([[&mid[0], &mid[0]], [&mid[1], &mid[1]]], 1, precision, dec)?;
        let normal = cross(&point[1], &point[2]);
        let centre = dot(&normal, &normal)
            .sqrt()
            .ok_or_else(|| invalid("area element at a cell midpoint is undefined"))?;
        let s = self.derivatives(cell, 3, precision, dec)?;
        let d = |order: (usize, usize)| &s[index(order)];
        let n = cross(d((1, 0)), d((0, 1)));
        let n_u = add(&cross(d((2, 0)), d((0, 1))), &cross(d((1, 0)), d((1, 1))));
        let n_v = add(&cross(d((1, 1)), d((0, 1))), &cross(d((1, 0)), d((0, 2))));
        let two = Ball::from_rational(&Rational::from_i64(2), precision);
        let n_uu = add(
            &add(
                &cross(d((3, 0)), d((0, 1))),
                &scale(&cross(d((2, 0)), d((1, 1))), &two),
            ),
            &cross(d((1, 0)), d((2, 1))),
        );
        let n_uv = add(
            &add(&cross(d((2, 1)), d((0, 1))), &cross(d((2, 0)), d((0, 2)))),
            &cross(d((1, 0)), d((1, 2))),
        );
        let n_vv = add(
            &add(
                &cross(d((1, 2)), d((0, 1))),
                &scale(&cross(d((1, 1)), d((0, 2))), &two),
            ),
            &cross(d((1, 0)), d((0, 3))),
        );
        let squared = dot(&n, &n);
        if squared.lo() <= Rational::zero() {
            if !coarse {
                return Ok(None);
            }
            // The area element is not separated from zero over this cell, so
            // no second-derivative bound divides by it. The area integral is
            // still between zero and the cell's measure times the element's
            // largest enclosure: a valid enclosure that splitting has already
            // narrowed as far as the cell budget allows.
            let f = squared.sqrt().expect("nonnegative enclosure");
            let measure = Ball::from_rational(&(&width[0] * &width[1]), precision);
            let slack = bounds(&Rational::zero(), &f.hi(), precision);
            return Ok(Some(&measure * &slack));
        }
        let f = squared.sqrt().expect("positive enclosure");
        let f3 = &(&f * &f) * &f;
        // f_ab = (N_a . N_b + N . N_ab) / f - (N . N_a)(N . N_b) / f^3.
        let second = |a: &Vector, b: &Vector, ab: &Vector| -> Result<Rational, Refusal> {
            let first = (&dot(a, b) + &dot(&n, ab))
                .checked_div(&f)
                .ok_or_else(|| invalid("area element is not separated from zero"))?;
            let correction = (&dot(&n, a) * &dot(&n, b))
                .checked_div(&f3)
                .ok_or_else(|| invalid("area element is not separated from zero"))?;
            Ok(magnitude(&(&first - &correction)))
        };
        let f_uu = second(&n_u, &n_u, &n_uu)?;
        let f_uv = second(&n_u, &n_v, &n_uv)?;
        let f_vv = second(&n_v, &n_v, &n_vv)?;
        let uu = &(&f_uu * &(&width[0] * &width[0])) / &Rational::from_i64(24);
        let uv = &(&(&f_uv * &width[0]) * &width[1]) / &Rational::from_i64(16);
        let vv = &(&f_vv * &(&width[1] * &width[1])) / &Rational::from_i64(24);
        let remainder = &(&uu + &uv) + &vv;
        let measure = Ball::from_rational(&(&width[0] * &width[1]), precision);
        let slack = bounds(&(&Rational::zero() - &remainder), &remainder, precision);
        Ok(Some(&measure * &(&centre + &slack)))
    }
}

impl BSplineSurface {
    /// Encloses the surface area of the whole patch, cutting each knot span
    /// into `cells` equal pieces along u and v. A cell whose area element is
    /// not separated from zero is halved a few times and then enclosed
    /// coarsely, so the result is always a valid enclosure.
    pub fn area_enclosure(
        &self,
        cells: [usize; 2],
        precision: u32,
        dec: &Decider,
    ) -> Result<Ball, Refusal> {
        let nets = Nets::new(self, precision, dec)?;
        self.area_with(&nets, cells, precision, dec)
    }

    /// Encloses the surface area of the whole patch, doubling the cells per
    /// span until the enclosure is at most `width` wide or more than
    /// `max_cells` cells would be needed; the last enclosure is returned either
    /// way, so the caller reads its width.
    pub fn certified_area(
        &self,
        width: &Rational,
        max_cells: usize,
        precision: u32,
        dec: &Decider,
    ) -> Result<Ball, Refusal> {
        let nets = Nets::new(self, precision, dec)?;
        let spans = [self.spans(0).count(), self.spans(1).count()];
        let mut cells = [1usize, 1usize];
        loop {
            let area = self.area_with(&nets, cells, precision, dec)?;
            let next = cells.map(|c| c * 2);
            if &(&area.hi() - &area.lo()) <= width
                || spans[0] * next[0] * spans[1] * next[1] > max_cells
            {
                return Ok(area);
            }
            cells = next;
        }
    }

    fn area_with(
        &self,
        nets: &Nets,
        cells: [usize; 2],
        precision: u32,
        dec: &Decider,
    ) -> Result<Ball, Refusal> {
        let pieces = |axis: usize| -> Vec<[Rational; 2]> {
            self.spans(axis)
                .flat_map(|[lo, hi]| {
                    let step = &(hi - lo) / &Rational::from_i64(cells[axis] as i64);
                    (0..cells[axis])
                        .map(|k| {
                            let a = lo + &(&step * &Rational::from_i64(k as i64));
                            let b = if k + 1 == cells[axis] {
                                hi.clone()
                            } else {
                                lo + &(&step * &Rational::from_i64(k as i64 + 1))
                            };
                            [a, b]
                        })
                        .collect::<Vec<_>>()
                })
                .collect()
        };
        let (us, vs) = (pieces(0), pieces(1));
        let cells: Vec<[[Rational; 2]; 2]> = us
            .iter()
            .flat_map(|u| vs.iter().map(move |v| [u.clone(), v.clone()]))
            .collect();
        let areas = enclose_cells(nets, &cells, precision, dec)?;
        // Summed in the cells' order, so the sum is the same on any number
        // of threads.
        Ok(areas
            .iter()
            .fold(zero(precision), |total, area| &total + area))
    }
}

// The area of one cell, halving where the area element is not yet separated
// from zero and enclosing a cell still unseparated at `MAX_SPLITS` between
// zero and its largest element. Depth-first in a fixed order, so the sum is
// the same on every run.
fn cell_area(
    nets: &Nets,
    cell: &[[Rational; 2]; 2],
    precision: u32,
    dec: &Decider,
) -> Result<Ball, Refusal> {
    let mut total = zero(precision);
    let mut pending = vec![(cell.clone(), 0u32)];
    while let Some(([u, v], depth)) = pending.pop() {
        dec.check_deadline()?;
        let last = depth >= MAX_SPLITS;
        if let Some(area) = nets.cell([[&u[0], &u[1]], [&v[0], &v[1]]], precision, dec, last)? {
            total = &total + &area;
            continue;
        }
        let half = |[a, b]: [Rational; 2]| {
            let m = &(&a + &b) / &Rational::from_i64(2);
            [[a, m.clone()], [m, b]]
        };
        for u in half(u) {
            for v in half(v.clone()) {
                pending.push(([u.clone(), v], depth + 1));
            }
        }
    }
    Ok(total)
}

// Every cell's area, in the cells' order. A cell's area depends on nothing
// but the cell and charges no work, so the kernel's helper threads take cells
// as they come and the writer reads the results in order; the first refusal
// in that order is the answer's.
fn enclose_cells(
    nets: &Nets,
    cells: &[[[Rational; 2]; 2]],
    precision: u32,
    dec: &Decider,
) -> Result<Vec<Ball>, Refusal> {
    let helpers = if cells.len() > 1 {
        catallax_numbers::Helpers::new(dec)
    } else {
        None
    };
    let Some(helpers) = helpers else {
        return cells
            .iter()
            .map(|cell| cell_area(nets, cell, precision, dec))
            .collect();
    };
    let slots: Vec<std::sync::OnceLock<Result<Ball, Refusal>>> = (0..cells.len())
        .map(|_| std::sync::OnceLock::new())
        .collect();
    let next = std::sync::atomic::AtomicUsize::new(0);
    let take = |dec: &Decider| {
        while !helpers.stopped() {
            let index = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let Some(cell) = cells.get(index) else {
                return;
            };
            let area = cell_area(nets, cell, precision, dec);
            if area.is_err() {
                // A refusal after the stop may be the stop's own; the writer
                // recomputes that cell below.
                if helpers.stopped() {
                    return;
                }
                let _ = slots[index].set(area);
                helpers.stop();
                return;
            }
            let _ = slots[index].set(area);
        }
    };
    std::thread::scope(|scope| {
        struct Stop<'h>(&'h catallax_numbers::Helpers);
        impl Drop for Stop<'_> {
            fn drop(&mut self) {
                self.0.stop();
            }
        }
        let _stop = Stop(&helpers);
        for _ in 0..helpers.count() {
            let spawned = std::thread::Builder::new()
                .name("catallax-spline-area".into())
                .stack_size(helpers.stack_size())
                .spawn_scoped(scope, || helpers.run(|worker| take(worker)));
            if spawned.is_err() {
                break;
            }
        }
        take(dec);
    });
    // A cell a stopped thread never reached is computed here, so the result
    // never depends on which thread stopped first.
    cells
        .iter()
        .zip(slots)
        .map(|(cell, slot)| match slot.into_inner() {
            Some(area) => area,
            None => cell_area(nets, cell, precision, dec),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Rational {
        Rational::from_i64(n)
    }

    fn pi() -> Ball {
        Ball::pi(256)
    }

    // A quarter circle of radius r about the origin in the plane z, as
    // rational quadratic controls and weights.
    fn arc(r: &Real, z: &Real) -> ([Point3; 3], [Real; 3]) {
        let zero = Real::zero();
        let half = Real::from_i64(2).sqrt() / Real::from_i64(2);
        (
            [
                Point3::new(r.clone(), zero.clone(), z.clone()),
                Point3::new(r.clone(), r.clone(), z.clone()),
                Point3::new(zero, r.clone(), z.clone()),
            ],
            [Real::one(), half, Real::one()],
        )
    }

    fn knots(degree: usize) -> Vec<Rational> {
        let mut knots = vec![q(0); degree + 1];
        knots.extend(vec![q(1); degree + 1]);
        knots
    }

    fn contains(area: &Ball, exact: &Ball) {
        assert!(
            area.lo() <= exact.lo() && exact.hi() <= area.hi(),
            "{area:?} misses {exact:?}"
        );
    }

    #[test]
    fn a_quarter_cylinder_encloses_its_area_and_narrows_with_the_cells() {
        let mut dec = Decider::default();
        let (r, length) = (Real::from_i64(2), Real::from_i64(3));
        let (bottom, w) = arc(&r, &Real::zero());
        let (top, _) = arc(&r, &length);
        let mut controls = Vec::new();
        let mut weights = Vec::new();
        for i in 0..3 {
            controls.extend([bottom[i].clone(), top[i].clone()]);
            weights.extend([w[i].clone(), w[i].clone()]);
        }
        let surface = BSplineSurface::new(
            [2, 1],
            [knots(2), knots(1)],
            [3, 2],
            controls,
            weights,
            &mut dec,
        )
        .unwrap();
        // r pi / 2 * length = 3 pi.
        let exact = &pi() * &Ball::from_rational(&q(3), 256);
        let coarse = surface.area_enclosure([1, 1], 128, &dec).unwrap();
        let fine = surface.area_enclosure([8, 8], 128, &dec).unwrap();
        contains(&coarse, &exact);
        contains(&fine, &exact);
        assert!(fine.rad() < coarse.rad(), "{fine:?} {coarse:?}");
        let target = Rational::from_ratio(1, 100).unwrap();
        let certified = surface.certified_area(&target, 1 << 14, 128, &dec).unwrap();
        contains(&certified, &exact);
        assert!(&certified.hi() - &certified.lo() <= target);
    }

    #[test]
    fn a_torus_patch_encloses_its_area() {
        let mut dec = Decider::default();
        // Tube quarter circle in the (radius, height) plane, swept by a
        // quarter turn about the z axis.
        let half = Real::from_i64(2).sqrt() / Real::from_i64(2);
        let tube = [
            (Real::from_i64(6), Real::zero(), Real::one()),
            (Real::from_i64(6), Real::one(), half.clone()),
            (Real::from_i64(5), Real::one(), Real::one()),
        ];
        let turn = [
            (Real::one(), Real::zero(), Real::one()),
            (Real::one(), Real::one(), half.clone()),
            (Real::zero(), Real::one(), Real::one()),
        ];
        let mut controls = Vec::new();
        let mut weights = Vec::new();
        for (cx, sx, wu) in &turn {
            for (radius, height, wv) in &tube {
                controls.push(Point3::new(radius * cx, radius * sx, height.clone()));
                weights.push(wu * wv);
            }
        }
        let surface = BSplineSurface::new(
            [2, 2],
            [knots(2), knots(2)],
            [3, 3],
            controls,
            weights,
            &mut dec,
        )
        .unwrap();
        // r (pi / 2) (R pi / 2 + r) with R = 5, r = 1.
        let quarter = &pi() * &Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), 256);
        let exact = &quarter
            * &(&(&quarter * &Ball::from_rational(&q(5), 256)) + &Ball::from_rational(&q(1), 256));
        let target = Rational::from_ratio(1, 100).unwrap();
        let area = surface.certified_area(&target, 1 << 14, 128, &dec).unwrap();
        contains(&area, &exact);
        assert!(&area.hi() - &area.lo() <= target);
    }
}
