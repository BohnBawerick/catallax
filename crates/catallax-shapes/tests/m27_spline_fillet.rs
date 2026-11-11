//! Opt-in M27 scope experiment. No production intersection entry calls this search.
use catallax_numbers::{Decider, MPoly, Rational, Refusal, with_operation_budget};

const N: usize = 6;
const M: usize = N - 1;
const ALLOWANCE: u64 = 20_000_000;
// Keep the M10 proof caps as well as its application allowance.
const MAX_DEPTH: u32 = 24;
const MAX_CELLS: usize = 16_384;
type Bounds = [[Rational; 2]; N];
type Interval = [Rational; 2];
type Matrix = Vec<Vec<Rational>>;

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn arithmetic(
    a: &Rational,
    b: &Rational,
    op: char,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, op == '*' || op == '/')?;
    Ok(match op {
        '+' => a + b,
        '-' => a - b,
        '*' => a * b,
        '/' => a / b,
        _ => unreachable!(),
    })
}
fn interval_add(a: &Interval, b: &Interval, dec: &mut Decider) -> Result<Interval, Refusal> {
    Ok([
        arithmetic(&a[0], &b[0], '+', dec)?,
        arithmetic(&a[1], &b[1], '+', dec)?,
    ])
}
fn interval_mul(a: &Interval, b: &Interval, dec: &mut Decider) -> Result<Interval, Refusal> {
    let mut values = Vec::new();
    for x in a {
        for y in b {
            values.push(arithmetic(x, y, '*', dec)?);
        }
    }
    Ok([
        values.iter().min().unwrap().clone(),
        values.iter().max().unwrap().clone(),
    ])
}
fn point(x: Rational) -> Interval {
    [x.clone(), x]
}
fn midpoint(b: &Interval, dec: &mut Decider) -> Result<Rational, Refusal> {
    arithmetic(&arithmetic(&b[0], &b[1], '+', dec)?, &q(2), '/', dec)
}

// Same first two polygon sides as LOFT-CASES' nonaffine parallel-section case.
// Parameter order is u,v,s,t,a,b. Both normals point out of the loft.
struct Definition {
    equations: [MPoly; M],
    first: [MPoly; 3],
    second: [MPoly; 3],
    na: [MPoly; 3],
    nb: [MPoly; 3],
}
fn definition() -> Definition {
    let [u, v, s, t, a, b] = std::array::from_fn(|i| MPoly::var(N, i));
    let c = |x| MPoly::constant(N, q(x));
    let first = [
        &(&c(-1) + &u.scale(&q(2))) + &(&v * &(&c(1) - &u)),
        &c(-1) + &(&u * &v),
        v.scale(&q(2)),
    ];
    let second = [
        &c(1) - &(&s * &t),
        &(&c(-1) + &s.scale(&q(2))) + &t,
        t.scale(&q(2)),
    ];
    let na = [
        v.scale(&q(2)),
        &v.scale(&q(2)) - &c(4),
        &u.scale(&q(2)) - &v,
    ];
    let nb = [c(4), t.scale(&q(2)), &s.scale(&q(2)) - &t];
    let mut equations = std::array::from_fn(|i| {
        if i < 3 {
            &(&(&first[i] - &second[i]) * &(&a * &b)).scale(&q(10))
                - &(&(&na[i] * &b) - &(&nb[i] * &a))
        } else {
            MPoly::zero(N)
        }
    });
    equations[3] = &(&a * &a) - &na.iter().fold(MPoly::zero(N), |sum, x| &sum + &(x * x));
    equations[4] = &(&b * &b) - &nb.iter().fold(MPoly::zero(N), |sum, x| &sum + &(x * x));
    Definition {
        equations,
        first,
        second,
        na,
        nb,
    }
}

#[derive(Clone)]
struct Net {
    count: [usize; N],
    values: Vec<Rational>,
}
impl Net {
    fn index(&self, mut flat: usize) -> [usize; N] {
        let mut at = [0; N];
        for axis in (0..N).rev() {
            at[axis] = flat % self.count[axis];
            flat /= self.count[axis];
        }
        at
    }
    fn stride(&self, axis: usize) -> usize {
        self.count[axis + 1..].iter().product()
    }
    // All degrees here are at most two. Convert each power on an absolute
    // interval to Bernstein coefficients, then take their tensor product.
    fn new(poly: &MPoly, bounds: &Bounds, dec: &mut Decider) -> Result<Self, Refusal> {
        let terms = poly.terms();
        let count = std::array::from_fn(|i| {
            1 + terms
                .iter()
                .map(|(m, _)| m.exponents()[i] as usize)
                .max()
                .unwrap_or(0)
        });
        assert!(count.iter().all(|&n| n <= 3));
        let mut net = Self {
            count,
            values: vec![q(0); count.iter().product()],
        };
        for (monomial, coefficient) in terms {
            let mut factors: Vec<Vec<Rational>> = Vec::new();
            for axis in 0..N {
                let [lo, hi] = &bounds[axis];
                let degree = count[axis] - 1;
                let power = monomial.exponents()[axis];
                let f = match power {
                    0 => vec![q(1); degree + 1],
                    1 if degree == 1 => vec![lo.clone(), hi.clone()],
                    1 => vec![lo.clone(), midpoint(&bounds[axis], dec)?, hi.clone()],
                    2 => vec![
                        arithmetic(lo, lo, '*', dec)?,
                        arithmetic(lo, hi, '*', dec)?,
                        arithmetic(hi, hi, '*', dec)?,
                    ],
                    _ => unreachable!(),
                };
                factors.push(f);
            }
            for k in 0..net.values.len() {
                let at = net.index(k);
                let mut term = coefficient.clone();
                for axis in 0..N {
                    if monomial.exponents()[axis] != 0 {
                        term = arithmetic(&term, &factors[axis][at[axis]], '*', dec)?;
                    }
                }
                net.values[k] = arithmetic(&net.values[k], &term, '+', dec)?;
            }
        }
        Ok(net)
    }
    fn range(&self, dec: &mut Decider) -> Result<Interval, Refusal> {
        dec.charge_work(2 * self.values.len() as u64)?;
        Ok([
            self.values.iter().min().unwrap().clone(),
            self.values.iter().max().unwrap().clone(),
        ])
    }
    fn derivative_range(
        &self,
        axis: usize,
        bounds: &Bounds,
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        if self.count[axis] == 1 {
            return Ok(point(q(0)));
        }
        let width = arithmetic(&bounds[axis][1], &bounds[axis][0], '-', dec)?;
        let factor = arithmetic(&q((self.count[axis] - 1) as i64), &width, '/', dec)?;
        let mut values = Vec::new();
        for k in 0..self.values.len() {
            if self.index(k)[axis] + 1 < self.count[axis] {
                let diff = arithmetic(
                    &self.values[k + self.stride(axis)],
                    &self.values[k],
                    '-',
                    dec,
                )?;
                values.push(arithmetic(&diff, &factor, '*', dec)?);
            }
        }
        dec.charge_work(2 * values.len() as u64)?;
        Ok([
            values.iter().min().unwrap().clone(),
            values.iter().max().unwrap().clone(),
        ])
    }
    fn halves(&self, axis: usize, dec: &mut Decider) -> Result<[Self; 2], Refusal> {
        let mut left = self.clone();
        let mut right = self.clone();
        let degree = self.count[axis] - 1;
        let stride = self.stride(axis);
        for k in 0..self.values.len() {
            if self.index(k)[axis] != 0 {
                continue;
            }
            let mut work = (0..=degree)
                .map(|i| self.values[k + i * stride].clone())
                .collect::<Vec<_>>();
            for level in 1..=degree {
                for i in 0..=degree - level {
                    work[i] = midpoint(&[work[i].clone(), work[i + 1].clone()], dec)?;
                }
                left.values[k + level * stride] = work[0].clone();
                right.values[k + (degree - level) * stride] = work[degree - level].clone();
            }
        }
        Ok([left, right])
    }
    fn centered_range(&self, free: usize, dec: &mut Decider) -> Result<Interval, Refusal> {
        let mut net = self.clone();
        for axis in 0..N {
            if axis == free {
                continue;
            }
            let left = net.halves(axis, dec)?[0].clone();
            let stride = net.stride(axis);
            let mut values = Vec::new();
            for k in 0..net.values.len() {
                if net.index(k)[axis] == 0 {
                    values.push(left.values[k + (net.count[axis] - 1) * stride].clone());
                }
            }
            net.count[axis] = 1;
            net.values = values;
        }
        net.range(dec)
    }
}

fn inverse(mut a: Matrix, dec: &mut Decider) -> Result<Option<Matrix>, Refusal> {
    let mut b: Matrix = (0..M)
        .map(|i| (0..M).map(|j| q(i64::from(i == j))).collect())
        .collect();
    for k in 0..M {
        dec.charge_work(M as u64)?;
        let Some(pivot) = (k..M).find(|&i| !a[i][k].is_zero()) else {
            return Ok(None);
        };
        a.swap(k, pivot);
        b.swap(k, pivot);
        let d = a[k][k].clone();
        for j in 0..M {
            a[k][j] = arithmetic(&a[k][j], &d, '/', dec)?;
            b[k][j] = arithmetic(&b[k][j], &d, '/', dec)?;
        }
        for i in 0..M {
            if i == k {
                continue;
            }
            let factor = a[i][k].clone();
            for j in 0..M {
                a[i][j] = arithmetic(
                    &a[i][j],
                    &arithmetic(&factor, &a[k][j], '*', dec)?,
                    '-',
                    dec,
                )?;
                b[i][j] = arithmetic(
                    &b[i][j],
                    &arithmetic(&factor, &b[k][j], '*', dec)?,
                    '-',
                    dec,
                )?;
            }
        }
    }
    Ok(Some(b))
}

#[derive(Clone)]
struct Cell {
    nets: Vec<Net>,
    bounds: Bounds,
    depth: u32,
}
// A parameter-dependent contraction proves one regular graph through the
// two free-axis faces, with no other zero in this cell. It needs no seed.
fn graph(cell: &Cell, dec: &mut Decider) -> Result<bool, Refusal> {
    let mut jacobian = Vec::new();
    for net in &cell.nets {
        let mut row = Vec::new();
        for axis in 0..N {
            row.push(net.derivative_range(axis, &cell.bounds, dec)?);
        }
        jacobian.push(row);
    }
    // Stable axis order, as in M10's fixed solved-axis preference.
    for free in (0..N).rev() {
        let solved: Vec<_> = (0..N).filter(|&i| i != free).collect();
        let mut mid_jac = Vec::new();
        for row in &jacobian {
            let mut mid_row = Vec::new();
            for &axis in &solved {
                mid_row.push(midpoint(&row[axis], dec)?);
            }
            mid_jac.push(mid_row);
        }
        let Some(preconditioner) = inverse(mid_jac, dec)? else {
            continue;
        };
        let mut residual = Vec::new();
        for net in &cell.nets {
            residual.push(net.centered_range(free, dec)?);
        }
        let mut certified = true;
        for (i, &axis) in solved.iter().enumerate() {
            let center = midpoint(&cell.bounds[axis], dec)?;
            let mut image = point(center.clone());
            for (k, f) in residual.iter().enumerate() {
                image = interval_add(
                    &image,
                    &interval_mul(&point(-&preconditioner[i][k]), f, dec)?,
                    dec,
                )?;
            }
            let mut norm = q(0);
            for (j, &other) in solved.iter().enumerate() {
                let mut r = point(q(i64::from(i == j)));
                for (k, row) in jacobian.iter().enumerate() {
                    r = interval_add(
                        &r,
                        &interval_mul(&point(-&preconditioner[i][k]), &row[other], dec)?,
                        dec,
                    )?;
                }
                let magnitude = r[0].abs().max(r[1].abs());
                norm = arithmetic(&norm, &magnitude, '+', dec)?;
                let mid = midpoint(&cell.bounds[other], dec)?;
                let delta = [
                    arithmetic(&cell.bounds[other][0], &mid, '-', dec)?,
                    arithmetic(&cell.bounds[other][1], &mid, '-', dec)?,
                ];
                image = interval_add(&image, &interval_mul(&r, &delta, dec)?, dec)?;
            }
            dec.charge_work(3)?;
            if norm >= q(1) || image[0] <= cell.bounds[axis][0] || image[1] >= cell.bounds[axis][1]
            {
                certified = false;
                break;
            }
        }
        if certified {
            return Ok(true);
        }
    }
    Ok(false)
}

#[derive(Debug, Default, PartialEq, Eq)]
struct Record {
    visited: usize,
    empty: usize,
    arcs: usize,
    pending: usize,
    depth: u32,
    setup_work: u64,
    used: u64,
    exhausted: bool,
    refusal: Option<(String, String)>,
    unresolved_box: Option<[[String; 2]; N]>,
}
fn run() -> Record {
    let mut record = Record::default();
    let (result, work) = with_operation_budget(ALLOWANCE, || -> Result<(), Refusal> {
        let mut dec = Decider::default();
        let equations = definition().equations;
        // ||Na|| >= 2, ||Nb|| >= 4, and both squared lengths <= 24.
        // These positive intervals exclude the extraneous negative roots.
        let bounds: Bounds = [
            [q(0), q(1)],
            [q(0), q(1)],
            [q(0), q(1)],
            [q(0), q(1)],
            [q(2), q(5)],
            [q(4), q(5)],
        ];
        let mut nets = Vec::new();
        for f in &equations {
            nets.push(Net::new(f, &bounds, &mut dec)?);
        }
        record.setup_work = catallax_numbers::current_operation_work().used;
        let mut pending = vec![Cell {
            nets,
            bounds,
            depth: 0,
        }];
        while let Some(cell) = pending.pop() {
            record.visited += 1;
            record.pending = pending.len() + 1;
            record.depth = record.depth.max(cell.depth);
            dec.check_deadline()?;
            let mut empty = false;
            for net in &cell.nets {
                let range = net.range(&mut dec)?;
                if range[0] > q(0) || range[1] < q(0) {
                    empty = true;
                    break;
                }
            }
            if empty {
                record.empty += 1;
                continue;
            }
            if graph(&cell, &mut dec)? {
                record.arcs += 1;
                continue;
            }
            let code = if cell.depth >= MAX_DEPTH {
                Some("experiment.m27-depth-limit")
            } else if record.empty + record.arcs + pending.len() + 2 > MAX_CELLS {
                Some("experiment.m27-cell-limit")
            } else {
                None
            };
            if let Some(code) = code {
                record.unresolved_box = Some(
                    cell.bounds
                        .each_ref()
                        .map(|b| b.each_ref().map(ToString::to_string)),
                );
                return Err(Refusal::new(
                    code,
                    "no empty-cell or regular-graph certificate at the unchanged M10 proof cap",
                ));
            }
            let mut axis = 0;
            let mut width = arithmetic(&cell.bounds[0][1], &cell.bounds[0][0], '-', &mut dec)?;
            for i in 1..N {
                let candidate = arithmetic(&cell.bounds[i][1], &cell.bounds[i][0], '-', &mut dec)?;
                if candidate > width {
                    axis = i;
                    width = candidate;
                }
            }
            let mid = midpoint(&cell.bounds[axis], &mut dec)?;
            let mut left = Cell {
                nets: Vec::new(),
                bounds: cell.bounds.clone(),
                depth: cell.depth + 1,
            };
            let mut right = left.clone();
            left.bounds[axis][1] = mid.clone();
            right.bounds[axis][0] = mid;
            for net in &cell.nets {
                let [a, b] = net.halves(axis, &mut dec)?;
                left.nets.push(a);
                right.nets.push(b);
            }
            pending.push(right);
            pending.push(left);
        }
        record.pending = 0;
        // A complete cell cover would still need exact endpoint joins.
        if record.arcs != 0 {
            return Err(Refusal::new(
                "experiment.m27-joins-unproved",
                "regular graph cells do not yet carry exact shared endpoint witnesses",
            ));
        }
        Ok(())
    });
    record.used = work.used;
    record.exhausted = work.exhausted;
    record.refusal = result.err().map(|r| (r.code.to_string(), r.detail));
    record
}

fn check_net_value(poly: &MPoly, net: &Net, bounds: &Bounds, sample: &[Rational; N]) {
    let at: [Rational; N] =
        std::array::from_fn(|i| (&sample[i] - &bounds[i][0]) / (&bounds[i][1] - &bounds[i][0]));
    let mut value = q(0);
    for (k, coefficient) in net.values.iter().enumerate() {
        let index = net.index(k);
        let mut term = coefficient.clone();
        for axis in 0..N {
            let t = &at[axis];
            let complement = &q(1) - t;
            let basis = match (net.count[axis], index[axis]) {
                (1, _) => q(1),
                (2, 0) => complement,
                (2, 1) => t.clone(),
                (3, 0) => &complement * &complement,
                (3, 1) => q(2) * (t * &complement),
                (3, 2) => t * t,
                _ => unreachable!(),
            };
            term = term * basis;
        }
        value = value + term;
    }
    assert_eq!(value, poly.eval(sample));
}

#[test]
#[ignore = "M27 scope experiment; excluded from ordinary tests and scoreboard"]
fn offset_spine_scope() {
    // Check the definition against independent bilinear interpolation of the
    // two existing loft sides, and the normals against their cross products.
    let Definition {
        equations,
        first,
        second,
        na,
        nb,
    } = definition();
    let cross = |p: &[MPoly; 3]| -> [MPoly; 3] {
        let du = p
            .each_ref()
            .map(|f| f.partial(if p == &first { 0 } else { 2 }));
        let dv = p
            .each_ref()
            .map(|f| f.partial(if p == &first { 1 } else { 3 }));
        std::array::from_fn(|i| {
            &(&du[(i + 1) % 3] * &dv[(i + 2) % 3]) - &(&du[(i + 2) % 3] * &dv[(i + 1) % 3])
        })
    };
    assert_eq!(cross(&first), na);
    assert_eq!(cross(&second), nb);
    let sample = [
        q(1) / q(3),
        q(2) / q(5),
        q(1) / q(4),
        q(3) / q(7),
        q(3),
        q(4),
    ];
    let lerp = |x: Rational, y: Rational, t: &Rational| &x * &(&q(1) - t) + &y * t;
    let a0 = [q(-1), q(-1), q(0)];
    let a1 = [q(1), q(-1), q(0)];
    let a2 = [q(0), q(-1), q(2)];
    let a3 = [q(1), q(0), q(2)];
    let b0 = a1.clone();
    let b1 = [q(1), q(1), q(0)];
    let b2 = a3.clone();
    let b3 = [q(0), q(2), q(2)];
    for i in 0..3 {
        assert_eq!(
            first[i].eval(&sample),
            lerp(
                lerp(a0[i].clone(), a1[i].clone(), &sample[0]),
                lerp(a2[i].clone(), a3[i].clone(), &sample[0]),
                &sample[1]
            )
        );
        assert_eq!(
            second[i].eval(&sample),
            lerp(
                lerp(b0[i].clone(), b1[i].clone(), &sample[2]),
                lerp(b2[i].clone(), b3[i].clone(), &sample[2]),
                &sample[3]
            )
        );
        let offset_difference = first[i].eval(&sample)
            - na[i].eval(&sample) / (&q(10) * &sample[4])
            - second[i].eval(&sample)
            + nb[i].eval(&sample) / (&q(10) * &sample[5]);
        assert_eq!(
            equations[i].eval(&sample),
            offset_difference * q(10) * sample[4].clone() * sample[5].clone()
        );
    }
    let (checks, work) = with_operation_budget(ALLOWANCE, || -> Result<(), Refusal> {
        let mut dec = Decider::default();
        let bounds = std::array::from_fn(|_| [q(0), q(1)]);
        let half = MPoly::constant(N, q(1) / q(2));
        let mut nets = Vec::new();
        for axis in 0..M {
            nets.push(Net::new(
                &(&MPoly::var(N, axis) - &half),
                &bounds,
                &mut dec,
            )?);
        }
        let mut cell = Cell {
            nets,
            bounds,
            depth: 0,
        };
        assert!(
            graph(&cell, &mut dec)?,
            "the full straight graph must certify"
        );
        cell.nets[0] = Net::new(&MPoly::constant(N, q(0)), &cell.bounds, &mut dec)?;
        assert!(
            !graph(&cell, &mut dec)?,
            "a rank-deficient sheet is not an arc"
        );
        let x = MPoly::var(N, 0);
        let square = &x * &x;
        let net = Net::new(&square, &cell.bounds, &mut dec)?;
        assert_eq!(net.values, vec![q(0), q(0), q(1)]);
        let [left, right] = net.halves(0, &mut dec)?;
        assert_eq!(left.values, vec![q(0), q(0), q(1) / q(4)]);
        assert_eq!(right.values, vec![q(1) / q(4), q(1) / q(2), q(1)]);
        assert_eq!(net.centered_range(5, &mut dec)?, point(q(1) / q(4)));
        let bounds = [
            [q(0), q(1)],
            [q(0), q(1)],
            [q(0), q(1)],
            [q(0), q(1)],
            [q(2), q(5)],
            [q(4), q(5)],
        ];
        for equation in &equations {
            let net = Net::new(equation, &bounds, &mut dec)?;
            check_net_value(equation, &net, &bounds, &sample);
            for axis in 0..N {
                let halves = net.halves(axis, &mut dec)?;
                let mut half_bounds = [bounds.clone(), bounds.clone()];
                let mid = (&bounds[axis][0] + &bounds[axis][1]) / q(2);
                half_bounds[0][axis][1] = mid.clone();
                half_bounds[1][axis][0] = mid;
                for half in 0..2 {
                    let center = std::array::from_fn(|i| {
                        (&half_bounds[half][i][0] + &half_bounds[half][i][1]) / q(2)
                    });
                    check_net_value(equation, &halves[half], &half_bounds[half], &center);
                }
            }
        }
        Ok(())
    });
    checks.unwrap();
    assert!(!work.exhausted);
    let expected = run();
    println!("M27 workers=1 {expected:?}");
    assert_eq!(
        expected.used, 2_910_656,
        "update the published scope if work changes"
    );
    assert_eq!(expected.setup_work, 24_550);
    assert!(!expected.exhausted);
    assert_eq!(
        (
            expected.visited,
            expected.empty,
            expected.arcs,
            expected.pending,
            expected.depth
        ),
        (102, 44, 0, 14, MAX_DEPTH)
    );
    assert_eq!(
        expected.refusal.as_ref().unwrap().0,
        "experiment.m27-depth-limit"
    );
    assert_eq!(
        expected.unresolved_box,
        Some(
            [
                ["7/8", "15/16"],
                ["7/16", "1/2"],
                ["0", "1/16"],
                ["3/8", "7/16"],
                ["53/16", "109/32"],
                ["4", "33/8"],
            ]
            .map(|b| b.map(str::to_owned))
        )
    );
    for workers in [2, 4] {
        std::thread::scope(|scope| {
            let jobs = (0..workers).map(|_| scope.spawn(run)).collect::<Vec<_>>();
            for job in jobs {
                assert_eq!(job.join().unwrap(), expected);
            }
        });
        println!("M27 workers={workers} identical independent requests");
    }
}
