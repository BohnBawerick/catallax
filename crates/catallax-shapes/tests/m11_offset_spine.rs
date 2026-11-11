//! Opt-in offset-spine research. No production fillet entry calls this code.
use catallax_numbers::{Decider, MPoly, Rational, Refusal, with_operation_budget};
use std::collections::BTreeMap;

const N: usize = 4;
const ALLOWANCE: u64 = 20_000_000;
const MAX_DEPTH: u32 = 24;
const MAX_CELLS: usize = 16_384;
type I = [Rational; 2];
type Bounds = [I; N];
type V = [I; 3];
type Result<T> = std::result::Result<T, Refusal>;

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn point(x: Rational) -> I {
    [x.clone(), x]
}
fn ar(a: &Rational, b: &Rational, op: char, d: &mut Decider) -> Result<Rational> {
    d.reserve_exact_arithmetic(a, b, op == '*' || op == '/')?;
    Ok(match op {
        '+' => a + b,
        '-' => a - b,
        '*' => a * b,
        '/' => a / b,
        _ => unreachable!(),
    })
}
fn mid(x: &I, d: &mut Decider) -> Result<Rational> {
    ar(&ar(&x[0], &x[1], '+', d)?, &q(2), '/', d)
}
fn add(a: &I, b: &I, d: &mut Decider) -> Result<I> {
    Ok([ar(&a[0], &b[0], '+', d)?, ar(&a[1], &b[1], '+', d)?])
}
fn neg(a: &I) -> I {
    [-&a[1], -&a[0]]
}
fn sub(a: &I, b: &I, d: &mut Decider) -> Result<I> {
    add(a, &neg(b), d)
}
fn mul(a: &I, b: &I, d: &mut Decider) -> Result<I> {
    let mut values = Vec::new();
    for x in a {
        for y in b {
            values.push(ar(x, y, '*', d)?);
        }
    }
    d.charge_work(8)?;
    Ok([
        values.iter().min().unwrap().clone(),
        values.iter().max().unwrap().clone(),
    ])
}
fn square(a: &I, d: &mut Decider) -> Result<I> {
    let x = ar(&a[0], &a[0], '*', d)?;
    let y = ar(&a[1], &a[1], '*', d)?;
    d.charge_work(3)?;
    Ok([
        if a[0] <= q(0) && a[1] >= q(0) {
            q(0)
        } else {
            x.clone().min(y.clone())
        },
        x.max(y),
    ])
}
fn div(a: &I, b: &I, d: &mut Decider) -> Result<I> {
    d.charge_work(1)?;
    assert!(b[0] > q(0));
    mul(
        a,
        &[ar(&q(1), &b[1], '/', d)?, ar(&q(1), &b[0], '/', d)?],
        d,
    )
}
fn sqrt(a: &I, d: &mut Decider) -> Result<I> {
    d.charge_work(1)?;
    assert!(a[0] >= q(0));
    let mut hi = q(1);
    while ar(&hi, &hi, '*', d)? < a[1] {
        hi = ar(&hi, &q(2), '*', d)?;
    }
    let mut answer = [q(0), hi.clone()];
    for side in 0..2 {
        if a[side].is_zero() {
            answer[side] = q(0);
            continue;
        }
        if ar(&hi, &hi, '*', d)? == a[side] {
            answer[side] = hi.clone();
            continue;
        }
        let mut bracket = [q(0), hi.clone()];
        for _ in 0..32 {
            let m = mid(&bracket, d)?;
            let sq = ar(&m, &m, '*', d)?;
            d.charge_work(1)?;
            if sq == a[side] {
                bracket = point(m);
                break;
            }
            if sq < a[side] {
                bracket[0] = m;
            } else {
                bracket[1] = m;
            }
        }
        answer[side] = bracket[side].clone();
    }
    Ok(answer)
}
fn dot(a: &V, b: &V, d: &mut Decider) -> Result<I> {
    let mut sum = point(q(0));
    for i in 0..3 {
        sum = add(&sum, &mul(&a[i], &b[i], d)?, d)?;
    }
    Ok(sum)
}
fn cross(a: &V, b: &V, d: &mut Decider) -> Result<V> {
    let mut out = std::array::from_fn(|_| point(q(0)));
    for i in 0..3 {
        out[i] = sub(
            &mul(&a[(i + 1) % 3], &b[(i + 2) % 3], d)?,
            &mul(&a[(i + 2) % 3], &b[(i + 1) % 3], d)?,
            d,
        )?;
    }
    Ok(out)
}
fn eval(p: &MPoly, b: &Bounds, d: &mut Decider) -> Result<I> {
    let mut sum = point(q(0));
    for (m, coefficient) in p.terms() {
        d.charge_work(1)?;
        let mut term = point(coefficient.clone());
        for (axis, bound) in b.iter().enumerate() {
            for _ in 0..m.exponents()[axis] {
                term = mul(&term, bound, d)?;
            }
        }
        sum = add(&sum, &term, d)?;
    }
    Ok(sum)
}
fn vector(p: &[MPoly; 3], b: &Bounds, d: &mut Decider) -> Result<V> {
    Ok([eval(&p[0], b, d)?, eval(&p[1], b, d)?, eval(&p[2], b, d)?])
}

struct Patch {
    p: [MPoly; 3],
    dp: [[MPoly; 3]; 2],
    normal: [MPoly; 3],
    dn: [[MPoly; 3]; 2],
}
impl Patch {
    fn new(p: [MPoly; 3], axes: [usize; 2]) -> Self {
        let dp = axes.map(|axis| p.each_ref().map(|x| x.partial(axis)));
        let normal: [MPoly; 3] = std::array::from_fn(|i| {
            &(&dp[0][(i + 1) % 3] * &dp[1][(i + 2) % 3])
                - &(&dp[0][(i + 2) % 3] * &dp[1][(i + 1) % 3])
        });
        let dn = axes.map(|axis| normal.each_ref().map(|x| x.partial(axis)));
        Self { p, dp, normal, dn }
    }
    fn offset(&self, b: &Bounds, d: &mut Decider) -> Result<(V, [V; 2])> {
        let p = vector(&self.p, b, d)?;
        let n = vector(&self.normal, b, d)?;
        let mut norm2 = point(q(0));
        for x in &n {
            norm2 = add(&norm2, &square(x, d)?, d)?;
        }
        if norm2[0] <= q(0) {
            return Err(Refusal::new(
                "experiment.m11-normal-range-unproved",
                "normal interval contains zero; this is not a singularity certificate",
            ));
        }
        let length = sqrt(&norm2, d)?;
        if length[0] <= q(0) {
            return Err(Refusal::new(
                "experiment.m11-normal-range-unproved",
                "positive normal length is unresolved at the enclosure precision",
            ));
        }
        let length3 = mul(&length, &norm2, d)?;
        let radius = point(Rational::from_ratio(1, 10).unwrap());
        let mut offset = p;
        for i in 0..3 {
            offset[i] = sub(&offset[i], &mul(&radius, &div(&n[i], &length, d)?, d)?, d)?;
        }
        let mut derivatives: [V; 2] = std::array::from_fn(|_| std::array::from_fn(|_| point(q(0))));
        for (axis, derivative) in derivatives.iter_mut().enumerate() {
            let dp = vector(&self.dp[axis], b, d)?;
            let dn = vector(&self.dn[axis], b, d)?;
            let ndn = dot(&n, &dn, d)?;
            for i in 0..3 {
                let correction = sub(
                    &div(&dn[i], &length, d)?,
                    &div(&mul(&n[i], &ndn, d)?, &length3, d)?,
                    d,
                )?;
                derivative[i] = sub(&dp[i], &mul(&radius, &correction, d)?, d)?;
            }
        }
        Ok((offset, derivatives))
    }
}
struct Case {
    name: &'static str,
    a: Patch,
    b: Patch,
}
fn fixture(index: usize) -> Case {
    let [u, v, s, t] = std::array::from_fn(|i| MPoly::var(N, i));
    let c = |n| MPoly::constant(N, q(n));
    let a = |u: &MPoly, v: &MPoly| {
        [
            &(&c(-1) + &u.scale(&q(2))) + &(v * &(&c(1) - u)),
            &c(-1) + &(u * v),
            v.scale(&q(2)),
        ]
    };
    let b = |u: &MPoly, v: &MPoly| {
        [
            &c(1) - &(u * v),
            &(&c(-1) + &u.scale(&q(2))) + v,
            v.scale(&q(2)),
        ]
    };
    let (name, mut first, mut second) = if index == 1 {
        (
            "loft-next-edge",
            b(&u, &v),
            [
                &(&c(1) - &s.scale(&q(2))) - &t,
                &c(1) + &(&t * &(&c(1) - &s.scale(&q(2)))),
                t.scale(&q(2)),
            ],
        )
    } else {
        (
            if index == 0 {
                "m27-loft-edge"
            } else {
                "biquadratic-edge"
            },
            a(&u, &v),
            b(&s, &t),
        )
    };
    if index == 2 {
        let bump = |u: &MPoly, v: &MPoly| {
            (&(&(u * &(&c(1) - u)) * v) * &(&c(1) - v)).scale(&(q(1) / q(5)))
        };
        first[2] = &first[2] + &bump(&u, &v);
        second[2] = &second[2] - &bump(&s, &t);
    }
    Case {
        name,
        a: Patch::new(first, [0, 1]),
        b: Patch::new(second, [2, 3]),
    }
}
fn offsets(case: &Case, bounds: &Bounds, d: &mut Decider) -> Result<(V, [V; N])> {
    let (a, da) = case.a.offset(bounds, d)?;
    let (b, db) = case.b.offset(bounds, d)?;
    Ok((
        [
            sub(&a[0], &b[0], d)?,
            sub(&a[1], &b[1], d)?,
            sub(&a[2], &b[2], d)?,
        ],
        [
            da[0].clone(),
            da[1].clone(),
            db[0].each_ref().map(neg),
            db[1].each_ref().map(neg),
        ],
    ))
}
fn inverse(mut a: Vec<Vec<Rational>>, d: &mut Decider) -> Result<Option<Vec<Vec<Rational>>>> {
    let mut out = (0..3)
        .map(|i| (0..3).map(|j| q(i64::from(i == j))).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for k in 0..3 {
        d.charge_work(3)?;
        let Some(pivot) = (k..3).find(|&i| !a[i][k].is_zero()) else {
            return Ok(None);
        };
        a.swap(k, pivot);
        out.swap(k, pivot);
        let factor = a[k][k].clone();
        for j in 0..3 {
            a[k][j] = ar(&a[k][j], &factor, '/', d)?;
            out[k][j] = ar(&out[k][j], &factor, '/', d)?;
        }
        for i in 0..3 {
            if i == k {
                continue;
            }
            let factor = a[i][k].clone();
            for j in 0..3 {
                a[i][j] = ar(&a[i][j], &ar(&factor, &a[k][j], '*', d)?, '-', d)?;
                out[i][j] = ar(&out[i][j], &ar(&factor, &out[k][j], '*', d)?, '-', d)?;
            }
        }
    }
    Ok(Some(out))
}
fn graph(case: &Case, bounds: &Bounds, jac: &[V; N], d: &mut Decider) -> Result<bool> {
    for free in (0..N).rev() {
        let solved = (0..N).filter(|&axis| axis != free).collect::<Vec<_>>();
        let mut j = vec![Vec::new(); 3];
        for &axis in &solved {
            for (row, interval) in j.iter_mut().zip(&jac[axis]) {
                row.push(mid(interval, d)?);
            }
        }
        let Some(p) = inverse(j, d)? else { continue };
        let mut centered = bounds.clone();
        for &axis in &solved {
            centered[axis] = point(mid(&bounds[axis], d)?);
        }
        let f = offsets(case, &centered, d)?.0;
        let mut certified = true;
        for (i, &axis) in solved.iter().enumerate() {
            let mut image = point(mid(&bounds[axis], d)?);
            for k in 0..3 {
                image = sub(&image, &mul(&point(p[i][k].clone()), &f[k], d)?, d)?;
            }
            let mut norm = q(0);
            for (j, &other) in solved.iter().enumerate() {
                let mut r = point(q(i64::from(i == j)));
                for k in 0..3 {
                    r = sub(&r, &mul(&point(p[i][k].clone()), &jac[other][k], d)?, d)?;
                }
                let magnitude = r[0].abs().max(r[1].abs());
                norm = ar(&norm, &magnitude, '+', d)?;
                let delta = sub(&bounds[other], &point(mid(&bounds[other], d)?), d)?;
                image = add(&image, &mul(&r, &delta, d)?, d)?;
            }
            d.charge_work(3)?;
            if norm >= q(1) || image[0] <= bounds[axis][0] || image[1] >= bounds[axis][1] {
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
struct Grid {
    positive: usize,
    negative: usize,
    unresolved: usize,
    samples: usize,
    ambiguous_nearest: usize,
}
fn quadrant(f: &[I; 2]) -> Option<i32> {
    if f[0][0] >= q(0) && f[1][0] > q(0) {
        Some(0)
    } else if f[0][1] < q(0) && f[1][0] >= q(0) {
        Some(1)
    } else if f[0][1] <= q(0) && f[1][1] < q(0) {
        Some(2)
    } else if f[0][0] > q(0) && f[1][1] <= q(0) {
        Some(3)
    } else {
        None
    }
}
fn grid_index(quadrants: &[Option<i32>]) -> Option<i32> {
    let mut sum = 0;
    for i in 0..4 {
        let (a, b) = (quadrants[i]?, quadrants[(i + 1) % 4]?);
        let delta = (b - a + 4) % 4;
        if delta == 2 {
            return None;
        }
        sum += if delta == 3 { -1 } else { delta };
    }
    assert_eq!(sum % 4, 0);
    Some(sum / 4)
}
// The nearest grid point and its field winding propose starts only.
fn candidates(case: &Case, d: &mut Decider) -> Result<Grid> {
    let mut target = Vec::new();
    for s in 0..=2 {
        for t in 0..=2 {
            let bounds = [
                point(q(0)),
                point(q(0)),
                point(q(s) / q(2)),
                point(q(t) / q(2)),
            ];
            let (p, dp) = case.b.offset(&bounds, d)?;
            target.push((p, cross(&dp[0], &dp[1], d)?));
        }
    }
    let mut cache = BTreeMap::new();
    let mut record = Grid::default();
    let mut pending = vec![([q(0), q(1)], [q(0), q(1)], 0)];
    while let Some((u, v, depth)) = pending.pop() {
        let corners = [
            (u[0].clone(), v[0].clone()),
            (u[1].clone(), v[0].clone()),
            (u[1].clone(), v[1].clone()),
            (u[0].clone(), v[1].clone()),
        ];
        let mut quadrants = Vec::new();
        for key in corners {
            if !cache.contains_key(&key) {
                let b = [
                    point(key.0.clone()),
                    point(key.1.clone()),
                    point(q(0)),
                    point(q(0)),
                ];
                let (p, dp) = case.a.offset(&b, d)?;
                let mut distances = Vec::new();
                for (other, _) in &target {
                    let mut sum = point(q(0));
                    for i in 0..3 {
                        sum = add(&sum, &square(&sub(&p[i], &other[i], d)?, d)?, d)?;
                    }
                    distances.push(sum);
                }
                let mut best = 0;
                for k in 1..target.len() {
                    if mid(&distances[k], d)? < mid(&distances[best], d)? {
                        best = k;
                    }
                }
                if distances
                    .iter()
                    .enumerate()
                    .any(|(k, x)| k != best && x[0] <= distances[best][1])
                {
                    record.ambiguous_nearest += 1;
                }
                let f = [
                    dot(&dp[0], &target[best].1, d)?,
                    dot(&dp[1], &target[best].1, d)?,
                ];
                cache.insert(key.clone(), quadrant(&f));
            }
            quadrants.push(cache[&key]);
        }
        d.charge_work(4)?;
        if let Some(index) = grid_index(&quadrants) {
            if index > 0 {
                record.positive += 1;
            } else if index < 0 {
                record.negative += 1;
            }
        } else if depth == 2 {
            record.unresolved += 1;
        } else {
            let um = mid(&u, d)?;
            let vm = mid(&v, d)?;
            for ub in [[u[0].clone(), um.clone()], [um, u[1].clone()]] {
                for vb in [[v[0].clone(), vm.clone()], [vm.clone(), v[1].clone()]] {
                    pending.push((ub.clone(), vb, depth + 1));
                }
            }
        }
    }
    record.samples = cache.len();
    Ok(record)
}
#[derive(Debug, Default, PartialEq, Eq)]
struct Record {
    grid: Grid,
    grid_work: u64,
    visited: usize,
    empty: usize,
    graphs: usize,
    pending: usize,
    depth: u32,
    used: u64,
    exhausted: bool,
    refusal: Option<(String, String)>,
    unresolved_box: Option<[[String; 2]; N]>,
    offset_gap_l1_nm: Option<String>,
    whole_face_deviation_nm: Option<String>,
}
fn run(index: usize, allowance: u64) -> Record {
    let case = fixture(index);
    let mut r = Record::default();
    let (result, work) = with_operation_budget(allowance, || -> Result<()> {
        let mut d = Decider::default();
        r.grid = candidates(&case, &mut d)?;
        r.grid_work = catallax_numbers::current_operation_work().used;
        let mut pending = vec![(std::array::from_fn(|_| [q(0), q(1)]), 0)];
        while let Some((bounds, depth)) = pending.pop() {
            d.charge_work(1)?;
            r.visited += 1;
            r.pending = pending.len() + 1;
            r.depth = r.depth.max(depth);
            r.unresolved_box = Some(
                bounds
                    .each_ref()
                    .map(|x| x.each_ref().map(ToString::to_string)),
            );
            r.offset_gap_l1_nm = None;
            let (f, jac) = offsets(&case, &bounds, &mut d)?;
            let mut gap = q(0);
            for x in &f {
                gap = ar(&gap, &x[0].abs().max(x[1].abs()), '+', &mut d)?;
            }
            r.offset_gap_l1_nm = Some(ar(&gap, &q(1_000_000), '*', &mut d)?.to_string());
            d.charge_work(3)?;
            if f.iter().any(|x| x[0] > q(0) || x[1] < q(0)) {
                r.empty += 1;
                continue;
            }
            if graph(&case, &bounds, &jac, &mut d)? {
                r.graphs += 1;
                continue;
            }
            if depth == MAX_DEPTH {
                return Err(Refusal::new(
                    "experiment.m11-offset-depth-limit",
                    "offset equality remains unresolved at the unchanged depth cap",
                ));
            }
            if r.empty + r.graphs + pending.len() + 2 > MAX_CELLS {
                return Err(Refusal::new(
                    "experiment.m11-offset-cell-limit",
                    "offset cover exceeds the unchanged proof cell cap",
                ));
            }
            let mut axis = 0;
            let mut width = ar(&bounds[0][1], &bounds[0][0], '-', &mut d)?;
            for (i, b) in bounds.iter().enumerate().skip(1) {
                let w = ar(&b[1], &b[0], '-', &mut d)?;
                if w > width {
                    axis = i;
                    width = w;
                }
            }
            let midpoint = mid(&bounds[axis], &mut d)?;
            let mut left = bounds.clone();
            let mut right = bounds;
            left[axis][1] = midpoint.clone();
            right[axis][0] = midpoint;
            pending.push((right, depth + 1));
            pending.push((left, depth + 1));
        }
        r.pending = 0;
        r.unresolved_box = None;
        if r.graphs > 0 {
            return Err(Refusal::new(
                "experiment.m11-offset-joins-unproved",
                "regular graph boxes need exact endpoints, domain clipping and joins",
            ));
        }
        Ok(())
    });
    r.used = work.used;
    r.exhausted = work.exhausted;
    r.refusal = result.err().map(|e| (e.code.to_string(), e.detail));
    r
}

#[test]
#[ignore = "M11 opt-in proof-route experiment; no production fillet changes"]
fn offset_distance_spine_scope() {
    // All three input pairs share an exact edge, including the curved control.
    for index in 0..3 {
        let case = fixture(index);
        for h in [q(0), q(1) / q(3), q(1)] {
            let uv = [q(1), h.clone(), q(0), h];
            for i in 0..3 {
                assert_eq!(case.a.p[i].eval(&uv), case.b.p[i].eval(&uv));
            }
        }
    }
    let expected = (0..3).map(|i| run(i, ALLOWANCE)).collect::<Vec<_>>();
    let published = [
        (5_431_053, 195, 91, 13, 47_321, 148_170),
        (4_518_992, 159, 71, 17, 48_939, 147_813),
        (10_168_139, 209, 98, 13, 78_593, 151_770),
    ];
    for (i, r) in expected.iter().enumerate() {
        println!("M11 workers=1 case={} {r:?}", fixture(i).name);
        assert!(r.used <= ALLOWANCE);
        assert!(r.refusal.is_some());
        assert!(r.whole_face_deviation_nm.is_none());
        let (used, visited, empty, pending, grid_work, gap_ceiling) = published[i];
        assert_eq!(
            (r.used, r.visited, r.empty, r.pending, r.grid_work),
            (used, visited, empty, pending, grid_work)
        );
        assert_eq!((r.graphs, r.depth, r.exhausted), (0, MAX_DEPTH, false));
        assert_eq!(
            r.refusal.as_ref().unwrap().0,
            "experiment.m11-offset-depth-limit"
        );
        assert_eq!(
            r.grid,
            Grid {
                samples: 4,
                ..Grid::default()
            }
        );
        let gap = Rational::parse(r.offset_gap_l1_nm.as_ref().unwrap()).unwrap();
        assert!(gap <= q(gap_ceiling) && gap > q(gap_ceiling - 1));
    }
    for workers in [2, 4] {
        let rows = std::thread::scope(|scope| {
            let jobs = (0..workers)
                .map(|worker| {
                    scope.spawn(move || {
                        (worker..3)
                            .step_by(workers)
                            .map(|i| (i, run(i, ALLOWANCE)))
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Vec<_>>();
            let mut rows = jobs
                .into_iter()
                .flat_map(|job| job.join().unwrap())
                .collect::<Vec<_>>();
            rows.sort_by_key(|(i, _)| *i);
            rows.into_iter().map(|(_, r)| r).collect::<Vec<_>>()
        });
        assert_eq!(rows, expected);
        println!("M11 workers={workers} identical ordered answers, work and evidence");
    }
}

#[test]
#[ignore = "M11 opt-in proof controls"]
fn interval_controls() {
    assert_eq!(grid_index(&[Some(0), Some(1), Some(2), Some(3)]), Some(1));
    assert_eq!(grid_index(&[Some(3), Some(2), Some(1), Some(0)]), Some(-1));
    assert_eq!(grid_index(&[Some(0), Some(0), Some(3), Some(3)]), Some(0));
    assert_eq!(grid_index(&[Some(0), Some(2), Some(2), Some(0)]), None);
    let (result, work) = with_operation_budget(ALLOWANCE, || -> Result<()> {
        let mut d = Decider::default();
        for n in [0, 1, 2, 3, 16] {
            let x = point(q(n));
            let root = sqrt(&x, &mut d)?;
            assert!(ar(&root[0], &root[0], '*', &mut d)? <= q(n));
            assert!(ar(&root[1], &root[1], '*', &mut d)? >= q(n));
        }
        assert_eq!(square(&[q(-2), q(3)], &mut d)?, [q(0), q(9)]);
        let case = fixture(0);
        let b = std::array::from_fn(|_| point(q(1) / q(2)));
        let (p, dp) = case.a.offset(&b, &mut d)?;
        assert!(p.iter().all(|x| x[0] <= x[1]));
        assert!(dp.iter().flatten().all(|x| x[0] <= x[1]));
        assert_eq!(case.a.normal[0].eval(&[q(0), q(1), q(0), q(0)]), q(2));
        let origin = std::array::from_fn(|_| point(q(0)));
        let (p, dp) = case.a.offset(&origin, &mut d)?;
        assert_eq!(p, [point(q(-1)), point(q(-9) / q(10)), point(q(0))]);
        assert_eq!(dp[0], [point(q(2)), point(q(0)), point(q(-1) / q(20))]);
        assert_eq!(
            dp[1],
            [point(q(19) / q(20)), point(q(0)), point(q(81) / q(40))]
        );
        let [u, v, s, t] = std::array::from_fn(|i| MPoly::var(N, i));
        let zero = MPoly::zero(N);
        let planes = Case {
            name: "regular-control",
            a: Patch::new([u.clone(), v.clone(), zero.clone()], [0, 1]),
            b: Patch::new([s.clone(), zero.clone(), t], [2, 3]),
        };
        let bounds = [[q(0), q(1)], [q(-1), q(1)], [q(-1), q(2)], [q(-1), q(1)]];
        let (_, jac) = offsets(&planes, &bounds, &mut d)?;
        assert!(graph(&planes, &bounds, &jac, &mut d)?);
        let coincident = Case {
            name: "rank-deficient-control",
            a: planes.a,
            b: Patch::new([s, MPoly::var(N, 3), zero], [2, 3]),
        };
        let (_, jac) = offsets(&coincident, &bounds, &mut d)?;
        assert!(!graph(&coincident, &bounds, &jac, &mut d)?);
        Ok(())
    });
    result.unwrap();
    assert!(!work.exhausted);
    let tiny = run(0, 1);
    assert!(tiny.exhausted);
    assert_eq!(tiny.refusal.unwrap().0, "numbers.operation-budget");
    assert!(tiny.whole_face_deviation_nm.is_none());
}
