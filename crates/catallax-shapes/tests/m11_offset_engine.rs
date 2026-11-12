//! M11's offset spine posed to M10's patch-pair engine. Opt-in research: no
//! production fillet calls this code.
//!
//! The engine meets two exact rational patches. An inward offset
//! `O = P - r N/|N|` is rational in the support's own parameters only when
//! `|N|^2` is a square, and on these supports it is not. On the two loft
//! edges the normal is linear in `(u, v)`, so `|N|^2` is quadratic in `u`
//! with a square leading coefficient `k^2`. The line `|N| = k u + tau` then
//! meets the quadric `w^2 = |N|^2` once more, which gives `u` and `|N|` as
//! exact rational functions of `(v, tau)`. With `tau` sheared along `v` and
//! scaled to `sigma` in `[0, 1]`, each offset becomes one exact rational
//! patch of degree 5 in `v` and 4 in `sigma`. The chart and its coverage
//! checks are symbolic, like M11's fixture construction, and are not charged.
//! The patch proof and any search are charged to one 20M request per edge.
//!
//! The biquadratic control has `|N|^2` of degree 4 in `u`, so this chart does
//! not exist for it. Each record names the step where posing stopped.

use catallax_numbers::{Decider, MPoly, Rational, Real, Refusal, with_operation_budget};
use catallax_shapes::{BSplineSurface, PatchPairIntersection, Point3, RationalPatch};
use std::num::NonZeroUsize;

type Q = Rational;
const ALLOWANCE: u64 = 20_000_000;
/// Polynomial variables: the support's first parameter, its second, and the
/// chart coordinate (`tau`, later `sigma`).
const X: usize = 0;
const Y: usize = 1;
const T: usize = 2;
const NV: usize = 3;

fn q(n: i64) -> Q {
    Q::from_i64(n)
}
fn r(text: &str) -> Q {
    Q::parse(text).unwrap()
}
fn k(n: Q) -> MPoly {
    MPoly::constant(NV, n)
}
fn var(i: usize) -> MPoly {
    MPoly::var(NV, i)
}
fn pow(p: &MPoly, n: usize) -> MPoly {
    (0..n).fold(k(q(1)), |a, _| &a * p)
}

/// `p` split by powers of `axis`: entry `i` multiplies `axis^i`.
fn by_power(p: &MPoly, axis: usize) -> Vec<MPoly> {
    let mut out: Vec<MPoly> = Vec::new();
    for (m, c) in p.terms() {
        let mut e = m.exponents().to_vec();
        let n = e[axis] as usize;
        e[axis] = 0;
        while out.len() <= n {
            out.push(MPoly::zero(NV));
        }
        out[n] = &out[n] + &MPoly::from_terms(NV, vec![(e, c)]);
    }
    out
}

/// `p(x = num/den) den^degree`, for `p` of degree at most `degree` in `x`.
fn homogenize(p: &MPoly, num: &MPoly, den: &MPoly, degree: usize) -> MPoly {
    let parts = by_power(p, X);
    assert!(parts.len() <= degree + 1);
    parts.iter().enumerate().fold(MPoly::zero(NV), |a, (i, f)| {
        &a + &(&(f * &pow(num, i)) * &pow(den, degree - i))
    })
}

fn binomial(n: usize, k: usize) -> Q {
    (0..k).fold(q(1), |out, i| {
        &(&out * &q((n - i) as i64)) / &q((i + 1) as i64)
    })
}

/// Degrees of `ps` in `(Y, T)`.
fn degrees(ps: &[&MPoly]) -> [usize; 2] {
    let mut d = [0, 0];
    for p in ps {
        for (m, _) in p.terms() {
            assert_eq!(m.exponents()[X], 0);
            d[0] = d[0].max(m.exponents()[Y] as usize);
            d[1] = d[1].max(m.exponents()[T] as usize);
        }
    }
    d
}

/// Bernstein coefficients of `p` in `(Y, T)` on `[0,1]^2`.
fn bernstein(p: &MPoly, degree: [usize; 2]) -> Vec<Vec<Q>> {
    let mut power = vec![vec![q(0); degree[1] + 1]; degree[0] + 1];
    for (m, c) in p.terms() {
        let e = m.exponents();
        power[e[Y] as usize][e[T] as usize] = c;
    }
    let mut out = vec![vec![q(0); degree[1] + 1]; degree[0] + 1];
    for (a, row) in out.iter_mut().enumerate() {
        for (b, entry) in row.iter_mut().enumerate() {
            for (i, line) in power.iter().enumerate().take(a + 1) {
                for (j, c) in line.iter().enumerate().take(b + 1) {
                    if c.is_zero() {
                        continue;
                    }
                    let f = &(&binomial(a, i) / &binomial(degree[0], i))
                        * &(&binomial(b, j) / &binomial(degree[1], j));
                    *entry = &*entry + &(&f * c);
                }
            }
        }
    }
    out
}

/// Whether `p` is positive on `[0,1]^2` in `(Y, T)`: every Bernstein
/// coefficient positive, after halving the domain up to 16 times.
fn positive(p: &MPoly) -> bool {
    fn halves(p: &MPoly, depth: u32) -> bool {
        let d = degrees(&[p]);
        let b = bernstein(p, d);
        if b.iter().flatten().all(|c| c > &q(0)) {
            return true;
        }
        let corners = [&b[0][0], &b[0][d[1]], &b[d[0]][0], &b[d[0]][d[1]]];
        if depth == 0 || corners.iter().any(|c| *c <= &q(0)) {
            return false;
        }
        let axis = if d[0] >= d[1] { Y } else { T };
        let half = r("1/2");
        [q(0), half.clone()].into_iter().all(|lo| {
            let map = &var(axis).scale(&half) + &k(lo);
            halves(&p.substitute(axis, &map), depth - 1)
        })
    }
    halves(p, 16)
}

fn cross(a: &[MPoly; 3], b: &[MPoly; 3]) -> [MPoly; 3] {
    std::array::from_fn(|i| {
        &(&a[(i + 1) % 3] * &b[(i + 2) % 3]) - &(&a[(i + 2) % 3] * &b[(i + 1) % 3])
    })
}

fn dot(a: &[MPoly; 3], b: &[MPoly; 3]) -> MPoly {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}

/// An exact rational chart of one inward offset: `tau = tau0 + slope y +
/// width sigma`, `x = (c(y) - tau^2) / (2 k tau - b(y))` and
/// `|N| = k x + tau`, where `|N|^2 = k^2 x^2 + b(y) x + c(y)`.
struct Chart {
    support: [MPoly; 3],
    normal: [MPoly; 3],
    k: Q,
    tau: [Q; 3],
    degree: [usize; 2],
    num: [MPoly; 3],
    den: MPoly,
}

fn no_chart(detail: String) -> Refusal {
    Refusal::new("experiment.m11-offset-no-rational-chart", detail)
}

impl Chart {
    fn new(support: &[MPoly; 3], radius: &Q) -> Result<Self, Refusal> {
        let px = support.each_ref().map(|f| f.partial(X));
        let py = support.each_ref().map(|f| f.partial(Y));
        let normal = cross(&px, &py);
        let nn = dot(&normal, &normal);
        let parts = by_power(&nn, X);
        if parts.len() != 3 || parts[2].total_degree() != Some(0) {
            return Err(no_chart(format!(
                "the squared normal length has degree {} in u; this chart needs degree 2 \
                 with a constant leading coefficient",
                parts.len() - 1
            )));
        }
        let a = parts[2].eval(&[q(0), q(0), q(0)]);
        let kk = (1..=64)
            .map(q)
            .find(|c| (c * c) == a)
            .ok_or_else(|| no_chart(format!("the leading coefficient {a} is not a square")))?;
        if support
            .iter()
            .chain(&normal)
            .any(|f| by_power(f, X).len() > 2)
        {
            return Err(no_chart(
                "the support or its normal is not linear in u".into(),
            ));
        }
        // Binary64 picks a sheared band round the image of the support domain,
        // `tau` from `|N(1, y)| - k` up to `|N(0, y)|`; exact checks follow.
        let image = |x: f64, y: f64| {
            let at = [
                Q::from_f64_exact(x).unwrap(),
                Q::from_f64_exact(y).unwrap(),
                q(0),
            ];
            nn.eval(&at).to_f64_nearest().sqrt() - kk.to_f64_nearest() * x
        };
        let samples: Vec<[f64; 3]> = (0..=1024)
            .map(|i| {
                let y = f64::from(i) / 1024.0;
                [y, image(1.0, y), image(0.0, y)]
            })
            .collect();
        let mut best: Option<(i64, i64, i64)> = None;
        for step in -128..=128i64 {
            let slope = step as f64 / 32.0;
            let lo = samples
                .iter()
                .map(|[y, l, _]| l - slope * y)
                .fold(f64::INFINITY, f64::min);
            let hi = samples
                .iter()
                .map(|[y, _, h]| h - slope * y)
                .fold(f64::NEG_INFINITY, f64::max);
            let base = ((lo - 1.0 / 32.0) * 32.0).floor() as i64;
            let width = ((hi + 1.0 / 32.0) * 32.0).ceil() as i64 - base;
            if best.is_none_or(|b| width < b.2) {
                best = Some((step, base, width));
            }
        }
        let (step, base, width) = best.unwrap();
        let slope = Q::from_ratio(step, 32).unwrap();
        let tau0 = Q::from_ratio(base, 32).unwrap();
        let width = Q::from_ratio(width, 32).unwrap();

        let tau = var(T);
        let xnum = &parts[0] - &(&tau * &tau);
        let xden = &tau.scale(&(&q(2) * &kk)) - &parts[1];
        let wnum = &xnum.scale(&kk) + &(&tau * &xden);
        let ph = support.each_ref().map(|f| homogenize(f, &xnum, &xden, 1));
        let nh = normal.each_ref().map(|f| homogenize(f, &xnum, &xden, 1));
        let lower = &k(tau0.clone()) + &var(Y).scale(&slope);
        let upper = &lower + &k(width.clone());
        let sigma = &lower + &var(T).scale(&width);
        // O = P - r N/|N| = (ph wnum - r xden nh) / (xden wnum).
        let num: [MPoly; 3] = std::array::from_fn(|i| {
            (&(&ph[i] * &wnum) - &(&xden * &nh[i]).scale(radius)).substitute(T, &sigma)
        });
        let den = (&xden * &wnum).substitute(T, &sigma);
        // `tau` falls strictly as `x` rises, from `|N(0, y)|` at `x = 0` to
        // `|N(1, y)| - k` at `x = 1`, once `|N|` never vanishes. So the band
        // holds the support domain when these are positive on `[0, 1]^2`.
        let shifted = &lower + &k(kk.clone());
        let checks = [
            (
                "the normal never vanishes",
                &parts[0].scale(&(&q(4) * &a)) - &(&parts[1] * &parts[1]),
            ),
            ("the upper line is positive", upper.clone()),
            (
                "u = 0 lies below the upper line",
                &(&upper * &upper) - &nn.substitute(X, &k(q(0))),
            ),
            ("the lower line plus k is positive", shifted.clone()),
            (
                "u = 1 lies above the lower line",
                &nn.substitute(X, &k(q(1))) - &(&shifted * &shifted),
            ),
            (
                "the chart denominator is positive",
                xden.substitute(T, &sigma),
            ),
            (
                "the chart takes the positive root",
                wnum.substitute(T, &sigma),
            ),
        ];
        for (name, p) in checks {
            if !positive(&p) {
                return Err(Refusal::new(
                    "experiment.m11-offset-chart-unproved",
                    format!("chart check failed: {name}"),
                ));
            }
        }
        let degree = degrees(&[&num[0], &num[1], &num[2], &den]);
        Ok(Self {
            support: support.clone(),
            normal,
            k: kk,
            tau: [tau0, slope, width],
            degree,
            num,
            den,
        })
    }

    /// The chart as one rational Bezier patch on `[0, 1]^2` in `(v, sigma)`.
    /// The patch proof is charged to `dec`.
    fn patch(&self, dec: &mut Decider) -> Result<RationalPatch, Refusal> {
        let degree = self.degree;
        let w = bernstein(&self.den, degree);
        let nb = self.num.each_ref().map(|f| bernstein(f, degree));
        let mut controls = Vec::new();
        let mut weights = Vec::new();
        for i in 0..=degree[0] {
            for j in 0..=degree[1] {
                let wij = &w[i][j];
                assert!(wij > &q(0), "chart weights are positive");
                controls.push(Point3::new(
                    (&nb[0][i][j] / wij).into(),
                    (&nb[1][i][j] / wij).into(),
                    (&nb[2][i][j] / wij).into(),
                ));
                weights.push(Real::from(wij.clone()));
            }
        }
        let knots = |d: usize| -> Vec<Q> { (0..2 * d + 2).map(|i| q(i64::from(i > d))).collect() };
        let carrier = BSplineSurface::new(
            degree,
            [knots(degree[0]), knots(degree[1])],
            [degree[0] + 1, degree[1] + 1],
            controls,
            weights,
            dec,
        )?;
        RationalPatch::new(carrier, dec)
    }

    fn summary(&self) -> String {
        let [tau0, slope, width] = &self.tau;
        format!(
            "k {}, tau = {tau0} + {slope} v + {width} sigma, degree {:?}",
            self.k, self.degree
        )
    }
}

/// M11's three support pairs, each support in `(X, Y)`; each offset is
/// charted on its own.
fn supports(index: usize) -> (&'static str, [MPoly; 3], [MPoly; 3]) {
    let one = k(q(1));
    let a = |u: &MPoly, v: &MPoly| -> [MPoly; 3] {
        [
            &(&k(q(-1)) + &u.scale(&q(2))) + &(v * &(&one - u)),
            &k(q(-1)) + &(u * v),
            v.scale(&q(2)),
        ]
    };
    let b = |u: &MPoly, v: &MPoly| -> [MPoly; 3] {
        [
            &one - &(u * v),
            &(&k(q(-1)) + &u.scale(&q(2))) + v,
            v.scale(&q(2)),
        ]
    };
    let c = |s: &MPoly, t: &MPoly| -> [MPoly; 3] {
        [
            &(&one - &s.scale(&q(2))) - t,
            &one + &(t * &(&one - &s.scale(&q(2)))),
            t.scale(&q(2)),
        ]
    };
    let (u, v) = (var(X), var(Y));
    match index {
        0 => ("m27-loft-edge", a(&u, &v), b(&u, &v)),
        1 => ("loft-next-edge", b(&u, &v), c(&u, &v)),
        _ => {
            let bump = (&(&(&u * &(&one - &u)) * &v) * &(&one - &v)).scale(&r("1/5"));
            let mut first = a(&u, &v);
            let mut second = b(&u, &v);
            first[2] = &first[2] + &bump;
            second[2] = &second[2] - &bump;
            ("biquadratic-edge", first, second)
        }
    }
}

/// What one edge's request proved, in a form that compares across thread
/// counts.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Record {
    name: &'static str,
    charts: Vec<String>,
    /// The offset (0 first, 1 second) whose posing refused, if one did.
    refused_offset: Option<usize>,
    answer: Result<String, (String, String)>,
    searched: bool,
    used: u64,
    exhausted: bool,
}

fn answer(found: &PatchPairIntersection) -> String {
    let closed = found.branches().iter().filter(|b| b.is_closed()).count();
    let s = found.search();
    let mut out = format!(
        "{} open, {closed} closed; {} arcs, {} cells, depth {}",
        found.branches().len() - closed,
        s.arcs,
        s.cells,
        s.depth
    );
    for branch in found.branches() {
        for arc in branch.arcs() {
            out.push_str(&format!(
                "; {:?} {:?} {:?} {:?} {:?}",
                arc.first_cell,
                arc.second_cell,
                arc.certificate,
                arc.ends[0].enclosure(),
                arc.ends[1].enclosure()
            ));
        }
    }
    out
}

fn run(index: usize) -> Record {
    let (name, first, second) = supports(index);
    let radius = r("1/10");
    let mut record = Record {
        name,
        charts: Vec::new(),
        refused_offset: None,
        answer: Ok(String::new()),
        searched: false,
        used: 0,
        exhausted: false,
    };
    let (result, work) = with_operation_budget(ALLOWANCE, || {
        let mut dec = Decider::default();
        let mut patches = Vec::new();
        for (side, support) in [&first, &second].into_iter().enumerate() {
            record.refused_offset = Some(side);
            let chart = Chart::new(support, &radius)?;
            record.charts.push(chart.summary());
            patches.push(chart.patch(&mut dec)?);
        }
        record.refused_offset = None;
        record.searched = true;
        patches[0].patch_meeting(&patches[1], &mut dec)
    });
    record.answer = match result {
        Ok(found) => Ok(answer(&found)),
        Err(refusal) => Err((refusal.code.to_string(), refusal.detail)),
    };
    record.used = work.used;
    record.exhausted = work.exhausted;
    record
}

/// Every chart point is an exact offset point: at rational `(v, sigma)` the
/// chart gives rational `u` and `w` with `w^2 = |N(u, v)|^2` exactly and
/// `w > 0`, and the patch point equals `P - r N / w` exactly.
#[test]
fn offset_charts_are_exact() {
    let radius = r("1/10");
    for index in 0..2 {
        let (name, first, second) = supports(index);
        for support in [first, second] {
            let chart = Chart::new(&support, &radius).unwrap_or_else(|e| panic!("{name}: {e:?}"));
            let [tau0, slope, width] = &chart.tau;
            let nn = dot(&chart.normal, &chart.normal);
            let parts = by_power(&nn, X);
            for (y, s) in [("0", "0"), ("1/3", "1/2"), ("1", "1"), ("2/7", "5/6")] {
                let (y, s) = (r(y), r(s));
                let tau = &(tau0 + &(slope * &y)) + &(width * &s);
                let at = [q(0), y.clone(), q(0)];
                let (b, c) = (parts[1].eval(&at), parts[0].eval(&at));
                let x = &(&c - &(&tau * &tau)) / &(&(&(&q(2) * &chart.k) * &tau) - &b);
                let w = &(&chart.k * &x) + &tau;
                let point = [x.clone(), y.clone(), q(0)];
                assert!(w > q(0));
                assert_eq!(&w * &w, nn.eval(&point));
                let chart_at = [q(0), y.clone(), s.clone()];
                let den = chart.den.eval(&chart_at);
                for i in 0..3 {
                    let offset = &chart.support[i].eval(&point)
                        - &(&(&radius * &chart.normal[i].eval(&point)) / &w);
                    assert_eq!(&chart.num[i].eval(&chart_at) / &den, offset, "{name}");
                }
            }
        }
    }
    let (_, first, _) = supports(2);
    let refusal = Chart::new(&first, &radius).err().expect("no chart");
    assert_eq!(refusal.code, "experiment.m11-offset-no-rational-chart");
}

/// The three M11 edges posed to the M10 engine under the unchanged 20M
/// allowance, with the same records at 1, 2 and 4 workers.
#[test]
#[ignore = "M11 opt-in: posing offset spines to the M10 engine"]
fn offset_spines_posed_to_the_pair_engine() {
    catallax_numbers::set_run_ahead_wait(true);
    let mut seen: Option<Vec<Record>> = None;
    for threads in [1, 2, 4] {
        catallax_numbers::set_threads(NonZeroUsize::new(threads));
        let records: Vec<Record> = (0..3).map(run).collect();
        for record in &records {
            println!("M11 workers={threads} {record:?}");
            assert!(record.used <= ALLOWANCE);
        }
        match &seen {
            Some(first) => assert_eq!(first, &records, "{threads} workers"),
            None => seen = Some(records),
        }
    }
    catallax_numbers::set_threads(None);
    catallax_numbers::set_run_ahead_wait(false);
    // Published: the first offset's patch proof alone spends the allowance on
    // both loft edges, and the biquadratic control has no chart. No search
    // starts.
    let published = [
        ("numbers.operation-budget", 20_000_000, true, 1),
        ("numbers.operation-budget", 19_999_978, true, 1),
        ("experiment.m11-offset-no-rational-chart", 0, false, 0),
    ];
    for (record, (code, used, exhausted, charts)) in seen.unwrap().iter().zip(published) {
        let refusal = record.answer.as_ref().expect_err("no edge is answered");
        assert_eq!(
            (
                refusal.0.as_str(),
                record.used,
                record.exhausted,
                record.charts.len()
            ),
            (code, used, exhausted, charts),
            "{}",
            record.name
        );
        assert_eq!(record.refused_offset, Some(0));
        assert!(!record.searched);
    }
}
