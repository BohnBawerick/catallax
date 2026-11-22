//! The field's own proof of a rolling ball's centre path, the spine, between
//! two whole rational patches: step 7.3 of 4.2a ([Checkpoint 2 decisions,
//! section 4.6], and the milestone note `docs/M42A-FIELD-FILLETS.md`).
//!
//! A ball of radius `r` touches support `A` at `S_A = X_A / W_A` and support
//! `B` at `S_B`, each from a chosen side `s` of the support's natural normal.
//! With `T_u = W X_u - W_u X` and `T_v = W X_v - W_v X` (each `W^2` times a
//! partial of `S`), the normal's numerator `M = T_u x T_v` is a positive
//! multiple of `S_u x S_v`. With an unknown length `L = |M|`, the centre `c`
//! satisfies, for each support,
//!
//! `L (W c - X) - s r W M = 0` and `L^2 - M . M = 0`,
//!
//! eight polynomial equations in nine values: the centre, the four contact
//! parameters and the two lengths. One parameter of `A`, `tau`, runs along the
//! spine and the other eight are unknown; `L > 0` picks the side. The B-rep's
//! fit (`catallax-shapes`, `spline_fillet`) eliminates the centre and the
//! lengths and solves for the three other contact parameters, so the two
//! proofs share no equations. Carrying `L` as an unknown is how the M10 offset
//! search keeps an offset polynomial (`docs/M11-OFFSET-PROOF.md`).
//!
//! The proof has three stages.
//!
//! 1. Trace, in binary64, with no authority: a grid on the three free contact
//!    parameters at the seed value of `tau`, Newton's method, then steps of
//!    `2^-6` in `tau` each way until a contact parameter leaves its support's
//!    square or `tau` reaches its own border. Those are the two ends.
//! 2. Ends: the system with the leaving parameter pinned at its border and
//!    `tau` unknown has one root, proved by the Krawczyk operator (Krawczyk
//!    1969; Moore, Kearfott and Cloud, "Introduction to Interval Analysis",
//!    2009, section 8.2) on a box round the traced end.
//! 3. Spans: the range of `tau` between the ends is cut into dyadic spans,
//!    eighths first. On each span the solution is interpolated at the six
//!    Chebyshev points by polynomials `p(s)` of degree five in the span's own
//!    variable `s` in `[-1, 1]` (Trefethen, "Approximation Theory and
//!    Approximation Practice", 2013, chapter 8: the error of an analytic
//!    function's interpolant falls geometrically with the degree). The
//!    unknowns are written `p(s) + y`, and the Krawczyk operator over the whole
//!    span, `K = -C G(0) + (I - C J) Y`, with `G(0)` the residual `F(p(s), s)`
//!    and `J` the Jacobian over every `s` and every `y` in `Y`, inside `Y`
//!    proves one root `y` in `Y` for each `tau` of the span (the parametric
//!    form: Neumaier, "Interval Methods for Systems of Equations", 1990,
//!    section 5.1; Martin, Goldsztejn, Granvilliers and Jermann, "Certified
//!    parallelotope continuation for one-manifolds", SINUM 51(6), 2013, the
//!    route the research note names). The residual and the Jacobian are
//!    computed on polynomials in `s` with interval coefficients, truncated
//!    above degree ten with the dropped terms bounded on `[-1, 1]` (a Taylor
//!    model in the sense of Makino and Berz, "Taylor models and other
//!    validated functional inclusion methods", IJPAM 4(4), 2003), so the
//!    residual is bounded by its fit error and not by the span's width. Each
//!    span's tube is the model curve of the centre widened by the centre part
//!    of `K`: its width `eps` is the largest such widening. A span whose proof
//!    fails or whose tube is wider than the target is halved, five times at
//!    most.
//!
//! Adjacent spans share a root at their common parameter. If a span-wide
//! enclosure is too wide for the next span's uniqueness box, a separate
//! point proof isolates the root there.
//! Each end root lies in its span's box, each contact parameter stays in its
//! support's square between the ends (a parameter that crosses its border at
//! an end is proved monotone there, from an enclosure of the spine's slope
//! `dx/dtau = -J_x^-1 F_tau`), and both lengths stay positive.
//!
//! OpenCASCADE (`BlendFunc_ConstRad` with `BRepBlend_Walking`: four surface
//! parameters on a section plane normal to a guide curve, solved by
//! `math_FunctionSetRoot`, 30 Newton steps, and marched with a deflection
//! test) and SMLib (`SmConstantRadiusFS`: equal offset-surface points in four
//! parameters, Newton steps marched by `SmFilletIntersector`, which checks each
//! span's Hermite curves against samples) trace the spine with Newton steps.
//! The field uses those steps to propose the spine, with
//! `math_FunctionSetRoot`'s 30 steps, and proves it with the Krawczyk operator.
//!
//! Polynomial evaluation, linear algebra and Newton proposals charge one
//! unit per binary64 or interval product, the limb rule's price for one-limb
//! operands (`exact_cost` in `catallax-numbers`). The proof is sequential.
//!
//! [Checkpoint 2 decisions, section 4.6]: ../../../docs/checkpoints/checkpoint-2-decisions.md

use std::cell::Cell;
use std::f64::consts::FRAC_1_SQRT_2;

use catallax_numbers::{Decider, Interval, Rational, Real, Refusal, Sign};
use catallax_shapes::{Pose, RationalPatch, Rotation, Surface};

pub mod sweep;

/// Degree of each span's model of the spine in the span's own variable.
pub const MODEL_DEGREE: usize = 5;
const TERMS: usize = MODEL_DEGREE + 1;
/// Highest power of `s` a span's polynomials keep, twice the model's degree,
/// so the product of two models is kept whole; higher terms are bounded on
/// `[-1, 1]` and folded into the constant term.
const ORDER: usize = 2 * MODEL_DEGREE;
/// The decisions' 100 nm cap in millimetres, rounded down.
pub const TUBE_CAP: f64 = 1e-4_f64.next_down();
/// Binary64 values of cos((2k+1)π/12), the roots of the degree-six
/// Chebyshev polynomial. Store the nodes so the proposal uses the same
/// numbers on every platform.
const NODES: [f64; TERMS] = [
    0.9659258262890683,
    FRAC_1_SQRT_2,
    0.25881904510252074,
    -0.25881904510252074,
    -FRAC_1_SQRT_2,
    -0.9659258262890683,
];
/// First spans are eighths of `tau`, as the B-rep fit starts from eight spans
/// (`docs/M11-OFFSET-PROOF.md`).
const FIRST_LEVEL: u32 = 3;
/// Spans are halved at most five times, the B-rep fit's limit, down to
/// `2^-8` of `tau`.
const LAST_LEVEL: u32 = 8;
/// The trace steps `2^-6` in `tau`, eight steps per first span.
const TRACE_STEP: f64 = 1.0 / 64.0;
/// The widest Krawczyk box in any unknown: `2^-14` mm, the largest power of
/// two below the 100 nm cap, so a proved centre lies inside the cap even
/// before the operator narrows it.
const BOX: f64 = 1.0 / 16_384.0;
/// Newton steps per solve, `math_FunctionSetRoot`'s count in OpenCASCADE's
/// `BRepBlend_Walking`.
const NEWTON_STEPS: usize = 30;
/// Newton stops once no unknown moves by more than `2^-44` of its size: well
/// below the operator's floor, and only a proposal, which the proof checks.
const NEWTON_STOP: f64 = 1.0 / 17_592_186_044_416.0;
/// Points per free parameter of the seed grid, as `PatchDistance`'s first
/// binary64 grid (`fast_nearest`).
const SEED_GRID: usize = 9;
/// Epsilon-inflations of a Krawczyk box before a proof gives up (Rump,
/// "Solving algebraic problems with high accuracy", 1983).
const INFLATIONS: usize = 8;

thread_local! {
    /// Products made since the last charge, on this thread.
    static PRODUCTS: Cell<u64> = const { Cell::new(0) };
}

fn count(products: u64) {
    PRODUCTS.with(|c| c.set(c.get().saturating_add(products)));
}

/// Charges the products made since the last charge.
fn charge(dec: &mut Decider) -> Result<(), Refusal> {
    let products = PRODUCTS.with(|c| c.replace(0));
    dec.charge_work(products)
}

fn refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new("field.spline-spine", detail.into())
}

fn definition_refusal(detail: &str) -> Refusal {
    Refusal::new("field.spline-spine-definition", detail)
}

fn zero() -> Interval {
    Interval::point(0.0).expect("zero")
}

fn point(value: f64) -> Interval {
    Interval::point(value).unwrap_or_else(Interval::entire)
}

fn magnitude(x: Interval) -> f64 {
    x.lo().abs().max(x.hi().abs())
}

/// The interval `[-m, m]`.
fn spread(m: f64) -> Interval {
    Interval::new(-m, m).unwrap_or_else(Interval::entire)
}

fn is_zero(x: Interval) -> bool {
    x.lo() == 0.0 && x.hi() == 0.0
}

fn product(a: Interval, b: Interval) -> Interval {
    count(1);
    a * b
}

/// Arithmetic the system is written in: binary64 for the trace, polynomials
/// in a span's variable with interval coefficients for the proofs.
trait Num: Clone {
    fn constant(value: Interval) -> Self;
    fn plus(&self, other: &Self) -> Self;
    fn minus(&self, other: &Self) -> Self;
    fn times(&self, other: &Self) -> Self;
    fn scaled(&self, k: Interval) -> Self;
}

impl Num for f64 {
    fn constant(value: Interval) -> Self {
        value.mid()
    }
    fn plus(&self, other: &Self) -> Self {
        self + other
    }
    fn minus(&self, other: &Self) -> Self {
        self - other
    }
    fn times(&self, other: &Self) -> Self {
        count(1);
        self * other
    }
    fn scaled(&self, k: Interval) -> Self {
        count(1);
        self * k.mid()
    }
}

impl Num for Interval {
    fn constant(value: Interval) -> Self {
        value
    }
    fn plus(&self, other: &Self) -> Self {
        *self + *other
    }
    fn minus(&self, other: &Self) -> Self {
        *self - *other
    }
    fn times(&self, other: &Self) -> Self {
        product(*self, *other)
    }
    fn scaled(&self, k: Interval) -> Self {
        product(*self, k)
    }
}

/// The value and first two coefficients of a directional Taylor expansion.
/// The product rule here gives the Hessian applied twice to one direction,
/// without forming a separate Hessian for each of the eight equations.
#[derive(Clone)]
struct Directional([Interval; 3]);

impl Num for Directional {
    fn constant(value: Interval) -> Self {
        Self([value, zero(), zero()])
    }
    fn plus(&self, other: &Self) -> Self {
        Self(add3(&self.0, &other.0))
    }
    fn minus(&self, other: &Self) -> Self {
        Self(sub3(&self.0, &other.0))
    }
    fn times(&self, other: &Self) -> Self {
        Self(std::array::from_fn(|k| {
            (0..=k).fold(zero(), |s, i| s + product(self.0[i], other.0[k - i]))
        }))
    }
    fn scaled(&self, k: Interval) -> Self {
        Self(mul3(&k, &self.0))
    }
}

/// A polynomial in `s` on `[-1, 1]` with interval coefficients, at most
/// `ORDER` in degree. It stands for every polynomial whose coefficients lie in
/// the intervals; a term dropped above `ORDER` is bounded by its coefficient's
/// magnitude on `[-1, 1]` and added to the constant term.
#[derive(Clone, Debug)]
struct Tm(Vec<Interval>);

impl Tm {
    /// Bounds on the polynomial over `[-1, 1]`.
    fn range(&self) -> Interval {
        self.0[1..]
            .iter()
            .fold(self.0[0], |out, c| out + spread(magnitude(*c)))
    }

    /// Bounds on the polynomial over `s`, by Horner's rule.
    fn at(&self, s: Interval) -> Interval {
        self.0
            .iter()
            .rev()
            .fold(zero(), |out, c| product(out, s) + *c)
    }
}

impl Num for Tm {
    fn constant(value: Interval) -> Self {
        Tm(vec![value])
    }
    fn plus(&self, other: &Self) -> Self {
        let (long, short) = if self.0.len() >= other.0.len() {
            (self, other)
        } else {
            (other, self)
        };
        let mut out = long.0.clone();
        for (o, c) in out.iter_mut().zip(&short.0) {
            *o = *o + *c;
        }
        Tm(out)
    }
    fn minus(&self, other: &Self) -> Self {
        self.plus(&Tm(other.0.iter().map(|c| -*c).collect()))
    }
    fn times(&self, other: &Self) -> Self {
        let length = (self.0.len() + other.0.len() - 1).min(ORDER + 1);
        let mut out = vec![zero(); length];
        let mut dropped = zero();
        for (i, a) in self.0.iter().enumerate() {
            if is_zero(*a) {
                continue;
            }
            for (j, b) in other.0.iter().enumerate() {
                let term = product(*a, *b);
                if i + j <= ORDER {
                    out[i + j] = out[i + j] + term;
                } else {
                    dropped = dropped + spread(magnitude(term));
                }
            }
        }
        out[0] = out[0] + dropped;
        Tm(out)
    }
    fn scaled(&self, k: Interval) -> Self {
        Tm(self.0.iter().map(|c| product(*c, k)).collect())
    }
}

fn add3<N: Num>(a: &[N; 3], b: &[N; 3]) -> [N; 3] {
    std::array::from_fn(|k| a[k].plus(&b[k]))
}

fn sub3<N: Num>(a: &[N; 3], b: &[N; 3]) -> [N; 3] {
    std::array::from_fn(|k| a[k].minus(&b[k]))
}

fn mul3<N: Num>(k: &N, a: &[N; 3]) -> [N; 3] {
    std::array::from_fn(|i| k.times(&a[i]))
}

fn cross<N: Num>(a: &[N; 3], b: &[N; 3]) -> [N; 3] {
    [
        a[1].times(&b[2]).minus(&a[2].times(&b[1])),
        a[2].times(&b[0]).minus(&a[0].times(&b[2])),
        a[0].times(&b[1]).minus(&a[1].times(&b[0])),
    ]
}

fn dot<N: Num>(a: &[N; 3], b: &[N; 3]) -> N {
    a[0].times(&b[0])
        .plus(&a[1].times(&b[1]))
        .plus(&a[2].times(&b[2]))
}

/// `W X' - W' X`, a partial of `S = X / W` times `W^2`.
fn tangent<N: Num>(w: &N, x_d: &[N; 3], w_d: &N, x: &[N; 3]) -> [N; 3] {
    sub3(&mul3(w, x_d), &mul3(w_d, x))
}

/// One support: the power coefficients of its homogeneous coordinates
/// `(X, Y, Z, W)` and of their partials on its span's unit square, and `s r`.
#[derive(Clone, Debug)]
struct Support {
    degrees: [usize; 2],
    /// `tables[kind][i * (n + 1) + j]` holds the coefficients of `u^i v^j` of
    /// the value (kind 0), `_u`, `_v`, `_uu`, `_uv` and `_vv` (kinds 1 to 5).
    tables: [Vec<[Interval; 4]>; 6],
    /// The side times the radius.
    side: Interval,
    /// Radius for measuring back to the original contact from the centre.
    contact_side: Interval,
}

fn binomial(n: usize, k: usize) -> Rational {
    (0..k).fold(Rational::one(), |out, i| {
        &(&out * &Rational::from_i64((n - i) as i64)) / &Rational::from_i64((i + 1) as i64)
    })
}

fn enclose(value: &Real, dec: &mut Decider) -> Result<Interval, Refusal> {
    Ok(match value.rational_ref() {
        Some(rational) => Interval::from_rational(rational),
        None => dec.enclose_budgeted_at(value, 128)?.to_interval(),
    })
}

impl Support {
    fn new(patch: &RationalPatch, side: Interval, dec: &mut Decider) -> Result<Self, Refusal> {
        let degrees = patch.carrier().degrees();
        if degrees.contains(&0) {
            return Err(definition_refusal(
                "a support of degree zero along a side has no area",
            ));
        }
        let mut spans = patch.bezier_spans(dec)?;
        let Some((_, controls)) = spans.next() else {
            return Err(definition_refusal("a support has no span"));
        };
        if spans.next().is_some() {
            return Err(definition_refusal(
                "a support of several spans is not supported; split it into whole patches",
            ));
        }
        let [m, n] = degrees;
        let rows = n + 1;
        let mut net = Vec::with_capacity(controls.len());
        for control in controls {
            let mut c = [zero(); 4];
            for (out, value) in c.iter_mut().zip(control) {
                *out = enclose(value, dec)?;
            }
            if c[3].lo() <= 0.0 {
                return Err(definition_refusal(
                    "a support weight is not proved positive",
                ));
            }
            net.push(c);
        }
        // `B_i^m(u) = sum_{a >= i} C(m, a) C(a, i) (-1)^(a - i) u^a`.
        let basis = degrees.map(|degree| {
            (0..=degree)
                .map(|a| {
                    let outer = binomial(degree, a);
                    (0..=a)
                        .map(|i| {
                            let term = &outer * &binomial(a, i);
                            let signed = if (a - i) % 2 == 0 { term } else { -term };
                            Interval::from_rational(&signed)
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        });
        let mut value = vec![[zero(); 4]; (m + 1) * rows];
        for a in 0..=m {
            for b in 0..=n {
                dec.charge_work(5 * (a as u64 + 1) * (b as u64 + 1))?;
                let mut out = [zero(); 4];
                for i in 0..=a {
                    for j in 0..=b {
                        let factor = basis[0][a][i] * basis[1][b][j];
                        for (o, c) in out.iter_mut().zip(net[i * rows + j]) {
                            *o = *o + factor * c;
                        }
                    }
                }
                value[a * rows + b] = out;
            }
        }
        let along_u =
            |table: &[[Interval; 4]], dec: &mut Decider| -> Result<Vec<[Interval; 4]>, Refusal> {
                let mut out = vec![[zero(); 4]; table.len()];
                for i in 0..m {
                    dec.charge_work(4 * rows as u64)?;
                    for j in 0..=n {
                        let k = point((i + 1) as f64);
                        out[i * rows + j] = table[(i + 1) * rows + j].map(|c| c * k);
                    }
                }
                Ok(out)
            };
        let along_v =
            |table: &[[Interval; 4]], dec: &mut Decider| -> Result<Vec<[Interval; 4]>, Refusal> {
                let mut out = vec![[zero(); 4]; table.len()];
                for i in 0..=m {
                    dec.charge_work(4 * n as u64)?;
                    for j in 0..n {
                        let k = point((j + 1) as f64);
                        out[i * rows + j] = table[i * rows + j + 1].map(|c| c * k);
                    }
                }
                Ok(out)
            };
        let u = along_u(&value, dec)?;
        let v = along_v(&value, dec)?;
        let uu = along_u(&u, dec)?;
        let uv = along_v(&u, dec)?;
        let vv = along_v(&v, dec)?;
        Ok(Self {
            degrees,
            tables: [value, u, v, uu, uv, vv],
            side,
            contact_side: side,
        })
    }

    /// The first `kinds` tables' values at `(u, v)`.
    fn values<N: Num>(&self, u: &N, v: &N, kinds: usize) -> Vec<[N; 4]> {
        let [m, n] = self.degrees;
        let rows = n + 1;
        let one = N::constant(point(1.0));
        let mut up = vec![one.clone(), u.clone()];
        for i in 2..=m {
            up.push(up[i - 1].times(u));
        }
        let mut vp = vec![one, v.clone()];
        for j in 2..=n {
            vp.push(vp[j - 1].times(v));
        }
        let mut monomials = Vec::with_capacity((m + 1) * rows);
        for ui in up.iter().take(m + 1) {
            for (j, vj) in vp.iter().take(n + 1).enumerate() {
                monomials.push(if j == 0 { ui.clone() } else { ui.times(vj) });
            }
        }
        let empty = N::constant(zero());
        (0..kinds)
            .map(|kind| {
                let mut out: [N; 4] = std::array::from_fn(|_| empty.clone());
                for (coefficients, monomial) in self.tables[kind].iter().zip(&monomials) {
                    for (o, c) in out.iter_mut().zip(coefficients) {
                        if !is_zero(*c) {
                            *o = o.plus(&monomial.scaled(*c));
                        }
                    }
                }
                out
            })
            .collect()
    }

    /// The support's four equations at centre `c`, parameters `(u, v)` and
    /// length `l`, and with `jacobian` their partials in `c`, `u`, `v` and `l`.
    fn equations<N: Num>(
        &self,
        c: &[N; 3],
        u: &N,
        v: &N,
        l: &N,
        jacobian: bool,
    ) -> ([N; 4], Option<[[N; 6]; 4]>) {
        let f = self.values(u, v, if jacobian { 6 } else { 3 });
        let split = |k: usize| -> ([N; 3], N) {
            let q = &f[k];
            ([q[0].clone(), q[1].clone(), q[2].clone()], q[3].clone())
        };
        let (x, w) = split(0);
        let (xu, wu) = split(1);
        let (xv, wv) = split(2);
        let tu = tangent(&w, &xu, &wu, &x);
        let tv = tangent(&w, &xv, &wv, &x);
        let m = cross(&tu, &tv);
        let gap = sub3(&mul3(&w, c), &x);
        let wm = mul3(&w, &m);
        let values = [
            l.times(&gap[0]).minus(&wm[0].scaled(self.side)),
            l.times(&gap[1]).minus(&wm[1].scaled(self.side)),
            l.times(&gap[2]).minus(&wm[2].scaled(self.side)),
            l.times(l).minus(&dot(&m, &m)),
        ];
        if !jacobian {
            return (values, None);
        }
        let (xuu, wuu) = split(3);
        let (xuv, wuv) = split(4);
        let (xvv, wvv) = split(5);
        // The partials of `T_u` and `T_v`.
        let tuu = tangent(&w, &xuu, &wuu, &x);
        let tuv = sub3(
            &add3(&mul3(&wv, &xu), &mul3(&w, &xuv)),
            &add3(&mul3(&wuv, &x), &mul3(&wu, &xv)),
        );
        let tvu = sub3(
            &add3(&mul3(&wu, &xv), &mul3(&w, &xuv)),
            &add3(&mul3(&wuv, &x), &mul3(&wv, &xu)),
        );
        let tvv = tangent(&w, &xvv, &wvv, &x);
        let mu = add3(&cross(&tuu, &tv), &cross(&tu, &tvu));
        let mv = add3(&cross(&tuv, &tv), &cross(&tu, &tvv));
        let lw = l.times(&w);
        let empty = N::constant(zero());
        let partial = |k: usize, wd: &N, xd: &[N; 3], md: &[N; 3]| -> N {
            l.times(&wd.times(&c[k]).minus(&xd[k]))
                .minus(&wd.times(&m[k]).plus(&w.times(&md[k])).scaled(self.side))
        };
        let row = |k: usize| -> [N; 6] {
            let mut out: [N; 6] = std::array::from_fn(|_| empty.clone());
            out[k] = lw.clone();
            out[3] = partial(k, &wu, &xu, &mu);
            out[4] = partial(k, &wv, &xv, &mv);
            out[5] = gap[k].clone();
            out
        };
        let two = point(-2.0);
        let last = [
            empty.clone(),
            empty.clone(),
            empty.clone(),
            dot(&m, &mu).scaled(two),
            dot(&m, &mv).scaled(two),
            l.scaled(point(2.0)),
        ];
        let rows = [row(0), row(1), row(2), last];
        (values, Some(rows))
    }

    /// The ball's centre over the support point `(u, v)` and `|M|`, in
    /// binary64; `None` where the normal vanishes.
    fn centre(&self, u: f64, v: f64) -> Option<([f64; 3], f64)> {
        let f = self.values(&u, &v, 3);
        let x = [f[0][0], f[0][1], f[0][2]];
        let tu = tangent(&f[0][3], &[f[1][0], f[1][1], f[1][2]], &f[1][3], &x);
        let tv = tangent(&f[0][3], &[f[2][0], f[2][1], f[2][2]], &f[2][3], &x);
        let m = cross(&tu, &tv);
        let length = dot(&m, &m).sqrt();
        if !length.is_finite() || length == 0.0 {
            return None;
        }
        let side = self.side.mid();
        count(3);
        Some((
            std::array::from_fn(|k| x[k] / f[0][3] + side * m[k] / length),
            length,
        ))
    }
}

/// Where in the nine values each quantity sits: the centre, `A`'s two
/// parameters, `B`'s two, then `L_A` and `L_B`.
const CENTRE: usize = 0;
const PARAMETERS: usize = 3;
const LENGTHS: usize = 7;

/// The eight unknowns when value `pinned` is held.
fn free(pinned: usize) -> [usize; 8] {
    let mut out = [0; 8];
    let mut k = 0;
    for index in 0..9 {
        if index != pinned {
            out[k] = index;
            k += 1;
        }
    }
    out
}

/// Gaussian elimination with partial pivoting, in binary64.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let size = b.len();
    for column in 0..size {
        let pivot =
            (column..size).max_by(|&i, &j| a[i][column].abs().total_cmp(&a[j][column].abs()))?;
        if a[pivot][column] == 0.0 || !a[pivot][column].is_finite() {
            return None;
        }
        a.swap(column, pivot);
        b.swap(column, pivot);
        for row in column + 1..size {
            let factor = a[row][column] / a[column][column];
            count((size - column + 1) as u64);
            for k in column..size {
                a[row][k] -= factor * a[column][k];
            }
            b[row] -= factor * b[column];
        }
    }
    let mut x = vec![0.0; size];
    for row in (0..size).rev() {
        let mut sum = b[row];
        for k in row + 1..size {
            sum -= a[row][k] * x[k];
        }
        count((size - row) as u64);
        x[row] = sum / a[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// The inverse of a square binary64 matrix, column by column.
fn inverse(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let size = a.len();
    let mut columns = Vec::with_capacity(size);
    for column in 0..size {
        let unit = (0..size)
            .map(|k| f64::from(u8::from(k == column)))
            .collect();
        columns.push(solve(a.to_vec(), unit)?);
    }
    Some(
        (0..size)
            .map(|row| (0..size).map(|column| columns[column][row]).collect())
            .collect(),
    )
}

/// A proved Krawczyk box: `K` inside `Y = [-radius, radius]`, the
/// preconditioner and `I - C J`.
struct Krawczyk {
    image: [Interval; 8],
    radius: [f64; 8],
    preconditioner: Vec<Vec<f64>>,
    contraction: [[Interval; 8]; 8],
    norm: f64,
}

/// The Krawczyk operator `K = -C g + (I - C J) Y` about zero, for residual
/// bounds `g` and the Jacobian `jacobian` over every point the unknowns can
/// take, grown by epsilon-inflation from the residual until `K` lies inside
/// `Y`. `scale` sizes each unknown's smallest box.
fn krawczyk(
    g: &[Interval; 8],
    jacobian: &[[Interval; 9]; 8],
    pinned: usize,
    scale: &[f64; 8],
) -> Option<Krawczyk> {
    let unknowns = free(pinned);
    let middle: Vec<Vec<f64>> = (0..8)
        .map(|row| unknowns.iter().map(|&k| jacobian[row][k].mid()).collect())
        .collect();
    let preconditioner = inverse(&middle)?;
    let step: [Interval; 8] = std::array::from_fn(|i| {
        -(0..8).fold(zero(), |sum, k| {
            sum + product(point(preconditioner[i][k]), g[k])
        })
    });
    let contraction: [[Interval; 8]; 8] = std::array::from_fn(|i| {
        std::array::from_fn(|c| {
            let unit = point(f64::from(u8::from(i == c)));
            unit - (0..8).fold(zero(), |sum, k| {
                sum + product(point(preconditioner[i][k]), jacobian[k][unknowns[c]])
            })
        })
    });
    // Four units in the last place of each unknown's size: room for the
    // operator's own outward rounding.
    let floor: [f64; 8] = std::array::from_fn(|i| 4.0 * f64::EPSILON * scale[i]);
    let mut radius: [f64; 8] = std::array::from_fn(|i| 2.0 * magnitude(step[i]) + floor[i]);
    for _ in 0..INFLATIONS {
        if radius.iter().any(|r| !r.is_finite() || *r > BOX) {
            return None;
        }
        let spread: [Interval; 8] = std::array::from_fn(|i| {
            (0..8).fold(zero(), |sum, c| {
                sum + product(contraction[i][c], spread(radius[c]))
            })
        });
        let image = std::array::from_fn(|i| step[i] + spread[i]);
        let norm = (0..8)
            .map(|i| (point(magnitude(spread[i])) / point(radius[i])).hi())
            .fold(0.0, f64::max);
        if norm < 1.0
            && image
                .iter()
                .zip(&radius)
                .all(|(k, r)| k.lo() > -r && k.hi() < *r)
        {
            return Some(Krawczyk {
                image,
                radius,
                preconditioner,
                contraction,
                norm,
            });
        }
        radius = std::array::from_fn(|i| 2.0 * magnitude(image[i]) + floor[i]);
    }
    None
}

/// Bounds on the spine's slope `d = dz / dtau` over a proved span, from
/// `J_x d = -F_tau`: `d = -C F_tau + (I - C J) d` for the span's
/// preconditioner `C`. A proved box `Y` gives `sum_c |I - C J|_ic y_c < y_i`
/// for its radii `y` (the operator's spread lies inside `Y`), so `I - C J`
/// contracts in the maximum norm weighted by `y`, and every solution has
/// `|d|_y <= |C F_tau|_y / (1 - |I - C J|_y)`, the contraction bound for
/// preconditioned interval linear systems (Neumaier, "Interval Methods for
/// Systems of Equations", 1990, chapter 4). The identity then narrows the
/// bounds, eight times at most.
fn solve_derivative(proof: &Krawczyk, load: [Interval; 8]) -> Option<[Interval; 8]> {
    let pull: [Interval; 8] = std::array::from_fn(|i| {
        -(0..8).fold(zero(), |sum, k| {
            sum + product(point(proof.preconditioner[i][k]), load[k])
        })
    });
    let y = &proof.radius;
    let largest = (0..8)
        .map(|i| (point(magnitude(pull[i])) / point(y[i])).hi())
        .fold(0.0, f64::max);
    let bound = (point(largest) / (point(1.0) - point(proof.norm))).hi();
    let mut d: [Interval; 8] =
        std::array::from_fn(|i| spread(product(point(bound), point(y[i])).hi()));
    for _ in 0..INFLATIONS {
        let next: [Interval; 8] = std::array::from_fn(|i| {
            (0..8).fold(pull[i], |sum, c| {
                sum + product(proof.contraction[i][c], d[c])
            })
        });
        for (out, n) in d.iter_mut().zip(next) {
            *out = Interval::new(out.lo().max(n.lo()), out.hi().min(n.hi()))?;
        }
    }
    Some(d)
}

/// The two supports, the radius and the side choices: the spine's system.
#[derive(Clone, Debug)]
struct System {
    supports: [Support; 2],
    axis: usize,
}

impl System {
    /// The pinned value of a span: `A`'s parameter `tau`.
    fn tau(&self) -> usize {
        PARAMETERS + self.axis
    }

    /// The contact parameters other than `tau`.
    fn contacts(&self) -> [usize; 3] {
        [PARAMETERS + 1 - self.axis, PARAMETERS + 2, PARAMETERS + 3]
    }

    /// The eight equations at the nine values, and with `jacobian` their
    /// partials in all nine.
    fn system<N: Num>(&self, z: &[N; 9], jacobian: bool) -> ([N; 8], Option<[[N; 9]; 8]>) {
        let c = [z[0].clone(), z[1].clone(), z[2].clone()];
        let mut values: [N; 8] = std::array::from_fn(|_| N::constant(zero()));
        let mut rows: Option<[[N; 9]; 8]> =
            jacobian.then(|| std::array::from_fn(|_| std::array::from_fn(|_| N::constant(zero()))));
        for (which, support) in self.supports.iter().enumerate() {
            let (u, v, l) = (
                PARAMETERS + 2 * which,
                PARAMETERS + 2 * which + 1,
                LENGTHS + which,
            );
            let (e, j) = support.equations(&c, &z[u], &z[v], &z[l], jacobian);
            for (k, value) in e.into_iter().enumerate() {
                values[4 * which + k] = value;
            }
            if let (Some(rows), Some(j)) = (rows.as_mut(), j) {
                for (k, partials) in j.into_iter().enumerate() {
                    let columns = [CENTRE, CENTRE + 1, CENTRE + 2, u, v, l];
                    for (column, partial) in columns.into_iter().zip(partials) {
                        rows[4 * which + k][column] = partial;
                    }
                }
            }
        }
        (values, rows)
    }

    /// Newton's method on the unknowns with value `pinned` held, from `z`.
    fn newton(&self, pinned: usize, mut z: [f64; 9]) -> Option<[f64; 9]> {
        let unknowns = free(pinned);
        for _ in 0..NEWTON_STEPS {
            let (f, j) = self.system(&z, true);
            let j = j.expect("a Jacobian");
            let a = (0..8)
                .map(|row| unknowns.iter().map(|&k| j[row][k]).collect())
                .collect();
            let step = solve(a, f.to_vec())?;
            let mut largest: f64 = 0.0;
            for (k, &index) in unknowns.iter().enumerate() {
                z[index] -= step[k];
                largest = largest.max(step[k].abs() / (1.0 + z[index].abs()));
            }
            if !z.iter().all(|v| v.is_finite()) {
                return None;
            }
            if largest <= NEWTON_STOP {
                return (z[LENGTHS] > 0.0 && z[LENGTHS + 1] > 0.0).then_some(z);
            }
        }
        None
    }

    /// The spine's slope `dz / dtau` at `z`, in binary64.
    fn slope(&self, z: &[f64; 9]) -> Option<[f64; 9]> {
        let pinned = self.tau();
        let unknowns = free(pinned);
        let (_, j) = self.system(z, true);
        let j = j.expect("a Jacobian");
        let a = (0..8)
            .map(|row| unknowns.iter().map(|&k| j[row][k]).collect())
            .collect();
        let b = (0..8).map(|row| -j[row][pinned]).collect();
        let d = solve(a, b)?;
        let mut out = [0.0; 9];
        out[pinned] = 1.0;
        for (k, &index) in unknowns.iter().enumerate() {
            out[index] = d[k];
        }
        Some(out)
    }

    /// A spine point at `tau`: the closest pair of ball centres over a grid
    /// on the three other contact parameters, then Newton's method.
    fn seed(&self, tau: f64) -> Option<[f64; 9]> {
        let grid = |i: usize| i as f64 / (SEED_GRID - 1) as f64;
        let mut a_params = [0.0; 2];
        a_params[self.axis] = tau;
        let mut first = Vec::new();
        for i in 0..SEED_GRID {
            a_params[1 - self.axis] = grid(i);
            if let Some(centre) = self.supports[0].centre(a_params[0], a_params[1]) {
                first.push((a_params, centre));
            }
        }
        let mut best: Option<(f64, [f64; 9])> = None;
        for i in 0..SEED_GRID {
            for j in 0..SEED_GRID {
                let Some((cb, lb)) = self.supports[1].centre(grid(i), grid(j)) else {
                    continue;
                };
                for (a, (ca, la)) in &first {
                    count(3);
                    let gap: f64 = (0..3)
                        .map(|k| {
                            let d = ca[k] - cb[k];
                            d * d
                        })
                        .sum();
                    if best.as_ref().is_none_or(|(g, _)| gap < *g) {
                        let c: [f64; 3] = std::array::from_fn(|k| (ca[k] + cb[k]) / 2.0);
                        best = Some((
                            gap,
                            [c[0], c[1], c[2], a[0], a[1], grid(i), grid(j), *la, lb],
                        ));
                    }
                }
            }
        }
        let (_, z) = best?;
        self.newton(self.tau(), z)
    }

    /// Follows the spine from `seed` both ways to its two ends. Returns the
    /// traced points in increasing `tau`, ends included, and the two ends.
    fn trace(
        &self,
        seed: [f64; 9],
        dec: &mut Decider,
    ) -> Result<(Vec<[f64; 9]>, [End; 2]), Refusal> {
        let tau = self.tau();
        let mut found: Vec<Option<End>> = vec![None, None];
        let mut points = vec![seed];
        for (side, direction) in [-1.0, 1.0].into_iter().enumerate() {
            let mut z = seed;
            loop {
                charge(dec)?;
                let here = z[tau];
                if (direction < 0.0 && here <= 0.0) || (direction > 0.0 && here >= 1.0) {
                    found[side] = Some(End::own(tau, z));
                    break;
                }
                let next = (here + direction * TRACE_STEP).clamp(0.0, 1.0);
                let lost = || refusal(format!("Newton's method lost the spine near tau = {here}"));
                let proposed = self.slope(&z);
                charge(dec)?;
                let d = proposed.ok_or_else(lost)?;
                count(9);
                let mut guess: [f64; 9] = std::array::from_fn(|k| z[k] + (next - here) * d[k]);
                guess[tau] = next;
                let proposed = self.newton(tau, guess);
                charge(dec)?;
                let ahead = proposed.ok_or_else(lost)?;
                if let Some(end) = self.leaves(&z, &ahead) {
                    points.push(end.z);
                    found[side] = Some(end);
                    break;
                }
                points.push(ahead);
                z = ahead;
            }
        }
        points.sort_by(|a, b| a[tau].total_cmp(&b[tau]));
        let [Some(low), Some(high)] = [found[0].take(), found[1].take()] else {
            unreachable!("both directions end")
        };
        Ok((points, [low, high]))
    }

    /// Where the spine between traced points `z` and `ahead` first leaves a
    /// support's square, if it does: the root with the leaving parameter at
    /// its border.
    fn leaves(&self, z: &[f64; 9], ahead: &[f64; 9]) -> Option<End> {
        let tau = self.tau();
        let mut best: Option<End> = None;
        for q in self.contacts() {
            let border = if ahead[q] < 0.0 {
                0.0
            } else if ahead[q] > 1.0 {
                1.0
            } else {
                continue;
            };
            let fraction = (border - z[q]) / (ahead[q] - z[q]);
            count(9);
            let mut guess: [f64; 9] = std::array::from_fn(|k| z[k] + fraction * (ahead[k] - z[k]));
            guess[q] = border;
            let Some(root) = self.newton(q, guess) else {
                continue;
            };
            let nearer = best
                .as_ref()
                .is_none_or(|b| (root[tau] - z[tau]).abs() < (b.z[tau] - z[tau]).abs());
            if nearer {
                best = Some(End {
                    pinned: q,
                    border,
                    z: root,
                });
            }
        }
        best.or_else(|| {
            // A parameter left, but no border root was found: the trace stops
            // at the last point inside, and the end proof refuses.
            self.contacts()
                .into_iter()
                .any(|q| !(0.0..=1.0).contains(&ahead[q]))
                .then_some(End {
                    pinned: self.contacts()[0],
                    border: f64::NAN,
                    z: *z,
                })
        })
    }

    /// Proves one root with value `pinned` held exactly at `z[pinned]`, in a
    /// box round `z`.
    fn prove_point(
        &self,
        pinned: usize,
        z: &[f64; 9],
        dec: &mut Decider,
    ) -> Result<Option<[Interval; 9]>, Refusal> {
        let unknowns = free(pinned);
        let at: [Tm; 9] = std::array::from_fn(|k| Tm::constant(point(z[k])));
        let (g, _) = self.system(&at, false);
        charge(dec)?;
        let mut around = at.clone();
        for &k in &unknowns {
            around[k] = Tm::constant(point(z[k]) + spread(BOX));
        }
        let (_, j) = self.system(&around, true);
        charge(dec)?;
        let g: [Interval; 8] = std::array::from_fn(|k| g[k].range());
        let j = j.expect("a Jacobian");
        let j: [[Interval; 9]; 8] =
            std::array::from_fn(|r| std::array::from_fn(|c| j[r][c].range()));
        let scale: [f64; 8] = std::array::from_fn(|i| 1.0 + z[unknowns[i]].abs());
        let proof = krawczyk(&g, &j, pinned, &scale);
        charge(dec)?;
        Ok(proof.map(|k| {
            let mut out: [Interval; 9] = std::array::from_fn(|i| point(z[i]));
            for (i, &index) in unknowns.iter().enumerate() {
                out[index] = point(z[index]) + k.image[i];
            }
            out
        }))
    }

    /// The traced spine at `t`, interpolated between traced points (and
    /// extended past the ends by the nearest two).
    fn guess(&self, points: &[[f64; 9]], t: f64) -> [f64; 9] {
        let tau = self.tau();
        let after = points.partition_point(|p| p[tau] < t);
        let (a, b) = match points.len() {
            0 | 1 => return points[0],
            len => {
                let i = after.clamp(1, len - 1);
                (&points[i - 1], &points[i])
            }
        };
        let span = b[tau] - a[tau];
        let fraction = if span == 0.0 {
            0.0
        } else {
            (t - a[tau]) / span
        };
        let mut out: [f64; 9] = std::array::from_fn(|k| a[k] + fraction * (b[k] - a[k]));
        out[tau] = t;
        out
    }

    /// One span's proof on `tau` in `[t0, t1]`, or `None` when the operator
    /// does not prove it.
    fn prove_span(
        &self,
        t: [f64; 2],
        points: &[[f64; 9]],
        dec: &mut Decider,
    ) -> Result<Option<SpineSpan>, Refusal> {
        let tau = self.tau();
        let unknowns = free(tau);
        let middle = (t[0] + t[1]) / 2.0;
        let half = (t[1] - t[0]) / 2.0;
        let mut solutions = Vec::with_capacity(TERMS);
        for s in NODES {
            let at = middle + half * s;
            let Some(z) = self.newton(tau, self.guess(points, at)) else {
                charge(dec)?;
                return Ok(None);
            };
            solutions.push(z);
        }
        let vandermonde: Vec<Vec<f64>> = NODES
            .iter()
            .map(|s| {
                let mut row = vec![1.0, *s];
                for k in 2..TERMS {
                    row.push(row[k - 1] * s);
                }
                count((TERMS - 2) as u64);
                row
            })
            .collect();
        let mut model = [[0.0; TERMS]; 8];
        for (i, &index) in unknowns.iter().enumerate() {
            let values = solutions.iter().map(|z| z[index]).collect();
            let Some(coefficients) = solve(vandermonde.clone(), values) else {
                charge(dec)?;
                return Ok(None);
            };
            model[i].copy_from_slice(&coefficients);
        }
        charge(dec)?;
        let models: [Tm; 8] = std::array::from_fn(|i| Tm(model[i].map(point).to_vec()));
        let mut inputs: [Tm; 9] = std::array::from_fn(|_| Tm::constant(zero()));
        inputs[tau] = Tm(vec![point(middle), point(half)]);
        for (i, &index) in unknowns.iter().enumerate() {
            inputs[index] = models[i].clone();
        }
        let (residual, _) = self.system(&inputs, false);
        charge(dec)?;
        let g: [Interval; 8] = std::array::from_fn(|k| residual[k].range());
        for &index in &unknowns {
            let mut widened = inputs[index].clone();
            widened.0[0] = widened.0[0] + spread(BOX);
            inputs[index] = widened;
        }
        let (_, j) = self.system(&inputs, true);
        charge(dec)?;
        let j = j.expect("a Jacobian");
        let j: [[Interval; 9]; 8] =
            std::array::from_fn(|r| std::array::from_fn(|c| j[r][c].range()));
        let scale: [f64; 8] = std::array::from_fn(|i| 1.0 + magnitude(models[i].range()));
        let Some(proof) = krawczyk(&g, &j, tau, &scale) else {
            charge(dec)?;
            return Ok(None);
        };
        let slope = solve_derivative(&proof, j.map(|row| row[tau]));
        // Twice differentiating F(z(tau)) = 0 gives the same linear
        // system with the directional second derivative of F as its load.
        let acceleration = slope.and_then(|d| {
            let mut z: [Directional; 9] = std::array::from_fn(|_| Directional::constant(zero()));
            z[tau] = Directional([Interval::new(t[0], t[1]).unwrap(), point(1.0), zero()]);
            for (i, index) in unknowns.into_iter().enumerate() {
                z[index] = Directional([models[i].range() + proof.image[i], d[i], zero()]);
            }
            let (g, _) = self.system(&z, false);
            let second = solve_derivative(&proof, g.map(|value| product(point(2.0), value.0[2])))?;
            Some([second[0], second[1], second[2]])
        });
        charge(dec)?;
        let eps = (0..3)
            .fold(zero(), |sum, i| {
                sum + point(magnitude(proof.image[i])).square()
            })
            .sqrt()
            .map_or(f64::INFINITY, |e| e.hi());
        // Scaling the uniqueness box preserves the weighted contraction,
        // and J already covers BOX in every coordinate. Retain the smaller
        // image for the tube, but allow an independently enclosed end root
        // to join in the larger box even when the model residual is tiny.
        let mut uniqueness = proof.radius;
        let mut largest = uniqueness.into_iter().fold(0.0, f64::max);
        while largest <= BOX / 2.0 {
            uniqueness = uniqueness.map(|r| r * 2.0);
            largest *= 2.0;
            count(9);
        }
        charge(dec)?;
        Ok(Some(SpineSpan {
            tau: t,
            model,
            offset: proof.image,
            radius: uniqueness,
            slope,
            acceleration,
            eps,
        }))
    }
}

/// A traced end: the root with value `pinned` at `border`.
#[derive(Clone, Debug)]
struct End {
    pinned: usize,
    border: f64,
    z: [f64; 9],
}

impl End {
    /// An end where `tau` itself reaches its border.
    fn own(tau: usize, z: [f64; 9]) -> Self {
        Self {
            pinned: tau,
            border: z[tau],
            z,
        }
    }
}

/// What the spine is: two supports, a radius, the side of each support the
/// centre lies on, the parameter of the first support that runs along the
/// spine, a seed value of it and the part's tolerance.
#[derive(Clone, Debug)]
pub struct SpineDefinition {
    /// The two supports, each one whole rational patch of one span.
    /// Parameters are the span's unit square, which is the patch's own
    /// parameters when its domain is the unit square.
    pub supports: [RationalPatch; 2],
    /// The ball's radius in millimetres.
    pub radius: Rational,
    /// The side of each support's natural normal `S_u x S_v` the centre lies on.
    pub sides: [Sign; 2],
    /// Which parameter of the first support runs along the spine, `tau`: 0 for
    /// `u`, 1 for `v`.
    pub axis: usize,
    /// A value of `tau` the spine passes, with both contacts inside their
    /// supports.
    pub seed: Rational,
    /// The part's tolerance in millimetres; the tube must be no wider than
    /// half of it, and no wider than 100 nm.
    pub tolerance: Rational,
}

/// One proved span: for every `tau` in `tau`, the spine's unknowns lie in
/// `model(s) + offset`, with `s = (tau - middle) / half`, and nowhere else in
/// `model(s) + [-radius, radius]`.
#[derive(Clone, Debug)]
pub struct SpineSpan {
    /// The span's range of `tau`, dyadic ends.
    pub tau: [f64; 2],
    /// Power coefficients in `s` of the eight unknowns, in the order centre,
    /// the free parameter of the first support, both parameters of the second,
    /// `L_A`, `L_B`.
    pub model: [[f64; TERMS]; 8],
    /// Where each unknown lies relative to its model.
    pub offset: [Interval; 8],
    /// The box in which the root is the only one.
    pub radius: [f64; 8],
    /// Bounds on `dz / dtau` over the span, where the operator contracts in
    /// the maximum norm.
    pub slope: Option<[Interval; 8]>,
    acceleration: Option<[Interval; 3]>,
    /// The tube's width on this span: the largest distance from the model's
    /// centre to the spine's centre, in millimetres.
    pub eps: f64,
}

impl SpineSpan {
    fn model_at(&self, i: usize, s: Interval) -> Interval {
        Tm(self.model[i].map(point).to_vec()).at(s)
    }

    /// The span's variable `s` for `tau`.
    fn variable(&self, tau: Interval) -> Interval {
        let middle = (self.tau[0] + self.tau[1]) / 2.0;
        let half = (self.tau[1] - self.tau[0]) / 2.0;
        (tau - point(middle)) / point(half)
    }

    /// Bounds on the unknowns over the span.
    fn ranges(&self) -> [Interval; 8] {
        std::array::from_fn(|i| Tm(self.model[i].map(point).to_vec()).range() + self.offset[i])
    }
}

/// A proved end: the spine point where a contact parameter, or `tau` itself,
/// is at its support's border.
#[derive(Clone, Debug)]
pub struct SpineEnd {
    /// Which of the four contact parameters is at its border: 0 and 1 for the
    /// first support's `u` and `v`, 2 and 3 for the second's.
    pub parameter: usize,
    /// The border value, 0 or 1.
    pub border: f64,
    /// Bounds on the nine values: centre, four parameters, two lengths.
    pub point: [Interval; 9],
}

/// A proved spine.
#[derive(Clone, Debug)]
pub struct Spine {
    spans: Vec<SpineSpan>,
    ends: [SpineEnd; 2],
    eps: f64,
    target: f64,
    system: System,
    radius: Interval,
    seed_root: [Interval; 9],
}

impl Spine {
    /// Proves the spine of `definition`: both ends, then every span, with a
    /// tube no wider than `min(100 nm, tolerance / 2)`.
    ///
    /// Refuses with `field.spline-spine-definition` for a definition it cannot
    /// pose, `field.spline-spine` where no spine is found or a proof fails,
    /// and `field.spline-spine-eps` with the width reached where every span is
    /// proved but some tube stays wider than the target.
    pub fn prove(definition: &SpineDefinition, dec: &mut Decider) -> Result<Self, Refusal> {
        Self::prove_offsets(definition, None, dec)
    }

    /// Proves a spline meeting a plane's exact rational parametrization.
    /// The selected planar chart is translated to the analytic offset,
    /// then enters the system with zero offset, as in M10's offset model.
    pub fn prove_spline_plane(
        definition: &SpineDefinition,
        plane_index: usize,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if plane_index >= definition.supports.len() {
            return Err(definition_refusal("the plane must be support zero or one"));
        }
        Self::prove_offsets(definition, Some(plane_index), dec)
    }

    fn prove_offsets(
        definition: &SpineDefinition,
        plane_index: Option<usize>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        PRODUCTS.with(|c| c.set(0));
        if definition.radius.sign() != Sign::Positive {
            return Err(definition_refusal("the radius must be positive"));
        }
        if definition.tolerance.sign() != Sign::Positive {
            return Err(definition_refusal("the tolerance must be positive"));
        }
        if definition.axis > 1 {
            return Err(definition_refusal("the axis must be 0 (u) or 1 (v)"));
        }
        if definition.seed.sign() == Sign::Negative || definition.seed > Rational::one() {
            return Err(definition_refusal("the seed must lie in [0, 1]"));
        }
        let radius = Interval::from_rational(&definition.radius);
        let mut supports = Vec::with_capacity(2);
        for (i, (patch, orientation)) in
            definition.supports.iter().zip(definition.sides).enumerate()
        {
            let side = match orientation {
                Sign::Positive => radius,
                Sign::Negative => -radius,
                Sign::Zero => return Err(definition_refusal("each side must be nonzero")),
            };
            let prepared = if plane_index == Some(i) {
                let surface = Surface::BSpline(patch.clone().into());
                let Some(Surface::Plane(plane)) = crate::flat_support(&surface) else {
                    return Err(definition_refusal(
                        "the analytic support must be a whole bilinear plane chart",
                    ));
                };
                let n = plane.normal;
                let d = match orientation {
                    Sign::Positive => definition.radius.clone(),
                    Sign::Negative => -&definition.radius,
                    Sign::Zero => unreachable!("the side was validated"),
                };
                let translation = if let Some(n) = n.as_rational()
                    && let Some(length) = crate::profile::sqrt(&n.norm_squared())
                {
                    n.scaled(&(&d / &length)).to_real()
                } else {
                    n.scaled(&(&Real::from(d) / &n.norm_squared().sqrt()))
                };
                let pose = Pose::new(Rotation::identity(), translation, dec)?;
                Support::new(&patch.transformed(&pose), zero(), dec).map(|mut p| {
                    p.contact_side = side;
                    p
                })
            } else {
                Support::new(patch, side, dec)
            };
            supports.push(prepared?);
        }
        let half_tolerance = &definition.tolerance / &Rational::from_i64(2);
        let target = TUBE_CAP.min(Interval::from_rational(&half_tolerance).lo());
        let system = System {
            supports: [supports[0].clone(), supports[1].clone()],
            axis: definition.axis,
        };
        let tau = system.tau();
        let seed_tau = Interval::from_rational(&definition.seed).mid();
        let proposed = system.seed(seed_tau);
        charge(dec)?;
        let seed = proposed.ok_or_else(|| refusal("no spine point was found at the seed"))?;
        if system
            .contacts()
            .iter()
            .any(|&q| !(0.0..=1.0).contains(&seed[q]))
        {
            return Err(refusal("the seed's contacts are not inside their supports"));
        }
        let (points, traced) = system.trace(seed, dec)?;
        let mut ends = Vec::with_capacity(2);
        for end in &traced {
            if end.border.is_nan() {
                return Err(refusal(format!(
                    "a contact leaves its support near tau = {} with no border root",
                    end.z[tau]
                )));
            }
            let proved = system
                .prove_point(end.pinned, &end.z, dec)?
                .ok_or_else(|| {
                    refusal(format!(
                        "the end near tau = {} is not proved by the Krawczyk operator",
                        end.z[tau]
                    ))
                })?;
            ends.push(SpineEnd {
                parameter: end.pinned - PARAMETERS,
                border: end.border,
                point: proved,
            });
        }
        let ends: [SpineEnd; 2] = [ends[0].clone(), ends[1].clone()];
        if !(ends[0].point[tau].hi() < seed_tau && seed_tau < ends[1].point[tau].lo()) {
            return Err(refusal("the border roots do not enclose the seed"));
        }
        let range = [ends[0].point[tau].lo(), ends[1].point[tau].hi()];
        let width = 0.5f64.powi(FIRST_LEVEL as i32);
        let first = (range[0] / width).floor().max(0.0) as i64;
        let last = (range[1] / width).ceil().min(1.0 / width) as i64;
        let mut spans = Vec::new();
        let mut widest_failure: Option<f64> = None;
        for j in first..last {
            let t = [j as f64 * width, (j + 1) as f64 * width];
            system.cover(
                t,
                FIRST_LEVEL,
                range,
                &points,
                target,
                &mut spans,
                &mut widest_failure,
                dec,
            )?;
        }
        if let Some(reached) = widest_failure {
            return Err(Refusal::new(
                "field.spline-spine-eps",
                format!(
                    "the tube reached {:.6e} nm after {} halvings of a span, above the target of {:.6e} nm",
                    nanometres(reached),
                    LAST_LEVEL - FIRST_LEVEL,
                    nanometres(target)
                ),
            ));
        }
        let mut spine = Self {
            eps: spans.iter().map(|s| s.eps).fold(0.0, f64::max),
            spans,
            ends,
            target,
            system,
            radius,
            seed_root: [zero(); 9],
        };
        spine.seed_root = spine.check(&spine.system, &seed, dec)?;
        charge(dec)?;
        Ok(spine)
    }

    /// The proofs that tie the spans together: joins, ends, contacts inside
    /// their supports, positive lengths.
    fn check(
        &self,
        system: &System,
        seed: &[f64; 9],
        dec: &mut Decider,
    ) -> Result<[Interval; 9], Refusal> {
        let tau = system.tau();
        let unknowns = free(tau);
        for pair in self.spans.windows(2) {
            let (a, b) = (&pair[0], &pair[1]);
            let inside = |from: &SpineSpan, s_from: f64, to: &SpineSpan, s_to: f64| {
                (0..8).all(|i| {
                    let gap = from.model_at(i, point(s_from)) + from.offset[i]
                        - to.model_at(i, point(s_to));
                    gap.lo() > -to.radius[i] && gap.hi() < to.radius[i]
                })
            };
            let join = a.tau[1];
            let refused = || refusal(format!("the spans meeting at tau = {join} do not join"));
            if join.to_bits() != b.tau[0].to_bits() {
                return Err(refused());
            }
            if inside(a, 1.0, b, -1.0) || inside(b, -1.0, a, 1.0) {
                continue;
            }
            // A span-wide enclosure can be too wide to fit in its neighbour,
            // even when their exact endpoint agrees. Isolate the root at the
            // shared parameter to tie both uniqueness boxes to that root.
            let mut guess = [join; 9];
            for (i, &index) in unknowns.iter().enumerate() {
                guess[index] = a.model_at(i, point(1.0)).mid();
            }
            let proposed = system.newton(tau, guess);
            charge(dec)?;
            let root = system
                .prove_point(tau, &proposed.ok_or_else(refused)?, dec)?
                .ok_or_else(refused)?;
            for (span, s) in [(a, 1.0), (b, -1.0)] {
                if !unknowns.iter().enumerate().all(|(i, &index)| {
                    let gap = root[index] - span.model_at(i, point(s));
                    gap.lo() > -span.radius[i] && gap.hi() < span.radius[i]
                }) {
                    return Err(refused());
                }
            }
        }
        let seed_root = system
            .prove_point(tau, seed, dec)?
            .ok_or_else(|| refusal("the seed root is not proved"))?;
        let roots = [
            ("end", &self.ends[0].point),
            ("end", &self.ends[1].point),
            ("seed", &seed_root),
        ];
        for (name, root) in roots {
            let at = root[tau];
            let mut held = false;
            for span in &self.spans {
                if at.hi() < span.tau[0] || at.lo() > span.tau[1] {
                    continue;
                }
                let s = span.variable(at);
                let inside = unknowns.iter().enumerate().all(|(i, &index)| {
                    let gap = root[index] - span.model_at(i, s);
                    gap.lo() > -span.radius[i] && gap.hi() < span.radius[i]
                });
                if !inside {
                    return Err(refusal(format!(
                        "the {name} near tau = {} is not the spine's point there",
                        at.mid()
                    )));
                }
                held = true;
            }
            if !held {
                return Err(refusal(format!("the {name} lies outside every span")));
            }
        }
        // Lengths stay positive, so each centre is on its chosen side.
        for span in &self.spans {
            let ranges = span.ranges();
            if ranges[6].lo() <= 0.0 || ranges[7].lo() <= 0.0 {
                return Err(refusal(format!(
                    "a normal may vanish on the span from tau = {}",
                    span.tau[0]
                )));
            }
        }
        // Contacts inside their squares. A parameter that crosses its border
        // at an end may leave the square on the spans nearest that end, where
        // its slope must have the sign that carries it inside.
        let mut excused = vec![[false; 3]; self.spans.len()];
        for (which, end) in self.ends.iter().enumerate() {
            let index = PARAMETERS + end.parameter;
            if index == tau {
                continue;
            }
            let Some(q) = system.contacts().iter().position(|&c| c == index) else {
                continue;
            };
            // Unknown order: centre, then the three contacts.
            let slot = 3 + q;
            // Inside the square means above the border 0 or below the border
            // 1; the spine lies above `tau` at the first end and below it at
            // the second.
            let into = if end.border == 0.0 { 1.0 } else { -1.0 };
            let along = if which == 0 { 1.0 } else { -1.0 };
            let order: Vec<usize> = if which == 0 {
                (0..self.spans.len()).collect()
            } else {
                (0..self.spans.len()).rev().collect()
            };
            for k in order {
                let span = &self.spans[k];
                let range = span.ranges()[slot];
                if range.lo() >= 0.0 && range.hi() <= 1.0 {
                    break;
                }
                let other = if end.border == 0.0 {
                    range.hi() <= 1.0
                } else {
                    range.lo() >= 0.0
                };
                let carried = span.slope.is_some_and(|d| {
                    let d = d[slot];
                    (into * along > 0.0 && d.lo() > 0.0) || (into * along < 0.0 && d.hi() < 0.0)
                });
                if !other || !carried {
                    return Err(refusal(format!(
                        "a contact is not proved inside its support near tau = {}",
                        span.tau[0]
                    )));
                }
                excused[k][q] = true;
            }
        }
        for (k, span) in self.spans.iter().enumerate() {
            let ranges = span.ranges();
            for q in 0..3 {
                let range = ranges[3 + q];
                if !excused[k][q] && (range.lo() < 0.0 || range.hi() > 1.0) {
                    return Err(refusal(format!(
                        "a contact is not proved inside its support on the span from tau = {}",
                        span.tau[0]
                    )));
                }
            }
        }
        Ok(seed_root)
    }

    /// The proved spans, in increasing `tau`.
    pub fn spans(&self) -> &[SpineSpan] {
        &self.spans
    }

    /// The two ends, the lower `tau` first.
    pub fn ends(&self) -> &[SpineEnd; 2] {
        &self.ends
    }

    /// The tube's width: the largest span width, in millimetres.
    pub fn eps(&self) -> f64 {
        self.eps
    }

    /// The width the tube had to reach, `min(100 nm, tolerance / 2)`.
    pub fn target(&self) -> f64 {
        self.target
    }

    /// Bounds on the nine values at `tau`: centre, the first support's `u`
    /// and `v`, the second's, `L_A` and `L_B`. `None` outside the spans.
    pub fn at(&self, tau: f64) -> Option<[Interval; 9]> {
        let span = self
            .spans
            .iter()
            .find(|s| s.tau[0] <= tau && tau <= s.tau[1])?;
        let pinned = self.system.tau();
        let s = span.variable(point(tau));
        let mut out = [point(tau); 9];
        for (i, &index) in free(pinned).iter().enumerate() {
            out[index] = span.model_at(i, s) + span.offset[i];
        }
        Some(out)
    }
}

impl System {
    /// Proves `tau` in `t` at the given level, halving where needed, and
    /// appends the proved spans in order. Halves that miss `range` are
    /// dropped.
    #[allow(clippy::too_many_arguments)]
    fn cover(
        &self,
        t: [f64; 2],
        level: u32,
        range: [f64; 2],
        points: &[[f64; 9]],
        target: f64,
        spans: &mut Vec<SpineSpan>,
        widest_failure: &mut Option<f64>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if t[1] < range[0] || t[0] > range[1] {
            return Ok(());
        }
        let proof = self.prove_span(t, points, dec)?;
        if let Some(span) = &proof
            && span.eps <= target
        {
            spans.push(span.clone());
            return Ok(());
        }
        if level >= LAST_LEVEL {
            return match proof {
                Some(span) => {
                    *widest_failure = Some(widest_failure.map_or(span.eps, |w| w.max(span.eps)));
                    spans.push(span);
                    Ok(())
                }
                None => Err(refusal(format!(
                    "the Krawczyk operator does not prove the span from tau = {} to {} after {} halvings",
                    t[0],
                    t[1],
                    LAST_LEVEL - FIRST_LEVEL
                ))),
            };
        }
        let middle = (t[0] + t[1]) / 2.0;
        for half in [[t[0], middle], [middle, t[1]]] {
            self.cover(
                half,
                level + 1,
                range,
                points,
                target,
                spans,
                widest_failure,
                dec,
            )?;
        }
        Ok(())
    }
}

/// A width in millimetres as nanometres, rounded up.
fn nanometres(width: f64) -> f64 {
    (point(width) * point(1e6)).hi()
}
