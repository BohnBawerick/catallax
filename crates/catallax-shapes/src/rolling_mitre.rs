//! The exact band a ball sweeps rolling over a mitre between two rounds of
//! one radius whose axes meet.
//!
//! Take unit round axes `a1`, `a2` meeting at `O`, with `g = a1.a2`, and the
//! frame `x = -(a1 + a2) / |a1 + a2|`, `y = a1 x a2 / |a1 x a2|`, `z = x x y`.
//! With `sa = sqrt((1 - g) / 2)` and `sc = sqrt((1 + g) / 2)`, a ball of radius
//! `rho` inside both rounds of radius `R` has its centre `k = R - rho` from
//! both axes. Those centres form the ellipse in the mitre plane
//!
//! ```text
//! c(t) = O + (k cos t / sa) x + k sin t y,      t in [0, pi/2],
//! ```
//!
//! and the contact rays toward the rounds are
//! `m1,2(t) = sa cos t x + sin t y -+ sc cos t z`. Both are square to `c'(t)`,
//! so the band is the shorter arc between them in the normal plane, a
//! constant-radius pipe. At `t = 0` that arc is a straight strip's end section
//! and at `t = pi/2` both rays are `y`: the band closes at the apex `O + R y`.
//!
//! With `tau = tan(t/2)`, `C = 1 - tau^2`, `S = 2 tau`, `W = 1 + tau^2` and
//! `G = sa^2 C^2 + S^2`, projecting the section circle from `-m1` and scaling
//! its chart by the arc's end gives, for `v` in `[0, 1]`,
//!
//! ```text
//! N_x = (k C / sa)(G + v^2 E) + rho ((G - v^2 E) sa C + 2 v sa sc^2 C^3)
//! N_y = k S (G + v^2 E) + rho ((G - v^2 E) S + 2 v sc^2 S C^2)
//! N_z = rho (2 v sc C G - (G - v^2 E) sc C)
//! D   = W (G + v^2 E),        E = sc^2 C^2,
//! ```
//!
//! and the point `O + (N_x x + N_y y + N_z z) / D`. That is an exact rational
//! patch of degree (6, 2): rows `v = 0` and `v = 1` are the contacts with the
//! first and second round, column `u = 0` is the end section and column `u = 1`
//! collapses to the apex. `G + E = W^2` and `W G` have positive power
//! coefficients, so every weight is positive and nothing is fitted.

use std::sync::Arc;

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};

use crate::{BSplineSurface, Circle, Ellipse, EnvelopeGrid, Isometry, Point3, Pose, Vector3, math};

mod query;

fn invalid(detail: &str) -> Refusal {
    Refusal::new("shapes.rolling-mitre", detail)
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &math::exact_sqrt(&v.norm_squared())))
}

// A polynomial in `tau` in the power basis. Real arithmetic keeps every
// node, so each coefficient that is structurally rational collapses to one
// leaf as it is made. Irrational factors join only after the basis change.
type Poly = Vec<Real>;

fn fold(value: Real) -> Real {
    value.as_rational().map_or(value, Real::from)
}
fn poly(values: &[i64]) -> Poly {
    values.iter().map(|v| Real::from_i64(*v)).collect()
}
fn times(a: &Poly, b: &Poly) -> Poly {
    let mut out = vec![Real::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] = fold(&out[i + j] + &(x * y));
        }
    }
    out
}
fn plus(a: &Poly, b: &Poly) -> Poly {
    (0..a.len().max(b.len()))
        .map(|i| match (a.get(i), b.get(i)) {
            (Some(x), Some(y)) => fold(x + y),
            (Some(x), None) => x.clone(),
            (None, Some(y)) => y.clone(),
            (None, None) => Real::zero(),
        })
        .collect()
}
fn scaled(a: &Poly, c: &Real) -> Poly {
    a.iter().map(|x| fold(x * c)).collect()
}
fn eval(a: &Poly, x: &Real) -> Real {
    a.iter()
        .rev()
        .fold(Real::zero(), |sum, c| fold(&(&sum * x) + c))
}

// A polynomial times a constant factor, and a sum of such terms.
#[derive(Clone)]
struct Term {
    base: Poly,
    factor: Real,
}
type Terms = Vec<Term>;

fn term(base: Poly, factor: Real) -> Terms {
    vec![Term {
        base,
        factor: fold(factor),
    }]
}

// The chart's numerators and denominator, each by power of `v`.
struct Chart {
    numerators: [[Terms; 3]; 3],
    denominator: [Poly; 3],
}

// The squares of `along` and `across` come directly from the axes' cosine,
// so they stay rational when it is.
#[derive(Clone, Copy)]
struct Shape<'a> {
    offset: &'a Real,
    radius: &'a Real,
    along: &'a Real,
    across: &'a Real,
    along_sq: &'a Real,
    across_sq: &'a Real,
}

fn chart(d: Shape) -> Chart {
    let (offset, radius, along, across) = (d.offset, d.radius, d.along, d.across);
    let c = poly(&[1, 0, -1]);
    let s = poly(&[0, 2]);
    let w = poly(&[1, 0, 1]);
    let c2 = times(&c, &c);
    let c3 = times(&c2, &c);
    let g = plus(&scaled(&c2, d.along_sq), &times(&s, &s));
    let e = scaled(&c2, d.across_sq);
    let two = Real::from_i64(2);
    let k_over = offset / along;
    let cg = times(&c, &g);
    let ce = times(&c, &e);
    let sg = times(&s, &g);
    let se = times(&s, &e);
    let x = [
        term(cg.clone(), &k_over + &(radius * along)),
        term(c3, &(&two * radius) * &(along * d.across_sq)),
        term(ce.clone(), &k_over - &(radius * along)),
    ];
    let y = [
        term(sg, offset + radius),
        term(times(&s, &c2), &(&two * radius) * d.across_sq),
        term(se, offset - radius),
    ];
    let z = [
        term(cg.clone(), -(radius * across)),
        term(cg, &(&two * radius) * across),
        term(ce, radius * across),
    ];
    Chart {
        numerators: [x, y, z],
        denominator: [times(&w, &g), Vec::new(), times(&w, &e)],
    }
}

// Bernstein coefficients of degree `n` of a power-basis polynomial.
fn bernstein(a: &Poly, n: usize) -> Vec<Real> {
    let binomial = |n: usize, k: usize| -> i64 {
        (0..k).fold(1, |value, i| value * (n - i) as i64 / (i as i64 + 1))
    };
    (0..=n)
        .map(|i| {
            (0..=i.min(a.len().saturating_sub(1))).fold(Real::zero(), |sum, j| {
                let factor = Real::from(
                    Rational::from_ratio(binomial(i, j), binomial(n, j)).expect("nonzero"),
                );
                fold(&sum + &(&factor * &a[j]))
            })
        })
        .collect()
}

// Each term's Bernstein coefficients times its factor, summed.
fn bernstein_terms(terms: &Terms, n: usize) -> Vec<Real> {
    let mut out: Vec<Option<Real>> = vec![None; n + 1];
    for t in terms {
        for (slot, coefficient) in out.iter_mut().zip(bernstein(&t.base, n)) {
            if coefficient.as_rational().is_some_and(|r| r.is_zero()) {
                continue;
            }
            let product = &coefficient * &t.factor;
            *slot = Some(match slot.take() {
                Some(sum) => &sum + &product,
                None => product,
            });
        }
    }
    out.into_iter()
        .map(|value| value.unwrap_or_else(Real::zero))
        .collect()
}

// Degree-two Bernstein coefficients in `v` of `p0 + p1 v + p2 v^2`.
fn across_v(p: &[Poly; 3]) -> [Poly; 3] {
    let half = Real::from(Rational::from_ratio(1, 2).expect("nonzero"));
    [
        p[0].clone(),
        plus(&p[0], &scaled(&p[1], &half)),
        plus(&plus(&p[0], &p[1]), &p[2]),
    ]
}
fn across_terms(p: &[Terms; 3]) -> [Terms; 3] {
    let half = Real::from(Rational::from_ratio(1, 2).expect("nonzero"));
    let halved: Terms = p[1]
        .iter()
        .map(|t| Term {
            base: t.base.clone(),
            factor: fold(&t.factor * &half),
        })
        .collect();
    [
        p[0].clone(),
        [p[0].clone(), halved].concat(),
        [p[0].clone(), p[1].clone(), p[2].clone()].concat(),
    ]
}

#[derive(Debug)]
struct Data {
    origin: Point3,
    axes: [Vector3; 2],
    frame: [Vector3; 3],
    round: Real,
    radius: Real,
    offset: Real,
    along: Real,
    across: Real,
    along_sq: Real,
    across_sq: Real,
    convex: bool,
    // A reflection leaves the frame left-handed: `z = -(x x y)`.
    reversing: bool,
    carrier: BSplineSurface,
}

/// The exact band of a ball narrower than two equal rounds, rolling over the
/// mitre where they meet from a straight strip's end to the mitre's apex.
#[derive(Clone, Debug)]
pub struct RollingMitre {
    data: Arc<Data>,
}

impl RollingMitre {
    /// Builds the band of a ball of `radius` between two rounds of radius
    /// `round` whose axes meet at `origin`. Each axis points from the corner
    /// along its round; `convex` puts the ball on the material side. Refuses
    /// parallel axes and a ball not strictly narrower than the rounds.
    pub fn new(
        origin: Point3,
        axes: [Vector3; 2],
        round: Real,
        radius: Real,
        convex: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&origin, dec)?;
        for axis in &axes {
            math::validate_vector(axis, dec)?;
            math::positive(&axis.norm_squared(), "shapes.rolling-mitre", dec)?;
        }
        math::validate(&round, dec)?;
        math::validate(&radius, dec)?;
        math::positive(&radius, "shapes.non-positive-radius", dec)?;
        let offset = &round - &radius;
        if math::sign(&offset, dec)? != Sign::Positive {
            return Err(invalid(
                "a ball as wide as the rounds does not roll over their mitre",
            ));
        }
        let axes = axes.each_ref().map(unit);
        let cosine = axes[0].dot(&axes[1]);
        let cross = axes[0].cross(&axes[1]);
        if math::sign(&cross.norm_squared(), dec)? != Sign::Positive {
            return Err(invalid("parallel round axes have no mitre"));
        }
        let up = unit(&cross);
        let out = unit(&axes[0].plus(&axes[1])).scaled(&Real::from_i64(-1));
        let side = out.cross(&up);
        let half = Real::from(Rational::from_ratio(1, 2).expect("nonzero"));
        let along_sq = fold(&(&Real::one() - &cosine) * &half);
        let across_sq = fold(&(&Real::one() + &cosine) * &half);
        let along = math::exact_sqrt(&along_sq);
        let across = math::exact_sqrt(&across_sq);
        let frame = [out, up, side];
        let shape = Shape {
            offset: &offset,
            radius: &radius,
            along: &along,
            across: &across,
            along_sq: &along_sq,
            across_sq: &across_sq,
        };
        let carrier = carrier(&origin, &frame, shape, dec)?;
        Ok(Self {
            data: Arc::new(Data {
                origin,
                axes,
                frame,
                round,
                radius,
                offset,
                along,
                across,
                along_sq,
                across_sq,
                convex,
                reversing: false,
                carrier,
            }),
        })
    }

    /// The exact rational carrier over `[0, 1] x [0, 1]`.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }
    /// Where the round axes meet.
    pub fn origin(&self) -> &Point3 {
        &self.data.origin
    }
    /// The unit round axes, each pointing from the corner along its round.
    pub fn axes(&self) -> &[Vector3; 2] {
        &self.data.axes
    }
    /// The frame `x`, `y`, `z` of the law.
    pub fn frame(&self) -> &[Vector3; 3] {
        &self.data.frame
    }
    /// The rounds' common radius.
    pub fn round(&self) -> &Real {
        &self.data.round
    }
    /// The ball's radius.
    pub fn radius(&self) -> &Real {
        &self.data.radius
    }
    /// `sqrt((1 - a1 . a2) / 2)` for the unit round axes `a1`, `a2`.
    pub fn along(&self) -> &Real {
        &self.data.along
    }
    /// `sqrt((1 + a1 . a2) / 2)` for the unit round axes `a1`, `a2`.
    pub fn across(&self) -> &Real {
        &self.data.across
    }
    /// Whether the ball is on the material side of both rounds.
    pub fn convex(&self) -> bool {
        self.data.convex
    }
    /// Whether a reflection transported the band. The chart's normal
    /// `S_u x S_v` then points toward the ball instead of away from it.
    pub fn reversing(&self) -> bool {
        self.data.reversing
    }
    fn local(&self, x: &Real, y: &Real, z: &Real) -> Point3 {
        let [a, b, c] = &self.data.frame;
        self.data
            .origin
            .plus(&a.scaled(x))
            .plus(&b.scaled(y))
            .plus(&c.scaled(z))
    }
    /// The ball centre at chart parameter `u = tan(t/2)`.
    pub fn centre(&self, u: &Rational) -> Point3 {
        let d = &self.data;
        let tau = Real::from(u.clone());
        let wide = &Real::one() + &tau.square();
        let cos = &(&Real::one() - &tau.square()) / &wide;
        let sin = &(&Real::from_i64(2) * &tau) / &wide;
        self.local(
            &(&(&d.offset * &cos) / &d.along),
            &(&d.offset * &sin),
            &Real::zero(),
        )
    }
    /// The exact point at chart parameters `u = tan(t/2)` and `v`, both in
    /// `[0, 1]`, from the law rather than the carrier.
    pub fn point(&self, u: &Rational, v: &Rational) -> Point3 {
        let law = chart(self.shape());
        let tau = Real::from(u.clone());
        let v = Real::from(v.clone());
        let powers = [Real::one(), v.clone(), v.square()];
        let weight = (0..3).fold(Real::zero(), |sum, i| {
            &sum + &(&powers[i] * &eval(&law.denominator[i], &tau))
        });
        let at = |p: &[Terms; 3]| {
            (0..3).fold(Real::zero(), |sum, i| {
                p[i].iter().fold(sum, |sum, t| {
                    &sum + &(&(&powers[i] * &eval(&t.base, &tau)) * &t.factor)
                })
            })
        };
        let [x, y, z] = law.numerators.each_ref().map(|n| &at(n) / &weight);
        self.local(&x, &y, &z)
    }
    fn shape(&self) -> Shape<'_> {
        let d = &self.data;
        Shape {
            offset: &d.offset,
            radius: &d.radius,
            along: &d.along,
            across: &d.across,
            along_sq: &d.along_sq,
            across_sq: &d.across_sq,
        }
    }
    /// The contact with round `i`, an ellipse centred where the axes meet:
    /// `O + cos t U + sin t R y` with `U = (k/sa + rho sa) x -+ rho sc z`. It
    /// is seamed at the end section's foot `O + U`, and its natural direction
    /// runs from there to the apex `O + R y` a quarter turn on.
    pub fn rail(&self, i: usize, dec: &mut Decider) -> Result<Ellipse, Refusal> {
        let d = &self.data;
        let [out, up, side] = &d.frame;
        let flip = if i == 0 { -1 } else { 1 };
        let major = out
            .scaled(&(&(&d.offset / &d.along) + &(&d.radius * &d.along)))
            .plus(&side.scaled(&(&Real::from_i64(flip) * &(&d.radius * &d.across))));
        Ellipse::new(
            d.origin.clone(),
            major.cross(up),
            major.clone(),
            major.norm_squared(),
            d.round.square(),
            dec,
        )
    }
    /// The circle of the band's end section at `t = 0`, seamed at the
    /// contact with the first round. Its natural direction from there to the
    /// second round's contact runs through the outward side: the band's end.
    pub fn section(&self, dec: &mut Decider) -> Result<Circle, Refusal> {
        let d = &self.data;
        let [out, _, side] = &d.frame;
        // `-y = x x z` in a right-handed frame. The cross product of the
        // transported axes keeps that turn after a reflection, where `-y`
        // alone would run the circle the other way round.
        Circle::new(
            self.local(&(&d.offset / &d.along), &Real::zero(), &Real::zero()),
            out.cross(side),
            d.radius.square(),
            out.scaled(&d.along).minus(&side.scaled(&d.across)),
            dec,
        )
    }
    /// The point where both contacts reach the mitre's end, `O + R y`.
    pub fn apex(&self) -> Point3 {
        self.data
            .origin
            .plus(&self.data.frame[1].scaled(&self.data.round))
    }
    /// A uniform display grid with a whole-cell chord certificate.
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<EnvelopeGrid, Refusal> {
        crate::rolling_pipe::grid::fit_net(
            self.carrier(),
            1,
            |u, v| Ok(self.point(u, v)),
            tolerance,
            max_points,
            dec,
        )
    }
    /// Rational bounds on one coordinate over the whole band. The centre's
    /// coordinate is `o + a cos t + b sin t` over the quarter turn, which
    /// stays between its ends and, where the turn reaches it, the amplitude
    /// `sqrt(a^2 + b^2)`. Every band point is within the ball's radius of
    /// its centre.
    pub fn axis_range(&self, axis: usize, dec: &Decider) -> Result<(Rational, Rational), Refusal> {
        const BITS: u32 = 64;
        let d = &self.data;
        let pick = |v: &Vector3| [&v.x, &v.y, &v.z][axis].clone();
        let a = &(&d.offset / &d.along) * &pick(&d.frame[0]);
        let b = &d.offset * &pick(&d.frame[1]);
        let amplitude = (&a.square() + &b.square()).sqrt();
        let origin = [&d.origin.x, &d.origin.y, &d.origin.z][axis].clone();
        let [a, b, amplitude, origin, radius] =
            [a, b, amplitude, origin, d.radius.clone()].map(|value| dec.enclose_at(&value, BITS));
        let (a, b, amplitude) = (a?, b?, amplitude?);
        let (origin, radius) = (origin?, radius?);
        let zero = Rational::zero();
        let mut low = a.lo().min(b.lo());
        let mut high = a.hi().max(b.hi());
        if a.lo() < zero && b.lo() < zero {
            low = low.min(-amplitude.hi());
        }
        if a.hi() > zero && b.hi() > zero {
            high = high.max(amplitude.hi());
        }
        Ok((
            &(&origin.lo() + &low) - &radius.hi(),
            &(&origin.hi() + &high) + &radius.hi(),
        ))
    }
    /// Every control and weight in world coordinates, enclosed in `u`-major order.
    pub fn enclosed_controls(&self, dec: &Decider) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "no rolling mitre control precision",
            ));
        }
        self.carrier()
            .controls()
            .iter()
            .zip(self.carrier().weights())
            .map(|(p, w)| Ok((math::point(p, bits, dec)?, dec.enclose_at(w, bits)?)))
            .collect()
    }
    /// Exact transport of the retained law and carrier.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }
    pub(crate) fn mapped(&self, map: &Isometry) -> Self {
        let d = &self.data;
        Self {
            data: Arc::new(Data {
                origin: map.apply_point(&d.origin),
                axes: d.axes.each_ref().map(|v| map.apply_vector(v)),
                frame: d.frame.each_ref().map(|v| map.apply_vector(v)),
                round: d.round.clone(),
                radius: d.radius.clone(),
                offset: d.offset.clone(),
                along: d.along.clone(),
                across: d.across.clone(),
                along_sq: d.along_sq.clone(),
                across_sq: d.across_sq.clone(),
                convex: d.convex,
                reversing: d.reversing != map.reversing(),
                carrier: d.carrier.mapped(map),
            }),
        }
    }
    /// Positive retained-definition identity for geometric query caches.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
}

fn carrier(
    origin: &Point3,
    frame: &[Vector3; 3],
    shape: Shape,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    const DEGREE: usize = 6;
    let law = chart(shape);
    let numerators = law.numerators.each_ref().map(|n| {
        across_terms(n)
            .each_ref()
            .map(|p| bernstein_terms(p, DEGREE))
    });
    let denominator = across_v(&law.denominator)
        .each_ref()
        .map(|p| bernstein(p, DEGREE));
    let mut controls = Vec::with_capacity((DEGREE + 1) * 3);
    let mut weights = Vec::with_capacity((DEGREE + 1) * 3);
    for i in 0..=DEGREE {
        for j in 0..3 {
            let weight = denominator[j][i].clone();
            if math::sign(&weight, dec)? != Sign::Positive {
                return Err(invalid("a rolling mitre weight is not positive"));
            }
            let [x, y, z] = numerators.each_ref().map(|n| &n[j][i] / &weight);
            controls.push(
                origin
                    .plus(&frame[0].scaled(&x))
                    .plus(&frame[1].scaled(&y))
                    .plus(&frame[2].scaled(&z)),
            );
            weights.push(weight);
        }
    }
    let ends = |n: usize| {
        let mut knots = vec![Rational::zero(); n + 1];
        knots.extend(vec![Rational::one(); n + 1]);
        knots
    };
    BSplineSurface::new(
        [DEGREE, 2],
        [ends(DEGREE), ends(2)],
        [DEGREE + 1, 3],
        controls,
        weights,
        dec,
    )
}
