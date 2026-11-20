//! Bend and twist (roadmap 3.4): a body sent through a space warp.
//!
//! The warps are Barr's global deformations (A. H. Barr, "Global and local
//! deformations of solid primitives", SIGGRAPH 1984), in the frames SMLib
//! gives them: `SmTwistVolume` turns each plane across the axis about it in
//! proportion to its distance along the axis, and `SmBendVolume` maps the
//! planes parallel to a neutral plane onto cylinders about the bend axis, so
//! that the neutral plane keeps its lengths. RECIPE.md section 4.20 defines
//! both as recipe operations.
//!
//! SMLib keeps a warped surface procedural (`SmSrfInVolume`: the source
//! surface composed with the volume map, evaluated through the chain rule)
//! and converts it to a B-spline only where an exact one exists. Catallax has
//! no procedural carrier, so it follows OpenCASCADE's conversion of a
//! procedural surface instead (`GeomConvert_ApproxSurface`): every face's
//! warped surface and every edge's warped curve is fitted in binary64 as a
//! cubic B-spline over the source's own parameters. A face therefore keeps
//! its source u/v chart, and every trim of the source stays a valid face
//! copy of the warped edge at the same edge parameter. The topology, keys and
//! face origins carry over unchanged. Numerical admission then measures every
//! edge gap, as for any other fitted body (rule T4); nothing here is a proof.

use crate::numerical_embedding::nurbs;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal};
use catallax_recipe::{Axis, Number, Operation};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Curve, Point3, RationalAtlas, SplineSupport, Surface,
};
use catallax_topology::Sense;
use catallax_topology::numerical::{
    CopyTag, NumericalDefinition, NumericalEdge, NumericalFace, NumericalNativeTrim, NumericalTrim,
    NumericalVertex,
};
use std::collections::BTreeSet;

const BITS: u32 = 128;
type V = [f64; 3];

fn refusal(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

fn bad(detail: impl Into<String>) -> Refusal {
    refusal("deform.definition", detail)
}

fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: V, s: f64) -> V {
    a.map(|x| x * s)
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn distance(a: V, b: V) -> f64 {
    let d = sub(a, b);
    dot(d, d).sqrt()
}

fn number(n: &Number) -> f64 {
    n.value().to_f64_nearest()
}

fn unit(v: V, what: &str) -> Result<V, Refusal> {
    let length = dot(v, v).sqrt();
    if !length.is_finite() || length == 0.0 {
        return Err(bad(format!("the {what} is zero or not finite")));
    }
    Ok(scale(v, 1.0 / length))
}

/// The unit direction of an exact vector. Dividing by its largest component
/// first puts every component in [-1, 1], one of them exactly 1 or -1, so no
/// component overflows or underflows in binary64 and the direction survives
/// rounding to within a unit in the last place.
fn exact_direction(v: &[Rational; 3], what: &str) -> Result<V, Refusal> {
    let largest =
        v.iter().map(Rational::abs).fold(
            Rational::zero(),
            |largest, x| {
                if x > largest { x } else { largest }
            },
        );
    let inverse = largest
        .recip()
        .ok_or_else(|| bad(format!("the {what} is zero")))?;
    unit(v.each_ref().map(|x| (x * &inverse).to_f64_nearest()), what)
}

fn exact(value: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(value)
        .ok_or_else(|| refusal("deform.fit", "a warped value is not finite"))
}

fn point(p: V) -> Result<Point3, Refusal> {
    Ok(Point3::new(
        exact(p[0])?.into(),
        exact(p[1])?.into(),
        exact(p[2])?.into(),
    ))
}

fn f64_of(value: &Real, dec: &Decider) -> Result<f64, Refusal> {
    Ok(dec.enclose_at(value, BITS)?.mid().to_f64_nearest())
}

fn ball_point(p: &[Ball; 3]) -> V {
    p.each_ref().map(|v| v.mid().to_f64_nearest())
}

/// The largest binary64 value not above `value`.
fn below(value: &Rational) -> f64 {
    let x = value.to_f64_nearest();
    if Rational::from_f64_exact(x).is_some_and(|q| &q > value) {
        x.next_down()
    } else {
        x
    }
}

/// The smallest binary64 value not below `value`.
fn above(value: &Rational) -> f64 {
    let x = value.to_f64_nearest();
    if Rational::from_f64_exact(x).is_some_and(|q| &q < value) {
        x.next_up()
    } else {
        x
    }
}

/// A space warp in binary64. Both keep orientation: the twist's Jacobian
/// determinant is one, the bend's is `u / radius`, positive on its domain.
#[derive(Clone, Debug)]
pub(crate) enum Warp {
    /// `SmTwistVolume`: a point at signed distance `w` along the axis from
    /// `origin` turns by `rate * w` radians about the axis.
    Twist { origin: V, axis: V, rate: f64 },
    /// `SmBendVolume`: with `u`, `v`, `w` a point's coordinates along `mid`,
    /// `binormal` and `axis` from `origin`, the point goes to the circle of
    /// radius `u` about the axis at angle `v / radius` from `mid`.
    Bend {
        origin: V,
        axis: V,
        mid: V,
        binormal: V,
        radius: f64,
    },
}

impl Warp {
    /// The warp a bend or twist step asks for, or `None` for any other step.
    pub(crate) fn of(op: &Operation) -> Result<Option<Self>, Refusal> {
        Ok(Some(match op {
            Operation::Twist {
                axis: Axis { point, direction },
                angle,
                length,
                ..
            } => {
                if length.value().sign() != catallax_numbers::Sign::Positive {
                    return Err(bad("a twist needs a positive length"));
                }
                let rate = number(angle).to_radians() / number(length);
                if !rate.is_finite() {
                    return Err(bad("the twist rate is not finite"));
                }
                Self::Twist {
                    origin: point.0.each_ref().map(number),
                    axis: exact_direction(
                        &direction.0.each_ref().map(|n| n.value().clone()),
                        "twist axis",
                    )?,
                    rate,
                }
            }
            Operation::Bend {
                axis: Axis { point, direction },
                toward,
                radius,
                ..
            } => {
                if radius.value().sign() != catallax_numbers::Sign::Positive {
                    return Err(bad("a bend needs a positive radius"));
                }
                // Decide "along the axis" on the exact inputs: in binary64 the
                // cross-axis part of a parallel `toward` is a rounding residue,
                // which `unit` would scale up into an arbitrary direction.
                let (a, t) = (&direction.0, &toward.0);
                let parallel = (0..3).all(|i| {
                    let (j, k) = ((i + 1) % 3, (i + 2) % 3);
                    let m = |x: &Number, y: &Number| x.value() * y.value();
                    (&m(&a[j], &t[k]) - &m(&a[k], &t[j])).is_zero()
                });
                if parallel {
                    return Err(bad(
                        "the bend direction runs along the axis (or is zero); it must point across it",
                    ));
                }
                // SMLib orients the bend by its axis and a mid direction
                // orthogonal to it; `toward` gives that direction up to its
                // component along the axis. That part is taken on the exact
                // inputs as well: |a|^2 t - (a.t) a is it times |a|^2. Rounding
                // the inputs first can leave only a residue pointing anywhere,
                // as for a = (3, 5, 7) and t = (3, 5, 7.0000000000000001),
                // whose part across is along (-21, -35, 34).
                let exact = |v: &[Number; 3]| v.each_ref().map(|n| n.value().clone());
                let (a, t) = (exact(a), exact(t));
                let inner = |x: &[Rational; 3], y: &[Rational; 3]| {
                    &(&(&x[0] * &y[0]) + &(&x[1] * &y[1])) + &(&x[2] * &y[2])
                };
                let (aa, at) = (inner(&a, &a), inner(&a, &t));
                let across: [Rational; 3] =
                    std::array::from_fn(|i| &(&t[i] * &aa) - &(&a[i] * &at));
                let axis = exact_direction(&a, "bend axis")?;
                let mid = exact_direction(&across, "bend direction across the axis")?;
                // One Gram-Schmidt step removes the rounding left between the two.
                let mid = unit(
                    sub(mid, scale(axis, dot(mid, axis))),
                    "bend direction across the axis",
                )?;
                Self::Bend {
                    origin: point.0.each_ref().map(number),
                    axis,
                    mid,
                    binormal: cross(axis, mid),
                    radius: number(radius),
                }
            }
            _ => return Ok(None),
        }))
    }

    /// The warped position of `p`.
    pub(crate) fn map(&self, p: V) -> V {
        match *self {
            Self::Twist { origin, axis, rate } => {
                let d = sub(p, origin);
                let w = dot(d, axis);
                let across = sub(d, scale(axis, w));
                let (sin, cos) = (rate * w).sin_cos();
                // Rodrigues' rotation of a vector orthogonal to the axis.
                add(
                    add(origin, scale(axis, w)),
                    add(scale(across, cos), scale(cross(axis, across), sin)),
                )
            }
            Self::Bend {
                origin,
                axis,
                mid,
                binormal,
                radius,
            } => {
                let d = sub(p, origin);
                let [u, v, w] = [dot(d, mid), dot(d, binormal), dot(d, axis)];
                let (sin, cos) = (v / radius).sin_cos();
                add(
                    origin,
                    add(
                        add(scale(mid, u * cos), scale(binormal, u * sin)),
                        scale(axis, w),
                    ),
                )
            }
        }
    }

    /// SMLib's bend domain: the half-space `u > 0` on the neutral plane's
    /// side of the axis, and less than one turn along the binormal, beyond
    /// which the map folds the body onto itself (`SmBendVolume.h`, notes 1
    /// and 2). A twist is one-to-one everywhere. On a polyhedral body both
    /// coordinates are linear, so the vertices hold their extremes.
    fn check_domain(&self, vertices: &[V]) -> Result<(), Refusal> {
        let Self::Bend {
            origin,
            mid,
            binormal,
            radius,
            ..
        } = *self
        else {
            return Ok(());
        };
        let mut v_range = [f64::INFINITY, f64::NEG_INFINITY];
        for p in vertices {
            let d = sub(*p, origin);
            let u = dot(d, mid);
            if u.is_nan() || u <= 0.0 {
                return Err(refusal(
                    "deform.bend-axis",
                    format!("the body reaches the bend axis or crosses it (distance {u:e} mm)"),
                ));
            }
            let v = dot(d, binormal);
            v_range = [v_range[0].min(v), v_range[1].max(v)];
        }
        let turn = (v_range[1] - v_range[0]) / radius;
        if !turn.is_finite() || turn >= std::f64::consts::TAU {
            return Err(refusal(
                "deform.bend-overlap",
                format!("the bent body would wrap {turn:.6} radians, a full turn or more"),
            ));
        }
        Ok(())
    }
}

/// A planar face chart in binary64: `origin + u x + v y`.
struct PlaneChart {
    origin: V,
    x: V,
    y: V,
}

impl PlaneChart {
    fn new(surface: &Surface, dec: &mut Decider) -> Result<Self, Refusal> {
        let zero = Ball::from_rational(&Rational::zero(), BITS);
        let origin = ball_point(&surface.evaluate(&zero, &zero, dec)?);
        let [x, y] = surface.derivatives(&zero, &zero, dec)?;
        Ok(Self {
            origin,
            x: ball_point(&x),
            y: ball_point(&y),
        })
    }
    fn point(&self, uv: [f64; 2]) -> V {
        add(self.origin, add(scale(self.x, uv[0]), scale(self.y, uv[1])))
    }
}

/// A straight edge in binary64: `origin + t dir`.
struct Line64 {
    origin: V,
    dir: V,
}

/// The spans of each fit, doubled from the first until the fit meets its
/// goal. These are the numerical embedding's sample densities: four spans
/// to start, and at most 256, because admission checks every span of every
/// loop pairwise (`numerical_embedding::DENSITIES` and `MAX_SPANS`).
const FIRST_SPANS: usize = 4;
const MAX_SPANS: usize = 256;
/// Each fit is checked at the eighths of every span, as the numerical
/// embedding's fits are; admission checks at its own, different locations.
const CHECKS_PER_SPAN: usize = 8;

fn nodes(range: [f64; 2], spans: usize) -> Vec<f64> {
    (0..=spans)
        .map(|i| {
            if i == spans {
                range[1]
            } else {
                range[0] + (range[1] - range[0]) * i as f64 / spans as f64
            }
        })
        .collect()
}

fn checks(range: [f64; 2], spans: usize) -> Vec<f64> {
    nodes(range, spans * CHECKS_PER_SPAN)
}

/// Cubic not-a-knot interpolation of `f` at equally spaced nodes, the
/// spans doubled until every check is within `goal`. Returns the knots, the
/// controls and the largest checked error.
fn fit_curve(
    range: [f64; 2],
    f: &impl Fn(f64) -> V,
    goal: f64,
    dec: &mut Decider,
) -> Result<(Vec<f64>, Vec<V>, f64), Refusal> {
    let mut spans = FIRST_SPANS;
    loop {
        let parameters = nodes(range, spans);
        dec.charge_work(parameters.len() as u64)?;
        let values: Vec<V> = parameters.iter().map(|t| f(*t)).collect();
        let (knots, controls) = nurbs::interpolate(&parameters, &values, dec)?;
        let mut error = 0.0_f64;
        let locations = checks(range, spans);
        dec.charge_work(locations.len() as u64)?;
        for t in locations {
            error = error.max(distance(nurbs::evaluate(&knots, &controls, t), f(t)));
        }
        if !error.is_finite() {
            return Err(refusal("deform.fit", "a warped edge is not finite"));
        }
        if error <= goal {
            return Ok((knots, controls, error));
        }
        spans *= 2;
        if spans > MAX_SPANS {
            return Err(refusal(
                "deform.fit",
                format!(
                    "a warped edge needs more than {MAX_SPANS} spans to come within {goal:e} mm; no tolerance was enlarged"
                ),
            ));
        }
    }
}

/// A tensor-product cubic B-spline in binary64, u-major controls.
struct Tensor {
    knots: [Vec<f64>; 2],
    controls: Vec<Vec<V>>,
}

impl Tensor {
    fn point(&self, uv: [f64; 2]) -> V {
        let column: Vec<V> = self
            .controls
            .iter()
            .map(|row| nurbs::evaluate(&self.knots[1], row, uv[1]))
            .collect();
        nurbs::evaluate(&self.knots[0], &column, uv[0])
    }
}

/// Tensor-product interpolation (The NURBS Book, A9.4): interpolate every
/// row along u, then every control column along v. The spans of each
/// direction double until its own error, measured where the other
/// direction sits on its nodes, is within half the goal, and the whole grid
/// of checks is within the goal.
fn fit_surface(
    region: [[f64; 2]; 2],
    f: &impl Fn([f64; 2]) -> V,
    goal: f64,
    dec: &mut Decider,
) -> Result<(Tensor, f64), Refusal> {
    let mut spans = [FIRST_SPANS; 2];
    loop {
        let u = nodes(region[0], spans[0]);
        let v = nodes(region[1], spans[1]);
        dec.charge_work((u.len() * v.len()) as u64)?;
        let mut columns = vec![Vec::with_capacity(u.len()); v.len()];
        for uu in &u {
            for (j, vv) in v.iter().enumerate() {
                columns[j].push(f([*uu, *vv]));
            }
        }
        // Along u for each v node: the controls of each v row.
        let mut knots_u = Vec::new();
        let mut rows: Vec<Vec<V>> = vec![Vec::with_capacity(v.len()); u.len()];
        for column in &columns {
            let (knots, controls) = nurbs::interpolate(&u, column, dec)?;
            knots_u = knots;
            for (i, control) in controls.into_iter().enumerate() {
                rows[i].push(control);
            }
        }
        // Along v for each u control.
        let mut knots_v = Vec::new();
        let mut controls = Vec::with_capacity(rows.len());
        for row in &rows {
            let (knots, fitted) = nurbs::interpolate(&v, row, dec)?;
            knots_v = knots;
            controls.push(fitted);
        }
        let tensor = Tensor {
            knots: [knots_u, knots_v],
            controls,
        };
        let cu = checks(region[0], spans[0]);
        let cv = checks(region[1], spans[1]);
        dec.charge_work((cu.len() * cv.len()) as u64)?;
        let mut along = [0.0_f64; 2];
        let mut whole = 0.0_f64;
        for (i, uu) in cu.iter().enumerate() {
            for (j, vv) in cv.iter().enumerate() {
                let error = distance(tensor.point([*uu, *vv]), f([*uu, *vv]));
                whole = whole.max(error);
                if j % CHECKS_PER_SPAN == 0 {
                    along[0] = along[0].max(error);
                }
                if i % CHECKS_PER_SPAN == 0 {
                    along[1] = along[1].max(error);
                }
            }
        }
        if !whole.is_finite() {
            return Err(refusal("deform.fit", "a warped face is not finite"));
        }
        if whole <= goal {
            return Ok((tensor, whole));
        }
        let mut refine = along.map(|error| error > goal / 2.0);
        if refine == [false; 2] {
            refine = [true; 2];
        }
        for axis in 0..2 {
            if refine[axis] {
                spans[axis] *= 2;
            }
        }
        if spans.iter().any(|s| *s > MAX_SPANS) {
            return Err(refusal(
                "deform.fit",
                format!(
                    "a warped face needs more than {MAX_SPANS} spans a side to come within {goal:e} mm; no tolerance was enlarged"
                ),
            ));
        }
    }
}

fn exact_knots(knots: &[f64]) -> Result<Vec<Rational>, Refusal> {
    knots.iter().map(|k| exact(*k)).collect()
}

fn curve(knots: &[f64], controls: &[V], dec: &mut Decider) -> Result<BSplineCurve, Refusal> {
    BSplineCurve::new(
        knots.len() - controls.len() - 1,
        exact_knots(knots)?,
        controls
            .iter()
            .map(|p| point(*p))
            .collect::<Result<_, _>>()?,
        vec![Real::one(); controls.len()],
        dec,
    )
}

fn surface(tensor: &Tensor, dec: &mut Decider) -> Result<Surface, Refusal> {
    let count = [tensor.controls.len(), tensor.controls[0].len()];
    let degrees = std::array::from_fn(|axis| tensor.knots[axis].len() - count[axis] - 1);
    let carrier = BSplineSurface::new(
        degrees,
        [
            exact_knots(&tensor.knots[0])?,
            exact_knots(&tensor.knots[1])?,
        ],
        count,
        tensor
            .controls
            .iter()
            .flatten()
            .map(|p| point(*p))
            .collect::<Result<_, _>>()?,
        vec![Real::one(); count[0] * count[1]],
        dec,
    )?;
    // The fitted rolling-ball strips are admitted the same way: exact span
    // coefficients, a proved regular chart cover, and one logical face.
    let [atlas] = <[_; 1]>::try_from(RationalAtlas::partition_with_polynomial_bounds(
        carrier,
        [false, false],
        dec,
    )?)
    .map_err(|_| {
        refusal(
            "deform.fit",
            "a warped face does not form one regular logical face",
        )
    })?;
    Ok(Surface::BSpline(SplineSupport::RationalAtlas(atlas)))
}

/// The u/v line of a source trim at the edge's own parameter: the face copy
/// the source admission measured, written as a native curve because an
/// atlas face needs one on every use.
fn uv_line(
    edge: &NumericalEdge,
    trim: &NumericalTrim,
    dec: &mut Decider,
) -> Result<BSplineCurve, Refusal> {
    let [first, second] = &trim.uv;
    let (at_start, at_end) = match trim.sense {
        Sense::Same => (first, second),
        Sense::Reversed => (second, first),
    };
    let [t0, t1] = &edge.parameters;
    let (lo, hi, a, b) = if t0 < t1 {
        (t0, t1, at_start, at_end)
    } else {
        (t1, t0, at_end, at_start)
    };
    let row =
        |uv: &[Rational; 2]| Point3::new(uv[0].clone().into(), uv[1].clone().into(), Real::zero());
    BSplineCurve::new(
        1,
        vec![lo.clone(), lo.clone(), hi.clone(), hi.clone()],
        vec![row(a), row(b)],
        vec![Real::one(), Real::one()],
        dec,
    )
}

/// Warps an embedded body. Every face of this first build is planar and
/// every edge straight; curved input is coverage follow-on item 3
/// (ROADMAP.md, "Widen 3.4 from a first build to general input").
pub(crate) fn warp(
    source: &NumericalDefinition,
    native_trims: &[NumericalNativeTrim],
    warp: &Warp,
    target_mm: f64,
    dec: &mut Decider,
) -> Result<
    (
        NumericalDefinition,
        Vec<NumericalNativeTrim>,
        BTreeSet<String>,
    ),
    Refusal,
> {
    if !native_trims.is_empty() {
        return Err(refusal(
            "deform.unsupported-edge",
            "bend and twist take straight edges on planar faces in their first build",
        ));
    }
    // A quarter of the target for each fit, as the numerical embedding's
    // fits meet a quarter of the tolerance: an edge use's gap combines the
    // edge's fit and its face's fit.
    let goal = target_mm / 4.0;
    let positions: Vec<V> = source
        .vertices
        .iter()
        .map(|v| v.position.each_ref().map(Rational::to_f64_nearest))
        .collect();
    warp.check_domain(&positions)?;
    let vertices = source
        .vertices
        .iter()
        .zip(&positions)
        .map(|(v, p)| {
            let q = warp.map(*p);
            Ok(NumericalVertex {
                key: v.key.clone(),
                position: [exact(q[0])?, exact(q[1])?, exact(q[2])?],
            })
        })
        .collect::<Result<Vec<_>, Refusal>>()?;
    let mut edges = Vec::with_capacity(source.edges.len());
    for edge in &source.edges {
        dec.charge_work(1)?;
        let Curve::Line(line) = &edge.curve else {
            return Err(refusal(
                "deform.unsupported-edge",
                format!(
                    "edge {} is curved; bend and twist take straight edges in their first build",
                    edge.key
                ),
            ));
        };
        let line = Line64 {
            origin: [
                f64_of(&line.origin.x, dec)?,
                f64_of(&line.origin.y, dec)?,
                f64_of(&line.origin.z, dec)?,
            ],
            dir: [
                f64_of(&line.dir.x, dec)?,
                f64_of(&line.dir.y, dec)?,
                f64_of(&line.dir.z, dec)?,
            ],
        };
        let [t0, t1] = &edge.parameters;
        let (lo, hi) = if t0 < t1 { (t0, t1) } else { (t1, t0) };
        let range = [below(lo), above(hi)];
        let f = |t: f64| warp.map(add(line.origin, scale(line.dir, t)));
        let (knots, controls, _) = fit_curve(range, &f, goal, dec)?;
        edges.push(NumericalEdge {
            key: edge.key.clone(),
            curve: Curve::BSpline(curve(&knots, &controls, dec)?),
            vertices: edge.vertices,
            parameters: edge.parameters.clone(),
        });
    }
    let mut faces = Vec::with_capacity(source.faces.len());
    let mut native = Vec::new();
    let mut sampled = BTreeSet::new();
    for face in &source.faces {
        dec.charge_work(1)?;
        if !matches!(face.surface, Surface::Plane(_)) {
            return Err(refusal(
                "deform.unsupported-face",
                format!(
                    "face {} is curved; bend and twist take planar faces in their first build",
                    face.key
                ),
            ));
        }
        let chart = PlaneChart::new(&face.surface, dec)?;
        let mut uvs = face.loops.iter().flatten().flat_map(|trim| &trim.uv);
        let first = uvs
            .next()
            .ok_or_else(|| bad(format!("face {} has no boundary", face.key)))?;
        let mut bounds = [[&first[0], &first[0]], [&first[1], &first[1]]];
        for uv in uvs {
            for axis in 0..2 {
                bounds[axis][0] = bounds[axis][0].min(&uv[axis]);
                bounds[axis][1] = bounds[axis][1].max(&uv[axis]);
            }
        }
        if bounds.iter().any(|[lo, hi]| lo >= hi) {
            return Err(bad(format!("face {} has a degenerate boundary", face.key)));
        }
        let region = bounds.map(|[lo, hi]| [below(lo), above(hi)]);
        let f = |uv: [f64; 2]| warp.map(chart.point(uv));
        let (tensor, _) = fit_surface(region, &f, goal, dec)?;
        let mut loops = Vec::with_capacity(face.loops.len());
        for (boundary, row) in face.loops.iter().enumerate() {
            let mut trims = Vec::with_capacity(row.len());
            for (use_index, trim) in row.iter().enumerate() {
                let edge = &source.edges[trim.edge];
                native.push(NumericalNativeTrim {
                    face: face.key.clone(),
                    boundary,
                    use_index,
                    edge: edge.key.clone(),
                    curve: uv_line(edge, trim, dec)?,
                    // Both copies were refitted: admission measures again.
                    tag: CopyTag::Independent,
                });
                trims.push(NumericalTrim {
                    edge: trim.edge,
                    sense: trim.sense,
                    uv: trim.uv.clone(),
                    tag: CopyTag::Independent,
                });
            }
            loops.push(trims);
        }
        sampled.insert(face.key.clone());
        faces.push(NumericalFace {
            key: face.key.clone(),
            surface: surface(&tensor, dec)?,
            sense: face.sense,
            loops,
        });
    }
    Ok((
        NumericalDefinition {
            name: source.name.clone(),
            vertices,
            edges,
            faces,
        },
        native,
        sampled,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn twist(degrees: f64, length: f64) -> Warp {
        Warp::Twist {
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            rate: degrees.to_radians() / length,
        }
    }

    fn bend(radius: f64) -> Warp {
        Warp::Bend {
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            mid: [0.0, 1.0, 0.0],
            binormal: [-1.0, 0.0, 0.0],
            radius,
        }
    }

    fn close(a: V, b: V) -> bool {
        distance(a, b) <= 1e-12 * (1.0 + dot(b, b).sqrt())
    }

    #[test]
    fn a_twist_turns_each_cross_section_and_keeps_the_axis_plane() {
        let warp = twist(90.0, 100.0);
        assert!(close(warp.map([5.0, 0.0, 0.0]), [5.0, 0.0, 0.0]));
        assert!(close(warp.map([5.0, 0.0, 100.0]), [0.0, 5.0, 100.0]));
        assert!(close(warp.map([5.0, 0.0, -100.0]), [0.0, -5.0, -100.0]));
        // SMLib: the inverse of a twist is the twist of opposite sense.
        let back = twist(-90.0, 100.0);
        let p = [3.0, -4.0, 37.0];
        assert!(close(back.map(warp.map(p)), p));
    }

    #[test]
    fn a_bend_keeps_the_neutral_plane_lengths_and_the_mid_plane() {
        let warp = bend(20.0);
        // The plane through the axis along `mid` stays put.
        assert!(close(warp.map([0.0, 7.0, 3.0]), [0.0, 7.0, 3.0]));
        // A neutral-plane point at arc length 20 * pi / 2 turns a quarter.
        let p = warp.map([-10.0 * std::f64::consts::PI, 20.0, 1.0]);
        assert!(close(p, [-20.0, 0.0, 1.0]), "{p:?}");
        assert!(warp.check_domain(&[[0.0, 1.0, 0.0]]).is_ok());
        assert_eq!(
            warp.check_domain(&[[0.0, -1.0, 0.0]]).unwrap_err().code,
            "deform.bend-axis"
        );
        let wrap = 20.0 * std::f64::consts::TAU;
        assert_eq!(
            warp.check_domain(&[[0.0, 20.0, 0.0], [-wrap, 20.0, 0.0]])
                .unwrap_err()
                .code,
            "deform.bend-overlap"
        );
    }

    /// The part of `toward` across the axis comes from the exact inputs. In
    /// binary64 both vectors below round to (3, 5, 7), whose difference is a
    /// residue near (0, 0, -1); the exact part is along (-21, -35, 34).
    #[test]
    fn a_nearly_parallel_bend_direction_keeps_its_exact_part_across() {
        let recipe = catallax_recipe::Recipe::parse(
            r#"{"format": "catallax-recipe/5", "name": "near", "result": "bent", "steps": [
                {"id": "plate", "op": "block", "size": ["1", "1", "1"],
                 "pose": {"rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"]}},
                {"id": "bent", "op": "bend", "body": "plate",
                 "axis": {"point": ["0", "0", "0"], "direction": ["3", "5", "7"]},
                 "toward": ["3", "5", "7.0000000000000001"], "radius": "20"}]}"#,
        )
        .unwrap();
        let Some(Warp::Bend {
            axis,
            mid,
            binormal,
            ..
        }) = Warp::of(&recipe.steps[1].op).unwrap()
        else {
            panic!("a bend step gives a bend");
        };
        let length = (21.0f64 * 21.0 + 35.0 * 35.0 + 34.0 * 34.0).sqrt();
        assert!(
            close(mid, [-21.0 / length, -35.0 / length, 34.0 / length]),
            "{mid:?}"
        );
        assert!(dot(mid, axis).abs() <= 1e-15);
        assert!(close(binormal, cross(axis, mid)));
        // The two frames disagree on which side of the axis (0, 0, -1) lies.
        assert!(dot([0.0, 0.0, -1.0], mid) < 0.0);
    }

    #[test]
    fn fits_double_their_spans_until_the_checks_pass() {
        let warp = twist(90.0, 100.0);
        let f = |t: f64| warp.map([5.0, 5.0, t]);
        let (knots, controls, error) = catallax_numbers::with_operation_budget(1_000_000, || {
            fit_curve([0.0, 100.0], &f, 1e-6, &mut Decider::default())
        })
        .0
        .unwrap();
        assert!(error <= 1e-6);
        // Not-a-knot interpolation has one control per node; the helix needed more.
        assert!(controls.len() > FIRST_SPANS + 1, "{}", controls.len());
        for t in [0.0, 12.5, 63.0, 100.0] {
            assert!(distance(nurbs::evaluate(&knots, &controls, t), f(t)) <= 1e-6);
        }
        let g = |uv: [f64; 2]| warp.map([uv[0], 5.0, uv[1]]);
        let (tensor, error) = catallax_numbers::with_operation_budget(10_000_000, || {
            fit_surface(
                [[-5.0, 5.0], [0.0, 100.0]],
                &g,
                1e-6,
                &mut Decider::default(),
            )
        })
        .0
        .unwrap();
        // The twist is linear across the axis, so u keeps its first spans.
        assert_eq!(tensor.controls.len(), FIRST_SPANS + 1);
        assert!(tensor.controls[0].len() > FIRST_SPANS + 1);
        assert!(error <= 1e-6);
    }
}
