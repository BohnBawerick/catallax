//! Section-law sweeps admitted and exported only as numerical B-reps.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{Loop, Number, Point3 as RecipePoint, Pose as RecipePose};
use catallax_shapes::{
    BSplineCurve, BSplineSurface, Curve, Plane, Point3, RationalPatch, SplineBoundary, Surface,
    Vector3,
};
use catallax_topology::{
    Sense,
    numerical::{
        self, NumericalBody, NumericalDefinition, NumericalEdge, NumericalFace, NumericalOptions,
        NumericalTrim, NumericalVertex,
    },
};
use serde::{Deserialize, Serialize};

mod recipe;
#[cfg(test)]
mod tests;
pub(crate) use recipe::{Context, path_needs_numerical};
pub use recipe::{RecipeCandidate, RecipeSwept, prepare_recipe, requires_numerical};
/// Persisted generating definition. Reopening rebuilds and rechecks the geometry.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SweepRecipe {
    /// Stable operation identity.
    pub name: String,
    /// Section boundaries in local XY. Nesting follows ordinary sketch semantics.
    pub loops: Vec<Loop>,
    /// Authored world-space sketch plane. Omission attaches local XY at the spine start.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plane: Option<RecipePose>,
    /// Spine, parameterized from zero to one.
    pub path: SweepPath,
    /// Fixed translation or rotation-minimising transport. A guide supplies its own frame.
    #[serde(default = "rotation_minimising")]
    pub orientation: catallax_recipe::SweepOrientation,
    /// Auxiliary spine. Correspondence uses equal normalized parameters, not arc length.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub guide: Option<SweepPath>,
    /// Initial section X direction, projected perpendicular to the tangent.
    pub reference: [Number; 3],
    /// Additional linear rotation about the moving tangent, in degrees.
    pub twist_degrees: Number,
    /// Positive final homothetic scale about the spine. This is not a draft angle.
    pub end_scale: Number,
    /// Admission and sampled fitting tolerance in millimetres, fixed before building.
    pub tolerance_mm: Number,
    /// Maximum number of cubic sections the caller permits.
    pub max_sections: usize,
}

/// Supported analytic and rational spine definitions.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum SweepPath {
    /// A directed line segment.
    Line {
        /// Initial spine point.
        start: RecipePoint,
        /// Terminal spine point.
        end: RecipePoint,
    },
    /// The ordinary recipe helix definition, with its starting point retained.
    Helix {
        /// Initial spine point.
        start: RecipePoint,
        /// Axis point in the starting radial plane.
        centre: RecipePoint,
        /// Nonzero axis direction.
        axis: catallax_recipe::Vector3,
        /// Signed rotation in degrees.
        angle: Number,
        /// Signed displacement along the unit axis.
        rise: Number,
    },
    /// A clamped, positive-weight rational B-spline.
    Spline {
        /// Polynomial degree.
        degree: usize,
        /// Clamped curve controls.
        control_points: Vec<RecipePoint>,
        /// One positive weight per control.
        weights: Vec<Number>,
        /// Expanded clamped knot vector.
        knots: Vec<Number>,
    },
}

/// Pointwise fitting checks. These do not claim a continuous or whole-body bound.
#[derive(Clone, Debug, Serialize)]
pub struct SweepChecks {
    /// Accepted cubic section count.
    pub sections: usize,
    /// Number of checked profile and section sample locations.
    pub checked_points: usize,
    /// Largest section-law residual at a checked binary64 sample.
    pub max_sample_error_mm: f64,
    /// Maximum source-profile sample displacement after numerical scalar fitting.
    pub max_profile_sample_error_mm: f64,
}

/// Checked numerical geometry and the definition needed for a parameter edit.
#[derive(Clone, Debug)]
pub struct Swept {
    recipe: SweepRecipe,
    body: NumericalBody,
    checks: SweepChecks,
}

impl Swept {
    /// Retained generating definition.
    pub fn recipe(&self) -> &SweepRecipe {
        &self.recipe
    }
    /// Numerically admitted B-rep and its admission evidence.
    pub fn body(&self) -> &NumericalBody {
        &self.body
    }
    /// Pointwise fitting diagnostics, separate from admission checks.
    pub fn checks(&self) -> &SweepChecks {
        &self.checks
    }
    /// Moves the admitted body to a numerical consumer without copying its tables.
    pub fn into_body(self) -> NumericalBody {
        self.body
    }

    /// Rebuilds transactionally. Previous admission evidence remains historical.
    pub fn edit(&self, recipe: SweepRecipe, dec: &mut Decider) -> Result<Self, Refusal> {
        let (definition, checks, _) = definition(&recipe, dec)?;
        let body = self.body.rebuild(definition, options(&recipe), dec)?;
        Ok(Self {
            recipe,
            body,
            checks,
        })
    }
}

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new("sweep.numerical-definition", detail)
}
fn q(x: i64) -> Rational {
    Rational::from_i64(x)
}
fn rational(x: f64) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(x).ok_or_else(|| bad("nonfinite sweep coordinate"))
}
fn real(x: f64) -> Result<Real, Refusal> {
    Ok(rational(x)?.into())
}
type V = [f64; 3];
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn mul(a: V, b: f64) -> V {
    a.map(|x| x * b)
}
fn dot(a: V, b: V) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn unit(a: V) -> Result<V, Refusal> {
    let length = libm::sqrt(dot(a, a));
    if !length.is_finite() || length <= 0.0 {
        return Err(bad("zero or nonfinite frame direction"));
    }
    Ok(mul(a, 1.0 / length))
}
fn number(n: &Number) -> f64 {
    n.value().to_f64_nearest()
}
fn recipe_point(p: &RecipePoint) -> V {
    p.0.each_ref().map(number)
}
fn point(p: V) -> Result<Point3, Refusal> {
    Ok(Point3::new(real(p[0])?, real(p[1])?, real(p[2])?))
}
fn vector(p: V) -> Result<Vector3, Refusal> {
    Ok(Vector3::new(real(p[0])?, real(p[1])?, real(p[2])?))
}
fn flip(s: Sense) -> Sense {
    crate::extrude::reverse(s)
}
fn rotation_minimising() -> catallax_recipe::SweepOrientation {
    catallax_recipe::SweepOrientation::RotationMinimising
}

fn options(recipe: &SweepRecipe) -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: recipe.tolerance_mm.value().clone(),
    }
}

#[derive(Clone)]
enum Path {
    Line(V, V),
    Helix {
        radius: f64,
        pitch: f64,
        turns: f64,
        origin: V,
        axes: [V; 3],
    },
    Spline(BSplineCurve),
}
impl Path {
    fn new(path: &SweepPath, dec: &mut Decider) -> Result<Self, Refusal> {
        match path {
            SweepPath::Line { start, end } => {
                let a = recipe_point(start);
                let delta = sub(recipe_point(end), a);
                unit(delta)?;
                Ok(Self::Line(a, delta))
            }
            SweepPath::Helix {
                start,
                centre,
                axis,
                angle,
                rise,
            } => {
                let as_point = |p: &RecipePoint| {
                    catallax_shapes::RPoint3::new(
                        p.0[0].value().clone(),
                        p.0[1].value().clone(),
                        p.0[2].value().clone(),
                    )
                };
                let a = as_point(start);
                let c = as_point(centre);
                let axis_r = catallax_shapes::RVector3::new(
                    axis.0[0].value().clone(),
                    axis.0[1].value().clone(),
                    axis.0[2].value().clone(),
                );
                let radius_r = a.minus(&c);
                if axis_r.is_zero()
                    || radius_r.is_zero()
                    || angle.value().is_zero()
                    || !axis_r.dot(&radius_r).is_zero()
                {
                    return Err(bad(
                        "a helix needs nonzero axis, radius and angle, with its start in the centre plane",
                    ));
                }
                let z = unit(axis.0.each_ref().map(number))?;
                let radial = sub(recipe_point(start), recipe_point(centre));
                let radius = libm::sqrt(dot(radial, radial));
                let x = unit(radial)?;
                let y = unit(cross(z, x))?;
                let turns = number(angle) / 360.0;
                Ok(Self::Helix {
                    radius,
                    pitch: number(rise) / turns,
                    turns,
                    origin: recipe_point(centre),
                    axes: [x, y, z],
                })
            }
            SweepPath::Spline {
                degree,
                control_points,
                weights,
                knots,
            } => {
                let points = control_points
                    .iter()
                    .map(|p| {
                        let [x, y, z] = p.0.each_ref().map(|n| Real::from(n.value().clone()));
                        Point3::new(x, y, z)
                    })
                    .collect();
                Ok(Self::Spline(BSplineCurve::new(
                    *degree,
                    knots.iter().map(|k| k.value().clone()).collect(),
                    points,
                    weights
                        .iter()
                        .map(|w| Real::from(w.value().clone()))
                        .collect(),
                    dec,
                )?))
            }
        }
    }
    fn at(&self, t: f64, dec: &mut Decider) -> Result<(V, V), Refusal> {
        dec.charge_work(1)?;
        match self {
            Self::Line(a, delta) => Ok((add(*a, mul(*delta, t)), *delta)),
            Self::Helix {
                radius,
                pitch,
                turns,
                origin,
                axes,
            } => {
                let omega = std::f64::consts::TAU * turns;
                let (s, c) = libm::sincos(omega * t);
                let local = [radius * c, radius * s, pitch * turns * t];
                let tangent = [-radius * omega * s, radius * omega * c, pitch * turns];
                let transform = |v: V| {
                    add(
                        add(mul(axes[0], v[0]), mul(axes[1], v[1])),
                        mul(axes[2], v[2]),
                    )
                };
                Ok((add(*origin, transform(local)), transform(tangent)))
            }
            Self::Spline(curve) => {
                let [lo, hi] = curve.domain();
                let parameter = lo + &((hi - lo) * rational(t)?);
                let jet = curve.jet_box([&parameter, &parameter], 128, dec)?;
                Ok((
                    jet.value.map(|v| v.mid().to_f64_nearest()),
                    jet.first.map(|v| v.mid().to_f64_nearest()),
                ))
            }
        }
    }
    fn minimum_sections(&self) -> usize {
        match self {
            Self::Helix { turns, .. } => (turns.abs() * 16.0).ceil() as usize,
            Self::Spline(curve) => curve.knots().len().max(16),
            Self::Line(..) => 1,
        }
    }
}

#[derive(Clone, Copy)]
struct Frame {
    origin: V,
    x: V,
    y: V,
    scale: f64,
}
impl Frame {
    fn at(
        path: &Path,
        guide: Option<&Path>,
        reference: V,
        initial: Option<[V; 2]>,
        normal_offset: f64,
        twist: f64,
        end_scale: f64,
        t: f64,
        orientation: catallax_recipe::SweepOrientation,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let base = |t: f64, dec: &mut Decider| -> Result<[V; 3], Refusal> {
            let (origin, tangent) = path.at(t, dec)?;
            let tangent = unit(tangent)?;
            if guide.is_none() {
                if let Path::Helix {
                    radius,
                    pitch,
                    turns,
                    axes,
                    ..
                } = path
                {
                    let omega = std::f64::consts::TAU * turns;
                    let (s, c) = libm::sincos(omega * t);
                    let radial = add(mul(axes[0], c), mul(axes[1], s));
                    let binormal = unit(cross(tangent, radial))?;
                    let rise = pitch * turns;
                    let speed = libm::sqrt(radius * radius * omega * omega + rise * rise);
                    let (spin_s, spin_c) = libm::sincos(-omega * rise / speed * t);
                    let x = add(mul(radial, spin_c), mul(binormal, spin_s));
                    return Ok([x, cross(tangent, x), tangent]);
                }
            }
            let direction = if let Some(guide) = guide {
                sub(guide.at(t, dec)?.0, origin)
            } else {
                reference
            };
            let x = unit(sub(direction, mul(tangent, dot(direction, tangent))))?;
            Ok([x, unit(cross(tangent, x))?, tangent])
        };
        let first = base(0.0, dec)?;
        let current = base(t, dec)?;
        let initial = initial.unwrap_or([first[0], first[1]]);
        let carry = |v: V| {
            if orientation == catallax_recipe::SweepOrientation::Fixed && guide.is_none() {
                return v;
            }
            let components = first.map(|axis| dot(v, axis));
            add(
                add(
                    mul(current[0], components[0]),
                    mul(current[1], components[1]),
                ),
                mul(current[2], components[2]),
            )
        };
        let (s, c) = libm::sincos(twist * t);
        let turn = |v: V| {
            add(
                add(mul(v, c), mul(cross(current[2], v), s)),
                mul(current[2], dot(v, current[2]) * (1.0 - c)),
            )
        };
        let x = turn(carry(initial[0]));
        let y = turn(carry(initial[1]));
        let scale = 1.0 + (end_scale - 1.0) * t;
        let origin = add(
            path.at(t, dec)?.0,
            mul(unit(cross(x, y))?, normal_offset * scale),
        );
        Ok(Self {
            origin,
            x,
            y,
            scale,
        })
    }
    fn apply(&self, p: V) -> V {
        add(
            self.origin,
            add(
                mul(self.x, p[0] * self.scale),
                mul(self.y, p[1] * self.scale),
            ),
        )
    }
    fn plane(&self, dec: &mut Decider) -> Result<Plane, Refusal> {
        let x = vector(self.x)?;
        let y = vector(self.y)?;
        Plane::new(point(self.origin)?, x.cross(&y), x, dec)
    }
    fn uv(&self, p: V) -> Result<[Rational; 2], Refusal> {
        let delta = sub(p, self.origin);
        Ok([rational(dot(delta, self.x))?, rational(dot(delta, self.y))?])
    }
}

struct Section {
    t: [Rational; 2],
    frames: [Frame; 4],
}
fn controls(values: [V; 4]) -> [V; 4] {
    let [a, b, c, d] = values;
    [
        a,
        std::array::from_fn(|i| (-5.0 * a[i] + 18.0 * b[i] - 9.0 * c[i] + 2.0 * d[i]) / 6.0),
        std::array::from_fn(|i| (2.0 * a[i] - 9.0 * b[i] + 18.0 * c[i] - 5.0 * d[i]) / 6.0),
        d,
    ]
}
fn cubic(p: &[V; 4], t: f64) -> V {
    let s = 1.0 - t;
    add(
        add(mul(p[0], s * s * s), mul(p[1], 3.0 * s * s * t)),
        add(mul(p[2], 3.0 * s * t * t), mul(p[3], t * t * t)),
    )
}

/// Builds an editable result without promoting fitted geometry into exact answers.
pub fn build(recipe: &SweepRecipe, dec: &mut Decider) -> Result<Swept, Refusal> {
    let (definition, checks, _) = definition(recipe, dec)?;
    let body = numerical::build(definition, options(recipe), dec)?;
    Ok(Swept {
        recipe: recipe.clone(),
        body,
        checks,
    })
}

fn definition(
    recipe: &SweepRecipe,
    dec: &mut Decider,
) -> Result<(NumericalDefinition, SweepChecks, Vec<catallax_recipe::Role>), Refusal> {
    if recipe.name.trim().is_empty()
        || recipe.tolerance_mm.value().sign() != Sign::Positive
        || recipe.end_scale.value().sign() != Sign::Positive
        || recipe.max_sections == 0
    {
        return Err(bad(
            "a name, positive tolerance, positive end scale and nonzero section allowance are required",
        ));
    }
    let profile = crate::Profile::new(
        recipe.plane.as_ref().unwrap_or(&RecipePose::identity()),
        &recipe.loops,
        dec,
    )?;
    if profile.is_empty() {
        return Err(bad("a sweep needs a nonempty profile"));
    }
    let path = Path::new(&recipe.path, dec)?;
    let guide = recipe
        .guide
        .as_ref()
        .map(|p| Path::new(p, dec))
        .transpose()?;
    let twist = number(&recipe.twist_degrees).to_radians();
    let end_scale = number(&recipe.end_scale);
    let mut initial = None;
    let mut normal_offset = 0.0;
    let mut anchor = [0.0; 3];
    if recipe.plane.is_some() {
        let mut coordinate = |x: &Real| -> Result<f64, Refusal> {
            Ok(dec.enclose_budgeted_at(x, 128)?.mid().to_f64_nearest())
        };
        let mut axes = [[0.0; 3]; 2];
        for (i, axis) in axes.iter_mut().enumerate() {
            let mut input = [0; 3];
            input[i] = 1;
            let v = profile.pose.rotation.apply(&Vector3::new(
                input[0].into(),
                input[1].into(),
                input[2].into(),
            ));
            *axis = [coordinate(&v.x)?, coordinate(&v.y)?, coordinate(&v.z)?];
        }
        let origin = [
            coordinate(&profile.pose.translation.x)?,
            coordinate(&profile.pose.translation.y)?,
            coordinate(&profile.pose.translation.z)?,
        ];
        anchor = path.at(0.0, dec)?.0;
        normal_offset = dot(sub(origin, anchor), unit(cross(axes[0], axes[1]))?);
        initial = Some(axes);
    }
    let reference=initial.map(|axes|axes[0])
        .unwrap_or_else(||recipe.reference.each_ref().map(number));
    let frame = |t, dec: &mut Decider| {
        Frame::at(
            &path,
            guide.as_ref(),
            reference,
            initial,
            normal_offset,
            twist,
            end_scale,
            t,
            recipe.orientation,
            dec,
        )
    };
    let mut charts = Vec::new();
    let mut sample_controls = Vec::new();
    let mut profile_error = 0.0_f64;
    let mut profile_samples = 0;
    let to_section = |world: V| {
        if let Some(axes) = initial {
            let delta = sub(world, anchor);
            [dot(delta, axes[0]), dot(delta, axes[1]), 0.0]
        } else {
            world
        }
    };
    for (lp_index, lp) in recipe.loops.iter().enumerate() {
        let start = lp.0.first().ok_or_else(|| bad("empty profile loop"))?.id();
        let (parts, ids, forward) = crate::loft::charts::build(&profile, lp, start, dec)?;
        for part in parts {
            let coordinates = part
                .curve
                .controls()
                .map(|p| {
                    let world = [
                        dec.enclose_budgeted_at(&p.x, 128)?.mid().to_f64_nearest(),
                        dec.enclose_budgeted_at(&p.y, 128)?.mid().to_f64_nearest(),
                        dec.enclose_budgeted_at(&p.z, 128)?.mid().to_f64_nearest(),
                    ];
                    Ok(to_section(world))
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            let fitted_weights = part
                .curve
                .weights()
                .map(|weight| Ok(Real::from(dec.enclose_budgeted_at(weight, 128)?.mid())))
                .collect::<Result<Vec<_>, Refusal>>()?;
            let local = BSplineCurve::new(
                part.curve.degree(),
                part.curve.knots().to_vec(),
                coordinates
                    .iter()
                    .map(|p| point(*p))
                    .collect::<Result<Vec<_>, _>>()?,
                fitted_weights,
                dec,
            )?;
            for i in 0..=8 {
                let t = &part.bounds[0]
                    + &((&part.bounds[1] - &part.bounds[0])
                        * Rational::from_ratio(i, 8).expect("fixed samples"));
                let source = part
                    .curve
                    .evaluate_box([&t, &t], 128, dec)?
                    .map(|b| b.mid().to_f64_nearest());
                let emitted = local
                    .evaluate_box([&t, &t], 128, dec)?
                    .map(|b| b.mid().to_f64_nearest());
                let delta = sub(to_section(source), emitted);
                profile_error = profile_error.max(libm::sqrt(dot(delta, delta)) * end_scale.max(1.0));
                profile_samples += 1;
            }
            sample_controls.extend(coordinates.iter().copied());
            let mid = (&part.bounds[0] + &part.bounds[1]).mul_pow2(-1);
            let segment = ids
                .iter()
                .enumerate()
                .find(|(i, _)| {
                    mid < Rational::from_ratio((i + 1) as i64, ids.len() as i64)
                        .expect("segment count")
                })
                .expect("authored interval")
                .1
                .clone();
            let fitted_chart = crate::loft::charts::Chart {
                bounds: part.bounds,
                curve: local,
            };
            charts.push((lp_index, fitted_chart, coordinates, forward, segment));
        }
    }
    let mut n = path
        .minimum_sections()
        .max((twist.abs() / (std::f64::consts::PI / 8.0)).ceil() as usize)
        .max(1);
    if let Some(g) = &guide {
        n = n.max(g.minimum_sections());
    }
    let tolerance = number(&recipe.tolerance_mm);
    if !profile_error.is_finite() || profile_error > tolerance / 4.0 {
        return Err(bad(
            "the numerical profile proposal does not meet the stated sample tolerance",
        ));
    }
    let (sections, checks) = loop {
        if n > recipe.max_sections {
            return Err(bad(
                "sampled fitting needs more than max_sections; no tolerance was enlarged",
            ));
        }
        dec.charge_work(n as u64)?;
        let mut sections = Vec::with_capacity(n);
        let mut checks = SweepChecks {
            sections: n,
            checked_points: profile_samples,
            max_sample_error_mm: 0.0,
            max_profile_sample_error_mm: profile_error,
        };
        for i in 0..n {
            let t = |j| (i as f64 + j as f64 / 3.0) / n as f64;
            let frames = [
                frame(t(0), dec)?,
                frame(t(1), dec)?,
                frame(t(2), dec)?,
                frame(t(3), dec)?,
            ];
            for p in &sample_controls {
                let fit = controls(frames.map(|f| f.apply(*p)));
                for j in 1..8 {
                    let local = j as f64 / 8.0;
                    let target = frame((i as f64 + local) / n as f64, dec)?.apply(*p);
                    let delta = sub(cubic(&fit, local), target);
                    let error = libm::sqrt(dot(delta, delta));
                    if !error.is_finite() {
                        return Err(bad("nonfinite sampled fitting error"));
                    }
                    checks.checked_points += 1;
                    checks.max_sample_error_mm = checks.max_sample_error_mm.max(error);
                }
            }
            sections.push(Section {
                t: [
                    Rational::from_ratio(i as i64, n as i64).expect("section count"),
                    Rational::from_ratio((i + 1) as i64, n as i64).expect("section count"),
                ],
                frames,
            });
        }
        if checks.max_sample_error_mm <= tolerance / 4.0 {
            break (sections, checks);
        }
        n = n
            .checked_mul(2)
            .ok_or_else(|| bad("section count overflow"))?;
    };
    let first_frame = frame(0.0, dec)?;
    let orientation = dot(cross(first_frame.x, first_frame.y), path.at(0.0, dec)?.1);
    if orientation == 0.0 {
        return Err(bad("the initial sweep direction lies in the profile plane"));
    }
    let positive = orientation > 0.0;
    let mut assembly = Assembly::new(recipe.name.clone(), recipe.loops.len(), positive);
    let linear = matches!(path, Path::Line(..)) && guide.is_none() && twist == 0.0;
    let degree = if linear { 1 } else { 3 };
    for (lp, chart, coordinates, forward, segment) in charts {
        for (station, section) in sections.iter().enumerate() {
            dec.charge_work(1)?;
            let mut points = Vec::with_capacity(coordinates.len() * (degree + 1));
            let mut weights = Vec::with_capacity(coordinates.len() * (degree + 1));
            for (i, p) in coordinates.iter().enumerate() {
                if linear {
                    for frame in [section.frames[0], section.frames[3]] {
                        let row=affine_row(frame,&[rational(p[0])?,rational(p[1])?])?;
                        let [x,y,z]=row;
                        points.push(Point3::new(x.into(),y.into(),z.into()));
                        weights.push(chart.curve.weight(i).clone());
                    }
                } else {
                    for p in exact_controls(section.frames, *p)? {
                        points.push(p);
                        weights.push(chart.curve.weight(i).clone());
                    }
                }
            }
            let mut vknots = vec![section.t[0].clone(); degree + 1];
            vknots.extend(vec![section.t[1].clone(); degree + 1]);
            let carrier = BSplineSurface::new(
                [chart.curve.degree(), degree],
                [chart.curve.knots().to_vec(), vknots],
                [coordinates.len(), degree + 1],
                points,
                weights,
                dec,
            )?;
            let span =
                catallax_numbers::OperationSpan::begin("brep", &recipe.name, "sweep-patch", dec);
            let patch = RationalPatch::new(carrier, dec);
            span.finish(dec, patch.as_ref().err(), None);
            let patch = patch?;
            assembly.side(
                lp,
                &segment,
                station,
                &patch,
                forward,
                station == 0,
                station + 1 == sections.len(),
                &section.frames[0],
                &section.frames[3],
            )?;
        }
    }
    assembly.caps(frame(0.0, dec)?, frame(1.0, dec)?, dec)?;
    Ok((assembly.definition, checks, assembly.roles))
}

fn affine_row(frame:Frame,p:&[Rational;2])->Result<[Rational;3],Refusal> {
    let scale=rational(frame.scale)?;
    let mut row=[q(0),q(0),q(0)];
    for (axis,value) in row.iter_mut().enumerate() {
        *value=rational(frame.origin[axis])? +
            &scale * &(&rational(frame.x[axis])? * &p[0] + &rational(frame.y[axis])? * &p[1]);
    }
    Ok(row)
}

fn exact_controls(frames: [Frame; 4], p: V) -> Result<[Point3; 4], Refusal> {
    let p=[rational(p[0])?,rational(p[1])?];
    let mut rows: [[Rational; 3]; 4] = std::array::from_fn(|_| std::array::from_fn(|_| q(0)));
    for (row,frame) in rows.iter_mut().zip(frames) {
        *row=affine_row(frame,&p)?;
    }
    let middle = |coefficients: [i64; 4]| -> [Rational; 3] {
        std::array::from_fn(|axis| {
            (0..4).fold(q(0), |sum, i| sum + &q(coefficients[i]) * &rows[i][axis]) / q(6)
        })
    };
    let values = [
        rows[0].clone(),
        middle([-5, 18, -9, 2]),
        middle([2, -9, 18, -5]),
        rows[3].clone(),
    ];
    Ok(values.map(|[x, y, z]| Point3::new(x.into(), y.into(), z.into())))
}

type VertexKey = (usize, Rational, Rational);
type EdgeKey = (usize, usize, Rational, Rational, Rational);
struct Assembly {
    definition: NumericalDefinition,
    vertices: BTreeMap<VertexKey, usize>,
    edges: BTreeMap<EdgeKey, usize>,
    caps: Vec<[Vec<(Rational, NumericalTrim)>; 2]>,
    positive: bool,
    roles: Vec<catallax_recipe::Role>,
}
impl Assembly {
    fn new(name: String, loops: usize, positive: bool) -> Self {
        Self {
            definition: NumericalDefinition {
                name,
                vertices: Vec::new(),
                edges: Vec::new(),
                faces: Vec::new(),
            },
            vertices: BTreeMap::new(),
            edges: BTreeMap::new(),
            caps: (0..loops).map(|_| [Vec::new(), Vec::new()]).collect(),
            positive,
            roles: Vec::new(),
        }
    }
    fn vertex(
        &mut self,
        lp: usize,
        u: &Rational,
        v: &Rational,
        p: &Point3,
    ) -> Result<usize, Refusal> {
        let u = if *u == q(1) { q(0) } else { u.clone() };
        let key = (lp, u, v.clone());
        if let Some(id) = self.vertices.get(&key) {
            return Ok(*id);
        }
        let coordinate = |x: &Real| {
            x.as_rational()
                .ok_or_else(|| bad("fitted control is not rational"))
        };
        let position = [coordinate(&p.x)?, coordinate(&p.y)?, coordinate(&p.z)?];
        let id = self.definition.vertices.len();
        self.definition.vertices.push(NumericalVertex {
            key: format!("vertex/{lp}/{}/{v}", key.1),
            position,
        });
        self.vertices.insert(key, id);
        Ok(id)
    }
    fn side(
        &mut self,
        lp: usize,
        segment: &catallax_recipe::ElementId,
        station: usize,
        patch: &RationalPatch,
        forward: bool,
        first: bool,
        last: bool,
        start: &Frame,
        end: &Frame,
    ) -> Result<(), Refusal> {
        let forward = forward == self.positive;
        let b = patch.domain();
        let carrier = patch.carrier();
        let [nu, nv] = carrier.control_count();
        let p = carrier.controls();
        let uv = [
            [b[0][0].clone(), b[1][0].clone()],
            [b[0][1].clone(), b[1][0].clone()],
            [b[0][1].clone(), b[1][1].clone()],
            [b[0][0].clone(), b[1][1].clone()],
        ];
        let indices = [0, (nu - 1) * nv, nu * nv - 1, nv - 1];
        let mut vertices = [0; 4];
        for i in 0..4 {
            vertices[i] = self.vertex(lp, &uv[i][0], &uv[i][1], &p[indices[i]])?;
        }
        let norm = |u: &Rational| if *u == q(1) { q(0) } else { u.clone() };
        let descriptions = [
            (
                SplineBoundary::VMin,
                (lp, 0, b[1][0].clone(), b[0][0].clone(), b[0][1].clone()),
                [0, 1],
                Sense::Same,
            ),
            (
                SplineBoundary::UMax,
                (lp, 1, norm(&b[0][1]), b[1][0].clone(), b[1][1].clone()),
                [1, 2],
                Sense::Same,
            ),
            (
                SplineBoundary::VMax,
                (lp, 0, b[1][1].clone(), b[0][0].clone(), b[0][1].clone()),
                [3, 2],
                Sense::Reversed,
            ),
            (
                SplineBoundary::UMin,
                (lp, 1, norm(&b[0][0]), b[1][0].clone(), b[1][1].clone()),
                [0, 3],
                Sense::Reversed,
            ),
        ];
        let mut trims = Vec::new();
        for (side, key, ends, mut sense) in descriptions {
            let edge = if let Some(id) = self.edges.get(&key) {
                *id
            } else {
                let id = self.definition.edges.len();
                self.definition.edges.push(NumericalEdge {
                    key: format!("edge/{lp}/{}/{}/{}/{}", key.1, key.2, key.3, key.4),
                    curve: Curve::BSpline(patch.boundary(side)),
                    vertices: ends.map(|i| vertices[i]),
                    parameters: [key.3.clone(), key.4.clone()],
                });
                self.edges.insert(key, id);
                id
            };
            if !forward {
                sense = flip(sense);
            }
            let mut ends = ends;
            if sense == Sense::Reversed {
                ends.reverse();
            }
            let trim = NumericalTrim {
                edge,
                sense,
                uv: ends.map(|i| uv[i].clone()),
            };
            if (first && side == SplineBoundary::VMin) || (last && side == SplineBoundary::VMax) {
                let cap = usize::from(side == SplineBoundary::VMax);
                let frame = if cap == 0 { start } else { end };
                let mut ends = ends;
                ends.reverse();
                let positions = ends.map(|i| {
                    let at = &p[indices[i]];
                    [
                        at.x.as_rational().expect("fitted control").to_f64_nearest(),
                        at.y.as_rational().expect("fitted control").to_f64_nearest(),
                        at.z.as_rational().expect("fitted control").to_f64_nearest(),
                    ]
                });
                self.caps[lp][cap].push((
                    b[0][0].clone(),
                    NumericalTrim {
                        edge,
                        sense: flip(sense),
                        uv: [frame.uv(positions[0])?, frame.uv(positions[1])?],
                    },
                ));
            }
            trims.push(trim);
        }
        if !forward {
            trims.reverse();
        }
        self.definition.faces.push(NumericalFace {
            key: format!("side/{lp}/{segment}/{station}/{}", b[0][0]),
            surface: Surface::BSpline(patch.clone().into()),
            sense: if forward {
                Sense::Same
            } else {
                Sense::Reversed
            },
            loops: vec![trims],
        });
        self.roles.push(catallax_recipe::Role::Side {
            segment: segment.clone(),
            path: None,
        });
        Ok(())
    }
    fn caps(&mut self, start: Frame, end: Frame, dec: &mut Decider) -> Result<(), Refusal> {
        for (cap, frame) in [start, end].iter().enumerate() {
            let mut loops = Vec::new();
            for pair in &mut self.caps {
                let mut rows = std::mem::take(&mut pair[cap]);
                rows.sort_by(|a, b| a.0.cmp(&b.0));
                let mut remaining: Vec<_> = rows.into_iter().map(|(_, t)| t).collect();
                let mut walk = Vec::new();
                let first = remaining.pop().ok_or_else(|| bad("missing cap boundary"))?;
                let finish = |t: &NumericalTrim| {
                    self.definition.edges[t.edge].vertices[usize::from(t.sense == Sense::Same)]
                };
                let mut at = finish(&first);
                walk.push(first);
                while !remaining.is_empty() {
                    let next = remaining
                        .iter()
                        .position(|t| {
                            self.definition.edges[t.edge].vertices
                                [usize::from(t.sense == Sense::Reversed)]
                                == at
                        })
                        .ok_or_else(|| bad("cap boundary does not form a connected walk"))?;
                    let t = remaining.remove(next);
                    at = finish(&t);
                    walk.push(t);
                }
                loops.push(walk);
            }
            self.definition.faces.push(NumericalFace {
                key: if cap == 0 { "start" } else { "end" }.into(),
                surface: Surface::Plane(frame.plane(dec)?),
                sense: if (cap == 0) == self.positive {
                    Sense::Reversed
                } else {
                    Sense::Same
                },
                loops,
            });
            self.roles.push(if cap == 0 {
                catallax_recipe::Role::Start
            } else {
                catallax_recipe::Role::End
            });
        }
        Ok(())
    }
}
