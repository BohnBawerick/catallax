//! Certified source-incidence measurements. Samples provide lower witnesses only.

use catallax_numbers::{Ball, Decider, Rational, Real, RealPoly, Refusal};
use catallax_shapes::{BSplineCurve, Point3};

use crate::exact::{dot, sub};
use crate::{Curve, Face, ImportRefusal, IncidenceGap, OrientedEdge, Placement, Surface, Vec3};

const BITS: u32 = 128;
const CELLS: usize = 4096;

fn ball(q: &Rational) -> Ball {
    Ball::from_rational(q, BITS)
}
fn range(a: &Rational, b: &Rational) -> Ball {
    ball(a).hull(&ball(b))
}
fn width() -> Rational {
    Rational::from_ratio(1, 1_000_000_000).unwrap()
}
fn error(detail: &str) -> Refusal {
    Refusal::new("import.incidence-measurement-unproven", detail)
}
fn constant(q: &Rational) -> RealPoly {
    RealPoly::constant(1, Real::from_rational(q.clone()))
}
fn choose(n: usize, k: usize) -> Rational {
    (0..k.min(n - k)).fold(Rational::one(), |v, i| {
        v * Rational::from_i64((n - i) as i64) / Rational::from_i64((i + 1) as i64)
    })
}

// Convert the univariate power polynomial to Bernstein form on this cell.
// Its control hull encloses every value, including stationary interior extrema.
fn polynomial_box(p: &RealPoly, t: &Ball, dec: &mut Decider) -> Result<Ball, Refusal> {
    let terms = p.terms();
    let degree = terms
        .iter()
        .map(|(m, _)| m.exponents()[0] as usize)
        .max()
        .unwrap_or(0);
    let mut coefficients = vec![Real::zero(); degree + 1];
    for (m, c) in terms {
        coefficients[m.exponents()[0] as usize] = c;
    }
    let a = t.lo();
    let h = &t.hi() - &a;
    let mut powers = vec![Real::zero(); degree + 1];
    for (j, value) in powers.iter_mut().enumerate() {
        dec.check_deadline()?;
        let hp = h
            .pow(i32::try_from(j).map_err(|_| error("polynomial degree cannot be represented"))?)
            .ok_or_else(|| error("polynomial power unproved"))?;
        for (k, c) in coefficients.iter().enumerate().skip(j) {
            let ap = a
                .pow(
                    i32::try_from(k - j)
                        .map_err(|_| error("polynomial degree cannot be represented"))?,
                )
                .ok_or_else(|| error("polynomial power unproved"))?;
            *value = &*value + &(c * &Real::from_rational(choose(k, j) * ap * hp.clone()));
        }
    }
    let mut result: Option<Ball> = None;
    for i in 0..=degree {
        let value = (0..=i).fold(Real::zero(), |v, j| {
            v + &powers[j] * &Real::from_rational(choose(i, j) / choose(degree, j))
        });
        let bound = dec.enclose_budgeted_at(&value, BITS)?;
        result = Some(result.map_or_else(|| bound.clone(), |old| old.hull(&bound)));
    }
    Ok(result.expect("nonempty polynomial basis"))
}
fn divide(a: &Ball, b: &Ball) -> Result<Ball, Refusal> {
    a.checked_div(b)
        .ok_or_else(|| error("distance denominator includes zero"))
}
fn root(a: &Ball) -> Result<Ball, Refusal> {
    if a.hi() < Rational::zero() {
        return Err(error("negative squared distance enclosure"));
    }
    range(&a.lo().max(Rational::zero()), &a.hi())
        .sqrt()
        .ok_or_else(|| error("distance square root unproved"))
}
fn absolute(a: &Ball) -> (Rational, Rational) {
    let lo = a.lo();
    let hi = a.hi();
    let lower = if lo <= Rational::zero() && hi >= Rational::zero() {
        Rational::zero()
    } else {
        lo.abs().min(hi.abs())
    };
    (lower, lo.abs().max(hi.abs()))
}

#[derive(Clone, Debug)]
struct Bounds {
    lower: Rational,
    upper: Rational,
    converged: bool,
}

pub(crate) fn target(refusal: &ImportRefusal) -> Option<(String, Option<String>, &'static str)> {
    let ImportRefusal::NotSolid { detail } = refusal else {
        return None;
    };
    let step = if detail.starts_with("spline endpoint incidence:") {
        "spline_endpoint_incidence"
    } else if detail.starts_with("spline curve-surface incidence:") {
        "spline_curve_surface_incidence"
    } else if detail.starts_with("a vertex of edge ") {
        "vertex_curve_incidence"
    } else if detail.starts_with("edge ") && detail.contains("does not lie exactly on face") {
        "edge_surface_incidence"
    } else {
        return None;
    };
    let names = detail
        .split_whitespace()
        .filter(|s| s.starts_with('#'))
        .collect::<Vec<_>>();
    Some((
        names.first()?.to_string(),
        names.get(1).map(|s| s.to_string()),
        step,
    ))
}

fn point(p: &Vec3) -> Point3 {
    Point3::new(
        Real::from_rational(p[0].clone()),
        Real::from_rational(p[1].clone()),
        Real::from_rational(p[2].clone()),
    )
}

fn point_curve(p: &Vec3, curve: &Curve, dec: &mut Decider) -> Result<Bounds, Refusal> {
    let value = match curve {
        Curve::Line(o, d) => {
            let v = sub(p, o);
            let d2 = dot(d, d);
            let t = dot(&v, d) / d2;
            let residual = std::array::from_fn(|i| &v[i] - &(&t * &d[i]));
            Real::from_rational(dot(&residual, &residual)).sqrt()
        }
        Curve::Circle(c, r) => {
            let v = sub(p, &c.origin);
            let h = dot(&v, &c.axis);
            let n2 = dot(&c.axis, &c.axis);
            let h2 = &h * &h / n2;
            let radial2 = &dot(&v, &v) - &h2;
            let radial = Real::from_rational(radial2).sqrt() - Real::from_rational(r.clone());
            (Real::from_rational(h2) + radial.square()).sqrt()
        }
        Curve::BSpline(s) => return point_spline(p, &s.geometry.to_curve(dec)?, dec),
        Curve::Ellipse(..) => return Err(error("ellipse endpoint distance is not covered")),
    };
    let enclosure = dec.enclose_budgeted_at(&value, BITS)?;
    Ok(Bounds {
        lower: enclosure.lo().max(Rational::zero()),
        upper: enclosure.hi(),
        converged: enclosure.hi() - enclosure.lo() <= width(),
    })
}

struct PointCell {
    span: usize,
    lo: Rational,
    hi: Rational,
    lower: Rational,
    upper: Rational,
}

fn point_cell(
    p: &Vec3,
    curve: &BSplineCurve,
    span: usize,
    lo: Rational,
    hi: Rational,
    dec: &mut Decider,
) -> Result<PointCell, Refusal> {
    let controls = curve.span_bezier_controls_on(span, [&lo, &hi], dec)?;
    let mut axes: [Option<Ball>; 3] = std::array::from_fn(|_| None);
    for control in controls {
        for i in 0..3 {
            let coordinate = dec.enclose_budgeted_at(&(&control[i] / &control[3]), BITS)?;
            axes[i] = Some(
                axes[i]
                    .take()
                    .map_or_else(|| coordinate.clone(), |old| old.hull(&coordinate)),
            );
        }
    }
    let mut d2 = ball(&Rational::zero());
    for i in 0..3 {
        d2 = d2 + (axes[i].as_ref().unwrap() - &ball(&p[i])).square();
    }
    let lower = root(&d2)?.lo().max(Rational::zero());
    let mid = (&lo + &hi).mul_pow2(-1);
    let witness = curve.point_at(&Real::from_rational(mid), dec)?;
    let distance = witness.minus(&point(p)).norm_squared().sqrt();
    let upper = dec.enclose_budgeted_at(&distance, BITS)?.hi();
    Ok(PointCell {
        span,
        lo,
        hi,
        lower,
        upper,
    })
}

fn point_spline(p: &Vec3, curve: &BSplineCurve, dec: &mut Decider) -> Result<Bounds, Refusal> {
    let mut cells = Vec::new();
    let mut upper: Option<Rational> = None;
    for (span, [a, b]) in curve.spans() {
        let cell = point_cell(p, curve, span, a.clone(), b.clone(), dec)?;
        upper = Some(upper.map_or_else(|| cell.upper.clone(), |old| old.min(cell.upper.clone())));
        cells.push(cell);
        for t in [a, b] {
            let value = curve
                .point_at(&Real::from_rational(t.clone()), dec)?
                .minus(&point(p))
                .norm_squared()
                .sqrt();
            upper = Some(
                upper
                    .unwrap()
                    .min(dec.enclose_budgeted_at(&value, BITS)?.hi()),
            );
        }
    }
    let mut upper = upper.ok_or_else(|| error("empty spline domain"))?;
    for _ in 0..CELLS {
        dec.check_deadline()?;
        cells.retain(|c| c.lower <= upper);
        let Some((index, cell)) = cells
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.lower.cmp(&b.lower))
        else {
            return Ok(Bounds {
                lower: upper.clone(),
                upper,
                converged: true,
            });
        };
        let lower = cell.lower.clone();
        if &upper - &lower <= width() {
            return Ok(Bounds {
                lower,
                upper,
                converged: true,
            });
        }
        let cell = cells.swap_remove(index);
        let mid = (&cell.lo + &cell.hi).mul_pow2(-1);
        for [a, b] in [[cell.lo, mid.clone()], [mid, cell.hi]] {
            let next = point_cell(p, curve, cell.span, a, b, dec)?;
            upper = upper.min(next.upper.clone());
            cells.push(next);
        }
    }
    let lower = cells
        .iter()
        .map(|c| c.lower.clone())
        .min()
        .unwrap_or_else(|| upper.clone())
        .min(upper.clone());
    Ok(Bounds {
        lower,
        upper,
        converged: false,
    })
}

fn frame(p: &Placement) -> [[Real; 3]; 2] {
    let axis = p.axis.clone().map(Real::from_rational);
    let reference = p.ref_dir.clone().map(Real::from_rational);
    let norm = |v: &[Real; 3]| {
        v.iter()
            .fold(Real::zero(), |sum, x| sum + x.square())
            .sqrt()
    };
    let a = norm(&axis);
    let r = norm(&reference);
    let n = axis.map(|v| &v / &a);
    let u = reference.map(|v| &v / &r);
    let v = [
        &n[1] * &u[2] - &n[2] * &u[1],
        &n[2] * &u[0] - &n[0] * &u[2],
        &n[0] * &u[1] - &n[1] * &u[0],
    ];
    [u, v]
}

struct Span {
    p: [RealPoly; 4],
    lo: Rational,
    hi: Rational,
}

fn spans(edge: &OrientedEdge, dec: &mut Decider) -> Result<(Vec<Span>, bool), Refusal> {
    let one = constant(&Rational::one());
    let t = RealPoly::var(1, 0);
    match &edge.curve {
        Curve::Line(o, d) => {
            let d2 = dot(d, d);
            let parameters =
                [&edge.start.point, &edge.end.point].map(|v| &dot(&sub(v, o), d) / &d2);
            let ends =
                parameters.map(|t| std::array::from_fn::<_, 3, _>(|i| &o[i] + &(&t * &d[i])));
            let p = std::array::from_fn(|i| {
                if i == 3 {
                    one.clone()
                } else {
                    constant(&ends[0][i]) + t.scale(&Real::from_rational(&ends[1][i] - &ends[0][i]))
                }
            });
            Ok((
                vec![Span {
                    p,
                    lo: Rational::zero(),
                    hi: Rational::one(),
                }],
                true,
            ))
        }
        Curve::BSpline(s) => {
            let curve = s.geometry.to_curve(dec)?;
            let clamped = s.geometry.clone().into_clamped(dec)?;
            let controls = clamped.controls();
            let ends = if edge.curve_same_sense {
                [&edge.start.point, &edge.end.point]
            } else {
                [&edge.end.point, &edge.start.point]
            };
            let complete = ends[0] == &controls[0] && ends[1] == &controls[controls.len() - 1];
            let spans = curve
                .spans()
                .map(|(span, [lo, hi])| {
                    Ok(Span {
                        p: curve.span_polynomials(span, dec)?,
                        lo: lo.clone(),
                        hi: hi.clone(),
                    })
                })
                .collect::<Result<_, Refusal>>()?;
            Ok((spans, complete))
        }
        Curve::Circle(c, r) => {
            let [u, v] = frame(c);
            let w = &one + &t.pow_u32(2);
            let x = &one - &t.pow_u32(2);
            let y = t.scale(&Real::from_i64(2));
            let mut spans = Vec::new();
            for quadrant in 0..4 {
                let (x, y) = match quadrant {
                    0 => (x.clone(), y.clone()),
                    1 => (-y.clone(), x.clone()),
                    2 => (-x.clone(), -y.clone()),
                    _ => (y.clone(), -x.clone()),
                };
                let p = std::array::from_fn(|i| {
                    if i == 3 {
                        w.clone()
                    } else {
                        w.scale(&Real::from_rational(c.origin[i].clone()))
                            + x.scale(&(&u[i] * &Real::from_rational(r.clone())))
                            + y.scale(&(&v[i] * &Real::from_rational(r.clone())))
                    }
                });
                spans.push(Span {
                    p,
                    lo: Rational::zero(),
                    hi: Rational::one(),
                });
            }
            Ok((spans, false))
        }
        Curve::Ellipse(..) => Err(error("ellipse/surface distance is not covered")),
    }
}

struct Metric {
    weight: RealPoly,
    height: RealPoly,
    squared: RealPoly,
    residual: RealPoly,
    axis_norm: Ball,
    radius: Ball,
    minor: Option<Ball>,
    slope: Option<(Ball, Ball)>,
    plane: bool,
}

impl Metric {
    fn new(p: &[RealPoly; 4], surface: &Surface) -> Result<Self, Refusal> {
        let (placement, radius) = match surface {
            Surface::Plane(p) => (p, Rational::zero()),
            Surface::Cylinder(p, r)
            | Surface::Sphere(p, r)
            | Surface::Cone45(p, r)
            | Surface::Cone(p, r, _) => (p, r.clone()),
            Surface::Torus(p, r, _) => (p, r.clone()),
            Surface::BSpline(_) => {
                return Err(error(
                    "spline surface distance needs a certified closest-point chart",
                ));
            }
        };
        let d: [RealPoly; 3] = std::array::from_fn(|i| {
            &p[i] - &p[3].scale(&Real::from_rational(placement.origin[i].clone()))
        });
        let height = (0..3).fold(RealPoly::zero(1), |sum, i| {
            sum + d[i].scale(&Real::from_rational(placement.axis[i].clone()))
        });
        let n2 = dot(&placement.axis, &placement.axis);
        let norm = ball(&n2).sqrt().ok_or_else(|| error("zero surface axis"))?;
        let full_squared = d
            .iter()
            .fold(RealPoly::zero(1), |sum, x| sum + x.pow_u32(2));
        let mut squared = full_squared.clone();
        if !matches!(surface, Surface::Sphere(..) | Surface::Plane(..)) {
            squared = squared
                - height
                    .pow_u32(2)
                    .scale(&Real::from_rational(&Rational::one() / &n2));
        }
        let residual = if let Surface::Torus(_, _, minor) = surface {
            let a = full_squared
                + p[3].pow_u32(2).scale(&Real::from_rational(
                    &(&radius * &radius) - &(minor * minor),
                ));
            a.pow_u32(2)
                - (squared.clone() * p[3].pow_u32(2)).scale(&Real::from_rational(
                    &(&radius * &radius) * &Rational::from_i64(4),
                ))
        } else {
            &squared
                - &p[3]
                    .pow_u32(2)
                    .scale(&Real::from_rational(&radius * &radius))
        };
        let slope = match surface {
            Surface::Cone45(..) => Some((
                ball(&Rational::one()),
                ball(&Rational::from_ratio(1, 2).unwrap()).sqrt().unwrap(),
            )),
            Surface::Cone(_, _, angle) => {
                let mut a = ball(&angle.value);
                if angle.degrees {
                    a = a * Ball::pi(BITS) * ball(&Rational::from_ratio(1, 180).unwrap());
                }
                let cosine = a.cos();
                Some((divide(&a.sin(), &cosine)?, cosine))
            }
            _ => None,
        };
        Ok(Self {
            weight: p[3].clone(),
            height,
            squared,
            residual,
            axis_norm: norm,
            radius: ball(&radius),
            minor: if let Surface::Torus(_, _, r) = surface {
                Some(ball(r))
            } else {
                None
            },
            slope,
            plane: matches!(surface, Surface::Plane(_)),
        })
    }

    fn evaluate(&self, t: &Ball, dec: &mut Decider) -> Result<(Rational, Rational), Refusal> {
        let w = polynomial_box(&self.weight, t, dec)?;
        let h = polynomial_box(&self.height, t, dec)?;
        if self.plane {
            return Ok(absolute(&divide(&h, &(&w * &self.axis_norm))?));
        }
        let squared = polynomial_box(&self.squared, t, dec)?;
        let radial = root(&squared)?;
        if let Some(minor) = &self.minor {
            let rho = divide(&radial, &w)?;
            let z = divide(&h, &(&w * &self.axis_norm))?;
            let d = root(&((&rho - &self.radius).square() + z.square()))?;
            let factor = (&rho + &self.radius).square() + z.square() - minor.square();
            let floor = self.radius.square() - minor.square();
            if floor.lo() <= Rational::zero() {
                return Err(error(
                    "torus distance factor needs a ring-torus certificate",
                ));
            }
            let factor = range(&factor.lo().max(floor.lo()), &factor.hi());
            // F = (D-r)(D+r)((rho+R)^2+z^2-r^2). The last factor is
            // at least R^2-r^2 on a ring torus, so division is certified.
            let residual = polynomial_box(&self.residual, t, dec)?;
            let denominator = w.square().square() * (&d + minor) * factor;
            return Ok(absolute(&divide(&residual, &denominator)?));
        }
        if let Some((slope, cosine)) = &self.slope {
            let rho = divide(&radial, &w)?;
            let z = divide(&h, &(&w * &self.axis_norm))?;
            let radius = &self.radius + &(&z * slope);
            // The double-cone distance is a lower bound for either nappe.
            // Equal-height matching on the declared radius gives an upper bound
            // only where that radius is nonnegative.
            if radius.lo() < Rational::zero() {
                return Err(error("cone distance needs a retained-nappe certificate"));
            }
            let residual = rho - radius;
            let (lo, hi) = absolute(&residual);
            return Ok(((&ball(&lo) * cosine).lo().max(Rational::zero()), hi));
        }
        let residual = polynomial_box(&self.residual, t, dec)?;
        Ok(absolute(&divide(
            &residual,
            &(&w * &(&radial + &(&self.radius * &w))),
        )?))
    }
}

fn edge_surface(
    edge: &OrientedEdge,
    face: &Face,
    dec: &mut Decider,
) -> Result<(Bounds, bool), Refusal> {
    let (spans, complete) = spans(edge, dec)?;
    let metrics = spans
        .iter()
        .map(|s| Metric::new(&s.p, &face.surface))
        .collect::<Result<Vec<_>, _>>()?;
    let mut cells = Vec::new();
    let mut lower = Rational::zero();
    for (index, span) in spans.iter().enumerate() {
        let (cell_lower, upper) = metrics[index].evaluate(&range(&span.lo, &span.hi), dec)?;
        if complete {
            for t in [&span.lo, &span.hi, &(&span.lo + &span.hi).mul_pow2(-1)] {
                lower = lower.max(metrics[index].evaluate(&ball(t), dec)?.0);
            }
        }
        cells.push((index, span.lo.clone(), span.hi.clone(), cell_lower, upper));
    }
    for _ in 0..CELLS {
        dec.check_deadline()?;
        let (index, cell) = cells
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.4.cmp(&b.4))
            .ok_or_else(|| error("empty edge cover"))?;
        let upper = cell.4.clone();
        let actual_lower = if complete {
            lower.clone()
        } else {
            cells.iter().map(|c| c.3.clone()).min().unwrap()
        };
        if &upper - &actual_lower <= width() {
            return Ok((
                Bounds {
                    lower: actual_lower,
                    upper,
                    converged: true,
                },
                complete,
            ));
        }
        // An enclosing carrier still proves an admission upper bound. Refining
        // its maximum cannot establish the maximum of an unknown subarc.
        if !complete && upper <= crate::gap::auto_cap() {
            return Ok((
                Bounds {
                    lower: actual_lower,
                    upper,
                    converged: false,
                },
                complete,
            ));
        }
        let (span, a, b, _, _) = cells.swap_remove(index);
        let mid = (&a + &b).mul_pow2(-1);
        if complete {
            lower = lower.max(metrics[span].evaluate(&ball(&mid), dec)?.0)
        }
        for [lo, hi] in [[a, mid.clone()], [mid, b]] {
            let (minimum, upper) = metrics[span].evaluate(&range(&lo, &hi), dec)?;
            cells.push((span, lo, hi, minimum, upper));
        }
    }
    let upper = cells.iter().map(|c| c.4.clone()).max().unwrap();
    if !complete {
        lower = cells.iter().map(|c| c.3.clone()).min().unwrap()
    }
    Ok((
        Bounds {
            lower,
            upper,
            converged: false,
        },
        complete,
    ))
}

/// Measures the edge named by a saved first-incidence refusal against its
/// immutable source geometry. Enclosing carriers are labelled explicitly.
pub fn measure_incidence(
    bytes: &[u8],
    refusal: &ImportRefusal,
) -> Result<IncidenceGap, ImportRefusal> {
    let (edge_name, face_name, step) = target(refusal).ok_or_else(|| ImportRefusal::Malformed {
        detail: "expected a saved incidence refusal".into(),
    })?;
    let file = crate::read_step(bytes)?;
    measure_faces(
        file.solids.iter().flat_map(|s| &s.faces),
        edge_name,
        face_name,
        step,
    )
}

pub(crate) fn measure_solid(
    solid: &crate::Solid,
    refusal: &ImportRefusal,
) -> Result<IncidenceGap, ImportRefusal> {
    let (edge_name, face_name, step) = target(refusal).ok_or_else(|| ImportRefusal::Malformed {
        detail: "expected a saved incidence refusal".into(),
    })?;
    measure_faces(solid.faces.iter(), edge_name, face_name, step)
}

fn measure_faces<'a>(
    faces: impl Iterator<Item = &'a Face>,
    edge_name: String,
    face_name: Option<String>,
    step: &str,
) -> Result<IncidenceGap, ImportRefusal> {
    let pair = faces
        .filter(|f| face_name.as_ref().is_none_or(|name| name == &f.name))
        .find_map(|face| {
            face.loops
                .iter()
                .flat_map(|b| &b.edges)
                .find(|e| e.edge == edge_name)
                .map(|e| (face, e))
        })
        .ok_or_else(|| ImportRefusal::Malformed {
            detail: "saved incidence entities are absent from source".into(),
        })?;
    let mut dec = Decider::default();
    let measured = (|| {
        if face_name.is_some() {
            return edge_surface(pair.1, pair.0, &mut dec);
        }
        let a = point_curve(&pair.1.start.point, &pair.1.curve, &mut dec)?;
        let b = point_curve(&pair.1.end.point, &pair.1.curve, &mut dec)?;
        Ok((
            Bounds {
                lower: a.lower.max(b.lower),
                upper: a.upper.max(b.upper),
                converged: a.converged && b.converged,
            },
            true,
        ))
    })()
    .map_err(|e: Refusal| {
        if catallax_numbers::current_operation_work().exhausted {
            ImportRefusal::ResourceLimit
        } else {
            ImportRefusal::Unsupported {
                entity: edge_name.clone(),
                type_name: format!(
                    "EDGE_CURVE (incidence measurement: {}: {})",
                    e.code, e.detail
                ),
            }
        }
    })?;
    let metric = if face_name.is_some() {
        "maximum curve-to-analytic-support distance"
    } else {
        "maximum endpoint distance to the retained source curve"
    };
    Ok(IncidenceGap {
        edge: edge_name,
        face: face_name,
        step: step.into(),
        metric: metric.into(),
        distance_lower_mm: crate::gap::decimal(&measured.0.lower),
        distance_upper_mm: crate::gap::decimal(&measured.0.upper),
        coverage: if measured.1 {
            "finite_source_feature"
        } else {
            "enclosing_retained_carrier"
        }
        .into(),
        converged: measured.0.converged,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_shapes::spline_tools::RationalSplineCurve;

    fn q(s: &str) -> Rational {
        Rational::parse(s).unwrap()
    }
    fn plane() -> Face {
        Face {
            name: "#1".into(),
            surface: Surface::Plane(Placement {
                origin: [q("0"), q("0"), q("0")],
                axis: [q("0"), q("0"), q("1")],
                ref_dir: [q("1"), q("0"), q("0")],
            }),
            same_sense: true,
            loops: Vec::new(),
        }
    }
    fn spline() -> OrientedEdge {
        let controls = [
            [q("0"), q("0"), q("0")],
            [q("1"), q("0"), q("0.002")],
            [q("2"), q("0"), q("0")],
        ];
        OrientedEdge {
            edge: "#2".into(),
            start: crate::Vertex {
                name: "#3".into(),
                point: controls[0].clone(),
            },
            end: crate::Vertex {
                name: "#4".into(),
                point: controls[2].clone(),
            },
            curve_same_sense: true,
            forward: true,
            curve: Curve::BSpline(crate::StepSplineCurve {
                geometry: RationalSplineCurve::new(
                    2,
                    vec![q("0"), q("0"), q("0"), q("1"), q("1"), q("1")],
                    controls.to_vec(),
                    vec![q("1"); 3],
                    &Decider::default(),
                )
                .unwrap(),
                closed: Some(false),
                self_intersect: Some(false),
            }),
        }
    }
    #[test]
    fn quadratic_midspan_miss_has_a_whole_curve_bound() {
        let (b, complete) = edge_surface(&spline(), &plane(), &mut Decider::default()).unwrap();
        assert!(complete && b.converged);
        assert!(b.lower <= q("0.001") && b.upper >= q("0.001"));
        assert!(b.upper - b.lower <= width());
    }
    #[test]
    fn endpoint_distance_minimizes_over_the_complete_spline() {
        let edge = spline();
        let mut dec = Decider::default();
        let b = point_curve(&[q("-0.000001"), q("0"), q("0")], &edge.curve, &mut dec).unwrap();
        assert!(b.converged && b.lower <= q("0.000001") && b.upper >= q("0.000001"));
    }
    #[test]
    fn analytic_endpoint_gap_is_a_distance_not_a_squared_residual() {
        let curve = Curve::Line([q("0"), q("0"), q("0")], [q("2"), q("0"), q("0")]);
        let b = point_curve(
            &[q("1"), q("0.00001"), q("0")],
            &curve,
            &mut Decider::default(),
        )
        .unwrap();
        assert!(b.lower <= q("0.00001") && b.upper >= q("0.00001"));
    }

    #[test]
    fn almost_incident_torus_parallel_has_a_certified_carrier_bound() {
        let placement = Placement {
            origin: [q("0"), q("0"), q("0")],
            axis: [q("0"), q("0"), q("1")],
            ref_dir: [q("1"), q("0"), q("0")],
        };
        let face = Face {
            name: "#1".into(),
            surface: Surface::Torus(placement.clone(), q("5"), q("1")),
            same_sense: true,
            loops: Vec::new(),
        };
        let mut edge = spline();
        edge.curve = Curve::Circle(placement, q("6.000001"));
        let (result, work) = catallax_numbers::with_operation_budget(20_000_000, || {
            edge_surface(&edge, &face, &mut Decider::default())
        });
        assert!(!work.exhausted);
        let (b, complete) = result.unwrap();
        assert!(!complete);
        assert!(b.lower <= q("0.000001") && b.upper >= q("0.000001"));
        assert!(b.upper <= q("0.01"));
    }
}
