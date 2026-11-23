//! Continuous distance bounds for the fitted sweep comparison.

use crate::Options;
use catallax_features::{BuildError, general_sweep, numerical_embedding};
use catallax_field::fitted_motion::{Chart, Motion, SectionBox};
use catallax_numbers::{Decider, Interval as I, OperationSpan, Rational, Real, Refusal};
use catallax_recipe::{Operation, Origin, Recipe, Role};
use catallax_shapes::{BSplineSurface, Surface};
use catallax_topology::numerical::{self, NumericalOptions};
use std::collections::BTreeMap;

type Vec3 = [I; 3];
type Domain = [Rational; 2];

fn station(t: &Rational) -> Domain {
    [t.clone(), t.clone()]
}

/// Distance bounds over a complete fitted sweep and its display mesh.
#[derive(Debug, PartialEq, Eq)]
pub struct FittedSweepReport {
    /// Number of fitted faces covered by the authored law.
    pub faces: usize,
    /// Number of positions in the certified display mesh.
    pub display_vertices: usize,
    /// Maximum authored-to-fit bound in either direction, in millimetres.
    pub source_upper_mm: Rational,
    /// Maximum fitted-surface-to-display bound, in millimetres.
    pub display_upper_mm: Rational,
    /// Maximum sum of source and display bounds for any one face.
    pub composed_upper_mm: Rational,
}

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}

fn number(n: i64) -> I {
    I::from_rational(&q(n))
}

fn correspondence() -> Refusal {
    Refusal::new(
        "compare.sweep-correspondence",
        "the fitted boundary has no complete authored chart correspondence",
    )
}

fn tolerance() -> Refusal {
    Refusal::new(
        "compare.sweep-source-distance",
        "the continuous source and display bound could not be certified at the requested tolerance",
    )
}

fn build_error(error: BuildError) -> Refusal {
    match error {
        BuildError::Geometry { refusal, .. } => refusal,
        other => Refusal::new("compare.sweep-recipe", other.to_string()),
    }
}

fn norm(v: Vec3) -> I {
    v.into_iter()
        .fold(number(0), |sum, x| sum + x.square())
        .sqrt()
        .expect("sum of squares")
}

fn difference(a: Vec3, b: Vec3) -> Vec3 {
    std::array::from_fn(|i| a[i] - b[i])
}

fn upper(x: I) -> Result<Rational, Refusal> {
    Rational::from_f64_exact(x.hi()).ok_or_else(tolerance)
}

fn rational_interval(real: &Real) -> Result<I, Refusal> {
    real.rational_ref()
        .map(I::from_rational)
        .ok_or_else(correspondence)
}

/// Compare an isolated fitted sweep under a single caller-owned work allowance.
/// Authored source and fitted chart bounds cover the full boundary in both
/// directions. Success certifies geometric distance at the requested accuracy.
pub fn check_fitted_sweep(
    recipe: &Recipe,
    modelling: &NumericalOptions,
    options: &Options,
    dec: &mut Decider,
) -> Result<FittedSweepReport, Refusal> {
    recipe
        .check()
        .map_err(|e| Refusal::new("compare.sweep-recipe", e.to_string()))?;
    let step = recipe.step(&recipe.result).ok_or_else(correspondence)?;
    if !matches!(step.op, Operation::Sweep { .. })
        || recipe.steps.iter().any(|s| {
            s.id != step.id && !matches!(s.op, Operation::Sketch { .. } | Operation::Path { .. })
        })
        || !recipe.materials.is_empty()
        || options.region.is_some()
        || options.accuracy <= q(0)
    {
        return Err(Refusal::new(
            "compare.sweep-scope",
            "fitted comparison requires one complete sweep with sketch and path inputs",
        ));
    }

    let (_, authored) = general_sweep::authored_definition(recipe, step, &modelling.tolerance_mm)?;
    let candidate =
        general_sweep::prepare_recipe(recipe, modelling.clone(), dec).map_err(build_error)?;
    let (_, (model, body)) = candidate.admit(
        |model, body, options, dec| {
            let (definition, trims) = numerical_embedding::definition_from_exact_model(
                &model,
                body,
                recipe.name.clone(),
                &options.tolerance_mm,
                dec,
            )?;
            let admitted =
                numerical::build_with_native_trims(definition, trims, options.clone(), dec)?;
            Ok((admitted, (model, body)))
        },
        dec,
    )?;

    let span = OperationSpan::begin("field", &recipe.result, "authored-sweep", dec);
    let source = Motion::new(&authored, dec);
    span.finish(dec, source.as_ref().err(), None);
    let source = source?;

    let span = OperationSpan::begin("brep", &recipe.result, "fitted-sweep-display", dec);
    let display = catallax_mesh::mesh_body(&model, body, &options.mesh, dec);
    span.finish(
        dec,
        display.as_ref().err(),
        display.as_ref().ok().map(|x| &x.marks),
    );
    let display = display?;
    if !display.marks.is_empty() || !display.value.unavailable_controls.is_empty() {
        return Err(Refusal::new(
            "compare.sweep-display",
            "the fitted display carries an unproved control",
        ));
    }
    let display = display.value;
    if display.faces.len() != model.faces_of_body(body).len() || display.positions.is_empty() {
        return Err(correspondence());
    }

    let span = OperationSpan::begin("compare", &recipe.result, "fitted-sweep-bound", dec);
    let result = (|| {
        let bound = |display: &Rational| {
            modelling
                .tolerance_mm
                .clone()
                .min(&options.accuracy - display)
        };
        let mut report = FittedSweepReport {
            faces: display.faces.len(),
            display_vertices: display.positions.len(),
            source_upper_mm: q(0),
            display_upper_mm: q(0),
            composed_upper_mm: q(0),
        };
        let mut add_face = |source: Rational, chord: &Rational| {
            if source > bound(chord) {
                return Err(tolerance());
            }
            report.composed_upper_mm = report.composed_upper_mm.clone().max(&source + chord);
            report.source_upper_mm = report.source_upper_mm.clone().max(source);
            report.display_upper_mm = report.display_upper_mm.clone().max(chord.clone());
            Ok(())
        };

        let mut sections = Sections {
            source: &source,
            cells: BTreeMap::new(),
        };
        let mut coverage = vec![Vec::<Domain>::new(); source.charts().len()];
        let mut caps = [None, None];
        let mut cap_error = [q(0), q(0)];
        for face in &display.faces {
            dec.charge_work(1)?;
            if face.triangles.is_empty() {
                return Err(correspondence());
            }
            let model_face = model.face(face.face);
            let surface = model.geom().surface(model_face.surface);
            let [Origin { step, role, .. }] = model.face_origins(face.face) else {
                return Err(correspondence());
            };
            if step != &recipe.result {
                return Err(correspondence());
            }
            match (role, surface) {
                (Role::Side { segment, .. }, Surface::BSpline(fit)) => {
                    let carrier = fit.carrier();
                    let [u, v] = carrier.domain();
                    let mut matches = source
                        .charts()
                        .iter()
                        .enumerate()
                        .filter(|(_, c)| &c.segment == segment && c.domain.each_ref() == u);
                    let (index, chart) = matches.next().ok_or_else(correspondence)?;
                    if matches.next().is_some() {
                        return Err(correspondence());
                    }
                    let v = v.map(Clone::clone);
                    let proof = FaceBound::new(chart, carrier, &v)?;
                    let (maximum, ends) =
                        proof.check(&mut sections, &bound(&face.max_chord_error), dec)?;
                    for k in 0..2 {
                        if v[k] == q(k as i64) {
                            cap_error[k] = cap_error[k].clone().max(ends[k].clone());
                        }
                    }
                    coverage[index].push(v);
                    add_face(maximum, &face.max_chord_error)?;
                }
                (Role::Start | Role::End, Surface::Plane(plane)) => {
                    let index = usize::from(matches!(role, Role::End));
                    if caps[index]
                        .replace((plane, &face.max_chord_error))
                        .is_some()
                    {
                        return Err(correspondence());
                    }
                }
                _ => return Err(correspondence()),
            }
        }
        for mut spans in coverage {
            spans.sort_unstable();
            if spans.first().map(|v| &v[0]) != Some(&q(0))
                || spans.last().map(|v| &v[1]) != Some(&q(1))
                || spans.windows(2).any(|pair| pair[0][1] != pair[1][0])
            {
                return Err(correspondence());
            }
        }

        // The fitter constructs the cap trim loops from the shared end rails.
        // Orthogonally project a fitted cap onto the source plane. The boundary
        // moves by at most the rail bound under this projection. Winding is
        // constant along the boundary homotopy away from its swept strip, so
        // the projected cap differs only within that bound of either outline.
        // The omitted height attains its extremes on the cap's boundary.
        for (i, cap) in caps.into_iter().enumerate() {
            let (plane, chord) = cap.ok_or_else(correspondence)?;
            let end = station(&q(i as i64));
            let normal = sections.at(&end, dec)?.normal();
            let [x, y, z] =
                [&plane.normal.x, &plane.normal.y, &plane.normal.z].map(rational_interval);
            let dot = normal[0] * x? + normal[1] * y? + normal[2] * z?;
            if dot.contains_zero() {
                return Err(correspondence());
            }
            let bound = I::from_rational(&cap_error[i]) * number(2).sqrt().expect("positive");
            add_face(upper(bound)?, chord)?;
        }
        Ok(report)
    })();
    span.finish(dec, result.as_ref().err(), None);
    result
}

struct Sections<'a> {
    source: &'a Motion,
    cells: BTreeMap<Domain, SectionBox>,
}

impl Sections<'_> {
    fn at(&mut self, range: &Domain, dec: &mut Decider) -> Result<&SectionBox, Refusal> {
        use std::collections::btree_map::Entry;
        // Keep room for two endpoints and one interior enclosure per
        // binary64 subdivision level. Eviction affects work, not the bound.
        if self.cells.len() >= 3 * f64::MANTISSA_DIGITS as usize {
            self.cells.clear();
        }
        Ok(match self.cells.entry(range.clone()) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => entry.insert(self.source.section(range.each_ref(), dec)?),
        })
    }
}

struct FaceBound<'a> {
    domain: &'a Domain,
    chart: &'a Chart,
    rails: Vec<Vec<Vec3>>,
    weight_error: I,
}

impl<'a> FaceBound<'a> {
    fn new(chart: &'a Chart, fit: &BSplineSurface, domain: &'a Domain) -> Result<Self, Refusal> {
        let [u, v] = fit.control_count();
        if fit.degrees() != [u - 1, v - 1]
            || u != chart.controls.len()
            || !matches!(v, 2 | 4)
            || domain[0] >= domain[1]
        {
            return Err(correspondence());
        }
        let mut delta = number(0);
        for (row, source) in fit.weights().chunks_exact(v).zip(&chart.weights) {
            if source.lo() <= 0.0
                || row
                    .iter()
                    .any(|w| w.rational_ref() != row[0].rational_ref())
            {
                return Err(correspondence());
            }
            let relative = (rational_interval(&row[0])? - *source) / *source;
            delta = delta.hull(&norm([relative, number(0), number(0)]));
        }
        if delta.hi() >= 1.0 {
            return Err(correspondence());
        }
        let mut diameter = number(0);
        for (i, &a) in chart.controls.iter().enumerate() {
            for &b in &chart.controls[i + 1..] {
                diameter = diameter.hull(&norm(difference(a, b)));
            }
        }
        // If a positive chart weight changes by a relative e_i, the convex
        // barycentric coordinate is scaled by (1+e_i)/(1+mean(e)). Subtract
        // the source mean before summing to bound the effect by the diameter.
        let weight_error = delta * diameter / (number(1) - delta);
        let rails = fit
            .controls()
            .chunks_exact(v)
            .map(|row| {
                row.iter()
                    .map(|p| {
                        Ok([
                            rational_interval(&p.x)?,
                            rational_interval(&p.y)?,
                            rational_interval(&p.z)?,
                        ])
                    })
                    .collect()
            })
            .collect::<Result<_, Refusal>>()?;
        Ok(Self {
            domain,
            chart,
            rails,
            weight_error,
        })
    }

    fn check(
        self,
        sections: &mut Sections<'_>,
        allowed: &Rational,
        dec: &mut Decider,
    ) -> Result<(Rational, [Rational; 2]), Refusal> {
        let mut maximum = q(0);
        let mut ends = [q(0), q(0)];
        for (rail, &control) in self.rails.into_iter().zip(&self.chart.controls) {
            for i in 0..2 {
                let source = sections
                    .at(&station(&self.domain[i]), dec)?
                    .control(control)
                    .0;
                let fitted = rail[if i == 0 { 0 } else { rail.len() - 1 }];
                ends[i] = ends[i]
                    .clone()
                    .max(upper(norm(difference(source, fitted)) + self.weight_error)?);
            }
            let mut pending = vec![(self.domain.clone(), rail, 0)];
            while let Some((range, curve, depth)) = pending.pop() {
                dec.charge_work(1)?;
                let a = sections.at(&station(&range[0]), dec)?.control(control).0;
                let b = sections.at(&station(&range[1]), dec)?.control(control).0;
                let a = norm(difference(a, curve[0]));
                let b = norm(difference(b, curve[curve.len() - 1]));
                let error = I::new(a.lo().max(b.lo()), a.hi().max(b.hi())).expect("ordered maxima");
                if Rational::from_f64_exact(error.lo()).is_some_and(|lo| lo > *allowed) {
                    return Err(tolerance());
                }
                match sections.at(&range, dec) {
                    Ok(section) => {
                        let (_, source_accel) = section.control(control);
                        let duration = I::from_rational(&(&range[1] - &range[0])).square();
                        let fit_accel = cubic_acceleration(&curve);
                        let accel = norm(difference(source_accel.map(|x| x * duration), fit_accel));
                        // The interpolation remainder on the unit interval
                        // is bounded by sup |E''| * max(s(1-s))/2.
                        let bound = upper(error + accel / number(8) + self.weight_error)?;
                        if bound <= *allowed {
                            maximum = maximum.max(bound);
                            continue;
                        }
                    }
                    Err(e) if e.code == "field.sweep-guide-enclosure" => (),
                    Err(e) => return Err(e),
                }
                // The interval law works on a binary64 enclosure grid.
                if depth >= f64::MANTISSA_DIGITS - 1 {
                    return Err(tolerance());
                }
                let mid = (&range[0] + &range[1]).mul_pow2(-1);
                let [left, right] = halve(curve);
                let [a, b] = range;
                pending.push(([mid.clone(), b], right, depth + 1));
                pending.push(([a, mid], left, depth + 1));
            }
        }
        Ok((maximum, ends))
    }
}

fn cubic_acceleration(curve: &[Vec3]) -> Vec3 {
    if let [a, b, c, d] = curve {
        std::array::from_fn(|i| {
            // Two Bernstein controls of the second derivative of a cubic.
            let first = (a[i] - number(2) * b[i] + c[i]) * number(6);
            let last = (b[i] - number(2) * c[i] + d[i]) * number(6);
            first.hull(&last)
        })
    } else {
        [number(0); 3]
    }
}

fn halve(mut curve: Vec<Vec3>) -> [Vec<Vec3>; 2] {
    let mut left = Vec::new();
    let mut right = Vec::new();
    while let Some(&last) = curve.last() {
        left.push(curve[0]);
        right.push(last);
        for i in 0..curve.len() - 1 {
            curve[i] = std::array::from_fn(|k| (curve[i][k] + curve[i + 1][k]) / number(2));
        }
        curve.pop();
    }
    right.reverse();
    [left, right]
}
