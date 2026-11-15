//! Separate approximate STEP artifacts for checked numerical bodies.
//! No exact answer, complete bound, model extraction or healing is provided.
use crate::graph::{Field, Graph, Id, boolean, refs, string, text};
use catallax_numbers::{
    Decider, NumericalEvidence, Rational, Real, Refusal, Sign, operation_budget_active,
    with_operation_budget,
};
use catallax_shapes::{Curve, Point3, Surface, Vector3};
use catallax_topology::{
    Sense,
    numerical::{DEFAULT_WORK, NumericalBody},
};

mod analytic;
mod curved;

/// Controls fixed before a numerical export. They do not change source tolerance.
#[derive(Clone, Debug)]
pub struct NumericalStepOptions {
    /// STEP product name, independent of the retained source name.
    pub product_name: String,
    /// Significant digits, from 1 through 17 in this first writer.
    pub significant_digits: u32,
    /// Maximum decimal error of each written length scalar, in millimetres.
    /// This is not a whole-body approximation bound.
    pub max_scalar_rounding_error_mm: Rational,
}
impl Default for NumericalStepOptions {
    fn default() -> Self {
        Self {
            product_name: "catallax numerical approximate".into(),
            significant_digits: 17,
            max_scalar_rounding_error_mm: Rational::from_ratio(1, 10_000_000).unwrap(),
        }
    }
}
/// Rounding of native spline parameters and weights, separate from lengths.
#[derive(Clone, Debug, Default)]
pub struct NumericalSplineRounding {
    /// Maximum error in a knot, finite curve parameter or patch UV coordinate,
    /// measured in that native parameter's units. No spatial bound is inferred.
    pub max_parameter_error: Rational,
    /// Maximum error of a dimensionless positive control weight.
    pub max_weight_error: Rational,
}
/// Source and file checks with separate scopes; always a numerical result.
#[derive(Clone, Debug)]
pub struct NumericalStepReport {
    /// Original generating input, before any explicit sewing or later edits.
    pub original_definition: catallax_topology::numerical::NumericalDefinition,
    /// Native parameter curves in the original generating input.
    pub original_native_trims: Vec<catallax_topology::numerical::NumericalNativeTrim>,
    /// Native parameter curves checked on the current source.
    pub source_native_trims: Vec<catallax_topology::numerical::NumericalNativeTrim>,
    /// Native parameter curves reconstructed from the emitted decimal values.
    pub emitted_native_trims: Vec<catallax_topology::numerical::NumericalNativeTrim>,
    /// Earlier checked states. Their checks apply only to their own definitions.
    pub source_history: Vec<catallax_topology::numerical::NumericalHistory>,
    /// Retained generating definition name.
    pub source_name: String,
    /// Immutable copy of the generating tables, trims and stable source keys.
    pub source_definition: catallax_topology::numerical::NumericalDefinition,
    /// N01/N02 evidence, including tolerance, actual samples and repair history.
    /// Sample maxima are never promoted into a complete bound.
    pub source_evidence: NumericalEvidence,
    /// N01/N02 checks on the emitted decimal-value model with retained topology.
    /// This is separate from the source checks and external reference validation.
    pub emitted_evidence: NumericalEvidence,
    /// Requested error limit on written length scalars only.
    pub scalar_rounding_limit_mm: Rational,
    /// Actual maximum rounding error among written length scalars only.
    pub max_scalar_rounding_error_mm: Rational,
    /// Maximum rounding error of normalized DIRECTION ratios, dimensionless.
    /// This is not an angular or spatial deviation bound.
    pub max_direction_ratio_rounding_error: Rational,
    /// Stored faces, explicit edges and vertices, with no inferred topology.
    pub counts: [usize; 3],
    /// Native spline rounding, when the curved writer is used.
    pub spline_rounding: Option<NumericalSplineRounding>,
    /// Checks performed by this writer, distinct from external reference checks.
    pub checks: Vec<String>,
}
/// An approximate file, with no implicit conversion to a proven/bounded result.
/// The private bytes were checked before this value was returned.
/// Raw reference import and edit acceptance require an independent check.
///
/// ```compile_fail
/// use catallax_step::{StepReport, numerical::NumericalStep};
/// use catallax_numbers::Answer;
/// fn promote(file: NumericalStep) -> Answer<StepReport> { file.into() }
/// ```
#[derive(Clone, Debug)]
pub struct NumericalStep {
    bytes: Vec<u8>,
    report: NumericalStepReport,
}
impl NumericalStep {
    /// Approximate STEP bytes; the file itself states the numerical category.
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
    /// Source provenance and scoped file checks, never an exact STEP report.
    pub fn report(&self) -> &NumericalStepReport {
        &self.report
    }
    /// Checks bytes reread from the actual output, including complete byte identity
    /// to this checked artifact and the seam syntax contract.
    /// This is not an OCCT import, healing, geometric certification or proof.
    pub fn check_written(&self, reread: &[u8], dec: &mut Decider) -> Result<(), Refusal> {
        guarded(|| {
            if reread.len() != self.bytes.len() {
                return Err(bad("file-changed", "reread file length differs"));
            }
            dec.charge_work(reread.len() as u64)?;
            if reread != self.bytes {
                return Err(bad("file-changed", "reread file bytes differ"));
            }
            crate::check_seams(reread)?;
            dec.check_deadline()?;
            Ok(())
        })
    }
}
fn bad(kind: &str, detail: impl Into<String>) -> Refusal {
    let code = match kind {
        "file-changed" => "step.numerical-file-changed",
        "file-number" => "step.numerical-file-number",
        "rounding-topology" => "step.numerical-rounding-topology",
        "rounding-limit" => "step.numerical-rounding-limit",
        "unsupported" => "step.numerical-unsupported",
        "options" => "step.numerical-options",
        "gap" => "step.numerical-gap",
        _ => unreachable!(),
    };
    Refusal::new(code, detail)
}
fn guarded<T>(run: impl FnOnce() -> Result<T, Refusal>) -> Result<T, Refusal> {
    if operation_budget_active() {
        return run();
    }
    let (result, work) = with_operation_budget(DEFAULT_WORK, run);
    if work.exhausted {
        Err(Refusal::new(
            "numbers.operation-budget",
            "numerical STEP exhausted its work allowance",
        ))
    } else {
        result
    }
}
fn rational_point(p: &Point3) -> Result<[Rational; 3], Refusal> {
    [&p.x, &p.y, &p.z]
        .map(|v| {
            v.as_rational().ok_or_else(|| {
                bad(
                    "unsupported",
                    "first numerical writer needs rational coordinates",
                )
            })
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| unreachable!())
}
fn rational_vector(p: &Vector3) -> Result<[Rational; 3], Refusal> {
    rational_point(&Point3::new(p.x.clone(), p.y.clone(), p.z.clone()))
}
fn axis(v: &[Rational; 3], unit: bool) -> Option<usize> {
    let nonzero: Vec<_> = v.iter().enumerate().filter(|(_, q)| !q.is_zero()).collect();
    if nonzero.len() != 1 || (unit && nonzero[0].1.abs() != Rational::one()) {
        return None;
    }
    Some(nonzero[0].0)
}
fn cross(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}
// Parse the selected STEP decimal token independently. Rational::parse needs a
// digit after '.', while STEP permits a trailing decimal point before E/end.
fn written(q: &Rational, options: &NumericalStepOptions) -> Result<Rational, Refusal> {
    let mut token = crate::number::rational(q, options.significant_digits)
        .text
        .replace(".E", ".0E");
    if token.ends_with('.') {
        token.push('0');
    }
    Rational::parse(&token).map_err(|_| bad("file-number", "written decimal cannot be reread"))
}
fn point_from(p: &[Rational; 3]) -> Point3 {
    Point3::new(
        p[0].clone().into(),
        p[1].clone().into(),
        p[2].clone().into(),
    )
}

fn arithmetic(
    a: &Rational,
    b: &Rational,
    op: char,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    dec.reserve_exact_arithmetic(a, b, matches!(op, '*' | '/'))?;
    Ok(match op {
        '+' => a + b,
        '-' => a - b,
        '*' => a * b,
        '/' => a / b,
        _ => unreachable!(),
    })
}
fn dot(a: &[Rational; 3], b: &[Rational; 3], dec: &mut Decider) -> Result<Rational, Refusal> {
    let mut result = Rational::zero();
    for i in 0..3 {
        let product = arithmetic(&a[i], &b[i], '*', dec)?;
        result = arithmetic(&result, &product, '+', dec)?;
    }
    Ok(result)
}
fn difference(
    a: &[Rational; 3],
    b: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    Ok([
        arithmetic(&a[0], &b[0], '-', dec)?,
        arithmetic(&a[1], &b[1], '-', dec)?,
        arithmetic(&a[2], &b[2], '-', dec)?,
    ])
}
fn charged_cross(
    a: &[Rational; 3],
    b: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 3], Refusal> {
    let mut result = std::array::from_fn(|_| Rational::zero());
    for (i, q) in result.iter_mut().enumerate() {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        let x = arithmetic(&a[j], &b[k], '*', dec)?;
        let y = arithmetic(&a[k], &b[j], '*', dec)?;
        *q = arithmetic(&x, &y, '-', dec)?;
    }
    Ok(result)
}
struct DirectionValues {
    values: [Rational; 3],
    tokens: [String; 3],
    error: Rational,
}
// STEP DIRECTION ratios are dimensionless and normalized by the reader. Scale
// them before decimal selection; never put their error in the length report.
fn direction_values(
    v: &[Rational; 3],
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<DirectionValues, Refusal> {
    if axis(v, false).is_some() {
        let values = v.each_ref().map(|q| {
            if q.is_zero() {
                Rational::zero()
            } else if q.sign() == Sign::Positive {
                Rational::one()
            } else {
                -Rational::one()
            }
        });
        let tokens = values.each_ref().map(|q| {
            if q.is_zero() {
                "0."
            } else if q.sign() == Sign::Positive {
                "1."
            } else {
                "-1."
            }
            .to_owned()
        });
        return Ok(DirectionValues {
            values,
            tokens,
            error: Rational::zero(),
        });
    }
    let mut scale = Rational::zero();
    for q in v {
        dec.reserve_exact_arithmetic(q, &scale, true)?;
        scale = scale.max(q.abs());
    }
    if scale.is_zero() {
        return Err(bad("rounding-topology", "direction vanished"));
    }
    let mut result = DirectionValues {
        values: std::array::from_fn(|_| Rational::zero()),
        tokens: std::array::from_fn(|_| String::new()),
        error: Rational::zero(),
    };
    for (i, q) in v.iter().enumerate() {
        let ratio = arithmetic(q, &scale, '/', dec)?;
        dec.charge_work(1)?;
        let number = crate::number::rational(&ratio, options.significant_digits);
        dec.reserve_exact_arithmetic(&result.error, &number.error, true)?;
        result.error = result.error.max(number.error);
        result.values[i] = written(&ratio, options)?;
        result.tokens[i] = number.text;
    }
    Ok(result)
}
fn inverse_uv(
    plane: &catallax_shapes::Plane,
    p: &[Rational; 3],
    dec: &mut Decider,
) -> Result<[Rational; 2], Refusal> {
    let origin = rational_point(&plane.origin)?;
    let x = rational_vector(&plane.ref_dir)?;
    let y = charged_cross(&rational_vector(&plane.normal)?, &x, dec)?;
    let delta = difference(p, &origin, dec)?;
    let mut uv = [Rational::zero(), Rational::zero()];
    for (i, axis) in [&x, &y].into_iter().enumerate() {
        let squared = dot(axis, axis, dec)?;
        if squared.sign() != Sign::Positive {
            return Err(bad("rounding-topology", "written plane frame vanished"));
        }
        let value = &Real::from(dot(&delta, axis, dec)?) / &Real::from(squared).sqrt();
        uv[i] = dec.enclose_budgeted_at(&value, 128)?.mid();
    }
    Ok(uv)
}
// Same IDs/loops, and only values actually selected for serialization. The UVs
// below are inverse coordinates on emitted planar frames, not inferred topology.
fn emitted_checks(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalEvidence, Refusal> {
    let mut d = body.definition().clone();
    d.name.push_str(" STEP decimal payload");
    for vertex in &mut d.vertices {
        dec.charge_work(1)?;
        for q in &mut vertex.position {
            *q = written(q, options)?;
        }
    }
    for edge in &mut d.edges {
        dec.charge_work(1)?;
        let Curve::Line(line) = &mut edge.curve else {
            return Err(bad("unsupported", "finite line required"));
        };
        let mut origin = rational_point(&line.origin)?;
        let mut dir = rational_vector(&line.dir)?;
        for q in &mut origin {
            *q = written(q, options)?;
        }
        if let Some(a) = axis(&dir, false) {
            let magnitude = written(&dir[a].abs(), options)?;
            dir[a] = if dir[a].sign() == Sign::Negative {
                -magnitude
            } else {
                magnitude
            };
            if dir[a].is_zero() {
                return Err(bad("rounding-topology", "line direction vanished"));
            }
            edge.parameters = std::array::from_fn(|i| {
                (&d.vertices[edge.vertices[i]].position[a] - &origin[a]) / dir[a].clone()
            });
        } else {
            // Unit magnitude in STEP and these retained rational ratios define
            // the same oriented line carrier. Its finite parameters are local
            // inverse coordinates; EDGE_CURVE keeps the explicit vertex IDs.
            dir = direction_values(&dir, options, dec)?.values;
            let squared = dot(&dir, &dir, dec)?;
            for i in 0..2 {
                let delta = difference(&d.vertices[edge.vertices[i]].position, &origin, dec)?;
                edge.parameters[i] = arithmetic(&dot(&delta, &dir, dec)?, &squared, '/', dec)?;
            }
        }
        line.origin = point_from(&origin);
        line.dir = Vector3::new(
            dir[0].clone().into(),
            dir[1].clone().into(),
            dir[2].clone().into(),
        );
    }
    for face in &mut d.faces {
        dec.charge_work(1)?;
        let Surface::Plane(plane) = &mut face.surface else {
            return Err(bad("unsupported", "plane required"));
        };
        let mut origin = rational_point(&plane.origin)?;
        for q in &mut origin {
            *q = written(q, options)?;
        }
        plane.origin = point_from(&origin);
        let mut x = rational_vector(&plane.ref_dir)?;
        let normal = rational_vector(&plane.normal)?;
        let coordinate_frame = axis(&normal, true).is_some() && axis(&x, true).is_some();
        if !coordinate_frame {
            let normal = direction_values(&normal, options, dec)?.values;
            let reference = direction_values(&x, options, dec)?.values;
            // gp_Ax2 projects the written reference onto the written normal's
            // perpendicular plane. Reconstruct that actual STEP frame.
            let y = charged_cross(&normal, &reference, dec)?;
            x = charged_cross(&y, &normal, dec)?;
            if dot(&x, &x, dec)?.is_zero() {
                return Err(bad("rounding-topology", "written plane axes are parallel"));
            }
            plane.normal = Vector3::new(
                normal[0].clone().into(),
                normal[1].clone().into(),
                normal[2].clone().into(),
            );
            plane.ref_dir = Vector3::new(
                x[0].clone().into(),
                x[1].clone().into(),
                x[2].clone().into(),
            );
        }
        let y = cross(&rational_vector(&plane.normal)?, &x);
        for trim in face.loops.iter_mut().flatten() {
            let edge = &d.edges[trim.edge];
            let endpoints = if trim.sense == Sense::Same {
                edge.vertices
            } else {
                [edge.vertices[1], edge.vertices[0]]
            };
            for (i, endpoint) in endpoints.into_iter().enumerate() {
                let p = &d.vertices[endpoint].position;
                trim.uv[i] = if coordinate_frame {
                    [&x, &y].map(|axis| {
                        (0..3).fold(Rational::zero(), |sum, k| {
                            sum + (&p[k] - &origin[k]) * axis[k].clone()
                        })
                    })
                } else {
                    inverse_uv(plane, p, dec)?
                };
            }
        }
    }
    Ok(
        catallax_topology::numerical::build(d, body.options().clone(), dec)?
            .evidence()
            .clone(),
    )
}

struct Writer<'a> {
    graph: Graph,
    options: &'a NumericalStepOptions,
    dec: &'a mut Decider,
    error: Rational,
    direction_error: Rational,
}
impl Writer<'_> {
    fn scalar(&mut self, q: &Rational) -> Result<Field, Refusal> {
        self.dec.charge_work(1)?;
        let n = crate::number::rational(q, self.options.significant_digits);
        if n.error > self.options.max_scalar_rounding_error_mm {
            return Err(bad(
                "rounding-limit",
                "written scalar exceeds its separate rounding limit",
            ));
        }
        self.error = self.error.clone().max(n.error);
        Ok(text(n.text))
    }
    fn point(&mut self, q: &[Rational; 3]) -> Result<Id, Refusal> {
        let xyz = q
            .iter()
            .map(|q| self.scalar(q))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self
            .graph
            .add("CARTESIAN_POINT", vec![text("''"), Field::List(xyz)]))
    }
    fn direction(&mut self, q: &[Rational; 3]) -> Result<Id, Refusal> {
        let values = direction_values(q, self.options, self.dec)?;
        self.direction_error = self.direction_error.clone().max(values.error);
        Ok(self.graph.add(
            "DIRECTION",
            vec![
                text("''"),
                Field::List(values.tokens.into_iter().map(text).collect()),
            ],
        ))
    }
}
/// Builds a separate approximate STEP artifact for a checked numerical body.
/// The planar line family retains its rational incidence and oblique-frame
/// checks. Native spline edges and rational patches recheck the decimal payload
/// numerically under the body's stated tolerance, with all physical loops.
/// Unsupported carriers refuse. No topology is inferred and no vertex is healed.
/// The source tolerance is metadata, not a STEP global uncertainty or body bound.
/// The caller's allowance covers all work; no global policy is changed.
pub fn build(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    guarded(|| {
        if analytic::applies(body) {
            analytic::build(body, options, dec)
        } else {
            build_inner(body, options, dec)
        }
    })
}
fn build_inner(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    if options.product_name.trim().is_empty()
        || !(1..=17).contains(&options.significant_digits)
        || options.max_scalar_rounding_error_mm.sign() != Sign::Positive
    {
        return Err(bad(
            "options",
            "named product, 1..17 digits and positive scalar rounding limit required",
        ));
    }
    let d = body.definition();
    let count = body.counts();
    if count.shells != 1 || count.lumps != 1 {
        return Err(bad("unsupported", "one retained closed shell required"));
    }
    if !body.native_trims().is_empty()
        || d.edges
            .iter()
            .any(|edge| !matches!(edge.curve, Curve::Line(_)))
        || d.faces
            .iter()
            .any(|face| !matches!(face.surface, Surface::Plane(_)))
    {
        return curved::build(body, options, dec);
    }
    let mut w = Writer {
        graph: Graph::default(),
        options,
        dec,
        error: Rational::zero(),
        direction_error: Rational::zero(),
    };
    let mut vertices = Vec::new();
    for vertex in &d.vertices {
        let point = w.point(&vertex.position)?;
        vertices.push(w.graph.add(
            "VERTEX_POINT",
            vec![text(string(&vertex.key)), point.into()],
        ));
    }
    let mut edges = Vec::new();
    for edge in &d.edges {
        w.dec.charge_work(1)?;
        let Curve::Line(line) = &edge.curve else {
            return Err(bad(
                "unsupported",
                "only retained finite line edges are supported",
            ));
        };
        let origin = rational_point(&line.origin)?;
        let dir = rational_vector(&line.dir)?;
        for (i, t) in edge.parameters.iter().enumerate() {
            let value = std::array::from_fn::<_, 3, _>(|k| &origin[k] + &(&dir[k] * t));
            if value != d.vertices[edge.vertices[i]].position {
                return Err(bad(
                    "gap",
                    "retained line endpoint and vertex differ; no export healing",
                ));
            }
        }
        let point = w.point(&origin)?;
        let direction = w.direction(&dir)?;
        let magnitude = if let Some(a) = axis(&dir, false) {
            w.scalar(&dir[a].abs())?
        } else {
            w.scalar(&Rational::one())?
        };
        let vector = w
            .graph
            .add("VECTOR", vec![text("''"), direction.into(), magnitude]);
        let line = w.graph.add(
            "LINE",
            vec![text(string(&edge.key)), point.into(), vector.into()],
        );
        edges.push(w.graph.add(
            "EDGE_CURVE",
            vec![
                text(string(&edge.key)),
                vertices[edge.vertices[0]].into(),
                vertices[edge.vertices[1]].into(),
                line.into(),
                boolean(edge.parameters[0] < edge.parameters[1]),
            ],
        ));
    }
    let mut faces = Vec::new();
    for face in &d.faces {
        w.dec.charge_work(1)?;
        let Surface::Plane(plane) = &face.surface else {
            return Err(bad(
                "unsupported",
                "numerical curved-face export is outside this bounded family",
            ));
        };
        let origin = rational_point(&plane.origin)?;
        let normal = rational_vector(&plane.normal)?;
        let x = rational_vector(&plane.ref_dir)?;
        let coordinate_frame = axis(&normal, true).is_some() && axis(&x, true).is_some();
        if coordinate_frame && axis(&normal, true) == axis(&x, true) {
            return Err(bad("unsupported", "plane axes must be orthogonal"));
        }
        let point = w.point(&origin)?;
        let n = w.direction(&normal)?;
        let reference = w.direction(&x)?;
        let placement = w.graph.add(
            "AXIS2_PLACEMENT_3D",
            vec![text("''"), point.into(), n.into(), reference.into()],
        );
        let surface = w
            .graph
            .add("PLANE", vec![text(string(&face.key)), placement.into()]);
        let mut bounds = Vec::new();
        for (index, boundary) in face.loops.iter().enumerate() {
            let mut uses = Vec::new();
            for trim in boundary {
                w.dec.charge_work(1)?;
                let edge = &d.edges[trim.edge];
                let endpoint = if trim.sense == Sense::Same {
                    edge.vertices
                } else {
                    [edge.vertices[1], edge.vertices[0]]
                };
                for endpoint in endpoint {
                    // STEP retains these world-space boundary vertices. Its
                    // implicit planar parameters need not equal the source's
                    // numerical UV midpoints; no point or carrier is moved.
                    let position = &d.vertices[endpoint].position;
                    let delta = difference(position, &origin, w.dec)?;
                    if !dot(&normal, &delta, w.dec)?.is_zero() {
                        return Err(bad(
                            "gap",
                            "retained boundary vertex is off its plane; no export healing",
                        ));
                    }
                }
                uses.push(w.graph.add(
                    "ORIENTED_EDGE",
                    vec![
                        text("''"),
                        text("*"),
                        text("*"),
                        edges[trim.edge].into(),
                        boolean(trim.sense == Sense::Same),
                    ],
                ));
            }
            let loop_ = w.graph.add("EDGE_LOOP", vec![text("''"), refs(uses)]);
            bounds.push(w.graph.add(
                if index == 0 {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                },
                vec![text("''"), loop_.into(), boolean(true)],
            ));
        }
        faces.push(w.graph.add(
            "ADVANCED_FACE",
            vec![
                text(string(&face.key)),
                refs(bounds),
                surface.into(),
                boolean(face.sense == Sense::Same),
            ],
        ));
    }
    let shell = w.graph.add("CLOSED_SHELL", vec![text("''"), refs(faces)]);
    let solid = w.graph.add(
        "MANIFOLD_SOLID_BREP",
        vec![text("'numerical approximate'"), shell.into()],
    );
    let emitted_evidence = emitted_checks(body, options, w.dec)?;
    w.finish(
        solid,
        body,
        emitted_evidence,
        None,
        vec![
            "source.retained-line-endpoints-and-plane-incidence".into(),
            "file.decimal-model-numerical-admission".into(),
            "file.canonical-byte-identity".into(),
            "file.planar-seam-contract".into(),
            "file.approximate-category-no-global-bound".into(),
        ],
    )
}

impl Writer<'_> {
    fn finish(
        self,
        solid: Id,
        body: &NumericalBody,
        emitted_evidence: NumericalEvidence,
        spline_rounding: Option<NumericalSplineRounding>,
        checks: Vec<String>,
    ) -> Result<NumericalStep, Refusal> {
        let options = self.options;
        let d = body.definition();
        let count = body.counts();
        let mut description = format!(
            "catallax numerical approximate;source={};source_tolerance_mm={};source_checks={};source_repairs={};complete_body_bound=unestablished;max_scalar_rounding_error_mm={}",
            d.name,
            body.options().tolerance_mm,
            body.evidence().checks().len(),
            body.evidence().repairs().len(),
            self.error
        );
        if !self.direction_error.is_zero() {
            description.push_str(&format!(
                ";max_direction_ratio_rounding_error={}",
                self.direction_error
            ));
        }
        if let Some(rounding) = &spline_rounding {
            description.push_str(&format!(
                ";max_native_parameter_rounding_error={};max_weight_rounding_error={}",
                rounding.max_parameter_error, rounding.max_weight_error
            ));
        }
        let mut file = self.graph.write(&[solid], |roots| {
            // Reuse AP242 product/unit framing, not its exact-model uncertainty claim.
            crate::export::approximate_header(&options.product_name, roots, &description)
        });
        file.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        let artifact = NumericalStep {
            bytes: file.into_bytes(),
            report: NumericalStepReport {
                original_definition: body.source_definition().clone(),
                original_native_trims: body.source_native_trims().to_vec(),
                source_native_trims: body.native_trims().to_vec(),
                emitted_native_trims: Vec::new(),
                source_history: body.history().to_vec(),
                source_name: d.name.clone(),
                source_definition: d.clone(),
                source_evidence: body.evidence().clone(),
                emitted_evidence,
                scalar_rounding_limit_mm: options.max_scalar_rounding_error_mm.clone(),
                max_scalar_rounding_error_mm: self.error,
                max_direction_ratio_rounding_error: self.direction_error,
                counts: [count.faces, count.edges, count.vertices],
                checks,
                spline_rounding,
            },
        };
        artifact.check_written(artifact.bytes(), self.dec)?;
        Ok(artifact)
    }
}
