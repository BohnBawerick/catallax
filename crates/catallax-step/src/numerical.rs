//! Separate approximate STEP artifacts for checked numerical planar bodies.
//! No exact answer, complete bound, model extraction or healing is provided.
use crate::graph::{Field, Graph, Id, boolean, refs, string, text};
use catallax_numbers::{
    Decider, NumericalEvidence, Rational, Refusal, Sign, operation_budget_active,
    with_operation_budget,
};
use catallax_shapes::{Curve, Point3, Surface, Vector3};
use catallax_topology::{
    Sense,
    numerical::{DEFAULT_WORK, NumericalBody},
};

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
/// Source and file checks with separate scopes; always a numerical result.
#[derive(Clone, Debug)]
pub struct NumericalStepReport {
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
    /// Stored faces, explicit edges and vertices, with no inferred topology.
    pub counts: [usize; 3],
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
    /// to this checked artifact and the existing planar seam syntax contract.
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
fn axis(v: &[Rational; 3], unit: bool) -> Result<usize, Refusal> {
    let nonzero: Vec<_> = v.iter().enumerate().filter(|(_, q)| !q.is_zero()).collect();
    if nonzero.len() != 1 || (unit && nonzero[0].1.abs() != Rational::one()) {
        return Err(bad(
            "unsupported",
            "first numerical writer needs axis-aligned lines and unit orthogonal plane frames",
        ));
    }
    Ok(nonzero[0].0)
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
        let a = axis(&dir, false)?;
        for q in &mut origin {
            *q = written(q, options)?;
        }
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
        let x = rational_vector(&plane.ref_dir)?;
        let y = cross(&rational_vector(&plane.normal)?, &x);
        for trim in face.loops.iter_mut().flatten() {
            let edge = &d.edges[trim.edge];
            let endpoints = if trim.sense == Sense::Same {
                edge.vertices
            } else {
                [edge.vertices[1], edge.vertices[0]]
            };
            trim.uv = std::array::from_fn(|i| {
                let p = &d.vertices[endpoints[i]].position;
                [&x, &y].map(|axis| {
                    (0..3).fold(Rational::zero(), |sum, k| {
                        sum + (&p[k] - &origin[k]) * axis[k].clone()
                    })
                })
            });
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
    fn direction(&mut self, q: &[Rational; 3]) -> Id {
        self.graph.add(
            "DIRECTION",
            vec![
                text("''"),
                Field::List(
                    q.iter()
                        .map(|q| {
                            text(if q.is_zero() {
                                "0."
                            } else if q.sign() == Sign::Positive {
                                "1."
                            } else {
                                "-1."
                            })
                        })
                        .collect(),
                ),
            ],
        )
    }
}
/// Builds a separate approximate numerical STEP artifact for N02's finite planar
/// linear family. Source endpoint/trim coordinates must agree exactly in retained
/// rational values; a merely tolerated gap refuses rather than being healed.
/// Curved carriers, non-axis frames and unsupported numerical bodies refuse.
/// Both loops and all shared source edge/vertex identities are serialized.
/// The source tolerance is metadata, not a STEP global uncertainty or body bound.
/// The caller's allowance covers all work; no global policy is changed.
pub fn build(
    body: &NumericalBody,
    options: &NumericalStepOptions,
    dec: &mut Decider,
) -> Result<NumericalStep, Refusal> {
    guarded(|| build_inner(body, options, dec))
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
    let mut w = Writer {
        graph: Graph::default(),
        options,
        dec,
        error: Rational::zero(),
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
        axis(&dir, false)?;
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
        let direction = w.direction(&dir);
        let magnitude = w.scalar(&dir.iter().find(|v| !v.is_zero()).unwrap().abs())?;
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
        if axis(&normal, true)? == axis(&x, true)? {
            return Err(bad("unsupported", "plane axes must be orthogonal"));
        }
        let y = cross(&normal, &x);
        let point = w.point(&origin)?;
        let n = w.direction(&normal);
        let reference = w.direction(&x);
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
                for i in 0..2 {
                    let uv = &trim.uv[i];
                    let value = std::array::from_fn::<_, 3, _>(|k| {
                        &origin[k] + &(&x[k] * &uv[0]) + (&y[k] * &uv[1])
                    });
                    if value != d.vertices[endpoint[i]].position {
                        return Err(bad(
                            "gap",
                            "retained trim and plane endpoint differ; no export healing",
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
    let description = format!(
        "catallax numerical approximate;source={};source_tolerance_mm={};source_checks={};source_repairs={};complete_body_bound=unestablished;max_scalar_rounding_error_mm={}",
        d.name,
        body.options().tolerance_mm,
        body.evidence().checks().len(),
        body.evidence().repairs().len(),
        w.error
    );
    let mut file = w.graph.write(&[solid], |roots| {
        // Reuse AP242 product/unit framing, not its exact-model uncertainty claim.
        crate::export::header(&options.product_name, roots, "0.")
            .lines()
            .filter(|line| !line.starts_with("#13="))
            .collect::<Vec<_>>()
            .join("\n")
            .replace("'catallax export'", &string(&description))
            .replace("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#13))", "")
            + "\n"
    });
    file.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    let emitted_evidence = emitted_checks(body, options, w.dec)?;
    let artifact = NumericalStep {
        bytes: file.into_bytes(),
        report: NumericalStepReport {
            source_name: d.name.clone(),
            source_definition: d.clone(),
            source_evidence: body.evidence().clone(),
            emitted_evidence,
            scalar_rounding_limit_mm: options.max_scalar_rounding_error_mm.clone(),
            max_scalar_rounding_error_mm: w.error,
            counts: [count.faces, count.edges, count.vertices],
            checks: vec![
                "source.retained-endpoint-and-trim-consistency".into(),
                "file.decimal-model-numerical-admission".into(),
                "file.canonical-byte-identity".into(),
                "file.planar-seam-contract".into(),
                "file.approximate-category-no-global-bound".into(),
            ],
        },
    };
    artifact.check_written(artifact.bytes(), w.dec)?;
    Ok(artifact)
}
