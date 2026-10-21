//! Catallax CLI and the legacy `exact-brep` exchange program.
//!
//! Reads question files, answers them with pieces 1 to 7, and writes the
//! manifest and STEP file defined in `docs/JUDGE-INTERFACE.md`. This crate is
//! not a piece; it holds no geometry of its own.

mod manifest;
mod question;

pub use manifest::{
    Bound, KernelInfo, Manifest, MarkEntry, Properties, RefusalEntry, ResultSummary, Stats, Status,
    StepStats, StepSummary, render_manifest,
};
pub use question::{
    ExportRequest, Limits, MeasureRequest, Operand, Question, QuestionError, ShapeSpec, Step,
    parse_question,
};

use std::{
    collections::BTreeMap,
    fs::{self, File},
    io,
    path::Path,
    time::{Duration, Instant},
};

use catallax_numbers::{BoundKind, Decider, Mark, Refusal};
use catallax_topology::{BodyId, BodyKind, Model};

/// The question format this build reads.
pub const QUESTION_FORMAT: &str = "exact-brep-question/2";

/// The manifest format this build writes.
pub const MANIFEST_FORMAT: &str = "exact-brep-manifest/2";

/// Answers a question and writes its STEP file when requested.
/// The caller writes the returned manifest with [`render_manifest`].
/// Output failures are I/O errors; geometric refusals are manifests.
pub fn answer(question: &Question, out_dir: &Path) -> io::Result<Manifest> {
    answer_since(question, out_dir, Instant::now())
}

/// Answers using the executable's start time for the deadline and wall time.
pub fn answer_since(question: &Question, out_dir: &Path, started: Instant) -> io::Result<Manifest> {
    prepare_output(out_dir)?;
    let mut manifest = failure(
        Some(question.id.clone()),
        Status::Refused,
        "cli.internal-error",
        "answer did not complete",
    );
    let mut model = Model::new();
    let mut dec = Decider::default();
    dec.set_deadline(started.checked_add(Duration::from_millis(question.limits.deadline_ms)));
    let result = evaluate(question, out_dir, &mut model, &mut dec, &mut manifest)
        .and_then(|()| dec.check_deadline().map_err(Failure::from));
    manifest.stats.wall_ms = milliseconds(started.elapsed());
    manifest.stats.decisions = dec.take_stats();
    match result {
        Ok(()) => {}
        Err(Failure::Refusal(refusal, step, stage)) => {
            if refusal.code == "numbers.deadline" {
                prepare_output(out_dir)?;
            }
            manifest.status = Status::Refused;
            manifest.exact = false;
            manifest.result = None;
            manifest.properties = None;
            manifest.step = None;
            manifest.refusal = Some(RefusalEntry {
                code: refusal.code.into(),
                step,
                stage,
                detail: refusal.detail,
            });
        }
        Err(Failure::Io(error)) => return Err(error),
    }
    Ok(manifest)
}

/// Creates the output directory and removes any stale STEP answer.
pub fn prepare_output(out_dir: &Path) -> io::Result<()> {
    fs::create_dir_all(out_dir)?;
    match fs::remove_file(out_dir.join("answer.step")) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error),
    }
}

/// Makes an invalid-question or internal-error manifest.
pub fn failure(
    question_id: Option<String>,
    status: Status,
    code: &str,
    detail: impl Into<String>,
) -> Manifest {
    Manifest {
        question_id,
        kernel: kernel_info(),
        status,
        exact: false,
        result: None,
        properties: None,
        step: None,
        marks: Vec::new(),
        refusal: Some(RefusalEntry {
            code: code.into(),
            step: None,
            stage: None,
            detail: detail.into(),
        }),
        stats: Stats::default(),
    }
}

fn kernel_info() -> KernelInfo {
    KernelInfo {
        name: "exact-brep".into(),
        version: env!("CARGO_PKG_VERSION").into(),
        commit: env!("CATALLAX_COMMIT").into(),
    }
}

enum Failure {
    Refusal(Refusal, Option<String>, Option<String>),
    Io(io::Error),
}

impl From<Refusal> for Failure {
    fn from(value: Refusal) -> Self {
        Self::Refusal(value, None, None)
    }
}

impl From<io::Error> for Failure {
    fn from(value: io::Error) -> Self {
        Self::Io(value)
    }
}

struct Built {
    body: BodyId,
    marks: Vec<MarkEntry>,
}

fn evaluate(
    question: &Question,
    out_dir: &Path,
    model: &mut Model,
    dec: &mut Decider,
    manifest: &mut Manifest,
) -> Result<(), Failure> {
    let mut bodies = BTreeMap::new();
    for (name, operand) in &question.operands {
        dec.check_deadline()?;
        let pose = &operand.pose;
        let body = match &operand.shape {
            ShapeSpec::Block { size } => catallax_topology::block(model, size.clone(), pose),
            ShapeSpec::Cylinder { radius, height } => {
                catallax_topology::cylinder(model, radius.clone(), height.clone(), pose)
            }
            ShapeSpec::Sphere { radius } => catallax_topology::sphere(model, radius.clone(), pose),
            ShapeSpec::Cone {
                base_radius,
                top_radius,
                height,
            } => catallax_topology::cone(
                model,
                base_radius.clone(),
                top_radius.clone(),
                height.clone(),
                pose,
            ),
            ShapeSpec::Torus { major, minor } => {
                catallax_topology::torus(model, major.clone(), minor.clone(), pose)
            }
        }?;
        let marks = model
            .body(body)
            .marks
            .iter()
            .map(|mark| mark_entry(mark, None))
            .collect();
        bodies.insert(name.as_str(), Built { body, marks });
    }
    for step in &question.steps {
        dec.check_deadline().map_err(|refusal| {
            Failure::Refusal(refusal, Some(step.id.clone()), Some("intersect".into()))
        })?;
        let a = &bodies[step.args[0].as_str()];
        let b = &bodies[step.args[1].as_str()];
        if a.body == b.body {
            let body = match step.op {
                catallax_boolean::BoolOp::Unite | catallax_boolean::BoolOp::Intersect => a.body,
                catallax_boolean::BoolOp::Subtract => model.add_body(catallax_topology::Body {
                    kind: BodyKind::Empty,
                    lumps: Vec::new(),
                    marks: model.body(a.body).marks.clone(),
                }),
                _ => panic!("unknown boolean operation"),
            };
            let marks = a.marks.clone();
            manifest.stats.per_step.push(StepStats {
                id: step.id.clone(),
                wall_ms: 0,
                decisions: 0,
                not_exact: 0,
            });
            bodies.insert(step.id.as_str(), Built { body, marks });
            continue;
        }
        let before = dec.stats().clone();
        let start = Instant::now();
        let result = catallax_boolean::boolean(model, a.body, b.body, step.op, dec);
        manifest.stats.per_step.push(StepStats {
            id: step.id.clone(),
            wall_ms: milliseconds(start.elapsed()),
            decisions: dec.stats().decisions - before.decisions,
            not_exact: dec.stats().not_exact - before.not_exact,
        });
        let answer = result.map_err(|error| {
            Failure::Refusal(
                error.refusal,
                Some(step.id.clone()),
                Some(stage_name(error.stage).into()),
            )
        })?;
        let added = answer
            .marks
            .len()
            .checked_sub(a.marks.len() + b.marks.len())
            .expect("boolean retains operand marks");
        let mut marks: Vec<_> = answer
            .marks
            .iter()
            .take(added)
            .map(|mark| mark_entry(mark, Some(&step.id)))
            .collect();
        marks.extend_from_slice(&a.marks);
        marks.extend_from_slice(&b.marks);
        bodies.insert(
            step.id.as_str(),
            Built {
                body: answer.body,
                marks,
            },
        );
    }
    let built = &bodies[question.result.as_str()];
    manifest.marks.clone_from(&built.marks);
    let request = &question.measure;
    let mut properties = Properties {
        volume: None,
        area: None,
        centroid: None,
    };
    if request.centroid {
        let answer = catallax_measure::mass_properties(model, built.body, &request.accuracy, dec)?;
        manifest
            .marks
            .extend(answer.marks.iter().map(|mark| mark_entry(mark, None)));
        if request.volume {
            properties.volume = Some(bound(&answer.value.volume, &request.accuracy));
        }
        if request.area {
            properties.area = Some(bound(&answer.value.area, &request.accuracy));
        }
        properties.centroid = Some(answer.value.centroid.as_ref().map(|values| {
            values
                .each_ref()
                .map(|value| bound(value, &request.accuracy))
        }));
    } else {
        if request.volume {
            let answer = catallax_measure::volume(model, built.body, &request.accuracy, dec)?;
            manifest
                .marks
                .extend(answer.marks.iter().map(|mark| mark_entry(mark, None)));
            properties.volume = Some(bound(&answer.value, &request.accuracy));
        }
        if request.area {
            let answer = catallax_measure::area(model, built.body, &request.accuracy, dec)?;
            manifest
                .marks
                .extend(answer.marks.iter().map(|mark| mark_entry(mark, None)));
            properties.area = Some(bound(&answer.value, &request.accuracy));
        }
    }
    mark_decimal_accuracy(&properties, &request.accuracy, &mut manifest.marks);
    let summary = result_summary(model, built.body, dec)?;
    if question.export.step && model.body(built.body).kind != BodyKind::Empty {
        let path = out_dir.join("answer.step");
        let mut file = File::create(&path)?;
        let exported = catallax_step::write_step(
            model,
            &[built.body],
            &question.export.options,
            dec,
            &mut file,
        );
        match exported {
            Ok(answer) => {
                let report = answer.value;
                manifest
                    .marks
                    .extend(answer.marks.iter().map(|mark| mark_entry(mark, None)));
                let max_deviation = report.max_deviation.to_decimal_bounds(40).1;
                let decimal_deviation = catallax_numbers::Rational::parse(&max_deviation)
                    .expect("decimal deviation parses");
                if report.max_deviation <= question.export.options.curve_tolerance
                    && decimal_deviation > question.export.options.curve_tolerance
                {
                    manifest.marks.push(mark_entry(
                        &Mark::new(
                            "step.tolerance-not-reached",
                            BoundKind::Distance,
                            decimal_deviation,
                            "outward-rounded deviation exceeds the requested tolerance",
                        ),
                        None,
                    ));
                }
                manifest.step = Some(StepSummary {
                    file: "answer.step".into(),
                    max_deviation,
                    approximated_curves: report.approximated_curves,
                    added_vertices: report.added_vertices,
                    added_edges: report.added_edges,
                });
            }
            Err(error) => {
                drop(file);
                fs::remove_file(&path)?;
                if error.code == "step.io-error" {
                    return Err(Failure::Io(io::Error::other(error.detail)));
                }
                return Err(error.into());
            }
        }
    }
    manifest.status = Status::Answered;
    manifest.exact = manifest.marks.is_empty();
    manifest.result = Some(summary);
    manifest.properties = Some(properties);
    manifest.refusal = None;
    Ok(())
}

fn result_summary(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<ResultSummary, Refusal> {
    let counts = model.counts(body);
    let mut faces_by_kind = BTreeMap::new();
    for (kind, count) in counts.faces_by_kind {
        let name = match kind {
            catallax_shapes::SurfaceKind::Plane => "plane",
            catallax_shapes::SurfaceKind::Cylinder => "cylinder",
            catallax_shapes::SurfaceKind::Cone => "cone",
            catallax_shapes::SurfaceKind::Sphere => "sphere",
            catallax_shapes::SurfaceKind::Torus => "torus",
            _ => panic!("unknown surface kind"),
        };
        faces_by_kind.insert(name.into(), count);
    }
    let mut edges_by_kind = BTreeMap::new();
    for edge in model.edges_of_body(body) {
        let name = match catallax_intersect::curve_form(model.geom(), model.edge(edge).curve, dec)?
        {
            catallax_intersect::CurveForm::Line { .. } => "line",
            catallax_intersect::CurveForm::Circle { .. } => "circle",
            catallax_intersect::CurveForm::Ellipse { .. } => "ellipse",
            catallax_intersect::CurveForm::Branch => "intersection",
            _ => panic!("unknown curve form"),
        };
        *edges_by_kind.entry(name.into()).or_insert(0) += 1;
    }
    Ok(ResultSummary {
        kind: if model.body(body).kind == BodyKind::Empty {
            "empty"
        } else {
            "solid"
        }
        .into(),
        lumps: counts.lumps,
        shells: counts.shells,
        faces: counts.faces,
        edges: counts.edges,
        ring_edges: counts.ring_edges,
        vertices: counts.vertices,
        faces_by_kind,
        edges_by_kind,
    })
}

fn bound(value: &catallax_measure::Enclosure, accuracy: &catallax_measure::Accuracy) -> Bound {
    let magnitude = if value.lo.sign() == value.hi.sign() {
        value.lo.abs().min(value.hi.abs())
    } else {
        catallax_numbers::Rational::zero()
    };
    let allowed = (&accuracy.relative * &magnitude).max(accuracy.absolute.clone());
    let mut digits = 40;
    loop {
        let lo = value.lo.to_decimal_bounds(digits).0;
        let hi = value.hi.to_decimal_bounds(digits).1;
        let rounded_lo = catallax_numbers::Rational::parse(&lo).expect("decimal bound parses");
        let rounded_hi = catallax_numbers::Rational::parse(&hi).expect("decimal bound parses");
        if &rounded_hi - &rounded_lo <= allowed || value.width() >= allowed {
            return Bound { lo, hi };
        }
        digits = digits.checked_mul(2).expect("decimal precision fits u32");
    }
}

fn mark_decimal_accuracy(
    properties: &Properties,
    accuracy: &catallax_measure::Accuracy,
    marks: &mut Vec<MarkEntry>,
) {
    let mut widest = catallax_numbers::Rational::zero();
    let centroid = properties.centroid.as_ref().and_then(Option::as_ref);
    for bound in properties
        .volume
        .iter()
        .chain(&properties.area)
        .chain(centroid.into_iter().flatten())
    {
        let lo = catallax_numbers::Rational::parse(&bound.lo).expect("decimal lower bound parses");
        let hi = catallax_numbers::Rational::parse(&bound.hi).expect("decimal upper bound parses");
        let magnitude = if lo.sign() == hi.sign() {
            lo.abs().min(hi.abs())
        } else {
            catallax_numbers::Rational::zero()
        };
        let width = hi - lo;
        if width > (&accuracy.relative * &magnitude).max(accuracy.absolute.clone()) {
            widest = widest.max(width);
        }
    }
    if !widest.is_zero()
        && !marks.iter().any(|mark| {
            mark.code == "measure.accuracy-not-reached"
                && catallax_numbers::Rational::parse(&mark.bound).expect("mark bound parses")
                    >= widest
        })
    {
        marks.push(mark_entry(
            &Mark::new(
                "measure.accuracy-not-reached",
                BoundKind::Value,
                widest,
                "decimal property enclosure exceeds the requested accuracy",
            ),
            None,
        ));
    }
}

fn mark_entry(mark: &Mark, step: Option<&str>) -> MarkEntry {
    MarkEntry {
        code: mark.code.into(),
        kind: match mark.kind {
            BoundKind::Distance => "distance",
            BoundKind::Angle => "angle",
            BoundKind::Value => "value",
        }
        .into(),
        bound: mark.bound.to_decimal_bounds(40).1,
        step: step.map(str::to_owned),
        detail: mark.detail.clone(),
    }
}

fn stage_name(stage: catallax_boolean::Stage) -> &'static str {
    match stage {
        catallax_boolean::Stage::Intersect => "intersect",
        catallax_boolean::Stage::Imprint => "imprint",
        catallax_boolean::Stage::Classify => "classify",
        catallax_boolean::Stage::Keep => "keep",
        catallax_boolean::Stage::Sew => "sew",
    }
}

fn milliseconds(duration: Duration) -> u64 {
    duration.as_millis().try_into().unwrap_or(u64::MAX)
}

/// The JSON object printed by `exact-brep about`.
pub fn about() -> String {
    #[derive(serde::Serialize)]
    struct About {
        #[serde(flatten)]
        kernel: KernelInfo,
        question_formats: [&'static str; 1],
        manifest_format: &'static str,
        shapes: [&'static str; 5],
        operations: [&'static str; 3],
        properties: [&'static str; 3],
    }
    serde_json::to_string(&About {
        kernel: kernel_info(),
        question_formats: [QUESTION_FORMAT],
        manifest_format: MANIFEST_FORMAT,
        shapes: ["block", "cylinder", "sphere", "cone", "torus"],
        operations: ["unite", "subtract", "intersect"],
        properties: ["volume", "area", "centroid"],
    })
    .expect("about fields serialize as JSON")
}
