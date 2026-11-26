//! Display comparison records for `tools/display-diff`, in verification builds.
//!
//! Compiled only with the `display-capture` feature. Each recorded call writes
//! one line to standard error, with tab-separated fields:
//!
//! `displaycapture`, [`FORMAT`], thread name, function, spline carriers,
//! options, work counters and the SHA-256 digest of the result.
//!
//! The digest covers an explicit encoding of what a caller can observe:
//! coordinates as binary64 bits, indices, exact rationals in lowest terms,
//! refusal codes and texts, and marks. It does not use `Debug` text, so a
//! renamed field or changed formatting leaves it unchanged. Raise [`FORMAT`]
//! whenever the encoding of any record changes.

use catallax_numbers::{
    Answer, BoundKind, BoundedAnswer, CertificateScope, EvidenceUnit, MarkLog, OperationResult,
    OperationWork, Rational, Refusal, current_operation_work, operation_budget_active,
};
use catallax_shapes::{SplineSupport, Surface};
use catallax_topology::{BodyId, Model};
use sha2::{Digest, Sha256};

use crate::{DisplayMesh, Mesh, MeshControl, MeshOptions, PreviewMesh, PreviewOptions, Triangle};

/// Version of the record line and of every result encoding.
pub const FORMAT: u32 = 1;

/// A running digest of one result.
pub struct Encoder(Sha256);

impl Encoder {
    /// Length-prefixed text.
    pub fn text(&mut self, text: &str) {
        self.count(text.len());
        self.0.update(text.as_bytes());
    }

    /// A count, index or tag.
    pub fn count(&mut self, value: usize) {
        self.0.update((value as u64).to_le_bytes());
    }

    /// A topological handle index.
    pub fn id(&mut self, value: u32) {
        self.0.update(u64::from(value).to_le_bytes());
    }

    /// A coordinate, as its binary64 bits.
    pub fn float(&mut self, value: f64) {
        self.0.update(value.to_bits().to_le_bytes());
    }

    /// An exact rational in lowest terms.
    pub fn rational(&mut self, value: &Rational) {
        self.text(&value.to_string());
    }

    /// An optional value, with its presence.
    pub fn optional<T: ?Sized>(&mut self, value: Option<&T>, encode: impl FnOnce(&mut Self, &T)) {
        match value {
            None => self.count(0),
            Some(value) => {
                self.count(1);
                encode(self, value);
            }
        }
    }

    /// Coordinates, with their count.
    pub fn positions(&mut self, positions: &[[f64; 3]]) {
        self.count(positions.len());
        for position in positions {
            for &value in position {
                self.float(value);
            }
        }
    }
}

/// A result with an explicit encoding for comparison records.
pub trait Observable {
    /// Writes everything a caller can observe of this value.
    fn encode(&self, out: &mut Encoder);
}

/// Options with an explicit text for the record key.
pub trait Describe {
    /// One line naming every option and its value.
    fn describe(&self) -> String;
}

fn optional_rational(value: Option<&Rational>) -> String {
    value.map_or_else(|| "-".to_owned(), Rational::to_string)
}

impl Describe for MeshOptions {
    fn describe(&self) -> String {
        format!(
            "chord_tolerance={} angular_deviation={} max_edge_length={} max_aspect_ratio={} max_vertices={}",
            self.chord_tolerance,
            optional_rational(self.angular_deviation.as_ref()),
            optional_rational(self.max_edge_length.as_ref()),
            optional_rational(self.max_aspect_ratio.as_ref()),
            self.max_vertices,
        )
    }
}

impl Describe for PreviewOptions {
    fn describe(&self) -> String {
        format!(
            "chord_tolerance={} angular_deviation={} max_vertices={}",
            self.chord_tolerance,
            optional_rational(self.angular_deviation.as_ref()),
            self.max_vertices,
        )
    }
}

impl Observable for Refusal {
    fn encode(&self, out: &mut Encoder) {
        out.text(self.code);
        out.text(&self.detail);
        out.count(usize::from(self.undecided().is_some()));
    }
}

impl Observable for MarkLog {
    fn encode(&self, out: &mut Encoder) {
        out.count(self.len());
        for mark in self.iter() {
            out.text(mark.code);
            out.text(&mark.detail);
            out.rational(&mark.bound);
            out.text(match mark.kind {
                BoundKind::Distance => "distance",
                BoundKind::Angle => "angle",
                BoundKind::Value => "value",
            });
        }
    }
}

impl<T: Observable> Observable for Result<T, Refusal> {
    fn encode(&self, out: &mut Encoder) {
        match self {
            Ok(value) => {
                out.text("ok");
                value.encode(out);
            }
            Err(refusal) => {
                out.text("refused");
                refusal.encode(out);
            }
        }
    }
}

impl<T: Observable> Observable for Answer<T> {
    fn encode(&self, out: &mut Encoder) {
        self.value.encode(out);
        self.marks.encode(out);
    }
}

fn unit(unit: EvidenceUnit) -> &'static str {
    match unit {
        EvidenceUnit::Millimetres => "millimetres",
        EvidenceUnit::Radians => "radians",
        EvidenceUnit::Scalar => "scalar",
        EvidenceUnit::SquareMillimetres => "square-millimetres",
        EvidenceUnit::CubicMillimetres => "cubic-millimetres",
        EvidenceUnit::Raw => "raw",
    }
}

impl<T: Observable> Observable for BoundedAnswer<T> {
    fn encode(&self, out: &mut Encoder) {
        self.value().encode(out);
        let evidence = self.evidence();
        let rule = evidence.rule();
        out.text(rule.case_kind);
        out.text(rule.name);
        out.text(rule.level);
        out.text(unit(rule.unit));
        out.rational(evidence.level());
        out.rational(evidence.measured());
        let definition = evidence.definition();
        out.text(&definition.revision);
        out.text(&definition.representation);
        out.text(&definition.target);
        match &definition.scope {
            CertificateScope::WholeSolid => out.text("whole-solid"),
            CertificateScope::Entities(entities) => {
                out.text("entities");
                out.count(entities.len());
                for entity in entities {
                    out.text(entity);
                }
            }
            CertificateScope::ParameterBox { chart, bounds } => {
                out.text("parameter-box");
                out.text(chart);
                out.count(bounds.len());
                for (lower, upper) in bounds {
                    out.rational(lower);
                    out.rational(upper);
                }
            }
            CertificateScope::Query(query) => {
                out.text("query");
                out.text(query);
            }
        }
    }
}

impl<T: Observable> Observable for OperationResult<T> {
    fn encode(&self, out: &mut Encoder) {
        match self {
            Self::Proven(answer) => {
                out.text("proven");
                answer.encode(out);
            }
            Self::Bounded(answer) => {
                out.text("bounded");
                answer.encode(out);
            }
            Self::Refused(refusal) => {
                out.text("refused");
                refusal.encode(out);
            }
        }
    }
}

fn triangles(out: &mut Encoder, triangles: &[Triangle]) {
    out.count(triangles.len());
    for triangle in triangles {
        for vertex in triangle.vertices {
            out.count(vertex);
        }
        out.id(triangle.face.index());
    }
}

impl Observable for Mesh {
    fn encode(&self, out: &mut Encoder) {
        out.positions(&self.positions);
        triangles(out, &self.triangles);
        out.count(self.faces.len());
        for face in &self.faces {
            out.id(face.face.index());
            out.count(face.triangles.start);
            out.count(face.triangles.end);
            out.rational(&face.max_chord_error);
            out.optional(face.max_aspect_ratio.as_ref(), Encoder::rational);
            out.optional(face.max_angular_deviation.as_ref(), Encoder::rational);
        }
        out.count(self.edges.len());
        for edge in &self.edges {
            out.id(edge.edge.index());
            out.count(edge.vertices.len());
            for &vertex in &edge.vertices {
                out.count(vertex);
            }
            out.rational(&edge.max_chord_error);
        }
        out.count(self.unavailable_controls.len());
        for control in &self.unavailable_controls {
            out.id(control.face.index());
            out.text(match control.control {
                MeshControl::AngularDeviation => "angular-deviation",
                MeshControl::EdgeLength => "edge-length",
                MeshControl::AspectRatio => "aspect-ratio",
            });
            out.text(control.reason);
        }
    }
}

impl Observable for PreviewMesh {
    fn encode(&self, out: &mut Encoder) {
        out.positions(&self.positions);
        triangles(out, &self.triangles);
        out.count(self.faces.len());
        for face in &self.faces {
            out.id(face.face.index());
            out.count(face.triangles.start);
            out.count(face.triangles.end);
        }
        out.count(self.edges.len());
        for edge in &self.edges {
            out.id(edge.edge.index());
            out.count(edge.vertices.len());
            for &vertex in &edge.vertices {
                out.count(vertex);
            }
        }
    }
}

impl Observable for DisplayMesh {
    fn encode(&self, out: &mut Encoder) {
        self.geometry.encode(out);
        out.text(self.measurement.method);
        out.rational(&self.measurement.maximum_deviation);
        out.count(self.measurement.covered_faces);
        self.source_marks.encode(out);
    }
}

/// Encodes a field display mesh. Each triangle gives its vertices and the
/// index of its field registry face.
pub fn field_mesh(
    out: &mut Encoder,
    positions: &[[f64; 3]],
    triangles: impl ExactSizeIterator<Item = ([usize; 3], usize)>,
    max_deviation: &Rational,
    cells: usize,
) {
    out.positions(positions);
    out.count(triangles.len());
    for (vertices, face) in triangles {
        for vertex in vertices {
            out.count(vertex);
        }
        out.count(face);
    }
    out.rational(max_deviation);
    out.count(cells);
}

/// Spline carrier kinds of a body's faces, sorted and comma-separated.
pub fn carriers(model: &Model, body: BodyId) -> String {
    let kinds = model
        .faces_of_body(body)
        .into_iter()
        .filter_map(
            |face| match model.geom().surface(model.face(face).surface) {
                Surface::BSpline(support) => Some(match support {
                    SplineSupport::BezierExtrusion(_) => "bezier",
                    SplineSupport::EllipticOffset(_) => "ellipse-offset",
                    SplineSupport::LinearFillet(_) => "linear-fillet",
                    SplineSupport::ObliqueDraft(_) => "oblique-draft",
                    SplineSupport::RationalAtlas(_) => "atlas",
                    SplineSupport::RationalPatch(_) => "rational-patch",
                    SplineSupport::RollingEnvelope(_) => "envelope",
                    SplineSupport::RollingMitre(_) => "mitre",
                    SplineSupport::RollingPipe(_) => "pipe",
                    SplineSupport::RollingSaddle(_) => "saddle",
                    SplineSupport::SplineFillet(_) => "spline-fillet",
                    _ => "other-spline",
                }),
                _ => None,
            },
        )
        .collect::<std::collections::BTreeSet<_>>();
    kinds.into_iter().collect::<Vec<_>>().join(",")
}

/// One display call being recorded.
pub struct Call {
    function: &'static str,
    carriers: String,
    before: Option<OperationWork>,
}

impl Call {
    /// Starts recording before the call does any work.
    pub fn start(function: &'static str, carriers: String) -> Self {
        Self {
            function,
            carriers,
            before: operation_budget_active().then(current_operation_work),
        }
    }

    /// Writes the record line for the finished call.
    pub fn finish(self, options: &impl Describe, result: &impl Observable) {
        let work = self.before.map_or_else(
            || "null".to_owned(),
            |before| {
                let after = current_operation_work();
                format!(
                    "[{}, {}, {}, {}]",
                    after.used.wrapping_sub(before.used),
                    after.refunded.wrapping_sub(before.refunded),
                    after.uncharged.wrapping_sub(before.uncharged),
                    u64::from(after.exhausted),
                )
            },
        );
        let mut out = Encoder(Sha256::new());
        result.encode(&mut out);
        let digest: String = out
            .0
            .finalize()
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect();
        let thread = std::thread::current();
        eprintln!(
            "\ndisplaycapture\t{FORMAT}\t{}\t{}\t{}\t{}\t{work}\t{digest}",
            thread.name().unwrap_or("unnamed"),
            self.function,
            self.carriers,
            options.describe(),
        );
    }
}
