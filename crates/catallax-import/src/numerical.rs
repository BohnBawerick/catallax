//! A separate numerical result for source-incidence repair. Exact bodies never
//! receive a numerical cast, and numerical bodies expose no exact model.

use std::collections::BTreeMap;

use catallax_numbers::{
    Decider, Rational, Real, Refusal, Sign, operation_budget_active, with_operation_budget,
};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Line, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, Model, Sense, numerical as n};

use crate::{ImportRefusal, RepairRecord, Solid, Vec3};

pub(crate) struct CheckedStageProbe<'a> {
    stage: &'static str,
    body: &'a str,
}

impl<'a> CheckedStageProbe<'a> {
    pub(crate) fn new(stage: &'static str, body: &'a str) -> Self {
        let work = catallax_numbers::current_operation_work();
        eprintln!(
            "CX_CHECKED_STAGE\tbegin\t{stage}\t{body}\t{}\t{}",
            work.used, work.exhausted
        );
        Self { stage, body }
    }
}

impl Drop for CheckedStageProbe<'_> {
    fn drop(&mut self) {
        let work = catallax_numbers::current_operation_work();
        eprintln!(
            "CX_CHECKED_STAGE\tend\t{}\t{}\t{}\t{}",
            self.stage, self.body, work.used, work.exhausted
        );
    }
}

mod analytic;
mod cones;

pub(crate) use analytic::constructed::sew_with_reconstruction;
pub use analytic::constructed::{SewnSurfaceBody, sew_constructed_surfaces};

/// Numerical reconstruction of finite parameter trims, retaining the immutable
/// source's shared vertices and analytic carriers. No geometric displacement or
/// complete source/body deviation bound is inferred from its sampled checks.
#[derive(Clone, Debug)]
pub struct NumericalSourceAdmission {
    /// The strict attempt's actual refusal, before numerical admission.
    pub strict_refusal: ImportRefusal,
    /// Fixed maximum residual at each completed numerical check, in millimetres.
    pub tolerance_mm: Rational,
    /// Versioned source interpretation and repair method.
    pub method: String,
    /// Explicit changes to source topology; source geometry is retained.
    pub topology_repairs: Vec<NumericalTopologyRepair>,
    /// Every within-solid source face pair, in source order. A tolerated
    /// incidence prerequisite is unresolved only for pairs using that face.
    /// Separate solids are not compared, matching the existing import policy.
    pub intersection_checks: Vec<NumericalSourcePairCheck>,
    /// One collective prepayment for optional curved source-pair proofs.
    pub optional_exact_work_reserved: u64,
    /// Actual work used inside that collective proof allowance.
    pub optional_exact_work_used: u64,
}

/// A derived cylinder seam joins two retained noncontractible source rims.
/// Its two coedges share one new 3D line and differ by one UV period.
#[derive(Clone, Debug, serde::Serialize)]
pub struct NumericalTopologyRepair {
    /// Source solid name.
    pub solid: String,
    /// Source cylindrical face name.
    pub face: String,
    /// Original complete circular rim edge names, in joined boundary order.
    pub source_rim_edges: [String; 2],
    /// Original rim anchor vertex names; no new vertex positions are introduced.
    pub source_vertices: [String; 2],
    /// Derived edge key, distinct from source entity names.
    pub derived_seam_edge: String,
    /// The two former one-edge loops become one bounded parameter cell.
    pub method: String,
}

/// Source intersection coverage retained separately from numerical samples.
#[derive(Clone, Debug, serde::Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum NumericalSourcePairCheck {
    /// The source checker passed this pair under its existing permitted contact
    /// rules. This does not assert separated boundaries.
    CheckedNoCrossing {
        /// Stable source solid name.
        solid: String,
        /// Stable source face names, in source order.
        faces: [String; 2],
        /// Work used by this pair check; shared preparation is recorded separately.
        exact_attempt_work_used: u64,
    },
    /// The source pair could not be checked exactly. This grants no separation.
    Unresolved {
        /// Stable source solid name.
        solid: String,
        /// Stable source face names, in source order.
        faces: [String; 2],
        /// The actual typed incidence prerequisite or undecided pair refusal.
        refusal: Box<ImportRefusal>,
        /// Work used by this pair check, or zero when no attempt completed.
        exact_attempt_work_used: u64,
    },
}

/// A source refusal or a numerical admission refusal. Neither contains partial geometry.
#[derive(Clone, Debug)]
pub enum ImportError {
    /// Recognition, source limits, repair or certified construction refused.
    Source(ImportRefusal),
    /// Numerical definition, tolerance or checked admission refused.
    Numerical(Refusal),
}

impl From<ImportRefusal> for ImportError {
    fn from(value: ImportRefusal) -> Self {
        Self::Source(value)
    }
}
impl From<Refusal> for ImportError {
    fn from(value: Refusal) -> Self {
        Self::Numerical(value)
    }
}

/// An import admitted by the existing certified path.
#[derive(Clone, Debug)]
pub struct CertifiedImport {
    model: Model,
    bodies: Vec<BodyId>,
    source_sha256: String,
}
impl CertifiedImport {
    /// The existing certified model; body marks retain their existing meaning.
    pub fn model(&self) -> &Model {
        &self.model
    }
    /// Source bodies, in source order.
    pub fn bodies(&self) -> &[BodyId] {
        &self.bodies
    }
    /// Actual immutable source digest.
    pub fn source_sha256(&self) -> &str {
        &self.source_sha256
    }
}

/// Numerical source bodies and the actual source-repair record.
/// Source repair history is separate from the numerical builder's own repairs.
#[derive(Clone, Debug)]
pub struct NumericalImport {
    bodies: Vec<(String, n::NumericalBody)>,
    record: RepairRecord,
    source_admission: Option<NumericalSourceAdmission>,
}
impl NumericalImport {
    /// Numerical bodies keyed by source solid. Zero residuals remain numerical.
    pub fn bodies(&self) -> &[(String, n::NumericalBody)] {
        &self.bodies
    }
    /// Complete source interpretation, source limits and source repair history.
    pub fn record(&self) -> &RepairRecord {
        &self.record
    }
    /// Actual numerical source interpretation, when finite analytic trims were
    /// rebuilt without an exact planar derivative. Body evidence lists every
    /// completed location and check; this record grants no exact incidence.
    pub fn source_admission(&self) -> Option<&NumericalSourceAdmission> {
        self.source_admission.as_ref()
    }
}

/// The result category is explicit and must be handled by the consumer.
#[derive(Clone, Debug)]
pub enum StepImport {
    /// An input already handled by the certified import path.
    Certified(Box<CertifiedImport>),
    /// A separately checked numerical B-rep, with no exact-model extraction.
    Numerical(Box<NumericalImport>),
}

fn unsupported(detail: &str) -> ImportError {
    Refusal::new("import.numerical-unsupported-family", detail).into()
}

/// Try the existing certified path, then numerically admit the bounded planar
/// source-incidence family. The raw source stays immutable.
///
/// The numerical path retains the recorded producer's actual source bounds and
/// repairs. It constructs a separate N01/N02 definition and rechecks it at the
/// caller's fixed tolerance. Every returned fallback stays numerical, including
/// zero sampled residuals and derivatives with cheap exact construction proofs.
/// No result is inserted into an exact caller model. The active request allowance
/// applies; without one, the whole query shares the unchanged 20M default.
pub fn import_step(
    bytes: &[u8],
    expected_sha256: &str,
    step: &StepId,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<StepImport, ImportError> {
    if operation_budget_active() {
        return import_inner(bytes, expected_sha256, step, options, dec);
    }
    let (result, work) = with_operation_budget(n::DEFAULT_WORK, || {
        import_inner(bytes, expected_sha256, step, options, dec)
    });
    if work.exhausted {
        return Err(ImportRefusal::ResourceLimit.into());
    }
    result
}

fn import_inner(
    bytes: &[u8],
    expected_sha256: &str,
    step: &StepId,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<StepImport, ImportError> {
    let _probe = CheckedStageProbe::new("numerical.file_import", "file");
    // Verification is part of this request, including a mismatch that would
    // otherwise return before any geometry work. Do not hash an unmetered slice.
    dec.charge_work(bytes.len() as u64).map_err(|refusal| {
        if catallax_numbers::current_operation_work().exhausted {
            ImportError::Source(ImportRefusal::ResourceLimit)
        } else {
            ImportError::Numerical(refusal)
        }
    })?;
    let mut model = Model::new();
    let strict_refusal = match crate::import_step_bodies(&mut model, bytes, expected_sha256, step) {
        Ok(bodies) => {
            return Ok(StepImport::Certified(Box::new(CertifiedImport {
                model,
                bodies,
                source_sha256: expected_sha256.to_owned(),
            })));
        }
        Err(error) if cones::eligible(bytes, &error, dec)? => {
            return analytic::import(bytes, expected_sha256, error, options, dec);
        }
        Err(error @ ImportRefusal::NotSolid { .. }) if matches!(&error, ImportRefusal::NotSolid { detail } if detail.contains("exactly")) => {
            error
        }
        Err(error) if analytic::eligible(&error) => {
            return analytic::import(bytes, expected_sha256, error, options, dec);
        }
        Err(error) => return Err(error.into()),
    };
    if options.tolerance_mm.sign() != Sign::Positive
        || options.tolerance_mm > crate::gap::auto_cap()
    {
        return Err(Refusal::new(
            "import.numerical-tolerance",
            "source tolerance must be positive and at most 0.01 mm",
        )
        .into());
    }
    let file = crate::read_step(bytes)?;
    if cones::in_file(&file, None, dec)? {
        return analytic::import(bytes, expected_sha256, strict_refusal, options, dec);
    }
    let imported = match crate::repair::apply_file(
        &mut Model::new(),
        file,
        expected_sha256,
        "source.step",
        &Default::default(),
        &[],
        Some(step),
        None,
    ) {
        Ok(imported) => imported,
        Err(error) if cones::eligible(bytes, &error, dec)? => {
            return analytic::import(bytes, expected_sha256, error, options, dec);
        }
        Err(error) if analytic::eligible(&error) => {
            return analytic::import(bytes, expected_sha256, error, options, dec);
        }
        Err(error) => return Err(error.into()),
    };
    if !imported.left_open.is_empty() {
        return Err(unsupported("the source repair left an open solid"));
    }
    for bound in &imported.record.face_tolerances {
        let upper = Rational::parse(&bound.deviation_upper_mm).expect("producer bound");
        if upper > options.tolerance_mm {
            return Err(Refusal::new(
                "import.numerical-source-limit",
                format!(
                    "source bound {} mm exceeds the requested {} mm",
                    bound.deviation_upper_mm, options.tolerance_mm
                ),
            )
            .into());
        }
    }
    let definition = imported.comparison.complete_definition().ok_or_else(|| {
        unsupported("the builder did not retain the complete repaired definition")
    })?;
    if !definition.surface_models.is_empty() {
        return Err(unsupported(
            "numerical continuation must retain every surface model and its surface body kind",
        ));
    }
    let mut bodies = Vec::new();
    for solid in &definition.solids {
        let definition = planar_definition(solid, dec)?;
        bodies.push((
            solid.name.clone(),
            n::build(definition, options.clone(), dec)?,
        ));
    }
    Ok(StepImport::Numerical(Box::new(NumericalImport {
        bodies,
        record: imported.record,
        source_admission: None,
    })))
}

// Only verified, retained source definitions from the opened-file transaction
// enter here. The existing byte entry and its planar-only policy are unchanged.
pub(crate) fn admit_retained_definition(
    file: &crate::StepFile,
    verified_sha256: &str,
    strict_refusal: ImportRefusal,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<StepImport, ImportError> {
    let _probe = CheckedStageProbe::new("numerical.retained_definition", "file");
    if options.tolerance_mm.sign() != Sign::Positive
        || options.tolerance_mm > crate::gap::auto_cap()
    {
        return Err(unsupported(
            "retained numerical tolerance is outside the import cap",
        ));
    }
    analytic::admit_definition(file, verified_sha256, strict_refusal, options, true, dec)
}

fn point(value: &Vec3) -> catallax_shapes::Point3 {
    RPoint3::new(value[0].clone(), value[1].clone(), value[2].clone()).to_real()
}
fn vector(value: &Vec3) -> catallax_shapes::Vector3 {
    RVector3::new(value[0].clone(), value[1].clone(), value[2].clone()).to_real()
}

fn vertex(
    vertices: &mut Vec<n::NumericalVertex>,
    indices: &mut BTreeMap<Vec3, usize>,
    source: &crate::Vertex,
    dec: &mut Decider,
) -> Result<usize, ImportError> {
    dec.charge_work(1)?;
    if let Some(index) = indices.get(&source.point) {
        return Ok(*index);
    }
    let index = vertices.len();
    indices.insert(source.point.clone(), index);
    vertices.push(n::NumericalVertex {
        key: source.name.clone(),
        position: source.point.clone(),
    });
    Ok(index)
}

fn planar_definition(
    solid: &Solid,
    dec: &mut Decider,
) -> Result<n::NumericalDefinition, ImportError> {
    let mut definition = n::NumericalDefinition {
        name: solid.name.clone(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for face in &solid.faces {
        dec.charge_work(1)?;
        let crate::Surface::Plane(plane) = &face.surface else {
            return Err(unsupported("this numerical family requires planar faces"));
        };
        let mut frame = plane.clone();
        let nonzero: Vec<_> = (0..3).filter(|i| !frame.axis[*i].is_zero()).collect();
        if nonzero.len() == 1 {
            // A coordinate plane has a cheap canonical unit frame. This changes
            // only its parameterization, not the retained world-space geometry.
            let axis = nonzero[0];
            let sign = frame.axis[axis].sign();
            frame.axis = std::array::from_fn(|_| Rational::zero());
            frame.axis[axis] = if sign == Sign::Positive {
                Rational::one()
            } else {
                -Rational::one()
            };
            frame.ref_dir = std::array::from_fn(|_| Rational::zero());
            frame.ref_dir[(axis + 1) % 3] = Rational::one();
        }
        let plane = &frame;
        // ElSLib::PlaneParameters uses the inverse orthonormal placement.
        // Keep the accepted plane's actual frame, including oblique directions;
        // retain 128-bit enclosure midpoints as numerical UV coordinates and let
        // N01 check their actual world-space residuals at the fixed tolerance.
        let y = crate::planar_arithmetic::cross(&plane.axis, &plane.ref_dir, dec)?;
        let axes = [&plane.ref_dir, &y];
        let squared = [
            crate::planar_arithmetic::dot(axes[0], axes[0], dec)?,
            crate::planar_arithmetic::dot(axes[1], axes[1], dec)?,
        ];
        if squared.iter().any(|q| q.sign() != Sign::Positive) {
            return Err(unsupported(
                "the retained repaired plane has a degenerate parameter frame",
            ));
        }
        let lengths = squared.map(|q| Real::from(q).sqrt());
        let uv = |p: &Vec3, dec: &mut Decider| -> Result<_, ImportError> {
            let delta = crate::planar_arithmetic::sub(p, &plane.origin, dec)?;
            let mut result = [Rational::zero(), Rational::zero()];
            for i in 0..2 {
                let projection = crate::planar_arithmetic::dot(&delta, axes[i], dec)?;
                let value = &Real::from(projection) / &lengths[i];
                result[i] = dec.enclose_budgeted_at(&value, 128)?.mid();
            }
            Ok(result)
        };
        let mut loops = Vec::new();
        let mut boundaries: Vec<_> = face.loops.iter().collect();
        boundaries.sort_by_key(|bound| !bound.outer);
        for bound in boundaries {
            if bound.vertex.is_some() {
                return Err(unsupported(
                    "vertex-loop faces are outside this numerical family",
                ));
            }
            for chain in crate::build::chains(&face.name, bound.orientation, &bound.edges)? {
                let mut trims = Vec::new();
                for (edge, forward) in chain {
                    if !matches!(edge.curve, crate::Curve::Line(..)) {
                        return Err(unsupported(
                            "this numerical family requires finite straight boundaries",
                        ));
                    }
                    let (a, b) = if forward {
                        (&edge.start, &edge.end)
                    } else {
                        (&edge.end, &edge.start)
                    };
                    let ai = vertex(&mut definition.vertices, &mut vertices, a, dec)?;
                    let bi = vertex(&mut definition.vertices, &mut vertices, b, dec)?;
                    let pair = [ai.min(bi), ai.max(bi)];
                    let index = if let Some(index) = edges.get(&pair) {
                        *index
                    } else {
                        let index = definition.edges.len();
                        let start = &definition.vertices[pair[0]].position;
                        let end = &definition.vertices[pair[1]].position;
                        let direction = crate::planar_arithmetic::sub(end, start, dec)?;
                        definition.edges.push(n::NumericalEdge {
                            key: edge.edge.clone(),
                            curve: Curve::Line(Line {
                                origin: point(start),
                                dir: vector(&direction),
                            }),
                            vertices: pair,
                            parameters: [Rational::zero(), Rational::one()],
                        });
                        edges.insert(pair, index);
                        index
                    };
                    trims.push(n::NumericalTrim {
                        edge: index,
                        sense: if ai == pair[0] {
                            Sense::Same
                        } else {
                            Sense::Reversed
                        },
                        uv: [uv(&a.point, dec)?, uv(&b.point, dec)?],
                    });
                }
                loops.push(trims);
            }
        }
        definition.faces.push(n::NumericalFace {
            key: face.name.clone(),
            surface: Surface::Plane(Plane {
                origin: point(&plane.origin),
                normal: vector(&plane.axis),
                ref_dir: vector(&plane.ref_dir),
            }),
            sense: if face.same_sense {
                Sense::Same
            } else {
                Sense::Reversed
            },
            loops,
        });
    }
    Ok(definition)
}
