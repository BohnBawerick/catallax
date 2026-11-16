//! Transactional direct edits of admitted imported solids.
//!
//! Entity keys stay scoped to the retained STEP SHA. Geometry construction can
//! be rational, but every resulting material body is separately numerical.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use catallax_numbers::{
    Decider, Rational, Refusal, operation_budget_active, with_operation_budget,
};
use catallax_topology::numerical as n;

use crate::numerical::{ImportError, SewnSurfaceBody, sew_constructed_surfaces};
use crate::planar_arithmetic as a;
use crate::tolerant::ToleranceOptions;
use crate::tolerant::opened::{BodyView, OpenedStep};
use crate::{Curve, Face, Placement, Surface, Vec3, Vertex};

mod removal;
mod support;

/// One requested material change, in millimetres and original entity keys.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Edit {
    /// Remove a four-sided face and extend its two side supports to meet.
    DeleteFace {
        /// Current source-scoped face key.
        face: String,
    },
    /// Remove a four-sided cylindrical blend between planar side supports.
    DeleteFillet {
        /// Current source-scoped cylindrical blend key.
        face: String,
    },
    /// Remove the connected coaxial cavity containing this inward wall.
    DeleteHole {
        /// One inward cylindrical wall of the connected cavity.
        face: String,
    },
    /// Translate a support and reconstruct its intersections with neighbours.
    MoveFace {
        /// Current source-scoped face key.
        face: String,
        /// Requested support translation, not a healing displacement.
        translation_mm: Vec3,
    },
    /// Positive distance follows the face's material-outward normal.
    OffsetFace {
        /// Current source-scoped face key.
        face: String,
        /// Signed material-outward normal displacement.
        distance_mm: Rational,
    },
    /// Replace the support, retaining and reconstructing its finite boundary.
    ReplaceFace {
        /// Current source-scoped face key.
        face: String,
        /// New geometric support; adjacent trims are reconstructed.
        surface: Surface,
    },
}

/// Accuracy of the changed supports and rational endpoint construction only.
/// This never upgrades the resulting body to exact material authority.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Construction {
    /// Rational support changes and exact intersections of retained supports.
    Rational,
    /// Upper bound on the requested parallel support displacement error.
    /// This is not a global source-to-result distance or topology certificate.
    BoundedOffset {
        /// Absolute whole-plane displacement error, in millimetres.
        error_upper_mm: Rational,
    },
}

/// Face correspondence for one completed edit. Failures publish no record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EditRecord {
    /// Request that produced this complete derivative.
    pub edit: Edit,
    /// Support construction accuracy, separate from body admission.
    pub construction: Construction,
    /// Retained face keys whose finite definitions changed.
    pub modified_faces: Vec<String>,
    /// Source-scoped face keys removed by this edit.
    pub deleted_faces: Vec<String>,
}

#[derive(Debug)]
struct Source {
    opened: Arc<OpenedStep>,
    body_index: usize,
    sha256: String,
    body_name: String,
    faces: Arc<[Face]>,
}

/// Retained admitted source and complete, separately checked edit derivatives.
/// Clone shares immutable geometry; editing never mutates an earlier result.
#[derive(Clone, Debug)]
pub struct DirectSolid {
    source: Arc<Source>,
    faces: Arc<[Face]>,
    material: Option<Arc<SewnSurfaceBody>>,
    history: Vec<EditRecord>,
}

impl DirectSolid {
    /// Start from the actual admitted imported solid, not a closed declaration.
    /// Keep its existing material and proof without rebuilding it from retained
    /// pre-repair declarations. Only a requested edit constructs new material.
    pub fn from_opened(
        opened: Arc<OpenedStep>,
        index: usize,
        dec: &mut Decider,
    ) -> Result<Self, ImportError> {
        if !operation_budget_active() {
            return with_operation_budget(n::DEFAULT_WORK, || {
                Self::from_opened(opened, index, dec)
            })
            .0;
        }
        match opened.body(index) {
            Some(BodyView::Certified { model, body })
                if model.body(body).kind == catallax_topology::BodyKind::Solid => {}
            Some(BodyView::Numerical(body))
                if body.kind() == catallax_topology::BodyKind::Solid => {}
            _ => {
                return Err(bad(
                    "direct.not-solid",
                    "direct edits require an admitted imported solid",
                ));
            }
        }
        let mut file = opened.source().body_definition(index, dec)?;
        if file.solids.len() != 1 || !file.surface_models.is_empty() {
            return Err(bad(
                "direct.body-kind",
                "one source solid with one shell is required",
            ));
        }
        let solid = file.solids.pop().expect("checked one solid");
        let faces: Arc<[Face]> = solid.faces.into();
        let source = Arc::new(Source {
            sha256: opened.source().source_sha256().to_owned(),
            body_name: solid.name,
            faces: faces.clone(),
            opened,
            body_index: index,
        });
        Ok(Self {
            source,
            faces,
            material: None,
            history: Vec::new(),
        })
    }

    /// SHA-256 of the original imported bytes, never a derivative fingerprint.
    pub fn source_sha256(&self) -> &str {
        &self.source.sha256
    }
    /// Original source-scoped solid key.
    pub fn source_body_name(&self) -> &str {
        &self.source.body_name
    }
    /// Original retained definitions, in millimetres, without edit mutation.
    pub fn source_faces(&self) -> &[Face] {
        &self.source.faces
    }
    /// Complete current derivative definitions.
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }
    /// Original admitted material before editing; numerical material afterwards.
    /// A repaired certified input retains its original repair marks and evidence.
    pub fn body(&self) -> BodyView<'_> {
        match &self.material {
            Some(material) => BodyView::Numerical(material.body()),
            None => self.source.opened.body(self.source.body_index).expect("retained admitted source"),
        }
    }
    /// Sewing, source-pair and finite reconstruction evidence for a derivative.
    /// The initial editor retains the original OpenedStep admission instead.
    pub fn admission(&self) -> Option<&SewnSurfaceBody> {
        self.material.as_deref()
    }
    /// Completed edits in order; refused requests never enter this history.
    pub fn history(&self) -> &[EditRecord] {
        &self.history
    }

    /// Construct privately, sew, rebuild finite charts, and admit the complete
    /// result before publishing it. A refusal leaves both source and input intact.
    pub fn apply(
        &self,
        edit: Edit,
        sewing: ToleranceOptions,
        options: n::NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, ImportError> {
        if !operation_budget_active() {
            return with_operation_budget(n::DEFAULT_WORK, || {
                self.apply(edit, sewing, options, dec)
            })
            .0;
        }
        charge_copy(&self.faces, dec)?;
        let mut faces = self.faces.to_vec();
        let construction = match &edit {
            Edit::DeleteFace { face } => {
                removal::face(&mut faces, face, false, &options.tolerance_mm, dec)?
            }
            Edit::DeleteFillet { face } => {
                removal::face(&mut faces, face, true, &options.tolerance_mm, dec)?
            }
            Edit::DeleteHole { face } => removal::hole(&mut faces, face, dec)?,
            Edit::MoveFace {
                face,
                translation_mm,
            } => support::move_face(&mut faces, face, translation_mm, dec)?,
            Edit::OffsetFace { face, distance_mm } => {
                support::offset(&mut faces, face, distance_mm, &options, dec)?
            }
            Edit::ReplaceFace { face, surface } => {
                support::replace(&mut faces, face, surface.clone(), dec)?
            }
        };
        let by_name: BTreeMap<_, _> = faces.iter().map(|f| (f.name.as_str(), f)).collect();
        let mut modified_faces = Vec::new();
        let mut deleted_faces = Vec::new();
        for previous in self.faces.iter() {
            dec.charge_work(1)?;
            match by_name.get(previous.name.as_str()) {
                None => deleted_faces.push(previous.name.clone()),
                Some(next) if previous != *next => {
                    if let (Some(old), Some(new)) =
                        (polygon_winding(previous, dec)?, polygon_winding(next, dec)?)
                    {
                        if old != new || new == catallax_numbers::Sign::Zero {
                            return Err(bad(
                                "direct.topology-change",
                                "the edit reverses or collapses an adjacent planar boundary",
                            ));
                        }
                    }
                    modified_faces.push(previous.name.clone());
                }
                _ => {}
            }
        }
        if modified_faces.is_empty() && deleted_faces.is_empty() {
            return Err(bad(
                "direct.no-change",
                "the request does not change the retained geometry",
            ));
        }
        let material =
            sew_constructed_surfaces(&self.source.body_name, &faces, sewing, options, dec)?;
        let mut history = self.history.clone();
        history.push(EditRecord {
            edit,
            construction,
            modified_faces,
            deleted_faces,
        });
        dec.check_deadline()?;
        Ok(Self {
            source: self.source.clone(),
            faces: faces.into(),
            material: Some(Arc::new(material)),
            history,
        })
    }
}

fn bad(code: &'static str, detail: &str) -> ImportError {
    Refusal::new(code, detail).into()
}

fn index(faces: &[Face], key: &str) -> Result<usize, ImportError> {
    faces.iter().position(|f| f.name == key).ok_or_else(|| {
        bad(
            "direct.face-key",
            "the selected face is absent from the current solid",
        )
    })
}

fn plane(surface: &Surface) -> Result<&Placement, ImportError> {
    match surface {
        Surface::Plane(p) => Ok(p),
        _ => Err(bad(
            "direct.neighbourhood",
            "this intersection requires planar adjacent supports",
        )),
    }
}

fn parallel(a0: &Vec3, b: &Vec3, dec: &mut Decider) -> Result<bool, ImportError> {
    Ok(crate::exact::is_zero(&a::cross(a0, b, dec)?))
}

fn triple(p: [&Placement; 3], dec: &mut Decider) -> Result<Vec3, ImportError> {
    let cross = [
        a::cross(&p[1].axis, &p[2].axis, dec)?,
        a::cross(&p[2].axis, &p[0].axis, dec)?,
        a::cross(&p[0].axis, &p[1].axis, dec)?,
    ];
    let det = a::dot(&p[0].axis, &cross[0], dec)?;
    if det.is_zero() {
        return Err(bad(
            "direct.intersection",
            "three supports do not have a unique finite intersection",
        ));
    }
    let mut result = [Rational::zero(), Rational::zero(), Rational::zero()];
    for i in 0..3 {
        let scale = a::arithmetic(&a::dot(&p[i].axis, &p[i].origin, dec)?, &det, '/', dec)?;
        result = a::offset(&result, &cross[i], &scale, '+', dec)?;
    }
    Ok(result)
}

fn owners(
    faces: &[Face],
    dec: &mut Decider,
) -> Result<BTreeMap<String, BTreeSet<usize>>, ImportError> {
    let mut result: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        for bound in &face.loops {
            for edge in &bound.edges {
                dec.charge_work(1)?;
                result.entry(edge.edge.clone()).or_default().insert(i);
            }
        }
    }
    Ok(result)
}

fn vertex_owners(
    faces: &[Face],
    dec: &mut Decider,
) -> Result<BTreeMap<String, BTreeSet<usize>>, ImportError> {
    let mut result: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
    for (i, face) in faces.iter().enumerate() {
        for bound in &face.loops {
            for edge in &bound.edges {
                dec.charge_work(1)?;
                result.entry(edge.start.name.clone()).or_default().insert(i);
                result.entry(edge.end.name.clone()).or_default().insert(i);
            }
        }
    }
    Ok(result)
}

fn charge_copy(faces: &[Face], dec: &mut Decider) -> Result<(), ImportError> {
    // Reserve topology copying before allocation; sewing separately meters the
    // complete carrier data and performs its own input/resource checks.
    for face in faces {
        dec.charge_work(1 + face.name.len() as u64)?;
        for boundary in &face.loops {
            dec.charge_work(1 + boundary.edges.len() as u64 * 64)?;
        }
    }
    Ok(())
}

/// Update every use of shared vertices and reconstruct incident straight edges.
fn update_vertices(
    faces: &mut [Face],
    vertices: &BTreeMap<String, Vertex>,
    circles: &BTreeMap<String, Curve>,
    dec: &mut Decider,
) -> Result<(), ImportError> {
    for face in faces {
        for boundary in &mut face.loops {
            if let Some(v) = &mut boundary.vertex {
                if let Some(next) = vertices.get(&v.name) {
                    *v = next.clone();
                }
            }
            for edge in &mut boundary.edges {
                dec.charge_work(1)?;
                let changed = vertices.contains_key(&edge.start.name)
                    || vertices.contains_key(&edge.end.name);
                if let Some(v) = vertices.get(&edge.start.name) {
                    edge.start = v.clone();
                }
                if let Some(v) = vertices.get(&edge.end.name) {
                    edge.end = v.clone();
                }
                if let Some(curve) = circles.get(&edge.edge) {
                    edge.curve = curve.clone();
                } else if changed {
                    if !matches!(edge.curve, Curve::Line(..)) {
                        return Err(bad(
                            "direct.boundary-curve",
                            "a moved curved boundary needs a supported reconstructed carrier",
                        ));
                    }
                    let mut direction = a::sub(&edge.end.point, &edge.start.point, dec)?;
                    if crate::exact::is_zero(&direction) {
                        return Err(bad(
                            "direct.topology-change",
                            "the edit collapses an incident edge",
                        ));
                    }
                    if !edge.curve_same_sense {
                        for component in &mut direction {
                            *component = -&*component;
                        }
                    }
                    edge.curve = Curve::Line(edge.start.point.clone(), direction);
                }
            }
        }
    }
    Ok(())
}

fn polygon_winding(
    face: &Face,
    dec: &mut Decider,
) -> Result<Option<catallax_numbers::Sign>, ImportError> {
    let Surface::Plane(plane) = &face.surface else {
        return Ok(None);
    };
    let Some(boundary) = face.loops.iter().find(|b| b.outer) else {
        return Ok(None);
    };
    if boundary
        .edges
        .iter()
        .any(|e| !matches!(e.curve, Curve::Line(..)))
    {
        return Ok(None);
    }
    let mut area = [Rational::zero(), Rational::zero(), Rational::zero()];
    for edge in &boundary.edges {
        let (start, end) = if edge.forward {
            (&edge.start.point, &edge.end.point)
        } else {
            (&edge.end.point, &edge.start.point)
        };
        let cross = a::cross(start, end, dec)?;
        area = a::offset(&area, &cross, &Rational::one(), '+', dec)?;
    }
    let value = a::dot(&area, &plane.axis, dec)?;
    Ok(Some(if face.same_sense {
        value.sign()
    } else {
        (-value).sign()
    }))
}
