//! Exact application queries on retained B-rep bodies, independent of display.
//!
//! Inputs and source handles belong to one model. Sections and contact patches
//! own a private model containing their new definitions. Persist origins rather
//! than handles across edits. Refusals keep the underlying decision and input
//! counts. Queries spend the caller's active allowance, or the application
//! 100M default when called alone. They never publish a partial answer.

use std::collections::BTreeMap;

use catallax_boolean::{BoolOp, FaceClass, boolean, classify, imprint, intersect};
use catallax_features::Built;
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{
    Answer, Answered, Decider, Outcome, Rational, Real, Refusal, StrictSign,
    operation_budget_active, with_operation_budget,
};
use catallax_recipe::StepId;
use catallax_shapes::{Circle, CurveId, PointId, RPoint3, RVector3, Surface, Vector3};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};

use crate::{
    ApplicationBody, DEFAULT_BUILD_WORK,
    curves::{CurveSpan, EdgeCurve, edge_curve},
};

/// Size of one queried operand, available even when the query refuses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BodyCounts {
    /// Model-local operand.
    pub body: BodyId,
    /// Number of faces.
    pub faces: usize,
    /// Number of non-planar faces.
    pub curved_faces: usize,
    /// Number of edges.
    pub edges: usize,
    /// Number of vertices.
    pub vertices: usize,
}

/// A query stopped without replacing exact geometry by display samples.
#[derive(Clone, Debug)]
pub struct QueryRefusal {
    /// Stable underlying code, explanation and any undecided numeric bound.
    pub refusal: Refusal,
    /// Operand sizes in argument order.
    pub operands: Vec<BodyCounts>,
}

impl std::fmt::Display for QueryRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}; operands: {:?}", self.refusal, self.operands)
    }
}
impl std::error::Error for QueryRefusal {}

/// A plane with exact rational coordinates. Normal length has no meaning.
#[derive(Clone, Debug)]
pub struct SectionPlane {
    /// A point on the plane, in millimetres.
    pub origin: RPoint3,
    /// Nonzero normal, in model coordinates.
    pub normal: RVector3,
}

/// One exact piece of a section: a cut across one face, or an edge of the
/// body that lies in the plane.
#[derive(Clone, Debug)]
pub struct SectionArc {
    /// Carrier in [`SectionCurves::model`]: the exact intersection of a
    /// face's surface with the plane, or the edge's own curve.
    pub curve: CurveId,
    /// Start along the carrier's natural direction; `None` with `end` for
    /// the whole closed carrier. An edge with one vertex on a closed carrier
    /// starts and ends there.
    pub start: Option<PointId>,
    /// End along the carrier's natural direction.
    pub end: Option<PointId>,
    /// Faces of the body the piece lies on, in `faces_of_body` order: the
    /// one face it cuts across, or the faces of the edge.
    pub faces: Vec<FaceId>,
    /// The edge, for an edge of the body that lies in the plane.
    pub edge: Option<EdgeId>,
}

/// Trimmed section definitions and isolated contact points.
#[derive(Clone, Debug)]
pub struct SectionCurves {
    /// Owns every curve and point handle below. The input model is untouched;
    /// edge and face ids are the input body's.
    pub model: Model,
    /// Edges in the plane in edge order, then cuts in face order, each face's
    /// cuts in carrier order and then along the carrier.
    pub curves: Vec<SectionArc>,
    /// Isolated contacts, not silently discarded as an empty section: tangent
    /// points, and vertices that touch the plane off every curve.
    pub points: Vec<PointId>,
}

/// Orthographic outline generators in model coordinates, before projection.
#[derive(Clone, Debug)]
pub struct Silhouette {
    /// Exact visible boundary edges of a proven convex polyhedron.
    pub edges: Vec<EdgeCurve>,
    /// Smooth contours of whole spheres, paired with their source face.
    pub circles: Vec<(FaceId, Circle)>,
}

/// One exact shared boundary patch, split at the other body's trims.
#[derive(Clone, Debug)]
pub struct ContactPatch {
    /// Face in [`ContactFaces::model`]. Origins retain its source identity.
    pub face: FaceId,
    /// `OnSame` for equally oriented shared boundaries; `OnOpposite` for
    /// opposing boundaries. These are area contacts, not edge or point touches.
    pub orientation: FaceClass,
}

/// Shared boundary patches on both operands, including partial-face contacts.
#[derive(Clone, Debug)]
pub struct ContactFaces {
    /// Owns the clipped faces and their exact boundary definitions.
    pub model: Model,
    /// Patches of the first operand.
    pub first: Vec<ContactPatch>,
    /// Patches of the second operand, without assuming a one-to-one pairing.
    pub second: Vec<ContactPatch>,
}

/// Proven relation between two closed B-rep solids. Contact and undecided
/// distances return a [`QueryRefusal`], never an inferred clear answer.
#[derive(Clone, Debug)]
pub enum BodyRelation {
    /// Strictly positive overlap volume, in cubic millimetres.
    Clash {
        /// Enclosure with a strictly positive lower endpoint.
        overlap: Enclosure,
    },
    /// Global minimum distance between disjoint solids, in millimetres.
    /// The lower endpoint is strictly positive.
    Clear {
        /// Enclosure of the minimum over both complete solids.
        gap: Enclosure,
    },
}

pub(crate) fn query<T>(
    model: &Model,
    bodies: &[BodyId],
    run: impl FnOnce() -> Answered<T>,
) -> Result<Answer<T>, QueryRefusal> {
    let result = if operation_budget_active() {
        run()
    } else {
        let (result, work) = with_operation_budget(DEFAULT_BUILD_WORK, run);
        if work.exhausted && result.is_ok() {
            Err(Refusal::new(
                "application.work-limit",
                "aggregate query work exhausted",
            ))
        } else {
            result
        }
    };
    result.map_err(|refusal| QueryRefusal {
        refusal,
        operands: bodies
            .iter()
            .map(|&body| {
                let faces = model.faces_of_body(body);
                BodyCounts {
                    body,
                    faces: faces.len(),
                    curved_faces: faces
                        .iter()
                        .filter(|&&face| {
                            !matches!(
                                model.geom().surface(model.face(face).surface),
                                Surface::Plane(_)
                            )
                        })
                        .count(),
                    edges: model.edges_of_body(body).len(),
                    vertices: model.vertices_of_body(body).len(),
                }
            })
            .collect(),
    })
}

pub(crate) fn sign(value: &Real, dec: &mut Decider) -> Result<i8, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(1),
        Outcome::Proven(StrictSign::Negative) => Ok(-1),
        Outcome::ProvenCoincident(_) => Ok(0),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    }
}

/// Section a body by a rational plane, face by face (`views/section.rs`,
/// OpenCASCADE's `BOPAlgo_Section` pattern against an unbounded plane). Each
/// face's surface meets the plane in exact curves, cut where the face
/// boundary crosses the plane and kept where the face holds them. A planar
/// face in the plane contributes its boundary edges; every edge of the body
/// that lies in the plane appears once, under its own id. Intersection and
/// membership refusals keep their own codes.
pub fn section_curves(
    model: &Model,
    body: BodyId,
    plane: &SectionPlane,
    dec: &mut Decider,
) -> Result<Answer<SectionCurves>, QueryRefusal> {
    query(model, &[body], || {
        dec.check_deadline()?;
        if plane.normal.is_zero() {
            return Err(Refusal::new("app.section-normal", "section normal is zero"));
        }
        crate::views::section::plane_section(model, body, &plane.origin, &plane.normal, dec)
    })
}

/// Exact orthographic outline for a single convex planar solid or a whole
/// sphere. Direction points toward the viewer. Every polyhedral vertex must
/// lie in every outward face halfspace; nonconvex and other curved bodies
/// refuse rather than returning hidden or incomplete contour edges.
pub fn silhouette(
    model: &Model,
    body: BodyId,
    direction: &RVector3,
    dec: &mut Decider,
) -> Result<Answer<Silhouette>, QueryRefusal> {
    query(model, &[body], || {
        dec.check_deadline()?;
        if direction.is_zero() {
            return Err(Refusal::new(
                "app.silhouette-direction",
                "view direction is zero",
            ));
        }
        let faces = model.faces_of_body(body);
        let empty = || Silhouette {
            edges: Vec::new(),
            circles: Vec::new(),
        };
        let marks = model.body(body).marks.clone();
        if faces.is_empty() {
            return Ok(Answer::with_marks(empty(), marks));
        }
        let lumps = &model.body(body).lumps;
        if lumps.len() != 1 || model.lump(lumps[0]).shells.len() != 1 {
            return Err(Refusal::new(
                "app.silhouette-occlusion",
                "outline needs one convex shell; occlusion between shells is unfinished",
            ));
        }
        let direction = direction.to_real();
        if faces.len() == 1 {
            let face = model.face(faces[0]);
            if face.loops.is_empty()
                && face.sense == Sense::Same
                && let Surface::Sphere(sphere) = model.geom().surface(face.surface)
            {
                let reference = perpendicular(&direction, dec)?;
                let circle = Circle {
                    center: sphere.center.clone(),
                    normal: direction,
                    radius_sq: sphere.radius.square(),
                    ref_dir: reference,
                };
                return Ok(Answer::with_marks(
                    Silhouette {
                        edges: Vec::new(),
                        circles: vec![(faces[0], circle)],
                    },
                    marks,
                ));
            }
        }
        let vertices = model.vertices_of_body(body);
        let mut uses: BTreeMap<EdgeId, Vec<i8>> = BTreeMap::new();
        for &face in &faces {
            dec.check_deadline()?;
            let stored = model.face(face);
            let Surface::Plane(plane) = model.geom().surface(stored.surface) else {
                return Err(Refusal::new(
                    "app.silhouette-surface",
                    "smooth contour trimming is unfinished for this surface",
                ));
            };
            let normal = plane
                .normal
                .scaled(&Real::from(if stored.sense == Sense::Same {
                    1
                } else {
                    -1
                }));
            for &vertex in &vertices {
                let point = model.geom().position(model.vertex(vertex).point);
                if sign(&point.minus(&plane.origin).dot(&normal), dec)? > 0 {
                    return Err(Refusal::new(
                        "app.silhouette-nonconvex",
                        "a vertex lies outside an outward face halfspace; visibility splitting is unfinished",
                    ));
                }
            }
            let facing = sign(&normal.dot(&direction), dec)?;
            for &lp in &stored.loops {
                let Loop::Coedges(coedges) = model.loop_(lp) else {
                    return Err(Refusal::new(
                        "app.silhouette-boundary",
                        "convex planar contour requires finite line boundaries",
                    ));
                };
                for &coedge in coedges {
                    uses.entry(model.coedge(coedge).edge)
                        .or_default()
                        .push(facing);
                }
            }
        }
        let mut result = empty();
        for (edge, facing) in uses {
            let curve = edge_curve(model, body, edge, dec)?;
            if !matches!(curve.span, CurveSpan::Line) || facing.len() != 2 {
                return Err(Refusal::new(
                    "app.silhouette-boundary",
                    "convex planar contour requires two face uses per finite line",
                ));
            }
            if facing.contains(&1) && facing.iter().any(|&f| f <= 0) {
                result.edges.push(curve);
            }
        }
        Ok(Answer::with_marks(result, marks))
    })
}

pub(crate) fn perpendicular(normal: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    for axis in [
        Vector3::new(Real::one(), Real::zero(), Real::zero()),
        Vector3::new(Real::zero(), Real::one(), Real::zero()),
    ] {
        let cross = normal.cross(&axis);
        if sign(&cross.norm_squared(), dec)? != 0 {
            return Ok(cross);
        }
    }
    Err(Refusal::new(
        "app.silhouette-direction",
        "view direction is zero",
    ))
}

/// Proved overlap volume by regularised exact intersection and the existing
/// volume integrator. Empty intersection is a successful zero, never an error.
/// Volume is an enclosure at the requested accuracy, with all marks retained.
pub fn overlap_volume(
    model: &Model,
    first: BodyId,
    second: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Result<Answer<Enclosure>, QueryRefusal> {
    query(model, &[first, second], || {
        dec.check_deadline()?;
        let mut scratch = model.clone();
        let result =
            boolean(&mut scratch, first, second, BoolOp::Intersect, dec).map_err(|e| e.refusal)?;
        let mut volume = catallax_measure::volume(&scratch, result.body, accuracy, dec)?;
        volume.marks.extend(&result.marks);
        Ok(volume)
    })
}

/// Classify retained, already placed bodies using the existing exact overlap
/// query, then the minimum distance between their closed boundaries
/// ([`catallax_measure::solid_gap`]). Zero overlap must be proved before that
/// query's noncrossing precondition holds; solids whose proven boxes are
/// strictly apart have zero overlap without the Boolean. Empty operands and
/// contact retain typed refusals. A gap whose lower bound is proven positive
/// is clear; an enclosure wider than the caller's accuracy carries the
/// `measure.accuracy-not-reached` mark, the bounded answer kind.
///
/// Where the overlap Boolean refuses or leaves the overlap unresolved, the
/// pair is decided from boundary samples instead: a boundary point strictly
/// inside the other solid proves a clash, boundaries proven apart with every
/// shell outside the other solid prove a clear pair, and one solid wholly
/// inside the other overlaps by its own volume. A clash proven from a sample
/// encloses the overlap between an inscribed ball and the boxes' overlap,
/// the bounded kind. A pair still undecided keeps the Boolean's refusal.
pub fn body_relation(
    model: &Model,
    first: BodyId,
    second: BodyId,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Result<Answer<BodyRelation>, QueryRefusal> {
    query(model, &[first, second], || {
        if accuracy.absolute < Rational::zero() || accuracy.relative < Rational::zero() {
            return Err(Refusal::new(
                "app.relation-accuracy",
                "accuracy must be nonnegative",
            ));
        }
        let first_faces = model.faces_of_body(first);
        let second_faces = model.faces_of_body(second);
        if first_faces.is_empty() || second_faces.is_empty() {
            return Err(Refusal::new(
                "app.relation-empty",
                "empty operand has no minimum gap",
            ));
        }
        for body in [first, second] {
            if let Some(mark) = model.body(body).marks.iter().find(|mark| {
                mark.kind != catallax_numbers::BoundKind::Value && !mark.bound.is_zero()
            }) {
                return Err(Refusal::new(
                    "app.relation-construction-bound",
                    format!(
                        "{} has {:?} bound {}; relation to the original definition is not certified",
                        mark.code,
                        mark.kind,
                        mark.bound.to_decimal_bounds(10).1
                    ),
                ));
            }
        }
        // Solids whose proven boxes are strictly apart share no point, so their
        // overlap is zero without a Boolean.
        let overlap = if boxes_apart(model, first, second, dec) {
            let mut marks = model.body(first).marks.clone();
            marks.extend(&model.body(second).marks);
            Answer::with_marks(
                Enclosure {
                    lo: Rational::zero(),
                    hi: Rational::zero(),
                },
                marks,
            )
        } else {
            // Where the Boolean refuses, the arrangement search answers.
            match overlap_volume(model, first, second, accuracy, dec) {
                Ok(overlap) => overlap,
                Err(refused) => {
                    return arrangement::relation(
                        model,
                        first,
                        second,
                        accuracy,
                        refused.refusal,
                        dec,
                    );
                }
            }
        };
        if overlap.value.lo > Rational::zero() {
            return Ok(Answer::with_marks(
                BodyRelation::Clash {
                    overlap: overlap.value,
                },
                overlap.marks,
            ));
        }
        if !overlap.value.hi.is_zero() {
            let refused = Refusal::new(
                "app.relation-overlap",
                format!("overlap volume is unresolved in {:?} mm^3", overlap.value),
            );
            return arrangement::relation(model, first, second, accuracy, refused, dec);
        }
        gap_relation(model, first, second, overlap.marks, None, accuracy, dec)
    })
}

/// Clear by the minimum gap between two solids proven not to overlap, or a
/// typed refusal for contact or a gap not proven positive. `floor` is a lower
/// bound of the gap proven by another route, such as faces proven apart.
fn gap_relation(
    model: &Model,
    first: BodyId,
    second: BodyId,
    mut marks: catallax_numbers::MarkLog,
    floor: Option<&Rational>,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Answered<BodyRelation> {
    let distance = catallax_measure::solid_gap(model, first, second, accuracy, dec)?;
    marks.extend(&distance.marks);
    let mut gap = distance.value;
    if let Some(floor) = floor
        && *floor > gap.lo
    {
        gap.lo = floor.clone().min(gap.hi.clone());
    }
    if gap.hi.is_zero() {
        return Err(Refusal::new(
            "assembly.contact",
            "proven minimum gap is 0 mm",
        ));
    }
    if gap.lo <= Rational::zero() {
        return Err(Refusal::new(
            "app.relation-gap",
            format!(
                "minimum gap is not proven positive: [{}, {}] mm",
                gap.lo.to_decimal_bounds(6).0,
                gap.hi.to_decimal_bounds(6).1
            ),
        ));
    }
    // A proven positive gap is a clear answer. An enclosure wider than
    // requested is the bounded kind and carries the accuracy mark.
    mark_width(
        &mut marks,
        &gap,
        accuracy,
        "minimum gap width exceeds the requested accuracy",
    );
    Ok(Answer::with_marks(BodyRelation::Clear { gap }, marks))
}

/// Adds the `measure.accuracy-not-reached` mark once when a positive
/// enclosure is wider than the caller's accuracy allows at its lower end.
fn mark_width(
    marks: &mut catallax_numbers::MarkLog,
    value: &Enclosure,
    accuracy: &Accuracy,
    detail: &str,
) {
    let allowed = accuracy
        .absolute
        .clone()
        .max(&accuracy.relative * &value.lo);
    if value.width() > allowed
        && !marks
            .iter()
            .any(|mark| mark.code == "measure.accuracy-not-reached")
    {
        marks.push(catallax_numbers::Mark::new(
            "measure.accuracy-not-reached",
            catallax_numbers::BoundKind::Value,
            value.width(),
            detail,
        ));
    }
}

/// Accuracy of the proven boxes that relation queries compare: 1e-6 mm, the
/// box accuracy of `catallax_measure::distance::PreparedFace`.
fn box_accuracy() -> Accuracy {
    Accuracy {
        relative: Rational::zero(),
        absolute: Rational::from_ratio(1, 1_000_000).unwrap(),
    }
}

/// Whether two solids' proven boxes are strictly apart along some axis.
fn boxes_apart(model: &Model, a: BodyId, b: BodyId, dec: &mut Decider) -> bool {
    let accuracy = box_accuracy();
    let [Some(x), Some(y)] = [a, b].map(|body| {
        catallax_measure::bounding_box(model, body, &accuracy, dec)
            .ok()
            .filter(|bounds| bounds.marks.is_empty())
            .and_then(|bounds| bounds.value)
    }) else {
        return false;
    };
    (0..3).any(|i| x.max[i] < y.min[i] || y.max[i] < x.min[i])
}

/// Exact shared area patches through intersect, imprint and classify. Edge and
/// point touches produce no contact faces. No overlap volume is needed, so a
/// measurement refusal does not hide a decidable boundary contact.
pub fn contact_faces(
    model: &Model,
    first: BodyId,
    second: BodyId,
    dec: &mut Decider,
) -> Result<Answer<ContactFaces>, QueryRefusal> {
    query(model, &[first, second], || {
        dec.check_deadline()?;
        let mut scratch = model.clone();
        let mut marks = model.body(first).marks.clone();
        marks.extend(&model.body(second).marks);
        if first == second {
            let patches = model
                .faces_of_body(first)
                .into_iter()
                .map(|face| ContactPatch {
                    face,
                    orientation: FaceClass::OnSame,
                })
                .collect::<Vec<_>>();
            return Ok(Answer::with_marks(
                ContactFaces {
                    model: scratch,
                    first: patches.clone(),
                    second: patches,
                },
                marks,
            ));
        }
        let graph = intersect(&mut scratch, first, second, dec)?;
        let split = imprint(&mut scratch, &graph.value, dec)?;
        let classes = classify(&mut scratch, &split.value, dec)?;
        marks.extend(&graph.marks);
        marks.extend(&split.marks);
        marks.extend(&classes.marks);
        let patches = |faces: Vec<FaceId>| {
            faces
                .into_iter()
                .filter_map(|face| {
                    let orientation = classes.value.class_of(face)?;
                    matches!(orientation, FaceClass::OnSame | FaceClass::OnOpposite)
                        .then_some(ContactPatch { face, orientation })
                })
                .collect()
        };
        let first = patches(split.value.pieces_of_first());
        let second = patches(split.value.pieces_of_second());
        Ok(Answer::with_marks(
            ContactFaces {
                model: scratch,
                first,
                second,
            },
            marks,
        ))
    })
}

/// Why a B-rep relation must not use `step`'s body, if it must not. A lattice
/// or gyroid step is field-only: the B-rep carries its region's envelope in
/// its place (PRD section 7, [`catallax_features::StandIn`]), and so does
/// every step that reads it, directly or through other steps. The part's
/// material lies only on the beams or the sheet, so a gap, overlap or contact
/// measured on the envelope is not the part's: the envelope can state a
/// smaller gap than the part has, or an overlap the part does not have. The
/// field kernel holds the material. Mass refuses the same envelope as
/// `assembly.mass-envelope` ([`crate::assembly::uniform_mass`]).
pub fn field_only_envelope(built: &Built, step: &StepId) -> Option<Refusal> {
    let envelopes = envelope_steps(built);
    envelopes
        .get(step)
        .map(|stand_in| envelope_refusal(step, stand_in))
}

/// Each step whose body the B-rep carries as a field-only envelope, with the
/// field-only step it reads. A step reads only steps before it, so one pass
/// in recipe order finds every step that reads one through others.
fn envelope_steps(built: &Built) -> BTreeMap<StepId, StepId> {
    let mut envelopes: BTreeMap<StepId, StepId> = built
        .stand_ins
        .iter()
        .map(|s| (s.step.clone(), s.step.clone()))
        .collect();
    if envelopes.is_empty() {
        return envelopes;
    }
    for step in &built.recipe().steps {
        let read = step
            .op
            .inputs()
            .into_iter()
            .find_map(|(input, _)| envelopes.get(input).cloned());
        if let Some(stand_in) = read {
            envelopes.entry(step.id.clone()).or_insert(stand_in);
        }
    }
    envelopes
}

fn envelope_refusal(step: &StepId, stand_in: &StepId) -> Refusal {
    Refusal::new(
        "app.relation-envelope",
        format!(
            "{step} carries the field-only step {stand_in} as its region's envelope, not its material"
        ),
    )
}

impl ApplicationBody {
    /// The refusal for a relation between this result and `other` when
    /// either body is a field-only envelope ([`field_only_envelope`]). A
    /// lattice step's body is its region's body, so the region refuses too.
    fn envelope_operand(&self, other: BodyId) -> Option<Refusal> {
        let built = &self.built;
        envelope_steps(built)
            .iter()
            .find(|(step, _)| {
                built
                    .bodies
                    .get(*step)
                    .is_some_and(|&body| body == built.result || body == other)
            })
            .map(|(step, stand_in)| envelope_refusal(step, stand_in))
    }
    /// Refuses `refusal` with both operands' counts, as the query would.
    fn refused<T>(&self, other: BodyId, refusal: Refusal) -> Result<Answer<T>, QueryRefusal> {
        query(&self.built.model, &[self.built.result, other], || {
            Err(refusal)
        })
    }
    /// Clash or proven minimum gap to another retained, placed body. A
    /// field-only envelope on either side refuses `app.relation-envelope`.
    pub fn body_relation(
        &self,
        other: BodyId,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Answer<BodyRelation>, QueryRefusal> {
        if let Some(refusal) = self.envelope_operand(other) {
            return self.refused(other, refusal);
        }
        body_relation(&self.built.model, self.built.result, other, accuracy, dec)
    }
    /// Exact section of this result, independent of its display.
    pub fn section_curves(
        &self,
        plane: &SectionPlane,
        dec: &mut Decider,
    ) -> Result<Answer<SectionCurves>, QueryRefusal> {
        section_curves(&self.built.model, self.built.result, plane, dec)
    }
    /// Exact outline of this result for a direction toward the viewer.
    pub fn silhouette(
        &self,
        direction: &RVector3,
        dec: &mut Decider,
    ) -> Result<Answer<Silhouette>, QueryRefusal> {
        silhouette(&self.built.model, self.built.result, direction, dec)
    }
    /// Overlap with another body retained in this build's model. A
    /// field-only envelope on either side refuses `app.relation-envelope`.
    pub fn overlap_volume(
        &self,
        other: BodyId,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Answer<Enclosure>, QueryRefusal> {
        if let Some(refusal) = self.envelope_operand(other) {
            return self.refused(other, refusal);
        }
        overlap_volume(&self.built.model, self.built.result, other, accuracy, dec)
    }
    /// Shared boundary patches with another body in this build's model. A
    /// field-only envelope on either side refuses `app.relation-envelope`.
    pub fn contact_faces(
        &self,
        other: BodyId,
        dec: &mut Decider,
    ) -> Result<Answer<ContactFaces>, QueryRefusal> {
        if let Some(refusal) = self.envelope_operand(other) {
            return self.refused(other, refusal);
        }
        contact_faces(&self.built.model, self.built.result, other, dec)
    }
}

mod arrangement;
#[cfg(test)]
mod tests;
