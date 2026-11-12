//! Exact application queries on retained B-rep bodies, independent of display.
//!
//! Inputs and source handles belong to one model. Sections and contact patches
//! own a private model containing their new definitions. Persist origins rather
//! than handles across edits. Refusals keep the underlying decision and input
//! counts. Queries spend the caller's active allowance, or the application
//! default when called alone. They never publish a partial answer.

use std::collections::BTreeMap;

use catallax_boolean::{BoolOp, FaceClass, SectionEdge, boolean, classify, imprint, intersect};
use catallax_measure::{Accuracy, Enclosure};
use catallax_numbers::{
    Answer, Answered, Decider, Outcome, Rational, Real, Refusal, StrictSign,
    operation_budget_active, with_operation_budget,
};
use catallax_shapes::{Circle, PointId, Pose, RPoint3, RVector3, Rotation, Surface, Vector3};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense, block};

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

/// Trimmed section definitions and isolated tangent points.
#[derive(Clone, Debug)]
pub struct SectionCurves {
    /// Owns every curve and point handle below. The input model is untouched.
    pub model: Model,
    /// Exact trimmed branches in deterministic graph order. Source faces are
    /// in `faces_of_first`; the second faces belong to the private cutting tool.
    pub curves: Vec<SectionEdge>,
    /// Isolated contacts, not silently discarded as an empty section.
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

fn query<T>(
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

fn sign(value: &Real, dec: &mut Decider) -> Result<i8, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(1),
        Outcome::Proven(StrictSign::Negative) => Ok(-1),
        Outcome::ProvenCoincident(_) => Ok(0),
        Outcome::NotExact(n) => Err(Refusal::from(n)),
        Outcome::Refused(r) => Err(r),
    }
}

/// Section a body by a rational plane using the Boolean engine's trimmed
/// intersection graph. A bounded cutting tool encloses the body's certified
/// box strictly on its other five sides. Coplanar faces refuse explicitly.
/// Bounding, intersection and tangent-point refusals retain their own codes.
pub fn section_curves(
    model: &Model,
    body: BodyId,
    plane: &SectionPlane,
    accuracy: &Accuracy,
    dec: &mut Decider,
) -> Result<Answer<SectionCurves>, QueryRefusal> {
    query(model, &[body], || {
        dec.check_deadline()?;
        if plane.normal.is_zero() {
            return Err(Refusal::new("app.section-normal", "section normal is zero"));
        }
        let bounds = catallax_measure::bounding_box(model, body, accuracy, dec)?;
        let mut scratch = model.clone();
        let mut marks = bounds.marks;
        marks.extend(&model.body(body).marks);
        let Some(bounds) = bounds.value else {
            return Ok(Answer::with_marks(
                SectionCurves {
                    model: scratch,
                    curves: Vec::new(),
                    points: Vec::new(),
                },
                marks,
            ));
        };
        let n = plane.normal.to_real();
        let origin = plane.origin.to_real();
        for face in model.faces_of_body(body) {
            if let Surface::Plane(p) = model.geom().surface(model.face(face).surface)
                && sign(&p.normal.cross(&n).norm_squared(), dec)? == 0
                && sign(&origin.minus(&p.origin).dot(&n), dec)? == 0
            {
                return Err(Refusal::new(
                    "app.section-coplanar",
                    "section contains a planar face; region boundary extraction is unfinished",
                ));
            }
        }
        let coordinates = [&plane.origin.x, &plane.origin.y, &plane.origin.z];
        let radius = (0..3).fold(Rational::one(), |sum, i| {
            sum + (&bounds.min[i] - coordinates[i])
                .abs()
                .max((&bounds.max[i] - coordinates[i]).abs())
        });
        let unit = n.scaled(&(Real::one() / n.norm_squared().sqrt()));
        let rotation = if plane.normal.x.is_zero() && plane.normal.y.is_zero() {
            if plane.normal.z.sign() == catallax_numbers::Sign::Negative {
                Rotation::from_quaternion(
                    Rational::zero(),
                    Rational::one(),
                    Rational::zero(),
                    Rational::zero(),
                )?
            } else {
                Rotation::identity()
            }
        } else {
            Rotation::from_real_quaternion(
                &Real::one() + &unit.z,
                -&unit.y,
                unit.x.clone(),
                Real::zero(),
                dec,
            )?
        };
        let offset = rotation.apply(&Vector3::new(
            (-&radius).into(),
            (-&radius).into(),
            Real::zero(),
        ));
        let pose = Pose::new(
            rotation,
            origin
                .plus(&offset)
                .minus(&catallax_shapes::Point3::origin()),
            dec,
        )?;
        let tool = block(
            &mut scratch,
            [&radius + &radius, &radius + &radius, radius],
            &pose,
        )?;
        let graph = intersect(&mut scratch, body, tool, dec)?;
        marks.extend(&graph.marks);
        let curves = graph.value.edges().to_vec();
        let points = graph.value.points().to_vec();
        Ok(Answer::with_marks(
            SectionCurves {
                model: scratch,
                curves,
                points,
            },
            marks,
        ))
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

fn perpendicular(normal: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
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

impl ApplicationBody {
    /// Exact section of this result, independent of its display.
    pub fn section_curves(
        &self,
        plane: &SectionPlane,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Answer<SectionCurves>, QueryRefusal> {
        section_curves(&self.built.model, self.built.result, plane, accuracy, dec)
    }
    /// Exact outline of this result for a direction toward the viewer.
    pub fn silhouette(
        &self,
        direction: &RVector3,
        dec: &mut Decider,
    ) -> Result<Answer<Silhouette>, QueryRefusal> {
        silhouette(&self.built.model, self.built.result, direction, dec)
    }
    /// Overlap with another body retained in this build's model.
    pub fn overlap_volume(
        &self,
        other: BodyId,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<Answer<Enclosure>, QueryRefusal> {
        overlap_volume(&self.built.model, self.built.result, other, accuracy, dec)
    }
    /// Shared boundary patches with another body in this build's model.
    pub fn contact_faces(
        &self,
        other: BodyId,
        dec: &mut Decider,
    ) -> Result<Answer<ContactFaces>, QueryRefusal> {
        contact_faces(&self.built.model, self.built.result, other, dec)
    }
}

#[cfg(test)]
mod tests;
