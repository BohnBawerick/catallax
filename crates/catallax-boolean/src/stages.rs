//! The five reusable stages.

use catallax_numbers::{Answer, Answered, Decider, MarkLog, Outcome, StrictSign};
use catallax_shapes::{CurveId, PointId};
use catallax_topology::{BodyId, EdgeId, FaceId, Model, Sense};
use std::collections::{BTreeMap, BTreeSet};

use crate::BoolOp;

/// An arc of a curve that lies on the boundaries of both bodies.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SectionEdge {
    /// The curve carrying the arc.
    pub curve: CurveId,
    /// Start of the arc; `None` with `end` for a whole closed curve.
    pub start: Option<PointId>,
    /// End of the arc; `None` with `start` for a whole closed curve.
    pub end: Option<PointId>,
    /// Faces of the first body the arc lies on, in `faces_of_body` order.
    pub faces_of_first: Vec<FaceId>,
    /// Faces of the second body the arc lies on, in `faces_of_body` order.
    pub faces_of_second: Vec<FaceId>,
}

/// Where two bodies' boundaries meet.
#[derive(Clone, Debug)]
pub struct IntersectionGraph {
    pub(crate) a: BodyId,
    pub(crate) b: BodyId,
    pub(crate) edges: Vec<SectionEdge>,
    pub(crate) points: Vec<PointId>,
    pub(crate) singular: BTreeSet<PointId>,
    pub(crate) marks: MarkLog,
}

impl IntersectionGraph {
    /// The first body.
    pub fn first(&self) -> BodyId {
        self.a
    }

    /// The second body.
    pub fn second(&self) -> BodyId {
        self.b
    }

    /// Section edges, in a deterministic order.
    pub fn edges(&self) -> &[SectionEdge] {
        &self.edges
    }

    /// Isolated meeting points (on no section edge), in a deterministic order.
    pub fn points(&self) -> &[PointId] {
        &self.points
    }
}

/// Both bodies split along their intersection graph.
#[derive(Clone, Debug)]
pub struct Imprinted {
    pub(crate) a: BodyId,
    pub(crate) b: BodyId,
    pub(crate) first: Vec<FaceId>,
    pub(crate) second: Vec<FaceId>,
    pub(crate) common: Vec<EdgeId>,
    pub(crate) singular: BTreeSet<PointId>,
    pub(crate) marks: MarkLog,
}

impl Imprinted {
    /// Face pieces of the first body.
    pub fn pieces_of_first(&self) -> Vec<FaceId> {
        self.first.clone()
    }

    /// Face pieces of the second body.
    pub fn pieces_of_second(&self) -> Vec<FaceId> {
        self.second.clone()
    }

    /// Edges shared by pieces of both bodies.
    pub fn common_edges(&self) -> Vec<EdgeId> {
        self.common.clone()
    }
}

/// Class of a face piece against the other body.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub enum FaceClass {
    /// Inside the other body.
    Inside,
    /// Outside the other body.
    Outside,
    /// On a face of the other body with the same normal.
    OnSame,
    /// On a face of the other body with the opposite normal.
    OnOpposite,
}

/// Classes of every face piece.
#[derive(Clone, Debug)]
pub struct Classification {
    classes: BTreeMap<FaceId, FaceClass>,
    coincident: BTreeMap<FaceId, FaceId>,
    marks: MarkLog,
}

impl Classification {
    /// The class of a face piece, if it is one.
    pub fn class_of(&self, face: FaceId) -> Option<FaceClass> {
        self.classes.get(&face).copied()
    }
}

/// Face pieces chosen for the result, with the sense to use.
#[derive(Clone, Debug)]
pub struct Selection {
    pub(crate) faces: Vec<(FaceId, Sense)>,
    pub(crate) coincident: BTreeMap<FaceId, FaceId>,
    pub(crate) intersection: bool,
    pub(crate) marks: MarkLog,
}

impl Selection {
    /// The chosen pieces in a fixed order.
    pub fn faces(&self) -> Vec<(FaceId, Sense)> {
        self.faces.clone()
    }
}

/// Stage 1: find where the boundaries of `a` and `b` meet.
pub fn intersect(
    model: &mut Model,
    a: BodyId,
    b: BodyId,
    dec: &mut Decider,
) -> Answered<IntersectionGraph> {
    crate::graph::intersect(model, a, b, dec).map_err(|mut refusal| {
        if refusal.code == "boolean.classify-ambiguous" {
            refusal.code = "boolean.imprint-ambiguous";
        }
        refusal
    })
}

/// Stage 2: split faces and edges of both bodies along the graph.
pub fn imprint(
    model: &mut Model,
    graph: &IntersectionGraph,
    dec: &mut Decider,
) -> Answered<Imprinted> {
    crate::imprinting::imprint(model, graph, dec).map_err(|mut refusal| {
        if refusal.code == "boolean.classify-ambiguous" {
            refusal.code = "boolean.imprint-ambiguous";
        }
        refusal
    })
}

/// Stage 3: label every face piece against the other body.
pub fn classify(
    model: &mut Model,
    imprinted: &Imprinted,
    dec: &mut Decider,
) -> Answered<Classification> {
    let mut classes = BTreeMap::new();
    let mut coincident = BTreeMap::new();
    for (pieces, other) in [
        (&imprinted.first, imprinted.b),
        (&imprinted.second, imprinted.a),
    ] {
        let other_faces = model.faces_of_body(other);
        for &piece in pieces {
            let face = model.face(piece).clone();
            let mut counterpart = None;
            let class = crate::points::find(model, piece, dec, |model, point, rational, dec| {
                for &other_face in &other_faces {
                    let opposite = model.face(other_face).clone();
                    if crate::geometry::coincident(
                        model
                            .geom()
                            .same_surface(face.surface, opposite.surface, dec),
                        "boolean.classify-ambiguous",
                    )? && crate::region::locate(model, other_face, point, dec)?
                        == crate::region::Location::Inside
                    {
                        let agree = match model.geom().normals_agree(
                            face.surface,
                            opposite.surface,
                            point,
                            dec,
                        ) {
                            Outcome::Proven(StrictSign::Positive) => true,
                            Outcome::Proven(StrictSign::Negative) => false,
                            Outcome::Refused(e) => return Err(e),
                            _ => {
                                return Err(crate::geometry::refusal(
                                    "boolean.classify-ambiguous",
                                    "coincident face normals could not be ordered",
                                ));
                            }
                        };
                        let same = agree == (face.sense == opposite.sense);
                        if !model.face_origins(piece).is_empty()
                            || !model.face_origins(other_face).is_empty()
                        {
                            counterpart = Some(other_face);
                        }
                        return Ok(Some(if same {
                            FaceClass::OnSame
                        } else {
                            FaceClass::OnOpposite
                        }));
                    }
                }
                match crate::classify_point(model, other, rational, dec) {
                    Outcome::Proven(crate::Containment::Inside) => Ok(Some(FaceClass::Inside)),
                    Outcome::Proven(crate::Containment::Outside) => Ok(Some(FaceClass::Outside)),
                    Outcome::ProvenCoincident(_) => Ok(None),
                    Outcome::NotExact(n) => Err(crate::geometry::refusal(
                        "boolean.classify-ambiguous",
                        n.detail,
                    )),
                    Outcome::Refused(e) => Err(e),
                }
            })?;
            classes.insert(piece, class);
            if let Some(other) = counterpart {
                coincident.insert(piece, other);
            }
        }
    }
    Ok(Answer::exact(Classification {
        classes,
        coincident,
        marks: MarkLog::new(),
    }))
}

/// Stage 4: choose pieces by the operation's table.
pub fn keep(classification: &Classification, imprinted: &Imprinted, op: BoolOp) -> Selection {
    let mut faces = Vec::new();
    for &piece in &imprinted.first {
        let class = classification
            .class_of(piece)
            .expect("classification covers every imprinted piece");
        let selected = match op {
            BoolOp::Unite => matches!(class, FaceClass::Outside | FaceClass::OnSame),
            BoolOp::Intersect => matches!(class, FaceClass::Inside | FaceClass::OnSame),
            BoolOp::Subtract => matches!(class, FaceClass::Outside | FaceClass::OnOpposite),
        };
        if selected {
            faces.push((piece, Sense::Same));
        }
    }
    for &piece in &imprinted.second {
        let class = classification
            .class_of(piece)
            .expect("classification covers every imprinted piece");
        let selected = match op {
            BoolOp::Unite => class == FaceClass::Outside,
            BoolOp::Intersect | BoolOp::Subtract => class == FaceClass::Inside,
        };
        if selected {
            faces.push((
                piece,
                if op == BoolOp::Subtract {
                    Sense::Reversed
                } else {
                    Sense::Same
                },
            ));
        }
    }
    let mut marks = imprinted.marks.clone();
    marks.extend(&classification.marks);
    let coincident = faces
        .iter()
        .filter_map(|(face, _)| {
            classification
                .coincident
                .get(face)
                .map(|other| (*face, *other))
        })
        .collect();
    Selection {
        faces,
        coincident,
        intersection: op == BoolOp::Intersect,
        marks,
    }
}

/// Stage 5: join the chosen pieces into a new body.
pub fn sew(
    model: &mut Model,
    imprinted: &Imprinted,
    selection: &Selection,
    dec: &mut Decider,
) -> Answered<BodyId> {
    crate::sewing::sew(model, imprinted, selection, dec)
}

/// Remove redundant same-support boundaries and degree-two vertices from a
/// structurally valid body. Preserve its proven lump/void partition and origins.
/// Unlike a boolean, this requires no second solid or surface classification.
/// Shell orientation and nesting must already be valid; this does not repair
/// self-intersections or classify material.
pub fn make_maximal(model: &mut Model, body: BodyId, dec: &mut Decider) -> Answered<BodyId> {
    crate::sewing::make_maximal(model, body, dec)
}
