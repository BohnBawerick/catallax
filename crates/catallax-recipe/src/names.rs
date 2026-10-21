//! Names: what made a face, the edges and vertices between faces, and references.

use std::cmp::Ordering;
use std::fmt;

use serde::de::{self, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

use crate::ids::{ElementId, StepId};
use crate::number::Point3;

/// What made a face: the step and the part of its operation that made it.
///
/// Origins are ordered by their compact canonical JSON text, compared byte by byte
/// ([`Origin::key`]). Edge and vertex names list their origins in that order.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Origin {
    /// The step that made the face.
    pub step: StepId,
    /// The part of the step's operation that made it.
    pub role: Role,
}

impl Origin {
    /// The compact canonical JSON text of this origin. It defines the canonical order.
    pub fn key(&self) -> String {
        serde_json::to_string(self).expect("an origin always serialises")
    }
}

impl Ord for Origin {
    fn cmp(&self, other: &Self) -> Ordering {
        if self == other {
            Ordering::Equal
        } else {
            self.key().cmp(&other.key())
        }
    }
}

impl PartialOrd for Origin {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// The part of an operation that made a face. `docs/RECIPE.md` lists which operation
/// makes which roles.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum Role {
    /// The starting profile region, or the lower cap of a stock solid.
    Start,
    /// The ending profile region, or the upper cap of a stock solid when present.
    End,
    /// A swept profile segment or the fixed side segment of a stock solid.
    Side {
        /// The profile segment, or a stock solid's fixed segment id.
        segment: ElementId,
        /// The path segment, for a sweep only.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<ElementId>,
    },
    /// A face of a hole.
    Hole {
        /// The hole placement.
        at: ElementId,
        /// Which face of the hole.
        part: HolePart,
    },
    /// A fillet face rounding an edge.
    Blend {
        /// The edge it rounds.
        edge: EdgeName,
    },
    /// A chamfer face bevelling an edge.
    Bevel {
        /// The edge it bevels.
        edge: EdgeName,
    },
    /// A patch where fillets or chamfers meet at a vertex.
    Corner {
        /// The vertex.
        vertex: VertexName,
    },
    /// The inner wall of a shell, at the thickness from a kept face.
    Inner {
        /// The face it is offset from.
        of: Box<Origin>,
    },
    /// A face at a constant distance from an edge, made by a shell or an offset.
    Round {
        /// The edge.
        edge: EdgeName,
    },
    /// A face at a constant distance from a vertex, made by a shell or an offset.
    RoundCorner {
        /// The vertex.
        vertex: VertexName,
    },
    /// The far face of a thickened face.
    Offset {
        /// The thickened face.
        of: Box<Origin>,
    },
    /// The side wall of a thickened face set along one of its boundary edges.
    Rim {
        /// The boundary edge.
        edge: EdgeName,
    },
    /// A face of a copy made by a pattern or a mirror.
    Copy {
        /// Which copy: `k` of a pattern, 1 for a mirror image.
        instance: u32,
        /// The face it copies.
        of: Box<Origin>,
    },
    /// The surface of a lattice's beams.
    Beams,
    /// The surface of a gyroid.
    Surface,
    /// The blend surface of a smooth blend.
    Smooth,
    /// A face of an imported STEP body.
    Imported {
        /// The face's entity name in the imported file, such as `"#1234"`.
        key: String,
    },
    /// The surface of an imported mesh.
    Mesh,
}

impl Role {
    /// The role's name, as written in a recipe.
    pub fn name(&self) -> &'static str {
        match self {
            Role::Start => "start",
            Role::End => "end",
            Role::Side { .. } => "side",
            Role::Hole { .. } => "hole",
            Role::Blend { .. } => "blend",
            Role::Bevel { .. } => "bevel",
            Role::Corner { .. } => "corner",
            Role::Inner { .. } => "inner",
            Role::Round { .. } => "round",
            Role::RoundCorner { .. } => "round_corner",
            Role::Offset { .. } => "offset",
            Role::Rim { .. } => "rim",
            Role::Copy { .. } => "copy",
            Role::Beams => "beams",
            Role::Surface => "surface",
            Role::Smooth => "smooth",
            Role::Imported { .. } => "imported",
            Role::Mesh => "mesh",
        }
    }
}

/// The faces of a hole, from the entry down.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum HolePart {
    /// The cylinder above the placement point, where material lies above it.
    Entry,
    /// The wall of a counterbore.
    CounterboreWall,
    /// The flat floor of a counterbore.
    CounterboreFloor,
    /// The cone of a countersink.
    Countersink,
    /// The main wall.
    Wall,
    /// The flat floor of a blind hole with a flat tip.
    Floor,
    /// The cone of a blind hole with a drill-point tip.
    Tip,
}

impl HolePart {
    /// The part's name, as written in a recipe.
    pub fn name(self) -> &'static str {
        match self {
            HolePart::Entry => "entry",
            HolePart::CounterboreWall => "counterbore_wall",
            HolePart::CounterboreFloor => "counterbore_floor",
            HolePart::Countersink => "countersink",
            HolePart::Wall => "wall",
            HolePart::Floor => "floor",
            HolePart::Tip => "tip",
        }
    }
}

/// Why a name is malformed.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct NameError(&'static str);

impl fmt::Display for NameError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}

impl std::error::Error for NameError {}

/// The name of an edge: the origins of the two faces it separates, in canonical order.
/// A face with several origins contributes the first of them.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(try_from = "RawEdgeName")]
pub struct EdgeName {
    between: Box<[Origin; 2]>,
}

impl EdgeName {
    /// The edge between faces of origins `a` and `b`, in either order.
    pub fn new(a: Origin, b: Origin) -> Result<Self, NameError> {
        match a.cmp(&b) {
            Ordering::Less => Ok(EdgeName {
                between: Box::new([a, b]),
            }),
            Ordering::Greater => Ok(EdgeName {
                between: Box::new([b, a]),
            }),
            Ordering::Equal => Err(NameError("an edge name needs two different origins")),
        }
    }

    /// The two origins, in canonical order.
    pub fn between(&self) -> &[Origin; 2] {
        &self.between
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawEdgeName {
    between: [Origin; 2],
}

impl TryFrom<RawEdgeName> for EdgeName {
    type Error = NameError;
    fn try_from(raw: RawEdgeName) -> Result<Self, NameError> {
        let [a, b] = raw.between;
        EdgeName::new(a, b)
    }
}

/// The name of a vertex: the origins of the faces that meet there, in canonical order.
/// A face with several origins contributes the first of them.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(try_from = "RawVertexName")]
pub struct VertexName {
    faces: Vec<Origin>,
}

impl VertexName {
    /// The vertex where faces of these origins meet, in any order. A singular vertex
    /// such as a cone apex has one origin. Needs at least one origin and no repeats.
    pub fn new(mut faces: Vec<Origin>) -> Result<Self, NameError> {
        faces.sort();
        if faces.is_empty() {
            return Err(NameError("a vertex name needs at least one origin"));
        }
        if faces.windows(2).any(|w| w[0] == w[1]) {
            return Err(NameError("a vertex name lists an origin twice"));
        }
        Ok(VertexName { faces })
    }

    /// The origins, in canonical order.
    pub fn faces(&self) -> &[Origin] {
        &self.faces
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct RawVertexName {
    faces: Vec<Origin>,
}

impl TryFrom<RawVertexName> for VertexName {
    type Error = NameError;
    fn try_from(raw: RawVertexName) -> Result<Self, NameError> {
        VertexName::new(raw.faces)
    }
}

/// Matches origins in a reference: every origin made by `step`, or with `role` only
/// those with that role.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OriginPattern {
    /// The step that made the face.
    pub step: StepId,
    /// The role; absent matches every role of the step.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub role: Option<Role>,
}

impl OriginPattern {
    /// True when `origin` matches.
    pub fn matches(&self, origin: &Origin) -> bool {
        self.step == origin.step && self.role.as_ref().is_none_or(|r| *r == origin.role)
    }
}

impl From<Origin> for OriginPattern {
    fn from(origin: Origin) -> Self {
        OriginPattern {
            step: origin.step,
            role: Some(origin.role),
        }
    }
}

/// Chooses faces of a step's input body by their origins.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FaceRef {
    /// Faces with an origin this matches.
    pub face: OriginPattern,
    /// Keeps only faces that share an edge with a face matching each of these.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub touching: Vec<OriginPattern>,
    /// When present, keeps only the face nearest this point, which must be strictly
    /// nearer than every other.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub near: Option<Point3>,
}

/// Chooses edges of a step's input body by the origins of the faces they separate.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EdgeRef {
    /// Edges between a face matching one of these and a face matching the other. The
    /// order matters only where an operation says so, such as a two-distance chamfer.
    pub between: [OriginPattern; 2],
    /// Keeps only edges with an end vertex on a face matching each of these.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub touching: Vec<OriginPattern>,
    /// When present, keeps only the edge nearest this point, which must be strictly
    /// nearer than every other.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub near: Option<Point3>,
}

/// The edges an edge operation acts on: written `"all"` or as a list of references.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum EdgeSelection {
    /// Every edge of the input body.
    All,
    /// The union of what these references choose.
    Refs(Vec<EdgeRef>),
}

impl Serialize for EdgeSelection {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            EdgeSelection::All => serializer.serialize_str("all"),
            EdgeSelection::Refs(refs) => refs.serialize(serializer),
        }
    }
}

impl<'de> Deserialize<'de> for EdgeSelection {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_any(SelectionVisitor)
    }
}

struct SelectionVisitor;

impl<'de> Visitor<'de> for SelectionVisitor {
    type Value = EdgeSelection;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("\"all\" or a list of edge references")
    }

    fn visit_str<E: de::Error>(self, text: &str) -> Result<EdgeSelection, E> {
        match text {
            "all" => Ok(EdgeSelection::All),
            _ => Err(E::invalid_value(de::Unexpected::Str(text), &self)),
        }
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<EdgeSelection, A::Error> {
        let mut refs = Vec::new();
        while let Some(r) = seq.next_element()? {
            refs.push(r);
        }
        Ok(EdgeSelection::Refs(refs))
    }
}
