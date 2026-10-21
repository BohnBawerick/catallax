//! The Catallax recipe language: its data model.
//!
//! A recipe is the feature tree both kernels build a part from. It is a list of steps;
//! each step runs one operation, reads the values of earlier steps by their ids, and
//! names faces and edges by the operations that made them, never by their order in a
//! result. This crate holds the types, their canonical JSON text (`catallax-recipe/3`)
//! and the checks that need no geometry. It evaluates nothing.
//!
//! The contract is `docs/RECIPE.md`.
//!
//! ```
//! use catallax_recipe::Recipe;
//!
//! let text = r#"{
//!   "format": "catallax-recipe/3",
//!   "name": "plate",
//!   "steps": [
//!     { "id": "outline", "op": "sketch",
//!       "plane": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] },
//!       "loops": [[
//!         { "kind": "line", "id": "bottom", "to": ["40", "0"] },
//!         { "kind": "line", "id": "right", "to": ["40", "30"] },
//!         { "kind": "line", "id": "top", "to": ["0", "30"] },
//!         { "kind": "line", "id": "left", "to": ["0", "0"] }
//!       ]] },
//!     { "id": "plate", "op": "extrude", "profile": "outline", "distance": "5.0" },
//!     { "id": "round", "op": "fillet", "body": "plate", "radius": "1/2",
//!       "edges": [{ "between": [
//!         { "step": "plate", "role": "end" },
//!         { "step": "plate", "role": { "side": { "segment": "top" } } }
//!       ] }] }
//!   ],
//!   "result": "round"
//! }"#;
//! let recipe = Recipe::parse(text).unwrap();
//! let canonical = recipe.to_text();
//! assert!(canonical.contains(r#""distance": "5""#));
//! assert!(canonical.contains(r#""radius": "0.5""#));
//! assert_eq!(Recipe::parse(&canonical).unwrap(), recipe);
//! ```

mod check;
mod ids;
mod names;
mod number;
mod operation;
mod recipe;
mod text;

pub use ids::{ElementId, IdError, StepId};
pub use names::{
    EdgeName, EdgeRef, EdgeSelection, FaceRef, HolePart, NameError, Origin, OriginPattern, Role,
    VertexName,
};
pub use number::{Axis, Number, Plane, Point2, Point3, Pose, Vector3};
pub use operation::{
    BooleanKind, FileRef, Grade, GyroidKind, HoleDepth, HoleHead, HolePlacement, HoleTip,
    LatticeCell, LoftSection, LoftShape, Loop, Operation, PathSegment, PatternLayout, Segment,
    SweepOrientation, ThickenSide, Turn, ValueKind,
};
pub use recipe::{FORMAT, Recipe, RecipeError, Step};
