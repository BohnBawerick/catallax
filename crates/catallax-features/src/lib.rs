//! Exact recipe-based features over the B-rep model.
//!
//! [`evaluate`] builds stock solids, sketches, extrusions, revolutions and booleans.
//! Every returned face has canonical recipe origins, including faces split or merged
//! by a boolean. The model remains compatible with `catallax-measure` and
//! `catallax-step`; the legacy CLI protocol is not extended.
//!
//! Exact analytic frames are rational in this implementation. Unsupported algebraic
//! placements, generated radii and elliptic cylinders have explicit refusal codes;
//! `docs/RECIPE.md` section 7.14 lists these unfinished legal-input cases.
//!
//! ```
//! use catallax_features::evaluate;
//! use catallax_numbers::Decider;
//! use catallax_recipe::Recipe;
//!
//! let recipe = Recipe::parse(r#"{
//!   "format":"catallax-recipe/3", "name":"tube",
//!   "steps":[
//!     {"id":"section","op":"sketch",
//!      "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
//!      "loops":[[{"kind":"circle","id":"wall","centre":["0","0"],"radius":"2"}]]},
//!     {"id":"tube","op":"extrude","profile":"section","distance":"3"}
//!   ], "result":"tube"
//! }"#)?;
//! let built = evaluate(&recipe, &mut Decider::default())?;
//! for face in built.model.faces_of_body(built.result) {
//!     println!("{:?}", built.model.face_origins(face));
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod evaluate;
mod extrude;
mod profile;
mod revolve;
mod stock;

pub use evaluate::{BuildError, Built, evaluate};
pub use extrude::extrude;
pub use profile::Profile;
pub use revolve::revolve;
