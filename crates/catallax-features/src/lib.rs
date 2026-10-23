//! Recipe-based B-rep features with retained exact definitions.
//!
//! [`evaluate`] builds stock solids, sketches, extrusions, revolutions, elementary
//! fillets, booleans, exact body patterns and mirrors.
//! Every returned face has canonical recipe origins, including faces split or merged
//! by a boolean. The model remains compatible with `catallax-measure` and
//! `catallax-step`; the legacy CLI protocol is not extended.
//!
//! Body-copy placements may be algebraic.
//! [`evaluate_certified`] additionally builds the checked inward elliptic-offset
//! prism family with genuine spline geometry and a complete geometry-bound record.
//! Its retained exact definition remains the target for later offsets and properties.
//! Unsupported legal geometry keeps an explicit coverage refusal; a surface enum
//! variant is not general feature or boolean coverage. `docs/RECIPE.md` section 7.14
//! records the remaining legal-input coverage gaps.
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

mod certified;
mod edge_blend;
mod evaluate;
mod extrude;
/// Exact elementary rolling-ball definitions retained beside their topology.
pub mod fillet;
mod fillet_build;
mod offset_fit;
mod placement;
mod profile;
mod revolve;
mod stock;

pub use certified::{BodyCertificate, CertificateKey, CertifiedBuilt, PrismTarget};
pub use evaluate::{BuildError, Built, evaluate, evaluate_certified};
pub use extrude::extrude;
pub use offset_fit::{EllipticOffsetFit, OffsetSpanBound, fit_elliptic_offset};
pub use placement::{copy_origin, mirror_placement, pattern_placement};
pub use profile::Profile;
pub use revolve::revolve;
pub use stock::evaluate_canonical_stock;
