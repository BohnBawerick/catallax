//! Opt-in unproved display. Certified display and exports remain separate calls.

use catallax_features::{Built, evaluate};
use catallax_mesh::{PreviewMesh, PreviewOptions, preview_mesh_body};
use catallax_numbers::Decider;
use catallax_recipe::Recipe;

use crate::{Error, Names, bounded_build, result_names};

/// Unproved display and picking names for one retained body.
#[derive(Clone, Debug)]
pub struct PreviewDisplay {
    /// No certificate, deviation bounds, or export guarantees.
    pub mesh: PreviewMesh,
    /// Face origins and edge/vertex names from the same model revision.
    pub names: Names,
}

/// A retained exact build with an explicitly unproved display.
#[derive(Clone, Debug)]
pub struct PreviewApplicationBody {
    /// Retained definitions. Use these for later certified mesh or export requests.
    pub built: Built,
    /// Unproved display and revision-local picking handles.
    pub preview: PreviewDisplay,
}

/// Requests a preview from an existing build without rebuilding or certifying it.
///
/// Charge this call to the caller's preview/display allowance, separately from
/// any later proof request. An active allowance always takes precedence; this
/// call never resets, duplicates, or enlarges it. Without one, the usual
/// [`crate::DEFAULT_BUILD_WORK`] application limit applies.
pub fn preview_body(
    built: &Built,
    options: &PreviewOptions,
    decider: &mut Decider,
) -> Result<PreviewDisplay, Error> {
    bounded_build(|| {
        let faults = built.model.check(built.result);
        if !faults.is_empty() {
            return Err(Error::InvalidBody(faults));
        }
        let mesh = preview_mesh_body(&built.model, built.result, options, decider)
            .map_err(Error::Display)?;
        Ok(PreviewDisplay {
            mesh,
            names: result_names(&built.model, built.result),
        })
    })
}

/// Builds retained exact definitions and an unproved preview, never a certified mesh.
/// Construction and preview share one allowance, including a caller-owned one.
/// Existing [`crate::build`] and [`crate::build_certified`] keep their contracts.
pub fn build_preview(
    recipe: &Recipe,
    options: &PreviewOptions,
    decider: &mut Decider,
) -> Result<PreviewApplicationBody, Error> {
    bounded_build(|| {
        let built = evaluate(recipe, decider).map_err(Error::Build)?;
        let preview = preview_body(&built, options, decider)?;
        Ok(PreviewApplicationBody { built, preview })
    })
}
