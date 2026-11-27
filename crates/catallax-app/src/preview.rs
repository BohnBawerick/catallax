//! Opt-in unproved display. Certified display and exports remain separate calls.

use std::collections::BTreeMap;

use catallax_features::{Built, evaluate};
use catallax_mesh::float::{
    DisplayDescription, DisplayOptions, describe_body, describe_numerical_body,
};
use catallax_mesh::{PreviewMesh, PreviewOptions, preview_mesh_body};
use catallax_numbers::Decider;
use catallax_recipe::{Origin, Recipe};

use crate::{Error, Names, NumericalApplicationBody, bounded_build, result_names};

pub use catallax_field::display::PreviewOptions as ConstructionOptions;

/// Limits for the first paint and its caller-scheduled refinement.
#[derive(Clone, Debug)]
pub struct ProgressiveOptions {
    /// Vertex and cell limits for both paints. Its tolerance is overridden.
    pub limits: ConstructionOptions,
    /// Portion of the preview allocation reserved for the 0.1 mm first paint.
    /// The remainder pays for the 0.01 mm refinement, including recompilation.
    pub first_paint_work: u64,
}

/// Unproved display and picking names for one retained body.
#[derive(Clone, Debug)]
pub struct PreviewDisplay {
    /// No certificate, deviation bounds, or export guarantees.
    pub mesh: PreviewMesh,
    /// Face origins and edge/vertex names from the same model revision.
    pub names: Names,
}

/// Unproved GPU preview input and picking names from the same retained build.
#[derive(Clone, Debug, PartialEq)]
pub struct DescribedDisplay {
    /// Shared edge samples, surface parameters and oriented trims.
    pub description: DisplayDescription,
    /// Face origins and edge/vertex names keyed by revision-local handles.
    pub names: Names,
}

/// Describe the retained body for GPU preview without rebuilding, meshing or
/// running proof. The CPU preview and certified paths are separate calls.
pub fn display_description(
    built: &Built,
    options: &DisplayOptions,
) -> Result<DescribedDisplay, Error> {
    Ok(DescribedDisplay {
        description: describe_body(&built.model, built.result, options).map_err(Error::Display)?,
        names: result_names(&built.model, built.result),
    })
}

/// Unproved GPU input for a fitted body and its stable face origins.
#[derive(Clone, Debug, PartialEq)]
pub struct DescribedNumericalDisplay {
    /// Edge and face keys belong to the admitted numerical definition.
    pub description: DisplayDescription<String, String>,
    /// Face origins keyed by the same stable face keys.
    pub face_origins: BTreeMap<String, Vec<Origin>>,
}

/// Describe a retained fitted B-rep, without repeating its admission checks.
pub fn numerical_display_description(
    built: &NumericalApplicationBody,
    options: &DisplayOptions,
) -> Result<DescribedNumericalDisplay, Error> {
    Ok(DescribedNumericalDisplay {
        description: describe_numerical_body(built.body(), options).map_err(Error::Display)?,
        face_origins: built.face_origins().clone(),
    })
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
/// 100M [`crate::DEFAULT_BUILD_WORK`] application limit applies.
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

/// Construction display from immutable definitions, before B-rep evaluation.
/// Registry indices are local to this preview and are never B-rep handles.
#[derive(Clone, Debug)]
pub struct ConstructionDisplay {
    /// Unproved triangles. Incomplete cell covers are refused before delivery.
    pub mesh: catallax_field::display::PreviewMesh,
    /// Origins of displayed faces only, keyed by preview registry index.
    pub origins: std::collections::BTreeMap<usize, catallax_recipe::Origin>,
}

/// Constructs a selectable display without first constructing an exact B-rep.
/// Unsupported field definitions and incomplete displays refuse explicitly.
/// This call shares the active allowance and supplies no export certificate.
pub fn construction_preview(
    recipe: &Recipe,
    options: &ConstructionOptions,
    decider: &mut Decider,
) -> Result<ConstructionDisplay, Error> {
    bounded_build(|| {
        let plan = catallax_field::Plan::compile_checked(recipe, decider)
            .map_err(Error::ConstructionPreview)?;
        let mesh = catallax_field::display::preview_plan(&plan, options, decider)
            .map_err(Error::ConstructionPreview)?;
        if mesh.unresolved_cells != 0 {
            return Err(Error::ConstructionPreview(catallax_numbers::Refusal::new(
                "application.preview-incomplete",
                format!(
                    "{} construction display cells remain unresolved",
                    mesh.unresolved_cells
                ),
            )));
        }
        let mut origins = std::collections::BTreeMap::new();
        for triangle in &mesh.triangles {
            let index = triangle.face.index();
            origins
                .entry(index)
                .or_insert_with(|| plan.faces()[index].origin.clone());
        }
        Ok(ConstructionDisplay { mesh, origins })
    })
}
