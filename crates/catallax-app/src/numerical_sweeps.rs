//! Numerical sweep recipes stay separate from exact application checkpoints.

use catallax_features::general_sweep::{self, RecipeSwept, SweepChecks};
use catallax_measure::{Accuracy, VolumeAreaCentroid};
use catallax_mesh::{Mesh, MeshOptions};
use catallax_numbers::{Decider, MarkLog, Rational, Refusal};
use catallax_recipe::{Origin, Recipe, StepId};
use catallax_topology::{
    BodyId, Model,
    numerical::{self, NumericalBody, NumericalOptions},
};
use std::collections::BTreeMap;

/// Numerical application geometry and its privately retained rational fit.
/// No exact build or evaluation checkpoint is exposed for the authored recipe.
pub struct SweepBody {
    checked: RecipeSwept,
    fitted: Model,
    fitted_body: BodyId,
}

/// Display geometry derived from the rational fit, not an exact original-body mesh.
pub struct SweepDisplay {
    /// Mesh of the fitted B-rep.
    pub fitted_mesh: Mesh,
    /// Display marks retained from the mesh producer.
    pub marks: MarkLog,
    /// The modelling tolerance, not a complete-body mesh error bound.
    pub source_tolerance_mm: Rational,
}

/// Enclosures about the rational fit only. The sweep's error to the authored
/// model remains numerical and is not included in these integration enclosures.
pub struct SweepProperties {
    /// Integration enclosures about the fitted B-rep.
    pub fitted: VolumeAreaCentroid,
    /// Integration marks of the fitted representation.
    pub marks: MarkLog,
    /// The modelling tolerance, not a bound to add to these enclosures.
    pub source_tolerance_mm: Rational,
}

impl SweepBody {
    /// Admitted numerical body.
    pub fn body(&self) -> &NumericalBody {
        &self.checked.body
    }
    /// Authored recipe needed to save and rebuild this result.
    pub fn recipe(&self) -> &Recipe {
        &self.checked.recipe
    }
    /// Pointwise fitting diagnostics per sweep operation.
    pub fn checks(&self) -> &BTreeMap<StepId, SweepChecks> {
        &self.checked.checks
    }
    /// Source operation origins, independent of numerical chart subdivisions.
    pub fn face_origins(&self) -> &BTreeMap<String, Vec<Origin>> {
        &self.checked.face_origins
    }

    /// A complete parameter rebuild; the old result is unchanged on refusal.
    pub fn edit(&self, recipe: &Recipe, dec: &mut Decider) -> Result<Self, super::Error> {
        build(recipe, self.body().options().clone(), dec)
    }

    /// Unproved preview of the retained fit, on the caller's preview allowance.
    pub fn preview(
        &self,
        options: &catallax_mesh::PreviewOptions,
        dec: &mut Decider,
    ) -> Result<catallax_mesh::PreviewMesh, Refusal> {
        catallax_mesh::preview_mesh_body(&self.fitted, self.fitted_body, options, dec)
    }

    /// Certified display of the fit. This does not certify its modelling error.
    pub fn display(
        &self,
        options: &MeshOptions,
        dec: &mut Decider,
    ) -> Result<SweepDisplay, Refusal> {
        let mesh = catallax_mesh::mesh_body(&self.fitted, self.fitted_body, options, dec)?;
        Ok(SweepDisplay {
            fitted_mesh: mesh.value,
            marks: mesh.marks,
            source_tolerance_mm: self.body().options().tolerance_mm.clone(),
        })
    }

    /// Measures the fitted approximation without returning an exact answer
    /// about the original helix or sweep definition.
    pub fn properties(
        &self,
        accuracy: &Accuracy,
        dec: &mut Decider,
    ) -> Result<SweepProperties, Refusal> {
        let answer =
            catallax_measure::volume_area_centroid(&self.fitted, self.fitted_body, accuracy, dec)?;
        Ok(SweepProperties {
            fitted: answer.value,
            marks: answer.marks,
            source_tolerance_mm: self.body().options().tolerance_mm.clone(),
        })
    }
}

/// Builds a recipe with helical sweeps under the caller's stated tolerance.
/// Reopening a saved ordinary recipe calls this same builder, not a cached fit.
pub fn build(
    recipe: &Recipe,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<SweepBody, super::Error> {
    super::bounded_build(|| {
        let candidate =
            general_sweep::prepare_recipe(recipe, options, dec).map_err(super::Error::Build)?;
        let (checked, (fitted, fitted_body)) = candidate
            .admit(
                |model, id, options, dec| {
                    let (definition, trims) =
                        catallax_features::numerical_embedding::definition_from_exact_model(
                            &model,
                            id,
                            recipe.name.clone(),
                            &options.tolerance_mm,
                            dec,
                        )?;
                    let body = numerical::build_with_native_trims(
                        definition,
                        trims,
                        options.clone(),
                        dec,
                    )?;
                    Ok((body, (model, id)))
                },
                dec,
            )
            .map_err(|refusal| {
                super::Error::Build(catallax_features::BuildError::Geometry {
                    step: recipe.result.clone(),
                    refusal,
                })
            })?;
        Ok(SweepBody {
            checked,
            fitted,
            fitted_body,
        })
    })
}
