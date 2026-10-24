//! Opt-in bounded geometry; ordinary exact application builds are unchanged.

use catallax_features::{CertifiedBuilt, evaluate_certified};
use catallax_mesh::{Mesh, MeshOptions, mesh_body};
use catallax_numbers::{Answer, Answered, Decider, Refusal};
use catallax_recipe::Recipe;

use crate::{Error, Names, result_names};

/// Geometry-bound certified result, picking names and optional display. The
/// fields stay private so a mesh or name table cannot accidentally be paired
/// with another revision's certificate. Numeric handles belong to this result.
#[derive(Clone, Debug)]
pub struct CertifiedApplicationBody {
    built: CertifiedBuilt,
    display: Option<Result<Answer<Mesh>, Refusal>>,
    names: Names,
}

impl CertifiedApplicationBody {
    /// Immutable model, retained original target and complete geometry record.
    pub fn built(&self) -> &CertifiedBuilt {
        &self.built
    }
    /// Requested model-relative display; absent means unrequested. A display
    /// refusal leaves the complete certified body and its names available.
    pub fn display(&self) -> Option<&Result<Answer<Mesh>, Refusal>> {
        self.display.as_ref()
    }
    /// Stable origins behind this revision's picking handles.
    pub fn names(&self) -> &Names {
        &self.names
    }

    /// Requested properties of the original exact target, not the fitted solid
    /// or display. A failed measurement does not discard the valid body.
    pub fn mass_properties(
        &self,
        accuracy: &catallax_measure::Accuracy,
        dec: &mut Decider,
    ) -> Answered<catallax_measure::MassProperties> {
        let target = self.built.result_certificate().target();
        if target.is_empty() {
            let built = self.built.built();
            catallax_measure::mass_properties(&built.model, built.result, accuracy, dec)
        } else {
            catallax_measure::elliptic_offset_prism_mass_properties(
                target.definition(),
                target.axial(),
                accuracy,
                dec,
            )
        }
    }

    /// Volume, area and centroid without computing second moments.
    pub fn volume_area_centroid(
        &self,
        accuracy: &catallax_measure::Accuracy,
        dec: &mut Decider,
    ) -> Answered<catallax_measure::VolumeAreaCentroid> {
        let target = self.built.result_certificate().target();
        if target.is_empty() {
            let built = self.built.built();
            catallax_measure::volume_area_centroid(&built.model, built.result, accuracy, dec)
        } else {
            catallax_measure::elliptic_offset_prism_volume_area_centroid(
                target.definition(),
                target.axial(),
                accuracy,
                dec,
            )
        }
    }
}

/// Builds the complete certified feature result and its names. A display is
/// computed only when requested; it never replaces the original target or
/// changes modelling settings. All exposed geometry has already passed the
/// producer's complete correspondence, local-limit and material-topology checks.
pub fn build_certified(
    recipe: &Recipe,
    display: Option<&MeshOptions>,
    dec: &mut Decider,
) -> Result<CertifiedApplicationBody, Error> {
    crate::bounded_build(|| {
        let built = evaluate_certified(recipe, dec).map_err(Error::Build)?;
        let represented = built.built();
        let display =
            display.map(|options| mesh_body(&represented.model, represented.result, options, dec));
        let names = result_names(&represented.model, represented.result);
        Ok(CertifiedApplicationBody {
            built,
            display,
            names,
        })
    })
}
