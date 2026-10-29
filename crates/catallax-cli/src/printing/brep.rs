//! Recipe construction for the shared body-based print adapter.
use super::Options;
use catallax_features::BuildError;
pub(super) use catallax_field::printing::brep::Output;
use catallax_numbers::{Decider, Refusal};
use catallax_recipe::Recipe;

fn build_error(error: BuildError) -> Refusal {
    match error {
        BuildError::Geometry { step, refusal } => {
            Refusal::new(refusal.code, format!("step {step}: {}", refusal.detail))
        }
        BuildError::Recipe(error) => Refusal::new("cli.invalid-recipe", error.to_string()),
    }
}

pub(super) fn mesh(
    recipe: &Recipe,
    options: &Options,
    certified: bool,
    dec: &mut Decider,
) -> Result<Output, Refusal> {
    if super::native_field(recipe) {
        return Err(Refusal::new(
            "print.brep-field-definition",
            "a field envelope is not a B-rep print definition; request the field route",
        ));
    }
    if certified {
        let built = catallax_features::evaluate_certified(recipe, dec).map_err(build_error)?;
        catallax_field::printing::brep::mesh_certified(&built, options, dec)
    } else {
        let built = catallax_features::evaluate(recipe, dec).map_err(build_error)?;
        catallax_field::printing::brep::mesh(&built, options, dec)
    }
}
