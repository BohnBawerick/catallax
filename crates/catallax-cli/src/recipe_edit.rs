//! Scalar edits rebuild the checked recipe; no internal setter bypasses names.

use catallax_recipe::{Number, Recipe};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ScalarEdit {
    pub(crate) step: String,
    pub(crate) parameter: String,
    pub(crate) value: Number,
}

impl ScalarEdit {
    pub(crate) fn apply(&self, recipe: &Recipe) -> Result<Recipe, String> {
        let mut value: Value =
            serde_json::from_str(&recipe.to_text()).map_err(|e| e.to_string())?;
        let steps = value["steps"].as_array_mut().unwrap();
        let step = steps
            .iter_mut()
            .find(|s| s["id"] == self.step)
            .ok_or_else(|| "edited step does not exist".to_owned())?;
        let parameter = step
            .get_mut(&self.parameter)
            .filter(|v| v.is_string())
            .ok_or_else(|| "edit must name an existing scalar parameter".to_owned())?;
        *parameter = json!(self.value);
        Recipe::parse(&value.to_string()).map_err(|e| e.to_string())
    }
}
