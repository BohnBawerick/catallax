//! Embedded, immutable STEP sources for recipe, field and comparison requests.

use catallax_features::ImportFiles;
use catallax_numbers::Refusal;
use serde_json::Value;
use sha2::{Digest, Sha256};

pub(crate) fn parse(value: &Value) -> Result<Option<(&Value, ImportFiles<'_>)>, Refusal> {
    if value["format"] != "catallax-import-recipe/1" {
        return Ok(None);
    }
    let bad = || {
        Refusal::new(
            "cli.import-sources",
            "import wrapper needs exactly format, recipe and files; files maps SHA-256 to STEP text",
        )
    };
    let object = value.as_object().ok_or_else(bad)?;
    if object.len() != 3 || !object.contains_key("recipe") {
        return Err(bad());
    }
    let mut files = ImportFiles::new();
    for (hash, text) in value["files"].as_object().ok_or_else(bad)? {
        let text = text.as_str().ok_or_else(bad)?;
        let actual = format!("{:x}", Sha256::digest(text.as_bytes()));
        if *hash != actual {
            return Err(Refusal::new(
                "import.hash-mismatch",
                format!("embedded STEP source expected {hash}, actual {actual}"),
            ));
        }
        files.insert(hash.clone(), text.as_bytes());
    }
    Ok(Some((&value["recipe"], files)))
}
