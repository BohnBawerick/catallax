//! Rebuild named entities after a parameter edit, independently of display work.

use std::collections::BTreeMap;

use catallax_features::{BuildError, Built};
use catallax_numbers::Decider;
use catallax_recipe::{EdgeName, Origin, Recipe};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct NameEdit {
    format: String,
    recipe: Recipe,
    edited_recipe: Recipe,
}

#[derive(Default)]
struct Names {
    face_count: usize,
    edges: Vec<EdgeName>,
    origins: BTreeMap<Origin, usize>,
    unnamed_faces: usize,
    unnamed_edges: usize,
}
fn names(built: &Built) -> Names {
    let mut result = Names::default();
    for face in built.model.faces_of_body(built.result) {
        let origins = built.model.face_origins(face);
        if origins.is_empty() {
            result.unnamed_faces += 1;
        }
        result.face_count += 1;
        for origin in origins {
            *result.origins.entry(origin.clone()).or_default() += 1;
        }
    }
    for edge in built.model.edges_of_body(built.result) {
        if let Some(name) = built.model.edge_name(built.result, edge) {
            result.edges.push(name);
        } else {
            result.unnamed_edges += 1;
        }
    }
    result
}
fn build(recipe: &Recipe) -> Result<Built, Value> {
    catallax_features::evaluate(recipe, &mut Decider::default()).map_err(|e| match e {
        BuildError::Geometry { step, refusal } => json!({"code":refusal.code,"detail":refusal.detail,"step":step,"operation":recipe.steps.iter().find(|s| s.id == step).map(|s| s.op.name())}),
        BuildError::Recipe(e) => json!({"code":"cli.invalid-recipe","detail":e.to_string()}),
    })
}

pub(crate) fn check(problem: &Value) -> Value {
    let input: NameEdit = match serde_json::from_value(problem.clone()) {
        Ok(input) => input,
        Err(error) => {
            return json!({"status":"invalid","code":"names.invalid-edit","detail":error.to_string()});
        }
    };
    if input.format != "catallax-name-edit/1" {
        return json!({"status":"invalid","code":"names.invalid-edit","detail":"expected catallax-name-edit/1"});
    }
    for recipe in [&input.recipe, &input.edited_recipe] {
        if let Err(error) = recipe.check() {
            return json!({"status":"invalid","code":"names.invalid-edit","detail":error.to_string()});
        }
    }
    // This adapter checks parameter edits, not deletion or replacement of steps.
    if input.recipe.result != input.edited_recipe.result
        || input.recipe.steps.len() != input.edited_recipe.steps.len()
        || input
            .recipe
            .steps
            .iter()
            .zip(&input.edited_recipe.steps)
            .any(|(a, b)| a.id != b.id || a.op.name() != b.op.name())
    {
        return json!({"status":"invalid","code":"names.not-parameter-edit","detail":"result, step identities and operation kinds must stay fixed"});
    }
    if input.recipe.steps == input.edited_recipe.steps {
        return json!({"status":"invalid","code":"names.no-edit"});
    }
    let before = build(&input.recipe);
    let after = build(&input.edited_recipe);
    let field_before = catallax_field::Plan::compile(&input.recipe);
    let field_after = catallax_field::Plan::compile(&input.edited_recipe);
    let mut waiting = Vec::new();
    for (kernel, phase, error) in [
        ("brep", "before", before.as_ref().err().cloned()),
        ("brep", "after", after.as_ref().err().cloned()),
        (
            "field",
            "before",
            field_before
                .as_ref()
                .err()
                .map(|e| json!({"code":e.code,"detail":e.detail})),
        ),
        (
            "field",
            "after",
            field_after
                .as_ref()
                .err()
                .map(|e| json!({"code":e.code,"detail":e.detail})),
        ),
    ] {
        if let Some(refusal) = error {
            waiting.push(json!({"kernel":kernel,"phase":phase,"refusal":refusal}));
        }
    }
    if !waiting.is_empty() {
        return json!({"status":"waiting","code":"names.build-required","waiting":waiting});
    }
    let before = before.unwrap();
    let after = after.unwrap();
    let old = names(&before);
    let new = names(&after);
    let missing_faces: Vec<_> = old
        .origins
        .iter()
        .filter_map(|(origin, count)| {
            let resolved = after.model.faces_with_origin(after.result, origin).len();
            (resolved == 0).then(|| json!({"name":origin,"before":count,"resolved":resolved}))
        })
        .collect();
    let missing_edges: Vec<_> = old
        .edges
        .iter()
        .filter_map(|name| {
            let resolved = after.model.edges_with_name(after.result, name).len();
            (resolved == 0).then(|| json!({"name":name,"resolved":resolved}))
        })
        .collect();
    let mut missing_field = Vec::new();
    for (phase, plan, origins) in [
        ("before", field_before.unwrap(), &old.origins),
        ("after", field_after.unwrap(), &new.origins),
    ] {
        for origin in origins.keys() {
            if plan.select(origin).is_none() {
                missing_field.push(json!({"phase":phase,"origin":origin}));
            }
        }
    }
    let unnamed = old.unnamed_faces + new.unnamed_faces + old.unnamed_edges + new.unnamed_edges;
    let passed = unnamed == 0
        && missing_faces.is_empty()
        && missing_edges.is_empty()
        && missing_field.is_empty();
    json!({"status":if passed {"passed"} else {"failed"},"code":if passed {"names.survived"} else {"names.lost-reference"},
        "brep":{"coverage":"all result face origins and physical edge names; split pieces and merged origin sets resolve by containment",
            "faces_checked":old.face_count,"edges_checked":old.edges.len(),
            "missing_faces":missing_faces,"missing_edges":missing_edges,
            "unnamed_before":{"faces":old.unnamed_faces,"edges":old.unnamed_edges},"unnamed_after":{"faces":new.unnamed_faces,"edges":new.unnamed_edges}},
        "field":{"coverage":"B-rep result origins resolve in each independently compiled definition registry; not a boundary certificate","missing":missing_field}})
}
