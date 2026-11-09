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
    result.edges.sort_by(|a, b| a.between().cmp(b.between()));
    result
}
fn build(recipe: &Recipe) -> Result<Built, Value> {
    catallax_features::evaluate(recipe, &mut Decider::default()).map_err(|e| match e {
        BuildError::Geometry { step, refusal } => json!({"code":refusal.code,"detail":refusal.detail,"step":step,"operation":recipe.steps.iter().find(|s| s.id == step).map(|s| s.op.name())}),
        BuildError::Recipe(e) => json!({"code":"cli.invalid-recipe","detail":e.to_string()}),
    })
}

fn summary(rows: &[Value]) -> Value {
    let count = |status: &str| rows.iter().filter(|row| row["status"] == status).count();
    json!({"survived":count("survived"),"renamed":count("renamed"),"lost":count("lost")})
}

fn brep_check(after: &Built, old: &Names, new: &Names) -> Value {
    let faces: Vec<_> = old
        .origins
        .iter()
        .map(|(origin, count)| {
            let resolved = after.model.faces_with_origin(after.result, origin).len();
            json!({"name":origin,"before":count,"resolved":resolved,
            "status":if resolved == 0 {"lost"} else {"survived"}})
        })
        .collect();
    let edges: Vec<_> = old.edges.chunk_by(|a, b| a == b).map(|group| {
        let name = &group[0];
        let resolved = after.model.edges_with_name(after.result, name).len();
        let published = new.edges.iter().filter(|candidate| *candidate == name).count();
        json!({"name":name,"before":group.len(),"resolved":resolved,"published_after":published,
            "status":if resolved == 0 {"lost"} else if published == 0 {"renamed"} else {"survived"}})
    }).collect();
    let missing_faces: Vec<_> = faces.iter().filter(|row| row["status"] == "lost").collect();
    let missing_edges: Vec<_> = edges.iter().filter(|row| row["status"] == "lost").collect();
    let renamed_edges: Vec<_> = edges
        .iter()
        .filter(|row| row["status"] == "renamed")
        .collect();
    let unnamed = old.unnamed_faces + new.unnamed_faces + old.unnamed_edges + new.unnamed_edges;
    let passed = unnamed == 0
        && missing_faces.is_empty()
        && missing_edges.is_empty()
        && renamed_edges.is_empty();
    json!({"status":if passed {"passed"} else {"failed"},
        "coverage":"all result face origins and physical edge construction names; splits and merged origin sets resolve by containment, never coordinates",
        "faces_checked":old.face_count,"edges_checked":old.edges.len(),
        "summary":{"faces":summary(&faces),"edges":summary(&edges)},
        "faces":faces,"edges":edges,
        "missing_faces":missing_faces,"missing_edges":missing_edges,"renamed_edges":renamed_edges,
        "unnamed_before":{"faces":old.unnamed_faces,"edges":old.unnamed_edges},
        "unnamed_after":{"faces":new.unnamed_faces,"edges":new.unnamed_edges}})
}

fn field_check(
    before: &catallax_field::Plan,
    after: &catallax_field::Plan,
    old: &Names,
    new: &Names,
) -> Value {
    let mut rows = Vec::new();
    for (phase, plan, origins) in [
        ("before", before, &old.origins),
        ("after", after, &new.origins),
    ] {
        for origin in origins.keys() {
            let resolved = plan
                .faces()
                .iter()
                .enumerate()
                .filter(|(index, face)| {
                    &face.origin == origin && plan.result_face(*index).is_some()
                })
                .count();
            rows.push(json!({"phase":phase,"name":origin,"resolved":resolved,
                "status":if resolved == 0 {"lost"} else {"survived"}}));
        }
    }
    let missing: Vec<_> = rows.iter().filter(|row| row["status"] == "lost").collect();
    json!({"status":if missing.is_empty() {"passed"} else {"failed"},
        "coverage":"B-rep result origins resolve to current independent field registry candidates; no field boundary or edge-incidence certificate",
        "summary":summary(&rows),"names":rows,"missing":missing})
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
    let (brep, field) = match (&before, &after) {
        (Ok(before), Ok(after)) => {
            let old = names(before);
            let new = names(after);
            let brep = brep_check(after, &old, &new);
            let field = match (&field_before, &field_after) {
                (Ok(before), Ok(after)) => field_check(before, after, &old, &new),
                _ => json!({"status":"waiting","code":"names.build-required"}),
            };
            (brep, field)
        }
        _ => (
            json!({"status":"waiting","code":"names.build-required"}),
            json!({"status":"waiting","code":"names.brep-result-required"}),
        ),
    };
    let status = brep["status"].clone();
    let code = match status.as_str() {
        Some("passed") => "names.survived",
        Some("failed") => "names.lost-reference",
        _ => "names.build-required",
    };
    json!({"status":status,"code":code,"gate":"brep","brep":brep,"field":field,"waiting":waiting})
}
