//! The recipe format: round trips, canonical text, and what the reader and the
//! structural check reject. Needs no kernel.

use catallax_recipe::{
    EdgeName, ElementId, Number, Operation, Origin, Recipe, RecipeError, Role, StepId, VertexName,
};
use serde_json::Value;

const EVERY_OPERATION: &str = include_str!("fixtures/every-operation.json");

fn every_operation() -> Recipe {
    Recipe::parse(EVERY_OPERATION).expect("the fixture is a valid recipe")
}

#[test]
fn every_operation_round_trips_to_the_same_bytes() {
    let recipe = every_operation();
    assert_eq!(
        recipe.to_text(),
        EVERY_OPERATION,
        "the fixture is in canonical form"
    );
    assert_eq!(Recipe::parse(&recipe.to_text()).unwrap(), recipe);
}

/// Writes `value` compactly with every object's keys in reverse order, every exact
/// number respelled, and every edge and vertex name listed backwards.
fn respell(value: &Value, key: Option<&str>, out: &mut String) {
    match value {
        Value::Object(map) => {
            out.push('{');
            for (i, (k, v)) in map.iter().rev().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                out.push_str(&serde_json::to_string(k).unwrap());
                out.push(':');
                let name = matches!(key, Some("edge" | "vertex"))
                    && matches!(k.as_str(), "between" | "faces");
                if let (true, Value::Array(items)) = (name, v) {
                    let reversed: Vec<Value> = items.iter().rev().cloned().collect();
                    respell(&Value::Array(reversed), Some(k), out);
                } else {
                    respell(v, Some(k), out);
                }
            }
            out.push('}');
        }
        Value::Array(items) => {
            out.push('[');
            for (i, v) in items.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                respell(v, key, out);
            }
            out.push(']');
        }
        Value::String(s) if Number::parse(s).is_ok() => {
            let (sign, body) = s.strip_prefix('-').map_or(("", s.as_str()), |b| ("-", b));
            let spelled = match body.split_once('/') {
                Some((p, q)) => format!("{sign}0{p}/0{q}"),
                None if body.contains('.') => format!("{sign}0{body}0e0"),
                None => format!("{sign}0{body}.0e0"),
            };
            assert_eq!(Number::parse(&spelled).unwrap(), Number::parse(s).unwrap());
            out.push_str(&serde_json::to_string(&spelled).unwrap());
        }
        other => out.push_str(&serde_json::to_string(other).unwrap()),
    }
}

#[test]
fn layout_key_order_and_spelling_do_not_change_the_text() {
    let value: Value = serde_json::from_str(EVERY_OPERATION).unwrap();
    let mut text = String::new();
    respell(&value, None, &mut text);
    assert_ne!(text, EVERY_OPERATION);
    assert!(text.contains("\"010.0e0\""), "numbers were respelled");
    let recipe = Recipe::parse(&text).unwrap();
    assert_eq!(recipe, every_operation());
    assert_eq!(recipe.to_text(), EVERY_OPERATION);
}

#[test]
fn writing_is_deterministic() {
    let a = every_operation();
    let b = a.clone();
    let text = a.to_text();
    for _ in 0..3 {
        assert_eq!(a.to_text(), text);
        assert_eq!(b.to_text(), text);
    }
    let rebuilt = Recipe::new(a.name.clone(), a.steps.clone(), a.result.clone());
    let rebuilt = Recipe {
        description: a.description.clone(),
        ..rebuilt
    };
    assert_eq!(rebuilt.to_text(), text);
}

#[test]
fn numbers_have_one_canonical_spelling() {
    for (text, canonical) in [
        ("10", "10"),
        ("10.0", "10"),
        ("007", "7"),
        ("1.10", "1.1"),
        ("-2.5e-3", "-0.0025"),
        ("1e3", "1000"),
        ("2/4", "0.5"),
        ("-6/8", "-0.75"),
        ("1/3", "1/3"),
        ("2/6", "1/3"),
        ("-0", "0"),
        ("0.000", "0"),
    ] {
        let n = Number::parse(text).unwrap();
        assert_eq!(n.canonical(), canonical, "{text}");
        assert_eq!(Number::parse(&n.canonical()).unwrap(), n);
        assert_eq!(
            serde_json::to_string(&n).unwrap(),
            format!("\"{canonical}\"")
        );
    }
    for bad in ["", "1.", ".5", "1/0", "0x10", "1e", "one"] {
        assert!(Number::parse(bad).is_err(), "{bad:?}");
    }
}

fn origin(step: &str, role: Role) -> Origin {
    Origin {
        step: StepId::new(step).unwrap(),
        role,
    }
}

fn side(step: &str, segment: &str) -> Origin {
    origin(
        step,
        Role::Side {
            segment: ElementId::new(segment).unwrap(),
            path: None,
        },
    )
}

#[test]
fn edge_and_vertex_names_are_in_canonical_order() {
    let (end, top) = (origin("base", Role::End), side("base", "top"));
    assert!(
        end < top,
        "a role without fields sorts first: {} {}",
        end.key(),
        top.key()
    );
    let ab = EdgeName::new(end.clone(), top.clone()).unwrap();
    let ba = EdgeName::new(top.clone(), end.clone()).unwrap();
    assert_eq!(ab, ba);
    assert_eq!(ab.between(), &[end.clone(), top.clone()]);
    assert!(EdgeName::new(end.clone(), end.clone()).is_err());
    let left = side("base", "left");
    let v = VertexName::new(vec![left.clone(), top.clone(), end.clone()]).unwrap();
    assert_eq!(v.faces(), [end.clone(), left.clone(), top.clone()]);
    assert!(VertexName::new(vec![end.clone(), top, end]).is_err());
}

/// A small valid recipe; tests break it one way at a time.
const PLATE: &str = r#"{
  "format": "catallax-recipe/3",
  "name": "plate",
  "steps": [
    { "id": "outline", "op": "sketch",
      "plane": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] },
      "loops": [[
        { "kind": "line", "id": "bottom", "to": ["40", "0"] },
        { "kind": "line", "id": "right", "to": ["40", "30"] },
        { "kind": "line", "id": "top", "to": ["0", "30"] },
        { "kind": "line", "id": "left", "to": ["0", "0"] }
      ]] },
    { "id": "plate", "op": "extrude", "profile": "outline", "distance": "5" },
    { "id": "round", "op": "fillet", "body": "plate", "radius": "1",
      "edges": [{ "between": [
        { "step": "plate", "role": "end" },
        { "step": "plate", "role": { "side": { "segment": "top" } } }
      ] }] }
  ],
  "result": "round"
}"#;

fn broken(from: &str, to: &str) -> Result<Recipe, RecipeError> {
    assert!(PLATE.contains(from), "{from}");
    Recipe::parse(&PLATE.replacen(from, to, 1))
}

fn code(result: Result<Recipe, RecipeError>) -> &'static str {
    match result {
        Ok(_) => "ok",
        Err(e) => {
            assert!(e.to_string().starts_with(e.code()), "{e}");
            e.code()
        }
    }
}

#[test]
fn the_reader_rejects_what_is_not_this_format() {
    assert_eq!(code(Recipe::parse(PLATE)), "ok");
    let cases = [
        // A JSON number where an exact number belongs.
        (r#""distance": "5""#, r#""distance": 5"#, "recipe.parse"),
        // Unknown keys, in a step and at the top.
        (
            r#""distance": "5""#,
            r#""distance": "5", "draft": "1""#,
            "recipe.parse",
        ),
        (
            r#""name": "plate","#,
            r#""name": "plate", "units": "mm","#,
            "recipe.parse",
        ),
        // Unknown operation, segment kind and role.
        (r#""op": "extrude""#, r#""op": "boss""#, "recipe.parse"),
        (
            r#""kind": "line", "id": "top""#,
            r#""kind": "spline", "id": "top""#,
            "recipe.parse",
        ),
        (r#""role": "end""#, r#""role": "lid""#, "recipe.parse"),
        // Ids outside the grammar.
        (r#""id": "plate""#, r#""id": "Plate""#, "recipe.parse"),
        (r#""id": "top""#, r#""id": "top edge""#, "recipe.parse"),
        // Another format.
        ("catallax-recipe/3", "catallax-recipe/5", "recipe.format"),
        ("catallax-recipe/3", "catallax-recipe/2", "recipe.format"),
        ("catallax-recipe/3", "catallax-recipe/1", "recipe.format"),
        (r#""format": "catallax-recipe/3","#, "", "recipe.format"),
    ];
    for (from, to, want) in cases {
        assert_eq!(code(broken(from, to)), want, "{from} -> {to}");
    }
    assert_eq!(code(Recipe::parse("{")), "recipe.parse");
}

#[test]
fn malformed_names_are_rejected() {
    let same = r#"{ "step": "plate", "role": { "blend": { "edge": { "between": [
        { "step": "plate", "role": "end" }, { "step": "plate", "role": "end" } ] } } } }"#;
    let text = PLATE.replacen(r#"{ "step": "plate", "role": "end" }"#, same, 1);
    assert_eq!(code(Recipe::parse(&text)), "recipe.parse");
    let short = r#"{ "step": "plate", "role": { "corner": { "vertex": { "faces": [] } } } }"#;
    let text = PLATE.replacen(r#"{ "step": "plate", "role": "end" }"#, short, 1);
    assert_eq!(code(Recipe::parse(&text)), "recipe.parse");
}

#[test]
fn the_check_rejects_inconsistent_recipes() {
    let cases = [
        (
            r#""id": "round""#,
            r#""id": "plate""#,
            "recipe.duplicate-step",
        ),
        (
            r#""profile": "outline""#,
            r#""profile": "sketch""#,
            "recipe.unknown-input",
        ),
        (
            r#""body": "plate""#,
            r#""body": "outline""#,
            "recipe.wrong-input",
        ),
        (
            r#""profile": "outline""#,
            r#""profile": "plate""#,
            "recipe.unknown-input",
        ),
        (
            r#""id": "right""#,
            r#""id": "bottom""#,
            "recipe.duplicate-element",
        ),
        (
            r#""segment": "top""#,
            r#""segment": "lid""#,
            "recipe.unknown-element",
        ),
        (
            r#""step": "plate", "role": "end""#,
            r#""step": "round", "role": "end""#,
            "recipe.unknown-origin",
        ),
        (
            r#""step": "plate", "role": "end""#,
            r#""step": "outline""#,
            "recipe.no-faces",
        ),
        (
            r#""step": "plate", "role": "end""#,
            r#""step": "plate", "role": "beams""#,
            "recipe.wrong-role",
        ),
        (
            r#""result": "round""#,
            r#""result": "outline""#,
            "recipe.result",
        ),
        (
            r#""result": "round""#,
            r#""result": "missing""#,
            "recipe.result",
        ),
    ];
    for (from, to, want) in cases {
        assert_eq!(code(broken(from, to)), want, "{from} -> {to}");
    }
}

/// Replaces the fillet step of [`PLATE`] with `step`.
fn with_last_step(step: &str) -> Result<Recipe, RecipeError> {
    let start = PLATE.find(r#"{ "id": "round""#).unwrap();
    let end = PLATE.find("\n  ],\n  \"result\"").unwrap();
    Recipe::parse(&format!("{}{step}{}", &PLATE[..start], &PLATE[end..]))
}

#[test]
fn the_check_knows_which_operation_makes_which_names() {
    let ok =
        r#"{ "id": "round", "op": "chamfer", "body": "plate", "distance": "1", "edges": "all" }"#;
    assert_eq!(code(with_last_step(ok)), "ok");
    let cases = [
        // A two-distance chamfer must say which face is first.
        (
            r#"{ "id": "round", "op": "chamfer", "body": "plate", "distance": "1", "other_distance": "2", "edges": "all" }"#,
            "recipe.asymmetric-all",
        ),
        // A boolean keeps its inputs' names and makes none.
        (
            r#"{ "id": "both", "op": "boolean", "kind": "union", "target": "plate", "tool": "plate" },
               { "id": "round", "op": "fillet", "body": "both", "radius": "1",
                 "edges": [{ "between": [{ "step": "both" }, { "step": "plate" }] }] }"#,
            "recipe.no-faces",
        ),
        // A blind hole with a flat tip has a floor and no cone.
        (
            r#"{ "id": "round", "op": "hole", "target": "plate", "diameter": "3", "depth": { "blind": "2" }, "tip": "flat",
                 "at": [{ "id": "h", "point": ["5", "5", "5"], "direction": ["0", "0", "-1"] }] },
               { "id": "r2", "op": "fillet", "body": "round", "radius": "1",
                 "edges": [{ "between": [{ "step": "round", "role": { "hole": { "at": "h", "part": "tip" } } }, { "step": "plate" }] }] }"#,
            "recipe.wrong-role",
        ),
        (
            r#"{ "id": "round", "op": "hole", "target": "plate", "diameter": "3", "depth": "through", "tip": "flat",
                 "at": [{ "id": "h", "point": ["5", "5", "5"], "direction": ["0", "0", "-1"] }] },
               { "id": "r2", "op": "fillet", "body": "round", "radius": "1",
                 "edges": [{ "between": [{ "step": "round", "role": { "hole": { "at": "g", "part": "wall" } } }, { "step": "plate" }] }] }"#,
            "recipe.unknown-element",
        ),
        // A pattern of three makes copies 1 and 2; a mirror makes copy 1.
        (
            r#"{ "id": "round", "op": "pattern", "body": "plate", "layout": { "linear": { "vector": ["50", "0", "0"], "count": 3 } } },
               { "id": "r2", "op": "fillet", "body": "round", "radius": "1",
                 "edges": [{ "between": [{ "step": "round", "role": { "copy": { "instance": 3, "of": { "step": "plate", "role": "end" } } } }, { "step": "plate" }] }] }"#,
            "recipe.no-such-copy",
        ),
        (
            r#"{ "id": "round", "op": "mirror", "body": "plate", "keep_original": true, "plane": { "point": ["0", "0", "0"], "normal": ["1", "0", "0"] } },
               { "id": "r2", "op": "fillet", "body": "round", "radius": "1",
                 "edges": [{ "between": [{ "step": "round", "role": { "copy": { "instance": 0, "of": { "step": "plate", "role": "end" } } } }, { "step": "plate" }] }] }"#,
            "recipe.no-such-copy",
        ),
        // A name inside a role names a face of the operation's input, so an earlier step.
        (
            r#"{ "id": "round", "op": "fillet", "body": "plate", "radius": "1", "edges": "all" },
               { "id": "r2", "op": "fillet", "body": "round", "radius": "1",
                 "edges": [{ "between": [{ "step": "round", "role": { "blend": { "edge": { "between": [
                   { "step": "round", "role": "end" }, { "step": "plate", "role": "end" } ] } } } }, { "step": "plate" }] }] }"#,
            "recipe.unknown-origin",
        ),
        // A loft needs two sections, each naming one start segment per loop.
        (
            r#"{ "id": "round", "op": "loft", "shape": "ruled", "sections": [{ "profile": "outline", "start": ["top"] }] }"#,
            "recipe.loft-sections",
        ),
        (
            r#"{ "id": "round", "op": "loft", "shape": "ruled",
                 "sections": [{ "profile": "outline", "start": ["top", "left"] }, { "profile": "outline", "start": ["top"] }] }"#,
            "recipe.loft-starts",
        ),
        (
            r#"{ "id": "round", "op": "import_mesh", "file": { "name": "scan.stl", "sha256": "ABC" }, "scale": "1",
                 "pose": { "rotation": ["1", "0", "0", "0"], "translation": ["0", "0", "0"] } }"#,
            "recipe.file-hash",
        ),
    ];
    for (step, want) in cases {
        assert_eq!(code(with_last_step(step)), want, "{step}");
    }
}

#[test]
fn steps_expose_their_inputs_and_elements() {
    let recipe = every_operation();
    let id = |s: &str| StepId::new(s).unwrap();
    let blob = recipe.step(&id("blob")).unwrap();
    let inputs: Vec<&str> = blob.op.inputs().iter().map(|(s, _)| s.as_str()).collect();
    assert_eq!(inputs, ["mirrored", "sponge"]);
    let outline = recipe.step(&id("outline")).unwrap();
    let elements: Vec<&str> = outline.op.elements().iter().map(|e| e.as_str()).collect();
    assert_eq!(
        elements,
        ["bottom", "right", "corner", "top", "left", "bore"]
    );
    assert!(matches!(outline.op, Operation::Sketch { .. }));
}

const BLOCK_AND_LATTICE: &str = r#"{"format":"catallax-recipe/3","name":"terminal lattice","steps":[{"id":"region","op":"block","size":["4","4","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"radius":{"constant":"0.25"}}],"result":"beams"}"#;

#[test]
fn a_lattice_that_nothing_reads_is_terminal() {
    let recipe = Recipe::parse(BLOCK_AND_LATTICE).unwrap();
    assert_eq!(recipe.check_lattice_terminal(), Ok(()));
}

#[test]
fn a_cad_step_after_a_lattice_is_not_terminal() {
    let text = BLOCK_AND_LATTICE
        .replace(r#""result":"beams""#, r#""result":"after""#)
        .replace(
            "}],\"result",
            r#"},{"id":"after","op":"offset","body":"beams","distance":"0.1"}],"result"#,
        );
    let recipe = Recipe::parse(&text).unwrap();
    let error = recipe.check_lattice_terminal().unwrap_err();
    assert_eq!(error.code(), "recipe.lattice-not-terminal");
    assert_eq!(
        error,
        RecipeError::LatticeNotTerminal {
            step: StepId::new("after").unwrap(),
            lattice: StepId::new("beams").unwrap(),
        }
    );
}

#[test]
fn a_union_of_a_lattice_is_not_terminal() {
    let error = every_operation().check_lattice_terminal().unwrap_err();
    assert!(matches!(
        error,
        RecipeError::LatticeNotTerminal { ref lattice, .. } if lattice.as_str() == "truss"
    ));
}
