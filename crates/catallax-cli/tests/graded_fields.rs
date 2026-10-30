//! V7d-K public print and stand-in requests, fixed by the independent V7d designer:
//! the six graded states, the wrapper and schema refusals, and the explicit exchange.
use catallax_cli::verification::{Request, answer};
use catallax_numbers::Rational;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    sync::atomic::{AtomicUsize, Ordering},
};

static NEXT: AtomicUsize = AtomicUsize::new(0);
const KERNEL_WORK: u64 = 99_999_999;

fn d0() -> Value {
    json!({"format":"catallax-recipe/4","name":"V7d-D0","steps":[
        {"id":"region","op":"block","size":["6","6","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-3","-3","20"]}},
        {"id":"field","op":"lattice","region":"region","cell":"cubic","size":"100",
         "frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
         "radius":{"linear":{"from":["0","0","20"],"to":["0","0","30"],"start":"1","end":"2"}}}],"result":"field"})
}
fn g0() -> Value {
    json!({"format":"catallax-recipe/4","name":"V7d-G0","steps":[
        {"id":"region","op":"block","size":["2","2","2"],"pose":{"rotation":["1","0","0","0"],"translation":["-1","-1","-1"]}},
        {"id":"field","op":"gyroid","region":"region","period":"10",
         "level":{"linear":{"from":["0","0","-1"],"to":["0","0","1"],"start":"-1/10","end":"1/10"}},
         "kind":"network","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"field"})
}
/// Designer section 4: the full print problem.
fn print(recipe: Value, tolerance: &str) -> Value {
    json!({"format":"catallax-print/1","recipe":recipe,"route":"field","tolerance":tolerance,
        "max_cells":2000000,"max_vertices":1000000})
}
/// Designer section 4: the full stand-in problem.
fn stand_in(recipe: Value) -> Value {
    json!({"format":"catallax-stand-in/1","recipe":recipe,"accuracy":{
        "volume":{"relative":"0","absolute":"1/100"},
        "centroid":{"relative":"0","absolute":"1/1000"},
        "inertia":{"relative":"0","absolute":"1/10"}},"max_cells":2000000})
}
fn request(source: &str, problem: Value) -> Request {
    Request {
        format: "catallax-request/1".into(),
        source: source.into(),
        problem,
        out: std::env::temp_dir().join(format!(
            "catallax-v7d-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        )),
        max_work: KERNEL_WORK,
    }
}
fn refusal(source: &str, problem: Value) -> (String, Request) {
    let r = request(source, problem);
    let response = answer(&r).unwrap();
    let code = response["manifest"]["refusal"]["code"]
        .as_str()
        .unwrap_or_else(|| panic!("answered: {response}"))
        .to_string();
    (code, r)
}
fn q(value: &Value) -> Rational {
    Rational::parse(value.as_str().unwrap()).unwrap()
}

/// Designer section 7: wrapper schema and recipe parse refusals happen before any work.
#[test]
fn wrappers_and_recipes_refuse_unknown_controls_and_functions() {
    let with = |problem: Value, field: &str, value: Value| {
        let mut problem = problem;
        problem[field] = value;
        problem
    };
    let detail = with(print(g0(), "1/10"), "minimum_detail", json!("1/10"));
    let work = with(print(g0(), "1/10"), "max_work", json!(1));
    let density = with(stand_in(g0()), "density", json!("1"));
    for (source, problem, expected) in [
        ("print", detail, "cli.invalid-print"),
        ("print", work, "cli.invalid-print"),
        ("recipe", density, "cli.invalid-stand-in"),
    ] {
        let (code, r) = refusal(source, problem);
        assert_eq!(code, expected);
        assert!(!r.out.join("print.stl").exists() && !r.out.join("answer.step").exists());
        fs::remove_dir_all(r.out).unwrap();
    }
    for (field, value) in [
        ("tolerance", json!("0")),
        ("max_cells", json!(0)),
        ("max_vertices", json!(3)),
    ] {
        let (code, r) = refusal("print", with(print(d0(), "1/10"), field, value));
        assert_eq!(code, "print.invalid-options", "{field}");
        fs::remove_dir_all(r.out).unwrap();
    }
    let radial = json!({"axis_distance":{"origin":["0","0","0"],"direction":["0","0","1"],"span":"1","near":"1/10","far":"1/5"}});
    let mut cases = Vec::new();
    for (base, key) in [(d0(), "radius"), (g0(), "level")] {
        let mut recipe = base;
        recipe["steps"][1][key] = radial.clone();
        cases.push(recipe);
    }
    let mut expression = g0();
    expression["steps"][1]["level"] = json!({"expression":"sin(z)"});
    cases.push(expression);
    let mut nonfinite = g0();
    nonfinite["steps"][1]["level"]["linear"]["end"] = json!("NaN");
    cases.push(nonfinite);
    for recipe in cases {
        let (code, r) = refusal("print", print(recipe, "1/10"));
        assert_eq!(code, "cli.invalid-recipe");
        fs::remove_dir_all(r.out).unwrap();
    }
}

/// Work-zero-L and Work-zero-G: an unmodified valid print with an outer allowance of
/// zero refuses at the first checkpoint with the outer allowance's code, with no geometry
/// and no file (designer ruling on the zero-work rows).
#[test]
fn a_zero_outer_allowance_refuses_before_any_geometry() {
    for recipe in [d0(), g0()] {
        let mut r = request("print", print(recipe, "1/10"));
        r.max_work = 0;
        let response = answer(&r).unwrap();
        let manifest = &response["manifest"];
        assert_eq!(manifest["status"], "refused");
        assert_eq!(
            manifest["refusal"]["code"], "numbers.operation-budget",
            "{response}"
        );
        let printing = &manifest["printing"];
        assert_eq!(printing["route"], "field");
        assert!(printing["file"].is_null() && printing["deviation_bound_mm"].is_null());
        assert!(printing.get("sha256").is_none_or(Value::is_null));
        assert!(manifest.get("files").is_none_or(Value::is_null));
        assert_eq!(response["work"]["exhausted"], true);
        assert_eq!(response["work"]["used"], 0);
        for name in ["print.stl", "display.obj", "answer.step"] {
            assert!(!r.out.join(name).exists(), "{name}");
        }
        fs::remove_dir_all(r.out).unwrap();
    }
}

/// Designer sections 4 and 8: graded prints answer on the field route with a labelled
/// file whose metadata and hash match the published bytes, within the kernel allowance.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d prints run in release mode")]
fn graded_prints_publish_labelled_field_files() {
    for (recipe, tolerance) in [(d0(), "1/10"), (d0(), "1/2"), (g0(), "1/10"), (g0(), "1/2")] {
        let r = request("print", print(recipe, tolerance));
        let response = answer(&r).unwrap();
        let manifest = &response["manifest"];
        assert_eq!(manifest["status"], "answered", "{response}");
        let printing = &manifest["printing"];
        assert_eq!(printing["route"], "field");
        assert_eq!(printing["metadata"], "catallax-print-stl/1");
        let bound = q(&printing["deviation_bound_mm"]);
        assert!(bound > Rational::zero() && bound <= Rational::parse(tolerance).unwrap());
        assert!(printing["cells"].as_u64().unwrap() <= 2_000_000);
        assert!(printing["vertices"].as_u64().unwrap() <= 1_000_000);
        let bytes = fs::read(r.out.join("print.stl")).unwrap();
        assert_eq!(&bytes[..8], b"CXPRINT1");
        assert_eq!(
            printing["sha256"].as_str().unwrap(),
            format!("{:x}", Sha256::digest(&bytes))
        );
        let work = &response["work"];
        assert_eq!(work["exhausted"], false);
        assert!(work["used"].as_u64().unwrap() < KERNEL_WORK);
        fs::remove_dir_all(r.out).unwrap();
    }
}

/// Designer section 8: an explicit stand-in of each family exports a massless surface
/// envelope and unit-density proxies whose exchanged bounds lie inside the field's.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d stand-ins run in release mode")]
fn graded_stand_ins_exchange_inside_the_field_enclosures() {
    for recipe in [d0(), g0()] {
        let r = request("recipe", stand_in(recipe));
        let response = answer(&r).unwrap();
        let manifest = &response["manifest"];
        assert_eq!(manifest["status"], "answered", "{response}");
        let stand_in = &manifest["stand_in"];
        assert_eq!(
            stand_in["representation"],
            "massless-surface-envelope-plus-unit-density-solids"
        );
        assert_eq!(stand_in["density_metadata"], false);
        assert_eq!(stand_in["within_field_enclosures"], true);
        assert_eq!(manifest["step"]["surface_envelopes"], 1);
        let field = &stand_in["field_properties"];
        let exchanged = &stand_in["exported_property_bounds"];
        let inside = |outer: &Value, inner: &Value| {
            q(&outer["lo"]) <= q(&inner["lo"]) && q(&inner["hi"]) <= q(&outer["hi"])
        };
        assert!(inside(&field["volume"], &exchanged["volume"]));
        for i in 0..3 {
            assert!(inside(&field["centroid"][i], &exchanged["centroid"][i]));
            for j in 0..3 {
                assert!(inside(&field["inertia"][i][j], &exchanged["inertia"][i][j]));
            }
        }
        let step = fs::read_to_string(r.out.join("answer.step")).unwrap();
        assert_eq!(step.matches("SHELL_BASED_SURFACE_MODEL").count(), 1);
        assert_eq!(
            step.matches("MANIFOLD_SOLID_BREP").count(),
            stand_in["mass_bodies"].as_u64().unwrap() as usize
        );
        assert!(response["work"]["used"].as_u64().unwrap() < KERNEL_WORK);
        fs::remove_dir_all(r.out).unwrap();
    }
}

/// Negative-envelope and plain exchange: a subtracted field has no outer envelope, and
/// a field recipe has no plain STEP boundary. Neither writes a file.
#[test]
fn field_envelopes_never_stand_in_for_a_boundary() {
    let mut negative = d0();
    negative["steps"].as_array_mut().unwrap().push(
        json!({"id":"negative","op":"boolean","kind":"subtract","target":"region","tool":"field"}),
    );
    negative["result"] = json!("negative");
    let (code, r) = refusal("recipe", stand_in(negative));
    assert_eq!(code, "field.stand-in-envelope-coverage");
    assert!(!r.out.join("answer.step").exists());
    fs::remove_dir_all(r.out).unwrap();
    for recipe in [d0(), g0()] {
        let (code, r) = refusal("recipe", recipe);
        assert_eq!(code, "field.stand-in-export");
        assert!(!r.out.join("answer.step").exists());
        fs::remove_dir_all(r.out).unwrap();
    }
}

/// Designer section 8: identical pinned requests in fresh request threads publish the
/// same bytes, metadata and deterministic work.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7d prints run in release mode")]
fn identical_requests_publish_identical_bytes_and_work() {
    let run = || {
        let r = request("print", print(g0(), "1/10"));
        let response = answer(&r).unwrap();
        let bytes = fs::read(r.out.join("print.stl")).unwrap();
        fs::remove_dir_all(r.out).unwrap();
        (
            response["manifest"]["printing"].clone(),
            response["work"].clone(),
            bytes,
        )
    };
    assert_eq!(run(), run());
}
