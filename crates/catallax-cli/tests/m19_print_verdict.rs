//! M19 release controls through the public print command and its artifact cleanup.
use serde_json::{Value, json};
use std::{
    fs,
    io::{BufRead, BufReader, Read, Write},
    process::{Command, Stdio},
};

fn without_clocks(value: &mut Value) {
    match value {
        Value::Object(map) => {
            map.remove("wall_ms");
            map.values_mut().for_each(without_clocks);
        }
        Value::Array(values) => values.iter_mut().for_each(without_clocks),
        _ => {}
    }
}

fn controls() -> [(&'static str, Value, &'static str); 4] {
    let p6: Value =
        serde_json::from_str(include_str!("../../../parts/P6-lattice-bracket.json")).unwrap();
    let sheet = json!({
        "format":"catallax-recipe/3", "name":"M19 sheet print control",
        "steps":[
            {"id":"box","op":"block","size":["4","1","1"],
             "pose":{"rotation":["1","0","0","0"],"translation":["-0.5","-0.5","-0.5"]}},
            {"id":"sheet","op":"gyroid","region":"box","period":"4",
             "level":{"constant":"0.25"},"kind":"sheet",
             "frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}
        ], "result":"sheet"
    });
    [
        ("P6-default", p6.clone(), "0.01"),
        ("P6-explicit-coarser", p6, "0.02"),
        ("sheet-coarse", sheet.clone(), "0.125"),
        ("sheet-default", sheet, "0.01"),
    ]
}

fn check(label: &str, cases: &[(&str, Value, &str)]) {
    let dir =
        std::env::temp_dir().join(format!("catallax-m19-print-{}-{label}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let mut reference = None;
    for workers in [1, 2, 4] {
        let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
            .env_remove("CATALLAX_THREADS")
            .args(["--threads", &workers.to_string(), "stream"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        let mut input = child.stdin.take().unwrap();
        let mut reader = BufReader::new(child.stdout.take().unwrap());
        let mut rows = Vec::new();
        for (name, recipe, accuracy) in cases {
            let out = dir.join(name);
            fs::create_dir_all(&out).unwrap();
            for file in ["print.stl", "display.obj", "answer.step"] {
                fs::write(out.join(file), b"stale artifact").unwrap();
            }
            writeln!(
                input,
                "{}",
                json!({
                    "format":"catallax-request/1", "source":"print", "max_work":1_000_000,
                    "problem":{"format":"catallax-print/1","recipe":recipe,
                        "tolerance":accuracy,"route":"field"}, "out":out
                })
            )
            .unwrap();
            input.flush().unwrap();
            let mut line = String::new();
            assert!(reader.read_line(&mut line).unwrap() > 0);
            let row: Value = serde_json::from_str(&line).unwrap();
            eprintln!(
                "{name} workers={workers} refusal={} work={}",
                row["manifest"]["refusal"], row["work"]
            );
            rows.push(row);
        }
        drop(input);
        let mut extra = String::new();
        reader.read_to_string(&mut extra).unwrap();
        assert!(child.wait().unwrap().success());
        assert!(extra.is_empty());
        assert_eq!(rows.len(), cases.len());
        let mut signature = Vec::new();
        for ((name, _, accuracy), mut row) in cases.iter().zip(rows) {
            let manifest = &row["manifest"];
            assert_eq!(manifest["status"], "refused", "{name}: {row}");
            let code = manifest["refusal"]["code"].as_str().unwrap();
            assert!(code.starts_with("field.print-") || code == "numbers.operation-budget");
            let printing = &manifest["printing"];
            assert_eq!(printing["route"], "field");
            let requested: catallax_numbers::Rational = accuracy.parse().unwrap();
            assert_eq!(printing["tolerance_mm"], requested.to_string());
            assert_eq!(printing["max_cells"], 2_000_000);
            assert_eq!(printing["max_vertices"], 1_000_000);
            assert_eq!(row["work"]["limit"], 1_000_000);
            assert_eq!(row["work"]["exhausted"], code == "numbers.operation-budget");
            assert!(printing["deviation_bound_mm"].is_null());
            assert!(printing["file"].is_null());
            for file in ["print.stl", "display.obj", "answer.step"] {
                assert!(!dir.join(name).join(file).exists(), "{name}: {file}");
            }
            let saved: Value =
                serde_json::from_slice(&fs::read(dir.join(name).join("manifest.json")).unwrap())
                    .unwrap();
            assert_eq!(&saved, manifest);
            without_clocks(&mut row);
            signature.push(row);
        }
        if let Some(expected) = &reference {
            assert_eq!(expected, &signature, "workers={workers}");
        } else {
            reference = Some(signature);
        }
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
#[ignore = "bounded P6 diagnostics; run explicitly with --ignored"]
fn bounded_p6_requests_keep_settings_and_refuse_identically_on_workers() {
    check("p6", &controls()[..2]);
}

#[test]
fn coarse_sheet_refusal_keeps_settings_cleans_files_and_matches_workers() {
    check("sheet-coarse", &controls()[2..3]);
}

#[test]
#[ignore = "bounded fine-sheet diagnostics; run explicitly with --ignored"]
fn fine_sheet_request_keeps_settings_cleans_files_and_matches_workers() {
    check("sheet-fine", &controls()[3..]);
}
