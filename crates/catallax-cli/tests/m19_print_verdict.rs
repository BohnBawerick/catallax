//! M19 release controls through the public print command and its artifact cleanup.
use serde_json::{Value, json};
use std::{
    fs,
    io::Write,
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

#[test]
#[ignore = "full-default print qualification; run explicitly with --ignored"]
fn full_p6_and_sheet_requests_keep_settings_and_refuse_identically_on_workers() {
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
    let cases = [
        ("P6-default", p6.clone(), "0.01"),
        ("P6-explicit-coarser", p6, "0.02"),
        ("sheet-coarse", sheet.clone(), "0.125"),
        ("sheet-default", sheet, "0.01"),
    ];
    let dir = std::env::temp_dir().join(format!("catallax-m19-print-{}", std::process::id()));
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
        for (name, recipe, accuracy) in &cases {
            let out = dir.join(name);
            fs::create_dir_all(&out).unwrap();
            for file in ["print.stl", "display.obj", "answer.step"] {
                fs::write(out.join(file), b"stale artifact").unwrap();
            }
            writeln!(
                input,
                "{}",
                json!({
                    "format":"catallax-request/1", "source":"print", "max_work":100_000_000,
                    "problem":{"format":"catallax-print/1","recipe":recipe,
                        "tolerance":accuracy,"route":"field"}, "out":out
                })
            )
            .unwrap();
        }
        drop(input);
        let output = child.wait_with_output().unwrap();
        assert!(output.status.success());
        let rows: Vec<Value> = String::from_utf8(output.stdout)
            .unwrap()
            .lines()
            .map(|line| serde_json::from_str(line).unwrap())
            .collect();
        assert_eq!(rows.len(), cases.len());
        let mut signature = Vec::new();
        for ((name, _, accuracy), mut row) in cases.iter().zip(rows) {
            let manifest = &row["manifest"];
            eprintln!(
                "{name} workers={workers} refusal={} work={}",
                manifest["refusal"], row["work"]
            );
            assert_eq!(manifest["status"], "refused", "{name}: {row}");
            assert!(
                manifest["refusal"]["code"]
                    .as_str()
                    .unwrap()
                    .starts_with("field.print-")
            );
            let printing = &manifest["printing"];
            assert_eq!(printing["route"], "field");
            let requested: catallax_numbers::Rational = accuracy.parse().unwrap();
            assert_eq!(printing["tolerance_mm"], requested.to_string());
            assert_eq!(printing["max_cells"], 2_000_000);
            assert_eq!(printing["max_vertices"], 1_000_000);
            assert_eq!(row["work"]["limit"], 100_000_000);
            assert_eq!(row["work"]["exhausted"], false);
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
