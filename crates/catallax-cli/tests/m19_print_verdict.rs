//! M19 release controls through the public print command and its artifact cleanup.

mod support;

use serde_json::{Value, json};
use std::{
    fs,
    io::{BufRead, BufReader, Read, Write},
    process::Stdio,
};

use support::bin_command;

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

fn check(label: &str, cases: &[(&str, Value, &str)], allowance: u64) {
    let dir =
        std::env::temp_dir().join(format!("catallax-m19-print-{}-{label}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let mut reference = None;
    for workers in [1, 2, 4] {
        let mut child = bin_command!("catallax-verify")
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
                    "format":"catallax-request/1", "source":"print", "max_work":allowance,
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
            assert_eq!(
                code,
                if name.starts_with("sheet") {
                    "field.print-detail-topology"
                } else {
                    "numbers.operation-budget"
                }
            );
            let printing = &manifest["printing"];
            assert_eq!(printing["route"], "field");
            let requested: catallax_numbers::Rational = accuracy.parse().unwrap();
            assert_eq!(printing["tolerance_mm"], requested.to_string());
            assert_eq!(printing["max_cells"], 2_000_000);
            assert_eq!(printing["max_vertices"], 1_000_000);
            assert_eq!(row["work"]["limit"], allowance);
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
    check("p6", &controls()[..2], 1_000_000);
}

#[test]
fn coarse_sheet_refusal_keeps_settings_cleans_files_and_matches_workers() {
    check("sheet-coarse", &controls()[2..3], 1_000_000);
}

#[test]
#[ignore = "bounded fine-sheet diagnostics; run explicitly with --ignored"]
fn fine_sheet_request_keeps_settings_cleans_files_and_matches_workers() {
    check("sheet-fine", &controls()[3..], 1_000_000);
}

#[test]
fn p6_allowance_refusal_keeps_settings_cleans_files_and_matches_workers() {
    check("p6-budget", &controls()[..2], 1000);
}

#[test]
fn fine_sheet_allowance_refusal_keeps_settings_cleans_files_and_matches_workers() {
    check("sheet-fine-budget", &controls()[3..], 1000);
}

fn stl_crossings(stl: &[u8], axis: usize, probe: [f64; 3]) -> Vec<(f64, bool)> {
    let [j, k] = [(axis + 1) % 3, (axis + 2) % 3];
    let mut hits = Vec::new();
    for facet in stl[84..].as_chunks::<50>().0 {
        let [a, b, c]: [[f64; 3]; 3] = std::array::from_fn(|vertex| {
            std::array::from_fn(|component| {
                let offset = 12 + 12 * vertex + 4 * component;
                f64::from(f32::from_le_bytes(
                    facet[offset..offset + 4].try_into().unwrap(),
                ))
            })
        });
        let [uj, uk] = [b[j] - a[j], b[k] - a[k]];
        let [vj, vk] = [c[j] - a[j], c[k] - a[k]];
        let det = uj * vk - uk * vj;
        if det == 0.0 {
            continue;
        }
        let [pj, pk] = [probe[j] - a[j], probe[k] - a[k]];
        let u = (pj * vk - pk * vj) / det;
        let v = (uj * pk - uk * pj) / det;
        if u >= 0.0 && v >= 0.0 && u + v <= 1.0 {
            hits.push((
                a[axis] + u * (b[axis] - a[axis]) + v * (c[axis] - a[axis]),
                det > 0.0,
            ));
        }
    }
    hits.sort_by(|a, b| a.0.total_cmp(&b.0));
    hits
}

fn assert_p6_material(stl: &[u8]) {
    // Each ray crosses two skin walls and the cavity between them, away from beams.
    for (axis, probe) in [(2, [4.431, 4.573, 0.0]), (1, [4.431, 0.0, 4.573])] {
        let hits = stl_crossings(stl, axis, probe);
        assert_eq!(hits.len(), 4, "cavity at {probe:?}: {hits:?}");
        for ((height, outward), (expected, sense)) in
            hits.iter()
                .zip([(0.0, false), (1.0, true), (2.0, false), (3.0, true)])
        {
            assert!(
                (height - expected).abs() <= 0.01,
                "cavity at {probe:?}: {hits:?}"
            );
            assert_eq!(*outward, sense, "cavity at {probe:?}: {hits:?}");
        }
    }
    assert!(
        stl_crossings(stl, 2, [6.03, 5.07, 0.0]).is_empty(),
        "blocked bore"
    );
    // At z = 1.93 go from the bore into the base and turn up the flange.
    for (axis, probe) in [(1, [6.03, 0.0, 1.93]), (2, [6.03, 1.93, 0.0])] {
        let hits = stl_crossings(stl, axis, probe);
        assert!(
            hits.iter().all(|(t, _)| !(1.93..=5.07).contains(t)),
            "blocked passage at {probe:?}: {hits:?}"
        );
    }
}

#[test]
#[cfg_attr(debug_assertions, ignore = "full P6 proof runs in release mode")]
fn full_p6_writes_the_same_material_certificate_and_stl_at_each_worker_count() {
    let (_, recipe, accuracy) = &controls()[0];
    let root =
        std::env::temp_dir().join(format!("catallax-p6-material-print-{}", std::process::id()));
    let mut reference = None;
    for workers in [1, 2, 4] {
        let dir = root.join(workers.to_string());
        let mut child = Command::new(env!("CARGO_BIN_EXE_catallax-verify"))
            .env_remove("CATALLAX_THREADS")
            .args(["--threads", &workers.to_string(), "stream"])
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .spawn()
            .unwrap();
        writeln!(
            child.stdin.take().unwrap(),
            "{}",
            json!({
                "format":"catallax-request/1", "source":"print", "max_work":100_000_000,
                "problem":{"format":"catallax-print/1","recipe":recipe,
                    "tolerance":accuracy,"route":"field"}, "out":dir
            })
        )
        .unwrap();
        let response = child.wait_with_output().unwrap();
        assert!(response.status.success());
        let mut row: Value = serde_json::from_slice(&response.stdout).unwrap();
        assert_eq!(row["work"]["limit"], 100_000_000);
        assert_eq!(row["work"]["target"], 20_000_000);
        assert_eq!(row["work"]["exhausted"], false);
        let manifest = &row["manifest"];
        assert_eq!(manifest["status"], "answered", "{row}");
        let record = &manifest["printing"];
        assert_eq!(record["route"], "field");
        assert_eq!(
            record["coverage"],
            "complete certified material isotopy of the field and encoded STL"
        );
        assert_eq!(record["tolerance_mm"], "1/100");
        assert_eq!(record["max_cells"], 2_000_000);
        assert_eq!(record["max_vertices"], 1_000_000);
        let bound: catallax_numbers::Rational = record["deviation_bound_mm"]
            .as_str()
            .unwrap()
            .parse()
            .unwrap();
        assert!(bound <= "1/100".parse().unwrap());
        let saved: Value =
            serde_json::from_slice(&fs::read(dir.join("manifest.json")).unwrap()).unwrap();
        assert_eq!(&saved, manifest);
        let stl = fs::read(dir.join("print.stl")).unwrap();
        assert_eq!(&stl[..12], b"CXPRINT1\x02\x01\x01\x01");
        assert_eq!(
            stl.len() as u64,
            84 + 50 * record["triangles"].as_u64().unwrap()
        );
        assert_p6_material(&stl);
        if reference.is_none() {
            eprintln!(
                "P6 work={} uncharged={} cells={} vertices={} triangles={} rounding={} bound={} sha256={}",
                row["work"]["used"],
                row["work"]["uncharged"],
                record["cells"],
                record["vertices"],
                record["triangles"],
                record["encoding_rounding_bound_mm"],
                record["deviation_bound_mm"],
                record["sha256"],
            );
        }
        without_clocks(&mut row);
        if let Some((expected_row, expected_stl)) = &reference {
            assert_eq!(expected_row, &row, "workers={workers}");
            assert_eq!(expected_stl, &stl, "workers={workers}");
        } else {
            reference = Some((row, stl));
        }
    }
    fs::remove_dir_all(root).unwrap();
}
