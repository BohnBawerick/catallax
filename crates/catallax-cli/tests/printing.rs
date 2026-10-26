//! Print files, not preview success, carry the accuracy and topology contract.
use catallax_cli::verification::{Request, answer};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
fn request(steps: Value, result: &str) -> Request {
    Request {
        format: "catallax-request/1".into(),
        source: "print".into(),
        problem: json!({"format":"catallax-print/1","recipe":{"format":"catallax-recipe/3","name":"print regression","steps":steps,"result":result},"tolerance":"0.01"}),
        out: std::env::temp_dir().join(format!(
            "catallax-brep-print-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        )),
        max_work: 100_000_000,
    }
}
fn block(id: &str, size: [&str; 3], translation: [&str; 3]) -> Value {
    json!({"id":id,"op":"block","size":size,"pose":{"rotation":["1","0","0","0"],"translation":translation}})
}
fn stl(path: PathBuf) -> (Vec<u8>, Vec<[[f64; 3]; 3]>) {
    let bytes = fs::read(path).unwrap();
    let count = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
    assert_eq!(bytes.len(), 84 + 50 * count);
    let triangles = bytes[84..]
        .chunks_exact(50)
        .map(|facet| {
            std::array::from_fn(|v| {
                std::array::from_fn(|axis| {
                    let offset = 12 + 12 * v + 4 * axis;
                    f64::from(f32::from_le_bytes(
                        facet[offset..offset + 4].try_into().unwrap(),
                    ))
                })
            })
        })
        .collect();
    (bytes, triangles)
}
#[test]
fn exact_box_uses_twelve_facets_with_lossless_settings_bound_to_file_bytes() {
    let mut r = request(
        json!([block("box", ["4", "4", "1"], ["0", "0", "0"])]),
        "box",
    );
    r.problem.as_object_mut().unwrap().remove("tolerance");
    // Unwrapped recipes retain the 0.01 mm default.
    r.problem = r.problem["recipe"].clone();
    let response = answer(&r).unwrap();
    assert_eq!(response["manifest"]["status"], "answered", "{response}");
    let p = &response["manifest"]["printing"];
    assert_eq!(p["route"], "brep");
    assert_eq!(p["tolerance_mm"], "1/100");
    assert_eq!(p["triangles"], 12);
    assert_eq!(p["vertices"], 8);
    assert_eq!(p["deviation_bound_mm"], "0");
    assert_eq!(p["cell_cap_applicable"], false);
    let (bytes, triangles) = stl(r.out.join("print.stl"));
    assert_eq!(&bytes[..8], b"CXPRINT1");
    assert_eq!(&bytes[8..12], &[1, 1, 1, 0]);
    assert_eq!(i64::from_le_bytes(bytes[12..20].try_into().unwrap()), 1);
    assert_eq!(u64::from_le_bytes(bytes[20..28].try_into().unwrap()), 100);
    assert_eq!(
        u64::from_le_bytes(bytes[28..36].try_into().unwrap()),
        2_000_000
    );
    assert_eq!(
        u64::from_le_bytes(bytes[36..44].try_into().unwrap()),
        1_000_000
    );
    assert_eq!(p["sha256"], format!("{:x}", Sha256::digest(&bytes)));
    let volume6: f64 = triangles
        .iter()
        .map(|[a, b, c]| {
            (0..3)
                .map(|i| a[i] * (b[(i + 1) % 3] * c[(i + 2) % 3] - b[(i + 2) % 3] * c[(i + 1) % 3]))
                .sum::<f64>()
        })
        .sum();
    assert_eq!(volume6, 96.0);
    fs::remove_dir_all(r.out).unwrap();
}
#[test]
fn decimal_coordinates_include_binary32_error_in_the_achieved_bound() {
    let r = request(
        json!([block("box", ["0.1", "0.2", "0.3"], ["0.1", "0.1", "0.1"])]),
        "box",
    );
    let response = answer(&r).unwrap();
    assert_eq!(response["manifest"]["status"], "answered", "{response}");
    let p = &response["manifest"]["printing"];
    let q = |key: &str| catallax_numbers::Rational::parse(p[key].as_str().unwrap()).unwrap();
    assert!(q("encoding_rounding_bound_mm") > catallax_numbers::Rational::zero());
    assert!(
        q("deviation_bound_mm") >= &q("tessellation_bound_mm") + &q("encoding_rounding_bound_mm")
    );
    assert!(q("deviation_bound_mm") <= q("tolerance_mm"));
    fs::remove_dir_all(r.out).unwrap();
}
#[test]
fn cavities_keep_inward_normals_and_disconnected_material_is_not_filtered() {
    let r = request(
        json!([
            block("outer",["4","4","4"],["0","0","0"]),block("void",["2","2","2"],["1","1","1"]),
            {"id":"hollow","op":"boolean","kind":"subtract","target":"outer","tool":"void"},
            block("tiny",["0.03125","0.03125","0.03125"],["6","0","0"]),
            {"id":"part","op":"boolean","kind":"union","target":"hollow","tool":"tiny"}
        ]),
        "part",
    );
    let response = answer(&r).unwrap();
    assert_eq!(response["manifest"]["status"], "answered", "{response}");
    let (_, triangles) = stl(r.out.join("print.stl"));
    let mut volumes = [0.0; 3];
    for [a, b, c] in triangles {
        let is_void = [a, b, c]
            .iter()
            .all(|p| p.iter().all(|x| (1.0..=3.0).contains(x)));
        let index = if a[0] >= 6.0 {
            2
        } else if is_void {
            1
        } else {
            0
        };
        volumes[index] += (0..3)
            .map(|i| a[i] * (b[(i + 1) % 3] * c[(i + 2) % 3] - b[(i + 2) % 3] * c[(i + 1) % 3]))
            .sum::<f64>()
            / 6.0;
    }
    assert!((volumes[0] - 64.0).abs() < 1e-10);
    assert!((volumes[1] + 8.0).abs() < 1e-10);
    assert!((volumes[2] - 0.03125_f64.powi(3)).abs() < 1e-10);
    fs::remove_dir_all(r.out).unwrap();
}
#[test]
fn binary32_sheet_collapse_refuses_without_changing_route_or_publishing_files() {
    let r = request(
        json!([block("box", ["0.000000001", "1", "1"], ["1", "0", "0"])]),
        "box",
    );
    let response = answer(&r).unwrap();
    assert_eq!(
        response["manifest"]["refusal"]["code"], "print.brep-encoding-topology",
        "{response}"
    );
    assert_eq!(response["manifest"]["printing"]["route"], "brep");
    assert!(!r.out.join("print.stl").exists());
    assert!(!r.out.join("display.obj").exists());
    fs::remove_dir_all(r.out).unwrap();
}
#[test]
fn vertex_and_work_caps_refuse_without_an_achieved_print_certificate() {
    let mut r = request(
        json!([block("box", ["4", "4", "1"], ["0", "0", "0"])]),
        "box",
    );
    r.problem["max_vertices"] = json!(4);
    let response = answer(&r).unwrap();
    assert_eq!(response["manifest"]["status"], "refused", "{response}");
    assert_eq!(response["manifest"]["printing"]["route"], "brep");
    assert!(response["manifest"]["printing"]["deviation_bound_mm"].is_null());
    r.max_work = 0;
    let response = answer(&r).unwrap();
    assert_eq!(
        response["manifest"]["refusal"]["code"],
        "numbers.operation-budget"
    );
    assert_eq!(response["manifest"]["printing"]["route"], "brep");
    assert!(!r.out.join("print.stl").exists());
    fs::remove_dir_all(r.out).unwrap();
}
#[test]
fn explicit_field_budget_is_not_retried_as_a_brep_print() {
    let mut r = request(
        json!([block("box", ["4", "4", "1"], ["0", "0", "0"])]),
        "box",
    );
    r.problem["route"] = json!("field");
    r.problem["max_cells"] = json!(1);
    let response = answer(&r).unwrap();
    assert_eq!(
        response["manifest"]["refusal"]["code"], "field.print-budget",
        "{response}"
    );
    assert_eq!(response["manifest"]["printing"]["route"], "field");
    assert_eq!(
        response["manifest"]["printing"]["cell_cap_applicable"],
        true
    );
    assert!(!r.out.join("print.stl").exists());
    fs::remove_dir_all(r.out).unwrap();
}
#[test]
fn an_unencodable_requested_rational_is_not_truncated_in_the_header() {
    let mut r = request(
        json!([block("box", ["4", "4", "1"], ["0", "0", "0"])]),
        "box",
    );
    r.problem["tolerance"] = json!("100000000000000000000");
    let response = answer(&r).unwrap();
    assert_eq!(
        response["manifest"]["refusal"]["code"], "print.metadata-range",
        "{response}"
    );
    assert!(!r.out.join("print.stl").exists());
    fs::remove_dir_all(r.out).unwrap();
}
