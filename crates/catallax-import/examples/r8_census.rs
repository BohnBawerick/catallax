//! R8 observation tool. Production import policy and geometry stay unchanged.

use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};

use catallax_import::{ImportRefusal, import_step_bodies, sha256_hex};
use catallax_numbers::with_operation_budget;
use catallax_recipe::StepId;
use catallax_topology::Model;
use serde_json::{Value, json};

// Reuse the production syntax parser, without stopping at geometry recognition.
#[path = "../src/part21.rs"]
mod part21;
use part21::{Entity, Param};

#[path = "r8_census/refusal_detail.rs"]
mod refusal_detail;

// The existing eager application build allowance, catallax_app::DEFAULT_BUILD_WORK.
const IMPORT_WORK: u64 = 20_000_000;

#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod common;

fn refs(params: &[Param], out: &mut Vec<u64>) {
    for param in params {
        match param {
            Param::Ref(id) => out.push(*id),
            Param::List(items) | Param::Typed(_, items) => refs(items, out),
            _ => {}
        }
    }
}

fn face(entity: &Entity) -> bool {
    entity
        .parts
        .iter()
        .any(|(kind, _)| matches!(kind.as_str(), "ADVANCED_FACE" | "FACE_SURFACE" | "FACE"))
}

fn root(entity: &Entity) -> bool {
    entity.parts.iter().any(|(kind, _)| {
        kind.ends_with("SHAPE_REPRESENTATION")
            || matches!(
                kind.as_str(),
                "MANIFOLD_SOLID_BREP"
                    | "BREP_WITH_VOIDS"
                    | "FACETED_BREP"
                    | "SHELL_BASED_SURFACE_MODEL"
                    | "CLOSED_SHELL"
                    | "OPEN_SHELL"
            )
    })
}

fn surface_kind(entity: &Entity) -> &'static str {
    if entity.part("RATIONAL_B_SPLINE_SURFACE").is_some() {
        return "rational_bspline";
    }
    if entity.parts.iter().any(|(kind, _)| {
        matches!(
            kind.as_str(),
            "B_SPLINE_SURFACE" | "B_SPLINE_SURFACE_WITH_KNOTS" | "BEZIER_SURFACE"
        )
    }) {
        return "bspline";
    }
    match entity.type_name() {
        "PLANE" => "plane",
        "CYLINDRICAL_SURFACE" => "cylinder",
        "CONICAL_SURFACE" => "cone",
        "SPHERICAL_SURFACE" => "sphere",
        "TOROIDAL_SURFACE" => "torus",
        _ => "other",
    }
}

struct DataSections {
    count: usize,
    unsupported: usize,
    scan_complete: bool,
}

impl DataSections {
    fn single_supported(&self) -> bool {
        self.scan_complete && self.count == 1 && self.unsupported == 0
    }
}

// The production parser reads the first literal DATA; only. Scan the entire
// file before making a completeness or absence claim about that graph.
fn data_sections(bytes: &[u8]) -> DataSections {
    let mut result = DataSections {
        count: 0,
        unsupported: 0,
        scan_complete: true,
    };
    let (mut quoted, mut comment) = (false, false);
    let (mut start, mut prefix, mut len) = (0, [0_u8; 5], 0);
    let mut i = 0;
    let record =
        |prefix: &[u8], start: usize, end: usize, terminated: bool, result: &mut DataSections| {
            if prefix.len() >= 4
                && prefix[..4].eq_ignore_ascii_case(b"DATA")
                && (prefix.len() == 4 || !(prefix[4].is_ascii_alphanumeric() || prefix[4] == b'_'))
            {
                result.count += 1;
                if !terminated || &bytes[start..end] != b"DATA" {
                    result.unsupported += 1;
                }
            }
        };
    while i < bytes.len() {
        if comment {
            if bytes[i..].starts_with(b"*/") {
                comment = false;
                i += 2;
            } else {
                i += 1;
            }
            continue;
        }
        let marker = bytes[i..].starts_with(b"DATA;") || bytes[i..].starts_with(b"ENDSEC;");
        if quoted {
            // A literal marker in a string can fool the parser's substring
            // search. Its graph is not a verified section observation.
            result.scan_complete &= !marker;
            if bytes[i] == b'\'' {
                if bytes.get(i + 1) == Some(&b'\'') {
                    i += 2;
                    continue;
                }
                quoted = false;
            }
            i += 1;
            continue;
        }
        if bytes[i..].starts_with(b"/*") {
            comment = true;
            i += 2;
            continue;
        }
        if marker && len != 0 {
            result.scan_complete = false;
        }
        if bytes[i] == b';' {
            record(&prefix[..len], start, i, true, &mut result);
            len = 0;
        } else if len != 0 || !bytes[i].is_ascii_whitespace() {
            if len == 0 {
                start = i;
            }
            if len < prefix.len() {
                prefix[len] = bytes[i];
                len += 1;
            }
            quoted = bytes[i] == b'\'';
        }
        i += 1;
    }
    if len != 0 {
        record(&prefix[..len], start, bytes.len(), false, &mut result);
    }
    result.scan_complete &= !quoted && !comment && len == 0;
    result
}

fn census(bytes: &[u8]) -> Value {
    let sections = data_sections(bytes);
    let graph = match part21::parse(bytes) {
        Ok(graph) => graph,
        Err(_) => {
            return json!({"stage": "faces", "source_parse": "malformed",
                "coverage_complete": false, "face_kinds": {},
                "all_face_kinds_buildable": null, "needs_spline_faces": null,
                "data_sections": sections.count, "unsupported_data_sections": sections.unsupported,
                "data_section_scan_complete": sections.scan_complete});
        }
    };
    let mut pending: Vec<_> = graph
        .iter()
        .filter(|(_, entity)| root(entity))
        .map(|(id, _)| *id)
        .collect();
    let roots = pending.len();
    let mut visited = BTreeSet::new();
    let mut missing = BTreeSet::new();
    while let Some(id) = pending.pop() {
        if !visited.insert(id) {
            continue;
        }
        if let Some(entity) = graph.get(&id) {
            for (_, params) in &entity.parts {
                refs(params, &mut pending);
            }
        } else {
            missing.insert(id);
        }
    }
    let mut kinds = BTreeMap::<&str, usize>::new();
    let mut unreferenced = 0;
    let mut unresolved_faces = 0;
    for (&id, entity) in &graph {
        if !face(entity) {
            continue;
        }
        if !visited.contains(&id) {
            unreferenced += 1;
            continue;
        }
        // A complex face splits inherited attributes between parts. Do not
        // guess its surface slot from the simple ADVANCED_FACE layout.
        let kind = if let [(ty, args)] = entity.parts.as_slice()
            && matches!(ty.as_str(), "ADVANCED_FACE" | "FACE_SURFACE")
            && let Some(Param::Ref(surface)) = args.get(2)
        {
            graph.get(surface).map(surface_kind)
        } else {
            None
        };
        if kind.is_none() {
            unresolved_faces += 1;
        }
        *kinds.entry(kind.unwrap_or("unresolved")).or_default() += 1;
    }
    let count_type = |kind: &str| {
        graph
            .values()
            .filter(|entity| entity.part(kind).is_some())
            .count()
    };
    let instances = count_type("MAPPED_ITEM")
        + count_type("NEXT_ASSEMBLY_USAGE_OCCURRENCE")
        + count_type("REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION");
    let mut unknown_items = 0;
    for entity in graph.values() {
        for (kind, args) in &entity.parts {
            if !kind.ends_with("SHAPE_REPRESENTATION") {
                continue;
            }
            if let Some(Param::List(items)) = args.get(1) {
                for item in items {
                    let known = match item {
                        Param::Ref(id) => graph.get(id).is_some_and(|entity| {
                            matches!(
                                entity.type_name(),
                                "AXIS2_PLACEMENT_3D"
                                    | "MANIFOLD_SOLID_BREP"
                                    | "BREP_WITH_VOIDS"
                                    | "FACETED_BREP"
                                    | "SHELL_BASED_SURFACE_MODEL"
                            )
                        }),
                        _ => false,
                    };
                    unknown_items += usize::from(!known);
                }
            } else {
                unknown_items += 1;
            }
        }
    }
    let faces: usize = kinds.values().sum();
    let complete = sections.single_supported()
        && roots > 0
        && faces > 0
        && missing.is_empty()
        && unreferenced == 0
        && unresolved_faces == 0
        && instances == 0
        && unknown_items == 0;
    let spline = kinds.get("bspline").copied().unwrap_or(0)
        + kinds.get("rational_bspline").copied().unwrap_or(0);
    let analytic = kinds
        .iter()
        .filter(|(kind, _)| matches!(**kind, "plane" | "cylinder" | "cone" | "sphere" | "torus"))
        .map(|(_, count)| count)
        .sum::<usize>();
    json!({"stage": "faces", "source_parse": "parsed", "coverage_complete": complete,
        "data_sections": sections.count, "unsupported_data_sections": sections.unsupported,
        "data_section_scan_complete": sections.scan_complete,
        "referenced_faces": faces, "face_kinds": kinds, "unreferenced_faces": unreferenced,
        "missing_references": missing.len(), "unresolved_faces": unresolved_faces,
        "assembly_or_mapping_entities": instances,
        "unknown_representation_items": unknown_items,
        "solid_entities": count_type("MANIFOLD_SOLID_BREP"),
        "all_face_kinds_buildable": complete.then_some(analytic == faces),
        "needs_spline_faces": if sections.scan_complete && spline > 0 { Some(true) } else {
            (complete && !kinds.contains_key("other")).then_some(false)
        },
        "fillet_count": null, "chamfer_count": null,
        "feature_count_status": "unknown_no_feature_history_in_step"})
}

fn refusal_kind(refusal: &ImportRefusal) -> String {
    match refusal {
        ImportRefusal::Unsupported { type_name, .. } => {
            let kind = [
                "ADVANCED_FACE",
                "EDGE_CURVE",
                "SEAM_CURVE",
                "SURFACE_CURVE",
                "PCURVE",
                "B_SPLINE",
                "BOUNDED_CURVE",
                "ELLIPSE",
                "CONICAL_SURFACE",
                "CYLINDRICAL_SURFACE",
                "SPHERICAL_SURFACE",
                "TOROIDAL_SURFACE",
                "PLANE",
                "BREP_WITH_VOIDS",
                "FACETED_BREP",
                "SHELL_BASED_SURFACE_MODEL",
                "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
                "ORIENTED_EDGE",
                "EDGE_LOOP",
                "FACE_BOUND",
                "VERTEX_LOOP",
                "TRIMMED_CURVE",
                "SURFACE_OF_REVOLUTION",
                "SURFACE_OF_LINEAR_EXTRUSION",
                "DEFINITIONAL_REPRESENTATION",
            ]
            .into_iter()
            .find(|kind| type_name.starts_with(kind))
            .unwrap_or("other_step_type");
            let suffix = if type_name.contains("spline solid certification") {
                ":spline_solid_certification"
            } else {
                ""
            };
            format!("unsupported:{kind}{suffix}")
        }
        _ => serde_json::to_value(refusal).unwrap()["code"]
            .as_str()
            .unwrap()
            .to_owned(),
    }
}

fn evaluate(bytes: &[u8]) -> Value {
    let mut model = Model::new();
    let step = StepId::new("r8-import").unwrap();
    let (result, work) = with_operation_budget(IMPORT_WORK, || {
        import_step_bodies(&mut model, bytes, &sha256_hex(bytes), &step)
    });
    let result = if work.exhausted {
        Err(ImportRefusal::ResourceLimit)
    } else {
        result
    };
    let mut observation = match result {
        Ok(bodies) => {
            let tolerances=bodies.iter().flat_map(|body|model.body(*body).marks.iter())
                .filter(|mark|mark.code.starts_with("import.incidence-"))
                .map(|mark|json!({"code":mark.code,"bound_mm":mark.bound.to_decimal_exact().expect("recorded decimal"),"definition":mark.detail})).collect::<Vec<_>>();
            json!({"stage": "import", "outcome": "admitted", "bodies": bodies.len(),"local_incidence_tolerances":tolerances})
        }
        Err(refusal) => {
            // Full diagnostics stay in the private stderr file, never the table.
            eprintln!("{}", serde_json::to_string(&refusal).unwrap());
            let code = serde_json::to_value(&refusal).unwrap()["code"].clone();
            json!({"stage": "import", "outcome": "refused", "refusal_code": code,
                "refusal_kind": refusal_kind(&refusal)})
        }
    };
    observation["work_used"] = json!(work.used);
    observation["work_limit"] = json!(IMPORT_WORK);
    observation
}

fn emit(value: Value) {
    println!("{value}");
    io::stdout().flush().unwrap();
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args_os().skip(1);
    let path = args.next().ok_or("expected one STEP input")?;
    if path == "--refusal-detail" {
        let bytes = std::fs::read(args.next().ok_or("expected STEP input")?)?;
        let saved = std::fs::read(args.next().ok_or("expected saved refusal JSON")?)?;
        let refusal: ImportRefusal = serde_json::from_slice(&saved)?;
        emit(refusal_detail::inspect(&bytes, &refusal)?);
        return Ok(());
    }
    let bytes = std::fs::read(path)?;
    // Checkpoint source counts before import, so a timeout cannot erase them.
    emit(census(&bytes));
    emit(evaluate(&bytes));
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_box_counts_and_imports_through_the_recipe_entry() {
        let bytes = common::box_step(common::BoxFixture::new(["1", "2", "3"]));
        let row = census(bytes.as_bytes());
        assert_eq!(row["coverage_complete"], true);
        assert_eq!(row["face_kinds"]["plane"], 6);
        assert_eq!(row["all_face_kinds_buildable"], true);
        assert_eq!(evaluate(bytes.as_bytes())["outcome"], "admitted");
    }

    #[test]
    fn counts_continue_past_an_import_refusal_and_ignore_unused_surfaces() {
        let bytes = b"ISO-10303-21; DATA;
            #1=MANIFOLD_SOLID_BREP('',#2); #2=CLOSED_SHELL('',(#3,#4));
            #3=ADVANCED_FACE('',(),#5,.T.); #4=ADVANCED_FACE('',(),#6,.T.);
            #5=PLANE($); #6=(BOUNDED_SURFACE() B_SPLINE_SURFACE()
            RATIONAL_B_SPLINE_SURFACE()); #7=TOROIDAL_SURFACE($); ENDSEC;";
        let row = census(bytes);
        assert_eq!(row["coverage_complete"], true);
        assert_eq!(row["referenced_faces"], 2);
        assert_eq!(row["face_kinds"]["rational_bspline"], 1);
        assert_eq!(row["all_face_kinds_buildable"], false);
        assert_eq!(row["needs_spline_faces"], true);
        assert_eq!(evaluate(bytes)["outcome"], "refused");
    }

    #[test]
    fn incomplete_or_empty_traversal_never_claims_whole_part_coverage() {
        for bytes in [
            b"DATA; #1=MANIFOLD_SOLID_BREP('',#2); ENDSEC;".as_slice(),
            b"DATA; #1=ADVANCED_FACE('',(),#2,.T.); #2=PLANE($); ENDSEC;",
            b"DATA; #1=CLOSED_SHELL('',()); ENDSEC;",
            b"not STEP",
        ] {
            let row = census(bytes);
            assert_eq!(row["coverage_complete"], false);
            assert!(row["all_face_kinds_buildable"].is_null());
        }
    }

    #[test]
    fn a_later_data_section_cannot_be_hidden_by_an_analytic_prefix() {
        let bytes = b"ISO-10303-21; HEADER; ENDSEC; DATA;
            #1=MANIFOLD_SOLID_BREP('',#2); #2=CLOSED_SHELL('',(#3));
            #3=ADVANCED_FACE('',(),#4,.T.); #4=PLANE($); ENDSEC;
            DATA; #101=MANIFOLD_SOLID_BREP('',#102);
            #102=CLOSED_SHELL('',(#103));
            #103=ADVANCED_FACE('',(),#104,.T.); #104=B_SPLINE_SURFACE();
            ENDSEC; END-ISO-10303-21;";
        let row = census(bytes);
        assert_eq!(row["face_kinds"]["plane"], 1);
        assert_eq!(row["data_sections"], 2);
        assert_eq!(row["unsupported_data_sections"], 0);
        assert_eq!(row["coverage_complete"], false);
        assert!(row["all_face_kinds_buildable"].is_null());
        assert!(row["needs_spline_faces"].is_null());
    }

    #[test]
    fn unsupported_data_declarations_also_make_absence_unknown() {
        for declaration in ["DATA('second',('schema'));", "data;", "DATA ;"] {
            let bytes = format!(
                "DATA; #1=CLOSED_SHELL('',(#2));
                #2=ADVANCED_FACE('',(),#3,.T.); #3=PLANE($); ENDSEC;
                {declaration} ENDSEC;"
            );
            let row = census(bytes.as_bytes());
            assert_eq!(row["data_sections"], 2);
            assert_eq!(row["unsupported_data_sections"], 1);
            assert_eq!(row["coverage_complete"], false);
            assert!(row["all_face_kinds_buildable"].is_null());
            assert!(row["needs_spline_faces"].is_null());
        }
    }

    #[test]
    fn section_scan_handles_comments_escaped_quotes_and_unclosed_text() {
        let sections = data_sections(
            b"/* DATA; ENDSEC; */ HEADER;
            FILE_DESCRIPTION(('a ''quoted'' name'),''); ENDSEC;
            /* a comment */ DATA; #1=PLANE('ordinary label'); ENDSEC; END-ISO-10303-21;",
        );
        assert!(sections.single_supported());
        for bytes in [
            b"DATA; /* unclosed".as_slice(),
            b"DATA; 'unclosed",
            b"DATA; DATA",
        ] {
            assert!(!data_sections(bytes).scan_complete);
        }
        let quoted =
            data_sections(b"HEADER; FILE_DESCRIPTION(('DATA;'),''); ENDSEC; DATA; ENDSEC;");
        assert_eq!(quoted.count, 1);
        assert!(!quoted.scan_complete);
    }
}
