//! Read-only reference audit for a previously recorded self-intersection pair.

use catallax_import::ImportRefusal;
use serde_json::{Value, json};

use super::part21::{Graph, Param};
use super::{data_sections, surface_kind};

fn args<'a>(graph: &'a Graph, id: u64, kind: &str) -> Result<&'a [Param], &'static str> {
    match graph.get(&id).map(|entity| entity.parts.as_slice()) {
        Some([(name, args)]) if name == kind => Ok(args),
        _ => Err("unexpected_entity_layout"),
    }
}

fn reference(params: &[Param], index: usize) -> Result<u64, &'static str> {
    match params.get(index) {
        Some(Param::Ref(id)) => Ok(*id),
        _ => Err("missing_entity_reference"),
    }
}

fn references(params: &[Param], index: usize) -> Result<Vec<u64>, &'static str> {
    match params.get(index) {
        Some(Param::List(items)) => items
            .iter()
            .map(|item| reference(std::slice::from_ref(item), 0))
            .collect(),
        _ => Err("missing_reference_list"),
    }
}

fn line(graph: &Graph, mut id: u64) -> Result<bool, &'static str> {
    let mut seen = std::collections::BTreeSet::new();
    loop {
        if !seen.insert(id) {
            return Err("cyclic_curve_reference");
        }
        let entity = graph.get(&id).ok_or("missing_curve_entity")?;
        match entity.type_name() {
            "LINE" => return Ok(true),
            "CIRCLE" | "ELLIPSE" => return Ok(false),
            "SURFACE_CURVE" | "SEAM_CURVE" => {
                id = reference(args(graph, id, entity.type_name())?, 1)?;
            }
            _ => return Err("unexpected_curve_in_self_intersection_pair"),
        }
    }
}

fn face(graph: &Graph, id: u64) -> Result<(&'static str, bool), &'static str> {
    let parameters = args(graph, id, "ADVANCED_FACE")?;
    let surface = graph
        .get(&reference(parameters, 2)?)
        .ok_or("missing_face_surface")?;
    let kind = surface_kind(surface);
    if !matches!(kind, "plane" | "cylinder" | "cone" | "sphere" | "torus") {
        return Err("unexpected_surface_in_self_intersection_pair");
    }
    if kind != "plane" {
        return Ok((kind, false));
    }
    let bounds = references(parameters, 1)?;
    if bounds.is_empty() {
        return Ok((kind, false));
    }
    for bound in bounds {
        let entity = graph.get(&bound).ok_or("missing_face_bound")?;
        if !matches!(entity.type_name(), "FACE_BOUND" | "FACE_OUTER_BOUND") {
            return Err("unexpected_face_bound");
        }
        let lp = reference(args(graph, bound, entity.type_name())?, 1)?;
        let loop_entity = graph.get(&lp).ok_or("missing_face_loop")?;
        if loop_entity.type_name() == "VERTEX_LOOP" {
            continue;
        }
        for oriented in references(args(graph, lp, "EDGE_LOOP")?, 1)? {
            let edge = reference(args(graph, oriented, "ORIENTED_EDGE")?, 3)?;
            let curve = reference(args(graph, edge, "EDGE_CURVE")?, 3)?;
            if !line(graph, curve)? {
                return Ok((kind, false));
            }
        }
    }
    Ok((kind, true))
}

pub(super) fn inspect(bytes: &[u8], refusal: &ImportRefusal) -> Result<Value, &'static str> {
    let ImportRefusal::SelfIntersecting { solid, faces } = refusal else {
        return Err("expected_saved_self_intersection_refusal");
    };
    if !data_sections(bytes).single_supported() {
        return Err("unverified_data_section_scope");
    }
    let graph = super::part21::parse(bytes).map_err(|_| "source_syntax_failure")?;
    let id = |name: &str| {
        name.strip_prefix('#')
            .and_then(|name| name.parse::<u64>().ok())
            .ok_or("invalid_saved_entity_id")
    };
    let shell = reference(args(&graph, id(solid)?, "MANIFOLD_SOLID_BREP")?, 1)?;
    let shell_faces = references(args(&graph, shell, "CLOSED_SHELL")?, 1)?;
    let pair = [id(&faces[0])?, id(&faces[1])?];
    if pair[0] == pair[1] || pair.iter().any(|id| !shell_faces.contains(id)) {
        return Err("saved_pair_not_in_saved_solid");
    }
    let a = face(&graph, pair[0])?;
    let b = face(&graph, pair[1])?;
    let mut kinds = [a.0, b.0];
    kinds.sort();
    // Only the dispatch predicate is reconstructed. No intersection proof
    // is repeated or suppressed, and no counterfactual import is attempted.
    let path = if a.1 && b.1 {
        "planar_polygon_pass"
    } else {
        "curved_or_circular_edge_pass"
    };
    Ok(
        json!({"import_stage": "build_solid.check_solid_before_topology",
        "intersection_path": path, "face_kinds": kinds,
        "source_reference_audit": true,
        "counterfactual_status": "unmeasured_first_refusal_only"}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_matches_decoded_box_carriers_and_checks_membership() {
        let bytes = crate::common::box_step(crate::common::BoxFixture::new(["1", "2", "3"]));
        let file = catallax_import::read_step(bytes.as_bytes()).unwrap();
        let solid = &file.solids[0];
        let mut saved = ImportRefusal::SelfIntersecting {
            solid: solid.name.clone(),
            faces: [solid.faces[0].name.clone(), solid.faces[1].name.clone()],
        };
        let detail = inspect(bytes.as_bytes(), &saved).unwrap();
        assert_eq!(detail["face_kinds"], json!(["plane", "plane"]));
        assert_eq!(detail["intersection_path"], "planar_polygon_pass");
        assert!(solid.faces.iter().all(|face| face.loops.iter().all(|lp| {
            lp.edges
                .iter()
                .all(|edge| matches!(edge.curve, catallax_import::Curve::Line(..)))
        })));
        if let ImportRefusal::SelfIntersecting { faces, .. } = &mut saved {
            faces[1] = "#9999999".to_owned();
        }
        assert_eq!(
            inspect(bytes.as_bytes(), &saved),
            Err("saved_pair_not_in_saved_solid")
        );
    }

    #[test]
    fn circular_bounds_send_a_plane_pair_to_the_curved_pass() {
        let bytes = b"DATA; #1=MANIFOLD_SOLID_BREP('',#2); #2=CLOSED_SHELL('',(#3,#4));
            #3=ADVANCED_FACE('',(#5),#9,.T.); #4=ADVANCED_FACE('',(#5),#9,.T.);
            #5=FACE_OUTER_BOUND('',#6,.T.); #6=EDGE_LOOP('',(#7));
            #7=ORIENTED_EDGE('',*,*,#8,.T.); #8=EDGE_CURVE('',#10,#10,#11,.T.);
            #9=PLANE($); #10=VERTEX_POINT($,$); #11=SEAM_CURVE('',#12,(),.CURVE_3D.);
            #12=SURFACE_CURVE('',#13,(),.CURVE_3D.); #13=CIRCLE($); ENDSEC;";
        let saved = ImportRefusal::SelfIntersecting {
            solid: "#1".to_owned(),
            faces: ["#3".to_owned(), "#4".to_owned()],
        };
        let detail = inspect(bytes, &saved).unwrap();
        assert_eq!(detail["face_kinds"], json!(["plane", "plane"]));
        assert_eq!(detail["intersection_path"], "curved_or_circular_edge_pass");
    }
}
