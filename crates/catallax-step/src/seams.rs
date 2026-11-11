//! The seam-curve exchange contract, checked on the finished file.
//!
//! Under [`crate::Seams::Required`] every face whose surface is not a `PLANE`
//! must place each of its edges on itself through a `PCURVE` of a
//! `SURFACE_CURVE` or `SEAM_CURVE`, and an edge that the face uses twice must
//! be a `SEAM_CURVE` with two such pcurves. A `VERTEX_LOOP` is admitted only on
//! a `SPHERICAL_SURFACE`, where the writer uses it for the reader's own pole.
//! The check reads the written entities, so it covers every encoder alike.

use std::collections::BTreeMap;

use catallax_numbers::Refusal;

/// One written entity: its type name (empty for a complex instance) and every
/// entity it references, in written order. Quoted text is skipped.
struct Entity<'a> {
    name: &'a str,
    refs: Vec<usize>,
}

fn entities(file: &str) -> BTreeMap<usize, Entity<'_>> {
    let data = file
        .split_once("DATA;\n")
        .map_or("", |(_, data)| data)
        .split("ENDSEC;")
        .next()
        .unwrap_or("");
    let mut out = BTreeMap::new();
    for line in data.split(";\n") {
        let Some((id, body)) = line
            .trim()
            .strip_prefix('#')
            .and_then(|l| l.split_once('='))
        else {
            continue;
        };
        let Ok(id) = id.parse() else { continue };
        let name = &body[..body.find('(').unwrap_or(body.len())];
        out.insert(
            id,
            Entity {
                name,
                refs: references(body),
            },
        );
    }
    out
}

fn references(body: &str) -> Vec<usize> {
    let mut refs = Vec::new();
    let mut quoted = false;
    for (at, c) in body.char_indices() {
        match c {
            '\'' => quoted = !quoted,
            '#' if !quoted => {
                let digits: String = body[at + 1..]
                    .chars()
                    .take_while(char::is_ascii_digit)
                    .collect();
                if let Ok(id) = digits.parse() {
                    refs.push(id);
                }
            }
            _ => {}
        }
    }
    refs
}

/// Refuses `step.seam-unsupported` when the file breaks the seam contract.
pub(crate) fn check(file: &str) -> Result<(), Refusal> {
    let entities = entities(file);
    let name = |id: &usize| entities.get(id).map_or("", |e| e.name);
    let refs = |id: &usize| entities.get(id).map_or(&[][..], |e| &e.refs[..]);
    // Pcurves of an edge's geometry that lie on `surface`.
    let pcurves_on = |geometry: &usize, surface: usize| {
        refs(geometry)
            .iter()
            .skip(1)
            .filter(|p| name(p) == "PCURVE" && refs(p).first() == Some(&surface))
            .count()
    };
    let mut curved = 0;
    let mut failed: BTreeMap<&str, usize> = BTreeMap::new();
    let mut first = None;
    for (id, face) in &entities {
        if face.name != "ADVANCED_FACE" {
            continue;
        }
        let Some((&surface, bounds)) = face.refs.split_last() else {
            continue;
        };
        let kind = name(&surface);
        if kind == "PLANE" {
            continue;
        }
        curved += 1;
        let mut edges = Vec::new();
        let mut ok = true;
        for bound in bounds {
            let Some(loop_) = refs(bound).first() else {
                ok = false;
                continue;
            };
            match name(loop_) {
                "VERTEX_LOOP" => ok &= kind == "SPHERICAL_SURFACE",
                "EDGE_LOOP" => edges.extend(
                    refs(loop_)
                        .iter()
                        .filter_map(|oriented| refs(oriented).last().copied()),
                ),
                _ => ok = false,
            }
        }
        for edge in &edges {
            let uses = edges.iter().filter(|other| *other == edge).count();
            let Some(geometry) = refs(edge).last() else {
                ok = false;
                continue;
            };
            ok &= match name(geometry) {
                "SEAM_CURVE" => pcurves_on(geometry, surface) >= uses,
                "SURFACE_CURVE" => uses == 1 && pcurves_on(geometry, surface) >= 1,
                _ => false,
            };
        }
        if !ok {
            *failed
                .entry(if kind.is_empty() { "complex" } else { kind })
                .or_default() += 1;
            first.get_or_insert(*id);
        }
    }
    let Some(first) = first else {
        return Ok(());
    };
    let total: usize = failed.values().sum();
    let kinds = failed
        .iter()
        .map(|(kind, count)| format!("{kind} {count}"))
        .collect::<Vec<_>>()
        .join(", ");
    Err(Refusal::new(
        "step.seam-unsupported",
        format!(
            "{total} of {curved} non-planar faces lack a pcurve on an edge or a SEAM_CURVE \
             at a seam ({kinds}); the first is entity #{first}"
        ),
    ))
}

#[cfg(test)]
mod tests {
    use super::check;

    const HEAD: &str = "ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n";

    fn file(body: &str) -> String {
        format!("{HEAD}{body}ENDSEC;\nEND-ISO-10303-21;\n")
    }

    // A cylinder band: two rings and a seam used twice, each with pcurves.
    const BAND: &str = "#1=ADVANCED_FACE('face_#1',(#2),#3,.T.);\n\
        #2=FACE_BOUND('',#4,.T.);\n\
        #3=CYLINDRICAL_SURFACE('',#99,1.);\n\
        #4=EDGE_LOOP('',(#5,#6,#7,#8));\n\
        #5=ORIENTED_EDGE('',*,*,#10,.T.);\n\
        #6=ORIENTED_EDGE('',*,*,#11,.T.);\n\
        #7=ORIENTED_EDGE('',*,*,#12,.F.);\n\
        #8=ORIENTED_EDGE('',*,*,#11,.F.);\n\
        #10=EDGE_CURVE('',#90,#90,#20,.T.);\n\
        #11=EDGE_CURVE('',#90,#91,#21,.T.);\n\
        #12=EDGE_CURVE('',#91,#91,#22,.T.);\n\
        #20=SURFACE_CURVE('',#80,(#30),.PCURVE_S1.);\n\
        #21=SEAM_CURVE('',#81,(#31,#32),.PCURVE_S1.);\n\
        #22=SURFACE_CURVE('',#82,(#33),.PCURVE_S1.);\n\
        #30=PCURVE('',#3,#70);\n\
        #31=PCURVE('',#3,#71);\n\
        #32=PCURVE('',#3,#72);\n\
        #33=PCURVE('',#3,#73);\n";

    #[test]
    fn a_seamed_band_passes() {
        assert!(check(&file(BAND)).is_ok());
    }

    #[test]
    fn a_bare_ring_on_a_cylinder_refuses() {
        let bare = BAND.replace(
            "#20=SURFACE_CURVE('',#80,(#30),.PCURVE_S1.);",
            "#20=CIRCLE('',#80,1.);",
        );
        let refusal = check(&file(&bare)).unwrap_err();
        assert_eq!(refusal.code, "step.seam-unsupported");
        assert!(refusal.detail.contains("1 of 1"), "{}", refusal.detail);
        assert!(refusal.detail.contains("CYLINDRICAL_SURFACE 1"));
    }

    #[test]
    fn a_seam_needs_both_pcurves_and_its_own_entity() {
        let one = BAND.replace("(#31,#32)", "(#31)");
        assert!(check(&file(&one)).is_err());
        let surface = BAND.replace("#21=SEAM_CURVE", "#21=SURFACE_CURVE");
        assert!(check(&file(&surface)).is_err());
    }

    #[test]
    fn a_pcurve_on_another_surface_does_not_count() {
        let other = BAND.replace("#30=PCURVE('',#3,", "#30=PCURVE('',#98,");
        assert!(check(&file(&other)).is_err());
    }

    #[test]
    fn planes_need_nothing_and_vertex_loops_only_on_spheres() {
        let plane = "#1=ADVANCED_FACE('a',(#2),#3,.T.);\n#2=FACE_BOUND('',#4,.T.);\n\
            #3=PLANE('',#9);\n#4=EDGE_LOOP('',(#5));\n#5=ORIENTED_EDGE('',*,*,#6,.T.);\n\
            #6=EDGE_CURVE('',#7,#7,#8,.T.);\n#8=CIRCLE('',#9,1.);\n";
        assert!(check(&file(plane)).is_ok());
        let pole = "#1=ADVANCED_FACE('a',(#2),#3,.T.);\n#2=FACE_BOUND('',#4,.T.);\n\
            #3=SPHERICAL_SURFACE('',#9,1.);\n#4=VERTEX_LOOP('',#5);\n";
        assert!(check(&file(pole)).is_ok());
        let apex = pole.replace(
            "SPHERICAL_SURFACE('',#9,1.)",
            "CONICAL_SURFACE('',#9,0.,1.)",
        );
        assert!(check(&file(&apex)).is_err());
    }

    #[test]
    fn quoted_text_is_not_a_reference() {
        let named = BAND.replace("'face_#1'", "'#30 #31'");
        assert!(check(&file(&named)).is_ok());
    }
}
