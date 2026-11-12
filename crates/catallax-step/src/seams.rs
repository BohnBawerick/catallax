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
struct Entity {
    name: String,
    refs: Vec<usize>,
}

/// Part 21 records: the text between semicolons outside strings and comments,
/// with comments removed and surrounding whitespace trimmed. Line endings and
/// layout are whitespace like any other.
fn records(file: &str) -> Result<Vec<String>, String> {
    let mut out = Vec::new();
    let mut current = String::new();
    let mut quoted = false;
    let mut chars = file.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '\'' => {
                quoted = !quoted;
                current.push(c);
            }
            '/' if !quoted && chars.peek() == Some(&'*') => {
                chars.next();
                let mut closed = false;
                while let Some(c) = chars.next() {
                    if c == '*' && chars.peek() == Some(&'/') {
                        chars.next();
                        closed = true;
                        break;
                    }
                }
                if !closed {
                    return Err("a comment is not closed".into());
                }
                current.push(' ');
            }
            ';' if !quoted => out.push(std::mem::take(&mut current).trim().to_string()),
            _ => current.push(c),
        }
    }
    if quoted {
        return Err("a string is not closed".into());
    }
    Ok(out)
}

/// The entities of every DATA section. A file with no DATA section, a record
/// there that is not an instance, or a section without its ENDSEC is refused:
/// an unread file is never taken as one that meets the contract.
fn entities(file: &str) -> Result<BTreeMap<usize, Entity>, String> {
    let mut out = BTreeMap::new();
    let mut data = false;
    let mut sections = 0;
    for record in records(file)? {
        if !data {
            let keyword = record.split(['(', ' ', '\t', '\r', '\n']).next();
            if keyword == Some("DATA") {
                data = true;
                sections += 1;
            }
            continue;
        }
        if record == "ENDSEC" {
            data = false;
            continue;
        }
        let instance = record.strip_prefix('#').and_then(|r| r.split_once('='));
        let Some((id, body)) = instance else {
            return Err(format!(
                "a DATA record is not an entity instance: `{record}`"
            ));
        };
        let id: usize = id
            .trim()
            .parse()
            .map_err(|_| format!("unreadable entity name `#{}`", id.trim()))?;
        let body = body.trim();
        let entity = Entity {
            name: body[..body.find('(').unwrap_or(body.len())]
                .trim()
                .to_string(),
            refs: references(body),
        };
        if out.insert(id, entity).is_some() {
            return Err(format!("entity #{id} is defined twice"));
        }
    }
    if sections == 0 {
        return Err("there is no DATA section".into());
    }
    if data {
        return Err("the DATA section has no ENDSEC".into());
    }
    Ok(out)
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
    let unread = |why: String| {
        Refusal::new(
            "step.seam-unsupported",
            format!("the seam check could not read the file: {why}"),
        )
    };
    let entities = entities(file).map_err(unread)?;
    if !entities.values().any(|e| e.name == "ADVANCED_FACE") {
        return Err(unread("it has no ADVANCED_FACE".into()));
    }
    let name = |id: &usize| entities.get(id).map_or("", |e| e.name.as_str());
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

    fn bare(band: &str) -> String {
        band.replace(
            "#20=SURFACE_CURVE('',#80,(#30),.PCURVE_S1.);",
            "#20=CIRCLE('',#80,1.);",
        )
    }

    #[test]
    fn crlf_files_are_read_not_passed_unread() {
        let crlf = |text: String| text.replace('\n', "\r\n");
        assert!(check(&crlf(file(BAND))).is_ok());
        let refusal = check(&crlf(file(&bare(BAND)))).unwrap_err();
        assert_eq!(refusal.code, "step.seam-unsupported");
        assert!(refusal.detail.contains("1 of 1"), "{}", refusal.detail);
    }

    #[test]
    fn other_layouts_are_read_alike() {
        // Records split across lines, several to a line, with spaces and comments.
        let layout = |band: &str| {
            let text = file(band)
                .replace(";\n#", "; /* next */ #")
                .replace("=", " =\n  ")
                .replace(",(", ",\n(");
            text.replace("DATA;", "DATA ;")
        };
        assert!(check(&layout(BAND)).is_ok(), "{}", layout(BAND));
        let refusal = check(&layout(&bare(BAND))).unwrap_err();
        assert!(refusal.detail.contains("1 of 1"), "{}", refusal.detail);
    }

    #[test]
    fn unread_files_refuse() {
        for text in [
            // No DATA section.
            BAND.to_string(),
            // A DATA section with no faces.
            file("#1=CARTESIAN_POINT('',(0.,0.,0.));\n"),
            // A record that is not an instance.
            file(&format!("{BAND}JUNK;\n")),
            // An entity defined twice.
            file(&format!("{BAND}#3=PLANE('',#9);\n")),
            // A string that is not closed.
            file(&BAND.replace("'face_#1'", "'face_#1")),
            // No ENDSEC.
            format!("{HEAD}{BAND}"),
        ] {
            let refusal = check(&text).unwrap_err();
            assert_eq!(refusal.code, "step.seam-unsupported");
            assert!(
                refusal.detail.contains("could not read"),
                "{}",
                refusal.detail
            );
        }
    }

    #[test]
    fn quoted_text_is_not_a_reference() {
        let named = BAND.replace("'face_#1'", "'#30 #31'");
        assert!(check(&file(&named)).is_ok());
    }
}
