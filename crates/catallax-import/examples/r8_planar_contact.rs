//! Audit saved planar refusal pairs without repairing or admitting a model.

use std::collections::BTreeSet;

use catallax_import::{
    Curve, Face, ImportRefusal, Solid, StepFile, Surface, check_self_intersection, read_step,
};
use catallax_numbers::Rational;
use serde_json::{Value, json};

fn off_plane_vertices(face: &Face) -> Result<usize, &'static str> {
    let Surface::Plane(p) = &face.surface else {
        return Err("saved face is not planar");
    };
    if face.loops.is_empty() {
        return Err("saved face has no polygon boundary");
    }
    let mut bad = BTreeSet::new();
    for edge in face.loops.iter().flat_map(|bound| &bound.edges) {
        if !matches!(edge.curve, Curve::Line(..)) {
            return Err("saved face has a curved boundary");
        }
        for v in [&edge.start, &edge.end] {
            let residual = (0..3).fold(Rational::zero(), |sum, i| {
                sum + &p.axis[i] * &(&v.point[i] - &p.origin[i])
            });
            if !residual.is_zero() {
                bad.insert(&v.name);
            }
        }
    }
    Ok(bad.len())
}

fn inspect(bytes: &[u8], saved: &ImportRefusal) -> Result<Value, Box<dyn std::error::Error>> {
    let ImportRefusal::SelfIntersecting { solid, faces } = saved else {
        return Err("expected a saved self-intersection refusal".into());
    };
    let file = read_step(bytes).map_err(|error| format!("{error:?}"))?;
    let source = file
        .solids
        .iter()
        .find(|s| s.name == *solid)
        .ok_or("missing saved solid")?;
    let find = |name: &str| {
        source
            .faces
            .iter()
            .find(|f| f.name == name)
            .ok_or("missing saved face")
    };
    let (a, b) = (find(&faces[0])?, find(&faces[1])?);
    if a.name == b.name {
        return Err("saved faces are identical".into());
    }
    let bad = [off_plane_vertices(a)?, off_plane_vertices(b)?];
    let edges = |f: &Face| {
        f.loops
            .iter()
            .flat_map(|bound| &bound.edges)
            .map(|e| e.edge.clone())
            .collect::<BTreeSet<_>>()
    };
    let vertices = |f: &Face| {
        f.loops
            .iter()
            .flat_map(|bound| &bound.edges)
            .flat_map(|e| [e.start.name.clone(), e.end.name.clone()])
            .collect::<BTreeSet<_>>()
    };
    let shared_edges = edges(a).intersection(&edges(b)).count();
    let shared_vertices = vertices(a).intersection(&vertices(b)).count();
    let classification = if bad.iter().any(|n| *n != 0) {
        "nonincident_boundary"
    } else {
        let pair = StepFile {
            mm_per_unit: file.mm_per_unit.clone(),
            solids: vec![Solid {
                name: source.name.clone(),
                faces: vec![a.clone(), b.clone()],
            }],
        };
        match check_self_intersection(&pair) {
            Err(ImportRefusal::SelfIntersecting { .. }) => "interior_crossing",
            Err(error) => return Err(format!("{error:?}").into()),
            Ok(()) if shared_edges > 0 => "boundary_contact",
            Ok(()) if shared_vertices > 0 => "shared_vertex_contact",
            Ok(()) => "no_interior_crossing",
        }
    };
    Ok(
        json!({"classification": classification, "shared_edges": shared_edges,
        "shared_vertices": shared_vertices, "off_plane_vertices": bad,
        "geometry_relation_verified": classification != "nonincident_boundary"}),
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err("expected STEP input and saved refusal JSON".into());
    }
    let bytes = std::fs::read(&args[0])?;
    let saved = serde_json::from_slice(&std::fs::read(&args[1])?)?;
    println!("{}", inspect(&bytes, &saved)?);
    Ok(())
}

#[cfg(test)]
#[path = "../tests/common/mod.rs"]
mod common;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_shared_edge_is_a_contact() {
        let bytes = common::box_step(common::BoxFixture::new(["1", "1", "1"]));
        let file = read_step(bytes.as_bytes()).unwrap();
        let solid = &file.solids[0];
        let saved = ImportRefusal::SelfIntersecting {
            solid: solid.name.clone(),
            faces: [solid.faces[0].name.clone(), solid.faces[2].name.clone()],
        };
        let row = inspect(bytes.as_bytes(), &saved).unwrap();
        assert_eq!(row["classification"], "boundary_contact");
        assert_eq!(row["shared_edges"], 1);
        assert_eq!(row["off_plane_vertices"], json!([0, 0]));
    }

    #[test]
    fn a_nonincident_boundary_does_not_prove_either_touching_or_crossing() {
        let bytes = common::box_step(common::BoxFixture::new(["1", "1", "1"]));
        let mut file = read_step(bytes.as_bytes()).unwrap();
        let face = &mut file.solids[0].faces[5];
        let Surface::Plane(p) = &mut face.surface else {
            unreachable!()
        };
        p.axis[0] = Rational::parse("1e-80").unwrap();
        assert_eq!(off_plane_vertices(face), Ok(2));
    }
}
