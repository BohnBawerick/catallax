//! Shared physical boundary rows for one exact oblique cylindrical cell.

use super::*;
use catallax_shapes::{BSplineCurve, SplineBoundary, SplineSupport};

pub(super) fn check_source(model: &Model, faces: &[FaceId]) -> Result<(), Refusal> {
    let refuse = || {
        Refusal::new(
            "draft.unsupported-oblique-arrangement",
            "oblique draft source needs a complete cylinder with two shared circular cap rings",
        )
    };
    if faces.len() != 3 {
        return Err(refuse());
    }
    let mut uses = BTreeMap::<EdgeId, usize>::new();
    for &f in faces {
        let face = model.face(f);
        let loops = match model.geom().surface(face.surface) {
            Surface::Cylinder(_) => 2,
            Surface::Plane(_) => 1,
            _ => return Err(refuse()),
        };
        if face.loops.len() != loops {
            return Err(refuse());
        }
        for lp in &face.loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                return Err(refuse());
            };
            let [c] = cs.as_slice() else {
                return Err(refuse());
            };
            let edge = model.coedge(*c).edge;
            let e = model.edge(edge);
            if e.start.is_some()
                || e.end.is_some()
                || !matches!(model.geom().curve(e.curve), Curve::Circle(_))
            {
                return Err(refuse());
            }
            *uses.entry(edge).or_default() += 1;
        }
    }
    if uses.len() != 2 || uses.values().any(|n| *n != 2) {
        return Err(refuse());
    }
    Ok(())
}

pub(super) fn assemble(
    model: &mut Model,
    cell: &ObliqueCell,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let mut faces = Vec::with_capacity(3);
    let mut paths = Vec::with_capacity(2);
    for i in 0..2 {
        let curve = BSplineCurve::boundary(
            cell.support.carrier(),
            [SplineBoundary::VMin, SplineBoundary::VMax][i],
        );
        let curve = model.geom_mut().add_curve(Curve::BSpline(curve));
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        let (plane, origins) = &cell.caps[i];
        let cap_sense = if i == 0 { Sense::Reversed } else { Sense::Same };
        faces.push(crate::extrude::face(
            model,
            Surface::Plane(plane.clone()),
            Sense::Same,
            vec![vec![(edge, cap_sense)]],
            None,
            origins.clone(),
        ));
        paths.push(vec![(edge, crate::extrude::reverse(cap_sense))]);
    }
    faces.push(crate::extrude::face(
        model,
        Surface::BSpline(SplineSupport::ObliqueDraft(cell.support.clone())),
        Sense::Same,
        paths,
        None,
        cell.origins.clone(),
    ));
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    catallax_topology::DraftedCylinder::check(model, body, dec)?;
    Ok(body)
}

pub(super) fn assemble_annular(
    model: &mut Model,
    cell: &ObliqueAnnularBound,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    cell.require_separated()?;
    let mut faces = Vec::with_capacity(4);
    let mut paths = [Vec::new(), Vec::new()];
    for i in 0..2 {
        let mut cap_paths = Vec::with_capacity(2);
        for (j, ring) in cell.cap_intersections()[i].1.iter().enumerate() {
            let curve = model.geom_mut().add_curve(Curve::BSpline(ring.clone()));
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            });
            let cap_sense = if i == j { Sense::Reversed } else { Sense::Same };
            cap_paths.push(vec![(edge, cap_sense)]);
            paths[j].push(vec![(edge, crate::extrude::reverse(cap_sense))]);
        }
        let (plane, origins) = &cell.walls()[0].caps[i];
        faces.push(crate::extrude::face(
            model,
            Surface::Plane(plane.clone()),
            Sense::Same,
            cap_paths,
            None,
            origins.clone(),
        ));
    }
    for (j, (wall, paths)) in cell.walls().iter().zip(paths).enumerate() {
        faces.push(crate::extrude::face(
            model,
            Surface::BSpline(SplineSupport::ObliqueDraft(wall.support.clone())),
            if j == 0 { Sense::Same } else { Sense::Reversed },
            paths,
            None,
            wall.origins.clone(),
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    });
    catallax_topology::DraftedCylinder::check(model, body, dec)?;
    Ok(body)
}
