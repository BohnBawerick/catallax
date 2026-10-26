//! Assign public origins only to the exact finite facets that own their tiles.

use crate::geometry;
use catallax_features::draft::{
    arrangement::Polyhedron,
    cell::{PlanarCell, PlanePatch},
};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Point3, Surface, Vector3};

pub(super) struct Publication {
    pub surfaces: Vec<Surface>,
    pub faces: Vec<usize>,
}

fn bad(detail: &str) -> Refusal {
    super::recovery::bad(detail)
}

/// Closed membership in an outward-oriented convex polygon. All points have
/// already been proved to share the selected carrier; zero is a closed edge.
fn contains(
    vertices: &[Point3],
    cycle: &[usize],
    point: &Point3,
    normal: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for i in 0..cycle.len() {
        let a = &vertices[cycle[i]];
        let edge = vertices[cycle[(i + 1) % cycle.len()]].minus(a);
        if geometry::sign(&edge.cross(&point.minus(a)).dot(normal), dec)? == Sign::Negative {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Mutual exact containment proves equality of closed convex facets, including
/// cycles that retain different redundant collinear vertices. A support and its
/// outward sense alone cannot prove that a disconnected coplanar patch is here.
fn same_facet(
    cell: &PlanarCell,
    facet: &[usize],
    component: &Polyhedron,
    tile: &[usize],
    normal: &Vector3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if facet.len() < 3 || tile.len() < 3 {
        return Err(bad("a physical boundary tile needs a positive-area facet"));
    }
    for i in facet {
        if !contains(&component.vertices, tile, &cell.vertices[*i], normal, dec)? {
            return Ok(false);
        }
    }
    for i in tile {
        if !contains(&cell.vertices, facet, &component.vertices[*i], normal, dec)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Every arrangement boundary tile must match exactly one finite material-term
/// facet. Internal partition facets acquire no names. The term halfspaces and
/// their order are unchanged; this proves registry ownership before publication.
pub(super) fn terms(
    component: &Polyhedron,
    dec: &mut Decider,
) -> Result<Vec<Publication>, Refusal> {
    let mut owners: Vec<Vec<Option<usize>>> = component
        .faces
        .iter()
        .map(|f| vec![None; f.tiles.len()])
        .collect();
    let mut result = Vec::with_capacity(component.terms.len());
    for (term_index, term) in component.terms.iter().enumerate() {
        let surfaces: Vec<_> = term
            .sides
            .iter()
            .map(|(index, same)| {
                let mut plane = component.planes[*index].clone();
                if !same {
                    plane.normal = plane.normal.scaled(&Real::from_i64(-1));
                }
                Surface::Plane(plane)
            })
            .collect();
        let supports = surfaces
            .iter()
            .map(|surface| {
                let Surface::Plane(plane) = surface else {
                    unreachable!()
                };
                PlanePatch {
                    plane: plane.clone(),
                    origins: Vec::new(),
                }
            })
            .collect();
        let cell = PlanarCell::new(supports, dec)?
            .ok_or_else(|| bad("a swept material term has no finite solid cell"))?;
        // Arrangement terms have one bounding facet per declared side. Refuse
        // an unproved correspondence instead of indexing a merged carrier.
        if cell.supports.len() != term.sides.len() {
            return Err(bad(
                "swept term facets have no unique carrier correspondence",
            ));
        }
        let mut faces = Vec::new();
        for (face_index, face) in component.faces.iter().enumerate() {
            let Some(side) = term
                .sides
                .iter()
                .position(|s| *s == (face.support, face.same_sense))
            else {
                continue;
            };
            let facet = &cell.cycles[side];
            if facet.is_empty() {
                continue;
            }
            let normal = &cell.supports[side].plane.normal;
            let mut matched = false;
            for (tile_index, tile) in face.tiles.iter().enumerate() {
                if !same_facet(&cell, facet, component, tile, normal, dec)? {
                    continue;
                }
                if owners[face_index][tile_index].replace(term_index).is_some() {
                    return Err(bad(
                        "a finite boundary tile has more than one material owner",
                    ));
                }
                matched = true;
            }
            if matched {
                faces.push(face_index);
            }
        }
        result.push(Publication { surfaces, faces });
    }
    if owners.iter().flatten().any(Option::is_none) {
        return Err(bad(
            "a finite boundary tile has no proved material-term owner",
        ));
    }
    Ok(result)
}
