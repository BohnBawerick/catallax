//! Planar distance tubes decomposed into face prisms, edge cylinders and vertex balls.
//! Regularised union trims hidden sectors; no displaced-support mitre is accepted.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Real, Refusal};
use catallax_recipe::{EdgeName, Origin, Role, StepId, ThickenSide, VertexName};
use catallax_shapes::{Cylinder, Sphere, Surface};

use crate::{cylinder_skin::CylinderSkin, polyhedron::ConvexCell};

/// An exact piece of a boundary's distance neighbourhood.
#[derive(Clone, Debug)]
pub enum DistancePiece {
    /// Normal extrusion of a complete planar face.
    Prism(ConvexCell),
    /// Finite tube about a straight edge. Its end disks are buried by vertex balls.
    Cylinder(CylinderSkin),
    /// Ball about a vertex, with its complete generating vertex name.
    Ball(Surface, Vec<Origin>),
}

impl ConvexCell {
    /// Builds the complete outward distance neighbourhood of a convex planar skin.
    /// Union these pieces with the input solid. Every piece is within `distance`
    /// of an input face, edge or vertex; nearest-point strata prove completeness.
    pub fn outward_pieces(
        &self,
        distance: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<DistancePiece>, Refusal> {
        self.pieces(None, distance, false, step, dec)
    }

    /// Constructs the distance tubes of the kept faces. Intersect their union
    /// with the input solid to obtain the shell, including trimmed-edge witnesses.
    pub fn shell_pieces(
        &self,
        kept: &[bool],
        thickness: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<DistancePiece>, Refusal> {
        self.pieces(Some(kept), thickness, true, step, dec)
    }

    fn pieces(
        &self,
        selected: Option<&[bool]>,
        distance: &Real,
        shell: bool,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<DistancePiece>, Refusal> {
        if selected.is_some_and(|s| s.len() != self.faces.len()) {
            return Err(Refusal::new(
                "distance.face-inventory",
                "selection and planar skin differ",
            ));
        }
        let selected = |i: usize| selected.is_none_or(|s| s[i]);
        let mut pieces = Vec::new();
        let mut edges: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
        let mut vertices: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, polygon) in self.polygons.iter().enumerate() {
            if polygon.is_empty() {
                continue;
            }
            if selected(i) {
                let mut prism = self.face_sweep(
                    i,
                    distance,
                    if shell {
                        ThickenSide::Back
                    } else {
                        ThickenSide::Front
                    },
                    step,
                    dec,
                )?;
                prism.faces[1].origins = self.faces[i]
                    .origins
                    .iter()
                    .map(|of| {
                        if shell {
                            Origin {
                                step: step.clone(),
                                role: Role::Inner {
                                    of: Box::new(of.clone()),
                                },
                            }
                        } else {
                            of.clone()
                        }
                    })
                    .collect();
                for face in &mut prism.faces[2..] {
                    face.origins.clear();
                }
                pieces.push(DistancePiece::Prism(prism));
            }
            for j in 0..polygon.len() {
                let a = polygon[j];
                let b = polygon[(j + 1) % polygon.len()];
                edges.entry((a.min(b), a.max(b))).or_default().push(i);
                vertices.entry(a).or_default().push(i);
            }
        }
        for ((a, b), incident) in &edges {
            let [first, second] = incident.as_slice() else {
                return Err(Refusal::new(
                    "distance.nonmanifold-envelope",
                    "an offset edge requires two incident input faces",
                ));
            };
            if !selected(*first) && !selected(*second) {
                continue;
            }
            let origin = |i: usize| {
                self.faces[i].origins.first().cloned().ok_or_else(|| {
                    Refusal::new(
                        "distance.missing-origin",
                        "an offset generator requires a named input face",
                    )
                })
            };
            let edge = EdgeName::new(origin(*first)?, origin(*second)?)
                .map_err(|e| Refusal::new("distance.edge-name", e.to_string()))?;
            let cylinder = Cylinder::new(
                self.vertices()[*a].clone(),
                self.vertices()[*b].minus(&self.vertices()[*a]),
                distance.clone(),
                self.faces[*first].plane.normal.clone(),
                dec,
            )?;
            pieces.push(DistancePiece::Cylinder(CylinderSkin {
                cylinder,
                ends: [self.vertices()[*a].clone(), self.vertices()[*b].clone()],
                origins: [
                    Vec::new(),
                    Vec::new(),
                    vec![Origin {
                        step: step.clone(),
                        role: Role::Round { edge },
                    }],
                ],
            }));
        }
        for (vertex, mut incident) in vertices {
            incident.sort_unstable();
            incident.dedup();
            if !incident.iter().any(|i| selected(*i)) {
                continue;
            }
            let mut origins: Vec<_> = incident
                .iter()
                .filter_map(|i| self.faces[*i].origins.first().cloned())
                .collect();
            origins.sort();
            origins.dedup();
            let name = VertexName::new(origins)
                .map_err(|e| Refusal::new("distance.vertex-name", e.to_string()))?;
            let face = &self.faces[incident[0]].plane;
            let sphere = Sphere::new(
                self.vertices()[vertex].clone(),
                distance.clone(),
                face.normal.clone(),
                face.ref_dir.clone(),
                dec,
            )?;
            pieces.push(DistancePiece::Ball(
                Surface::Sphere(sphere),
                vec![Origin {
                    step: step.clone(),
                    role: Role::RoundCorner { vertex: name },
                }],
            ));
        }
        Ok(pieces)
    }
}

impl crate::draft::arrangement::Skin {
    /// Exact two-sided distance neighbourhoods of kept physical polygonal faces.
    /// Finite face prisms, complete boundary-edge tubes and vertex balls overlap
    /// as sets. Their regularised union resolves all meeting envelopes.
    pub fn shell_pieces(
        &self,
        kept: &[bool],
        thickness: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<DistancePiece>, Refusal> {
        use crate::polyhedron::PlaneFace;
        use catallax_shapes::{Plane, distance::PlanarBoundary};
        if kept.len() != self.patches.len() {
            return Err(Refusal::new(
                "distance.face-inventory",
                "selection and physical planar skin differ",
            ));
        }
        let mut pieces = Vec::new();
        let mut edges: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
        let mut vertices: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (i, patch) in self.patches.iter().enumerate() {
            for cycle in &patch.loops {
                for j in 0..cycle.len() {
                    let (a, b) = (cycle[j], cycle[(j + 1) % cycle.len()]);
                    edges.entry((a.min(b), a.max(b))).or_default().push(i);
                    vertices.entry(a).or_default().push(i);
                }
            }
            if !kept[i] {
                continue;
            }
            let boundary =
                PlanarBoundary::polygon(patch.plane.clone(), &self.vertices, &patch.loops);
            let offset = patch
                .plane
                .normal
                .scaled(&(thickness / &patch.plane.normal.norm_squared().sqrt()));
            let origins: Vec<_> = patch
                .origins
                .iter()
                .map(|of| Origin {
                    step: step.clone(),
                    role: Role::Inner {
                        of: Box::new(of.clone()),
                    },
                })
                .collect();
            for tile in boundary.tiles(dec)? {
                let mut faces: Vec<_> = [-1, 1]
                    .into_iter()
                    .map(|side| {
                        let sign = Real::from_i64(side);
                        PlaneFace {
                            plane: Plane {
                                origin: patch.plane.origin.plus(&offset.scaled(&sign)),
                                normal: patch.plane.normal.scaled(&sign),
                                ref_dir: patch.plane.ref_dir.clone(),
                            },
                            origins: origins.clone(),
                        }
                    })
                    .collect();
                for j in 0..tile.len() {
                    let edge = tile[(j + 1) % tile.len()].minus(&tile[j]);
                    faces.push(PlaneFace {
                        plane: Plane {
                            origin: tile[j].clone(),
                            normal: edge.cross(&patch.plane.normal),
                            ref_dir: patch.plane.normal.clone(),
                        },
                        origins: Vec::new(),
                    });
                }
                let prism = ConvexCell::new(faces, dec)?;
                if !prism.is_empty() {
                    pieces.push(DistancePiece::Prism(prism));
                }
            }
        }
        for ((a, b), mut incident) in edges {
            incident.sort_unstable();
            incident.dedup();
            let [first, second] = incident.as_slice() else {
                return Err(Refusal::new(
                    "distance.nonmanifold-envelope",
                    "a polygonal shell edge needs exactly two physical input faces",
                ));
            };
            if !kept[*first] && !kept[*second] {
                continue;
            }
            let mut origins = Vec::new();
            for a in &self.patches[*first].origins {
                for b in &self.patches[*second].origins {
                    let edge = EdgeName::new(a.clone(), b.clone())
                        .map_err(|e| Refusal::new("distance.edge-name", e.to_string()))?;
                    origins.push(Origin {
                        step: step.clone(),
                        role: Role::Round { edge },
                    });
                }
            }
            origins.sort();
            origins.dedup();
            let cylinder = Cylinder::new(
                self.vertices[a].clone(),
                self.vertices[b].minus(&self.vertices[a]),
                thickness.clone(),
                self.patches[*first].plane.normal.clone(),
                dec,
            )?;
            pieces.push(DistancePiece::Cylinder(CylinderSkin {
                cylinder,
                ends: [self.vertices[a].clone(), self.vertices[b].clone()],
                origins: [Vec::new(), Vec::new(), origins],
            }));
        }
        for (vertex, mut incident) in vertices {
            incident.sort_unstable();
            incident.dedup();
            if !incident.iter().any(|i| kept[*i]) || incident.len() < 3 {
                continue;
            }
            let mut names: Vec<_> = incident
                .iter()
                .flat_map(|i| self.patches[*i].origins.iter().cloned())
                .collect();
            names.sort();
            names.dedup();
            let name = VertexName::new(names)
                .map_err(|e| Refusal::new("distance.vertex-name", e.to_string()))?;
            let plane = &self.patches[incident[0]].plane;
            let sphere = Sphere::new(
                self.vertices[vertex].clone(),
                thickness.clone(),
                plane.normal.clone(),
                plane.ref_dir.clone(),
                dec,
            )?;
            pieces.push(DistancePiece::Ball(
                Surface::Sphere(sphere),
                vec![Origin {
                    step: step.clone(),
                    role: Role::RoundCorner { vertex: name },
                }],
            ));
        }
        Ok(pieces)
    }
}
