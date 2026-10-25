//! Independent finite source witnesses from primitive field domains.
//!
//! Dilation distributes over union. After proving the source attachment and
//! clearances, finite face images and complete edge/vertex tubes cover it.
//! Normal pieces overlap only in the proven source interior. The union trims
//! buried disks and rings; no unbounded support replaces a finite source patch.
use crate::{Instruction, Plan, geometry};
use catallax_features::{
    cylinder_join::{CylinderJoin, unsupported},
    cylinder_skin::CylinderSkin,
    polyhedron::{ConvexCell, PlaneFace},
};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_recipe::{BooleanKind, FaceRef, Origin, StepId, ThickenSide};
use catallax_shapes::Surface;
impl Plan {
    pub(super) fn cylinder_join(
        &self,
        node: usize,
        dec: &mut Decider,
    ) -> Result<Option<(CylinderJoin, [usize; 2])>, Refusal> {
        let Instruction::Boolean {
            kind: BooleanKind::Union,
            a,
            b,
        } = self.instructions[node]
        else {
            return Ok(None);
        };
        let mut parts = Vec::new();
        for node in [a, b] {
            let Instruction::Leaf(leaf) = self.instructions[node] else {
                return Ok(None);
            };
            // A primitive carries its complete disk/band or rectangle domain.
            // Arbitrary intersections and feature leaves need separate witnesses.
            if !matches!(self.leaves[leaf], geometry::Leaf::Primitive(_)) {
                return Ok(None);
            }
            let named = self.leaves[leaf]
                .constraints()
                .iter()
                .map(|c| {
                    let origins = self
                        .faces
                        .iter()
                        .filter(|f| f.leaf == leaf && f.surface.same_definition(c.surface()))
                        .map(|f| f.origin.clone())
                        .collect::<Vec<Origin>>();
                    (c.surface().clone(), origins)
                })
                .collect::<Vec<_>>();
            parts.push(named);
        }
        let mut cylinders = Vec::new();
        let mut plate = None;
        for (index, part) in parts.iter().enumerate() {
            if let Some(skin) = CylinderSkin::new(part, dec)? {
                cylinders.push(skin);
            } else if part.len() == 6 && part.iter().all(|(s, _)| matches!(s, Surface::Plane(_))) {
                let planes = part
                    .iter()
                    .map(|(s, origins)| {
                        let Surface::Plane(plane) = s else {
                            unreachable!()
                        };
                        PlaneFace {
                            plane: plane.clone(),
                            origins: origins.clone(),
                        }
                    })
                    .collect();
                plate = Some((index, ConvexCell::new(planes, dec)?));
            } else {
                return Ok(None);
            }
        }
        let join = if cylinders.len() == 2 {
            CylinderJoin::stepped(cylinders.remove(0), cylinders.remove(0), dec)?
        } else if cylinders.len() == 1 {
            let Some((_, cell)) = plate else {
                return Ok(None);
            };
            CylinderJoin::boss(cell, cylinders.remove(0), dec)?
        } else {
            return Ok(None);
        };
        Ok(Some((join, [a, b])))
    }
    fn join_cylinder_leaf(
        &mut self,
        skin: CylinderSkin,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut surfaces = skin.surfaces();
        for surface in &mut surfaces {
            reduce_surface(surface, dec)?;
        }
        let mut faces = Vec::new();
        for (surface, origins) in surfaces.iter().zip(skin.origins) {
            for origin in origins {
                faces.push(geometry::face(surface.clone(), origin, dec)?);
            }
        }
        self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(surfaces.into(), Some(bound.clone()))?),
            faces,
            dec,
        )
    }
    #[allow(clippy::too_many_arguments)]
    pub(super) fn cylinder_join_thicken(
        &mut self,
        join: &CylinderJoin,
        refs: &[FaceRef],
        t: &Real,
        side: ThickenSide,
        step: &StepId,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let index = join.select(refs)?;
        if let Some((prism, hole)) = join.plate_sweep(index, t, side, step, dec)? {
            let mut surfaces = Vec::new();
            let mut faces = Vec::new();
            for f in prism.faces {
                let surface = Surface::Plane(f.plane);
                for origin in f.origins {
                    faces.push(geometry::face(surface.clone(), origin, dec)?);
                }
                surfaces.push(surface);
            }
            let prism = self.leaf(
                geometry::Leaf::Cell(geometry::Cell::new(surfaces, Some(bound.clone()))?),
                faces,
                dec,
            )?;
            let hole = self.join_cylinder_leaf(hole, bound, dec)?;
            return self.boolean(BooleanKind::Subtract, prism, hole, dec);
        }
        if refs.is_empty() {
            return Err(unsupported(
                "thicken needs one complete selected join patch",
            ));
        }
        let band = join.face_band(index, t, side, step, dec)?;
        let outer = self.join_cylinder_leaf(band.outer, bound, dec)?;
        if let Some(inner) = band.inner {
            let inner = self.join_cylinder_leaf(inner, bound, dec)?;
            self.boolean(BooleanKind::Subtract, outer, inner, dec)
        } else {
            Ok(outer)
        }
    }
}

impl Plan {
    fn join_named_cell(
        &mut self,
        mut surfaces: Vec<Surface>,
        mut named: Vec<(Surface, Vec<Origin>)>,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        for surface in &mut surfaces {
            reduce_surface(surface, dec)?;
        }
        for (surface, _) in &mut named {
            reduce_surface(surface, dec)?;
        }
        let mut faces = Vec::new();
        for (surface, origins) in named {
            for origin in origins {
                faces.push(geometry::face(surface.clone(), origin, dec)?);
            }
        }
        self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(surfaces, Some(bound.clone()))?),
            faces,
            dec,
        )
    }
    pub(super) fn cylinder_join_offset(
        &mut self,
        mut node: usize,
        join: &CylinderJoin,
        d: &Real,
        step: &StepId,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        use catallax_features::polyhedral_offset::DistancePiece;
        use catallax_recipe::{EdgeName, Role};
        use catallax_shapes::Torus;
        join.clearance(d, dec)?;
        if let Some((cell, _)) = &join.plate {
            for piece in cell.outward_pieces(d, step, dec)? {
                let piece = match piece {
                    DistancePiece::Prism(mut prism) => {
                        // A rectangular face's inward normal extension stays
                        // inside the original plate. This exact overlap makes
                        // its artificial source plane an interior point of a
                        // leaf, without changing the distance neighbourhood.
                        let base = Surface::Plane(prism.faces[0].plane.clone())
                            .parallel_support(d, dec)?
                            .unwrap();
                        let Surface::Plane(base) = base else {
                            unreachable!()
                        };
                        prism.faces[0].plane = base;
                        let named = vec![(
                            Surface::Plane(prism.faces[1].plane.clone()),
                            prism.faces[1].origins.clone(),
                        )];
                        self.join_named_cell(
                            prism
                                .faces
                                .into_iter()
                                .map(|f| Surface::Plane(f.plane))
                                .collect(),
                            named,
                            bound,
                            dec,
                        )?
                    }
                    DistancePiece::Cylinder(skin) => self.join_cylinder_leaf(skin, bound, dec)?,
                    DistancePiece::Ball(surface, origins) => self.join_named_cell(
                        vec![surface.clone()],
                        vec![(surface, origins)],
                        bound,
                        dec,
                    )?,
                };
                node = self.boolean(BooleanKind::Union, node, piece, dec)?;
            }
        }
        for (_, source) in &join.bands {
            let mut outer = source.clone();
            let Some(Surface::Cylinder(cylinder)) =
                Surface::Cylinder(source.cylinder.clone()).parallel_support(d, dec)?
            else {
                unreachable!()
            };
            outer.cylinder = cylinder;
            outer.origins[0].clear();
            outer.origins[1].clear();
            let core = self.join_cylinder_leaf(outer, bound, dec)?;
            node = self.boolean(BooleanKind::Union, node, core, dec)?;
            let axial = source
                .cylinder
                .axis
                .scaled(&(d / &source.cylinder.axis.norm_squared().sqrt()));
            for end in 0..2 {
                let mut slab = source.clone();
                slab.ends = [
                    source.ends[end].plus(&axial.scaled(&Real::from_i64(-1))),
                    source.ends[end].plus(&axial),
                ];
                slab.origins = [Vec::new(), Vec::new(), Vec::new()];
                slab.origins[end] = source.origins[end].clone();
                let slab = self.join_cylinder_leaf(slab, bound, dec)?;
                node = self.boolean(BooleanKind::Union, node, slab, dec)?;
                let edge = EdgeName::new(
                    source.origins[end]
                        .first()
                        .cloned()
                        .ok_or_else(|| unsupported("a circular join needs a named cap"))?,
                    source.origins[2]
                        .first()
                        .cloned()
                        .ok_or_else(|| unsupported("a circular join needs a named wall"))?,
                )
                .map_err(|e| unsupported(&e.to_string()))?;
                let radius = source.cylinder.radius();
                let radius = if let Some(q) = dec.rational_value(&radius)? {
                    Real::from(q)
                } else {
                    radius
                };
                let surface = Surface::Torus(Torus::new(
                    source.ends[end].clone(),
                    source.cylinder.axis.clone(),
                    radius,
                    d.clone(),
                    source.cylinder.ref_dir.clone(),
                    dec,
                )?);
                let round = self.join_named_cell(
                    vec![surface.clone()],
                    vec![(
                        surface,
                        vec![Origin {
                            step: step.clone(),
                            role: Role::Round { edge },
                        }],
                    )],
                    bound,
                    dec,
                )?;
                node = self.boolean(BooleanKind::Union, node, round, dec)?;
            }
        }
        Ok(node)
    }
}

// Canonicalize only expressions proved rational. No approximation or new
// algebraic solver is introduced; irrational definitions remain unchanged.
fn reduce_surface(surface: &mut Surface, dec: &mut Decider) -> Result<(), Refusal> {
    fn value(x: &mut Real, dec: &mut Decider) -> Result<(), Refusal> {
        if x.rational_ref().is_none()
            && let Some(q) = dec.rational_value(x)?
        {
            *x = q.into();
        }
        Ok(())
    }
    match surface {
        Surface::Plane(p) => {
            for x in [
                &mut p.origin.x,
                &mut p.origin.y,
                &mut p.origin.z,
                &mut p.normal.x,
                &mut p.normal.y,
                &mut p.normal.z,
                &mut p.ref_dir.x,
                &mut p.ref_dir.y,
                &mut p.ref_dir.z,
            ] {
                value(x, dec)?;
            }
        }
        Surface::Cylinder(c) => {
            for x in [
                &mut c.origin.x,
                &mut c.origin.y,
                &mut c.origin.z,
                &mut c.axis.x,
                &mut c.axis.y,
                &mut c.axis.z,
                &mut c.ref_dir.x,
                &mut c.ref_dir.y,
                &mut c.ref_dir.z,
                &mut c.radius_sq,
            ] {
                value(x, dec)?;
            }
        }
        Surface::Torus(t) => {
            for x in [
                &mut t.center.x,
                &mut t.center.y,
                &mut t.center.z,
                &mut t.axis.x,
                &mut t.axis.y,
                &mut t.axis.z,
                &mut t.ref_dir.x,
                &mut t.ref_dir.y,
                &mut t.ref_dir.z,
                &mut t.major,
                &mut t.minor,
            ] {
                value(x, dec)?;
            }
        }
        _ => (),
    }
    Ok(())
}
