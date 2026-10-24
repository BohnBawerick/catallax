//! Distance-feature reductions from the field's own constraints and face registry.

use catallax_features::{
    polyhedron::{ConvexCell, PlaneFace},
    radial::RadialBand,
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, FaceRef, Origin, StepId, ThickenSide};
use catallax_shapes::Surface;

use crate::{Instruction, Plan, geometry};

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.unsupported-envelope", detail)
}
fn positive(t: &Real, dec: &mut Decider) -> Result<(), Refusal> {
    if geometry::sign(t, dec)? != Sign::Positive {
        return Err(Refusal::new(
            "distance.non-positive-thickness",
            "thickness must be positive",
        ));
    }
    Ok(())
}
impl Plan {
    fn radial_skin(&self, node: usize) -> Option<(Surface, Vec<Origin>)> {
        let Instruction::Leaf(leaf) = self.instructions[node] else {
            return None;
        };
        if matches!(self.leaves[leaf], geometry::Leaf::Draft(_)) {
            return None;
        }
        let [constraint] = self.leaves[leaf].constraints() else {
            return None;
        };
        let surface = constraint.surface();
        if !matches!(surface, Surface::Sphere(_) | Surface::Torus(_)) {
            return None;
        }
        let mut origins: Vec<_> = self
            .faces
            .iter()
            .filter(|f| f.leaf == leaf)
            .map(|f| f.origin.clone())
            .collect();
        origins.sort();
        origins.dedup();
        Some((surface.clone(), origins))
    }
    fn cylinder_skin(
        &self,
        node: usize,
        dec: &mut Decider,
    ) -> Result<Option<catallax_features::cylinder_skin::CylinderSkin>, Refusal> {
        let Instruction::Leaf(leaf) = self.instructions[node] else {
            return Ok(None);
        };
        if matches!(
            self.leaves[leaf],
            geometry::Leaf::Feature(_) | geometry::Leaf::Draft(_)
        ) {
            return Ok(None);
        }
        let surfaces: Vec<_> = self.leaves[leaf]
            .constraints()
            .iter()
            .map(|c| {
                let origins = self
                    .faces
                    .iter()
                    .filter(|f| f.leaf == leaf && f.surface.same_definition(c.surface()))
                    .map(|f| f.origin.clone())
                    .collect();
                (c.surface().clone(), origins)
            })
            .collect();
        catallax_features::cylinder_skin::CylinderSkin::new(&surfaces, dec)
    }
    fn cylinder_leaf(
        &mut self,
        skin: catallax_features::cylinder_skin::CylinderSkin,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let surfaces = skin.surfaces();
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
    fn cone_skin(
        &self,
        node: usize,
        dec: &mut Decider,
    ) -> Result<Option<catallax_features::cone_skin::ConeSkin>, Refusal> {
        fn has_cone(plan: &Plan, node: usize) -> bool {
            match plan.instructions[node] {
                Instruction::Leaf(leaf) => plan.leaves[leaf]
                    .constraints()
                    .iter()
                    .any(|c| matches!(c.surface(), Surface::Cone(_))),
                Instruction::Boolean {
                    kind: BooleanKind::Intersect,
                    a,
                    b,
                } => has_cone(plan, a) || has_cone(plan, b),
                _ => false,
            }
        }
        if !has_cone(self, node) {
            return Ok(None);
        }
        let mut stack = vec![node];
        let mut surfaces = Vec::new();
        while let Some(node) = stack.pop() {
            match self.instructions[node] {
                Instruction::Boolean {
                    kind: BooleanKind::Intersect,
                    a,
                    b,
                } => {
                    stack.push(b);
                    stack.push(a);
                }
                Instruction::Leaf(leaf)
                    if !matches!(
                        self.leaves[leaf],
                        geometry::Leaf::Feature(_) | geometry::Leaf::Draft(_)
                    ) =>
                {
                    for c in self.leaves[leaf].constraints() {
                        let origins = self
                            .faces
                            .iter()
                            .filter(|f| f.leaf == leaf && f.surface.same_definition(c.surface()))
                            .map(|f| f.origin.clone())
                            .collect();
                        surfaces.push((c.surface().clone(), origins));
                    }
                }
                _ => return Ok(None),
            }
        }
        catallax_features::cone_skin::ConeSkin::new(&surfaces, dec)
    }
    fn cone_leaf(
        &mut self,
        skin: catallax_features::cone_skin::ConeSkin,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut faces = Vec::new();
        for f in &skin.faces {
            for origin in &f.origins {
                faces.push(geometry::face(f.surface.clone(), origin.clone(), dec)?);
            }
        }
        self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(
                skin.constraints(),
                Some(bound.clone()),
            )?),
            faces,
            dec,
        )
    }
    fn cone_sweep(
        &mut self,
        sweep: &catallax_features::cone_sweep::ConeSweep,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        use catallax_numbers::Sign;
        use catallax_shapes::{Cylinder, Plane};
        // Strict auxiliary bounds have no exposed patch or physical origin.
        let mut radius = Real::zero();
        let mut lo = sweep.vertices[0][1].clone();
        let mut hi = lo.clone();
        for p in &sweep.vertices {
            if geometry::sign(&(&p[0] - &radius), dec)? == Sign::Positive {
                radius = p[0].clone();
            }
            if geometry::sign(&(&p[1] - &lo), dec)? == Sign::Negative {
                lo = p[1].clone();
            }
            if geometry::sign(&(&p[1] - &hi), dec)? == Sign::Positive {
                hi = p[1].clone();
            }
        }
        let bounds = vec![
            Surface::Cylinder(Cylinder {
                origin: sweep.origin.clone(),
                axis: sweep.axis.clone(),
                radius_sq: (&radius + &Real::one()).square(),
                ref_dir: sweep.reference.clone(),
            }),
            Surface::Plane(Plane {
                origin: sweep.origin.plus(&sweep.axis.scaled(&(&lo - &Real::one()))),
                normal: sweep.axis.scaled(&Real::from_i64(-1)),
                ref_dir: sweep.reference.clone(),
            }),
            Surface::Plane(Plane {
                origin: sweep.origin.plus(&sweep.axis.scaled(&(&hi + &Real::one()))),
                normal: sweep.axis.clone(),
                ref_dir: sweep.reference.clone(),
            }),
        ];
        let mut node = self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(bounds, Some(bound.clone()))?),
            Vec::new(),
            dec,
        )?;
        for f in &sweep.faces {
            let mut surfaces = vec![f.surface.clone()];
            if let Surface::Cone(c) = &f.surface {
                surfaces.push(Surface::Plane(Plane {
                    origin: c.apex.clone(),
                    normal: c.axis.scaled(&Real::from_i64(-1)),
                    ref_dir: c.ref_dir.clone(),
                }));
            }
            let faces = f
                .origins
                .iter()
                .map(|o| geometry::face(f.surface.clone(), o.clone(), dec))
                .collect::<Result<_, _>>()?;
            let leaf = self.leaf(
                geometry::Leaf::Cell(geometry::Cell::new(surfaces, Some(bound.clone()))?),
                faces,
                dec,
            )?;
            node = self.boolean(
                if f.reversed {
                    BooleanKind::Subtract
                } else {
                    BooleanKind::Intersect
                },
                node,
                leaf,
                dec,
            )?;
        }
        Ok(node)
    }
    fn planar_skin(&self, node: usize, dec: &mut Decider) -> Result<ConvexCell, Refusal> {
        let mut stack = vec![node];
        let mut faces = Vec::new();
        while let Some(node) = stack.pop() {
            match self.instructions[node] {
                Instruction::Boolean {
                    kind: BooleanKind::Intersect,
                    a,
                    b,
                } => {
                    stack.push(b);
                    stack.push(a);
                }
                Instruction::Leaf(leaf) => {
                    if matches!(self.leaves[leaf], geometry::Leaf::Draft(_)) {
                        return Err(bad(
                            "a draft region needs its retained trimmed boundary before distance construction",
                        ));
                    }
                    if let geometry::Leaf::Feature(feature) = &self.leaves[leaf] {
                        let Some(planes) = feature.convex_planes(dec)? else {
                            return Err(bad(
                                "a nonconvex profile needs the complete trimmed-face envelope arrangement",
                            ));
                        };
                        for (plane, c) in planes.into_iter().zip(feature.constraints.iter()) {
                            let origins = self
                                .faces
                                .iter()
                                .filter(|f| {
                                    f.leaf == leaf && f.surface.same_definition(c.surface())
                                })
                                .map(|f| f.origin.clone())
                                .collect();
                            faces.push(PlaneFace { plane, origins });
                        }
                        continue;
                    }
                    for c in self.leaves[leaf].constraints() {
                        let Surface::Plane(plane) = c.surface() else {
                            return Err(bad(
                                "curved trimmed skins require their offset and edge-envelope arrangement",
                            ));
                        };
                        let origins = self
                            .faces
                            .iter()
                            .filter(|f| f.leaf == leaf && f.surface.same_definition(c.surface()))
                            .map(|f| f.origin.clone())
                            .collect();
                        faces.push(PlaneFace {
                            plane: plane.clone(),
                            origins,
                        });
                    }
                }
                _ => {
                    return Err(bad(
                        "a nonconvex field needs the complete trimmed-face envelope arrangement",
                    ));
                }
            }
        }
        ConvexCell::new(faces, dec)
    }
    fn distance_cell(
        &mut self,
        cell: ConvexCell,
        bound: Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if cell.is_empty() {
            return Ok(0);
        }
        let mut faces = Vec::new();
        let mut surfaces = Vec::new();
        for f in cell.faces {
            let surface = Surface::Plane(f.plane);
            for origin in f.origins {
                faces.push(geometry::face(surface.clone(), origin, dec)?);
            }
            surfaces.push(surface);
        }
        self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(surfaces, Some(bound))?),
            faces,
            dec,
        )
    }
    fn radial_leaf(
        &mut self,
        surface: Surface,
        origins: Vec<Origin>,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let faces = origins
            .into_iter()
            .map(|o| geometry::face(surface.clone(), o, dec))
            .collect::<Result<_, _>>()?;
        self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(vec![surface], Some(bound.clone()))?),
            faces,
            dec,
        )
    }
    fn radial_band(
        &mut self,
        band: RadialBand,
        bound: Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let outer = self.radial_leaf(band.outer.0, band.outer.1, &bound, dec)?;
        if let Some((surface, origins)) = band.inner {
            let inner = self.radial_leaf(surface, origins, &bound, dec)?;
            self.boolean(BooleanKind::Subtract, outer, inner, dec)
        } else {
            Ok(outer)
        }
    }
    fn unite_distance_pieces(
        &mut self,
        mut result: usize,
        pieces: Vec<catallax_features::polyhedral_offset::DistancePiece>,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        use catallax_features::polyhedral_offset::DistancePiece;
        for piece in pieces {
            let piece = match piece {
                DistancePiece::Prism(cell) => self.distance_cell(cell, bound.clone(), dec)?,
                DistancePiece::Cylinder(skin) => self.cylinder_leaf(skin, bound, dec)?,
                DistancePiece::Ball(surface, origins) => {
                    self.radial_leaf(surface, origins, bound, dec)?
                }
            };
            result = self.boolean(BooleanKind::Union, result, piece, dec)?;
        }
        Ok(result)
    }
    pub(crate) fn distance_offset(
        &mut self,
        node: usize,
        d: &Rational,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if d.is_zero() || node == 0 {
            return Ok(node);
        }
        let bound = self.node_radius(node, dec)? + d.abs();
        if let Some((surface, origins)) = self.radial_skin(node) {
            return match surface.parallel_support(&d.clone().into(), dec)? {
                Some(surface) => self.radial_leaf(surface, origins, &bound, dec),
                None => Ok(0),
            };
        }
        if let Some(skin) = self.cylinder_skin(node, dec)? {
            if d > &Rational::zero() {
                return Err(bad(
                    "outward cylinder offsets require toroidal cap-edge joins",
                ));
            }
            return match skin.eroded(&Real::from(-d), dec)? {
                Some(skin) => self.cylinder_leaf(skin, &bound, dec),
                None => Ok(0),
            };
        }
        if let Some(skin) = self.cone_skin(node, dec)? {
            if d > &Rational::zero() {
                return Err(bad(
                    "outward conical offsets need circular-edge rolls and tip sectors",
                ));
            }
            return match skin.eroded(&Real::from(-d), dec)? {
                Some(skin) => self.cone_leaf(skin, &bound, dec),
                None => Ok(0),
            };
        }
        let cell = self.planar_skin(node, dec)?;
        if d > &Rational::zero() {
            return self.unite_distance_pieces(
                node,
                cell.outward_pieces(&d.clone().into(), step, dec)?,
                &bound,
                dec,
            );
        }
        self.distance_cell(cell.eroded(&Real::from(-d), dec)?, bound, dec)
    }
    pub(crate) fn distance_shell(
        &mut self,
        node: usize,
        open: &[FaceRef],
        t: &Rational,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        positive(&t.clone().into(), dec)?;
        if node == 0 {
            if open.is_empty() {
                return Ok(0);
            }
            return Err(Refusal::new(
                "reference.lost-face",
                "an empty body has no open face",
            ));
        }
        let bound = self.node_radius(node, dec)?;
        if let Some((surface, origins)) = self.radial_skin(node) {
            for reference in open {
                if !reference.touching.is_empty()
                    || !origins.iter().any(|o| reference.face.matches(o))
                {
                    return Err(Refusal::new(
                        "reference.lost-face",
                        format!("no closed input face resolves {reference:?}"),
                    ));
                }
            }
            if !open.is_empty() {
                return Ok(0);
            }
            return self.radial_band(
                RadialBand::shell(&surface, &origins, &t.clone().into(), step, dec)?,
                bound,
                dec,
            );
        }
        if let Some(skin) = self.cylinder_skin(node, dec)? {
            let selected = skin.select(open)?;
            if selected.iter().all(|v| *v) {
                return Ok(0);
            }
            let Some(cavity) = skin.cavity(&selected, &t.clone().into(), Some(step), dec)? else {
                return Ok(node);
            };
            let cavity = self.cylinder_leaf(cavity, &bound, dec)?;
            return self.boolean(BooleanKind::Subtract, node, cavity, dec);
        }
        if let Some(skin) = self.cone_skin(node, dec)? {
            let selected = skin.select(open)?;
            if selected.iter().all(|v| *v) {
                return Ok(0);
            }
            let cavity = match skin.cavity(&selected, &t.clone().into(), step, dec) {
                Ok(Some(cavity)) => cavity,
                Ok(None) => return Ok(node),
                Err(e) => return Err(e),
            };
            let cavity = self.cone_leaf(cavity, &bound, dec)?;
            return self.boolean(BooleanKind::Subtract, node, cavity, dec);
        }
        let cell = self.planar_skin(node, dec)?;
        let selected = cell.select(open, dec)?;
        if selected.iter().all(|v| *v) {
            return Ok(0);
        }
        let cavity = match cell.cavity(&selected, &t.clone().into(), step, dec) {
            Ok(cavity) => cavity,
            Err(e) if e.code == "distance.planar-domain" => {
                let kept: Vec<_> = selected.iter().map(|open| !open).collect();
                let walls = self.unite_distance_pieces(
                    0,
                    cell.shell_pieces(&kept, &t.clone().into(), step, dec)?,
                    &bound,
                    dec,
                )?;
                return self.boolean(BooleanKind::Intersect, node, walls, dec);
            }
            Err(e) => return Err(e),
        };
        if cavity.is_empty() {
            return Ok(node);
        }
        let cavity = self.distance_cell(cavity, bound, dec)?;
        self.boolean(BooleanKind::Subtract, node, cavity, dec)
    }
    pub(crate) fn distance_thicken(
        &mut self,
        node: usize,
        refs: &[FaceRef],
        t: &Rational,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        positive(&t.clone().into(), dec)?;
        if refs.is_empty() || node == 0 {
            return Err(Refusal::new(
                "reference.lost-face",
                "thicken needs a surviving selected face",
            ));
        }
        let bound = &self.node_radius(node, dec)? + t;
        if let Some((surface, origins)) = self.radial_skin(node) {
            for reference in refs {
                if !reference.touching.is_empty()
                    || !origins.iter().any(|o| reference.face.matches(o))
                {
                    return Err(Refusal::new(
                        "reference.lost-face",
                        format!("no closed input face resolves {reference:?}"),
                    ));
                }
            }
            return self.radial_band(
                RadialBand::thicken(&surface, &origins, &t.clone().into(), side, step, dec)?,
                bound,
                dec,
            );
        }
        if let Some(skin) = self.cylinder_skin(node, dec)? {
            let selected = skin.select(refs)?;
            let indices: Vec<_> = selected
                .iter()
                .enumerate()
                .filter_map(|(i, v)| v.then_some(i))
                .collect();
            let [index] = indices.as_slice() else {
                return Err(bad(
                    "selected cylinder joins require their projection normal sectors",
                ));
            };
            let band = skin.thicken(*index, &t.clone().into(), side, step, dec)?;
            let outer = self.cylinder_leaf(band.outer, &bound, dec)?;
            return if let Some(inner) = band.inner {
                let inner = self.cylinder_leaf(inner, &bound, dec)?;
                self.boolean(BooleanKind::Subtract, outer, inner, dec)
            } else {
                Ok(outer)
            };
        }
        if let Some(skin) = self.cone_skin(node, dec)? {
            let selected = skin.select(refs)?;
            let indices: Vec<_> = selected
                .iter()
                .enumerate()
                .filter_map(|(i, v)| v.then_some(i))
                .collect();
            let [index] = indices.as_slice() else {
                return Err(bad(
                    "sharp selected conical joins require their projection normal sectors",
                ));
            };
            let sweep = skin.face_sweep(*index, &t.clone().into(), side, step, dec)?;
            return self.cone_sweep(&sweep, &bound, dec);
        }
        let cell = self.planar_skin(node, dec)?;
        let selected = cell.select(refs, dec)?;
        let indices: Vec<_> = selected
            .iter()
            .enumerate()
            .filter_map(|(i, v)| v.then_some(i))
            .collect();
        let [index] = indices.as_slice() else {
            return Err(bad(
                "sharp selected joins require their projection normal sectors",
            ));
        };
        self.distance_cell(
            cell.face_sweep(*index, &t.clone().into(), side, step, dec)?,
            bound,
            dec,
        )
    }
}
