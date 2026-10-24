//! Independent nearest-face sweeps and complete edge tubes for circular cells.

use catallax_features::circular_skin::CircularSkin;
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, Origin, StepId};
use catallax_shapes::{Cone, Curve, Cylinder, Plane, Point3, Surface};

use crate::{Instruction, Plan, geometry};

struct SphereClip<'a> {
    center: &'a Point3,
    radius: Real,
    axial: Real,
    lower: bool,
}

impl Plan {
    pub(super) fn circular_skin(
        &self,
        node: usize,
        dec: &mut Decider,
    ) -> Result<Option<CircularSkin>, Refusal> {
        fn has_cylinder(plan: &Plan, node: usize) -> bool {
            match plan.instructions[node] {
                Instruction::Leaf(leaf) => plan.leaves[leaf]
                    .constraints()
                    .iter()
                    .any(|c| matches!(c.surface(), Surface::Cylinder(_))),
                Instruction::Boolean {
                    kind: BooleanKind::Intersect,
                    a,
                    b,
                } => has_cylinder(plan, a) || has_cylinder(plan, b),
                _ => false,
            }
        }
        if !has_cylinder(self, node) {
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
                        geometry::Leaf::Feature(_)
                            | geometry::Leaf::Draft(_)
                            | geometry::Leaf::Lattice(_)
                    ) =>
                {
                    for constraint in self.leaves[leaf].constraints() {
                        let surface = constraint.surface();
                        let origins = self
                            .faces
                            .iter()
                            .filter(|f| f.leaf == leaf && f.surface.same_definition(surface))
                            .map(|f| f.origin.clone())
                            .collect();
                        surfaces.push((surface.clone(), origins));
                    }
                }
                _ => return Ok(None),
            }
        }
        CircularSkin::new(&surfaces, dec)
    }

    fn circular_leaf(
        &mut self,
        surfaces: Vec<Surface>,
        named: &[(Surface, Vec<Origin>)],
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut faces = Vec::new();
        for (surface, origins) in named {
            for origin in origins {
                faces.push(geometry::face(surface.clone(), origin.clone(), dec)?);
            }
        }
        self.leaf(
            geometry::Leaf::Cell(geometry::Cell::new(surfaces, Some(bound.clone()))?),
            faces,
            dec,
        )
    }

    pub(super) fn circular_cell(
        &mut self,
        skin: CircularSkin,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let named: Vec<_> = skin
            .faces
            .into_iter()
            .map(|f| (f.surface, f.origins))
            .collect();
        let surfaces = named.iter().map(|(s, _)| s.clone()).collect();
        self.circular_leaf(surfaces, &named, bound, dec)
    }

    // The meridian endpoint of a sphere band sweeps a cone, not a plane.
    // Keep the positive nappe explicitly, then choose the correct half-sector.
    fn circular_sphere_clip(
        &mut self,
        node: usize,
        skin: &CircularSkin,
        clip: SphereClip<'_>,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let SphereClip {
            center,
            radius,
            axial,
            lower,
        } = clip;
        if geometry::sign(&radius, dec)? == Sign::Zero {
            return Ok(node);
        }
        let side = geometry::sign(&axial, dec)?;
        if side == Sign::Zero {
            let plane = Surface::Plane(Plane {
                origin: center.clone(),
                normal: skin
                    .axis
                    .scaled(&Real::from_i64(if lower { -1 } else { 1 })),
                ref_dir: skin.reference.clone(),
            });
            let clip = self.circular_leaf(vec![plane], &[], bound, dec)?;
            return self.boolean(BooleanKind::Intersect, node, clip, dec);
        }
        let positive = side == Sign::Positive;
        let axis = skin
            .axis
            .scaled(&Real::from_i64(if positive { 1 } else { -1 }));
        let cone = Surface::Cone(Cone {
            apex: center.clone(),
            axis: axis.clone(),
            slope: &radius / &if positive { axial } else { -axial },
            ref_dir: skin.reference.clone(),
        });
        let nappe = Surface::Plane(Plane {
            origin: center.clone(),
            normal: axis.scaled(&Real::from_i64(-1)),
            ref_dir: skin.reference.clone(),
        });
        let clip = self.circular_leaf(vec![cone, nappe], &[], bound, dec)?;
        self.boolean(
            if positive == lower {
                BooleanKind::Intersect
            } else {
                BooleanKind::Subtract
            },
            node,
            clip,
            dec,
        )
    }

    pub(super) fn circular_offset(
        &mut self,
        mut node: usize,
        skin: &CircularSkin,
        distance: &Real,
        step: &StepId,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        // Validate every carrier before adding any generated field instruction.
        let rounds = skin
            .boundaries
            .iter()
            .map(|b| skin.edge_round(b, distance, step, dec))
            .collect::<Result<Vec<_>, _>>()?;
        for (face, f) in skin.faces.iter().enumerate() {
            let shifted = f
                .surface
                .parallel_support(distance, dec)?
                .expect("positive convex offset support");
            let named = [(shifted.clone(), f.origins.clone())];
            let sweep = match &f.surface {
                Surface::Plane(plane) => {
                    let boundary = skin
                        .boundaries
                        .iter()
                        .find(|b| b.faces.contains(&face))
                        .expect("a disk has its circular boundary");
                    let Curve::Circle(circle) = &boundary.curve else {
                        unreachable!()
                    };
                    let cylinder = Surface::Cylinder(Cylinder {
                        origin: circle.center.clone(),
                        axis: skin.axis.clone(),
                        radius_sq: circle.radius_sq.clone(),
                        ref_dir: skin.reference.clone(),
                    });
                    let source = Surface::Plane(Plane {
                        origin: plane.origin.clone(),
                        normal: plane.normal.scaled(&Real::from_i64(-1)),
                        ref_dir: plane.ref_dir.clone(),
                    });
                    self.circular_leaf(vec![shifted, source, cylinder], &named, bound, dec)?
                }
                Surface::Cylinder(_) => {
                    let band = skin.bands.iter().find(|b| b.face == face).unwrap();
                    let lower = Surface::Plane(Plane {
                        origin: skin.origin.plus(&skin.axis.scaled(&band.heights[0])),
                        normal: skin.axis.scaled(&Real::from_i64(-1)),
                        ref_dir: skin.reference.clone(),
                    });
                    let upper = Surface::Plane(Plane {
                        origin: skin.origin.plus(&skin.axis.scaled(&band.heights[1])),
                        normal: skin.axis.clone(),
                        ref_dir: skin.reference.clone(),
                    });
                    let outer =
                        self.circular_leaf(vec![shifted, lower, upper], &named, bound, dec)?;
                    let inner = self.circular_leaf(vec![f.surface.clone()], &[], bound, dec)?;
                    self.boolean(BooleanKind::Subtract, outer, inner, dec)?
                }
                Surface::Sphere(sphere) => {
                    let outer = self.circular_leaf(vec![shifted], &named, bound, dec)?;
                    let inner = self.circular_leaf(vec![f.surface.clone()], &[], bound, dec)?;
                    let mut sweep = self.boolean(BooleanKind::Subtract, outer, inner, dec)?;
                    let band = skin.bands.iter().find(|b| b.face == face).unwrap();
                    for (i, height) in band.heights.iter().enumerate() {
                        let circle = skin.circle(face, height, dec)?;
                        let radius = circle.radius_sq.sqrt();
                        let axial = circle.center.minus(&sphere.center).dot(&skin.axis);
                        let clip = SphereClip {
                            center: &sphere.center,
                            radius,
                            axial,
                            lower: i == 0,
                        };
                        sweep = self.circular_sphere_clip(sweep, skin, clip, bound, dec)?;
                    }
                    sweep
                }
                _ => unreachable!(),
            };
            node = self.boolean(BooleanKind::Union, node, sweep, dec)?;
        }
        for round in rounds.into_iter().flatten() {
            let named = [(round.surface.clone(), round.origins)];
            let tube = self.circular_leaf(vec![round.surface], &named, bound, dec)?;
            node = self.boolean(BooleanKind::Union, node, tube, dec)?;
        }
        Ok(node)
    }
}
