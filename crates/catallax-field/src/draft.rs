//! Analytic draft cells evaluated independently of B-rep topology.

use crate::{Face, Instruction, Location, Plan, geometry, metric::PRECISION};
use catallax_features::draft::{
    DraftDefinition,
    arrangement::Polyhedron,
    bored::{self, Bore},
    bossed::{Boss, BossedCell},
    capped::CappedCylinder,
    cell::{PlanarCell, PlanePatch},
    near,
    round::RoundCell,
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Origin, StepId};
use catallax_shapes::{Cone, Cylinder, Isometry, Plane, Point3, RPoint3, Surface, Vector3};
use std::collections::BTreeSet;

#[path = "draft_annular.rs"]
mod annular;

/// The cone a recipe tool removes with its bore, kept unchanged by the draft.
pub(crate) struct PreservedHead {
    pub(crate) walls: Vec<usize>,
    pub(crate) cone: Cone,
}

/// One bore a recipe tool removes, kept unchanged by the draft.
pub(crate) struct PreservedBore {
    pub(crate) walls: Vec<usize>,
    pub(crate) cylinder: Cylinder,
    pub(crate) head: Option<PreservedHead>,
}

/// The recipe's bore tools, removed again after their target is drafted.
pub(crate) struct PreservedBores {
    pub(crate) tools: Vec<usize>,
    pub(crate) bores: Vec<PreservedBore>,
}

/// One boss a recipe tool adds, kept unchanged by the draft.
pub(crate) struct PreservedBoss {
    pub(crate) tool: usize,
    pub(crate) boss: Boss,
    pub(crate) wall: Vec<usize>,
    pub(crate) top: Vec<usize>,
    /// The target's patch carrying the boss.
    pub(crate) shoulder: usize,
    /// Whether the preserved bore, newer than the boss, cuts through it.
    pub(crate) through: bool,
}

/// Curved source material the draft leaves in place around its planar target.
/// At least one part is present.
pub(crate) struct Preserved {
    pub(crate) bores: Option<PreservedBores>,
    pub(crate) boss: Option<Box<PreservedBoss>>,
}

#[derive(Clone, Debug)]
pub(crate) enum Cell {
    Annular(Box<catallax_features::draft::annular::AnnularCell>),
    Planar(PlanarCell),
    Round(Box<RoundCell>),
    Polyhedral(Box<Polyhedron>),
}

impl Cell {
    fn contains_interior(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Annular(c) => c.contains_interior(p, dec),
            Self::Planar(c) => c.contains_interior(p, dec),
            Self::Round(c) => c.contains_interior(p, dec),
            Self::Polyhedral(c) => c.contains_interior(p, dec),
        }
    }
    fn anchored(&self, definition: &DraftDefinition, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Annular(c) => c.anchored(&definition.neutral, dec),
            Self::Planar(c) => c.anchored(&definition.neutral, dec),
            Self::Round(c) => c.anchored(&definition.neutral, dec),
            Self::Polyhedral(c) => c.anchored(&definition.neutral, dec),
        }
    }
}

#[derive(Clone, Debug)]
enum Fast {
    Annular {
        origin: [f64; 3],
        axis: [f64; 3],
        heights: [f64; 2],
        radii: [[f64; 2]; 2],
    },
    Planar {
        vertices: Vec<[f64; 3]>,
        planes: Vec<([f64; 3], [f64; 3])>,
    },
    Round {
        origin: [f64; 3],
        axis: [f64; 3],
        heights: [f64; 2],
        radii: [f64; 2],
    },
}

#[derive(Clone, Debug)]
pub(crate) struct Region {
    pub(crate) cell: Cell,
    pub(crate) constraints: Vec<geometry::Constraint>,
    edges: Vec<(usize, usize)>,
    fast: Fast,
    bounds: Option<crate::samples::Bounds>,
}

fn point(p: &Point3) -> [f64; 3] {
    [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()]
}
fn vector(v: &Vector3) -> [f64; 3] {
    [v.x.to_f64(), v.y.to_f64(), v.z.to_f64()]
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + a[2] * b[2]
}
fn minimum(best: &mut Option<Real>, candidate: Real, dec: &mut Decider) -> Result<(), Refusal> {
    if let Some(previous) = best {
        if geometry::sign(&(&candidate - previous), dec)? == Sign::Negative {
            *previous = candidate;
        }
    } else {
        *best = Some(candidate);
    }
    Ok(())
}
fn clamp(value: Real, dec: &mut Decider) -> Result<Real, Refusal> {
    Ok(if geometry::sign(&value, dec)? == Sign::Negative {
        Real::zero()
    } else if geometry::sign(&(&value - &Real::one()), dec)? == Sign::Positive {
        Real::one()
    } else {
        value
    })
}

impl Region {
    /// Resolve a current binding to the nearest complete cached definition name.
    /// Peeling to the original creator would conflate older copy instances.
    pub(crate) fn source_origin<'a>(&self, origin: &'a Origin) -> Option<&'a Origin> {
        let mut current = origin;
        loop {
            let present = match &self.cell {
                Cell::Annular(c) => {
                    c.outer.origins.contains(current)
                        || c.inner.origins.contains(current)
                        || c.outer
                            .caps
                            .iter()
                            .flatten()
                            .any(|(_, names)| names.contains(current))
                }
                Cell::Planar(cell) => cell
                    .supports
                    .iter()
                    .any(|support| support.origins.contains(current)),
                Cell::Polyhedral(cell) => {
                    cell.faces.iter().any(|face| face.origins.contains(current))
                }
                Cell::Round(cell) => {
                    cell.origins.contains(current)
                        || cell
                            .caps
                            .iter()
                            .flatten()
                            .any(|(_, names)| names.contains(current))
                }
            };
            if present {
                return Some(current);
            }
            let catallax_recipe::Role::Copy { of, .. } = &current.role else {
                return None;
            };
            current = of;
        }
    }

    /// A strict interior point of a published boundary tile has its stored
    /// outward normal. Tile edges still require the general ray proof.
    pub(crate) fn interior_normal(
        &self,
        origin: &Origin,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Option<Vector3>, Refusal> {
        let Cell::Polyhedral(cell) = &self.cell else {
            return Ok(None);
        };
        let Some(source) = self.source_origin(origin) else {
            return Ok(None);
        };
        let p = Point3::new(
            Real::from(p[0].clone()),
            Real::from(p[1].clone()),
            Real::from(p[2].clone()),
        );
        for face in &cell.faces {
            if !face.origins.contains(source) {
                continue;
            }
            let plane = cell.plane(face);
            if geometry::sign(&p.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero {
                continue;
            }
            for tile in &face.tiles {
                let mut interior = true;
                for i in 0..tile.len() {
                    let a = &cell.vertices[tile[i]];
                    let edge = cell.vertices[tile[(i + 1) % tile.len()]].minus(a);
                    if geometry::sign(&edge.cross(&p.minus(a)).dot(&plane.normal), dec)?
                        != Sign::Positive
                    {
                        interior = false;
                        break;
                    }
                }
                if interior {
                    return Ok(Some(plane.normal));
                }
            }
        }
        Ok(None)
    }

    /// Exact membership for algebraic chart points and their one-sided germs.
    pub(crate) fn classify(
        &self,
        point: &Point3,
        direction: Option<&Vector3>,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        if let Cell::Polyhedral(cell) = &self.cell {
            return Ok(if !cell.contains(point, direction, dec)? {
                Location::Outside
            } else if cell.on_boundary(point, direction, dec)? {
                Location::Boundary
            } else {
                Location::Inside
            });
        }
        let point = [point.x.clone(), point.y.clone(), point.z.clone()];
        let mut boundary = false;
        for constraint in &self.constraints {
            let sign = if let Some(direction) = direction {
                let mut sign = Sign::Zero;
                for coefficient in geometry::point_ray_coefficients(constraint, &point, direction) {
                    sign = geometry::sign(&coefficient, dec)?;
                    if sign != Sign::Zero {
                        break;
                    }
                }
                sign
            } else {
                geometry::sign(&constraint.polynomial.eval_real(&point), dec)?
            };
            match sign {
                Sign::Positive => return Ok(Location::Outside),
                Sign::Zero => boundary = true,
                Sign::Negative => (),
            }
        }
        Ok(if boundary {
            Location::Boundary
        } else {
            Location::Inside
        })
    }

    pub(crate) fn on_face_point(
        &self,
        origin: &Origin,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let Some(source) = self.source_origin(origin) else {
            return Ok(false);
        };
        if let Cell::Polyhedral(cell) = &self.cell {
            for face in &cell.faces {
                if face.origins.contains(source) && cell.on_face(face, point, None, dec)? {
                    return Ok(true);
                }
            }
            return Ok(false);
        }
        Ok(true)
    }

    fn new(cell: Cell, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut edges = BTreeSet::new();
        let mut include = |cycles: &[Vec<usize>]| {
            for cycle in cycles {
                for i in 0..cycle.len() {
                    let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                    edges.insert((a.min(b), a.max(b)));
                }
            }
        };
        match &cell {
            Cell::Annular(_) => (),
            Cell::Planar(c) => include(&c.cycles),
            Cell::Polyhedral(c) => {
                for face in &c.faces {
                    include(&face.loops);
                }
            }
            Cell::Round(_) => (),
        }
        let fast = Self::fast_geometry(&cell);
        let mut result = Self {
            cell,
            constraints: Vec::new(),
            edges: edges.into_iter().collect(),
            fast,
            bounds: None,
        };
        result.refresh_constraints()?;
        result.bounds = Some(result.display_box(dec)?.map(|b| [b.lo(), b.hi()]));
        Ok(result)
    }
    fn fast_geometry(cell: &Cell) -> Fast {
        match cell {
            Cell::Annular(c) => Fast::Annular {
                origin: point(&c.outer.origin),
                axis: vector(&c.outer.axis),
                heights: c.outer.heights.each_ref().map(Real::to_f64),
                radii: [
                    c.outer.radii.each_ref().map(Real::to_f64),
                    c.inner.radii.each_ref().map(Real::to_f64),
                ],
            },
            Cell::Planar(cell) => Fast::Planar {
                vertices: cell.vertices.iter().map(point).collect(),
                planes: cell
                    .supports
                    .iter()
                    .map(|s| (point(&s.plane.origin), vector(&s.plane.normal)))
                    .collect(),
            },
            Cell::Polyhedral(c) => Fast::Planar {
                vertices: c.vertices.iter().map(point).collect(),
                planes: c
                    .planes
                    .iter()
                    .map(|p| (point(&p.origin), vector(&p.normal)))
                    .collect(),
            },
            Cell::Round(cell) => Fast::Round {
                origin: point(&cell.origin),
                axis: vector(&cell.axis),
                heights: cell.heights.each_ref().map(Real::to_f64),
                radii: cell.radii.each_ref().map(Real::to_f64),
            },
        }
    }
    fn refresh_constraints(&mut self) -> Result<(), Refusal> {
        self.constraints.clear();
        match &self.cell {
            Cell::Annular(c) => {
                self.constraints
                    .push(geometry::constraint(c.outer.surface.clone())?);
                self.constraints
                    .push(geometry::constraint(c.inner.surface.clone())?.reversed());
                for (plane, _) in c.outer.caps.iter().flatten() {
                    self.constraints
                        .push(geometry::constraint(Surface::Plane(plane.clone()))?);
                }
                // Positive endpoint radii and finite caps already prove that
                // both cone graphs stay on their chosen nappes throughout.
            }
            Cell::Planar(cell) => {
                for support in &cell.supports {
                    self.constraints
                        .push(geometry::constraint(Surface::Plane(support.plane.clone()))?);
                }
            }
            Cell::Polyhedral(cell) => {
                for face in &cell.faces {
                    self.constraints
                        .push(geometry::constraint(Surface::Plane(cell.plane(face)))?);
                }
            }
            Cell::Round(cell) => {
                self.constraints
                    .push(geometry::constraint(cell.surface.clone())?);
                for (plane, _) in cell.caps.iter().flatten() {
                    self.constraints
                        .push(geometry::constraint(Surface::Plane(plane.clone()))?);
                }
                if let Surface::Cone(cone) = &cell.surface {
                    self.constraints
                        .push(geometry::constraint(Surface::Plane(Plane {
                            origin: cone.apex.clone(),
                            normal: cone.axis.scaled(&Real::from_i64(-1)),
                            ref_dir: cone.ref_dir.clone(),
                        }))?);
                }
            }
        }
        Ok(())
    }
    pub(crate) fn display_box(&self, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        let mut bounds: Option<[Ball; 3]> = None;
        let mut include = |p: &Point3, radius: &Real| -> Result<(), Refusal> {
            let r = geometry::enclose(radius, 128, dec)?.hi();
            let xyz = [
                geometry::enclose(&p.x, 128, dec)?,
                geometry::enclose(&p.y, 128, dec)?,
                geometry::enclose(&p.z, 128, dec)?,
            ];
            let local = std::array::from_fn(|i| {
                Ball::from_rational(&(&xyz[i].lo() - &r), 128)
                    .hull(&Ball::from_rational(&(&xyz[i].hi() + &r), 128))
            });
            if let Some(bounds) = &mut bounds {
                for i in 0..3 {
                    bounds[i] = bounds[i].hull(&local[i]);
                }
            } else {
                bounds = Some(local);
            }
            Ok(())
        };
        match &self.cell {
            Cell::Annular(c) => {
                for i in 0..2 {
                    include(&c.outer.centre(i), &c.outer.radii[i])?;
                }
            }
            Cell::Planar(cell) => {
                for p in &cell.vertices {
                    include(p, &Real::zero())?;
                }
            }
            Cell::Polyhedral(c) => {
                for p in &c.vertices {
                    include(p, &Real::zero())?;
                }
            }
            Cell::Round(cell) => {
                for i in 0..2 {
                    include(&cell.centre(i), &cell.radii[i])?;
                }
            }
        }
        Ok(bounds.expect("published draft cells have a nonempty bounded boundary"))
    }
    pub(crate) fn display_radius(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        Ok(self
            .display_box(dec)?
            .iter()
            .fold(Rational::zero(), |sum, b| {
                sum + b.lo().abs().max(b.hi().abs())
            }))
    }
    pub(crate) fn sample_bounds(&self) -> Result<crate::samples::Bounds, Refusal> {
        self.bounds
            .clone()
            .ok_or_else(crate::samples::placement_error)
    }
    pub(crate) fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        match &mut self.cell {
            Cell::Annular(c) => {
                for wall in [&mut c.outer, &mut c.inner] {
                    annular::transform_wall(wall, map);
                }
            }
            Cell::Planar(cell) => {
                for p in &mut cell.vertices {
                    *p = map.apply_point(p);
                }
                for s in &mut cell.supports {
                    let Surface::Plane(p) = map.surface(&Surface::Plane(s.plane.clone())) else {
                        unreachable!()
                    };
                    s.plane = p;
                }
            }
            Cell::Polyhedral(c) => {
                for point in &mut c.vertices {
                    *point = map.apply_point(point);
                }
                for plane in &mut c.planes {
                    let Surface::Plane(p) = map.surface(&Surface::Plane(plane.clone())) else {
                        unreachable!()
                    };
                    *plane = p;
                }
                for term in &mut c.terms {
                    term.witness = map.apply_point(&term.witness);
                }
                if map.reversing() {
                    for face in &mut c.faces {
                        for cycle in face.loops.iter_mut().chain(&mut face.tiles) {
                            cycle.reverse();
                        }
                    }
                }
            }
            Cell::Round(cell) => {
                cell.origin = map.apply_point(&cell.origin);
                cell.axis = map.apply_vector(&cell.axis);
                cell.reference = map.apply_vector(&cell.reference);
                cell.surface = map.surface(&cell.surface);
                for (p, _) in cell.caps.iter_mut().flatten() {
                    let Surface::Plane(transformed) = map.surface(&Surface::Plane(p.clone()))
                    else {
                        unreachable!()
                    };
                    *p = transformed;
                }
            }
        }
        self.bounds = self
            .bounds
            .as_ref()
            .and_then(|bounds| crate::samples::transform(bounds, map).ok());
        self.fast = Self::fast_geometry(&self.cell);
        self.refresh_constraints()
    }
    pub(crate) fn location(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        if matches!(self.cell, Cell::Polyhedral(_)) {
            let p = Point3::new(
                p[0].clone().into(),
                p[1].clone().into(),
                p[2].clone().into(),
            );
            return self.classify(&p, None, dec);
        }
        let mut boundary = false;
        for constraint in &self.constraints {
            match crate::polynomial_sign(&constraint.polynomial, p, dec)? {
                Sign::Positive => return Ok(Location::Outside),
                Sign::Zero => boundary = true,
                Sign::Negative => (),
            }
        }
        Ok(if boundary {
            Location::Boundary
        } else {
            Location::Inside
        })
    }
    pub(crate) fn proven_distance(
        &self,
        p: &RPoint3,
        location: Location,
        dec: &mut Decider,
    ) -> Result<Ball, Refusal> {
        if location == Location::Boundary {
            return Ok(Ball::from_rational(&Rational::zero(), PRECISION));
        }
        let p = p.to_real();
        let mut best = None;
        match &self.cell {
            Cell::Annular(c) => {
                best = Some(annular::distance_squared(c, &p, dec)?);
            }
            Cell::Planar(cell) => {
                for (index, (support, cycle)) in cell.supports.iter().zip(&cell.cycles).enumerate()
                {
                    if cycle.is_empty() {
                        continue;
                    }
                    let plane = &support.plane;
                    let signed = p.minus(&plane.origin).dot(&plane.normal);
                    let projected = p.plus(
                        &plane
                            .normal
                            .scaled(&(-&signed / plane.normal.norm_squared())),
                    );
                    let mut inside = true;
                    for (other, support) in cell.supports.iter().enumerate() {
                        // The foot lies on its own support by construction.
                        if other != index
                            && geometry::sign(
                                &projected
                                    .minus(&support.plane.origin)
                                    .dot(&support.plane.normal),
                                dec,
                            )? == Sign::Positive
                        {
                            inside = false;
                            break;
                        }
                    }
                    if inside {
                        minimum(
                            &mut best,
                            signed.square() / plane.normal.norm_squared(),
                            dec,
                        )?;
                    }
                }
                for (a, b) in &self.edges {
                    let a = &cell.vertices[*a];
                    let b = &cell.vertices[*b];
                    let direction = b.minus(a);
                    let t = clamp(p.minus(a).dot(&direction) / direction.norm_squared(), dec)?;
                    let foot = a.plus(&direction.scaled(&t));
                    minimum(&mut best, p.minus(&foot).norm_squared(), dec)?;
                }
            }
            Cell::Polyhedral(c) => {
                for face in &c.faces {
                    let plane = &c.planes[face.support];
                    let signed = p.minus(&plane.origin).dot(&plane.normal);
                    let projected = p.plus(
                        &plane
                            .normal
                            .scaled(&(-&signed / plane.normal.norm_squared())),
                    );
                    if c.on_face(face, &projected, None, dec)? {
                        minimum(
                            &mut best,
                            signed.square() / plane.normal.norm_squared(),
                            dec,
                        )?;
                    }
                }
                for (a, b) in &self.edges {
                    let a = &c.vertices[*a];
                    let direction = c.vertices[*b].minus(a);
                    let t = clamp(p.minus(a).dot(&direction) / direction.norm_squared(), dec)?;
                    minimum(
                        &mut best,
                        p.minus(&a.plus(&direction.scaled(&t))).norm_squared(),
                        dec,
                    )?;
                }
            }
            Cell::Round(cell) => {
                let d = p.minus(&cell.origin);
                let h = d.dot(&cell.axis);
                let radius = d.cross(&cell.axis).norm_squared().sqrt();
                let dr = &cell.radii[1] - &cell.radii[0];
                let dh = &cell.heights[1] - &cell.heights[0];
                let numerator = &((&radius - &cell.radii[0]) * dr.clone())
                    + &((&h - &cell.heights[0]) * dh.clone());
                let t = clamp(numerator / (dr.square() + dh.square()), dec)?;
                let rr = &radius - &(&cell.radii[0] + &(&t * &dr));
                let hh = &h - &(&cell.heights[0] + &(&t * &dh));
                minimum(&mut best, rr.square() + hh.square(), dec)?;
                for i in 0..2 {
                    if cell.caps[i].is_none() {
                        continue;
                    }
                    let radial = &radius - &cell.radii[i];
                    let radial = if geometry::sign(&radial, dec)? == Sign::Positive {
                        radial
                    } else {
                        Real::zero()
                    };
                    minimum(
                        &mut best,
                        radial.square() + (&h - &cell.heights[i]).square(),
                        dec,
                    )?;
                }
            }
        }
        let distance = geometry::enclose(&best.expect("nonempty boundary").sqrt(), PRECISION, dec)?;
        Ok(if location == Location::Inside {
            -distance
        } else {
            distance
        })
    }
    pub(crate) fn fast_distance(&self, p: [f64; 3]) -> f64 {
        let mut best = f64::INFINITY;
        let inside = match (&self.cell, &self.fast) {
            (
                Cell::Annular(_),
                Fast::Annular {
                    origin,
                    axis,
                    heights,
                    radii,
                },
            ) => {
                let (distance, inside) =
                    annular::fast_distance(p, *origin, *axis, *heights, *radii);
                best = distance;
                inside
            }
            (Cell::Planar(_), Fast::Planar { vertices, planes }) => {
                let inside = planes
                    .iter()
                    .all(|(origin, normal)| dot(sub(p, *origin), *normal) <= 0.0);
                for (origin, normal) in planes {
                    let value = dot(sub(p, *origin), *normal);
                    let norm = dot(*normal, *normal);
                    let foot = std::array::from_fn(|i| p[i] - normal[i] * (value / norm));
                    if planes.iter().all(|(o, n)| dot(sub(foot, *o), *n) <= 0.0) {
                        best = best.min(value * value / norm);
                    }
                }
                for (a, b) in &self.edges {
                    let direction = sub(vertices[*b], vertices[*a]);
                    let t = (dot(sub(p, vertices[*a]), direction) / dot(direction, direction))
                        .clamp(0.0, 1.0);
                    let delta =
                        std::array::from_fn(|i| p[i] - (vertices[*a][i] + t * direction[i]));
                    best = best.min(dot(delta, delta));
                }
                inside
            }
            (Cell::Polyhedral(c), Fast::Planar { vertices, planes }) => {
                for face in &c.faces {
                    let (origin, normal) = planes[face.support];
                    let value = dot(sub(p, origin), normal);
                    let norm = dot(normal, normal);
                    let foot = std::array::from_fn(|i| p[i] - normal[i] * (value / norm));
                    let inside = face.tiles.iter().any(|tile| {
                        tile.iter().enumerate().all(|(i, a)| {
                            let edge = sub(vertices[tile[(i + 1) % tile.len()]], vertices[*a]);
                            let d = sub(foot, vertices[*a]);
                            let cross = [
                                edge[1] * d[2] - edge[2] * d[1],
                                edge[2] * d[0] - edge[0] * d[2],
                                edge[0] * d[1] - edge[1] * d[0],
                            ];
                            dot(cross, normal) * if face.same_sense { 1.0 } else { -1.0 } >= 0.0
                        })
                    });
                    if inside {
                        best = best.min(value * value / norm);
                    }
                }
                for (a, b) in &self.edges {
                    let direction = sub(vertices[*b], vertices[*a]);
                    let t = (dot(sub(p, vertices[*a]), direction) / dot(direction, direction))
                        .clamp(0.0, 1.0);
                    let delta =
                        std::array::from_fn(|i| p[i] - (vertices[*a][i] + t * direction[i]));
                    best = best.min(dot(delta, delta));
                }
                c.terms.iter().any(|term| {
                    term.sides.iter().all(|(i, same)| {
                        let (o, n) = planes[*i];
                        dot(sub(p, o), n) * if *same { 1.0 } else { -1.0 } <= 0.0
                    })
                })
            }
            (
                Cell::Round(cell),
                Fast::Round {
                    origin,
                    axis,
                    heights,
                    radii,
                },
            ) => {
                let delta = sub(p, *origin);
                let h = dot(delta, *axis);
                let radial = sub(delta, axis.map(|a| a * h));
                let radius = dot(radial, radial).sqrt();
                let dr = radii[1] - radii[0];
                let dh = heights[1] - heights[0];
                let t = (((radius - radii[0]) * dr + (h - heights[0]) * dh) / (dr * dr + dh * dh))
                    .clamp(0.0, 1.0);
                best = (radius - (radii[0] + t * dr)).powi(2) + (h - (heights[0] + t * dh)).powi(2);
                for i in 0..2 {
                    if cell.caps[i].is_some() {
                        best = best
                            .min((radius - radii[i]).max(0.0).powi(2) + (h - heights[i]).powi(2));
                    }
                }
                h >= heights[0]
                    && h <= heights[1]
                    && radius <= radii[0] + ((h - heights[0]) / dh) * dr
            }
            _ => unreachable!("draft cache matches its exact geometry"),
        };
        if inside { -best.sqrt() } else { best.sqrt() }
    }
}
enum Candidate {
    Region(Cell),
    Capped(Box<CappedCylinder>),
    Bossed(Box<BossedCell>),
}

impl Candidate {
    fn contains_interior(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Region(cell) => cell.contains_interior(p, dec),
            Self::Capped(cell) => cell.contains_interior(p, dec),
            Self::Bossed(cell) => cell.contains_interior(p, dec),
        }
    }
    fn anchored(&self, definition: &DraftDefinition, dec: &mut Decider) -> Result<bool, Refusal> {
        match self {
            Self::Region(cell) => cell.anchored(definition, dec),
            Self::Capped(cell) => cell.anchored(&definition.neutral, dec),
            Self::Bossed(cell) => cell.anchored(&definition.neutral, dec),
        }
    }
}

struct Candidates {
    components: Vec<Candidate>,
    source_faces: Vec<Vec<usize>>,
}

fn physical_sources(
    skin: &catallax_features::draft::arrangement::Skin,
    faces: &[Face],
    source_faces: &[usize],
) -> Vec<Vec<usize>> {
    skin.patches
        .iter()
        .map(|patch| {
            source_faces
                .iter()
                .copied()
                .filter(|f| patch.origins.contains(&faces[*f].origin))
                .collect()
        })
        .collect()
}

/// Cached geometry keeps its definition names. A whole-body copy adds a common
/// outer copy prefix; use the nearest matching suffix, not the original creator,
/// so previously distinct copy instances never acquire one another's names.
pub(crate) fn remap_skin_origins(
    skin: &mut catallax_features::draft::arrangement::Skin,
    public: &[&Origin],
) -> Result<(), Refusal> {
    let mut bindings = std::collections::BTreeMap::<&Origin, (usize, BTreeSet<&Origin>)>::new();
    for origin in public {
        let mut current = *origin;
        let mut depth = 0;
        loop {
            let binding = bindings
                .entry(current)
                .or_insert_with(|| (depth, BTreeSet::new()));
            if depth < binding.0 {
                binding.0 = depth;
                binding.1.clear();
            }
            if depth == binding.0 {
                binding.1.insert(*origin);
            }
            let catallax_recipe::Role::Copy { of, .. } = &current.role else {
                break;
            };
            current = of;
            depth += 1;
        }
    }
    for patch in &mut skin.patches {
        let mut names = BTreeSet::new();
        for origin in &patch.origins {
            let binding = bindings.get(origin).ok_or_else(|| {
                Refusal::new(
                    "draft.missing-source-origin",
                    "cached source patch has no current named binding",
                )
            })?;
            names.extend(binding.1.iter().map(|name| (**name).clone()));
        }
        patch.origins = names.into_iter().collect();
    }
    Ok(())
}

fn selected_patches(
    skin: &catallax_features::draft::arrangement::Skin,
    faces: &[Face],
    source_faces: &[usize],
    selected: &mut BTreeSet<usize>,
) -> BTreeSet<usize> {
    let patches: BTreeSet<_> = skin
        .patches
        .iter()
        .enumerate()
        .filter_map(|(i, patch)| {
            selected
                .iter()
                .any(|f| patch.origins.contains(&faces[*f].origin))
                .then_some(i)
        })
        .collect();
    bind_patches(skin, &patches, faces, source_faces, selected);
    patches
}

fn bind_patches(
    skin: &catallax_features::draft::arrangement::Skin,
    patches: &BTreeSet<usize>,
    faces: &[Face],
    source_faces: &[usize],
    selected: &mut BTreeSet<usize>,
) {
    for i in patches {
        for f in source_faces {
            if skin.patches[*i].origins.contains(&faces[*f].origin) {
                selected.insert(*f);
            }
        }
    }
}

impl Plan {
    /// Chooses the uniquely nearest closed trimmed face of a plane-capped axial
    /// cylinder or of a coaxial annulus, then every alias of that physical face.
    fn nearest_axial(
        &self,
        source_faces: &[usize],
        matches: &[usize],
        near: &catallax_recipe::Point3,
        dec: &mut Decider,
    ) -> Result<Vec<usize>, Refusal> {
        let unsupported = || {
            Refusal::new(
                "draft.unsupported-reference",
                "curved field references need exact patch adjacency and distance",
            )
        };
        // One entry per physical cylinder and per physical cap: aliases of a
        // definition are the same face, two walls of different radii are not.
        let mut cylinders: Vec<(usize, &Cylinder)> = Vec::new();
        let mut caps = Vec::new();
        for f in source_faces {
            match &self.faces[*f].surface {
                Surface::Cylinder(c) => {
                    if !cylinders.iter().any(|(g, _)| {
                        self.faces[*g]
                            .support
                            .same_definition(&self.faces[*f].support)
                    }) {
                        cylinders.push((*f, c));
                    }
                }
                Surface::Plane(p) => {
                    if !caps.iter().any(|(g, _): &(usize, &Plane)| {
                        self.faces[*g]
                            .support
                            .same_definition(&self.faces[*f].support)
                    }) {
                        caps.push((*f, p));
                    }
                }
                _ => return Err(unsupported()),
            }
        }
        let [(_, first), (_, second)] = caps.as_slice() else {
            return Err(unsupported());
        };
        let walls = match cylinders.as_slice() {
            [_] => None,
            [(_, a), (_, b)] => Some(near::annulus_walls(a, b, dec)?),
            _ => return Err(unsupported()),
        };
        let [x, y, z] = near.0.each_ref().map(|n| n.value().clone().into());
        let point = Point3::new(x, y, z);
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for f in matches {
            match groups.iter_mut().find(|group| {
                self.faces[group[0]]
                    .support
                    .same_definition(&self.faces[*f].support)
            }) {
                Some(group) => group.push(*f),
                None => groups.push(vec![*f]),
            }
        }
        let mut distances = Vec::with_capacity(groups.len());
        for group in &groups {
            distances.push(match &self.faces[group[0]].surface {
                Surface::Plane(cap) => match walls {
                    Some(walls) => near::annular_cap_distance_squared(walls, cap, &point, dec)?,
                    None => near::cap_distance_squared(cylinders[0].1, cap, &point, dec)?,
                },
                Surface::Cylinder(wall) => {
                    near::wall_distance_squared(wall, [first, second], &point, dec)?
                }
                _ => return Err(unsupported()),
            });
        }
        Ok(groups.swap_remove(near::unique_nearest(&distances, dec)?))
    }

    /// Returns distinct connected material interiors, not a convex decomposition.
    fn draft_cells(
        &self,
        source: (
            Option<usize>,
            Option<(
                &catallax_features::draft::arrangement::Skin,
                &BTreeSet<usize>,
            )>,
        ),
        source_faces: &[usize],
        selected: &mut BTreeSet<usize>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Candidates, Refusal> {
        let (leaf, boundary) = source;
        if let Some((skin, patches)) = boundary {
            if skin.patches.is_empty() {
                return Ok(Candidates {
                    components: Vec::new(),
                    source_faces: Vec::new(),
                });
            }
            if skin.is_convex(dec)? {
                let mut supports = Vec::with_capacity(skin.patches.len());
                for (index, patch) in skin.patches.iter().enumerate() {
                    let plane = if patches.contains(&index) {
                        let Surface::Plane(plane) = definition
                            .support(&Surface::Plane(patch.plane.clone()), true, dec)?
                            .surface
                        else {
                            unreachable!("plane draft remains planar")
                        };
                        plane
                    } else {
                        patch.plane.clone()
                    };
                    supports.push(PlanePatch {
                        plane,
                        origins: patch.origins.clone(),
                    });
                }
                let cell = PlanarCell::new(supports, dec)?.ok_or_else(|| {
                    Refusal::new(
                        "draft.reversed-solid",
                        "replaced boundary has no full-dimensional material cell",
                    )
                })?;
                return Ok(Candidates {
                    components: vec![Candidate::Region(Cell::Planar(cell))],
                    source_faces: physical_sources(skin, &self.faces, source_faces),
                });
            }
            return Ok(Candidates {
                components: skin
                    .draft(patches, definition, dec)?
                    .into_iter()
                    .map(|cell| Candidate::Region(Cell::Polyhedral(Box::new(cell))))
                    .collect(),
                source_faces: physical_sources(skin, &self.faces, source_faces),
            });
        }
        let leaf = leaf.expect("a source is either a leaf or an extracted boundary");
        let planar = source_faces
            .iter()
            .all(|f| matches!(self.faces[*f].surface, Surface::Plane(_)));
        if planar {
            let convex = match self.leaves[leaf].certify_planar_convexity(dec) {
                Ok(()) => true,
                Err(e) if e.code == "draft.unsupported-arrangement" => false,
                Err(e) => return Err(e),
            };
            if !convex {
                let origins: Vec<_> = source_faces
                    .iter()
                    .map(|f| &self.faces[*f].origin)
                    .collect();
                let skin = match &self.leaves[leaf] {
                    geometry::Leaf::Feature(f) => f.planar_skin(&origins, None, dec)?,
                    geometry::Leaf::Draft(region) => match &region.cell {
                        Cell::Polyhedral(c) => {
                            let mut skin = c.skin();
                            remap_skin_origins(&mut skin, &origins)?;
                            skin
                        }
                        _ => {
                            return Err(Refusal::new(
                                "draft.unsupported-arrangement",
                                "source has no planar boundary definition",
                            ));
                        }
                    },
                    _ => {
                        return Err(Refusal::new(
                            "draft.unsupported-arrangement",
                            "source has no planar boundary definition",
                        ));
                    }
                };
                let patches = selected_patches(&skin, &self.faces, source_faces, selected);
                return Ok(Candidates {
                    components: skin
                        .draft(&patches, definition, dec)?
                        .into_iter()
                        .map(|c| Candidate::Region(Cell::Polyhedral(Box::new(c))))
                        .collect(),
                    source_faces: physical_sources(&skin, &self.faces, source_faces),
                });
            }
        }
        // A convex carrier has only one physical patch. Its aliases all move
        // together; coplanar disconnected patches took the skin path above.
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for f in source_faces {
            if let Some(group) = groups.iter_mut().find(|group| {
                self.faces[group[0]]
                    .support
                    .same_definition(&self.faces[*f].support)
            }) {
                group.push(*f);
            } else {
                groups.push(vec![*f]);
            }
        }
        for group in &groups {
            if group.iter().any(|f| selected.contains(f)) {
                selected.extend(group.iter().copied());
            }
        }
        if let Some(cell) = self.annular_cell(leaf, &groups, selected, definition, dec)? {
            return Ok(Candidates {
                components: vec![Candidate::Region(Cell::Annular(Box::new(cell)))],
                source_faces: groups,
            });
        }
        let mut planes = Vec::new();
        let mut cylinder = None;
        let witness = self.leaves[leaf].draft_witness(dec)?;
        for group in &groups {
            let f = &group[0];
            let face = &self.faces[*f];
            let origins: Vec<_> = group
                .iter()
                .map(|f| self.faces[*f].origin.clone())
                .collect();
            let mut surface = face.surface.clone();
            if let Surface::Plane(plane) = &mut surface {
                let side = geometry::sign(&witness.minus(&plane.origin).dot(&plane.normal), dec)?;
                if side == Sign::Zero {
                    return Err(Refusal::new(
                        "draft.unsupported-arrangement",
                        "a source support passes through its interior witness",
                    ));
                }
                if side == Sign::Positive {
                    plane.normal = plane.normal.scaled(&Real::from_i64(-1));
                }
            }
            if selected.contains(f) {
                surface = definition.support(&surface, true, dec)?.surface;
            }
            match surface {
                Surface::Plane(plane) => planes.push(PlanePatch { plane, origins }),
                Surface::Cylinder(_) | Surface::Cone(_) => {
                    let Surface::Cylinder(c) = &face.surface else {
                        return Err(Refusal::new(
                            "draft.unsupported-surface",
                            "selected field support is not an axial cylinder",
                        ));
                    };
                    if cylinder.is_some() {
                        return Err(Refusal::new(
                            "draft.unsupported-arrangement",
                            "mixed curved supports need oriented curved cells",
                        ));
                    }
                    cylinder = Some((c.clone(), origins, selected.contains(f)));
                }
                _ => {
                    return Err(Refusal::new(
                        "draft.unsupported-arrangement",
                        "field draft requires analytic oriented cells",
                    ));
                }
            }
        }
        let cell = if let Some((c, origins, move_side)) = cylinder {
            let caps = planes
                .into_iter()
                .map(|p| (p.plane, p.origins))
                .collect::<Vec<_>>()
                .try_into()
                .map_err(|_| {
                    Refusal::new(
                        "draft.unsupported-round-cap",
                        "axial field draft needs two complete cap planes",
                    )
                })?;
            if move_side {
                Candidate::Region(Cell::Round(Box::new(RoundCell::new(
                    &c, origins, caps, definition, dec,
                )?)))
            } else {
                Candidate::Capped(Box::new(CappedCylinder::new(c, origins, caps, dec)?))
            }
        } else {
            Candidate::Region(Cell::Planar(PlanarCell::new(planes, dec)?.ok_or_else(
                || {
                    Refusal::new(
                        "draft.reversed-solid",
                        "replaced boundary has no full-dimensional material cell",
                    )
                },
            )?))
        };
        Ok(Candidates {
            components: vec![cell],
            source_faces: groups,
        })
    }

    pub(crate) fn draft(
        &mut self,
        source: (
            usize,
            &[usize],
            Option<&catallax_features::draft::arrangement::Skin>,
        ),
        refs: &[FaceRef],
        definition: &DraftDefinition,
        (keep_regions, preserved): (Option<&[catallax_recipe::Point3]>, Option<&Preserved>),
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let (source, source_faces, boundary) = source;
        let actual_faces = boundary.map(|skin| {
            source_faces
                .iter()
                .copied()
                .filter(|f| {
                    skin.patches
                        .iter()
                        .any(|patch| patch.origins.contains(&self.faces[*f].origin))
                })
                .collect::<Vec<_>>()
        });
        let source_faces = actual_faces.as_deref().unwrap_or(source_faces);
        let leaf = match self.instructions[source] {
            Instruction::Leaf(leaf) => Some(leaf),
            _ => None,
        };
        let mut selected = BTreeSet::new();
        let patches = if let Some(skin) = boundary {
            let patches = skin.resolve_faces(refs, dec)?;
            bind_patches(skin, &patches, &self.faces, source_faces, &mut selected);
            Some(patches)
        } else {
            for reference in refs {
                let matches: Vec<_> = source_faces
                    .iter()
                    .copied()
                    .filter(|f| reference.face.matches(&self.faces[*f].origin))
                    .collect();
                if matches.is_empty() {
                    return Err(Refusal::new(
                        "recipe.missing-reference",
                        format!("no input face resolves {reference:?}"),
                    ));
                }
                if !reference.touching.is_empty() {
                    return Err(Refusal::new(
                        "draft.unsupported-reference",
                        "curved field references need exact patch adjacency and distance",
                    ));
                }
                match &reference.near {
                    Some(near) if matches.len() > 1 => {
                        selected.extend(self.nearest_axial(source_faces, &matches, near, dec)?)
                    }
                    _ => selected.extend(matches),
                }
            }
            None
        };
        let Candidates {
            mut components,
            source_faces: mut physical_sources,
        } = self.draft_cells(
            (leaf, boundary.zip(patches.as_ref())),
            source_faces,
            &mut selected,
            definition,
            dec,
        )?;
        if let Some(preserved) = preserved {
            let [Candidate::Region(Cell::Planar(cell))] = components.as_slice() else {
                return Err(crate::draft_source::failure(
                    "preserved source material needs a convex planar envelope",
                ));
            };
            let bores = preserved.bores.iter().flat_map(|bores| &bores.bores);
            match &preserved.boss {
                None => {
                    let bores: Vec<_> = bores.collect();
                    let cylinders: Vec<_> = bores.iter().map(|bore| &bore.cylinder).collect();
                    let planes = cell.supports.iter().map(|support| &support.plane);
                    let caps = bored::certify_bores(planes.clone(), &cylinders, dec)
                        .map_err(crate::draft_source::source_refusal)?;
                    for bore in bores {
                        if let Some(head) = &bore.head {
                            bored::certify_head(
                                planes.clone(),
                                &bore.cylinder,
                                &head.cone,
                                caps,
                                dec,
                            )
                            .map_err(crate::draft_source::source_refusal)?;
                        }
                    }
                }
                Some(boss) => {
                    let bores = bores
                        .map(|bore| Bore {
                            cylinder: bore.cylinder.clone(),
                            origins: bore
                                .walls
                                .iter()
                                .map(|wall| self.faces[*wall].origin.clone())
                                .collect(),
                        })
                        .collect();
                    let arrange = if boss.through {
                        BossedCell::through
                    } else {
                        BossedCell::new
                    };
                    let cell = arrange(cell.clone(), boss.boss.clone(), bores, definition, dec)
                        .map_err(crate::draft_source::source_refusal)?;
                    components = vec![Candidate::Bossed(Box::new(cell))];
                }
            }
        }
        let mut retained = vec![false; components.len()];
        if let Some(points) = keep_regions {
            for p in points {
                let [x, y, z] = p.0.each_ref().map(|n| n.value().clone().into());
                let point = Point3::new(x, y, z);
                let tools = preserved.and_then(|preserved| preserved.bores.as_ref());
                for tool in tools.iter().flat_map(|bores| &bores.tools) {
                    if self.register_location(*tool, &point, dec)? != Location::Outside {
                        return Err(Refusal::new(
                            "draft.invalid-keep-region",
                            "saved region point is on a boundary or outside all candidate material",
                        ));
                    }
                }
                let mut found = false;
                for (i, cell) in components.iter().enumerate() {
                    if cell.contains_interior(&point, dec)? {
                        retained[i] = true;
                        found = true;
                    }
                }
                if !found {
                    return Err(Refusal::new(
                        "draft.invalid-keep-region",
                        "saved region point is on a boundary or outside all candidate material",
                    ));
                }
            }
        } else {
            for (i, cell) in components.iter().enumerate() {
                retained[i] = cell.anchored(definition, dec)?;
            }
        }
        let mut surviving = BTreeSet::new();
        let exact_distance_components =
            components.iter().all(|c| matches!(c, Candidate::Region(_)));
        let mut result = 0;
        for (cell, keep) in components.into_iter().zip(retained) {
            if !keep {
                continue;
            }
            let mut faces = Vec::<Face>::new();
            let mut add = |surface: Surface, origins: &[Origin]| -> Result<(), Refusal> {
                for origin in origins {
                    faces.push(geometry::face(surface.clone(), origin.clone(), dec)?);
                }
                Ok(())
            };
            let planar = match &cell {
                Candidate::Region(Cell::Planar(c)) => Some(c),
                Candidate::Bossed(c) => Some(&c.outer),
                _ => None,
            };
            for (s, cycle) in planar.iter().flat_map(|c| c.supports.iter().zip(&c.cycles)) {
                if !cycle.is_empty() {
                    add(Surface::Plane(s.plane.clone()), &s.origins)?;
                }
            }
            match &cell {
                Candidate::Region(Cell::Annular(c)) => {
                    add(c.outer.surface.clone(), &c.outer.origins)?;
                    add(c.inner.surface.clone(), &c.inner.origins)?;
                    for (p, origins) in c.outer.caps.iter().flatten() {
                        add(Surface::Plane(p.clone()), origins)?;
                    }
                }
                Candidate::Region(Cell::Planar(_)) | Candidate::Bossed(_) => (),
                Candidate::Region(Cell::Round(c)) => {
                    add(c.surface.clone(), &c.origins)?;
                    for (p, origins) in c.caps.iter().flatten() {
                        add(Surface::Plane(p.clone()), origins)?;
                    }
                }
                Candidate::Region(Cell::Polyhedral(c)) => {
                    for face in &c.faces {
                        add(Surface::Plane(c.plane(face)), &face.origins)?;
                    }
                }
                Candidate::Capped(c) => {
                    add(Surface::Cylinder(c.cylinder.clone()), &c.origins)?;
                    for (plane, origins) in &c.caps {
                        add(Surface::Plane(plane.clone()), origins)?;
                    }
                }
            }
            surviving.extend(faces.iter().map(|f| f.origin.clone()));
            let leaf = match cell {
                Candidate::Region(cell) => geometry::Leaf::Draft(Region::new(cell, dec)?),
                Candidate::Capped(cell) => {
                    geometry::Leaf::Cell(geometry::Cell::from_capped(cell, dec)?)
                }
                // The tool stays in the tape and adds the boss back.
                Candidate::Bossed(cell) => {
                    geometry::Leaf::Draft(Region::new(Cell::Planar(cell.outer), dec)?)
                }
            };
            let register = self.leaf(leaf, faces, dec)?;
            result = if result == 0 {
                register
            } else {
                let union =
                    self.boolean(catallax_recipe::BooleanKind::Union, result, register, dec)?;
                if exact_distance_components {
                    self.separated_unions.insert(union);
                }
                union
            };
        }
        let by: Vec<_> = selected
            .iter()
            .map(|f| self.faces[*f].origin.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        // Walls survive with any retained material. Each is one physical face.
        if result == 0 {
            if let Some(bores) = preserved.and_then(|preserved| preserved.bores.as_ref()) {
                physical_sources.extend(bores.bores.iter().map(|bore| bore.walls.clone()));
                let heads = bores.bores.iter().flat_map(|bore| &bore.head);
                physical_sources.extend(heads.map(|head| head.walls.clone()));
            }
            if let Some(boss) = preserved.and_then(|preserved| preserved.boss.as_ref()) {
                physical_sources.extend([boss.wall.clone(), boss.top.clone()]);
            }
        }
        for group in physical_sources {
            if !group
                .iter()
                .all(|f| !surviving.contains(&self.faces[*f].origin))
            {
                continue;
            }
            let source = group
                .iter()
                .map(|f| self.faces[*f].origin.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            self.consumed.push(catallax_features::ConsumedFace {
                step: step.clone(),
                source,
                by: by.clone(),
            });
        }
        Ok(result)
    }
}
