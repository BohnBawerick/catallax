//! Analytic draft cells evaluated independently of B-rep topology.

use crate::{Face, Instruction, Location, Plan, geometry, metric::PRECISION};
use catallax_features::draft::{
    DraftDefinition,
    cell::{PlanarCell, PlanePatch},
    round::RoundCell,
};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Origin, StepId};
use catallax_shapes::{Isometry, Plane, Point3, RPoint3, Surface, Vector3};
use std::collections::BTreeSet;

#[derive(Clone, Debug)]
pub(crate) enum Cell {
    Planar(PlanarCell),
    Round(Box<RoundCell>),
}

#[derive(Clone, Debug)]
enum Fast {
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
    fn new(cell: Cell, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut edges = BTreeSet::new();
        if let Cell::Planar(cell) = &cell {
            for cycle in &cell.cycles {
                for i in 0..cycle.len() {
                    let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                    edges.insert((a.min(b), a.max(b)));
                }
            }
        }
        let fast = Self::fast_geometry(&cell);
        let mut result = Self {
            cell,
            constraints: Vec::new(),
            edges: edges.into_iter().collect(),
            fast,
            bounds: None,
        };
        result.refresh_constraints();
        result.bounds = Some(result.display_box(dec)?.map(|b| [b.lo(), b.hi()]));
        Ok(result)
    }
    fn fast_geometry(cell: &Cell) -> Fast {
        match cell {
            Cell::Planar(cell) => Fast::Planar {
                vertices: cell.vertices.iter().map(point).collect(),
                planes: cell
                    .supports
                    .iter()
                    .map(|s| (point(&s.plane.origin), vector(&s.plane.normal)))
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
    fn refresh_constraints(&mut self) {
        self.constraints.clear();
        match &self.cell {
            Cell::Planar(cell) => self.constraints.extend(
                cell.supports
                    .iter()
                    .map(|s| geometry::constraint(Surface::Plane(s.plane.clone()))),
            ),
            Cell::Round(cell) => {
                self.constraints
                    .push(geometry::constraint(cell.surface.clone()));
                for (plane, _) in cell.caps.iter().flatten() {
                    self.constraints
                        .push(geometry::constraint(Surface::Plane(plane.clone())));
                }
                if let Surface::Cone(cone) = &cell.surface {
                    self.constraints
                        .push(geometry::constraint(Surface::Plane(Plane {
                            origin: cone.apex.clone(),
                            normal: cone.axis.scaled(&Real::from_i64(-1)),
                            ref_dir: cone.ref_dir.clone(),
                        })));
                }
            }
        }
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
            Cell::Planar(cell) => {
                for p in &cell.vertices {
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
    pub(crate) fn transform(&mut self, map: &Isometry) {
        match &mut self.cell {
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
        self.refresh_constraints();
    }
    pub(crate) fn location(
        &self,
        p: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
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
            Cell::Planar(cell) => {
                for (support, cycle) in cell.supports.iter().zip(&cell.cycles) {
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
                    for support in &cell.supports {
                        if geometry::sign(
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

impl Plan {
    pub(crate) fn draft(
        &mut self,
        source: (usize, &[usize]),
        refs: &[FaceRef],
        definition: &DraftDefinition,
        keep_regions: Option<&[catallax_recipe::Point3]>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let (source, source_faces) = source;
        let Instruction::Leaf(leaf) = self.instructions[source] else {
            if source == 0 && refs.is_empty() {
                return Ok(0);
            }
            return Err(Refusal::new(
                "draft.unsupported-arrangement",
                "field draft needs an oriented multi-cell arrangement for Boolean inputs",
            ));
        };
        let mut selected = BTreeSet::new();
        for reference in refs {
            let matches: Vec<_> = source_faces
                .iter()
                .copied()
                .filter(|f| reference.face.matches(&self.faces[*f].origin))
                .collect();
            if matches.is_empty() {
                return Err(Refusal::new(
                    "draft.lost-face",
                    format!("no input face resolves {reference:?}"),
                ));
            }
            if !reference.touching.is_empty() || reference.near.is_some() {
                return Err(Refusal::new(
                    "draft.unsupported-reference",
                    "filtered field draft references need exact adjacency and trimmed-face distances",
                ));
            }
            selected.extend(matches);
        }
        // A convex leaf can publish several origins for one merged support,
        // such as adjacent collinear sketch segments. Replace the whole patch.
        for face in source_faces {
            if selected.iter().any(|s| {
                self.faces[*face]
                    .support
                    .same_definition(&self.faces[*s].support)
            }) {
                selected.insert(*face);
            }
        }
        let mut planes = Vec::new();
        let mut cylinder = None;
        let witness = self.leaves[leaf].draft_witness()?;
        for f in source_faces {
            let face = &self.faces[*f];
            let mut surface = face.surface.clone();
            if let Surface::Plane(plane) = &mut surface {
                let side = geometry::sign(&witness.minus(&plane.origin).dot(&plane.normal), dec)?;
                if side == Sign::Zero {
                    return Err(Refusal::new(
                        "draft.unsupported-arrangement",
                        "a source support passes through the proposed interior witness",
                    ));
                }
                if side == Sign::Positive {
                    plane.normal = plane.normal.scaled(&Real::from_i64(-1));
                }
            }
            if selected.contains(f) {
                // The plane above is already outward oriented.
                let replacement = definition.support(&surface, true, dec)?;
                surface = replacement.surface;
            }
            match surface {
                Surface::Plane(plane) => planes.push(PlanePatch {
                    plane,
                    origins: vec![face.origin.clone()],
                }),
                Surface::Cylinder(_) | Surface::Cone(_) => {
                    let original = &face.surface;
                    let Surface::Cylinder(c) = original else {
                        return Err(Refusal::new(
                            "draft.unsupported-surface",
                            "selected field support is not an axial cylinder",
                        ));
                    };
                    if !selected.contains(f) || cylinder.is_some() {
                        return Err(Refusal::new(
                            "draft.unsupported-arrangement",
                            "mixed curved support arrangements need additional cells",
                        ));
                    }
                    cylinder = Some((c.clone(), vec![face.origin.clone()]));
                }
                _ => {
                    return Err(Refusal::new(
                        "draft.unsupported-arrangement",
                        "field draft supports require analytic oriented cells",
                    ));
                }
            }
        }
        let cell = if let Some((c, origins)) = cylinder {
            let caps: Vec<_> = planes.into_iter().map(|p| (p.plane, p.origins)).collect();
            let caps = caps.try_into().map_err(|_| {
                Refusal::new(
                    "draft.unsupported-round-cap",
                    "axial field draft needs two complete cap planes",
                )
            })?;
            Cell::Round(Box::new(RoundCell::new(
                &c, origins, caps, definition, dec,
            )?))
        } else {
            self.leaves[leaf].certify_planar_convexity(dec)?;
            Cell::Planar(PlanarCell::new(planes, dec)?.ok_or_else(|| {
                Refusal::new(
                    "draft.reversed-solid",
                    "replaced boundary has no full-dimensional material cell",
                )
            })?)
        };
        let retained = if let Some(points) = keep_regions {
            for p in points {
                let [x, y, z] = p.0.each_ref().map(|n| n.value().clone().into());
                let p = Point3::new(x, y, z);
                let inside = match &cell {
                    Cell::Planar(c) => c.contains_interior(&p, dec)?,
                    Cell::Round(c) => c.contains_interior(&p, dec)?,
                };
                if !inside {
                    return Err(Refusal::new(
                        "draft.invalid-keep-region",
                        "saved region point is on a boundary or outside all candidate material",
                    ));
                }
            }
            !points.is_empty()
        } else {
            match &cell {
                Cell::Planar(c) => c.anchored(&definition.neutral, dec)?,
                Cell::Round(c) => c.anchored(&definition.neutral, dec)?,
            }
        };
        let mut surviving = BTreeSet::new();
        if retained {
            match &cell {
                Cell::Planar(cell) => {
                    for (support, cycle) in cell.supports.iter().zip(&cell.cycles) {
                        if !cycle.is_empty() {
                            surviving.extend(support.origins.iter().cloned());
                        }
                    }
                }
                Cell::Round(cell) => {
                    surviving.extend(cell.origins.iter().cloned());
                    for (_, origins) in cell.caps.iter().flatten() {
                        surviving.extend(origins.iter().cloned());
                    }
                }
            }
        }
        let by: Vec<_> = selected
            .iter()
            .map(|f| self.faces[*f].origin.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let mut removed: Vec<_> = source_faces
            .iter()
            .filter_map(|f| {
                let origin = &self.faces[*f].origin;
                (!surviving.contains(origin)).then(|| vec![origin.clone()])
            })
            .collect();
        removed.sort();
        removed.dedup();
        for source in removed {
            self.consumed.push(catallax_features::ConsumedFace {
                step: step.clone(),
                source,
                by: by.clone(),
            });
        }
        if !retained {
            return Ok(0);
        }
        let mut faces = Vec::<Face>::new();
        let mut add = |surface: Surface, origins: &[Origin]| -> Result<(), Refusal> {
            for origin in origins {
                faces.push(geometry::face(surface.clone(), origin.clone(), dec)?);
            }
            Ok(())
        };
        match &cell {
            Cell::Planar(cell) => {
                for (s, cycle) in cell.supports.iter().zip(&cell.cycles) {
                    if !cycle.is_empty() {
                        add(Surface::Plane(s.plane.clone()), &s.origins)?;
                    }
                }
            }
            Cell::Round(cell) => {
                add(cell.surface.clone(), &cell.origins)?;
                for (p, origins) in cell.caps.iter().flatten() {
                    add(Surface::Plane(p.clone()), origins)?;
                }
            }
        }
        let region = Region::new(cell, dec)?;
        self.leaf(geometry::Leaf::Draft(region), faces, dec)
    }
}
