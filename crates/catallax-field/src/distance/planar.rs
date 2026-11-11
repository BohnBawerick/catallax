//! Field-native planar section arrangements. Source and historical Boolean
//! germs decide the closed trim, never a B-rep, display mesh or residual sample.

use super::NamedPlanarRegion;
use crate::{
    Direction, FaceId, FaceQuery, FaceState, Instruction, Location, Plan, geometry,
    ownership::Domain,
};
use catallax_intersect::{
    Component, CurveForm, CurveSurface, SurfacePair, arc_point, curve_form,
    intersect_curve_surface, intersect_surfaces, same_curve, tangent_at,
};
use catallax_numbers::{Decider, Outcome, Rational, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{
    Circle, CurveId, GeomStore, Plane, Point3, PointDef, PointId, RPoint3, Surface, Vector3,
    distance::{PlanarBoundary, PlanarEdge},
};
use std::{
    collections::BTreeMap,
    sync::{Arc, OnceLock},
};

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.unsupported-planar-trim", detail)
}

/// A positive infinitesimal ray represented by a certified rational step. Every
/// source polynomial keeps its first nonzero Taylor sign throughout this step.
pub(super) fn ray_sample<'a>(
    constraints: impl Iterator<Item = &'a geometry::Constraint>,
    point: &Point3,
    direction: &Vector3,
    dec: &mut Decider,
) -> Result<Point3, Refusal> {
    let xyz = [point.x.clone(), point.y.clone(), point.z.clone()];
    let mut epsilon = Rational::one();
    for constraint in constraints {
        let coefficients = geometry::point_ray_coefficients(constraint, &xyz, direction);
        for (i, coefficient) in coefficients.iter().enumerate() {
            let sign = geometry::sign(coefficient, dec)?;
            if sign == Sign::Zero {
                continue;
            }
            let magnitude = if sign == Sign::Negative {
                -coefficient
            } else {
                coefficient.clone()
            };
            let tail = coefficients[i + 1..].iter().fold(Real::one(), |sum, c| {
                sum + c
                    .rational_ref()
                    .map(|q| Real::from(q.abs()))
                    .unwrap_or_else(|| c.square() + Real::one())
            });
            let bound = &magnitude / &(&magnitude + &tail);
            let mut bits = 128u32;
            loop {
                let lower = dec.enclose_budgeted_at(&bound, bits)?.lo();
                if lower > Rational::zero() {
                    epsilon = epsilon.min(lower.mul_pow2(-1));
                    break;
                }
                bits = bits.checked_mul(2).ok_or_else(|| {
                    Refusal::new(
                        "numbers.refinement-budget",
                        "planar germ precision overflow",
                    )
                })?;
            }
            break;
        }
    }
    Ok(point.plus(&direction.scaled(&Real::from(epsilon))))
}

struct Scene<'a> {
    plan: &'a Plan,
    active: Vec<bool>,
    leaves: Vec<usize>,
}

impl<'a> Scene<'a> {
    fn new(plan: &'a Plan) -> Self {
        let mut active = vec![false; plan.instructions.len()];
        active[plan.root] = true;
        for domain in &plan.domains {
            if let Domain::Boundary { node, .. } = domain {
                active[*node] = true;
            }
        }
        for i in (0..active.len()).rev() {
            if active[i]
                && let Instruction::Boolean { a, b, .. } = plan.instructions[i]
            {
                active[a] = true;
                active[b] = true;
            }
        }
        let leaves = plan
            .instructions
            .iter()
            .enumerate()
            .filter_map(|(i, node)| match node {
                Instruction::Leaf(leaf) if active[i] => Some(*leaf),
                _ => None,
            })
            .collect();
        Self {
            plan,
            active,
            leaves,
        }
    }
}

impl Scene<'_> {
    fn germ(
        &self,
        point: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<bool>, Refusal> {
        let mut values = vec![false; self.active.len()];
        for (node, instruction) in self
            .plan
            .instructions
            .iter()
            .enumerate()
            .filter(|(i, _)| self.active[*i])
        {
            values[node] = match *instruction {
                Instruction::Empty => false,
                Instruction::Leaf(leaf) => {
                    match super::ray(&self.plan.leaves[leaf], point, direction, dec)? {
                        Location::Inside => true,
                        Location::Outside => false,
                        _ => return Err(bad("a planar section has no strict source germ")),
                    }
                }
                Instruction::Boolean { kind, a, b } => match kind {
                    catallax_recipe::BooleanKind::Union => values[a] || values[b],
                    catallax_recipe::BooleanKind::Intersect => values[a] && values[b],
                    catallax_recipe::BooleanKind::Subtract => values[a] && !values[b],
                },
            };
        }
        Ok(values)
    }

    fn memberships(
        &self,
        faces: &[usize],
        point: &Point3,
        normal: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<bool>, Refusal> {
        let plus = self.germ(point, normal, dec)?;
        let minus = self.germ(point, &normal.scaled(&Real::from_i64(-1)), dec)?;
        if plus[self.plan.root] == minus[self.plan.root] {
            return Ok(vec![false; faces.len()]);
        }
        let mut domains = Vec::with_capacity(self.plan.domains.len());
        for domain in &self.plan.domains {
            domains.push(match *domain {
                Domain::Source => true,
                Domain::Absent => false,
                Domain::Boundary { node, parent } => domains[parent] && plus[node] != minus[node],
                Domain::Either { a, b } => domains[a] || domains[b],
            });
        }
        let mut result = Vec::with_capacity(faces.len());
        for &id in faces {
            let face = &self.plan.faces[id];
            let register = self.plan.leaf_registers[face.leaf].expect("live source leaf");
            let mut member = domains[face.domain] && plus[register] != minus[register];
            if member {
                member = match &self.plan.leaves[face.leaf] {
                    geometry::Leaf::Feature(feature) => {
                        feature.on_face_point(self.plan.definition_role(id), point, dec)?
                    }
                    geometry::Leaf::Draft(region) => {
                        region.on_face_point(&face.origin, point, dec)?
                    }
                    _ => true,
                };
            }
            result.push(member);
        }
        Ok(result)
    }

    fn sample(
        &self,
        point: &Point3,
        direction: &Vector3,
        trims: &[geometry::Constraint],
        dec: &mut Decider,
    ) -> Result<Point3, Refusal> {
        ray_sample(
            self.leaves
                .iter()
                .flat_map(|i| self.plan.leaves[*i].constraints())
                .chain(trims),
            point,
            direction,
            dec,
        )
    }
}

fn intern(
    geom: &mut GeomStore,
    points: &mut Vec<PointId>,
    point: Point3,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for id in points.iter() {
        if geometry::sign(&geom.position(*id).minus(&point).norm_squared(), dec)? == Sign::Zero {
            return Ok(());
        }
    }
    points.push(geom.add_point(PointDef::Explicit(point)));
    Ok(())
}

fn sort_points(
    geom: &GeomStore,
    points: &mut [PointId],
    form: &CurveForm,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    for i in 1..points.len() {
        let mut j = i;
        while j > 0 {
            let a = geom.position(points[j - 1]);
            let b = geom.position(points[j]);
            let after = match form {
                CurveForm::Line { dir, .. } => {
                    geometry::sign(&a.minus(&b).dot(dir), dec)? == Sign::Positive
                }
                CurveForm::Circle {
                    center,
                    normal,
                    ref_dir,
                    ..
                } => {
                    let a = a.minus(center);
                    let b = b.minus(center);
                    let half = |p: &Vector3, dec: &mut Decider| -> Result<bool, Refusal> {
                        let y = geometry::sign(&ref_dir.cross(p).dot(normal), dec)?;
                        Ok(y == Sign::Negative
                            || (y == Sign::Zero
                                && geometry::sign(&ref_dir.dot(p), dec)? == Sign::Negative))
                    };
                    let (ah, bh) = (half(&a, dec)?, half(&b, dec)?);
                    if ah != bh {
                        ah
                    } else {
                        geometry::sign(&a.cross(&b).dot(normal), dec)? == Sign::Negative
                    }
                }
                _ => {
                    return Err(bad(
                        "only line and circular-arc trims are in this distance domain",
                    ));
                }
            };
            if !after {
                break;
            }
            points.swap(j - 1, j);
            j -= 1;
        }
    }
    Ok(())
}

fn line_bounds(
    origin: &Point3,
    dir: &Vector3,
    radius: &Rational,
    dec: &mut Decider,
) -> Result<Option<[Point3; 2]>, Refusal> {
    let mut interval: Option<[Real; 2]> = None;
    for (x, d) in [
        (&origin.x, &dir.x),
        (&origin.y, &dir.y),
        (&origin.z, &dir.z),
    ] {
        let sign = geometry::sign(d, dec)?;
        if sign == Sign::Zero {
            if geometry::sign(&(x - &Real::from(radius.clone())), dec)? == Sign::Positive
                || geometry::sign(&(x + &Real::from(radius.clone())), dec)? == Sign::Negative
            {
                return Ok(None);
            }
            continue;
        }
        let mut range = [
            &(&Real::from(-radius) - x) / d,
            &(&Real::from(radius.clone()) - x) / d,
        ];
        if sign == Sign::Negative {
            range.swap(0, 1);
        }
        if let Some(current) = &mut interval {
            if geometry::sign(&(&range[0] - &current[0]), dec)? == Sign::Positive {
                current[0] = range[0].clone();
            }
            if geometry::sign(&(&range[1] - &current[1]), dec)? == Sign::Negative {
                current[1] = range[1].clone();
            }
        } else {
            interval = Some(range);
        }
    }
    let Some([a, b]) = interval else {
        return Err(bad("a planar section line has zero direction"));
    };
    if geometry::sign(&(&b - &a), dec)? != Sign::Positive {
        return Ok(None);
    }
    Ok(Some([
        origin.plus(&dir.scaled(&a)),
        origin.plus(&dir.scaled(&b)),
    ]))
}

fn append_named_trims(
    cuts: &mut Vec<geometry::Constraint>,
    own: &Plane,
    patch: &Plane,
    vertices: &[Point3],
    loops: &[Vec<usize>],
    dec: &mut Decider,
) -> Result<u8, Refusal> {
    if geometry::sign(&own.normal.cross(&patch.normal).norm_squared(), dec)? != Sign::Zero
        || geometry::sign(&patch.origin.minus(&own.origin).dot(&own.normal), dec)? != Sign::Zero
    {
        return Ok(0);
    }
    let orientation = match geometry::sign(&patch.normal.dot(&own.normal), dec)? {
        Sign::Positive => 1,
        Sign::Negative => 2,
        Sign::Zero => return Err(bad("a named planar source has no material-side normal")),
    };
    for cycle in loops {
        for i in 0..cycle.len() {
            let a = &vertices[cycle[i]];
            let edge = vertices[cycle[(i + 1) % cycle.len()]].minus(a);
            cuts.push(geometry::constraint(Surface::Plane(Plane {
                origin: a.clone(),
                normal: edge.cross(&own.normal),
                ref_dir: edge,
            }))?);
        }
    }
    Ok(orientation)
}

struct SectionGroup {
    faces: Vec<usize>,
    boundaries: OnceLock<Vec<Option<Arc<PlanarBoundary>>>>,
}

/// Reuse source germs and support groups, constructing only requested sections.
pub(super) struct Sections<'a> {
    scene: Scene<'a>,
    radius: Rational,
    groups: Vec<SectionGroup>,
    group_of: Vec<Option<(usize, usize)>>,
}

impl Plan {
    fn named_planar_region(
        &self,
        id: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<NamedPlanarRegion>, Refusal> {
        let face = &self.faces[id.0];
        let Surface::Plane(own) = &face.surface else {
            return Ok(None);
        };
        let mut trims = Vec::new();
        let mut source_orientation = 0;
        match &self.leaves[face.leaf] {
            geometry::Leaf::Feature(feature) => {
                let origins: Vec<_> = self
                    .faces
                    .iter()
                    .enumerate()
                    .filter(|(_, other)| other.leaf == face.leaf)
                    .map(|(index, other)| Origin {
                        step: other.origin.step.clone(),
                        role: self.definition_role(index).clone(),
                    })
                    .collect();
                let defining = Origin {
                    step: face.origin.step.clone(),
                    role: self.definition_role(id.0).clone(),
                };
                let names = origins.iter().collect::<Vec<_>>();
                let skin = match feature.planar_skin(&names, None, dec) {
                    Ok(skin) => skin,
                    Err(error)
                        if matches!(
                            error.code,
                            "draft.unsupported-arrangement" | "draft.missing-source-origin"
                        ) =>
                    {
                        return Ok(None);
                    }
                    Err(error) => return Err(error),
                };
                for patch in &skin.patches {
                    if patch.origins.contains(&defining) {
                        source_orientation |= append_named_trims(
                            &mut trims,
                            own,
                            &patch.plane,
                            &skin.vertices,
                            &patch.loops,
                            dec,
                        )?;
                    }
                }
            }
            geometry::Leaf::Draft(region) => {
                let crate::draft::Cell::Polyhedral(poly) = &region.cell else {
                    return Ok(None);
                };
                let Some(source) = region.source_origin(&face.origin) else {
                    return Ok(None);
                };
                for patch in &poly.faces {
                    if patch.origins.contains(source) {
                        source_orientation |= append_named_trims(
                            &mut trims,
                            own,
                            &poly.plane(patch),
                            &poly.vertices,
                            &patch.tiles,
                            dec,
                        )?;
                    }
                }
            }
            _ => return Ok(None),
        }
        if source_orientation == 0 {
            return Ok(None);
        }
        let scene = Scene::new(self);
        if !scene.leaves.iter().all(|leaf| match &self.leaves[*leaf] {
            geometry::Leaf::Primitive(_) => true,
            geometry::Leaf::Cell(cell) => cell.capped_chart().is_none(),
            geometry::Leaf::Feature(feature) => feature.complete_planar_carriers(),
            geometry::Leaf::Draft(region) => {
                matches!(&region.cell, crate::draft::Cell::Polyhedral(_))
            }
            geometry::Leaf::Lattice(_) | geometry::Leaf::Gyroid(_) => false,
        }) {
            return Ok(None);
        }
        let radius = match self.node_radius(self.root, dec) {
            Ok(radius) => radius + Rational::one(),
            Err(error) if error.code == "field.unbounded-target" => return Ok(None),
            Err(error) => return Err(error),
        };
        let boundaries = match scene.boundaries(&[id.0], &radius, &trims, dec) {
            Ok(boundaries) => boundaries,
            Err(error) if error.code == "distance.unsupported-planar-trim" => return Ok(None),
            Err(error) => return Err(error),
        };
        Ok(boundaries
            .into_iter()
            .next()
            .flatten()
            .map(|boundary| NamedPlanarRegion {
                boundary,
                source_orientation,
            }))
    }

    pub(super) fn named_planar_query(
        &self,
        id: FaceId,
        point: &RPoint3,
        dec: &mut Decider,
    ) -> Result<Option<FaceQuery>, Refusal> {
        let face = &self.faces[id.0];
        let Surface::Plane(plane) = &face.surface else {
            return Ok(None);
        };
        match &self.leaves[face.leaf] {
            geometry::Leaf::Feature(feature) if feature.complete_planar_carriers() => (),
            geometry::Leaf::Draft(region)
                if matches!(&region.cell, crate::draft::Cell::Polyhedral(_)) => {}
            _ => return Ok(None),
        }
        let regions = self.planar_memberships.0.get_or_init(Default::default);
        let region = {
            let mut regions = regions
                .lock()
                .expect("named planar membership initialization panicked");
            if let Some(region) = regions.get(&id.0) {
                region.clone()
            } else {
                // Refusals leave the cache retryable under a later proof policy.
                let region = self.named_planar_region(id, dec)?;
                regions.insert(id.0, region.clone());
                region
            }
        };
        let Some(region) = region else {
            return Ok(None);
        };
        if !region.boundary.contains(&point.to_real(), dec)? {
            return Ok(None);
        }
        let orientation = match self.leaf_polarities[face.leaf] {
            1 => region.source_orientation,
            2 => ((region.source_orientation & 1) << 1) | ((region.source_orientation & 2) >> 1),
            _ => 3,
        };
        let direction = match orientation {
            1 => crate::unit_direction(plane.normal.clone(), dec)?,
            2 => crate::unit_direction(plane.normal.scaled(&Real::from_i64(-1)), dec)?,
            _ => Direction::Ambiguous,
        };
        Ok(Some(FaceQuery {
            state: FaceState::OnFace,
            direction,
        }))
    }

    pub(super) fn prepare_planar_boundaries(
        &self,
        dec: &mut Decider,
    ) -> Result<Sections<'_>, Refusal> {
        let scene = Scene::new(self);
        let radius = self.node_radius(self.root, dec)? + Rational::one();
        let mut groups: Vec<SectionGroup> = Vec::new();
        let mut group_of = vec![None; self.faces.len()];
        for (face, definition) in self.faces.iter().enumerate() {
            if !definition.in_result
                || self.leaf_registers[definition.leaf].is_none()
                || !matches!(definition.surface, Surface::Plane(_))
            {
                continue;
            }
            let mut group = None;
            for (i, members) in groups.iter().enumerate() {
                if definition
                    .support
                    .equals(&self.faces[members.faces[0]].support, dec)?
                {
                    group = Some(i);
                    break;
                }
            }
            let i = if let Some(i) = group {
                i
            } else {
                groups.push(SectionGroup {
                    faces: Vec::with_capacity(1),
                    boundaries: OnceLock::new(),
                });
                groups.len() - 1
            };
            group_of[face] = Some((i, groups[i].faces.len()));
            groups[i].faces.push(face);
        }
        Ok(Sections {
            scene,
            radius,
            groups,
            group_of,
        })
    }
}

impl Sections<'_> {
    pub(super) fn get(&self, face: usize) -> Option<&Option<Arc<PlanarBoundary>>> {
        let (group, member) = self.group_of[face]?;
        Some(&self.groups[group].boundaries.get()?[member])
    }

    // Prepared's initialization mutex serializes writes to the section caches.
    pub(super) fn prepare(&self, face: usize, dec: &mut Decider) -> Result<(), Refusal> {
        let (group, _) = self.group_of[face].expect("distance face is a live plane");
        let group = &self.groups[group];
        if group.boundaries.get().is_none() {
            let boundaries = match self.scene.boundaries(&group.faces, &self.radius, &[], dec) {
                Ok(boundaries) => boundaries,
                Err(error) if error.code == "distance.unsupported-planar-trim" => {
                    vec![None; group.faces.len()]
                }
                Err(error) => return Err(error),
            };
            group
                .boundaries
                .set(boundaries)
                .expect("a planar section initializes under one lock");
        }
        Ok(())
    }
}

impl Scene<'_> {
    fn boundaries(
        &self,
        faces: &[usize],
        radius: &Rational,
        trims: &[geometry::Constraint],
        dec: &mut Decider,
    ) -> Result<Vec<Option<Arc<PlanarBoundary>>>, Refusal> {
        let Surface::Plane(plane) = &self.plan.faces[faces[0]].surface else {
            unreachable!()
        };
        let mut geom = GeomStore::new();
        let own = geom.add_surface(Surface::Plane(plane.clone()));
        let surfaces: Vec<_> = self
            .leaves
            .iter()
            .flat_map(|leaf| self.plan.leaves[*leaf].constraints())
            .chain(trims)
            .map(|c| geom.add_surface(c.surface().clone()))
            .collect();
        let mut curves: Vec<(CurveId, CurveForm)> = Vec::new();
        for &surface in &surfaces {
            let answer = intersect_surfaces(&mut geom, own, surface, dec)?;
            if !answer.marks.is_empty() {
                return Err(bad("a planar section has marked intersections"));
            }
            let SurfacePair::Meets(components) = answer.value else {
                continue;
            };
            for component in components {
                let Component::Curve { curve, .. } = component else {
                    continue;
                };
                let form = curve_form(&geom, curve, dec)?;
                if !matches!(form, CurveForm::Line { .. } | CurveForm::Circle { .. }) {
                    return Err(bad(
                        "a planar face needs conic or spline trim metrics outside this slice",
                    ));
                }
                let mut duplicate = false;
                for (other, _) in &curves {
                    match same_curve(&geom, curve, *other, dec) {
                        Outcome::ProvenCoincident(_) => {
                            duplicate = true;
                            break;
                        }
                        Outcome::Proven(_) => (),
                        Outcome::Refused(error) => return Err(error),
                        Outcome::NotExact(n) => {
                            return Err(
                                bad("a planar section curve identity is unresolved").keeping(n)
                            );
                        }
                    }
                }
                if !duplicate {
                    curves.push((curve, form));
                }
            }
        }
        let mut vertices = Vec::new();
        let mut edges = Vec::new();
        let mut selected = vec![Vec::new(); faces.len()];
        for (curve, form) in curves {
            let mut points = Vec::new();
            let closed = matches!(form, CurveForm::Circle { .. });
            let bound = match &form {
                CurveForm::Line { origin, dir } => {
                    let Some(ends) = line_bounds(origin, dir, radius, dec)? else {
                        continue;
                    };
                    for point in &ends {
                        intern(&mut geom, &mut points, point.clone(), dec)?;
                    }
                    Some(ends)
                }
                CurveForm::Circle {
                    center,
                    normal: _,
                    radius_sq,
                    ref_dir,
                } => {
                    let offset = ref_dir.scaled(&(radius_sq / &ref_dir.norm_squared()).sqrt());
                    intern(&mut geom, &mut points, center.plus(&offset), dec)?;
                    intern(
                        &mut geom,
                        &mut points,
                        center.plus(&offset.scaled(&Real::from_i64(-1))),
                        dec,
                    )?;
                    None
                }
                _ => unreachable!(),
            };
            for &surface in &surfaces {
                let answer = intersect_curve_surface(&mut geom, curve, surface, dec)?;
                if !answer.marks.is_empty() {
                    return Err(bad("a planar trim split has marked intersections"));
                }
                if let CurveSurface::Meets(hits) = answer.value {
                    for hit in hits {
                        let point = geom.position(hit.point);
                        if let (Some([a, b]), CurveForm::Line { dir, .. }) = (&bound, &form)
                            && (geometry::sign(&point.minus(a).dot(dir), dec)? == Sign::Negative
                                || geometry::sign(&b.minus(&point).dot(dir), dec)?
                                    == Sign::Negative)
                        {
                            continue;
                        }
                        intern(&mut geom, &mut points, point, dec)?;
                    }
                }
            }
            sort_points(&geom, &mut points, &form, dec)?;
            let count = if closed {
                points.len()
            } else {
                points.len() - 1
            };
            for i in 0..count {
                let a = points[i];
                let b = points[(i + 1) % points.len()];
                let answer = arc_point(&geom, curve, Some(a), Some(b), dec)?;
                if !answer.marks.is_empty() {
                    return Err(bad("a planar trim interval has no exact interior point"));
                }
                let middle = answer.value;
                let middle_id = geom.add_point(PointDef::Explicit(middle.clone()));
                let tangent = tangent_at(&geom, curve, middle_id, dec)?[0].clone();
                let across = plane.normal.cross(&tangent);
                let left = self.sample(&middle, &across, trims, dec)?;
                let right =
                    self.sample(&middle, &across.scaled(&Real::from_i64(-1)), trims, dec)?;
                let left = self.memberships(faces, &left, &plane.normal, dec)?;
                let right = self.memberships(faces, &right, &plane.normal, dec)?;
                if !left.iter().zip(&right).any(|(a, b)| a != b) {
                    continue;
                }
                let pair = [a, b].map(|id| {
                    let index = vertices.len();
                    vertices.push(geom.position(id));
                    index
                });
                let edge = match &form {
                    CurveForm::Line { .. } => PlanarEdge::Segment(pair),
                    CurveForm::Circle {
                        center,
                        normal,
                        radius_sq,
                        ref_dir,
                    } => PlanarEdge::Arc {
                        circle: Circle {
                            center: center.clone(),
                            normal: normal.clone(),
                            radius_sq: radius_sq.clone(),
                            ref_dir: ref_dir.clone(),
                        },
                        ends: Some(pair),
                    },
                    _ => unreachable!(),
                };
                let index = edges.len();
                edges.push(edge);
                for ((a, b), selected) in left.into_iter().zip(right).zip(&mut selected) {
                    if a != b {
                        selected.push(index);
                    }
                }
            }
        }
        let mut shared = BTreeMap::new();
        Ok(selected
            .into_iter()
            .map(|selected| {
                if selected.is_empty() {
                    return None;
                }
                Some(
                    shared
                        .entry(selected.clone())
                        .or_insert_with(|| {
                            Arc::new(PlanarBoundary {
                                plane: plane.clone(),
                                vertices: vertices.clone(),
                                edges: selected.into_iter().map(|i| edges[i].clone()).collect(),
                            })
                        })
                        .clone(),
                )
            })
            .collect())
    }
}
