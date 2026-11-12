//! Trimmed distances from the field's own recipe predicates, without B-rep input.

use crate::{FaceId, Instruction, Location, Plan, geometry};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{
    Point3, Vector3,
    distance::{self, Bounds, Domain, Membership, Options},
};
use std::{
    collections::BTreeMap,
    sync::{Arc, Mutex, OnceLock},
};

mod contact;
mod nearest;
pub(crate) mod patch;
mod planar;
pub use nearest::{Ambiguity, Candidate, NearestOptions, Owner, Ownership};

#[derive(Clone, Debug)]
struct NamedPlanarRegion {
    boundary: Arc<distance::PlanarBoundary>,
    source_orientation: u8,
}

type NamedPlanarRegions = BTreeMap<usize, Option<NamedPlanarRegion>>;

#[derive(Clone, Debug, Default)]
pub(super) struct PlanarMembershipCache(OnceLock<Arc<Mutex<NamedPlanarRegions>>>);

fn coordinates(p: &Point3) -> [Real; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
fn cone_witness(
    surface: &catallax_shapes::Surface,
    target: &Point3,
    step: &Rational,
    offset: i64,
) -> Result<Option<Point3>, Refusal> {
    let catallax_shapes::Surface::Cone(s) = surface else {
        return Ok(None);
    };
    let (Some(axis), Some(reference), Some(slope)) = (
        s.axis.as_rational(),
        s.ref_dir.as_rational(),
        s.slope.rational_ref(),
    ) else {
        return Ok(None);
    };
    let Some(axis_len) = crate::profile::sqrt(&axis.norm_squared()) else {
        return Ok(None);
    };
    let Some(ref_len) = crate::profile::sqrt(&reference.norm_squared()) else {
        return Ok(None);
    };
    let axis = axis.scaled(&axis_len.recip().unwrap());
    let u = reference.scaled(&ref_len.recip().unwrap());
    let v = axis.cross(&u);
    let d = target.minus(&s.apex);
    let approx = |v: &catallax_shapes::RVector3| {
        let mid = d.dot(&v.clone().to_real()).approx().mid();
        Rational::from_f64_exact(mid).map(|q| Ball::from_rational(&q, 14).mid())
    };
    let Some(mut height) = approx(&axis) else {
        return Ok(None);
    };
    height = match offset {
        0 => height,
        1 => &height - step,
        2 => &height + step,
        _ => return Ok(None),
    };
    if height <= Rational::zero() {
        return Ok(None);
    }
    let (Some(x), Some(y)) = (approx(&u), approx(&v)) else {
        return Ok(None);
    };
    let r_xy = Ball::from_rational(&(&x * &x + &y * &y), 16)
        .sqrt()
        .unwrap()
        .mid();
    let [cos, sin] = if (&r_xy + &x).is_zero() {
        [-Rational::one(), Rational::zero()]
    } else {
        let t = Ball::from_rational(&(&y / &(&r_xy + &x)), 14).mid();
        let t2 = &t * &t;
        let denom = Rational::one() + t2.clone();
        [
            (Rational::one() - t2) / denom.clone(),
            Rational::from_i64(2) * t / denom,
        ]
    };
    let radius = slope * &height;
    let offset = axis
        .scaled(&height)
        .plus(&u.scaled(&(&radius * &cos)))
        .plus(&v.scaled(&(&radius * &sin)));
    Ok(Some(s.apex.plus(&offset.to_real())))
}
fn location(leaf: &geometry::Leaf, p: &Point3, dec: &mut Decider) -> Result<Location, Refusal> {
    if let Some(p) = p.as_rational() {
        return leaf.location(&[p.x, p.y, p.z], dec);
    }
    match leaf {
        geometry::Leaf::Feature(f) => f.classify(p, dec),
        geometry::Leaf::Draft(region) => region.classify(p, None, dec),
        geometry::Leaf::Lattice(_) | geometry::Leaf::Gyroid(_) | geometry::Leaf::Mesh(_) => {
            Ok(Location::Indeterminate)
        }
        geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_) => {
            let mut boundary = false;
            let xyz = coordinates(p);
            let box_xyz = [
                geometry::enclose(&xyz[0], 128, dec)?,
                geometry::enclose(&xyz[1], 128, dec)?,
                geometry::enclose(&xyz[2], 128, dec)?,
            ];
            for c in leaf.constraints() {
                let fast = c.polynomial.eval_ball(&box_xyz, dec)?;
                if fast.lo() > Rational::zero() {
                    return Ok(Location::Outside);
                }
                if fast.hi() < Rational::zero() {
                    continue;
                }
                let value = match c.surface() {
                    catallax_shapes::Surface::Cone(s) if !c.surface().has_rational_definition() => {
                        let d = p.minus(&s.apex);
                        d.cross(&s.axis).norm_squared() - s.slope.square() * d.dot(&s.axis).square()
                    }
                    catallax_shapes::Surface::Plane(s)
                        if !c.surface().has_rational_definition() =>
                    {
                        p.minus(&s.origin).dot(&s.normal)
                    }
                    _ => c.polynomial.eval_real(&xyz),
                };
                match geometry::sign(&value, dec)? {
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
    }
}
fn ray(
    leaf: &geometry::Leaf,
    p: &Point3,
    d: &Vector3,
    dec: &mut Decider,
) -> Result<Location, Refusal> {
    if let geometry::Leaf::Draft(region) = leaf {
        return region.classify(p, Some(d), dec);
    }
    if let Some(p) = p.as_rational() {
        return leaf.ray(&[p.x, p.y, p.z], d, dec);
    }
    if let geometry::Leaf::Feature(feature) = leaf {
        return feature.classify(
            &planar::ray_sample(leaf.constraints().iter(), p, d, dec)?,
            dec,
        );
    }
    if matches!(
        leaf,
        geometry::Leaf::Lattice(_) | geometry::Leaf::Gyroid(_) | geometry::Leaf::Mesh(_)
    ) {
        return Ok(Location::Indeterminate);
    }
    let xyz = coordinates(p);
    let mut boundary = false;
    let mut epsilon = Rational::one();
    let coefficients = crate::workers::map(leaf.constraints(), |c| {
        geometry::point_ray_coefficients(c, &xyz, d)
    });
    for coefficients in coefficients {
        let mut found = false;
        for (i, c) in coefficients.iter().enumerate() {
            let fast = geometry::enclose(c, 128, dec)?;
            let sign = if fast.lo() > Rational::zero() {
                Sign::Positive
            } else if fast.hi() < Rational::zero() {
                Sign::Negative
            } else {
                geometry::sign(c, dec)?
            };
            if sign == Sign::Zero {
                continue;
            }
            found = true;
            if matches!(leaf, geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_)) {
                if sign == Sign::Positive {
                    return Ok(Location::Outside);
                }
                break;
            }
            let lo_mag = fast.lo().abs().min(fast.hi().abs());
            if !fast.contains_zero() && lo_mag > Rational::zero() {
                let mut tail_hi = Rational::one();
                for tail_c in &coefficients[i + 1..] {
                    let tb = geometry::enclose(tail_c, 128, dec)?;
                    tail_hi = tail_hi + tb.lo().abs().max(tb.hi().abs()) + Rational::one();
                }
                let step = (&lo_mag / &(&lo_mag + &tail_hi)).mul_pow2(-1);
                if step > Rational::zero() {
                    epsilon = epsilon.min(step);
                    break;
                }
            }
            let magnitude = if sign == Sign::Negative {
                -c
            } else {
                c.clone()
            };
            let denominator = coefficients[i + 1..]
                .iter()
                .fold(&magnitude + &Real::one(), |sum, c| {
                    sum + c.square() + Real::one()
                });
            let bound = magnitude / denominator;
            let mut bits = 128;
            loop {
                let lo = dec.enclose_budgeted_at(&bound, bits)?.lo();
                if lo > Rational::zero() {
                    epsilon = epsilon.min(lo);
                    break;
                }
                bits = bits.checked_mul(2).ok_or_else(|| {
                    Refusal::new(
                        "numbers.refinement-budget",
                        "distance ray separation precision overflow",
                    )
                })?;
            }
            break;
        }
        boundary |= !found;
    }
    if matches!(leaf, geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_)) {
        Ok(if boundary {
            Location::Boundary
        } else {
            Location::Inside
        })
    } else {
        location(leaf, &p.plus(&d.scaled(&Real::from(epsilon))), dec)
    }
}

// A physical patch is the union of its convex tiles; each tile is a conjunction.
type PlanarTrim = Vec<Vec<geometry::Constraint>>;

/// Prepared field-native distance domains. Reuse across point queries on a plan.
/// Face handles are valid only for the plan used to construct this object.
pub struct Prepared<'a> {
    plan: &'a Plan,
    domains: Vec<Option<Domain>>,
    coincident: Vec<Vec<Vec<bool>>>,
    planar_trims: Vec<Vec<PlanarTrim>>,
    planar: OnceLock<planar::Sections<'a>>,
    planar_init: Mutex<()>,
}
impl Plan {
    /// Prepares named chart covers only when every live leaf has registered faces.
    /// Use this entry point for whole-boundary observations. The named-face API
    /// [`Self::prepare_distances`] alone does not establish this coverage.
    pub fn prepare_boundary_distances(&self, dec: &mut Decider) -> Result<Prepared<'_>, Refusal> {
        // Match the solver coverage rule: aliases can register faces for shared
        // geometry, and compacted-away leaves do not belong to the live result.
        let mut named = vec![false; self.leaves.len()];
        for face in &self.faces {
            named[face.leaf] = true;
        }
        if (0..self.leaves.len()).any(|i| self.leaf_registers[i].is_some() && !named[i]) {
            return Err(Refusal::new(
                "solver.distance-coverage",
                "a live field leaf has no named boundary faces; whole-boundary distance is unproved",
            ));
        }
        self.prepare_distances(dec)
    }

    /// Prepares compact chart covers from recipe bounds, never from a B-rep mesh.
    pub fn prepare_distances(&self, dec: &mut Decider) -> Result<Prepared<'_>, Refusal> {
        let mut bounds: Vec<Option<Rational>> = Vec::with_capacity(self.instructions.len());
        for instruction in &self.instructions {
            let bound = match *instruction {
                Instruction::Empty => Some(Rational::zero()),
                Instruction::Leaf(i) => self.leaves[i].display_radius(dec)?,
                Instruction::Boolean { kind, a, b } => match kind {
                    catallax_recipe::BooleanKind::Subtract => bounds[a].clone(),
                    catallax_recipe::BooleanKind::Union => match (&bounds[a], &bounds[b]) {
                        (Some(a), Some(b)) => Some(a.max(b).clone()),
                        _ => None,
                    },
                    catallax_recipe::BooleanKind::Intersect => match (&bounds[a], &bounds[b]) {
                        (Some(a), Some(b)) => Some(a.min(b).clone()),
                        (a, b) => a.clone().or_else(|| b.clone()),
                    },
                },
            };
            bounds.push(bound);
        }
        let radius = bounds[self.root].as_ref().ok_or_else(|| {
            Refusal::new(
                "field.distance-unbounded",
                "result has no certified compact recipe bound",
            )
        })?;
        let bound = Ball::from_rational(&(-radius), 128).hull(&Ball::from_rational(radius, 128));
        let mut domains = Vec::with_capacity(self.faces.len());
        let mut coincident = Vec::with_capacity(self.faces.len());
        let mut constraint_supports = Vec::with_capacity(self.leaves.len());
        let mut planar_trims = Vec::with_capacity(self.leaves.len());
        let mut leaf_boxes = Vec::with_capacity(self.leaves.len());
        for (i, leaf) in self.leaves.iter().enumerate() {
            let mut supports = Vec::new();
            if self.leaf_registers[i].is_some() {
                supports.reserve(leaf.constraints().len());
                for c in leaf.constraints() {
                    supports.push(geometry::support_geometry(&c.polynomial, dec)?.0);
                }
                leaf_boxes.push(leaf.display_box(dec)?);
            } else {
                leaf_boxes.push(None);
            }
            let mut trims = Vec::new();
            if self.leaf_registers[i].is_some()
                && let geometry::Leaf::Draft(region) = leaf
                && let crate::draft::Cell::Polyhedral(poly) = &region.cell
            {
                for patch in &poly.faces {
                    let normal = poly.plane(patch).normal;
                    let mut tiles = Vec::with_capacity(patch.tiles.len());
                    for tile in &patch.tiles {
                        let mut sides = Vec::with_capacity(tile.len());
                        for j in 0..tile.len() {
                            let a = &poly.vertices[tile[j]];
                            let edge = poly.vertices[tile[(j + 1) % tile.len()]].minus(a);
                            sides.push(geometry::constraint(catallax_shapes::Surface::Plane(
                                catallax_shapes::Plane {
                                    origin: a.clone(),
                                    normal: edge.cross(&normal),
                                    ref_dir: edge,
                                },
                            ))?);
                        }
                        tiles.push(sides);
                    }
                    trims.push(tiles);
                }
            }
            planar_trims.push(trims);
            constraint_supports.push(supports);
        }
        for face in &self.faces {
            if !face.in_result || self.leaf_registers[face.leaf].is_none() {
                domains.push(None);
                coincident.push(Vec::new());
                continue;
            }
            let mut face_bounds = std::array::from_fn(|_| bound.clone());
            if let Some(lb) = &leaf_boxes[face.leaf] {
                for k in 0..3 {
                    let lo = face_bounds[k].lo().max(lb[k].lo());
                    let hi = face_bounds[k].hi().min(lb[k].hi());
                    if lo <= hi {
                        face_bounds[k] =
                            Ball::from_rational(&lo, 128).hull(&Ball::from_rational(&hi, 128));
                    }
                }
            }
            domains.push(Some(Domain::new(&face.surface, face_bounds, dec)?));
            let mut same = Vec::with_capacity(self.leaves.len());
            for supports in &constraint_supports {
                let mut constraints = Vec::with_capacity(supports.len());
                for support in supports {
                    constraints.push(face.support.equals(support, dec)?);
                }
                same.push(constraints);
            }
            coincident.push(same);
        }
        Ok(Prepared {
            plan: self,
            domains,
            coincident,
            planar_trims,
            planar: OnceLock::new(),
            planar_init: Mutex::new(()),
        })
    }
}
impl Prepared<'_> {
    pub(crate) fn belongs_to(&self, plan: &Plan) -> bool {
        std::ptr::eq(self.plan, plan)
    }

    fn planar_boundary(
        &self,
        id: usize,
        dec: &mut Decider,
    ) -> Result<Option<&distance::PlanarBoundary>, Refusal> {
        if let Some(boundary) = self.planar.get().and_then(|sections| sections.get(id)) {
            return Ok(boundary.as_deref());
        }
        let _initializing = self
            .planar_init
            .lock()
            .expect("planar distance initialization panicked");
        if self.planar.get().is_none() {
            let sections = self.plan.prepare_planar_boundaries(dec)?;
            assert!(
                self.planar.set(sections).is_ok(),
                "planar distance sections initialize under one lock"
            );
        }
        let sections = self.planar.get().expect("planar sections were initialized");
        sections.prepare(id, dec)?;
        Ok(sections
            .get(id)
            .expect("planar section was initialized")
            .as_deref())
    }
    fn member(&self, id: usize, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        if let Some(planar) = self
            .planar
            .get()
            .and_then(|sections| sections.get(id))
            .and_then(Option::as_deref)
        {
            return planar.contains(p, dec);
        }
        let face = &self.plan.faces[id];
        let leaf = &self.plan.leaves[face.leaf];
        if location(leaf, p, dec)? != Location::Boundary {
            return Ok(false);
        }
        if let geometry::Leaf::Feature(f) = leaf
            && !f.on_face_point(self.plan.definition_role(id), p, dec)?
        {
            return Ok(false);
        }
        if let geometry::Leaf::Draft(region) = leaf
            && !region.on_face_point(&face.origin, p, dec)?
        {
            return Ok(false);
        }
        if matches!(self.plan.instructions[self.plan.root], Instruction::Leaf(i) if i == face.leaf)
        {
            return Ok(true);
        }
        let leaf_locations = self
            .plan
            .leaves
            .iter()
            .enumerate()
            .map(|(i, leaf)| {
                if i == face.leaf {
                    Ok(Location::Boundary)
                } else if self.plan.leaf_registers[i].is_none() {
                    Ok(Location::Outside)
                } else {
                    location(leaf, p, dec)
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut tape = Vec::with_capacity(self.plan.instructions.len());
        if leaf_locations.iter().enumerate().all(|(i, loc)| {
            i == face.leaf
                || self.plan.leaf_registers[i].is_none()
                || matches!(loc, Location::Inside | Location::Outside)
        }) {
            // Only this source boundary can change membership locally. Its two
            // material sides exist by the source-face predicate, even at a
            // closed edge. Evaluate both Boolean assignments without rebuilding
            // radial roots and ray-germ separation bounds for the same source.
            let mut result = [Location::Indeterminate; 2];
            for (j, side) in [Location::Inside, Location::Outside]
                .into_iter()
                .enumerate()
            {
                tape.clear();
                for instruction in &self.plan.instructions {
                    let loc = match *instruction {
                        Instruction::Empty => Location::Outside,
                        Instruction::Leaf(i) => {
                            if i == face.leaf {
                                side
                            } else {
                                leaf_locations[i]
                            }
                        }
                        Instruction::Boolean { kind, a, b } => {
                            crate::combine(kind, tape[a], tape[b])
                        }
                    };
                    tape.push(loc);
                }
                result[j] = tape[self.plan.root];
            }
            return Ok(matches!(
                result,
                [Location::Inside, Location::Outside] | [Location::Outside, Location::Inside]
            ));
        }
        let xyz = coordinates(p);
        // Other active supports need a closed-trim incidence proof. Use nearby
        // open-face witnesses instead, so a buried tangent face cannot donate a
        // spurious upper bound at its isolated contact with the result boundary.
        for (i, other) in self.plan.leaves.iter().enumerate() {
            if self.plan.leaf_registers[i].is_none() {
                continue;
            }
            for (j, c) in other.constraints().iter().enumerate() {
                if !self.coincident[id][i][j]
                    && geometry::sign(&c.polynomial.eval_real(&xyz), dec)? == Sign::Zero
                {
                    return Ok(false);
                }
            }
        }
        let n = Vector3::new(
            face.gradient[0].eval_real(&xyz),
            face.gradient[1].eval_real(&xyz),
            face.gradient[2].eval_real(&xyz),
        );
        let mut inside = false;
        let mut outside = false;
        tape.clear();
        // Opposite normal germs prove both material and complement arbitrarily
        // near this admitted open-face point of the regularised result.
        for d in [n.clone(), n.scaled(&Real::from_i64(-1))] {
            if geometry::sign(&d.norm_squared(), dec)? == Sign::Zero {
                continue;
            }
            tape.clear();
            for instruction in &self.plan.instructions {
                let loc = match *instruction {
                    Instruction::Empty => Location::Outside,
                    Instruction::Leaf(i) => {
                        let loc = leaf_locations[i];
                        if matches!(loc, Location::Inside | Location::Outside) {
                            loc
                        } else {
                            ray(&self.plan.leaves[i], p, &d, dec)?
                        }
                    }
                    Instruction::Boolean { kind, a, b } => crate::combine(kind, tape[a], tape[b]),
                };
                tape.push(loc);
            }
            inside |= tape[self.plan.root] == Location::Inside;
            outside |= tape[self.plan.root] == Location::Outside;
            if inside && outside {
                return Ok(true);
            }
        }
        Ok(false)
    }
    fn classify(
        &self,
        id: usize,
        p: &Point3,
        cover: &[Ball; 3],
        dec: &mut Decider,
    ) -> Result<Membership, Refusal> {
        let face = &self.plan.faces[id];
        if matches!(
            &self.plan.leaves[face.leaf],
            geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_)
        ) {
            for (j, c) in self.plan.leaves[face.leaf].constraints().iter().enumerate() {
                if !self.coincident[id][face.leaf][j]
                    && c.polynomial.eval_ball(cover, dec)?.lo() > Rational::zero()
                {
                    return Ok(Membership {
                        witness: false,
                        excluded: true,
                    });
                }
            }
        }
        let witness = self.member(id, p, dec)?;
        if witness {
            return Ok(Membership {
                witness: true,
                excluded: false,
            });
        }
        if let geometry::Leaf::Draft(region) = &self.plan.leaves[face.leaf]
            && let crate::draft::Cell::Polyhedral(poly) = &region.cell
        {
            let source = region.source_origin(&face.origin);
            let mut possible = false;
            // Region constraints enumerate physical faces in this same order.
            for ((patch, tiles), &same_carrier) in poly
                .faces
                .iter()
                .zip(&self.planar_trims[face.leaf])
                .zip(&self.coincident[id][face.leaf])
            {
                if !same_carrier || source.is_none_or(|source| !patch.origins.contains(source)) {
                    continue;
                }
                for tile in tiles {
                    let mut separated = false;
                    for side in tile {
                        if side.polynomial.eval_ball(cover, dec)?.lo() > Rational::zero() {
                            separated = true;
                            break;
                        }
                    }
                    if !separated {
                        possible = true;
                        break;
                    }
                }
                if possible {
                    break;
                }
            }
            if !possible {
                return Ok(Membership {
                    witness: false,
                    excluded: true,
                });
            }
        }
        let mut stable = true;
        for (i, leaf) in self.plan.leaves.iter().enumerate() {
            if self.plan.leaf_registers[i].is_none() {
                continue;
            }
            for (j, c) in leaf.constraints().iter().enumerate() {
                if self.coincident[id][i][j] {
                    continue;
                }
                if c.polynomial.eval_ball(cover, dec)?.contains_zero() {
                    stable = false;
                    break;
                }
            }
            if !stable {
                break;
            }
        }
        Ok(Membership {
            witness: false,
            excluded: stable,
        })
    }
    fn domain(&self, id: FaceId) -> Result<&Domain, Refusal> {
        self.domains
            .get(id.0)
            .and_then(Option::as_ref)
            .ok_or_else(|| {
                Refusal::new("field.face-id", "distance face is absent from this result")
            })
    }
    /// Searches for a certified boundary witness without constructing trim regions.
    /// A missing upper means no witness was found, not that the face is empty.
    /// The lower bound is zero; only an admitted physical point supplies an upper.
    pub fn point_to_face_witness(
        &self,
        id: FaceId,
        point: &Point3,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        let domain = self.domain(id)?;
        if width <= &Rational::zero() {
            return Err(Refusal::new(
                "shapes.distance-domain",
                "positive width required",
            ));
        }
        self.support_witness(id, domain, point, width, true, dec)
    }

    fn support_witness(
        &self,
        id: FaceId,
        domain: &Domain,
        point: &Point3,
        width: &Rational,
        perturb: bool,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        let mut candidate = Bounds {
            lower: Rational::zero(),
            upper: None,
            witnesses: None,
            converged: false,
            cells: 0,
        };
        let step = width.mul_pow2(-3);
        for offset in 0..if perturb { 28 } else { 1 } {
            dec.check_deadline()?;
            let target = if offset == 0 {
                point.clone()
            } else {
                let n = offset - 1;
                point.plus(&Vector3::new(
                    Real::from(&Rational::from_i64(n % 3 - 1) * &step),
                    Real::from(&Rational::from_i64(n / 3 % 3 - 1) * &step),
                    Real::from(&Rational::from_i64(n / 9 - 1) * &step),
                ))
            };
            // Reuse the rational analytic support witness, not a B-rep face or
            // mesher. Rational frames avoid introducing needless radial radicals.
            let p = match catallax_boolean::surface_witness(
                domain.surface(),
                &target,
                &width.mul_pow2(-20),
                dec,
            )? {
                Some(p) => p.to_real(),
                None => match cone_witness(domain.surface(), point, &step, offset)? {
                    Some(p) => p,
                    None => break,
                },
            };
            let squared = point.minus(&p).norm_squared();
            let distance = if let Some(q) = squared.rational_ref() {
                Ball::from_rational(q, 128)
            } else {
                dec.enclose_budgeted_at(&squared, 128)?
            }
            .sqrt()
            .unwrap();
            if offset > 0 && distance.lo() > width.mul_pow2(1) {
                break;
            }
            if self.member(id.0, &p, dec)?
                && candidate.upper.as_ref().is_none_or(|u| distance.hi() < *u)
            {
                candidate.upper = Some(distance.hi());
                candidate.witnesses = Some([point.clone(), p]);
            }
            if candidate.upper.as_ref().is_some_and(|u| u <= width) {
                candidate.converged = true;
                return Ok(candidate);
            }
            if offset == 0 && distance.lo() > *width {
                break;
            }
        }
        Ok(candidate)
    }

    /// Certified point-to-face distance, including the face's closed trim loops.
    pub fn point_to_face(
        &self,
        id: FaceId,
        point: &Point3,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        let domain = self.domain(id)?;
        if options.width <= Rational::zero() || options.max_cells == 0 {
            return Err(Refusal::new(
                "shapes.distance-domain",
                "positive width and cell count required",
            ));
        }
        if let Some(planar) = self
            .planar
            .get()
            .and_then(|sections| sections.get(id.0))
            .and_then(Option::as_deref)
        {
            return planar.point_distance(point, options, dec);
        }
        let candidate = self.support_witness(
            id,
            domain,
            point,
            &options.width,
            options.max_cells == 1,
            dec,
        )?;
        if candidate.converged {
            return Ok(candidate);
        }
        // An admitted near witness already encloses the minimum to the requested
        // width. Reconstruct complete finite regions only when that proof fails.
        if matches!(domain.surface(), catallax_shapes::Surface::Plane(_))
            && let Some(planar) = self.planar_boundary(id.0, dec)?
        {
            return planar.point_distance(point, options, dec);
        }
        if options.max_cells == 1 {
            return Ok(candidate);
        }
        let mut result = distance::point_distance(
            domain,
            point,
            options,
            &mut |p: &Point3, b: &[Ball; 3], d: &mut Decider| self.classify(id.0, p, b, d),
            dec,
        )?;
        if let Some(upper) = candidate.upper
            && result.upper.as_ref().is_none_or(|u| upper < *u)
        {
            result.upper = Some(upper);
            result.witnesses = candidate.witnesses;
            result.converged = result
                .upper
                .as_ref()
                .is_some_and(|u| u - &result.lower <= options.width);
        }
        Ok(result)
    }
    /// Certified minimum distance between selected trimmed field faces.
    pub fn face_to_face(
        &self,
        a: FaceId,
        b: FaceId,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        let a_domain = self.domain(a)?;
        let b_domain = self.domain(b)?;
        if matches!(a_domain.surface(), catallax_shapes::Surface::Plane(_))
            && matches!(b_domain.surface(), catallax_shapes::Surface::Plane(_))
            && let (Some(first), Some(second)) = (
                self.planar_boundary(a.0, dec)?,
                self.planar_boundary(b.0, dec)?,
            )
            && let Some(answer) = first.face_distance(second, options, dec)?
        {
            return Ok(answer);
        }
        distance::face_distance(
            a_domain,
            b_domain,
            options,
            &mut |p: &Point3, c: &[Ball; 3], d: &mut Decider| self.classify(a.0, p, c, d),
            &mut |p: &Point3, c: &[Ball; 3], d: &mut Decider| self.classify(b.0, p, c, d),
            dec,
        )
    }
    /// Bounds distance to the union of all surviving named faces. A missing
    /// witness remains unbounded; empty or buried faces never supply an upper.
    pub fn point_to_boundary(
        &self,
        point: &Point3,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Bounds, Refusal> {
        if options.width <= Rational::zero() || options.max_cells == 0 {
            return Err(Refusal::new(
                "shapes.distance-domain",
                "positive width and cell count required",
            ));
        }
        let mut result = Bounds {
            lower: Rational::zero(),
            upper: None,
            witnesses: None,
            converged: false,
            cells: 0,
        };
        let mut remaining = options.max_cells;
        let mut face_lower: Vec<Option<Rational>> = vec![None; self.domains.len()];
        let mut face_settled = vec![false; self.domains.len()];
        for (i, domain) in self.domains.iter().enumerate() {
            let Some(domain) = domain else {
                face_settled[i] = true;
                continue;
            };
            if matches!(domain.surface(), catallax_shapes::Surface::Plane(_))
                && let Some(planar) = self.planar_boundary(i, dec)?
            {
                let bounds = planar.point_distance(point, options, dec)?;
                face_lower[i] = Some(bounds.lower);
                face_settled[i] = bounds.converged;
                result.cells += bounds.cells;
                if let Some(upper) = bounds.upper
                    && result.upper.as_ref().is_none_or(|u| upper < *u)
                {
                    result.upper = Some(upper);
                    result.witnesses = bounds.witnesses;
                }
                if result.upper.as_ref().is_some_and(|u| u <= &options.width) {
                    result.converged = true;
                    return Ok(result);
                }
            }
        }
        for stage in [2usize, 16, 128, 1024, options.max_cells] {
            let cap = stage.min(options.max_cells);
            for (i, domain) in self.domains.iter().enumerate() {
                if domain.is_none() || face_settled[i] {
                    continue;
                }
                if let (Some(lo), Some(best)) = (&face_lower[i], &result.upper)
                    && lo >= best
                {
                    face_settled[i] = true;
                    continue;
                }
                if remaining == 0 {
                    break;
                }
                let refinement = Options {
                    width: options.width.clone(),
                    max_cells: cap.min(remaining),
                };
                let bounds = self.point_to_face(FaceId(i), point, &refinement, dec)?;
                result.cells += bounds.cells;
                remaining = remaining.saturating_sub(bounds.cells);
                let next_lo = face_lower[i]
                    .take()
                    .map_or(bounds.lower.clone(), |old| old.max(bounds.lower));
                face_settled[i] =
                    bounds.converged || result.upper.as_ref().is_some_and(|best| &next_lo >= best);
                face_lower[i] = Some(next_lo);
                if let Some(upper) = bounds.upper
                    && result.upper.as_ref().is_none_or(|u| upper < *u)
                {
                    result.upper = Some(upper);
                    result.witnesses = bounds.witnesses;
                }
                if result.upper.as_ref().is_some_and(|u| u <= &options.width) {
                    result.converged = true;
                    return Ok(result);
                }
            }
            if face_settled.iter().all(|&s| s) || remaining == 0 || cap >= options.max_cells {
                break;
            }
        }
        if (0..self.domains.len()).any(|i| self.domains[i].is_some() && face_lower[i].is_none()) {
            result.lower = Rational::zero();
        } else {
            let min_lo = face_lower
                .into_iter()
                .flatten()
                .min()
                .unwrap_or_else(Rational::zero);
            result.lower = result
                .upper
                .as_ref()
                .map_or(min_lo.clone(), |u| min_lo.min(u.clone()));
        }
        result.converged = result
            .upper
            .as_ref()
            .is_some_and(|u| u - &result.lower <= options.width);
        Ok(result)
    }
}
