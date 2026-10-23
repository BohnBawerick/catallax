//! Trimmed distances from the field's own recipe predicates, without B-rep input.

use crate::{FaceId, Instruction, Location, Plan, geometry};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{
    Point3, Vector3,
    distance::{self, Bounds, Domain, Membership, Options},
};

fn coordinates(p: &Point3) -> [Real; 3] {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}
fn location(leaf: &geometry::Leaf, p: &Point3, dec: &mut Decider) -> Result<Location, Refusal> {
    match leaf {
        geometry::Leaf::Feature(f) => f.classify(p, dec),
        geometry::Leaf::Primitive(_) => {
            let mut boundary = false;
            let p = coordinates(p);
            for c in leaf.constraints() {
                match geometry::sign(&c.polynomial.eval_real(&p), dec)? {
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
    if let Some(p) = p.as_rational() {
        return leaf.ray(&[p.x, p.y, p.z], d, dec);
    }
    let xyz = coordinates(p);
    let mut boundary = false;
    let mut epsilon = Rational::one();
    for c in leaf.constraints() {
        let coefficients = geometry::point_ray_coefficients(c, &xyz, d);
        let mut found = false;
        for (i, c) in coefficients.iter().enumerate() {
            let sign = geometry::sign(c, dec)?;
            if sign == Sign::Zero {
                continue;
            }
            found = true;
            if matches!(leaf, geometry::Leaf::Primitive(_)) {
                if sign == Sign::Positive {
                    return Ok(Location::Outside);
                }
                break;
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
    if matches!(leaf, geometry::Leaf::Primitive(_)) {
        Ok(if boundary {
            Location::Boundary
        } else {
            Location::Inside
        })
    } else {
        location(leaf, &p.plus(&d.scaled(&Real::from(epsilon))), dec)
    }
}

/// Prepared field-native distance domains. Reuse across point queries on a plan.
/// Face handles are valid only for the plan used to construct this object.
pub struct Prepared<'a> {
    plan: &'a Plan,
    domains: Vec<Option<Domain>>,
    coincident: Vec<Vec<Vec<bool>>>,
}
impl Plan {
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
        for (i, leaf) in self.leaves.iter().enumerate() {
            let mut supports = Vec::new();
            if self.leaf_registers[i].is_some() {
                supports.reserve(leaf.constraints().len());
                for c in leaf.constraints() {
                    supports.push(geometry::support_geometry(&c.polynomial, dec)?.0);
                }
            }
            constraint_supports.push(supports);
        }
        for face in &self.faces {
            if self.leaf_registers[face.leaf].is_none() {
                domains.push(None);
                coincident.push(Vec::new());
                continue;
            }
            domains.push(Some(Domain::new(
                &face.surface,
                std::array::from_fn(|_| bound.clone()),
                dec,
            )?));
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
        })
    }
}
impl Prepared<'_> {
    fn member(&self, id: usize, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let face = &self.plan.faces[id];
        let leaf = &self.plan.leaves[face.leaf];
        if location(leaf, p, dec)? != Location::Boundary {
            return Ok(false);
        }
        if let geometry::Leaf::Feature(f) = leaf
            && !f.on_face_point(&face.origin.role, p, dec)?
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
        let witness = self.member(id, p, dec)?;
        if witness {
            return Ok(Membership {
                witness: true,
                excluded: false,
            });
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
            .ok_or_else(|| Refusal::new("field.face-id", "distance face is absent from this plan"))
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
        let mut candidate = Bounds {
            lower: Rational::zero(),
            upper: None,
            witnesses: None,
            converged: false,
            cells: 0,
        };
        let step = options.width.mul_pow2(-3);
        for offset in 0..if options.max_cells == 1 { 28 } else { 1 } {
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
            let Some(p) = catallax_boolean::surface_witness(
                domain.surface(),
                &target,
                &options.width.mul_pow2(-20),
                dec,
            )?
            else {
                break;
            };
            let p = p.to_real();
            let squared = point.minus(&p).norm_squared();
            let distance = if let Some(q) = squared.rational_ref() {
                Ball::from_rational(q, 128)
            } else {
                dec.enclose_budgeted_at(&squared, 128)?
            }
            .sqrt()
            .unwrap();
            if self.member(id.0, &p, dec)?
                && candidate.upper.as_ref().is_none_or(|u| distance.hi() < *u)
            {
                candidate.upper = Some(distance.hi());
                candidate.witnesses = Some([point.clone(), p]);
            }
            if candidate
                .upper
                .as_ref()
                .is_some_and(|u| u <= &options.width)
            {
                candidate.converged = true;
                return Ok(candidate);
            }
            if offset == 0 && distance.lo() > options.width {
                break;
            }
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
        distance::face_distance(
            self.domain(a)?,
            self.domain(b)?,
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
        let probe = Options {
            width: options.width.clone(),
            max_cells: 1,
        };
        let mut remaining = options.max_cells;
        for pass in 0..2 {
            let mut lower = None;
            for (i, domain) in self.domains.iter().enumerate() {
                if domain.is_none() {
                    continue;
                }
                if pass == 1 && remaining == 0 {
                    result.lower = Rational::zero();
                    return Ok(result);
                }
                let refinement = Options {
                    width: options.width.clone(),
                    max_cells: remaining,
                };
                let bounds = self.point_to_face(
                    FaceId(i),
                    point,
                    if pass == 0 { &probe } else { &refinement },
                    dec,
                )?;
                lower = Some(lower.map_or(bounds.lower.clone(), |x: Rational| x.min(bounds.lower)));
                result.cells += bounds.cells;
                remaining = remaining.saturating_sub(bounds.cells);
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
            result.lower = lower.unwrap_or_else(Rational::zero);
        }
        result.converged = result
            .upper
            .as_ref()
            .is_some_and(|u| u - &result.lower <= options.width);
        Ok(result)
    }
}
