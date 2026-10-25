//! Closed named-patch incidence at certified support contacts.
//!
//! Two independent smooth supports, or an exact plane arrangement, admit local
//! polynomial sign coordinates. At a tangent seam, transverse trim dominance
//! or strictly separated second-order graphs can prove an open source patch.
//! Every proof retains source/result and historical Boolean ownership.
//! Tangent-plane samples never assign signs to tangent supports. Singular or
//! unproved contacts remain unresolved.

use crate::{
    Direction, FaceId, FaceQuery, FaceState, Instruction, Location, Plan, Workspace, geometry,
};
use catallax_numbers::{Decider, Rational, Real, RealPoly, Refusal, Sign};
use catallax_shapes::{RPoint3, RVector3, Surface, Torus, Vector3};

fn unknown() -> FaceQuery {
    FaceQuery {
        state: FaceState::Indeterminate,
        direction: Direction::Ambiguous,
    }
}

fn positive_lower(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    if let Some(q) = value.rational_ref() {
        return Ok(q.clone());
    }
    let mut bits = crate::metric::PRECISION;
    loop {
        let lo = geometry::enclose(value, bits, dec)?.lo();
        if lo > Rational::zero() {
            return Ok(lo);
        }
        bits = bits.checked_mul(2).ok_or_else(|| {
            Refusal::new(
                "numbers.refinement-budget",
                "contact-sector separation exceeds the precision index",
            )
        })?;
    }
}

fn magnitude(value: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    Ok(if geometry::sign(value, dec)? == Sign::Negative {
        -value
    } else {
        value.clone()
    })
}

// Choose a positive rational e for which every nonzero leading linear term
// retains its sign in a + e*b and a - e*b. Zero leading terms use b instead.
fn epsilon(
    normals: &[Vector3],
    a: &Vector3,
    b: &Vector3,
    dec: &mut Decider,
) -> Result<Option<Rational>, Refusal> {
    let mut result = Rational::one();
    for normal in normals {
        dec.check_deadline()?;
        let leading = normal.dot(a);
        let tail = normal.dot(b);
        if geometry::sign(&leading, dec)? == Sign::Zero {
            if geometry::sign(&tail, dec)? == Sign::Zero {
                return Ok(None);
            }
            continue;
        }
        let leading = magnitude(&leading, dec)?;
        let tail = if let Some(q) = tail.rational_ref() {
            Real::from(q.abs())
        } else {
            tail.square() + Real::one()
        };
        let bound = &leading / &(&leading + &tail + Real::one());
        result = result.min(positive_lower(&bound, dec)?);
    }
    Ok(Some(result))
}

fn tangent_sectors(
    normal: &Vector3,
    others: &[Vector3],
    dec: &mut Decider,
) -> Result<Vec<Vector3>, Refusal> {
    if others.is_empty() {
        for axis in 0..3 {
            let [x, y, z] =
                std::array::from_fn(|i| Rational::from_i64(if axis == i { 1 } else { 0 }));
            let tangent = normal.cross(&RVector3::new(x, y, z).to_real());
            if geometry::sign(&tangent.norm_squared(), dec)? != Sign::Zero {
                return Ok(vec![tangent]);
            }
        }
        return Ok(Vec::new());
    }
    let mut result = Vec::new();
    // Every open sector of the central line arrangement touches an end of one
    // of its lines. Both infinitesimal sides at both ends visit all sectors.
    for cut in others {
        let line = normal.cross(cut);
        let across = normal.cross(&line);
        for sign in [-1, 1] {
            let line = line.scaled(&Real::from_i64(sign));
            let Some(e) = epsilon(others, &line, &across, dec)? else {
                continue;
            };
            for sign in [-1, 1] {
                result.push(
                    line.plus(&across.scaled(&Real::from(e.clone() * Rational::from_i64(sign)))),
                );
            }
        }
    }
    Ok(result)
}

enum LocalConstraint {
    Fixed(Sign),
    Own(Sign),
    Transverse(Vector3),
    TangentBound { sign: Sign, escape: Vector3 },
}
enum LocalLeaf {
    Fixed(Location),
    Cell(Vec<LocalConstraint>),
    Unknown,
}
impl LocalLeaf {
    fn side(&self, tangent: &Vector3, plus: bool, dec: &mut Decider) -> Result<Location, Refusal> {
        let Self::Cell(constraints) = self else {
            return Ok(match self {
                Self::Fixed(location) => *location,
                _ => Location::Indeterminate,
            });
        };
        let mut unknown = false;
        for constraint in constraints {
            let sign = match constraint {
                LocalConstraint::Fixed(sign) => *sign,
                LocalConstraint::Own(sign) => {
                    if plus {
                        *sign
                    } else {
                        match sign {
                            Sign::Positive => Sign::Negative,
                            Sign::Negative => Sign::Positive,
                            Sign::Zero => Sign::Zero,
                        }
                    }
                }
                LocalConstraint::Transverse(normal) => geometry::sign(&normal.dot(tangent), dec)?,
                LocalConstraint::TangentBound { sign, escape } => {
                    if geometry::sign(&escape.dot(tangent), dec)? == Sign::Zero {
                        Sign::Zero
                    } else {
                        *sign
                    }
                }
            };
            match sign {
                Sign::Positive => return Ok(Location::Outside),
                Sign::Zero => unknown = true,
                Sign::Negative => (),
            }
        }
        Ok(if unknown {
            Location::Indeterminate
        } else {
            Location::Inside
        })
    }
}

fn flip_polarity(polarity: u8) -> u8 {
    ((polarity & 1) << 1) | ((polarity & 2) >> 1)
}

// A complete circular tube lies in these three convex parallel supports.
// Equality is a single contact circle. Leaving that circle meridionally gives
// a strict sign on either carrier, unlike an unproved tangent-plane sample.
fn tube_enclosure(
    tube: &Torus,
    convex: &Surface,
    point: &[Real; 3],
    dec: &mut Decider,
) -> Result<Option<Vector3>, Refusal> {
    if geometry::sign(&tube.minor, dec)? != Sign::Positive
        || geometry::sign(&(&tube.major - &tube.minor), dec)? != Sign::Positive
    {
        return Ok(None);
    }
    let zero = |value: &Real, dec: &mut Decider| {
        Ok::<_, Refusal>(geometry::sign(value, dec)? == Sign::Zero)
    };
    Ok(match convex {
        Surface::Cylinder(c)
            if zero(&c.axis.cross(&tube.axis).norm_squared(), dec)?
                && zero(
                    &c.origin
                        .minus(&tube.center)
                        .cross(&tube.axis)
                        .norm_squared(),
                    dec,
                )?
                && zero(&(&c.radius() - &(&tube.major + &tube.minor)), dec)? =>
        {
            Some(tube.axis.clone())
        }
        Surface::Sphere(s)
            if zero(
                &s.center
                    .minus(&tube.center)
                    .cross(&tube.axis)
                    .norm_squared(),
                dec,
            )? && zero(
                &(&s.radius
                    - &(&(&tube.major.square() + &s.center.minus(&tube.center).norm_squared())
                        .sqrt()
                        + &tube.minor)),
                dec,
            )? =>
        {
            Some(tube.axis.clone())
        }
        Surface::Plane(p)
            if zero(&p.normal.cross(&tube.axis).norm_squared(), dec)?
                && zero(
                    &(&tube.center.minus(&p.origin).dot(&p.normal)
                        + &(&tube.minor * &p.normal.norm_squared().sqrt())),
                    dec,
                )? =>
        {
            let q = Vector3::new(
                &point[0] - &tube.center.x,
                &point[1] - &tube.center.y,
                &point[2] - &tube.center.z,
            );
            Some(
                q.minus(
                    &tube
                        .axis
                        .scaled(&(&q.dot(&tube.axis) / &tube.axis.norm_squared())),
                ),
            )
        }
        _ => None,
    })
}

fn tangent_bound(
    own: &Surface,
    other: &Surface,
    point: &[Real; 3],
    dec: &mut Decider,
) -> Result<Option<(Sign, Vector3)>, Refusal> {
    if let Surface::Torus(tube) = own
        && let Some(escape) = tube_enclosure(tube, other, point, dec)?
    {
        return Ok(Some((Sign::Negative, escape)));
    }
    if let Surface::Torus(tube) = other
        && let Some(escape) = tube_enclosure(tube, own, point, dec)?
    {
        return Ok(Some((Sign::Positive, escape)));
    }
    Ok(None)
}

impl Plan {
    /// Tests named-patch closure using transverse sectors, trim dominance or
    /// separated second-order graphs. Contacts without a surviving-patch proof
    /// remain indeterminate; mixed source polarity retains an ambiguous normal.
    /// Draft tile domains retain the ordinary named-face proof unchanged.
    pub fn closed_face_query_checked(
        &self,
        id: FaceId,
        point: &RPoint3,
        storage: &mut Workspace,
        dec: &mut Decider,
    ) -> Result<FaceQuery, Refusal> {
        let ordinary = self.face_query_checked(id, point, storage, dec)?;
        // Draft origins may partition coplanar tiles. Their existing named
        // domain proof remains authoritative: carrier sectors alone cannot
        // establish which tile approaches a closed contact.
        if ordinary.state != FaceState::Indeterminate
            || matches!(
                &self.leaves[self.faces[id.0].leaf],
                geometry::Leaf::Draft(_)
            )
        {
            return Ok(ordinary);
        }
        self.contact_face_query(id, point, dec)
    }

    fn dominated_contact(
        &self,
        id: FaceId,
        p: &[Rational; 3],
        active: &[(RealPoly, Vector3)],
        dec: &mut Decider,
    ) -> Result<Option<FaceQuery>, Refusal> {
        let face = &self.faces[id.0];
        let geometry::Leaf::Cell(source) = &self.leaves[face.leaf] else {
            return Ok(None);
        };
        if source.capped_chart().is_some() {
            return Ok(None);
        }
        // A source leaf used with both Boolean polarities can own oppositely
        // oriented patches. Do not infer one direction from a single witness.
        let mut polarities = Vec::with_capacity(self.instructions.len());
        for instruction in &self.instructions {
            polarities.push(match *instruction {
                Instruction::Empty => 0,
                Instruction::Leaf(leaf) => u8::from(leaf == face.leaf),
                Instruction::Boolean { kind, a, b } => {
                    polarities[a]
                        | if kind == catallax_recipe::BooleanKind::Subtract {
                            flip_polarity(polarities[b])
                        } else {
                            polarities[b]
                        }
                }
            });
        }
        if !matches!(polarities[self.root], 1 | 2) {
            return Ok(None);
        }
        let normal = crate::normal(face, p);
        let mut cuts = Vec::new();
        for (_, other) in active {
            if geometry::sign(&normal.cross(other).norm_squared(), dec)? != Sign::Zero {
                cuts.push(other.clone());
            }
        }
        let tangents = tangent_sectors(&normal, &cuts, dec)?;
        let xyz = p.clone().map(Real::from);
        let mut leaves = Vec::with_capacity(self.leaves.len());
        for (i, leaf) in self.leaves.iter().enumerate() {
            if self.leaf_registers[i].is_none() {
                leaves.push(LocalLeaf::Unknown);
                continue;
            }
            let location = leaf.location(p, dec)?;
            if matches!(location, Location::Inside | Location::Outside) {
                leaves.push(LocalLeaf::Fixed(location));
                continue;
            }
            let geometry::Leaf::Cell(cell) = leaf else {
                leaves.push(LocalLeaf::Unknown);
                continue;
            };
            if cell.capped_chart().is_some() {
                leaves.push(LocalLeaf::Unknown);
                continue;
            }
            let mut constraints = Vec::new();
            for constraint in leaf.constraints() {
                let at = crate::polynomial_sign(&constraint.polynomial, p, dec)?;
                if at != Sign::Zero {
                    constraints.push(LocalConstraint::Fixed(at));
                    continue;
                }
                let (support, gradient) = geometry::support_geometry(&constraint.polynomial, dec)?;
                let [x, y, z] = gradient.map(|g| g.eval_real(&xyz));
                let gradient = Vector3::new(x, y, z);
                constraints.push(if support.equals(&face.support, dec)? {
                    LocalConstraint::Own(geometry::sign(&gradient.dot(&normal), dec)?)
                } else if let Some((sign, escape)) =
                    tangent_bound(&face.surface, constraint.surface(), &xyz, dec)?
                {
                    LocalConstraint::TangentBound { sign, escape }
                } else {
                    LocalConstraint::Transverse(gradient)
                });
            }
            leaves.push(LocalLeaf::Cell(constraints));
        }
        let mut plus = Vec::with_capacity(self.instructions.len());
        let mut minus = Vec::with_capacity(self.instructions.len());
        let mut domains = Vec::new();
        for tangent in tangents {
            plus.clear();
            minus.clear();
            for instruction in &self.instructions {
                let pair = match *instruction {
                    Instruction::Empty => [Location::Outside; 2],
                    Instruction::Leaf(leaf) => [
                        leaves[leaf].side(&tangent, true, dec)?,
                        leaves[leaf].side(&tangent, false, dec)?,
                    ],
                    Instruction::Boolean { kind, a, b } => [
                        crate::combine(kind, plus[a], plus[b]),
                        crate::combine(kind, minus[a], minus[b]),
                    ],
                };
                plus.push(pair[0]);
                minus.push(pair[1]);
            }
            let source = self.leaf_registers[face.leaf].expect("retained source leaf");
            if !matches!(
                (plus[source], minus[source]),
                (Location::Inside, Location::Outside) | (Location::Outside, Location::Inside)
            ) || self.domain_state(face.domain, &plus, &minus, &mut domains) != FaceState::OnFace
            {
                continue;
            }
            let outward = match (plus[self.root], minus[self.root]) {
                (Location::Outside, Location::Inside) => normal.clone(),
                (Location::Inside, Location::Outside) => normal.scaled(&Real::from_i64(-1)),
                _ => continue,
            };
            // On the smooth source carrier q(t)=p+t*T+O(t²), every transverse
            // sign is fixed for small positive t. Unknown tangent signs were
            // allowed both values in the Boolean evaluation above. Thus an
            // open source patch, not merely its tangent ray, reaches this point.
            return Ok(Some(FaceQuery {
                state: FaceState::OnFace,
                direction: crate::unit_direction(outward, dec)?,
            }));
        }
        Ok(None)
    }

    fn contact_face_query(
        &self,
        id: FaceId,
        point: &RPoint3,
        dec: &mut Decider,
    ) -> Result<FaceQuery, Refusal> {
        let face = &self.faces[id.0]; // The ordinary query validated this handle.
        let p = [point.x.clone(), point.y.clone(), point.z.clone()];
        let source = &self.leaves[face.leaf];
        if source.location(&p, dec)? != Location::Boundary
            || source.on_face(&face.origin, self.definition_role(id.0), &p, dec)? != Some(true)
        {
            return Ok(unknown());
        }
        // A carrier shared by different source roles can have authored trim
        // partitions absent from the solid's support arrangement. Never borrow
        // another partition's open patch as an incidence proof for this name.
        for (j, other) in self.faces.iter().enumerate() {
            if j != id.0
                && other.leaf == face.leaf
                && self.definition_role(j) != self.definition_role(id.0)
                && crate::same_support(face, other, dec)?
            {
                return Ok(unknown());
            }
        }
        let mut active: Vec<(RealPoly, Vector3)> = Vec::new();
        let mut planar = true;
        for (i, leaf) in self.leaves.iter().enumerate() {
            if self.leaf_registers[i].is_none() {
                continue;
            }
            for constraint in leaf.constraints() {
                dec.check_deadline()?;
                if crate::polynomial_sign(&constraint.polynomial, &p, dec)? != Sign::Zero {
                    continue;
                }
                let (support, gradient) = geometry::support_geometry(&constraint.polynomial, dec)?;
                let mut duplicate = false;
                for (old, _) in &active {
                    duplicate |= support.equals(old, dec)?;
                }
                if duplicate {
                    continue;
                }
                planar &= constraint
                    .polynomial
                    .terms()
                    .iter()
                    .all(|(m, _)| m.exponents().iter().sum::<u32>() <= 1);
                let xyz = p.clone().map(Real::from);
                let [x, y, z] = gradient.map(|g| g.eval_real(&xyz));
                let normal = Vector3::new(x, y, z);
                if geometry::sign(&normal.norm_squared(), dec)? == Sign::Zero {
                    return Ok(unknown());
                }
                active.push((support, normal));
            }
        }
        if active.len() > 2 && !planar {
            return self.tangent_contact_query(id, &p, &active, dec);
        }
        let mut own = None;
        for (i, (support, normal)) in active.iter().enumerate() {
            if support.equals(&face.support, dec)? {
                own = Some(i);
            }
            for (_, other) in &active[i + 1..] {
                if geometry::sign(&normal.cross(other).norm_squared(), dec)? == Sign::Zero {
                    return self.tangent_contact_query(id, &p, &active, dec);
                }
            }
        }
        let Some(own) = own else {
            return Ok(unknown());
        };
        let n = crate::normal(face, &p);
        let others: Vec<_> = active
            .into_iter()
            .enumerate()
            .filter_map(|(i, (_, n))| (i != own).then_some(n))
            .collect();
        let sectors = tangent_sectors(&n, &others, dec)?;
        if sectors.is_empty() {
            return Ok(unknown());
        }
        let mut orientation = 0;
        let mut incomplete = false;
        let mut tape = Vec::new();
        let mut opposite = Vec::new();
        let mut domains = Vec::new();
        for tangent in sectors {
            dec.check_deadline()?;
            let Some(e) = epsilon(&others, &tangent, &n, dec)? else {
                incomplete = true;
                continue;
            };
            let offset = n.scaled(&Real::from(e));
            let directions = [tangent.plus(&offset), tangent.minus(&offset)];
            let source_sides = [
                source.ray(&p, &directions[0], dec)?,
                source.ray(&p, &directions[1], dec)?,
            ];
            if matches!(
                source_sides,
                [Location::Inside, Location::Inside] | [Location::Outside, Location::Outside]
            ) {
                continue;
            }
            if !matches!(
                source_sides,
                [Location::Inside, Location::Outside] | [Location::Outside, Location::Inside]
            ) {
                incomplete = true;
                continue;
            }
            let sides = [
                self.ray(&p, &directions[0], &mut tape, dec)?,
                self.ray(&p, &directions[1], &mut opposite, dec)?,
            ];
            match self.domain_state(face.domain, &tape, &opposite, &mut domains) {
                FaceState::NotOnFace => continue,
                FaceState::Indeterminate => {
                    incomplete = true;
                    continue;
                }
                FaceState::OnFace => (),
            }
            match sides {
                [Location::Outside, Location::Inside] => orientation |= 1,
                [Location::Inside, Location::Outside] => orientation |= 2,
                [Location::Inside, Location::Inside] | [Location::Outside, Location::Outside] => (),
                _ => incomplete = true,
            }
        }
        if orientation == 0 {
            return Ok(if incomplete {
                unknown()
            } else {
                FaceQuery {
                    state: FaceState::NotOnFace,
                    direction: Direction::Undefined,
                }
            });
        }
        let direction = if incomplete || orientation == 3 {
            Direction::Ambiguous
        } else {
            let n = if orientation == 1 {
                n
            } else {
                n.scaled(&Real::from_i64(-1))
            };
            crate::unit_direction(n, dec)?
        };
        Ok(FaceQuery {
            state: FaceState::OnFace,
            direction,
        })
    }

    // Try transverse-trim dominance before the more expensive curved germs.
    // A second-order blow-up separates a round from its tangent planar
    // continuation. Strict signs of every other support certify an open patch
    // between the two paths by the implicit function theorem. An isolated
    // buried tangency has no such source-and-result boundary witness.
    fn tangent_contact_query(
        &self,
        id: FaceId,
        p: &[Rational; 3],
        active: &[(RealPoly, Vector3)],
        dec: &mut Decider,
    ) -> Result<FaceQuery, Refusal> {
        if let Some(query) = self.dominated_contact(id, p, active, dec)? {
            return Ok(query);
        }

        let face = &self.faces[id.0];
        let n = crate::normal(face, p);
        let xyz = p.clone().map(Real::from);
        let mut others = Vec::new();
        let mut graphs = Vec::new();
        let mut own = None;
        for (support, normal) in active {
            if geometry::sign(&normal.cross(&n).norm_squared(), dec)? != Sign::Zero {
                others.push(normal.clone());
            } else {
                let gradient = std::array::from_fn::<_, 3, _>(|i| support.partial(i));
                let g = Vector3::new(
                    gradient[0].eval_real(&xyz),
                    gradient[1].eval_real(&xyz),
                    gradient[2].eval_real(&xyz),
                );
                if support.equals(&face.support, dec)? {
                    own = Some(graphs.len());
                }
                graphs.push((gradient, g.dot(&n)));
            }
        }
        let Some(own) = own else {
            return Ok(unknown());
        };
        let mut sectors = tangent_sectors(&n, &others, dec)?;
        if others.is_empty() {
            sectors.clear();
            for axis in 0..3 {
                let [x, y, z] =
                    std::array::from_fn(|i| Rational::from_i64(if i == axis { 1 } else { 0 }));
                let tangent = n.cross(&RVector3::new(x, y, z).to_real());
                if geometry::sign(&tangent.norm_squared(), dec)? != Sign::Zero {
                    sectors.push(tangent);
                }
            }
        }
        let mut orientation = 0;
        let mut tape = Vec::new();
        let mut opposite = Vec::new();
        let mut domains = Vec::new();
        let mut heights = Vec::with_capacity(graphs.len());
        for tangent in sectors {
            let d = [&tangent.x, &tangent.y, &tangent.z];
            heights.clear();
            for (gradient, slope) in &graphs {
                let mut quadratic = Real::zero();
                for i in 0..3 {
                    for j in 0..3 {
                        quadratic =
                            quadratic + d[i] * d[j] * gradient[i].partial(j).eval_real(&xyz);
                    }
                }
                heights.push(-quadratic / (slope * &Real::from_i64(2)));
            }
            let mut e = Real::one();
            let mut separated = true;
            for (i, height) in heights.iter().enumerate() {
                if i == own {
                    continue;
                }
                let difference = height - &heights[own];
                if geometry::sign(&difference, dec)? == Sign::Zero {
                    separated = false;
                    break;
                }
                let distance = magnitude(&difference, dec)? / Real::from_i64(3);
                if geometry::sign(&(&distance - &e), dec)? == Sign::Negative {
                    e = distance;
                }
            }
            if !separated {
                continue;
            }
            let plus = n.scaled(&(&heights[own] + &e));
            let minus = n.scaled(&(&heights[own] - &e));
            let sides = [
                self.parabolic_ray(p, &tangent, &plus, &mut tape, dec)?,
                self.parabolic_ray(p, &tangent, &minus, &mut opposite, dec)?,
            ];
            let slot = self.leaf_registers[face.leaf].expect("a live contact source");
            if !matches!(
                [tape[slot], opposite[slot]],
                [Location::Inside, Location::Outside] | [Location::Outside, Location::Inside]
            ) {
                continue;
            }
            if self.domain_state(face.domain, &tape, &opposite, &mut domains) != FaceState::OnFace {
                continue;
            }
            match sides {
                [Location::Outside, Location::Inside] => orientation |= 1,
                [Location::Inside, Location::Outside] => orientation |= 2,
                _ => (),
            }
        }
        if orientation == 0 {
            return Ok(unknown());
        }
        // A single-polarity source cannot reverse orientation on an unvisited
        // tangent sector. Mixed-polarity reuse retains an ambiguous normal.
        let mut polarity = vec![0u8; self.instructions.len()];
        polarity[self.root] = 1;
        for node in (0..self.instructions.len()).rev() {
            if let crate::Instruction::Boolean { kind, a, b } = self.instructions[node] {
                let sign = polarity[node];
                polarity[a] |= sign;
                polarity[b] |= if kind == catallax_recipe::BooleanKind::Subtract {
                    ((sign & 1) << 1) | ((sign & 2) >> 1)
                } else {
                    sign
                };
            }
        }
        let slot = self.leaf_registers[face.leaf].expect("a live contact source");
        let direction = if orientation == polarity[slot] && orientation != 3 {
            crate::unit_direction(
                if orientation == 1 {
                    n
                } else {
                    n.scaled(&Real::from_i64(-1))
                },
                dec,
            )?
        } else {
            Direction::Ambiguous
        };
        Ok(FaceQuery {
            state: FaceState::OnFace,
            direction,
        })
    }

    fn parabolic_ray(
        &self,
        p: &[Rational; 3],
        tangent: &Vector3,
        quadratic: &Vector3,
        tape: &mut Vec<Location>,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        tape.clear();
        for instruction in &self.instructions {
            let location = match *instruction {
                crate::Instruction::Empty => Location::Outside,
                crate::Instruction::Leaf(i) => {
                    let leaf = &self.leaves[i];
                    let location = leaf.location(p, dec)?;
                    if matches!(location, Location::Inside | Location::Outside) {
                        location
                    } else if matches!(leaf, geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_))
                    {
                        geometry::parabolic_constraints(
                            leaf.constraints(),
                            p,
                            tangent,
                            quadratic,
                            dec,
                        )?
                    } else {
                        Location::Indeterminate
                    }
                }
                crate::Instruction::Boolean { kind, a, b } => {
                    crate::combine(kind, tape[a], tape[b])
                }
            };
            tape.push(location);
        }
        Ok(tape[self.root])
    }
}
