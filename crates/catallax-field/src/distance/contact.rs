//! Closed named-patch incidence at certified support contacts.
//!
//! Transverse supports admit local polynomial sign coordinates. At a tangent
//! seam, a transverse trim can also make every other tangent support irrelevant.
//! In that case universal Boolean bounds prove an open source patch approaches
//! the contact. Tangent supports are never assigned signs from tangent samples.

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
    /// Tests named-patch closure using transverse contact sectors or a source
    /// patch whose transverse trims dominate all unresolved tangent supports.
    /// Singularities and contacts lacking either proof remain indeterminate.
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
            return Ok(self
                .dominated_contact(id, &p, &active, dec)?
                .unwrap_or_else(unknown));
        }
        let mut own = None;
        for (i, (support, normal)) in active.iter().enumerate() {
            if support.equals(&face.support, dec)? {
                own = Some(i);
            }
            for (_, other) in &active[i + 1..] {
                if geometry::sign(&normal.cross(other).norm_squared(), dec)? == Sign::Zero {
                    // Distinct tangent supports may enclose thin curved sectors
                    // invisible in their common tangent plane.
                    return Ok(self
                        .dominated_contact(id, &p, &active, dec)?
                        .unwrap_or_else(unknown));
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
}
