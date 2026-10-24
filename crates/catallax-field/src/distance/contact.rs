//! Closed named-patch incidence at certified transverse support contacts.
//!
//! Two independent smooth supports have local coordinates given by their
//! defining polynomials. Their open sign sectors therefore exhaust a
//! neighbourhood, including curved approaches to the contact. The same is
//! true for an exact plane arrangement. Tangencies and larger curved fans are
//! deliberately left unresolved rather than inferred from tangent samples.

use crate::{Direction, FaceId, FaceQuery, FaceState, Location, Plan, Workspace, geometry};
use catallax_numbers::{Decider, Rational, Real, RealPoly, Refusal, Sign};
use catallax_shapes::{RPoint3, RVector3, Vector3};

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

impl Plan {
    /// Tests the closure of a named surviving face, refining the ordinary face
    /// query with exact transverse-contact sectors. Tangencies, singularities,
    /// ambiguous source partitions and larger curved fans remain indeterminate.
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
            return Ok(unknown());
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
                    return Ok(unknown());
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
                self.ray(&p, &directions[1], &mut tape, dec)?,
            ];
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
