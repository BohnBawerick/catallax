//! Cap-hole distance bounds for the separated through-bore shell family.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use catallax_features::{
    draft::bored::{axis_point, certify_bores},
    polyhedron::ConvexCell,
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::BooleanKind;
use catallax_shapes::{
    Cylinder, Plane, Point3, Surface,
    distance::{Bounds, Options},
};

use super::Prepared;
use crate::{Instruction, Plan, geometry};

type Caps = BTreeMap<usize, Arc<[Cylinder]>>;

struct Skin {
    base: usize,
    cell: ConvexCell,
    bores: Vec<Cylinder>,
    caps: [usize; 2],
}

fn coplanar(a: &Plane, b: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(
        geometry::sign(&a.normal.cross(&b.normal).norm_squared(), dec)? == Sign::Zero
            && geometry::sign(&a.origin.minus(&b.origin).dot(&b.normal), dec)? == Sign::Zero,
    )
}

impl Plan {
    fn bored_distance_skin(
        &self,
        mut node: usize,
        dec: &mut Decider,
    ) -> Result<Option<Skin>, Refusal> {
        let mut tools = Vec::new();
        while let Instruction::Boolean {
            kind: BooleanKind::Subtract,
            a,
            b,
        } = self.instructions[node]
        {
            tools.push(b);
            node = a;
        }
        if tools.is_empty() {
            return Ok(None);
        }
        let cell = match self.planar_skin(node, dec) {
            Ok(cell) if !cell.is_empty() => cell,
            Ok(_) => return Ok(None),
            Err(error) if error.code.starts_with("distance.") => return Ok(None),
            Err(error) => return Err(error),
        };
        let mut bores = Vec::new();
        for tool in tools {
            match self.distance_tool_bores(tool, &cell, dec) {
                Ok(cylinders) => bores.extend(cylinders),
                Err(error) if error.code.starts_with("draft.") => return Ok(None),
                Err(error) => return Err(error),
            }
        }
        let cylinders: Vec<_> = bores.iter().collect();
        let caps = match certify_bores(cell.faces.iter().map(|face| &face.plane), &cylinders, dec) {
            Ok(caps) => caps,
            Err(error) if error.code.starts_with("draft.") => return Ok(None),
            Err(error) => return Err(error),
        };
        Ok(Some(Skin {
            base: node,
            cell,
            bores,
            caps,
        }))
    }

    pub(super) fn prepare_bored_caps(&self, dec: &mut Decider) -> Result<Caps, Refusal> {
        let mut caps = Caps::new();
        let Instruction::Boolean {
            kind: BooleanKind::Subtract,
            a,
            b,
        } = self.instructions[self.root]
        else {
            return Ok(caps);
        };
        let (Some(outer), Some(inner)) = (
            self.bored_distance_skin(a, dec)?,
            self.bored_distance_skin(b, dec)?,
        ) else {
            return Ok(caps);
        };
        if outer.bores.len() != inner.bores.len() {
            return Ok(caps);
        }
        // The cavity lies within the original hull. Every original hole lies
        // within its corresponding enlarged hole, so the cavity adds no material.
        for p in inner.cell.vertices() {
            for face in &outer.cell.faces {
                if geometry::sign(&p.minus(&face.plane.origin).dot(&face.plane.normal), dec)?
                    == Sign::Positive
                {
                    return Ok(caps);
                }
            }
        }
        for original in &outer.bores {
            let mut matched = false;
            for expanded in &inner.bores {
                let (a, b) = (original, expanded);
                if geometry::sign(&a.axis.cross(&b.axis).norm_squared(), dec)? == Sign::Zero
                    && geometry::sign(
                        &a.origin.minus(&b.origin).cross(&a.axis).norm_squared(),
                        dec,
                    )? == Sign::Zero
                    && geometry::sign(&(&b.radius_sq - &a.radius_sq), dec)? != Sign::Negative
                {
                    matched = true;
                    break;
                }
            }
            if !matched {
                return Ok(caps);
            }
        }
        for (skin, other) in [(&outer, &inner), (&inner, &outer)] {
            let mut leaves = BTreeSet::new();
            let mut stack = vec![skin.base];
            while let Some(node) = stack.pop() {
                match self.instructions[node] {
                    Instruction::Leaf(leaf) => {
                        leaves.insert(leaf);
                    }
                    Instruction::Boolean { a, b, .. } => stack.extend([a, b]),
                    Instruction::Empty => (),
                }
            }
            let holes: Arc<[Cylinder]> = skin.bores.clone().into();
            for cap in skin.caps {
                let plane = &skin.cell.faces[cap].plane;
                let mut open = false;
                for other_cap in other.caps {
                    open |= coplanar(plane, &other.cell.faces[other_cap].plane, dec)?;
                }
                // An open cap is a rim with a polygonal cavity opening, not the
                // polygon-minus-separated-circles domain certified here.
                if open {
                    continue;
                }
                for (id, face) in self.faces.iter().enumerate() {
                    if face.in_result
                        && leaves.contains(&face.leaf)
                        && let Surface::Plane(support) = &face.surface
                        && coplanar(plane, support, dec)?
                    {
                        caps.insert(id, holes.clone());
                    }
                }
            }
        }
        Ok(caps)
    }
}

impl Prepared<'_> {
    pub(super) fn bored_cap_distance(
        &self,
        id: usize,
        point: &Point3,
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Option<Bounds>, Refusal> {
        let Some(holes) = self.bored_caps.get(&id) else {
            return Ok(None);
        };
        let Surface::Plane(plane) = &self.plan.faces[id].surface else {
            unreachable!()
        };
        let projection = point.plus(&plane.normal.scaled(
            &(-point.minus(&plane.origin).dot(&plane.normal) / plane.normal.norm_squared()),
        ));
        for hole in holes.iter() {
            let center = axis_point(hole, plane);
            let radial = projection.minus(&center);
            let rho_sq = radial.norm_squared();
            if geometry::sign(&(&rho_sq - &hole.radius_sq), dec)? == Sign::Positive {
                continue;
            }
            let rho = rho_sq.sqrt();
            let radius = hole.radius_sq.sqrt();
            let distance =
                (point.minus(&projection).norm_squared() + (&radius - &rho).square()).sqrt();
            let mut bits = 128;
            let lower = loop {
                let bound = dec.enclose_budgeted_at(&distance, bits)?;
                if bound.hi() - bound.lo() <= options.width.mul_pow2(-2) {
                    break bound.lo().max(Rational::zero());
                }
                bits = bits.checked_mul(2).ok_or_else(|| {
                    Refusal::new(
                        "numbers.refinement-budget",
                        "cap-hole distance precision overflow",
                    )
                })?;
            };
            let direction = if geometry::sign(&rho_sq, dec)? == Sign::Zero {
                hole.ref_dir
                    .scaled(&(Real::one() / hole.ref_dir.norm_squared().sqrt()))
            } else {
                radial.scaled(&(&Real::one() / &rho))
            };
            let mut step = options.width.mul_pow2(-3);
            // A point just outside the complete hole circle supplies an admitted
            // upper witness. Strict wall/hole clearance guarantees an open annulus;
            // no point in the hole can lower the exact circular lower certificate.
            for _ in 0..options.max_cells {
                let witness =
                    center.plus(&direction.scaled(&(&radius + &Real::from(step.clone()))));
                if self.member(id, &witness, dec)? {
                    let upper = dec
                        .enclose_budgeted_at(&point.minus(&witness).norm_squared().sqrt(), bits)?
                        .hi();
                    if &upper - &lower <= options.width {
                        return Ok(Some(Bounds {
                            lower,
                            upper: Some(upper),
                            witnesses: Some([point.clone(), witness]),
                            converged: true,
                            cells: 0,
                        }));
                    }
                }
                step = step.mul_pow2(-1);
            }
            return Ok(None);
        }
        Ok(None)
    }
}
