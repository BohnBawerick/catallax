//! One unchanged raised cylinder on a stationary cap of a convex planar draft cell.
//!
//! Both kernels share this definition math. Neither passes topology through it.

use crate::draft::{
    DraftDefinition,
    bored::{self, Bore, BoredCell},
    cell::PlanarCell,
    sign,
};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cylinder, Plane, Point3};

/// One finite outward cylinder and its end disk, both unchanged by the draft.
#[derive(Clone, Debug)]
pub struct Boss {
    /// Unchanged wall carrier.
    pub cylinder: Cylinder,
    /// Complete source names of the wall.
    pub origins: Vec<Origin>,
    /// Unchanged end disk's carrier, outward from material.
    pub top: Plane,
    /// Complete source names of the end disk.
    pub top_origins: Vec<Origin>,
}

/// A convex planar envelope with one boss joined across a whole disk of a cap.
#[derive(Clone, Debug)]
pub struct BossedCell {
    /// Drafted base material, whose every point lies below the shoulder.
    pub outer: PlanarCell,
    /// The raised material above that cap.
    pub boss: Boss,
    /// Support index of the cap carrying the boss.
    pub shoulder: usize,
    /// The bores through `outer`, each strictly apart from the boss unless `through`.
    pub bored: Option<BoredCell>,
    /// Whether the bore runs on through the boss, strictly inside its footprint, to
    /// open on its end instead of the shoulder.
    pub through: bool,
}

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-arrangement", detail)
}

fn value(plane: &Plane, point: &Point3) -> Real {
    point.minus(&plane.origin).dot(&plane.normal)
}

/// Finds the cap carrying the boss and proves the boss's footprint strictly
/// inside every other halfspace at that cap. Only the footprint matters: a
/// convex envelope lies below its cap, so the planes' extensions above it bound
/// nothing and may cross the raised material.
pub fn certify_footprint<'a>(
    planes: impl Iterator<Item = &'a Plane> + Clone,
    boss: &Boss,
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    let Boss { cylinder, top, .. } = boss;
    if sign(&top.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero {
        return Err(unsupported("a preserved boss needs an axial end disk"));
    }
    let mut shoulder = None;
    for (index, plane) in planes.clone().enumerate() {
        if sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)? == Sign::Zero
            && sign(&plane.normal.dot(&top.normal), dec)? == Sign::Positive
            && shoulder.replace(index).is_some()
        {
            return Err(unsupported("a preserved boss needs exactly one shoulder"));
        }
    }
    let shoulder = shoulder.ok_or_else(|| unsupported("a preserved boss needs a shoulder"))?;
    let cap = planes.clone().nth(shoulder).expect("indexed above");
    let centre = bored::axis_point(cylinder, cap);
    if sign(&value(cap, &bored::axis_point(cylinder, top)), dec)? != Sign::Positive {
        return Err(unsupported("a preserved boss needs height above its cap"));
    }
    for (index, wall) in planes.enumerate() {
        if index != shoulder && !bored::ring_clear(wall, cylinder, &centre, dec)? {
            return Err(unsupported(
                "the boss footprint touches or crosses the edge of its cap",
            ));
        }
    }
    Ok(shoulder)
}

impl BossedCell {
    /// Certifies the boss and any bores against the replaced envelope before any
    /// publication. Each may clear every wall and still meet the other, so the two
    /// footprints are proved apart here.
    pub fn new(
        outer: PlanarCell,
        boss: Boss,
        bores: Vec<Bore>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Self::arranged(outer, boss, bores, false, definition, dec)
    }

    /// The same certification for a bore that runs on through the boss. Its footprint
    /// is proved strictly inside the boss's, so the bore never meets the boss's wall
    /// and the shoulder, whose ring stays outside it, needs no second exclusion.
    pub fn through(
        outer: PlanarCell,
        boss: Boss,
        bores: Vec<Bore>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Self::arranged(outer, boss, bores, true, definition, dec)
    }

    fn arranged(
        outer: PlanarCell,
        boss: Boss,
        bores: Vec<Bore>,
        through: bool,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        bored::check_axial(&boss.cylinder, definition, dec)?;
        let shoulder = certify_footprint(outer.supports.iter().map(|s| &s.plane), &boss, dec)?;
        let bored = if bores.is_empty() {
            None
        } else {
            let cell = BoredCell::new(outer.clone(), bores, definition, dec)?;
            for bore in &cell.bores {
                if through {
                    bored::certify_nested(&boss.cylinder, &bore.cylinder, dec)?;
                } else {
                    bored::certify_disjoint(&boss.cylinder, &bore.cylinder, dec)?;
                }
            }
            Some(cell)
        };
        Ok(Self {
            outer,
            boss,
            shoulder,
            bored,
            through,
        })
    }

    /// Exact strict material membership. The joined disk is interior: the base
    /// lies under it and the boss over it. A bore running through keeps cutting
    /// the boss above that disk.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let base = match &self.bored {
            Some(cell) => cell.contains_interior(point, dec)?,
            None => self.outer.contains_interior(point, dec)?,
        };
        if base {
            return Ok(true);
        }
        let cylinder = &self.boss.cylinder;
        let radial = point.minus(&cylinder.origin).cross(&cylinder.axis);
        let excess = radial.norm_squared() - &cylinder.radius_sq * &cylinder.axis.norm_squared();
        if sign(&excess, dec)? != Sign::Negative
            || sign(
                &value(&self.outer.supports[self.shoulder].plane, point),
                dec,
            )? == Sign::Negative
            || sign(&value(&self.boss.top, point), dec)? != Sign::Negative
        {
            return Ok(false);
        }
        match &self.bored {
            Some(cell) if self.through => cell.clear_of_bores(point, dec),
            _ => Ok(true),
        }
    }

    /// A relatively open neutral section, through the base or across the boss.
    /// The neutral plane lies across the axis, so it cuts or holds a whole disk.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.outer.anchored(neutral, dec)?
            || sign(
                &value(&self.outer.supports[self.shoulder].plane, &neutral.origin),
                dec,
            )? == Sign::Positive
                && sign(&value(&self.boss.top, &neutral.origin), dec)? != Sign::Positive)
    }
}
