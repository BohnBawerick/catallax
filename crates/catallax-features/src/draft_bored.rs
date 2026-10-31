//! Unchanged pull-parallel bores strictly inside a convex planar draft cell.
//!
//! Both kernels share this definition math. Neither passes topology through it.

use crate::draft::{DraftDefinition, cell::PlanarCell, sign};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cone, Cylinder, Plane, Point3, Vector3};

/// The finite cone that opens a countersunk bore on one cap. It shares one complete
/// ring with the bore's wall, so the hole has no ledge between them.
#[derive(Clone, Debug)]
pub struct Head {
    /// Unchanged cone carrier, on the nappe that widens toward its cap.
    pub cone: Cone,
    /// Complete source names of the cone.
    pub origins: Vec<Origin>,
}

/// One finite bore's unchanged wall.
#[derive(Clone, Debug)]
pub struct Bore {
    /// Unchanged bore wall carrier.
    pub cylinder: Cylinder,
    /// Complete source names of the bore wall.
    pub origins: Vec<Origin>,
}

/// A convex planar envelope with finite bores removed from cap to cap.
/// Every bore's support and both ring trims never change.
#[derive(Clone, Debug)]
pub struct BoredCell {
    /// Filled outer material, before any bore is removed.
    pub outer: PlanarCell,
    /// Mutually separated bores, all between the same two caps.
    pub bores: Vec<Bore>,
    /// The two cap support indices, in ascending order.
    pub caps: [usize; 2],
    /// The cone of a countersunk bore, which stands alone, and the cap it opens on.
    pub head: Option<(Head, usize)>,
}

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-arrangement", detail)
}

fn value(plane: &Plane, point: &Point3) -> Real {
    point.minus(&plane.origin).dot(&plane.normal)
}

/// The cylinder axis point on a plane perpendicular to that axis.
pub fn axis_point(cylinder: &Cylinder, plane: &Plane) -> Point3 {
    let height = -&(value(plane, &cylinder.origin) / plane.normal.dot(&cylinder.axis));
    cylinder.origin.plus(&cylinder.axis.scaled(&height))
}

/// Whether a circle about the axis at `centre` lies strictly inside a halfspace.
/// The circular maximum is the axis value plus radius times the normal's radial
/// length, compared on squares so the proof needs no root.
pub fn ring_clear(
    wall: &Plane,
    cylinder: &Cylinder,
    centre: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    circle_clear(wall, &cylinder.axis, &cylinder.radius_sq, centre, dec)
}

fn circle_clear(
    wall: &Plane,
    axis: &Vector3,
    radius_sq: &Real,
    centre: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let along = wall.normal.dot(axis);
    let radial = wall.normal.norm_squared() - along.square() / axis.norm_squared();
    let offset = value(wall, centre);
    Ok(sign(&offset, dec)? == Sign::Negative
        && sign(&(offset.square() - radius_sq * &radial), dec)? == Sign::Positive)
}

/// An unselected wall has no neutral curve, but this bounded family keeps its
/// axial mouths parallel to the neutral section: pull, axis and normal agree.
pub fn check_axial(
    cylinder: &Cylinder,
    definition: &DraftDefinition,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let parallel = |a: &Vector3, dec: &mut Decider| {
        Ok::<_, Refusal>(sign(&a.cross(&definition.pull).norm_squared(), dec)? == Sign::Zero)
    };
    if !parallel(&cylinder.axis, dec)? || !parallel(&definition.neutral.normal, dec)? {
        return Err(unsupported(
            "a preserved wall needs its axis and the neutral normal along pull",
        ));
    }
    Ok(())
}

/// Finds the two opposed caps and proves the finite cylinder between them is
/// strictly inside every other halfspace. A wall value is affine along the axis,
/// so strict inequalities at both rings cover the whole height.
pub fn certify_clear<'a>(
    planes: impl Iterator<Item = &'a Plane> + Clone,
    cylinder: &Cylinder,
    dec: &mut Decider,
) -> Result<[usize; 2], Refusal> {
    let mut lower = None;
    let mut upper = None;
    for (index, plane) in planes.clone().enumerate() {
        if sign(&plane.normal.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero {
            continue;
        }
        let side = if sign(&plane.normal.dot(&cylinder.axis), dec)? == Sign::Positive {
            &mut upper
        } else {
            &mut lower
        };
        if side.replace(index).is_some() {
            return Err(unsupported(
                "a preserved bore needs exactly two opposed caps",
            ));
        }
    }
    let (Some(lower), Some(upper)) = (lower, upper) else {
        return Err(unsupported(
            "a preserved bore needs exactly two opposed caps",
        ));
    };
    let caps = [lower, upper];
    let plane = |index: usize| planes.clone().nth(index).expect("indexed above");
    let centres = caps.map(|i| axis_point(cylinder, plane(i)));
    if sign(&centres[1].minus(&centres[0]).dot(&cylinder.axis), dec)? != Sign::Positive {
        return Err(unsupported("preserved bore caps have no axial interval"));
    }
    for (index, wall) in planes.enumerate() {
        if caps.contains(&index) {
            continue;
        }
        for centre in &centres {
            if !ring_clear(wall, cylinder, centre, dec)? {
                return Err(unsupported(
                    "the bore touches or crosses a wall of its planar envelope",
                ));
            }
        }
    }
    Ok(caps)
}

/// The axis point where a cone reaches the radius of a coaxial cylinder, on its nappe.
pub fn seam(cylinder: &Cylinder, cone: &Cone) -> Point3 {
    let distance = cylinder.radius_sq.sqrt() / (&cone.slope * &cone.axis.norm_squared().sqrt());
    cone.apex.plus(&cone.axis.scaled(&distance))
}

/// Proves a cone opens this bore on one of its two caps, which it returns: the two
/// are coaxial, and their shared ring lies strictly between the caps at a positive
/// radius. Every wall value is affine along each band and `certify_clear` has
/// cleared the cylinder, so only the wider mouth ring is left to prove strict. A
/// cone's extension past either ring stays out of the proof: no material meets it.
pub fn certify_head<'a>(
    planes: impl Iterator<Item = &'a Plane> + Clone,
    cylinder: &Cylinder,
    cone: &Cone,
    caps: [usize; 2],
    dec: &mut Decider,
) -> Result<usize, Refusal> {
    if sign(&cone.axis.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero
        || sign(
            &cone
                .apex
                .minus(&cylinder.origin)
                .cross(&cylinder.axis)
                .norm_squared(),
            dec,
        )? != Sign::Zero
    {
        return Err(unsupported(
            "a countersink needs a cone coaxial with its bore",
        ));
    }
    let plane = |index: usize| planes.clone().nth(index).expect("indexed caps");
    let reach = |cap: usize| axial_reach(cone, &axis_point(cylinder, plane(cap)));
    let (mouth, far) = if sign(&(reach(caps[0]) - reach(caps[1])), dec)? == Sign::Positive {
        (caps[0], caps[1])
    } else {
        (caps[1], caps[0])
    };
    let mouth_sq = nappe_radius_sq(cone, &reach(mouth));
    let beyond = |radius_sq: Real| radius_sq - cylinder.radius_sq.clone();
    if sign(&reach(mouth), dec)? != Sign::Positive
        || sign(&beyond(mouth_sq.clone()), dec)? != Sign::Positive
        || (sign(&reach(far), dec)? == Sign::Positive
            && sign(&beyond(nappe_radius_sq(cone, &reach(far))), dec)? != Sign::Negative)
    {
        return Err(unsupported(
            "a countersink's shared ring must lie between the bore's caps",
        ));
    }
    let centre = axis_point(cylinder, plane(mouth));
    for (index, wall) in planes.clone().enumerate() {
        if !caps.contains(&index) && !circle_clear(wall, &cylinder.axis, &mouth_sq, &centre, dec)? {
            return Err(unsupported(
                "the countersink touches or crosses a wall of its planar envelope",
            ));
        }
    }
    Ok(mouth)
}

/// The signed axial distance from a cone's apex, scaled by its axis's length.
fn axial_reach(cone: &Cone, point: &Point3) -> Real {
    point.minus(&cone.apex).dot(&cone.axis)
}

/// The squared nappe radius at a scaled axial distance.
fn nappe_radius_sq(cone: &Cone, reach: &Real) -> Real {
    (&cone.slope * reach).square() / cone.axis.norm_squared()
}

/// Whether a point lies in or on a cone's retained nappe, however far it runs.
fn in_nappe(cone: &Cone, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    let reach = axial_reach(cone, point);
    if sign(&reach, dec)? != Sign::Positive {
        return Ok(false);
    }
    let radial =
        point.minus(&cone.apex).cross(&cone.axis).norm_squared() / cone.axis.norm_squared();
    Ok(sign(&(radial - nappe_radius_sq(cone, &reach)), dec)? != Sign::Positive)
}

/// The squared distance between two parallel axes.
fn axis_offset(a: &Cylinder, b: &Cylinder) -> Real {
    b.origin.minus(&a.origin).cross(&a.axis).norm_squared() / a.axis.norm_squared()
}

/// Proves two walls strictly separated. Parallel walls are nearest along their
/// common normal, so the axis distance `d` must exceed the radii's sum. Squaring
/// `d^2-a^2-b^2 > 2ab` twice keeps the proof exact without a root. An exact
/// contact refuses, and a gap of any positive size stands unchanged.
pub fn certify_disjoint(a: &Cylinder, b: &Cylinder, dec: &mut Decider) -> Result<(), Refusal> {
    if sign(&a.axis.cross(&b.axis).norm_squared(), dec)? != Sign::Zero {
        return Err(unsupported("preserved walls need parallel axes"));
    }
    let spare = axis_offset(a, b) - a.radius_sq.clone() - b.radius_sq.clone();
    let reach = Real::from_i64(4) * (&a.radius_sq * &b.radius_sq);
    if sign(&spare, dec)? != Sign::Positive
        || sign(&(spare.square() - reach), dec)? != Sign::Positive
    {
        return Err(unsupported("preserved walls touch or cross each other"));
    }
    Ok(())
}

/// Whether `inner`'s axis lies strictly inside the wall of `outer`. A wall about such
/// an axis can only nest there, and one about any other axis can only stand apart, so
/// this names which of the two proofs a pair of parallel walls has to pass.
pub fn axis_inside(outer: &Cylinder, inner: &Cylinder, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(&(axis_offset(outer, inner) - outer.radius_sq.clone()), dec)? == Sign::Negative)
}

/// Proves `inner` strictly inside `outer`. Parallel walls are nearest along their
/// common normal, so the axis distance `d` must fall short of the radii's difference.
/// With `outer` the larger, squaring `a^2+b^2-d^2 > 2ab` twice keeps the proof exact
/// without a root. An exact contact refuses, and a wall of any thickness stands.
pub fn certify_nested(
    outer: &Cylinder,
    inner: &Cylinder,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if sign(&outer.axis.cross(&inner.axis).norm_squared(), dec)? != Sign::Zero {
        return Err(unsupported("preserved walls need parallel axes"));
    }
    let (a, b) = (&outer.radius_sq, &inner.radius_sq);
    let spare = a.clone() + b.clone() - axis_offset(outer, inner);
    let reach = Real::from_i64(4) * (a * b);
    if sign(&(a - b), dec)? != Sign::Positive
        || sign(&spare, dec)? != Sign::Positive
        || sign(&(spare.square() - reach), dec)? != Sign::Positive
    {
        return Err(unsupported(
            "a nested wall touches or crosses its outer wall",
        ));
    }
    Ok(())
}

/// Proves every bore clear of every other halfspace, between the same two caps
/// and strictly apart from each other. Opposite axis charts of one direction
/// agree, because their caps are compared as a set. Returns those caps ascending.
pub fn certify_bores<'a>(
    planes: impl Iterator<Item = &'a Plane> + Clone,
    cylinders: &[&Cylinder],
    dec: &mut Decider,
) -> Result<[usize; 2], Refusal> {
    let mut shared = None;
    for (index, cylinder) in cylinders.iter().enumerate() {
        let mut caps = certify_clear(planes.clone(), cylinder, dec)?;
        caps.sort_unstable();
        if shared.replace(caps).is_some_and(|other| other != caps) {
            return Err(unsupported("preserved bores need the same two caps"));
        }
        for other in &cylinders[..index] {
            certify_disjoint(other, cylinder, dec)?;
        }
    }
    shared.ok_or_else(|| unsupported("a preserved bore needs a wall"))
}

impl BoredCell {
    /// Certifies the bores against the replaced envelope before any publication.
    pub fn new(
        outer: PlanarCell,
        bores: Vec<Bore>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        for bore in &bores {
            check_axial(&bore.cylinder, definition, dec)?;
        }
        let cylinders: Vec<_> = bores.iter().map(|bore| &bore.cylinder).collect();
        let caps = certify_bores(outer.supports.iter().map(|s| &s.plane), &cylinders, dec)?;
        Ok(Self {
            outer,
            bores,
            caps,
            head: None,
        })
    }

    /// The same certification for one bore with its cone.
    pub fn countersunk(
        outer: PlanarCell,
        bore: Bore,
        head: Head,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut cell = Self::new(outer, vec![bore], definition, dec)?;
        let planes = cell.outer.supports.iter().map(|s| &s.plane);
        let cylinder = &cell.bores[0].cylinder;
        let mouth = certify_head(planes, cylinder, &head.cone, cell.caps, dec)?;
        cell.head = Some((head, mouth));
        Ok(cell)
    }

    /// Exact strict material membership, excluding every bore and boundary.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.outer.contains_interior(point, dec)? && self.clear_of_bores(point, dec)?)
    }

    /// Whether a point lies strictly outside every bore's complete cylinder and cone.
    pub fn clear_of_bores(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        for Bore { cylinder, .. } in &self.bores {
            let radial = point.minus(&cylinder.origin).cross(&cylinder.axis);
            let excess =
                radial.norm_squared() - &cylinder.radius_sq * &cylinder.axis.norm_squared();
            if sign(&excess, dec)? != Sign::Positive {
                return Ok(false);
            }
        }
        match &self.head {
            Some((head, _)) => Ok(!in_nappe(&head.cone, point, dec)?),
            None => Ok(true),
        }
    }
}
