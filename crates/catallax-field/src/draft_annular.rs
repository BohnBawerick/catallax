//! Annular recognition over independent field carriers and source membership.

use super::*;
use catallax_features::draft::annular::{AnnularCell, AnnularWall};

impl Plan {
    pub(super) fn annular_cell(
        &self,
        leaf: usize,
        groups: &[Vec<usize>],
        selected: &BTreeSet<usize>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Option<AnnularCell>, Refusal> {
        if groups
            .iter()
            .filter(|g| matches!(self.faces[g[0]].surface, Surface::Cylinder(_)))
            .count()
            != 2
        {
            return Ok(None);
        }
        let unsupported = || {
            Refusal::new(
                "draft.unsupported-arrangement",
                "field annular construction needs one complete coaxial annulus and perpendicular caps",
            )
        };
        let geometry::Leaf::Feature(feature) = &self.leaves[leaf] else {
            return Err(unsupported());
        };
        if groups.len() != 4 || self.leaves[leaf].constraints().len() != 4 {
            return Err(unsupported());
        }
        let mut walls = Vec::new();
        let mut caps = Vec::new();
        for group in groups {
            let face = &self.faces[group[0]];
            let origins = group
                .iter()
                .map(|f| self.faces[*f].origin.clone())
                .collect();
            match &face.surface {
                Surface::Cylinder(c) => walls.push(AnnularWall {
                    cylinder: c.clone(),
                    origins,
                    same_sense: false,
                    selected: selected.contains(&group[0]),
                }),
                Surface::Plane(p) => {
                    if geometry::sign(&p.normal.cross(&definition.pull).norm_squared(), dec)?
                        != Sign::Zero
                    {
                        return Err(unsupported());
                    }
                    if selected.contains(&group[0]) {
                        definition.support(&face.surface, true, dec)?;
                        return Err(unsupported());
                    }
                    caps.push((p.clone(), origins));
                }
                _ => return Err(unsupported()),
            }
        }
        if walls.len() != 2 || caps.len() != 2 {
            return Err(unsupported());
        }
        match geometry::sign(
            &(&walls[0].cylinder.radius_sq - &walls[1].cylinder.radius_sq),
            dec,
        )? {
            Sign::Negative => walls.swap(0, 1),
            Sign::Zero => return Err(unsupported()),
            Sign::Positive => (),
        }
        walls[0].same_sense = true;
        let span = caps[1]
            .0
            .origin
            .minus(&caps[0].0.origin)
            .dot(&definition.pull);
        match geometry::sign(&span, dec)? {
            Sign::Negative => caps.swap(0, 1),
            Sign::Zero => return Err(unsupported()),
            Sign::Positive => (),
        }
        for (i, (plane, _)) in caps.iter_mut().enumerate() {
            let expected = if i == 0 {
                Sign::Negative
            } else {
                Sign::Positive
            };
            if geometry::sign(&plane.normal.dot(&definition.pull), dec)? != expected {
                plane.normal = plane.normal.scaled(&Real::from_i64(-1));
            }
        }
        let walls: [AnnularWall; 2] = walls.try_into().map_err(|_| unsupported())?;
        let caps: [(Plane, Vec<Origin>); 2] = caps.try_into().map_err(|_| unsupported())?;
        let mut fixed = walls.clone();
        for wall in &mut fixed {
            wall.selected = false;
        }
        let original = AnnularCell::new(fixed, caps.clone(), definition, dec)?;
        let middle = (&original.outer.heights[0] + &original.outer.heights[1]) / Real::from_i64(2);
        let centre = original
            .outer
            .origin
            .plus(&original.outer.axis.scaled(&middle));
        let radius = (&original.outer.radii[0] + &original.inner.radii[0]) / Real::from_i64(2);
        let radial = original
            .outer
            .reference
            .scaled(&(&radius / &original.outer.reference.norm_squared().sqrt()));
        // A validated finite feature with only these four complete coaxial
        // carriers has no angular or additional radial boundary. Its connected
        // annular interior and axial hole have constant source membership.
        // Classify their strict representatives to retain the material sense;
        // carrier radii alone must never turn a cavity into a filled cylinder.
        if feature.classify(&centre.plus(&radial), dec)? != Location::Inside
            || feature.classify(&centre, dec)? != Location::Outside
        {
            return Err(unsupported());
        }
        AnnularCell::new(walls, caps, definition, dec).map(Some)
    }
}

pub(super) fn transform_wall(wall: &mut RoundCell, map: &Isometry) {
    wall.origin = map.apply_point(&wall.origin);
    wall.axis = map.apply_vector(&wall.axis);
    wall.reference = map.apply_vector(&wall.reference);
    wall.surface = map.surface(&wall.surface);
    for (plane, _) in wall.caps.iter_mut().flatten() {
        let Surface::Plane(transformed) = map.surface(&Surface::Plane(plane.clone())) else {
            unreachable!("isometries preserve planes")
        };
        *plane = transformed;
    }
}

pub(super) fn distance_squared(
    cell: &AnnularCell,
    p: &Point3,
    dec: &mut Decider,
) -> Result<Real, Refusal> {
    let delta = p.minus(&cell.outer.origin);
    let h = delta.dot(&cell.outer.axis);
    let radius = delta.cross(&cell.outer.axis).norm_squared().sqrt();
    let mut best = None;
    for wall in [&cell.outer, &cell.inner] {
        let dr = &wall.radii[1] - &wall.radii[0];
        let dh = &wall.heights[1] - &wall.heights[0];
        let t = clamp(
            (&(&radius - &wall.radii[0]) * &dr + &(&h - &wall.heights[0]) * &dh)
                / (dr.square() + dh.square()),
            dec,
        )?;
        let rr = &radius - &(&wall.radii[0] + &(&t * &dr));
        let hh = &h - &(&wall.heights[0] + &(&t * &dh));
        minimum(&mut best, rr.square() + hh.square(), dec)?;
    }
    for i in 0..2 {
        let radial = if geometry::sign(&(&radius - &cell.inner.radii[i]), dec)? == Sign::Negative {
            &cell.inner.radii[i] - &radius
        } else if geometry::sign(&(&radius - &cell.outer.radii[i]), dec)? == Sign::Positive {
            &radius - &cell.outer.radii[i]
        } else {
            Real::zero()
        };
        minimum(
            &mut best,
            radial.square() + (&h - &cell.outer.heights[i]).square(),
            dec,
        )?;
    }
    Ok(best.expect("four finite annular boundaries"))
}

pub(super) fn fast_distance(
    p: [f64; 3],
    origin: [f64; 3],
    axis: [f64; 3],
    heights: [f64; 2],
    radii: [[f64; 2]; 2],
) -> (f64, bool) {
    let delta = sub(p, origin);
    let h = dot(delta, axis);
    let radial = sub(delta, axis.map(|a| a * h));
    let radius = dot(radial, radial).sqrt();
    let dh = heights[1] - heights[0];
    let mut best = f64::INFINITY;
    for wall in radii {
        let dr = wall[1] - wall[0];
        let t = (((radius - wall[0]) * dr + (h - heights[0]) * dh) / (dr * dr + dh * dh))
            .clamp(0.0, 1.0);
        best = best.min((radius - wall[0] - t * dr).powi(2) + (h - heights[0] - t * dh).powi(2));
    }
    for i in 0..2 {
        let radial = (radii[1][i] - radius).max(0.0).max(radius - radii[0][i]);
        best = best.min(radial * radial + (h - heights[i]).powi(2));
    }
    let t = (h - heights[0]) / dh;
    let outer = radii[0][0] + t * (radii[0][1] - radii[0][0]);
    let inner = radii[1][0] + t * (radii[1][1] - radii[1][0]);
    (
        best,
        h >= heights[0] && h <= heights[1] && radius >= inner && radius <= outer,
    )
}
