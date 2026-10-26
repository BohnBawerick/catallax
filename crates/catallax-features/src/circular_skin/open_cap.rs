//! Finite projection proof for a cylinder with one spherical cap and open base.

use super::{CircularBand, CircularSkin, bad, reduced, sign, zero};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Curve, Surface};

fn range(band: &CircularBand, direction: &Real, reversed: bool) -> [Real; 2] {
    let [a, b] = if reversed {
        [&band.heights[1], &band.heights[0]]
    } else {
        [&band.heights[0], &band.heights[1]]
    };
    [a * direction, b * direction]
}

/// In coordinates pointing from the open base to the pole, require
/// `0 <= b < h`, `h² = R²-a²`, `R > a > 0`. The source is
/// `b <= u`, `r <= a`, `q² = r²+u² <= R²`. Its kept faces are the
/// cylinder `b <= u <= h` and sphere cap `h <= u <= R`.
///
/// Write `dc=a-r`, `ds=R-q`. These are lower bounds for distances to
/// the respective finite faces. If `dc <= ds`, then `q-r <= R-a`;
/// `u > h` would imply `q-r > R-a` for `r <= a`, so the cylindrical
/// projection lies in its finite band. If `ds <= dc`, then
/// `r <= q-(R-a) <= aq/R` because `q <= R`. With `u >= 0` this gives
/// `Ru >= hq`, so the spherical projection lies in the finite cap.
/// `q=0` cannot minimise `ds` because `a<R`. Thus the minimum always
/// has a finite witness, also at ties, the join and the closed open-rim
/// circle. Consequently `d(p,K)=min(dc,ds)` over the WHOLE source,
/// and its cavity is exactly the intersection of the original base
/// halfspace and the two inward parallels. No edge-tube approximation
/// or untrimmed-support assertion is used outside these hypotheses.
pub(super) fn prove(skin: &CircularSkin, open: &[bool], dec: &mut Decider) -> Result<(), Refusal> {
    let unsupported = || {
        bad("open curved circular cells need their kept-face edge-distance envelope")
    };
    if skin.faces.len() != 3
        || skin.bands.len() != 2
        || skin.boundaries.len() != 2
        || open.iter().filter(|v| **v).count() != 1
    {
        return Err(unsupported());
    }
    let base = open.iter().position(|v| *v).unwrap();
    let Surface::Plane(plane) = &skin.faces[base].surface else {
        return Err(unsupported());
    };
    if !zero(&plane.normal.cross(&skin.axis).norm_squared(), dec)? {
        return Err(unsupported());
    }
    let reversed = match sign(&plane.normal.dot(&skin.axis), dec)? {
        Sign::Negative => false,
        Sign::Positive => true,
        Sign::Zero => return Err(unsupported()),
    };
    let direction = Real::from_i64(if reversed { -1 } else { 1 });
    let cylinder = skin.bands.iter().find_map(|band| {
        if let Surface::Cylinder(cylinder) = &skin.faces[band.face].surface {
            Some((band, cylinder))
        } else {
            None
        }
    });
    let sphere = skin.bands.iter().find_map(|band| {
        if let Surface::Sphere(sphere) = &skin.faces[band.face].surface {
            Some((band, sphere))
        } else {
            None
        }
    });
    let (Some((wall, cylinder)), Some((cap, sphere))) = (cylinder, sphere) else {
        return Err(unsupported());
    };
    let center = &sphere.center.minus(&skin.origin).dot(&skin.axis) * &direction;
    let base_height = &plane.origin.minus(&skin.origin).dot(&skin.axis) * &direction;
    let h_sq = &sphere.radius.square() - &cylinder.radius_sq;
    if sign(&cylinder.radius_sq, dec)? != Sign::Positive
        || sign(&sphere.radius, dec)? != Sign::Positive
        || sign(&h_sq, dec)? != Sign::Positive
        || sign(&(&base_height - &center), dec)? == Sign::Negative
    {
        return Err(unsupported());
    }
    let join = &center + &reduced(h_sq.sqrt(), dec)?;
    let pole = &center + &sphere.radius;
    if sign(&(&join - &base_height), dec)? != Sign::Positive {
        return Err(unsupported());
    }
    let wall_range = range(wall, &direction, reversed);
    let cap_range = range(cap, &direction, reversed);
    for (actual, expected) in [
        (&wall_range[0], &base_height),
        (&wall_range[1], &join),
        (&cap_range[0], &join),
        (&cap_range[1], &pole),
    ] {
        if !zero(&(actual - expected), dec)? {
            return Err(unsupported());
        }
    }
    // The B-rep caller separately validates complete source loops and material
    // senses; the field extracts its own source supports. Check this written
    // definition's physical section inventory before applying the lemma.
    for (faces, height) in [
        ([base, wall.face], base_height),
        ([wall.face, cap.face], join),
    ] {
        let mut found = false;
        for boundary in &skin.boundaries {
            if !faces.iter().all(|f| boundary.faces.contains(f)) {
                continue;
            }
            let Curve::Circle(circle) = &boundary.curve else {
                return Err(unsupported());
            };
            let expected_center = skin.origin.plus(&skin.axis.scaled(&(&height * &direction)));
            if !zero(&circle.center.minus(&expected_center).norm_squared(), dec)?
                || !zero(&(&circle.radius_sq - &cylinder.radius_sq), dec)?
                || !zero(&circle.normal.cross(&skin.axis).norm_squared(), dec)?
            {
                return Err(unsupported());
            }
            found = true;
        }
        if !found {
            return Err(unsupported());
        }
    }
    Ok(())
}
