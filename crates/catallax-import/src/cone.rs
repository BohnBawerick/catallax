//! Cones at a semi-angle other than 45 degrees.
//!
//! The slope of such a cone, tan(semi-angle), is irrational for every angle a
//! file can write, so no decimal rim circle away from the placement plane lies
//! exactly on the declared cone. The importer builds the cone through the face's
//! two rims instead, or through its rim and apex. That cone has a rational slope.
//! Its whole-face deviation from the declared cone is certified, accepted only
//! within the automatic cap and recorded as the face's local tolerance.

use catallax_numbers::{Ball, Rational, Sign};

use crate::build::{chains, seam_points};
use crate::exact::{cross, dot, is_zero, sub, unit};
use crate::gap::{auto_cap, decimal};
use crate::{Curve, Face, FaceTolerance, ImportRefusal, SemiAngle, Solid, Surface, Vec3};

/// Working precision of the deviation enclosure, in bits.
const BITS: u32 = 256;
/// Decimal places of the reported bounds.
const SCALE: i64 = 1_000_000_000_000;

/// The cone built for a face, with its certified deviation from the declared cone.
#[derive(Clone, Debug)]
pub(crate) struct RimCone {
    /// Apex.
    pub apex: Vec3,
    /// Unit axis pointing into the nappe.
    pub axis: Vec3,
    /// Radius gained per unit axial distance, positive.
    pub slope: Rational,
    /// Lower and upper decimal bounds on the whole-face deviation, in mm.
    pub deviation: (Rational, Rational),
    /// Conservative nearest-feature size squared, in mm^2.
    pub feature2: Rational,
}

fn unsupported(face: &Face, why: &str) -> ImportRefusal {
    ImportRefusal::Unsupported {
        entity: face.name.clone(),
        type_name: format!("CONICAL_SURFACE ({why})"),
    }
}

fn not_solid(detail: String) -> ImportRefusal {
    ImportRefusal::NotSolid { detail }
}

/// Encloses tan(semi-angle), refusing an angle outside (0, 90) degrees.
fn tangent(face: &Face, angle: &SemiAngle) -> Result<(Ball, Ball), ImportRefusal> {
    let value = Ball::from_rational(&angle.value, BITS);
    let radians = if angle.degrees {
        value * Ball::pi(BITS) * Ball::from_rational(&Rational::from_ratio(1, 180).unwrap(), BITS)
    } else {
        value
    };
    let quarter = Ball::pi(BITS) * Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), BITS);
    let acute = radians.lo().sign() == Sign::Positive && radians.hi() < quarter.lo();
    if !acute {
        return Err(unsupported(
            face,
            "semi-angle is not strictly between 0 and 90 degrees",
        ));
    }
    let cos = radians.cos();
    let tan = radians
        .sin()
        .checked_div(&cos)
        .ok_or_else(|| unsupported(face, "semi-angle tangent could not be enclosed"))?;
    Ok((tan, cos))
}

/// Rounds an enclosure outwards to 1e-12 mm.
fn outward(lo: &Rational, hi: &Rational) -> (Rational, Rational) {
    let s = Rational::from_i64(SCALE);
    let down = &Rational::from_integer((lo * &s).floor()) / &s;
    let up = &Rational::from_integer((hi * &s).ceil()) / &s;
    (down.max(Rational::zero()), up)
}

/// The cone through the face's rims, or `None` when the face is not on such a cone.
pub(crate) fn rim_cone(face: &Face) -> Result<Option<RimCone>, ImportRefusal> {
    let Surface::Cone(p, radius, angle) = &face.surface else {
        return Ok(None);
    };
    let axis = unit(&p.axis).ok_or_else(|| unsupported(face, "axis length is not rational"))?;
    let height = |point: &Vec3| dot(&sub(point, &p.origin), &axis);
    let on_axis = |point: &Vec3| is_zero(&cross(&sub(point, &p.origin), &axis));
    // Each end is (axial height, radius); an apex has radius zero.
    let mut ends: Vec<(Rational, Rational)> = Vec::new();
    for bound in &face.loops {
        if let Some(v) = &bound.vertex {
            if !on_axis(&v.point) {
                return Err(unsupported(face, "vertex loop is not on the axis"));
            }
            ends.push((height(&v.point), Rational::zero()));
            continue;
        }
        for chain in chains(&face.name, bound.orientation, &bound.edges)? {
            let [(edge, _)] = chain.as_slice() else {
                return Err(unsupported(face, "bound is not one complete rim circle"));
            };
            let Curve::Circle(c, r) = &edge.curve else {
                return Err(unsupported(face, "bound is not one complete rim circle"));
            };
            if edge.start.point != edge.end.point
                || !is_zero(&cross(&c.axis, &axis))
                || !on_axis(&c.origin)
            {
                return Err(unsupported(
                    face,
                    "rim is not a complete circle about the axis",
                ));
            }
            ends.push((height(&c.origin), r.clone()));
        }
    }
    if ends.len() == 1 && !face.loops.iter().any(|b| b.vertex.is_some()) {
        // One rim and a seam that ends on the axis: a pointed cone.
        let apexes: Vec<_> = seam_points(face)
            .into_iter()
            .filter(|p| on_axis(p) && height(p) != ends[0].0)
            .collect();
        if let [apex] = apexes.as_slice() {
            ends.push((height(apex), Rational::zero()));
        }
    }
    let [a, b] = ends.as_slice() else {
        return Err(unsupported(
            face,
            "bounds other than two coaxial rims, or one rim and the apex",
        ));
    };
    let (mut lo, mut hi) = (a.clone(), b.clone());
    if lo.0 > hi.0 {
        std::mem::swap(&mut lo, &mut hi);
    }
    let span = &hi.0 - &lo.0;
    if span.is_zero() || lo.1 == hi.1 {
        return Err(unsupported(face, "rims do not span a cone"));
    }
    let signed = &(&hi.1 - &lo.1) / &span;
    let apex_height = &lo.0 - &(&lo.1 / &signed);
    let at = |t: &Rational| -> Vec3 { std::array::from_fn(|i| &p.origin[i] + &(&axis[i] * t)) };
    let into = if signed.sign() == Sign::Positive {
        axis.clone()
    } else {
        axis.clone().map(|c| -c)
    };

    // The declared radius at height h is radius + h tan(angle). The difference
    // from the built cone is linear in h, so its largest size over the face is
    // at an end. Matching points at equal height and angle give the upper bound.
    // An end's distance to the whole declared double cone is the radial miss
    // times cos(angle), a lower bound on the face's deviation.
    let (tan, cos) = tangent(face, angle)?;
    let mut upper = Rational::zero();
    let mut lower = Rational::zero();
    for (h, r) in [&lo, &hi] {
        let miss =
            Ball::from_rational(&(radius - r), BITS) + tan.clone() * Ball::from_rational(h, BITS);
        let size_hi = miss.lo().abs().max(miss.hi().abs());
        upper = upper.max(size_hi);
        if !miss.contains_zero() {
            let size_lo = miss.lo().abs().min(miss.hi().abs());
            lower = lower.max((Ball::from_rational(&size_lo, BITS) * cos.clone()).lo());
        }
    }
    let deviation = outward(&lower, &upper);
    let cap = auto_cap();
    if deviation.1 > cap {
        return Err(not_solid(format!(
            "conical face {}: its rims lie {} to {} mm off the declared cone, beyond the {} mm automatic cap",
            face.name,
            decimal(&deviation.0),
            decimal(&deviation.1),
            decimal(&cap)
        )));
    }
    let mut feature2 = &span * &span;
    for (_, r) in [&lo, &hi] {
        if r.sign() == Sign::Positive {
            feature2 = feature2.min(r * r);
        }
    }
    if &(&deviation.1 * &deviation.1) * &Rational::from_i64(100) > feature2 {
        return Err(not_solid(format!(
            "conical face {}: a {} mm deviation is not below a tenth of its nearest feature",
            face.name,
            decimal(&deviation.1)
        )));
    }
    Ok(Some(RimCone {
        apex: at(&apex_height),
        axis: into,
        slope: signed.abs(),
        deviation,
        feature2,
    }))
}

/// Local tolerances of the rim-defined cones of one solid, in face order.
pub(crate) fn tolerances(solid: &Solid) -> Result<Vec<FaceTolerance>, ImportRefusal> {
    let mut out = Vec::new();
    for face in &solid.faces {
        let Some(cone) = rim_cone(face)? else {
            continue;
        };
        if cone.deviation.1.is_zero() {
            continue;
        }
        out.push(FaceTolerance {
            solid: solid.name.clone(),
            face: face.name.clone(),
            deviation_lower_mm: decimal(&cone.deviation.0),
            deviation_upper_mm: decimal(&cone.deviation.1),
            approved_tolerance_mm: decimal(&auto_cap()),
            feature_size_squared_mm2: cone.feature2.to_string(),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FaceLoop, OrientedEdge, Placement, Vertex};

    fn q(s: &str) -> Rational {
        Rational::parse(s).unwrap()
    }
    fn rim(z: &str, radius: &str) -> FaceLoop {
        let point = [q(radius), q("0"), q(z)];
        let vertex = Vertex {
            name: format!("#v{z}"),
            point,
        };
        FaceLoop {
            outer: true,
            orientation: true,
            edges: vec![OrientedEdge {
                edge: format!("#e{z}"),
                start: vertex.clone(),
                end: vertex,
                curve_same_sense: true,
                curve: Curve::Circle(
                    Placement {
                        origin: [q("0"), q("0"), q(z)],
                        axis: [q("0"), q("0"), q("1")],
                        ref_dir: [q("1"), q("0"), q("0")],
                    },
                    q(radius),
                ),
                forward: true,
            }],
            vertex: None,
        }
    }
    fn face(angle: &str, degrees: bool, loops: Vec<FaceLoop>) -> Face {
        Face {
            name: "#f".into(),
            surface: Surface::Cone(
                Placement {
                    origin: [q("0"), q("0"), q("0")],
                    axis: [q("0"), q("0"), q("1")],
                    ref_dir: [q("1"), q("0"), q("0")],
                },
                q("1"),
                SemiAngle {
                    value: q(angle),
                    degrees,
                },
            ),
            same_sense: true,
            loops,
        }
    }

    #[test]
    fn a_thirty_degree_frustum_builds_through_its_rims_with_a_certified_deviation() {
        // tan 30 = 0.57735026918962..., so radius 1 + 2 tan 30 = 2.15470053837925...
        let f = face("30", true, vec![rim("0", "1"), rim("2", "2.154700538")]);
        let cone = rim_cone(&f).unwrap().unwrap();
        assert_eq!(cone.slope, q("0.577350269"));
        assert!(cone.deviation.0 > q("3e-10") && cone.deviation.1 < q("4e-10"));
        assert!(cone.deviation.0 <= cone.deviation.1);
        // The apex lies 1/slope below the lower rim.
        assert_eq!(cone.apex[2], -(q("1") / q("0.577350269")));
        // The same cone written in radians: pi/6 to fifteen places.
        let f = face(
            "0.523598775598299",
            false,
            vec![rim("0", "1"), rim("2", "2.154700538")],
        );
        assert_eq!(rim_cone(&f).unwrap().unwrap().slope, q("0.577350269"));
    }

    #[test]
    fn rims_beyond_the_cap_are_refused_with_numbers() {
        let f = face("30", true, vec![rim("0", "1"), rim("2", "2.17")]);
        let ImportRefusal::NotSolid { detail } = rim_cone(&f).unwrap_err() else {
            panic!("expected not_solid")
        };
        assert!(
            detail.contains("0.015") && detail.contains("0.01 mm"),
            "{detail}"
        );
    }

    #[test]
    fn other_bounds_and_angles_are_unsupported() {
        for f in [
            face("30", true, vec![rim("0", "1")]),
            face("30", true, vec![rim("0", "1"), rim("2", "1")]),
            face("90", true, vec![rim("0", "1"), rim("2", "2")]),
            face("1.6", false, vec![rim("0", "1"), rim("2", "2")]),
        ] {
            assert!(matches!(
                rim_cone(&f),
                Err(ImportRefusal::Unsupported { .. })
            ));
        }
    }
}
