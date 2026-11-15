//! Cones at a semi-angle other than 45 degrees.
//!
//! The slope of such a cone, tan(semi-angle), is irrational for every angle a
//! file can write, so no decimal rim circle away from the placement plane lies
//! exactly on the declared cone. The importer builds the cone through the face's
//! two rims instead, or through its rim and apex. Its slope per axis multiple is rational.
//! Its whole-face deviation from the declared cone is certified, accepted only
//! within the automatic cap and recorded as the face's local tolerance.

use std::collections::BTreeMap;

use catallax_numbers::{Ball, Rational, Real, Sign};

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
    /// Nonzero axis pointing into the nappe.
    pub axis: Vec3,
    /// Radius gained per axis multiple, positive.
    pub slope: Rational,
    /// Squared axis length.
    pub axis_norm2: Rational,
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
pub(crate) fn tangent(face: &Face, angle: &SemiAngle) -> Result<(Ball, Ball), ImportRefusal> {
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
    let axis = unit(&p.axis).unwrap_or_else(|| p.axis.clone());
    let axis_norm2 = dot(&axis, &axis);
    if axis_norm2.is_zero() {
        return Err(unsupported(face, "axis is zero"));
    }
    let (tan, cos) = tangent(face, angle)?;
    let norm = Ball::from_rational(&axis_norm2, BITS)
        .sqrt()
        .expect("positive axis norm");
    let height = |point: &Vec3| &dot(&sub(point, &p.origin), &axis) / &axis_norm2;
    let on_axis = |point: &Vec3| is_zero(&cross(&sub(point, &p.origin), &axis));
    // Each end is (axial height, radius); an apex has radius zero.
    let mut ends: Vec<(Rational, Rational)> = Vec::new();
    let mut heights = Vec::new();
    let mut rim_witnesses = true;
    for bound in &face.loops {
        if let Some(v) = &bound.vertex {
            rim_witnesses &= on_axis(&v.point);
            ends.push((height(&v.point), Rational::zero()));
            continue;
        }
        for chain in chains(&face.name, bound.orientation, &bound.edges)? {
            for (edge, _) in chain {
                heights.extend([height(&edge.start.point), height(&edge.end.point)]);
                match &edge.curve {
                    Curve::Circle(c, r) => {
                        let h = height(&c.origin);
                        let parallel = is_zero(&cross(&c.axis, &axis));
                        rim_witnesses &= parallel && on_axis(&c.origin);
                        if !parallel {
                            let tilt = cross(&c.axis, &axis);
                            let reach2 = &(r * r) * &dot(&tilt, &tilt)
                                / (&(&axis_norm2 * &axis_norm2) * &dot(&c.axis, &c.axis));
                            let reach = Ball::from_rational(&reach2, BITS)
                                .sqrt()
                                .expect("circle projection")
                                .hi();
                            heights.extend([&h - &reach, &h + &reach]);
                        }
                        heights.push(h.clone());
                        ends.push((h, r.clone()));
                    }
                    Curve::Line(..) => {}
                    Curve::BSpline(c) => {
                        // Positive weights keep the entire curve inside this hull.
                        heights.extend(c.geometry.controls().iter().map(height));
                    }
                    Curve::Ellipse(c, major, minor) => {
                        let y = cross(&c.axis, &c.ref_dir);
                        let xh = dot(&c.ref_dir, &axis);
                        let yh = dot(&y, &axis);
                        let x2 = &(&(major * major) * &(&xh * &xh)) / &dot(&c.ref_dir, &c.ref_dir);
                        let y2 = &(&(minor * minor) * &(&yh * &yh)) / &dot(&y, &y);
                        let reach2 = &(&x2 + &y2) / &(&axis_norm2 * &axis_norm2);
                        let reach = Ball::from_rational(&reach2, BITS)
                            .sqrt()
                            .expect("ellipse projection")
                            .hi();
                        let h = height(&c.origin);
                        heights.extend([&h - &reach, &h + &reach]);
                    }
                }
            }
        }
    }
    ends.sort();
    ends.dedup();
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
    if ends.len() == 1 {
        // A pointed sector closes two distinct generators at a shared vertex,
        // without a complete rim, repeated seam or explicit vertex loop.
        let mut candidates = Vec::new();
        let mut generators = BTreeMap::<Vec3, Vec<&Vec3>>::new();
        let mut dec = catallax_numbers::Decider::default();
        for edge in face.loops.iter().flat_map(|bound| &bound.edges) {
            dec.charge_work(1)
                .map_err(|_| ImportRefusal::ResourceLimit)?;
            if let Curve::Line(_, direction) = &edge.curve {
                for vertex in [&edge.start, &edge.end] {
                    generators
                        .entry(vertex.point.clone())
                        .or_default()
                        .push(direction);
                }
            }
        }
        let cap2 = &auto_cap() * &auto_cap();
        for (point, directions) in generators {
            let h = height(&point);
            let radial = cross(&sub(&point, &p.origin), &axis);
            let distance2 = &dot(&radial, &radial) / &axis_norm2;
            if h == ends[0].0 || distance2 > cap2 {
                continue;
            }
            if directions
                .iter()
                .skip(1)
                .any(|d| !is_zero(&cross(directions[0], d)))
            {
                candidates.push(point);
            }
        }
        if let [apex] = candidates.as_slice() {
            rim_witnesses &= on_axis(apex);
            ends.push((height(apex), Rational::zero()));
        }
    }
    if ends.len() < 2 {
        let (Some(low), Some(high)) = (heights.iter().min(), heights.iter().max()) else {
            return Err(unsupported(face, "bounds have no finite axial extent"));
        };
        if low == high {
            return Err(unsupported(face, "bounds have no finite axial extent"));
        }
        // A trim need not contain two parallels. Round the declared slope to
        // 17 significant decimal digits and certify it over the full trim hull.
        let decimal_slope = (tan.clone() * norm.clone()).mid().to_decimal_bounds(17).0;
        let slope = Rational::parse(&decimal_slope).expect("decimal slope");
        let (base_height, base_radius) = ends
            .first()
            .cloned()
            .unwrap_or_else(|| (Rational::zero(), radius.clone()));
        ends = [low, high]
            .map(|h| (h.clone(), &base_radius + &(&slope * &(h - &base_height))))
            .into();
        rim_witnesses = false;
    }
    let (mut lo, mut hi) = (ends.first().unwrap().clone(), ends.last().unwrap().clone());
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
    // from the built cone is linear in h. The trim hull bounds the complete face
    // height, so its largest size is at a hull end. Matching points at equal height and angle give the upper bound.
    // In a meridian plane the declared double cone is the two lines
    // rho = +-(radius + h tan(angle)), whatever the sign of the declared radius
    // at h. An end at (h, r) lies |r -+ (radius + h tan)| cos(angle) from them,
    // and the smaller of the two is its distance to the double cone, a lower
    // bound on the face's deviation. Where the declared radius is nonnegative
    // the nearer line is the declared nappe's own.
    heights.extend([lo.0.clone(), hi.0.clone()]);
    let domain_lo = heights.iter().min().expect("rim heights").clone();
    let domain_hi = heights.iter().max().expect("rim heights").clone();
    let fitted = |h: &Rational| &lo.1 + &(&signed * &(h - &lo.0));
    if fitted(&domain_lo).sign() == Sign::Negative || fitted(&domain_hi).sign() == Sign::Negative {
        return Err(unsupported(face, "trim domain crosses the fitted apex"));
    }
    let mut upper = Rational::zero();
    let mut lower = Rational::zero();
    for h in [&domain_lo, &domain_hi] {
        let r = fitted(h);
        let declared = Ball::from_rational(radius, BITS)
            + tan.clone() * norm.clone() * Ball::from_rational(h, BITS);
        let miss = declared.clone() - Ball::from_rational(&r, BITS);
        upper = upper.max(miss.lo().abs().max(miss.hi().abs()));
    }
    for (h, r) in [&lo, &hi] {
        let declared = Ball::from_rational(radius, BITS)
            + tan.clone() * norm.clone() * Ball::from_rational(h, BITS);
        let miss = declared.clone() - Ball::from_rational(r, BITS);
        let mirror = declared + Ball::from_rational(r, BITS);
        let least = |b: &Ball| {
            if b.contains_zero() {
                Rational::zero()
            } else {
                b.lo().abs().min(b.hi().abs())
            }
        };
        let size_lo = least(&miss).min(least(&mirror));
        lower = lower.max((Ball::from_rational(&size_lo, BITS) * cos.clone()).lo());
    }
    let deviation = outward(
        &if rim_witnesses {
            lower
        } else {
            Rational::zero()
        },
        &upper,
    );
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
    let mut feature2 = &(&span * &span) * &axis_norm2;
    for (_, r) in &ends {
        if r.sign() == Sign::Positive {
            feature2 = feature2.min(r * r);
        }
    }
    // A partial rim can bound a sliver much narrower than its radius or height.
    for edge in face.loops.iter().flat_map(|bound| &bound.edges) {
        let chord = sub(&edge.end.point, &edge.start.point);
        let length2 = dot(&chord, &chord);
        if !length2.is_zero() {
            feature2 = feature2.min(length2);
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
        axis_norm2,
        deviation,
        feature2,
    }))
}

impl RimCone {
    /// Physical slope for the native unit-height parameterization.
    pub fn physical_slope(&self) -> Real {
        Real::from(self.slope.clone()) / Real::from(self.axis_norm2.clone()).sqrt()
    }

    /// Axial coordinate in multiples of the retained axis.
    pub fn height(&self, p: &Vec3) -> Rational {
        &dot(&sub(p, &self.apex), &self.axis) / &self.axis_norm2
    }

    /// A complete straight generator, without testing only its endpoints.
    pub fn contains_line(&self, origin: &Vec3, direction: &Vec3) -> bool {
        let delta = sub(origin, &self.apex);
        let axial = dot(direction, &self.axis);
        let radial = cross(direction, &self.axis);
        is_zero(&cross(&delta, direction))
            && &dot(&radial, &radial) * &self.axis_norm2
                == &(&self.slope * &self.slope) * &(&axial * &axial)
    }
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

    fn sector() -> Face {
        let mut lower = rim("0", "1").edges.remove(0);
        let mut upper = rim("2", "2.154700538").edges.remove(0);
        let a = lower.start.clone();
        let mut b = a.clone();
        b.name = "#b".into();
        b.point = [q("0"), q("1"), q("0")];
        let d = upper.start.clone();
        let mut c = d.clone();
        c.name = "#c".into();
        c.point = [q("0"), q("2.154700538"), q("2")];
        lower.end = b.clone();
        upper.start = c.clone();
        upper.curve_same_sense = false;
        let line = |name: &str, start: Vertex, end: Vertex| OrientedEdge {
            edge: name.into(),
            curve: Curve::Line(start.point.clone(), sub(&end.point, &start.point)),
            start,
            end,
            curve_same_sense: true,
            forward: true,
        };
        face(
            "30",
            true,
            vec![FaceLoop {
                outer: true,
                orientation: true,
                vertex: None,
                edges: vec![lower, line("#bc", b, c), upper, line("#da", d, a)],
            }],
        )
    }

    #[test]
    fn a_pointed_sector_uses_the_shared_generator_vertex_as_its_apex() {
        let mut f = sector();
        let apex = Vertex {
            name: "#apex".into(),
            point: [q("0"), q("0"), q("-1.732050808")],
        };
        let lower = f.loops[0].edges[0].clone();
        let generator = |name: &str, start: Vertex, end: Vertex| OrientedEdge {
            edge: name.into(),
            curve: Curve::Line(start.point.clone(), sub(&end.point, &start.point)),
            start,
            end,
            curve_same_sense: true,
            forward: true,
        };
        f.loops[0].edges = vec![
            lower.clone(),
            generator("#to-apex", lower.end.clone(), apex.clone()),
            generator("#from-apex", apex.clone(), lower.start.clone()),
        ];
        let cone = rim_cone(&f).unwrap().unwrap();
        assert_eq!(cone.apex, apex.point);
        assert!(cone.deviation.1 < q("1e-9"));
        for edge in &f.loops[0].edges {
            if let Curve::Line(origin, direction) = &edge.curve {
                assert!(cone.contains_line(origin, direction));
            }
        }
    }

    #[test]
    fn one_parallel_and_a_finite_trim_hull_determine_a_certified_cone() {
        let mut f = sector();
        let edge = &mut f.loops[0].edges[2];
        edge.curve = Curve::Line(
            edge.start.point.clone(),
            sub(&edge.end.point, &edge.start.point),
        );
        let cone = rim_cone(&f).unwrap().unwrap();
        assert_eq!(cone.deviation.0, q("0"));
        assert!(cone.deviation.1 <= q("1e-12"));
        assert!(cone.slope > q("0.577350269189") && cone.slope < q("0.577350269190"));
    }

    #[test]
    fn partial_rims_and_complete_generators_define_the_same_fitted_cone() {
        let f = sector();
        let cone = rim_cone(&f).unwrap().unwrap();
        assert_eq!(cone.slope, q("0.577350269"));
        assert_eq!(cone.deviation, (q("0.000000000328"), q("0.000000000380")));
        for edge in &f.loops[0].edges {
            if let Curve::Line(origin, direction) = &edge.curve {
                assert!(cone.contains_line(origin, direction));
                let mut shifted = origin.clone();
                shifted[0] = &shifted[0] + &q("1");
                assert!(!cone.contains_line(&shifted, direction));
            }
        }
    }

    #[test]
    fn raw_axis_keeps_an_exact_apex_and_certified_physical_slope() {
        let mut f = face(
            "30",
            true,
            vec![rim("0", "1"), rim("2", "2.632993161855452")],
        );
        let Surface::Cone(p, ..) = &mut f.surface else {
            unreachable!()
        };
        p.axis = [q("1"), q("0"), q("1")];
        p.ref_dir = [q("0"), q("1"), q("0")];
        for edge in f.loops.iter_mut().flat_map(|b| &mut b.edges) {
            for v in [&mut edge.start, &mut edge.end] {
                v.point = [v.point[2].clone(), v.point[0].clone(), v.point[2].clone()];
            }
            let Curve::Circle(p, _) = &mut edge.curve else {
                unreachable!()
            };
            p.origin = [p.origin[2].clone(), q("0"), p.origin[2].clone()];
            p.axis = [q("1"), q("0"), q("1")];
            p.ref_dir = [q("0"), q("1"), q("0")];
        }
        let cone = rim_cone(&f).unwrap().unwrap();
        assert_eq!(cone.axis_norm2, q("2"));
        assert_eq!(cone.apex[0], cone.apex[2]);
        assert!(cone.deviation.1 <= q("1e-12"));
        let physical = catallax_numbers::Decider::default()
            .enclose_at(&cone.physical_slope(), 128)
            .unwrap();
        assert!(physical.lo() > q("0.577350269189") && physical.hi() < q("0.577350269190"));
        assert_eq!(
            cone.height(&[q("2"), q("0"), q("2")]) * cone.slope,
            q("2.632993161855452")
        );
    }

    #[test]
    fn spline_trim_hull_cannot_hide_a_large_cone_deviation() {
        let mut f = sector();
        // This fitted slope differs from the declared slope enough at the remote control.
        let Surface::Cone(_, _, angle) = &mut f.surface else {
            unreachable!()
        };
        angle.value = q("30.001");
        let edge = &mut f.loops[0].edges[1];
        edge.curve = Curve::BSpline(crate::StepSplineCurve {
            geometry: catallax_shapes::spline_tools::RationalSplineCurve::new(
                2,
                vec![q("0"), q("0"), q("0"), q("1"), q("1"), q("1")],
                vec![
                    edge.start.point.clone(),
                    [q("0"), q("1"), q("1000")],
                    edge.end.point.clone(),
                ],
                vec![q("1"); 3],
                &catallax_numbers::Decider::default(),
            )
            .unwrap(),
            closed: Some(false),
            self_intersect: Some(false),
        });
        assert!(
            matches!(rim_cone(&f), Err(ImportRefusal::NotSolid { detail }) if detail.contains("automatic cap"))
        );
    }

    #[test]
    fn a_narrow_partial_rim_uses_its_chord_as_a_feature_limit() {
        let mut f = sector();
        let t = q("1e-10");
        let denominator = &q("1") + &(&t * &t);
        let x = &(&q("1") - &(&t * &t)) / &denominator;
        let y = &(&q("2") * &t) / &denominator;
        for edge in &mut f.loops[0].edges {
            for v in [&mut edge.start, &mut edge.end] {
                if v.name == "#b" || v.name == "#c" {
                    let radius = v.point[1].clone();
                    v.point[0] = &radius * &x;
                    v.point[1] = &radius * &y;
                }
            }
        }
        assert!(
            matches!(rim_cone(&f), Err(ImportRefusal::NotSolid { detail }) if detail.contains("nearest feature"))
        );
    }

    #[test]
    fn partial_rims_preserve_the_neighbouring_feature_limit() {
        let mut f = sector();
        let Surface::Cone(_, radius, angle) = &mut f.surface else {
            unreachable!()
        };
        *radius = q("0.000001");
        angle.value = q("35");
        for edge in &mut f.loops[0].edges {
            for v in [&mut edge.start, &mut edge.end] {
                v.point = v.point.clone().map(|c| c * q("0.000001"));
            }
            match &mut edge.curve {
                Curve::Circle(p, r) => {
                    p.origin = p.origin.clone().map(|c| c * q("0.000001"));
                    *r = &*r * &q("0.000001");
                }
                Curve::Line(o, d) => {
                    *o = o.clone().map(|c| c * q("0.000001"));
                    *d = d.clone().map(|c| c * q("0.000001"));
                }
                _ => unreachable!(),
            }
        }
        assert!(
            matches!(rim_cone(&f), Err(ImportRefusal::NotSolid { detail }) if detail.contains("nearest feature"))
        );
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
    fn a_frustum_past_the_declared_apex_reports_its_distance_to_the_nearer_nappe() {
        // The declared radius 1 + h tan 30 is negative at both rims, so the rims
        // are nearest the opposite nappe: at h = -4, |1 - 1.3094...| cos 30.
        let f = face("30", true, vec![rim("-4", "1"), rim("-5", "2")]);
        let ImportRefusal::NotSolid { detail } = rim_cone(&f).unwrap_err() else {
            panic!("expected not_solid")
        };
        let bounds = detail.split("lie ").nth(1).unwrap();
        let (lower, rest) = bounds.split_once(" to ").unwrap();
        let upper = rest.split_once(' ').unwrap().0;
        let (lower, upper) = (q(lower), q(upper));
        // True distance of the h = -4 rim: (4 tan 30 - 2) cos 30 = 0.26794919...
        assert!(lower > q("0.2679") && lower < q("0.26795"), "{detail}");
        // Equal-height matching gives at most 2 + 5 tan 30 - 1 = 3.886...
        assert!(upper > q("3.886") && upper < q("3.8868"), "{detail}");
    }

    #[test]
    fn accepted_bounds_are_unchanged_where_the_declared_radius_is_positive() {
        let f = face("30", true, vec![rim("0", "1"), rim("2", "2.154700538")]);
        let cone = rim_cone(&f).unwrap().unwrap();
        // (2 tan 30 - 1.154700538) cos 30 and 2 tan 30 - 1.154700538, rounded outward.
        assert_eq!(cone.deviation, (q("0.000000000328"), q("0.000000000380")));
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
