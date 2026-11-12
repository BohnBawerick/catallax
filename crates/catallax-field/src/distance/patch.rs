//! Certified surface-area, centroid and normal-weighted integration on surviving
//! named field faces (planar arrangements with line/arc trims and axial curved
//! bands/sectors), without constructing or consulting a B-rep mesh.

use super::{
    Prepared,
    planar::{PlanarFaceSection, intern, ray_sample, sort_points},
};
use crate::{FaceId, Instruction, Location, geometry};
use catallax_intersect::{CurveForm, CurveSurface, arc_point, intersect_curve_surface};
use catallax_measure::Enclosure;
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{
    Circle, Cone, Curve, Cylinder, GeomStore, Point3, Sphere, Surface, Torus, Vector3,
    distance::PlanarEdge,
};

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}

fn ball_i64(v: i64, bits: u32) -> Ball {
    Ball::from_rational(&q(v), bits)
}

fn real_half() -> Real {
    Real::from_rational(Rational::from_ratio(1, 2).unwrap())
}

fn refusal(code: &'static str, detail: impl Into<String>) -> Refusal {
    Refusal::new(code, detail)
}

fn exact_sqrt(value: &Real) -> Real {
    if let Some(r) = value.rational_ref()
        && let Some(root) = crate::profile::sqrt(r)
    {
        return Real::from(root);
    }
    value.sqrt()
}

fn unit_vec(v: &Vector3) -> Vector3 {
    v.scaled(&(Real::one() / exact_sqrt(&v.norm_squared())))
}

fn enclose_vec(v: &Vector3, bits: u32, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        geometry::enclose(&v.x, bits, dec)?,
        geometry::enclose(&v.y, bits, dec)?,
        geometry::enclose(&v.z, bits, dec)?,
    ])
}

fn enclose_point(p: &Point3, bits: u32, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    Ok([
        geometry::enclose(&p.x, bits, dec)?,
        geometry::enclose(&p.y, bits, dec)?,
        geometry::enclose(&p.z, bits, dec)?,
    ])
}

fn ball_cross(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

fn ball_add(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] + &b[i])
}

fn ball_scale(a: &[Ball; 3], s: &Ball) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] * s)
}

fn to_enclosure(b: &Ball) -> Enclosure {
    Enclosure {
        lo: b.lo(),
        hi: b.hi(),
    }
}

/// Counter-clockwise angle in `(0, 2*pi]` from `pa - center` to `pb - center`
/// about `unit_axis`.
fn ccw_arc_angle_ball(
    center: &Point3,
    unit_axis: &Vector3,
    pa: &Point3,
    pb: &Point3,
    bits: u32,
    dec: &mut Decider,
) -> Result<Ball, Refusal> {
    let va = pa.minus(center);
    let vb = pb.minus(center);
    let cross = va.cross(&vb).dot(unit_axis);
    let dot = va.dot(&vb);
    let pi = Ball::pi_checked(bits, dec)?;
    let half_pi = &pi * &Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), bits);
    let two_pi = &pi + &pi;
    let c_sign = geometry::sign(&cross, dec)?;
    let d_sign = geometry::sign(&dot, dec)?;
    match (c_sign, d_sign) {
        (Sign::Zero, Sign::Negative) => Ok(pi),
        (Sign::Zero, Sign::Positive) => Ok(two_pi),
        (Sign::Zero, Sign::Zero) => Err(refusal(
            "solver.patch-coverage",
            "degenerate zero-radius arc in face boundary",
        )),
        (Sign::Positive, Sign::Zero) => Ok(half_pi),
        (Sign::Negative, Sign::Zero) => Ok(&pi + &half_pi),
        (c_s, d_s) => {
            let c_abs = if c_s == Sign::Negative { -cross } else { cross };
            let d_abs = if d_s == Sign::Negative { -dot } else { dot };
            let c_ball = geometry::enclose(&c_abs, bits, dec)?;
            let d_ball = geometry::enclose(&d_abs, bits, dec)?;
            let acute = if geometry::sign(&(&c_abs - &d_abs), dec)? != Sign::Positive {
                let ratio = c_ball.checked_div(&d_ball).ok_or_else(|| {
                    refusal(
                        "solver.patch-coverage",
                        "arc angle tangent division by zero",
                    )
                })?;
                ratio.atan_checked(dec)?
            } else {
                let ratio = d_ball.checked_div(&c_ball).ok_or_else(|| {
                    refusal(
                        "solver.patch-coverage",
                        "arc angle cotangent division by zero",
                    )
                })?;
                &half_pi - &ratio.atan_checked(dec)?
            };
            Ok(match (c_s, d_s) {
                (Sign::Positive, Sign::Positive) => acute,
                (Sign::Positive, Sign::Negative) => &pi - &acute,
                (Sign::Negative, Sign::Negative) => &pi + &acute,
                (Sign::Negative, Sign::Positive) => &two_pi - &acute,
                _ => unreachable!(),
            })
        }
    }
}

#[derive(Clone, Debug)]
struct PlanarArcTerm {
    circle: Circle,
    unit_circle_normal: Vector3,
    ends: Option<[Point3; 2]>,
    sigma_sec: i8,
    rc: Vector3,
}

#[derive(Clone, Debug)]
struct PlanarPieceTerm {
    origin: Point3,
    outward_normal: Vector3,
    twice_poly_area: Real,
    six_poly_moment: Vector3,
    arcs: Vec<PlanarArcTerm>,
}

#[derive(Clone, Debug)]
struct CylinderArcTerm {
    center: Point3,
    unit_axis: Vector3,
    radius: Real,
    radius_sq: Real,
    height: Real,
    mid_center: Point3,
    ends: Option<[Point3; 2]>,
    chord_cross_axis: Vector3,
    outward_sign: i8,
}

#[derive(Clone, Debug)]
struct ConeBandTerm {
    apex: Point3,
    unit_axis: Vector3,
    slope: Real,
    slant: Real,
    h_lo: Real,
    h_hi: Real,
    outward_sign: i8,
}

#[derive(Clone, Debug)]
struct SphereZoneTerm {
    center: Point3,
    unit_axis: Vector3,
    radius: Real,
    h_lo: Real,
    h_hi: Real,
    outward_sign: i8,
}

#[derive(Clone, Debug)]
struct TorusFullTerm {
    center: Point3,
    major: Real,
    minor: Real,
}

#[derive(Clone, Debug)]
enum PatchTerm {
    Planar(PlanarPieceTerm),
    Cylinder(CylinderArcTerm),
    Cone(ConeBandTerm),
    Sphere(SphereZoneTerm),
    Torus(TorusFullTerm),
}

/// Certified surface integral over all surviving patches of a named face.
#[derive(Clone, Debug)]
pub(crate) struct FaceIntegral {
    pub faces: Vec<FaceId>,
    pub pieces: usize,
    pub area: Enclosure,
    pub first_moment: [Enclosure; 3],
    pub centroid: [Enclosure; 3],
    pub normal_integral: [Enclosure; 3],
    pub normal_moment: [Enclosure; 3],
    pub unit_normal: Option<[Enclosure; 3]>,
}

type EvaluatedTerm = (Ball, [Ball; 3], [Ball; 3], [Ball; 3]);

impl PatchTerm {
    fn evaluate(&self, bits: u32, dec: &mut Decider) -> Result<EvaluatedTerm, Refusal> {
        let half = Ball::from_rational(&Rational::from_ratio(1, 2).unwrap(), bits);
        let sixth = Ball::from_rational(&Rational::from_ratio(1, 6).unwrap(), bits);
        let third = Ball::from_rational(&Rational::from_ratio(1, 3).unwrap(), bits);
        match self {
            Self::Planar(term) => {
                let mut twice_area = geometry::enclose(&term.twice_poly_area, bits, dec)?;
                let mut six_moment = enclose_vec(&term.six_poly_moment, bits, dec)?;
                for arc in &term.arcs {
                    let dtheta = match &arc.ends {
                        Some([pa, pb]) => ccw_arc_angle_ball(
                            &arc.circle.center,
                            &arc.unit_circle_normal,
                            pa,
                            pb,
                            bits,
                            dec,
                        )?,
                        None => {
                            let pi = Ball::pi_checked(bits, dec)?;
                            &pi + &pi
                        }
                    };
                    let r2 = geometry::enclose(&arc.circle.radius_sq, bits, dec)?;
                    let sign = ball_i64(i64::from(arc.sigma_sec), bits);
                    let twice_sec = &(&sign * &r2) * &dtheta;
                    twice_area = &twice_area + &twice_sec;
                    let rc_ball = enclose_vec(&arc.rc, bits, dec)?;
                    let sec_moment6 = ball_scale(&rc_ball, &(&twice_sec * &ball_i64(3, bits)));
                    six_moment = ball_add(&six_moment, &sec_moment6);
                }
                let area = &twice_area * &half;
                let p0 = enclose_point(&term.origin, bits, dec)?;
                let local_moment = ball_scale(&six_moment, &sixth);
                let first_moment = ball_add(&ball_scale(&p0, &area), &local_moment);
                let n = enclose_vec(&term.outward_normal, bits, dec)?;
                let normal_integral = ball_scale(&n, &area);
                let normal_moment = ball_cross(&first_moment, &n);
                Ok((area, first_moment, normal_integral, normal_moment))
            }
            Self::Cylinder(term) => {
                let dtheta = match &term.ends {
                    Some([pa, pb]) => {
                        ccw_arc_angle_ball(&term.center, &term.unit_axis, pa, pb, bits, dec)?
                    }
                    None => {
                        let pi = Ball::pi_checked(bits, dec)?;
                        &pi + &pi
                    }
                };
                let r = geometry::enclose(&term.radius, bits, dec)?;
                let r2 = geometry::enclose(&term.radius_sq, bits, dec)?;
                let h = geometry::enclose(&term.height, bits, dec)?;
                let area = &(&r * &h) * &dtheta;
                let mid_c = enclose_point(&term.mid_center, bits, dec)?;
                let cross = enclose_vec(&term.chord_cross_axis, bits, dec)?;
                let sign = ball_i64(i64::from(term.outward_sign), bits);
                let normal_integral = ball_scale(&cross, &(&h * &sign));
                let radial_moment = ball_scale(&cross, &(&h * &r));
                let first_moment = ball_add(&ball_scale(&mid_c, &area), &radial_moment);
                let normal_moment = ball_cross(&mid_c, &normal_integral);
                let _ = r2;
                Ok((area, first_moment, normal_integral, normal_moment))
            }
            Self::Cone(term) => {
                let pi = Ball::pi_checked(bits, dec)?;
                let two_pi = &pi + &pi;
                let k = geometry::enclose(&term.slope, bits, dec)?;
                let slant = geometry::enclose(&term.slant, bits, dec)?;
                let h0 = geometry::enclose(&term.h_lo, bits, dec)?;
                let h1 = geometry::enclose(&term.h_hi, bits, dec)?;
                let dh2 = &h1.square() - &h0.square();
                let dh3 = &h1.pow_u32(3) - &h0.pow_u32(3);
                let area = &(&(&(&k * &slant) * &dh2) * &half) * &two_pi;
                let apex = enclose_point(&term.apex, bits, dec)?;
                let axis = enclose_vec(&term.unit_axis, bits, dec)?;
                let axial_moment = &(&(&(&k * &slant) * &dh3) * &third) * &two_pi;
                let first_moment =
                    ball_add(&ball_scale(&apex, &area), &ball_scale(&axis, &axial_moment));
                let sign = ball_i64(-i64::from(term.outward_sign), bits);
                let n_axial = &(&(&k.square() * &dh2) * &half) * &(&two_pi * &sign);
                let normal_integral = ball_scale(&axis, &n_axial);
                let normal_moment = ball_cross(&apex, &normal_integral);
                Ok((area, first_moment, normal_integral, normal_moment))
            }
            Self::Sphere(term) => {
                let pi = Ball::pi_checked(bits, dec)?;
                let two_pi = &pi + &pi;
                let r = geometry::enclose(&term.radius, bits, dec)?;
                let h0 = geometry::enclose(&term.h_lo, bits, dec)?;
                let h1 = geometry::enclose(&term.h_hi, bits, dec)?;
                let dh = &h1 - &h0;
                let h_mid = &(&h0 + &h1) * &half;
                let area = &(&two_pi * &r) * &dh;
                let center = enclose_point(&term.center, bits, dec)?;
                let axis = enclose_vec(&term.unit_axis, bits, dec)?;
                let mid_pt = ball_add(&center, &ball_scale(&axis, &h_mid));
                let first_moment = ball_scale(&mid_pt, &area);
                let sign = ball_i64(i64::from(term.outward_sign), bits);
                let dh2 = &h1.square() - &h0.square();
                let n_axial = &(&pi * &dh2) * &sign;
                let normal_integral = ball_scale(&axis, &n_axial);
                let normal_moment = ball_cross(&center, &normal_integral);
                Ok((area, first_moment, normal_integral, normal_moment))
            }
            Self::Torus(term) => {
                let pi = Ball::pi_checked(bits, dec)?;
                let four = ball_i64(4, bits);
                let major = geometry::enclose(&term.major, bits, dec)?;
                let minor = geometry::enclose(&term.minor, bits, dec)?;
                let area = &(&(&four * &pi.square()) * &major) * &minor;
                let center = enclose_point(&term.center, bits, dec)?;
                let first_moment = ball_scale(&center, &area);
                let zero = [ball_i64(0, bits), ball_i64(0, bits), ball_i64(0, bits)];
                Ok((area, first_moment, zero.clone(), zero))
            }
        }
    }
}

impl Prepared<'_> {
    fn planar_face_section(
        &self,
        id: usize,
        dec: &mut Decider,
    ) -> Result<&PlanarFaceSection, Refusal> {
        if let Some(section) = self.planar.get().and_then(|sections| sections.section(id)) {
            return Ok(section);
        }
        let _initializing = self
            .planar_init
            .lock()
            .expect("planar distance initialization panicked");
        if self.planar.get().is_none() {
            let sections = self.plan.prepare_planar_boundaries(dec)?;
            assert!(
                self.planar.set(sections).is_ok(),
                "planar distance sections initialize under one lock"
            );
        }
        let sections = self.planar.get().expect("planar sections were initialized");
        sections.prepare(id, dec)?;
        Ok(sections
            .section(id)
            .expect("planar section was initialized"))
    }

    fn planar_patch_term(
        &self,
        section: &PlanarFaceSection,
        dec: &mut Decider,
    ) -> Result<Option<PlanarPieceTerm>, Refusal> {
        if !section.supported {
            return Err(refusal(
                "solver.patch-coverage",
                "planar face has unsupported conic or spline trim curves",
            ));
        }
        let Some(boundary) = &section.boundary else {
            return Ok(None);
        };
        let plane = &boundary.plane;
        let unit_n = unit_vec(&plane.normal);
        let outward_sign = if section.outward_sign == 0 {
            1
        } else {
            section.outward_sign
        };
        let outward_normal = unit_n.scaled(&Real::from_i64(i64::from(outward_sign)));
        let p0 = &plane.origin;
        let mut twice_poly_area = Real::zero();
        let mut six_poly_moment = Vector3::zero();
        let mut arcs = Vec::new();
        for (edge, &edge_sign) in boundary.edges.iter().zip(&section.edge_signs) {
            let sign = if edge_sign == 0 { 1 } else { edge_sign };
            match edge {
                PlanarEdge::Segment([a_idx, b_idx]) => {
                    let (i0, i1) = if sign >= 0 {
                        (*a_idx, *b_idx)
                    } else {
                        (*b_idx, *a_idx)
                    };
                    let r0 = boundary.vertices[i0].minus(p0);
                    let r1 = boundary.vertices[i1].minus(p0);
                    let da2 = r0.cross(&r1).dot(&unit_n);
                    twice_poly_area = &twice_poly_area + &da2;
                    six_poly_moment = six_poly_moment.plus(&r0.plus(&r1).scaled(&da2));
                }
                PlanarEdge::Arc { circle, ends } => {
                    let unit_c = unit_vec(&circle.normal);
                    let sc = if geometry::sign(&unit_c.dot(&unit_n), dec)? == Sign::Negative {
                        -1i8
                    } else {
                        1i8
                    };
                    let sigma_sec = sign * sc;
                    let rc = circle.center.minus(p0);
                    if let Some([a_idx, b_idx]) = ends {
                        let (i0, i1) = if sign >= 0 {
                            (*a_idx, *b_idx)
                        } else {
                            (*b_idx, *a_idx)
                        };
                        let r0 = boundary.vertices[i0].minus(p0);
                        let r1 = boundary.vertices[i1].minus(p0);
                        let da2_1 = r0.cross(&rc).dot(&unit_n);
                        twice_poly_area = &twice_poly_area + &da2_1;
                        six_poly_moment = six_poly_moment.plus(&r0.plus(&rc).scaled(&da2_1));
                        let da2_2 = rc.cross(&r1).dot(&unit_n);
                        twice_poly_area = &twice_poly_area + &da2_2;
                        six_poly_moment = six_poly_moment.plus(&rc.plus(&r1).scaled(&da2_2));
                        let chord_cross = r1
                            .minus(&r0)
                            .cross(&unit_n)
                            .scaled(&(&circle.radius_sq * &Real::from_i64(2)));
                        six_poly_moment = six_poly_moment.plus(&chord_cross);
                        arcs.push(PlanarArcTerm {
                            circle: circle.clone(),
                            unit_circle_normal: unit_c,
                            ends: Some([
                                boundary.vertices[*a_idx].clone(),
                                boundary.vertices[*b_idx].clone(),
                            ]),
                            sigma_sec,
                            rc,
                        });
                    } else {
                        arcs.push(PlanarArcTerm {
                            circle: circle.clone(),
                            unit_circle_normal: unit_c,
                            ends: None,
                            sigma_sec,
                            rc,
                        });
                    }
                }
            }
        }
        Ok(Some(PlanarPieceTerm {
            origin: p0.clone(),
            outward_normal,
            twice_poly_area,
            six_poly_moment,
            arcs,
        }))
    }

    fn surface_outward_sign(
        &self,
        point: &Point3,
        outward_candidate: &Vector3,
        dec: &mut Decider,
    ) -> Result<i8, Refusal> {
        let all_constraints = self
            .plan
            .leaves
            .iter()
            .enumerate()
            .filter(|(i, _)| self.plan.leaf_registers[*i].is_some())
            .flat_map(|(_, leaf)| leaf.constraints());
        let plus_pt = ray_sample(all_constraints.clone(), point, outward_candidate, dec)?;
        let minus_pt = ray_sample(
            all_constraints,
            point,
            &outward_candidate.scaled(&Real::from_i64(-1)),
            dec,
        )?;
        let eval_root = |pt: &Point3, dec: &mut Decider| -> Result<Location, Refusal> {
            let mut tape = Vec::with_capacity(self.plan.instructions.len());
            for instruction in &self.plan.instructions {
                let loc = match *instruction {
                    Instruction::Empty => Location::Outside,
                    Instruction::Leaf(i) => super::location(&self.plan.leaves[i], pt, dec)?,
                    Instruction::Boolean { kind, a, b } => crate::combine(kind, tape[a], tape[b]),
                };
                tape.push(loc);
            }
            Ok(tape[self.plan.root])
        };
        let plus = eval_root(&plus_pt, dec)?;
        let minus = eval_root(&minus_pt, dec)?;
        match (minus, plus) {
            (Location::Inside, Location::Outside) => Ok(1),
            (Location::Outside, Location::Inside) => Ok(-1),
            _ => Ok(if self.plan.leaf_polarities[0] == 2 {
                -1
            } else {
                1
            }),
        }
    }

    fn collect_axial_heights(
        &self,
        origin: &Point3,
        unit_axis: &Vector3,
        min_h: &Real,
        max_h: &Real,
        dec: &mut Decider,
    ) -> Result<Vec<Real>, Refusal> {
        let mut heights = vec![min_h.clone(), max_h.clone()];
        for (i, leaf) in self.plan.leaves.iter().enumerate() {
            if self.plan.leaf_registers[i].is_none() {
                continue;
            }
            for c in leaf.constraints() {
                if let Surface::Plane(plane) = c.surface()
                    && geometry::sign(&plane.normal.cross(unit_axis).norm_squared(), dec)?
                        == Sign::Zero
                {
                    let h = plane.origin.minus(origin).dot(unit_axis);
                    if geometry::sign(&(&h - min_h), dec)? != Sign::Negative
                        && geometry::sign(&(max_h - &h), dec)? != Sign::Negative
                    {
                        heights.push(h);
                    }
                }
            }
        }
        // Insertion sort and deduplication by exact sign.
        let mut unique: Vec<Real> = Vec::new();
        for h in heights {
            let mut pos = unique.len();
            let mut dup = false;
            for (idx, existing) in unique.iter().enumerate() {
                match geometry::sign(&(&h - existing), dec)? {
                    Sign::Zero => {
                        dup = true;
                        break;
                    }
                    Sign::Negative => {
                        pos = idx;
                        break;
                    }
                    Sign::Positive => (),
                }
            }
            if !dup {
                unique.insert(pos, h);
            }
        }
        Ok(unique)
    }

    fn cylinder_patch_terms(
        &self,
        id: usize,
        cyl: &Cylinder,
        dec: &mut Decider,
    ) -> Result<(Vec<PatchTerm>, usize), Refusal> {
        let unit_axis = unit_vec(&cyl.axis);
        let unit_ref = unit_vec(&cyl.ref_dir);
        let unit_v = unit_axis.cross(&unit_ref);
        let radius = cyl.radius();
        let leaf_box = self.plan.leaves[self.plan.faces[id].leaf]
            .display_box(dec)?
            .ok_or_else(|| refusal("solver.patch-coverage", "unbounded cylinder leaf"))?;
        let origin_box = enclose_point(&cyl.origin, 128, dec)?;
        let axis_box = enclose_vec(&unit_axis, 128, dec)?;
        let h_ball = (0..3).fold(ball_i64(0, 128), |sum, k| {
            &sum + &(&(&leaf_box[k] - &origin_box[k]) * &axis_box[k])
        });
        let min_h = Real::from(h_ball.lo() - q(1));
        let max_h = Real::from(h_ball.hi() + q(1));
        let heights = self.collect_axial_heights(&cyl.origin, &unit_axis, &min_h, &max_h, dec)?;
        let mut terms = Vec::new();
        let mut pieces = 0usize;
        let mut prev_arcs: Vec<[Point3; 2]> = Vec::new();
        let mut prev_full = false;
        let r_ball = geometry::enclose(&radius, 128, dec)?;
        for pair in heights.windows(2) {
            let (h0, h1) = (&pair[0], &pair[1]);
            let span = h1 - h0;
            if geometry::sign(&span, dec)? != Sign::Positive {
                continue;
            }
            let h_mid = (h0 + h1) * real_half();
            let c_mid = cyl.origin.plus(&unit_axis.scaled(&h_mid));
            let c0 = cyl.origin.plus(&unit_axis.scaled(h0));
            let c1 = cyl.origin.plus(&unit_axis.scaled(h1));
            let b0 = enclose_point(&c0, 128, dec)?;
            let b1 = enclose_point(&c1, 128, dec)?;
            let band_box: [Ball; 3] = std::array::from_fn(|k| {
                let hull = b0[k].hull(&b1[k]);
                (&hull - &r_ball).hull(&(&hull + &r_ball))
            });
            // Verify every active constraint intersecting this band is either
            // perpendicular to the axis (already split at heights) or parallel
            // to the axis (plane or parallel cylinder).
            let mut geom = GeomStore::new();
            let circle = Circle {
                center: c_mid.clone(),
                normal: unit_axis.clone(),
                radius_sq: cyl.radius_sq.clone(),
                ref_dir: unit_ref.clone(),
            };
            let curve = geom.add_curve(Curve::Circle(circle.clone()));
            let form = CurveForm::Circle {
                center: c_mid.clone(),
                normal: unit_axis.clone(),
                radius_sq: cyl.radius_sq.clone(),
                ref_dir: unit_ref.clone(),
            };
            let mut cut_surfaces = Vec::new();
            for (leaf_idx, leaf) in self.plan.leaves.iter().enumerate() {
                if self.plan.leaf_registers[leaf_idx].is_none() {
                    continue;
                }
                if let Some(lb) = leaf.display_box(dec)?
                    && (0..3)
                        .any(|k| lb[k].hi() < band_box[k].lo() || band_box[k].hi() < lb[k].lo())
                {
                    continue;
                }
                if matches!(leaf, geometry::Leaf::Primitive(_) | geometry::Leaf::Cell(_)) {
                    let mut outside = false;
                    for c in leaf.constraints() {
                        if c.polynomial.eval_ball(&band_box, dec)?.lo() > q(0) {
                            outside = true;
                            break;
                        }
                    }
                    if outside {
                        continue;
                    }
                }
                for (c_idx, c) in leaf.constraints().iter().enumerate() {
                    if self.coincident[id][leaf_idx][c_idx] {
                        continue;
                    }
                    if !c.polynomial.eval_ball(&band_box, dec)?.contains_zero() {
                        continue;
                    }
                    match c.surface() {
                        Surface::Plane(p) => {
                            if geometry::sign(&p.normal.cross(&unit_axis).norm_squared(), dec)?
                                == Sign::Zero
                            {
                                continue;
                            }
                            if geometry::sign(&p.normal.dot(&unit_axis), dec)? != Sign::Zero {
                                return Err(refusal(
                                    "solver.patch-coverage",
                                    "cylinder face cut by an oblique plane requires general curved surface quadrature",
                                ));
                            }
                            cut_surfaces.push(geom.add_surface(c.surface().clone()));
                        }
                        Surface::Cylinder(other_cyl) => {
                            if geometry::sign(
                                &other_cyl.axis.cross(&unit_axis).norm_squared(),
                                dec,
                            )? != Sign::Zero
                            {
                                return Err(refusal(
                                    "solver.patch-coverage",
                                    "cylinder face cut by a non-parallel cylinder requires general curved surface quadrature",
                                ));
                            }
                            cut_surfaces.push(geom.add_surface(c.surface().clone()));
                        }
                        _ => {
                            return Err(refusal(
                                "solver.patch-coverage",
                                "cylinder face cut by a non-prismatic surface requires general curved surface quadrature",
                            ));
                        }
                    }
                }
            }
            let offset = unit_ref.scaled(&radius);
            let mut points = Vec::new();
            intern(&mut geom, &mut points, c_mid.plus(&offset), dec)?;
            intern(
                &mut geom,
                &mut points,
                c_mid.plus(&offset.scaled(&Real::from_i64(-1))),
                dec,
            )?;
            let mut has_real_cuts = false;
            for surface in cut_surfaces {
                let answer = intersect_curve_surface(&mut geom, curve, surface, dec)?;
                if !answer.marks.is_empty() {
                    return Err(refusal(
                        "solver.patch-coverage",
                        "cylinder section trim has marked intersections",
                    ));
                }
                if let CurveSurface::Meets(hits) = answer.value {
                    if !hits.is_empty() {
                        has_real_cuts = true;
                    }
                    for hit in hits {
                        let pt = geom.position(hit.point);
                        intern(&mut geom, &mut points, pt, dec)?;
                    }
                }
            }
            sort_points(&geom, &mut points, &form, dec)?;
            let mut surviving_arcs: Vec<([Point3; 2], Point3, i8)> = Vec::new();
            for i in 0..points.len() {
                let a = points[i];
                let b = points[(i + 1) % points.len()];
                let mid = arc_point(&geom, curve, Some(a), Some(b), dec)?.value;
                if self.member(id, &mid, dec)? {
                    let radial = mid.minus(&c_mid);
                    let outward = self.surface_outward_sign(&mid, &radial, dec)?;
                    surviving_arcs.push(([geom.position(a), geom.position(b)], mid, outward));
                }
            }
            if surviving_arcs.is_empty() {
                prev_arcs.clear();
                prev_full = false;
                continue;
            }
            if !has_real_cuts && surviving_arcs.len() == points.len() {
                let outward = surviving_arcs[0].2;
                let border_pt = c0.plus(&offset);
                if !(prev_full && self.member(id, &border_pt, dec)?) {
                    pieces += 1;
                }
                prev_full = true;
                prev_arcs.clear();
                terms.push(PatchTerm::Cylinder(CylinderArcTerm {
                    center: c_mid.clone(),
                    unit_axis: unit_axis.clone(),
                    radius: radius.clone(),
                    radius_sq: cyl.radius_sq.clone(),
                    height: span,
                    mid_center: c_mid,
                    ends: None,
                    chord_cross_axis: Vector3::zero(),
                    outward_sign: outward,
                }));
                continue;
            }
            // Coalesce consecutive sub-arcs that were only split by artificial
            // diameter points when those diameter points are interior to the arc.
            let mut merged: Vec<([Point3; 2], Point3, i8)> = Vec::new();
            for (ends, mid, outward) in surviving_arcs {
                if let Some((prev_ends, _, prev_out)) = merged.last_mut()
                    && *prev_out == outward
                    && geometry::sign(&prev_ends[1].minus(&ends[0]).norm_squared(), dec)?
                        == Sign::Zero
                    && self.member(id, &ends[0], dec)?
                {
                    prev_ends[1] = ends[1].clone();
                } else {
                    merged.push((ends, mid, outward));
                }
            }
            if merged.len() > 1
                && let (Some(first), Some(last)) = (merged.first().cloned(), merged.last().cloned())
                && first.2 == last.2
                && geometry::sign(&last.0[1].minus(&first.0[0]).norm_squared(), dec)? == Sign::Zero
                && self.member(id, &first.0[0], dec)?
            {
                let last_idx = merged.len() - 1;
                merged[0].0[0] = merged[last_idx].0[0].clone();
                merged.pop();
            }
            if merged.len() == 1
                && geometry::sign(&merged[0].0[0].minus(&merged[0].0[1]).norm_squared(), dec)?
                    == Sign::Zero
            {
                let outward = merged[0].2;
                let border_pt = c0.plus(&offset);
                if !(prev_full && self.member(id, &border_pt, dec)?) {
                    pieces += 1;
                }
                prev_full = true;
                prev_arcs.clear();
                terms.push(PatchTerm::Cylinder(CylinderArcTerm {
                    center: c_mid.clone(),
                    unit_axis: unit_axis.clone(),
                    radius: radius.clone(),
                    radius_sq: cyl.radius_sq.clone(),
                    height: span,
                    mid_center: c_mid,
                    ends: None,
                    chord_cross_axis: Vector3::zero(),
                    outward_sign: outward,
                }));
                continue;
            }
            prev_full = false;
            let mut current_arcs = Vec::with_capacity(merged.len());
            for (ends, mid, outward) in merged {
                let d_mid = mid.minus(&c_mid);
                let border_mid = c0.plus(&d_mid);
                let continues = prev_arcs.iter().any(|prev| {
                    let da = prev[0].minus(&ends[0]);
                    let db = prev[1].minus(&ends[1]);
                    let da_perp = da.minus(&unit_axis.scaled(&da.dot(&unit_axis)));
                    let db_perp = db.minus(&unit_axis.scaled(&db.dot(&unit_axis)));
                    geometry::sign(&da_perp.norm_squared(), dec).ok() == Some(Sign::Zero)
                        && geometry::sign(&db_perp.norm_squared(), dec).ok() == Some(Sign::Zero)
                }) && self.member(id, &border_mid, dec)?;
                if !continues {
                    pieces += 1;
                }
                let chord = ends[1].minus(&ends[0]);
                let chord_cross_axis = chord.cross(&unit_axis);
                current_arcs.push(ends.clone());
                terms.push(PatchTerm::Cylinder(CylinderArcTerm {
                    center: c_mid.clone(),
                    unit_axis: unit_axis.clone(),
                    radius: radius.clone(),
                    radius_sq: cyl.radius_sq.clone(),
                    height: span.clone(),
                    mid_center: c_mid.clone(),
                    ends: Some(ends),
                    chord_cross_axis,
                    outward_sign: outward,
                }));
            }
            prev_arcs = current_arcs;
            let _ = unit_v;
        }
        Ok((terms, pieces))
    }

    fn cone_patch_terms(
        &self,
        id: usize,
        cone: &Cone,
        dec: &mut Decider,
    ) -> Result<(Vec<PatchTerm>, usize), Refusal> {
        let unit_axis = unit_vec(&cone.axis);
        let unit_ref = unit_vec(&cone.ref_dir);
        let slant = exact_sqrt(&(Real::one() + cone.slope.square()));
        let leaf_box = self.plan.leaves[self.plan.faces[id].leaf]
            .display_box(dec)?
            .ok_or_else(|| refusal("solver.patch-coverage", "unbounded cone leaf"))?;
        let apex_box = enclose_point(&cone.apex, 128, dec)?;
        let axis_box = enclose_vec(&unit_axis, 128, dec)?;
        let h_ball = (0..3).fold(ball_i64(0, 128), |sum, k| {
            &sum + &(&(&leaf_box[k] - &apex_box[k]) * &axis_box[k])
        });
        let min_h = Real::from(h_ball.lo().max(q(0)));
        let max_h = Real::from((h_ball.hi() + q(1)).max(q(1)));
        let heights = self.collect_axial_heights(&cone.apex, &unit_axis, &min_h, &max_h, dec)?;
        let mut terms = Vec::new();
        let mut pieces = 0usize;
        let mut prev_active = false;
        for pair in heights.windows(2) {
            let (h0, h1) = (&pair[0], &pair[1]);
            if geometry::sign(&(h1 - h0), dec)? != Sign::Positive
                || geometry::sign(h0, dec)? == Sign::Negative
            {
                continue;
            }
            let h_mid = (h0 + h1) * real_half();
            let r_mid = &cone.slope * &h_mid;
            let c_mid = cone.apex.plus(&unit_axis.scaled(&h_mid));
            let p1 = c_mid.plus(&unit_ref.scaled(&r_mid));
            let p2 = c_mid.plus(&unit_ref.scaled(&(-&r_mid)));
            if !(self.member(id, &p1, dec)? && self.member(id, &p2, dec)?) {
                prev_active = false;
                continue;
            }
            let r_max = &cone.slope * h1;
            let r_ball = geometry::enclose(&r_max, 128, dec)?;
            let b0 = enclose_point(&cone.apex.plus(&unit_axis.scaled(h0)), 128, dec)?;
            let b1 = enclose_point(&cone.apex.plus(&unit_axis.scaled(h1)), 128, dec)?;
            let band_box: [Ball; 3] = std::array::from_fn(|k| {
                let hull = b0[k].hull(&b1[k]);
                (&hull - &r_ball).hull(&(&hull + &r_ball))
            });
            for (leaf_idx, leaf) in self.plan.leaves.iter().enumerate() {
                if self.plan.leaf_registers[leaf_idx].is_none() {
                    continue;
                }
                if let Some(lb) = leaf.display_box(dec)?
                    && (0..3)
                        .any(|k| lb[k].hi() < band_box[k].lo() || band_box[k].hi() < lb[k].lo())
                {
                    continue;
                }
                for (c_idx, c) in leaf.constraints().iter().enumerate() {
                    if self.coincident[id][leaf_idx][c_idx] {
                        continue;
                    }
                    if !c.polynomial.eval_ball(&band_box, dec)?.contains_zero() {
                        continue;
                    }
                    if let Surface::Plane(p) = c.surface()
                        && geometry::sign(&p.normal.cross(&unit_axis).norm_squared(), dec)?
                            == Sign::Zero
                    {
                        continue;
                    }
                    return Err(refusal(
                        "solver.patch-coverage",
                        "non-coaxial cone trim requires general curved surface quadrature",
                    ));
                }
            }
            let normal_cand = unit_ref.minus(&unit_axis.scaled(&cone.slope));
            let outward = self.surface_outward_sign(&p1, &normal_cand, dec)?;
            if !prev_active {
                pieces += 1;
            }
            prev_active = true;
            terms.push(PatchTerm::Cone(ConeBandTerm {
                apex: cone.apex.clone(),
                unit_axis: unit_axis.clone(),
                slope: cone.slope.clone(),
                slant: slant.clone(),
                h_lo: h0.clone(),
                h_hi: h1.clone(),
                outward_sign: outward,
            }));
        }
        Ok((terms, pieces))
    }

    fn sphere_patch_terms(
        &self,
        id: usize,
        sphere: &Sphere,
        dec: &mut Decider,
    ) -> Result<(Vec<PatchTerm>, usize), Refusal> {
        let unit_axis = unit_vec(&sphere.axis);
        let unit_ref = unit_vec(&sphere.ref_dir);
        let min_h = -&sphere.radius;
        let max_h = sphere.radius.clone();
        let heights =
            self.collect_axial_heights(&sphere.center, &unit_axis, &min_h, &max_h, dec)?;
        let mut terms = Vec::new();
        let mut pieces = 0usize;
        let mut prev_active = false;
        let r_ball = geometry::enclose(&sphere.radius, 128, dec)?;
        for pair in heights.windows(2) {
            let (h0, h1) = (&pair[0], &pair[1]);
            if geometry::sign(&(h1 - h0), dec)? != Sign::Positive {
                continue;
            }
            let h_mid = (h0 + h1) * real_half();
            let r_mid = exact_sqrt(&(sphere.radius.square() - h_mid.square()));
            let c_mid = sphere.center.plus(&unit_axis.scaled(&h_mid));
            let p1 = c_mid.plus(&unit_ref.scaled(&r_mid));
            let p2 = c_mid.plus(&unit_ref.scaled(&(-&r_mid)));
            if !(self.member(id, &p1, dec)? && self.member(id, &p2, dec)?) {
                prev_active = false;
                continue;
            }
            let b0 = enclose_point(&sphere.center.plus(&unit_axis.scaled(h0)), 128, dec)?;
            let b1 = enclose_point(&sphere.center.plus(&unit_axis.scaled(h1)), 128, dec)?;
            let band_box: [Ball; 3] = std::array::from_fn(|k| {
                let hull = b0[k].hull(&b1[k]);
                (&hull - &r_ball).hull(&(&hull + &r_ball))
            });
            for (leaf_idx, leaf) in self.plan.leaves.iter().enumerate() {
                if self.plan.leaf_registers[leaf_idx].is_none() {
                    continue;
                }
                if let Some(lb) = leaf.display_box(dec)?
                    && (0..3)
                        .any(|k| lb[k].hi() < band_box[k].lo() || band_box[k].hi() < lb[k].lo())
                {
                    continue;
                }
                for (c_idx, c) in leaf.constraints().iter().enumerate() {
                    if self.coincident[id][leaf_idx][c_idx] {
                        continue;
                    }
                    if !c.polynomial.eval_ball(&band_box, dec)?.contains_zero() {
                        continue;
                    }
                    if let Surface::Plane(p) = c.surface()
                        && geometry::sign(&p.normal.cross(&unit_axis).norm_squared(), dec)?
                            == Sign::Zero
                    {
                        continue;
                    }
                    return Err(refusal(
                        "solver.patch-coverage",
                        "non-zonal sphere trim requires general curved surface quadrature",
                    ));
                }
            }
            let radial = p1.minus(&sphere.center);
            let outward = self.surface_outward_sign(&p1, &radial, dec)?;
            if !prev_active {
                pieces += 1;
            }
            prev_active = true;
            terms.push(PatchTerm::Sphere(SphereZoneTerm {
                center: sphere.center.clone(),
                unit_axis: unit_axis.clone(),
                radius: sphere.radius.clone(),
                h_lo: h0.clone(),
                h_hi: h1.clone(),
                outward_sign: outward,
            }));
        }
        Ok((terms, pieces))
    }

    fn torus_patch_terms(
        &self,
        id: usize,
        torus: &Torus,
        dec: &mut Decider,
    ) -> Result<(Vec<PatchTerm>, usize), Refusal> {
        let unit_ref = unit_vec(&torus.ref_dir);
        let p_outer = torus
            .center
            .plus(&unit_ref.scaled(&(&torus.major + &torus.minor)));
        let p_inner = torus
            .center
            .plus(&unit_ref.scaled(&(&torus.major - &torus.minor)));
        if !(self.member(id, &p_outer, dec)? && self.member(id, &p_inner, dec)?) {
            return Err(refusal(
                "solver.patch-coverage",
                "partially trimmed torus patch requires general curved surface quadrature",
            ));
        }
        let r_outer = geometry::enclose(&(&torus.major + &torus.minor), 128, dec)?;
        let center_box = enclose_point(&torus.center, 128, dec)?;
        let torus_box: [Ball; 3] =
            std::array::from_fn(|k| (&center_box[k] - &r_outer).hull(&(&center_box[k] + &r_outer)));
        for (leaf_idx, leaf) in self.plan.leaves.iter().enumerate() {
            if self.plan.leaf_registers[leaf_idx].is_none() {
                continue;
            }
            if let Some(lb) = leaf.display_box(dec)?
                && (0..3).any(|k| lb[k].hi() < torus_box[k].lo() || torus_box[k].hi() < lb[k].lo())
            {
                continue;
            }
            for (c_idx, c) in leaf.constraints().iter().enumerate() {
                if self.coincident[id][leaf_idx][c_idx] {
                    continue;
                }
                if c.polynomial.eval_ball(&torus_box, dec)?.contains_zero() {
                    return Err(refusal(
                        "solver.patch-coverage",
                        "trimmed torus patch requires general curved surface quadrature",
                    ));
                }
            }
        }
        Ok((
            vec![PatchTerm::Torus(TorusFullTerm {
                center: torus.center.clone(),
                major: torus.major.clone(),
                minor: torus.minor.clone(),
            })],
            1,
        ))
    }

    pub(crate) fn integrate_origin(
        &self,
        origin: &Origin,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<FaceIntegral, Refusal> {
        if width <= &Rational::zero() {
            return Err(refusal(
                "solver.patch-options",
                "patch integration requires a positive enclosure width",
            ));
        }
        if !self.plan.faces.iter().any(|f| &f.origin == origin) {
            return Err(refusal(
                "solver.face-name",
                "name is absent from this revision",
            ));
        }
        let candidates: Vec<usize> = self
            .plan
            .faces
            .iter()
            .enumerate()
            .filter_map(|(i, f)| (&f.origin == origin && self.domains[i].is_some()).then_some(i))
            .collect();
        if candidates.is_empty() {
            return Err(refusal(
                "solver.patch-membership",
                "selected face has no surviving boundary area in this result",
            ));
        }
        let mut surviving_faces = Vec::new();
        let mut terms = Vec::new();
        let mut total_pieces = 0usize;
        for id in candidates {
            match &self.plan.faces[id].surface {
                Surface::Plane(_) => {
                    let section = self.planar_face_section(id, dec)?;
                    if let Some(term) = self.planar_patch_term(section, dec)? {
                        surviving_faces.push(FaceId(id));
                        total_pieces += 1;
                        terms.push(PatchTerm::Planar(term));
                    }
                }
                Surface::Cylinder(cyl) => {
                    let (cyl_terms, count) = self.cylinder_patch_terms(id, cyl, dec)?;
                    if !cyl_terms.is_empty() {
                        surviving_faces.push(FaceId(id));
                        total_pieces += count;
                        terms.extend(cyl_terms);
                    }
                }
                Surface::Cone(cone) => {
                    let (cone_terms, count) = self.cone_patch_terms(id, cone, dec)?;
                    if !cone_terms.is_empty() {
                        surviving_faces.push(FaceId(id));
                        total_pieces += count;
                        terms.extend(cone_terms);
                    }
                }
                Surface::Sphere(sphere) => {
                    let (sphere_terms, count) = self.sphere_patch_terms(id, sphere, dec)?;
                    if !sphere_terms.is_empty() {
                        surviving_faces.push(FaceId(id));
                        total_pieces += count;
                        terms.extend(sphere_terms);
                    }
                }
                Surface::Torus(torus) => {
                    let (torus_terms, count) = self.torus_patch_terms(id, torus, dec)?;
                    if !torus_terms.is_empty() {
                        surviving_faces.push(FaceId(id));
                        total_pieces += count;
                        terms.extend(torus_terms);
                    }
                }
                _ => {
                    return Err(refusal(
                        "solver.patch-coverage",
                        "named face carrier is outside analytical patch integration coverage",
                    ));
                }
            }
        }
        if terms.is_empty() {
            return Err(refusal(
                "solver.patch-membership",
                "selected face has no surviving boundary area in this result",
            ));
        }
        // Fast exact rational check when all terms are straight-edged planar polygons.
        if terms
            .iter()
            .all(|t| matches!(t, PatchTerm::Planar(p) if p.arcs.is_empty()))
        {
            let mut total_area_q = Some(q(0));
            let mut total_moment_q = Some([q(0), q(0), q(0)]);
            let mut common_normal_q: Option<[Rational; 3]> = None;
            let mut same_normal = true;
            for term in &terms {
                let PatchTerm::Planar(p) = term else {
                    unreachable!()
                };
                if let (Some(a2), Some(m6), Some(p0), Some(n)) = (
                    p.twice_poly_area.as_rational(),
                    p.six_poly_moment.as_rational(),
                    p.origin.as_rational(),
                    p.outward_normal.as_rational(),
                ) {
                    let area = a2.mul_pow2(-1);
                    let m_local = [m6.x / q(6), m6.y / q(6), m6.z / q(6)];
                    let p0_arr = [p0.x, p0.y, p0.z];
                    let n_arr = [n.x, n.y, n.z];
                    if let Some(prev) = &common_normal_q {
                        if prev != &n_arr {
                            same_normal = false;
                        }
                    } else {
                        common_normal_q = Some(n_arr);
                    }
                    if let (Some(ta), Some(tm)) = (&mut total_area_q, &mut total_moment_q) {
                        *ta = &*ta + &area;
                        for k in 0..3 {
                            tm[k] = &tm[k] + &(&(&area * &p0_arr[k]) + &m_local[k]);
                        }
                    }
                } else {
                    total_area_q = None;
                    break;
                }
            }
            if let (Some(area_q), Some(moment_q), Some(n_q)) =
                (total_area_q, total_moment_q, common_normal_q)
                && same_normal
                && area_q > q(0)
            {
                let exact = |r: Rational| Enclosure {
                    lo: r.clone(),
                    hi: r,
                };
                let centroid_q: [Rational; 3] = std::array::from_fn(|k| &moment_q[k] / &area_q);
                let normal_int_q: [Rational; 3] = std::array::from_fn(|k| &area_q * &n_q[k]);
                let normal_mom_q: [Rational; 3] = [
                    &(&moment_q[1] * &n_q[2]) - &(&moment_q[2] * &n_q[1]),
                    &(&moment_q[2] * &n_q[0]) - &(&moment_q[0] * &n_q[2]),
                    &(&moment_q[0] * &n_q[1]) - &(&moment_q[1] * &n_q[0]),
                ];
                return Ok(FaceIntegral {
                    faces: surviving_faces,
                    pieces: total_pieces,
                    area: exact(area_q),
                    first_moment: moment_q.map(exact),
                    centroid: centroid_q.map(exact),
                    normal_integral: normal_int_q.map(exact),
                    normal_moment: normal_mom_q.map(exact),
                    unit_normal: Some(n_q.map(exact)),
                });
            }
        }
        let max_bits = dec.policy().max_precision_bits.max(128);
        let mut bits = 128u32.min(max_bits);
        loop {
            dec.check_deadline()?;
            let mut area = ball_i64(0, bits);
            let mut first_moment = [ball_i64(0, bits), ball_i64(0, bits), ball_i64(0, bits)];
            let mut normal_integral = first_moment.clone();
            let mut normal_moment = first_moment.clone();
            for term in &terms {
                let (a, fm, ni, nm) = term.evaluate(bits, dec)?;
                area = &area + &a;
                first_moment = ball_add(&first_moment, &fm);
                normal_integral = ball_add(&normal_integral, &ni);
                normal_moment = ball_add(&normal_moment, &nm);
            }
            if area.lo() > q(0) {
                let centroid: [Ball; 3] =
                    std::array::from_fn(|k| first_moment[k].checked_div(&area).unwrap());
                let area_enc = to_enclosure(&area);
                let centroid_enc = centroid.each_ref().map(to_enclosure);
                let normal_enc = normal_integral.each_ref().map(to_enclosure);
                let nm_enc = normal_moment.each_ref().map(to_enclosure);
                let unit_normal = if let [PatchTerm::Planar(p)] = terms.as_slice() {
                    Some(enclose_vec(&p.outward_normal, bits, dec)?.map(|b| to_enclosure(&b)))
                } else {
                    None
                };
                if area_enc.width() <= *width
                    && centroid_enc.iter().all(|c| c.width() <= *width)
                    && normal_enc.iter().all(|n| n.width() <= *width)
                {
                    return Ok(FaceIntegral {
                        faces: surviving_faces,
                        pieces: total_pieces,
                        area: area_enc,
                        first_moment: first_moment.map(|b| to_enclosure(&b)),
                        centroid: centroid_enc,
                        normal_integral: normal_enc,
                        normal_moment: nm_enc,
                        unit_normal,
                    });
                }
            }
            if bits >= max_bits {
                return Err(refusal(
                    "solver.patch-accuracy",
                    format!(
                        "face patch integral did not reach width {width} at {bits} bits (attained area=[{}, {}])",
                        area.lo(),
                        area.hi()
                    ),
                ));
            }
            bits = bits.saturating_mul(2).min(max_bits);
        }
    }
}
