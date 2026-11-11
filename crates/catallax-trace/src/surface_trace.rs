// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surfaces and planes of
// `spline.rs`; statuses become `Fault`, the intersection point is a plain
// struct and the intersector's tracing state the `Intersector`. Paths that
// need code from later parts of step 8 refuse with
// `SolverRefusal::NotPorted`: the step back from a surface singularity
// (`StepBackFromSingularity`), a short range on a plane's side line, and a
// side range that needs SMLib's general curve drop. Long side ranges become
// curves through `DropAndTrimCurve`.
// Translated from SMLib SmSurfaceIntersector.cpp (SmTsectPnt,
// SmTsectPnt::ClassifyPointParams, SmSurfaceIntersector::ComputePointValues),
// SmAdvSurfaceIntersector.cpp (ComputeSingularityPoint, MakeDupin,
// sm_AddSingDir) and
// SmSurface.cpp (smsurf_EvaluateGeometric, smsurf_ComputeSurfacePointType);
// rust-brep has no port.
// Indexed loops follow SMLib's per-surface loops.
#![allow(clippy::needless_range_loop)]
//! The intersection points of step 8: a pair of surface parameters with the
//! surfaces' derivatives and normals, the meeting curve's point and unit
//! tangent, and the tangent projected into each parameter domain.

use crate::Fault;
use crate::conic::{Conic, ConicKind};
use crate::curve_surface::{drop_vector, is_singularity};
use crate::extent::{Extent1, Extent2, PointType};
use crate::global_solver::{SolverOperation, SolverRefusal};
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI, ZONE_TOL_3D, max, min};
use crate::spline::BSplineCurve;
use crate::surface_curves::IntersectionCurve;
use crate::surface_eval::{evaluate, evaluate_normal};
use crate::surface_intersector::Intersector;
use crate::surface_surface::local_surface_solve;
use crate::vector::{Vector2, Vector3};

/// `SmIntersectionPointType`, in SMLib's order (the tracer compares them).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum PointKind {
    /// `SM_IP_UNDEFINED`.
    Undefined,
    /// `SM_IP_UNKNOWN`.
    Unknown,
    /// `SM_IP_CROSSING`: the normals are not parallel.
    Crossing,
    /// `SM_IP_TANGENT_POINT`: the surfaces touch at a point.
    TangentPoint,
    /// `SM_IP_TANGENT_CURVE`: the point lies on a tangent curve.
    TangentCurve,
    /// `SM_IP_SINGULARITY`: four or six branches meet.
    Singularity,
    /// `SM_IP_COINCIDENCE`: the second fundamental forms match.
    Coincidence,
    /// `SM_IP_TOUCHING`.
    Touching,
}

/// `SmTsectPnt`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TsectPoint {
    /// `m_ePointType`.
    pub kind: PointKind,
    /// `m_eUVParamType`: each surface's u and v against its domain.
    pub param_types: [[PointType; 2]; 2],
    /// `UVPos`.
    pub uv: [Vector2; 2],
    /// `UVDeriv`: the curve tangent in each parameter domain.
    pub uv_deriv: [Vector2; 2],
    /// `m_vSurfacePV`: SMLib's evaluation matrices, cell `k * 3 + l` for
    /// the derivative `S^(k,l)`, with the unit normal in cell 8.
    pub surface: [[Vector3; 9]; 2],
    /// `CrvPos`.
    pub curve_point: Vector3,
    /// `CrvDeriv`.
    pub curve_tangent: Vector3,
    /// `m_dTangentPlaneAngleRad`.
    pub tangent_plane_angle: f64,
    /// `m_dDeviation`.
    pub deviation: f64,
    /// `m_dCurveParameter`: arc length along the traced curve.
    pub curve_parameter: f64,
}

impl Default for TsectPoint {
    /// `SmTsectPnt()`: unknown and unclassified, its vectors undefined
    /// (`SM_UNDEF_DOUBLE`) and its numbers zero.
    fn default() -> Self {
        let undefined = Vector3::default();
        Self {
            kind: PointKind::Unknown,
            param_types: [[PointType::Outside; 2]; 2],
            uv: [Vector2::default(); 2],
            uv_deriv: [Vector2::default(); 2],
            surface: [[undefined; 9]; 2],
            curve_point: undefined,
            curve_tangent: undefined,
            tangent_plane_angle: 0.0,
            deviation: 0.0,
            curve_parameter: 0.0,
        }
    }
}

impl TsectPoint {
    /// `SmTsectPnt::IsBounded`: a parameter lies on a side of its domain or
    /// outside it. (SMLib's unknown class, which the port's unclassified
    /// default reads as outside, never reaches a caller.)
    pub fn is_bounded(&self) -> bool {
        self.param_types
            .iter()
            .flatten()
            .any(|t| *t != PointType::Inside)
    }

    /// `SrfPos`.
    pub fn position(&self, s: usize) -> Vector3 {
        self.surface[s][0]
    }
    /// `SrfDu`.
    pub fn du(&self, s: usize) -> Vector3 {
        self.surface[s][3]
    }
    /// `SrfDv`.
    pub fn dv(&self, s: usize) -> Vector3 {
        self.surface[s][1]
    }
    /// `SrfNorm`.
    pub fn normal(&self, s: usize) -> Vector3 {
        self.surface[s][8]
    }

    /// `ClassifyPointParams`: each parameter against its domain, within the
    /// 3D tolerance mapped by the first derivative (at most a tenth of the
    /// interval).
    fn classify(&mut self, uvs: [Vector2; 2], domains: [Extent2; 2], tolerance: f64) {
        for s in 0..2 {
            let (u, v) = (domains[s].u_interval(), domains[s].v_interval());
            let mut u_size = self.du(s).length();
            if u_size < EFF_ZERO {
                u_size = EFF_ZERO;
            }
            let mut v_size = self.dv(s).length();
            if v_size < EFF_ZERO {
                v_size = EFF_ZERO;
            }
            let mut tol_u = min(u.length() / 10.0, tolerance / max(0.99, u_size));
            let mut tol_v = min(v.length() / 10.0, tolerance / max(0.99, v_size));
            if tol_u < EFF_ZERO {
                tol_u = EFF_ZERO;
            }
            if tol_v < EFF_ZERO {
                tol_v = EFF_ZERO;
            }
            self.param_types[s] = [
                u.classify_point(uvs[s].x, tol_u),
                v.classify_point(uvs[s].y, tol_v),
            ];
        }
    }

    /// Whether the curve tangent leaves some side of a domain it lies on
    /// (NaN derivatives count as leaving, as in SMLib's negated tests).
    #[allow(clippy::neg_cmp_op_on_partial_ord)]
    fn points_outward(&self) -> bool {
        (0..2).any(|s| {
            let [tu, tv] = self.param_types[s];
            let d = self.uv_deriv[s];
            (tu == PointType::Start && !(d.x >= -EFF_ZERO))
                || (tu == PointType::End && !(d.x <= EFF_ZERO))
                || (tv == PointType::Start && !(d.y >= -EFF_ZERO))
                || (tv == PointType::End && !(d.y <= EFF_ZERO))
        })
    }
}

fn singular(intersector: &Intersector<'_>, s: usize, uv: Vector2) -> bool {
    // SMLib's IsSingularity reads a failed test as singular.
    is_singularity(intersector.surfaces[s], uv, EFF_ZERO, ZONE_TOL_3D, false)
        .map_or(true, |found| found.is_some())
}

/// `SmLocalSurfaceType`, from `smsurf_ComputeSurfacePointType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LocalSurface {
    Planar,
    Elliptical,
    Hyperbolic,
    Cylindrical,
}

/// `smsurf_ComputeSurfacePointType`: how the iso-curves bend against the
/// normal at a point, from the parametric second derivatives. SMLib reads
/// only `N . S_uu` and `N . S_vv`, so a twisted point with only a mixed term
/// (`z = u v`, a saddle) passes for planar and its meeting with a plane is
/// called a tangent point. The port also requires `N . S_uv` to vanish for a
/// planar point and calls such a point hyperbolic (deliberately unlike
/// SMLib).
fn surface_point_type(normal: Vector3, duu: Vector3, dvv: Vector3, duv: Vector3) -> LocalSurface {
    let tolerance = EFF_ZERO_SQRT * 10.0 * (1.0 + normal.length() + duu.length() + dvv.length());
    // 0 on, 1 above, 2 below.
    let side = |dot: f64| {
        if dot.abs() < tolerance {
            0
        } else if dot >= tolerance {
            1
        } else {
            2
        }
    };
    let (along_u, along_v) = (side(normal.dot(duu)), side(normal.dot(dvv)));
    if along_u == 0 && along_v == 0 {
        let twist_tolerance = EFF_ZERO_SQRT * 10.0 * (1.0 + normal.length() + duv.length());
        if normal.dot(duv).abs() < twist_tolerance {
            LocalSurface::Planar
        } else {
            LocalSurface::Hyperbolic
        }
    } else if along_u == 0 || along_v == 0 {
        LocalSurface::Cylindrical
    } else if along_u == along_v {
        LocalSurface::Elliptical
    } else {
        LocalSurface::Hyperbolic
    }
}

/// The principal curvatures and directions of `smsurf_EvaluateGeometric`.
struct Principal {
    k1: f64,
    k2: f64,
    v1: Vector3,
    v2: Vector3,
}

/// `smsurf_EvaluateGeometric` on an `SmTsectPnt` matrix: the principal
/// curvatures (`k1` the larger) and unit directions, with `v1 x v2` along
/// the normal. Like SMLib it replaces a vanishing first derivative in the
/// matrix by one perpendicular to the other and the normal.
fn evaluate_geometric(cells: &mut [Vector3; 9]) -> Result<Principal, Fault> {
    let mut normal = cells[3].cross(cells[1]);
    if normal.length_squared() < EFF_ZERO_SQ {
        normal = cells[8];
        if normal.length_squared() < EFF_ZERO_SQ {
            return Err(Fault::Degenerate);
        }
    }
    normal.unitize()?;
    if cells[3].length_squared() < EFF_ZERO_SQ {
        cells[3] = cells[1].cross(normal);
    }
    if cells[1].length_squared() < EFF_ZERO_SQ {
        cells[1] = normal.cross(cells[3]);
    }
    let (du, dv) = (cells[3], cells[1]);
    let e = du.dot(du);
    let f = du.dot(dv);
    let g = dv.dot(dv);
    let l = cells[6].dot(normal);
    let m = cells[4].dot(normal);
    let n = cells[2].dot(normal);
    let numer1 = l * n - m * m;
    let denom1 = e * g - f * f;
    if denom1.abs() <= EFF_ZERO * numer1.abs() {
        return Err(Fault::Degenerate);
    }
    let gaussian = numer1 / denom1;
    let numer2 = e * n - 2.0 * f * m + g * l;
    let denom2 = 2.0 * (e * g - f * f);
    if denom2.abs() <= EFF_ZERO * numer2.abs() {
        return Err(Fault::Degenerate);
    }
    let mean = numer2 / denom2;
    let mut spread = mean * mean - gaussian;
    if spread < EFF_ZERO_SQ {
        spread = 0.0;
    }
    let root = crate::scalar::sqrt(spread);
    let (k1, k2) = (mean + root, mean - root);
    if (k1 - k2).abs() < EFF_ZERO_SQRT * (1.0 + k1.abs() + k2.abs()) {
        let mut v1 = du;
        v1.unitize()?;
        let mut v2 = v1.cross(dv).cross(v1);
        v2.unitize()?;
        return Ok(Principal { k1, k2, v1, v2 });
    }
    let (count, roots) =
        crate::geometry::solve_quadratic([e * m - f * l, e * n - g * l, f * n - g * m], EFF_ZERO);
    let (mut v1, mut v2);
    if count < 1 {
        if du.max_dimension() > dv.max_dimension() {
            v1 = du;
            v2 = normal.cross(v1);
        } else {
            v1 = dv;
            v2 = v1.cross(normal);
        }
    } else {
        v1 = du + roots[0] * dv;
        if v1.length_squared() < EFF_ZERO_SQ {
            v1 = du;
        }
        v2 = normal.cross(v1);
    }
    let (wu, wv) = (v1.dot(du), v1.dot(dv));
    let a = (g * wu - f * wv) / denom1;
    let b = (-f * wu + e * wv) / denom1;
    let length = e * a * a + 2.0 * f * a * b + g * b * b;
    if length > EFF_ZERO_SQ {
        let normal_curvature = (l * a * a + 2.0 * m * a * b + n * b * b) / length;
        if (normal_curvature - k2).abs() < (normal_curvature - k1).abs() {
            std::mem::swap(&mut v1, &mut v2);
        }
    }
    if v1.length_squared() > EFF_ZERO_SQ {
        v1.unitize()?;
    }
    if v2.length_squared() > EFF_ZERO_SQ {
        v2.unitize()?;
    }
    if v1.cross(v2).dot(normal) < 0.0 {
        v2 = -v2;
    }
    Ok(Principal { k1, k2, v1, v2 })
}

/// The Dupin indicatrix of `SmAdvSurfaceIntersector::MakeDupin`: the
/// principal curvatures signed against a reference normal (one negligible
/// beside the other made zero) with their directions, and the indicatrix
/// conics where the surface bends towards the normal (`above`) or away.
struct Dupin {
    above: Vec<Conic>,
    below: Vec<Conic>,
    k1: f64,
    k2: f64,
    v1: Vector3,
    v2: Vector3,
}

/// `MakeDupin` on an `SmTsectPnt` matrix, which it updates in place like
/// [`evaluate_geometric`].
fn make_dupin(cells: &mut [Vector3; 9], reference: Vector3) -> Result<Dupin, Fault> {
    let principal = evaluate_geometric(cells)?;
    let (mut k1, mut k2) = (principal.k1, principal.k2);
    if reference.dot(cells[8]) < 0.0 {
        k1 = -k1;
        k2 = -k2;
    }
    let mut dupin = Dupin {
        above: Vec::new(),
        below: Vec::new(),
        k1,
        k2,
        v1: principal.v1,
        v2: principal.v2,
    };
    if k1.abs() < EFF_ZERO && k2.abs() < EFF_ZERO {
        return Ok(dupin);
    }
    if k1.abs() < k2.abs() / 1000.0 {
        k1 = 0.0;
    }
    if k2.abs() < k1.abs() / 1000.0 {
        k2 = 0.0;
    }
    dupin.k1 = k1;
    dupin.k2 = k2;
    let radius = |k: f64| 1.0 / crate::scalar::sqrt(k.abs());
    let place = |dupin: &mut Dupin, k: f64, conic: Conic| {
        if k > 0.0 {
            dupin.above.push(conic);
        } else {
            dupin.below.push(conic);
        }
    };
    let turned = |kind: ConicKind, a: f64, b: f64, angle: f64| {
        let mut conic = Conic::new(kind, a, b);
        conic.transform_2d(angle, Vector2::new(0.0, 0.0));
        conic
    };
    if k1.abs() < EFF_ZERO || k2.abs() < EFF_ZERO {
        let (k, angle) = if k1.abs() < EFF_ZERO {
            (k2, 0.0)
        } else {
            (k1, PI / 2.0)
        };
        place(
            &mut dupin,
            k,
            turned(ConicKind::Line, radius(k), 0.0, angle),
        );
        place(
            &mut dupin,
            k,
            turned(ConicKind::Line, radius(k), 0.0, PI + angle),
        );
    } else if k1 * k2 > 0.0 {
        place(
            &mut dupin,
            k1,
            Conic::new(ConicKind::Ellipse, radius(k1), radius(k2)),
        );
    } else if k1 * k2 < 0.0 {
        place(
            &mut dupin,
            k1,
            Conic::new(ConicKind::Hyperbola, radius(k1), radius(k2)),
        );
        place(
            &mut dupin,
            k1,
            turned(ConicKind::Hyperbola, radius(k1), radius(k2), PI),
        );
        place(
            &mut dupin,
            k2,
            turned(ConicKind::Hyperbola, radius(k2), radius(k1), PI / 2.0),
        );
        // SMLib builds the last branch from the second curvature alone.
        place(
            &mut dupin,
            k2,
            turned(ConicKind::Hyperbola, radius(k2), radius(k2), PI + PI / 2.0),
        );
    }
    Ok(dupin)
}

fn angle_between_normals(point: &TsectPoint) -> Result<f64, Fault> {
    let mut angle = point.normal(0).angle_between(point.normal(1))?;
    if angle > PI / 2.0 {
        angle = PI - angle;
    }
    Ok(angle)
}

impl Intersector<'_> {
    /// `SmAdvSurfaceIntersector::ComputeSingularityPoint`: the kind of a
    /// point where the normals are parallel, and the direction to trace.
    /// With one surface planar there, the other's principal curvatures
    /// decide: a tangent curve along a flat principal direction, a
    /// singularity along a saddle's asymptotic directions, or a tangent
    /// point. Otherwise the Dupin indicatrices of the two surfaces decide:
    /// matching curvatures give coincidence or tangent curves along the
    /// shared principal directions, and the indicatrices' meetings give the
    /// directions of the branches (none: a tangent point; two: a tangent
    /// curve; more: a singularity, where the branch inside both domains, or
    /// nearest the previous direction, is taken).
    #[allow(clippy::too_many_lines)]
    fn singularity_point(
        &self,
        point: &mut TsectPoint,
        previous: Option<&TsectPoint>,
        curve_end: Option<Vector3>,
        trace_direction: f64,
    ) -> Result<(), SolverRefusal> {
        if let Some(previous) = previous {
            point.curve_tangent = previous.curve_tangent;
        }
        let kinds = [0, 1].map(|s| {
            surface_point_type(
                point.normal(s),
                point.surface[s][6],
                point.surface[s][2],
                point.surface[s][4],
            )
        });
        // The reference direction: the previous point's (or the traced
        // curve's last) tangent projected into the first tangent plane.
        let reference = previous
            .map(|p| p.curve_tangent)
            .or(curve_end)
            .map(|tangent| {
                let end = point.curve_point + tangent;
                let on_plane = crate::geometry::point_project_to_plane(
                    end,
                    point.curve_point,
                    point.normal(0),
                );
                let mut direction = on_plane - point.curve_point;
                if direction.unitize().is_ok() {
                    direction
                } else {
                    tangent
                }
            });
        let angle_tolerance = self.angle_tolerance;
        // `ALIGN_DIR_WITH_PREVIOUS_POINT`; a failed angle leaves SMLib's 0.
        let align = |point: &mut TsectPoint| {
            if let Some(reference) = reference {
                let angle = point.curve_tangent.angle_between(reference).unwrap_or(0.0);
                if (angle - PI).abs() < min(PI / 30.0, 5.0 * angle_tolerance) {
                    point.curve_tangent = -point.curve_tangent;
                }
            }
        };
        let tangent_curve = |point: &mut TsectPoint, direction: Vector3| {
            point.kind = PointKind::TangentCurve;
            point.curve_tangent = direction * trace_direction;
            align(point);
        };
        // The four principal directions of matching indicatrices: the one
        // nearest the reference, or without one a coincidence along the
        // flatter direction.
        let four = |point: &mut TsectPoint, v1: Vector3, v2: Vector3, k1: f64, k2: f64| {
            let directions = [v1, -v1, v2, -v2];
            point.curve_tangent = directions[0];
            if let Some(reference) = reference {
                let mut best = reference.angle_between(directions[0])?;
                for &direction in &directions[1..] {
                    let test = reference.angle_between(direction)?;
                    if test < best {
                        best = test;
                        point.curve_tangent = direction;
                    }
                }
                point.kind = PointKind::Singularity;
            } else {
                point.curve_tangent = if k2.abs() < k1.abs() { v2 } else { v1 } * trace_direction;
                point.kind = PointKind::Coincidence;
            }
            Ok::<(), Fault>(())
        };
        let mut directions: Vec<Vector3> = Vec::new();
        let planar = if kinds[0] == LocalSurface::Planar {
            Some(0)
        } else if kinds[1] == LocalSurface::Planar {
            Some(1)
        } else {
            None
        };
        if let Some(planar) = planar {
            let other = 1 - planar;
            let mut principal = evaluate_geometric(&mut point.surface[other])?;
            let opposite = point.normal(other).dot(point.normal(planar)) < EFF_ZERO;
            if opposite {
                principal.v1 = -principal.v1;
                principal.v2 = -principal.v2;
            }
            let Principal { k1, k2, v1, v2 } = principal;
            if k1.abs() < EFF_ZERO_SQRT && k2.abs() > EFF_ZERO_SQRT {
                tangent_curve(point, v1);
                return Ok(());
            }
            if k2.abs() < EFF_ZERO_SQRT && k1.abs() > EFF_ZERO_SQRT {
                tangent_curve(point, v2);
                return Ok(());
            }
            #[allow(clippy::neg_cmp_op_on_partial_ord)]
            let saddle = k1 * k2 < 0.0;
            if !saddle {
                point.kind = PointKind::TangentPoint;
                return Ok(());
            }
            point.kind = PointKind::Singularity;
            let y = 10000.0;
            let mut x2 = (1.0 - k1 * y * y) / k2;
            if x2 < 0.0 {
                x2 = 0.0;
            }
            let x = crate::scalar::sqrt(x2);
            let mut out = y * v1 + x * v2;
            let mut inward = y * v1 - x * v2;
            if opposite {
                std::mem::swap(&mut out, &mut inward);
            }
            // SMLib ignores these statuses.
            let _ = out.unitize();
            let _ = inward.unitize();
            directions = vec![
                out * trace_direction,
                -out * trace_direction,
                inward * trace_direction,
                -inward * trace_direction,
            ];
        } else {
            let reference_normal = point.normal(0);
            let first = make_dupin(&mut point.surface[0], reference_normal)?;
            let mut second = make_dupin(&mut point.surface[1], reference_normal)?;
            let (k1, k2, v1, v2) = (first.k1, first.k2, first.v1, first.v2);
            let (k1_2, k2_2, v1_2, v2_2) = (second.k1, second.k2, second.v1, second.v2);
            let near = |a: f64, b: f64| (a - b).abs() < EFF_ZERO_SQRT;
            // The first surface is umbilic there.
            if (k1 - k2).abs() < EFF_ZERO * (1.0 + k1.abs()) {
                if near(k1, k1_2) {
                    if near(k2, k2_2) {
                        point.kind = PointKind::Coincidence;
                    } else {
                        tangent_curve(point, v1_2);
                    }
                    return Ok(());
                }
                if near(k2, k2_2) {
                    tangent_curve(point, v2_2);
                    return Ok(());
                }
            }
            // The second surface is umbilic there.
            if (k1_2 - k2_2).abs() < EFF_ZERO * (1.0 + k1_2.abs()) {
                if near(k1, k1_2) {
                    if near(k2, k2_2) {
                        point.kind = PointKind::Coincidence;
                    } else {
                        tangent_curve(point, v1);
                    }
                    return Ok(());
                }
                if near(k2, k2_2) {
                    tangent_curve(point, v2);
                    return Ok(());
                }
            }
            // The first principal directions agree.
            if v1.cross(v1_2).length() < EFF_ZERO_SQRT * 100.0 {
                if near(k1, k1_2) {
                    if near(k2, k2_2) {
                        four(point, v1, v2, k1, k2)?;
                    } else {
                        tangent_curve(point, v1);
                    }
                    return Ok(());
                }
                if near(k2, k2_2) {
                    tangent_curve(point, v2);
                    return Ok(());
                }
            }
            // The first direction of one is the second of the other.
            if v1.cross(v2_2).length() < EFF_ZERO_SQRT * 100.0 {
                if near(k1, k2_2) {
                    if near(k2, k1_2) {
                        four(point, v1, v2, k1, k2)?;
                    } else {
                        tangent_curve(point, v1);
                    }
                    return Ok(());
                }
                if near(k2, k1_2) {
                    tangent_curve(point, v2);
                    return Ok(());
                }
            }
            let angle = crate::surface_curves::ccw_angle_between(point.normal(0), v1, v1_2)?;
            for conic in second.above.iter_mut().chain(second.below.iter_mut()) {
                conic.transform_2d(angle, Vector2::new(0.0, 0.0));
            }
            for (ones, twos) in [(&first.above, &second.above), (&first.below, &second.below)] {
                for one in ones {
                    for two in twos {
                        for (t, _, _) in one.intersect_conic(two, EFF_ZERO_SQRT * 100.0) {
                            let p = one.evaluate_point(t);
                            let mut p = Vector3::new(p.x, p.y, 0.0);
                            p.unitize()?;
                            let mut direction = p.x * v1 + p.y * v2;
                            direction.unitize()?;
                            // `sm_AddSingDir`.
                            if !directions.iter().any(|d| d.dot(direction) > 0.9) {
                                directions.push(direction);
                            }
                        }
                    }
                }
            }
        }
        match directions.len() {
            0 => {
                point.kind = PointKind::TangentPoint;
                Ok(())
            }
            2 => {
                point.kind = PointKind::TangentCurve;
                point.curve_tangent = directions[0];
                if let Some(reference) = reference {
                    let first = reference.angle_between(directions[0])?;
                    let second = reference.angle_between(directions[1])?;
                    if second < first {
                        point.curve_tangent = directions[1];
                    }
                }
                align(point);
                Ok(())
            }
            n if n > 2 => {
                point.kind = PointKind::Singularity;
                point.curve_tangent = directions[0];
                if curve_end.is_none() {
                    let mut inside = Vec::new();
                    let tolerance = EFF_ZERO_SQRT / 100_000.0;
                    for (index, &direction) in directions.iter().enumerate() {
                        let along =
                            [0, 1].map(|s| drop_vector(point.du(s), point.dv(s), direction));
                        let along = [along[0]?, along[1]?];
                        let probe = |s: usize| point.uv[s] + EFF_ZERO_SQRT * along[s];
                        if self.domains[0].contains_point(probe(0), tolerance)
                            && self.domains[1].contains_point(probe(1), tolerance)
                        {
                            inside.push(index);
                        }
                    }
                    match (reference, inside.first()) {
                        (Some(reference), Some(&first)) => {
                            // SMLib leaves the first branch found when the
                            // first inside branch is the nearest; the port
                            // takes that inside branch.
                            point.curve_tangent = directions[first];
                            let mut best = reference.angle_between(directions[first])?;
                            for &index in &inside[1..] {
                                let test = reference.angle_between(directions[index])?;
                                if test < best {
                                    point.curve_tangent = directions[index];
                                    best = test;
                                }
                            }
                        }
                        (Some(reference), None) => point.curve_tangent = reference,
                        (None, Some(&first)) => point.curve_tangent = directions[first],
                        (None, None) => {}
                    }
                } else {
                    let reference = reference.ok_or(Fault::Degenerate)?;
                    let mut best = reference.angle_between(directions[0])?;
                    for &direction in &directions[1..] {
                        for candidate in [direction, -direction] {
                            let test = reference.angle_between(candidate)?;
                            if test < best {
                                best = test;
                                point.curve_tangent = candidate;
                            }
                        }
                    }
                }
                Ok(())
            }
            _ => Err(Fault::Degenerate.into()),
        }
    }

    /// `ComputePointValues`: evaluates both surfaces at `uvs` and fills the
    /// intersection point: positions, derivatives and unit normals, the
    /// parameter classes, the curve point and unit tangent along the
    /// normals' cross product (turned into the domains, or along the
    /// previous point's direction), its projections into the domains and the
    /// angle between the normals; near parallel normals after `previous` it
    /// refines the point by a normal solve.
    pub fn compute_point_values(
        &self,
        uvs: [Vector2; 2],
        point: &mut TsectPoint,
        previous: Option<&TsectPoint>,
        curve_end: Option<Vector3>,
        step_size: Option<&mut f64>,
        trace_direction: f64,
    ) -> Result<(), SolverRefusal> {
        for s in 0..2 {
            point.uv[s] = uvs[s];
            let cells = evaluate(
                self.surfaces[s],
                uvs[s],
                2,
                2,
                [true, true],
                true,
                true,
                true,
            )?;
            for (slot, value) in point.surface[s].iter_mut().zip(&cells) {
                *slot = *value;
            }
            let mut normal = point.du(s).cross(point.dv(s));
            if normal.length_squared() < EFF_ZERO_SQ {
                normal = evaluate_normal(self.surfaces[s], uvs[s], [true, true]);
                if normal.length_squared() < EFF_ZERO_SQ {
                    return Err(Fault::Degenerate.into());
                }
            }
            let _ = normal.unitize();
            point.surface[s][8] = normal;
        }
        point.classify(uvs, self.domains, self.approx_tolerance);
        point.curve_point = (point.position(0) + point.position(1)).divided(2.0);
        point.deviation = point.curve_point.distance(point.position(0));
        if point.kind != PointKind::Singularity {
            point.kind = PointKind::Crossing;
            point.curve_tangent = point.normal(0).cross(point.normal(1)) * trace_direction;
            point.tangent_plane_angle = angle_between_normals(point)?;
            match previous {
                None => {
                    for s in 0..2 {
                        point.uv_deriv[s] =
                            drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
                    }
                    if point.points_outward() {
                        point.uv_deriv = [-point.uv_deriv[0], -point.uv_deriv[1]];
                        point.curve_tangent = -point.curve_tangent;
                    }
                }
                Some(previous) => {
                    let same = if previous.kind == PointKind::Crossing {
                        let last = 0.0
                            < previous
                                .curve_tangent
                                .dot(previous.normal(0).cross(previous.normal(1)));
                        let this = 0.0
                            < point
                                .curve_tangent
                                .dot(point.normal(0).cross(point.normal(1)));
                        last == this
                    } else {
                        0.0 < previous.curve_tangent.dot(point.curve_tangent)
                    };
                    if !same {
                        point.curve_tangent = -point.curve_tangent;
                    }
                }
            }
            if point.tangent_plane_angle < self.start_parallel_tolerance
                || point.curve_tangent.length() < self.start_parallel_tolerance
                || previous.is_some_and(|p| p.kind == PointKind::TangentCurve)
            {
                match previous {
                    Some(previous) if singular(self, 0, uvs[0]) || singular(self, 1, uvs[1]) => {
                        point.curve_tangent = point
                            .normal(0)
                            .cross(previous.curve_tangent)
                            .cross(point.normal(0));
                        point.kind = PointKind::TangentCurve;
                    }
                    _ => {
                        // SMLib fails the whole point when this fails.
                        self.singularity_point(point, previous, curve_end, trace_direction)
                            .map_err(|refusal| match refusal {
                                SolverRefusal::Failed(_) => Fault::Degenerate.into(),
                                other => other,
                            })?;
                        if matches!(point.kind, PointKind::TangentPoint | PointKind::Coincidence) {
                            return Ok(());
                        }
                    }
                }
            }
        }
        let tangent_length = point.curve_tangent.length();
        if tangent_length < self.start_parallel_tolerance {
            point.kind = PointKind::TangentPoint;
            return Ok(());
        }
        point.curve_tangent /= tangent_length;
        if let (Some(_), Some(_)) = (previous, step_size) {
            for s in 0..2 {
                if (point.du(s).length_squared() < EFF_ZERO_SQRT
                    || point.dv(s).length_squared() < EFF_ZERO_SQRT)
                    && singular(self, s, uvs[s])
                {
                    // `StepBackFromSingularity`.
                    return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
                }
            }
        }
        for s in 0..2 {
            point.uv_deriv[s] = drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
        }
        point.tangent_plane_angle = angle_between_normals(point)?;
        if let Some(previous) = previous
            && point.kind == PointKind::Crossing
            && point.tangent_plane_angle < 3.0 * PI / 180.0
        {
            let mut found = None;
            if previous.kind != PointKind::Unknown && previous.kind != PointKind::Singularity {
                let grow = |s: usize| {
                    let further = point.uv[s] + (point.uv[s] - previous.uv[s]);
                    Extent2::new(previous.uv[s], previous.uv[s]).union(Extent2::new(
                        self.domains[s].clamp_point(further),
                        self.domains[s].clamp_point(further),
                    ))
                };
                let (domain1, domain2) = (grow(0), grow(1));
                found = local_surface_solve(
                    self.surfaces[0],
                    domain1,
                    self.surfaces[1],
                    domain2,
                    SolverOperation::Normalize,
                    domain1.evaluate(0.5, 0.5),
                    domain2.evaluate(0.5, 0.5),
                )?;
            }
            if let Some(pair) = found {
                point.kind = PointKind::Unknown;
                self.compute_point_values(
                    [pair.uv, pair.other_uv],
                    point,
                    None,
                    curve_end,
                    None,
                    trace_direction,
                )?;
                if previous.curve_tangent.dot(point.curve_tangent) < 0.0 {
                    point.curve_tangent = -point.curve_tangent;
                    for s in 0..2 {
                        point.uv_deriv[s] =
                            drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
                    }
                }
            } else if point.kind == PointKind::Crossing
                && point.tangent_plane_angle < self.start_parallel_tolerance
            {
                point.curve_tangent = previous.curve_tangent;
                let _ = point.curve_tangent.unitize();
                for s in 0..2 {
                    point.uv_deriv[s] = drop_vector(point.du(s), point.dv(s), point.curve_tangent)?;
                }
                return Ok(());
            }
        }
        Ok(())
    }
}

/// What [`Intersector::compute_next_point`] reports besides the point:
/// SMLib's `rbFoundGoodPoint`, `rdDeviationFound`, `rdAngleFoundRad`,
/// `rbClipped` and `rbBoundaryHit`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct StepOutcome {
    /// The step met the angle and distance tolerances.
    pub good: bool,
    /// The largest distance found along the span.
    pub deviation: f64,
    /// Whether the span test wrote `deviation`; SMLib's caller otherwise
    /// keeps its previous value.
    pub deviation_written: bool,
    /// The turn of the curve tangent over the span.
    pub angle: f64,
    /// A domain side shortened the step.
    pub clipped: bool,
    /// The step stopped on a domain side.
    pub boundary_hit: bool,
}

impl Intersector<'_> {
    /// `ComputeStepSize`: the old step, or a hundred tolerances to start.
    pub fn compute_step_size(&self, old_step: f64) -> f64 {
        if crate::scalar::is_zero(old_step) {
            100.0 * self.approx_tolerance
        } else {
            old_step
        }
    }

    /// `RefinePoint`: the meeting near `uvs` (looser for tangent and
    /// singular points); a tangent point is not refined and SMLib then reads
    /// an untouched solution, the parameters `(0, 0)` clamped into each
    /// domain.
    pub fn refine_point(
        &self,
        uvs: [Vector2; 2],
        kind: PointKind,
        plane_normal: Option<Vector3>,
        previous_uv: Option<Vector2>,
        boundary_testing: bool,
    ) -> Result<Option<[Vector2; 2]>, SolverRefusal> {
        let tolerance = match kind {
            PointKind::TangentCurve | PointKind::Singularity => self.approx_tolerance * 10.00,
            PointKind::Crossing => self.approx_tolerance,
            PointKind::TangentPoint => {
                let zero = Vector2::new(0.0, 0.0);
                return Ok(Some([
                    self.domains[0].clamp_point(zero),
                    self.domains[1].clamp_point(zero),
                ]));
            }
            _ => return Ok(None),
        };
        let Some(pair) = crate::surface_surface::local_surface_intersect(
            self.surfaces[0],
            self.domains[0],
            self.surfaces[1],
            self.domains[1],
            tolerance,
            uvs[0],
            uvs[1],
            plane_normal,
            previous_uv,
            boundary_testing,
        )?
        else {
            return Ok(None);
        };
        Ok(Some([
            self.domains[0].clamp_point(pair.uv),
            self.domains[1].clamp_point(pair.other_uv),
        ]))
    }

    /// `TestSpanAccuracy`: the span from `start` to `next` is good when the
    /// tangent turns less than the angle tolerance and, at five points of
    /// the Hermite curves through both ends, the two surfaces and the 3D
    /// curve stay within the tolerance of one another. A failing span
    /// demotes a singular `next` to crossing.
    pub fn test_span_accuracy(
        &self,
        start: &TsectPoint,
        next: &mut TsectPoint,
    ) -> Result<(bool, Option<f64>, f64), SolverRefusal> {
        use crate::hermite::Hermite;
        if next.curve_tangent.length() < self.stepping_parallel_tolerance {
            return Err(Fault::Degenerate.into());
        }
        let step = next.curve_parameter - start.curve_parameter;
        let angle = start.curve_tangent.angle_between(next.curve_tangent)?;
        let demote = |next: &mut TsectPoint| {
            if next.kind == PointKind::Singularity {
                next.kind = PointKind::Crossing;
            }
        };
        if angle > self.angle_tolerance && step > EFF_ZERO * 100.0 {
            demote(next);
            return Ok((false, None, angle));
        }
        let uv_curves = [0, 1].map(|s| {
            Hermite::planar(
                start.uv[s],
                start.uv_deriv[s] * step,
                next.uv[s],
                next.uv_deriv[s] * step,
            )
        });
        let curve = Hermite::new(
            start.curve_point,
            start.curve_tangent * step,
            next.curve_point,
            next.curve_tangent * step,
        );
        let tolerance_sq = self.approx_tolerance * self.approx_tolerance;
        let mut max_sq: f64 = 0.0;
        for i in 0..5 {
            let t = (i as f64 + 1.0) / (5.0 + 1.0);
            let on = |s: usize| -> Result<Vector3, SolverRefusal> {
                let q = uv_curves[s].evaluate_point(t);
                let uv = self.domains[s].clamp_point(Vector2::new(q.x, q.y));
                Ok(crate::surface_eval::evaluate_point(self.surfaces[s], uv)?)
            };
            let (p0, p1) = (on(0)?, on(1)?);
            let c = curve.evaluate_point(t);
            let distance_sq = crate::scalar::max(
                crate::scalar::max(p0.distance_squared(p1), p0.distance_squared(c)),
                p1.distance_squared(c),
            );
            if distance_sq > max_sq {
                max_sq = distance_sq;
            }
            if distance_sq > tolerance_sq && step > EFF_ZERO * 100.0 {
                demote(next);
                return Ok((false, Some(crate::scalar::sqrt(max_sq)), angle));
            }
        }
        Ok((true, Some(crate::scalar::sqrt(max_sq)), angle))
    }

    /// `ComputeNextPoint`: steps `step` along the previous point's tangent
    /// (clipped to the domains), refines the meeting there, computes its
    /// values and tests the span.
    pub fn compute_next_point(
        &self,
        previous: &TsectPoint,
        step: f64,
        next: &mut TsectPoint,
        curve_end: Option<Vector3>,
        trace_direction: f64,
    ) -> Result<StepOutcome, SolverRefusal> {
        let mut outcome = StepOutcome::default();
        if !self.domains[0].contains_point(previous.uv[0], EFF_ZERO)
            || !self.domains[1].contains_point(previous.uv[1], EFF_ZERO)
        {
            return Err(Fault::Degenerate.into());
        }
        let unclipped = step;
        let mut step = self.domains[0].clip_line(previous.uv[0], previous.uv_deriv[0], step);
        step = self.domains[1].clip_line(previous.uv[1], previous.uv_deriv[1], step);
        outcome.clipped = step < unclipped - EFF_ZERO * (1.0 + unclipped);
        if crate::scalar::is_zero(step) {
            outcome.boundary_hit = true;
            return Ok(outcome);
        }
        let guess = [0, 1]
            .map(|s| self.domains[s].clamp_point(previous.uv[s] + step * previous.uv_deriv[s]));
        let refine = |boundary: bool| {
            self.refine_point(
                guess,
                previous.kind,
                Some(previous.curve_tangent),
                Some(previous.uv[0]),
                boundary,
            )
        };
        // SMLib treats a failed refinement like a missed one.
        let Some(mut uvs) = refine(true).or_else(|r| match r {
            SolverRefusal::Failed(_) => Ok(None),
            other => Err(other),
        })?
        else {
            return Ok(outcome);
        };
        let scale = |step: &mut f64, uvs: &[Vector2; 2]| {
            let guess_step = (guess[0] - previous.uv[0]).length();
            let actual_step = (uvs[0] - previous.uv[0]).length();
            if guess_step > EFF_ZERO * actual_step {
                *step *= actual_step / guess_step;
            }
            (guess_step, actual_step)
        };
        let (guess_step, actual_step) = scale(&mut step, &uvs);
        if crate::scalar::is_zero(step) || actual_step < guess_step / 100.0 {
            let Some(again) = refine(false).or_else(|r| match r {
                SolverRefusal::Failed(_) => Ok(None),
                other => Err(other),
            })?
            else {
                return Ok(outcome);
            };
            uvs = again;
            scale(&mut step, &uvs);
            if crate::scalar::is_zero(step) {
                outcome.boundary_hit = true;
                return Ok(outcome);
            }
        }
        match self.compute_point_values(
            uvs,
            next,
            Some(previous),
            curve_end,
            Some(&mut step),
            trace_direction,
        ) {
            Ok(()) => {}
            // SMLib reads any failure here as a boundary hit.
            Err(SolverRefusal::Failed(_)) => {
                outcome.boundary_hit = true;
                return Ok(outcome);
            }
            Err(refusal) => return Err(refusal),
        }
        if next.kind == PointKind::TangentPoint || next.kind == PointKind::Coincidence {
            return Ok(outcome);
        }
        let along = crate::geometry::line_closest_point(
            Vector3::new(previous.uv[0].x, previous.uv[0].y, 0.0),
            Vector3::new(previous.uv_deriv[0].x, previous.uv_deriv[0].y, 0.0),
            Vector3::new(next.uv[0].x, next.uv[0].y, 0.0),
        );
        if along < 0.0 {
            if next.kind == PointKind::Singularity {
                next.kind = PointKind::Crossing;
            }
            return Ok(outcome);
        }
        next.curve_parameter =
            previous.curve_parameter + next.curve_point.distance(previous.curve_point);
        let (good, deviation, angle) = self.test_span_accuracy(previous, next)?;
        outcome.good = good;
        if let Some(deviation) = deviation {
            outcome.deviation = deviation;
            outcome.deviation_written = true;
        }
        outcome.angle = angle;
        Ok(outcome)
    }
}

impl crate::curve_solve::PointSolveCurve for crate::hermite::Hermite {
    fn evaluate(
        &self,
        t: f64,
        derivatives: usize,
        _from_left: bool,
        nonzero_tangent: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        crate::hermite::Hermite::evaluate(self, t, derivatives, nonzero_tangent)
    }
    fn evaluate_point(&self, t: f64) -> Result<Vector3, Fault> {
        Ok(crate::hermite::Hermite::evaluate_point(self, t))
    }
    /// `SmCurve::IsDegenerate` over the natural interval.
    fn is_degenerate(&self, tolerance: f64) -> bool {
        let interval = self.natural_interval();
        let p1 = crate::hermite::Hermite::evaluate_point(self, interval.min);
        let p2 = crate::hermite::Hermite::evaluate_point(self, interval.max);
        let scaled_zero = EFF_ZERO * (1.0 + p1.max_dimension());
        let tolerance_sq = max(tolerance * tolerance, scaled_zero * scaled_zero);
        if p1.distance_squared(p2) >= tolerance_sq {
            return false;
        }
        let lsq = interval.length() * interval.length();
        for t in [interval.min, interval.max, interval.evaluate(0.4567)] {
            let Ok(pv) = crate::hermite::Hermite::evaluate(self, t, 1, false) else {
                return false;
            };
            if pv[1].length_squared() * lsq > tolerance_sq {
                return false;
            }
            if t == interval.evaluate(0.4567) && p1.distance_squared(pv[0]) > tolerance_sq {
                return false;
            }
        }
        true
    }
}

/// `sm_IsPointOnCurve`: whether `test` lies on the Hermite span from `prev`
/// to `next`, within `tolerance` grown for shallow meetings, and either
/// against the span's direction or inside it within half the span.
pub fn is_point_on_span(
    test: &TsectPoint,
    prev: &TsectPoint,
    next: &TsectPoint,
    tolerance: f64,
) -> Result<bool, SolverRefusal> {
    use crate::hermite::Hermite;
    let mut tolerance = tolerance;
    let angle = max(prev.tangent_plane_angle, next.tangent_plane_angle);
    if angle > 0.0 {
        let mut tangent = (angle / 2.0).tan();
        if tangent < 0.001 {
            tangent = 0.001;
        }
        tolerance /= tangent;
    }
    let delta = next.curve_parameter - prev.curve_parameter;
    let (p1, p2) = (prev.curve_point, next.curve_point);
    let span = Hermite::new(
        p1,
        prev.curve_tangent * delta,
        p2,
        next.curve_tangent * delta,
    );
    let length = (p2 - p1).max_dimension();
    let (span_box, _) = span.bounding_boxes();
    let mut point_box = crate::extent::Extent3::singleton(test.curve_point);
    point_box.expand_absolute(tolerance * 2.0 + length / 4.0);
    if span_box.are_disjoint(point_box, 0.0) {
        return Ok(false);
    }
    for s in 0..2 {
        let uv_span = Hermite::planar(
            prev.uv[s],
            prev.uv_deriv[s] * delta,
            next.uv[s],
            next.uv_deriv[s] * delta,
        );
        let (mut uv_box, _) = uv_span.bounding_boxes();
        let diagonal = (uv_box.max - uv_box.min).length();
        for axis in 0..2 {
            let (lo, hi) = if axis == 0 {
                (uv_box.min.x, uv_box.max.x)
            } else {
                (uv_box.min.y, uv_box.max.y)
            };
            let mut grow = hi - lo;
            if grow < diagonal / 100.0 {
                grow = diagonal / 100.0;
            }
            if axis == 0 {
                uv_box.min.x = lo - grow;
                uv_box.max.x = hi + grow;
            } else {
                uv_box.min.y = lo - grow;
                uv_box.max.y = hi + grow;
            }
        }
        let uv = Vector3::new(test.uv[s].x, test.uv[s].y, 0.0);
        if !uv_box.contains_point(uv, 0.0) {
            return Ok(false);
        }
    }
    let interval = span.natural_interval();
    let Some(found) = crate::curve_solve::local_point_solve(
        &span,
        interval,
        SolverOperation::Minimize,
        test.curve_point,
        interval.evaluate(0.5),
    )?
    else {
        return Ok(false);
    };
    if found.distance <= tolerance {
        let pv = span.evaluate(found.parameter, 1, true)?;
        if test.curve_tangent.dot(pv[1]) < 0.0 {
            return Ok(true);
        }
    }
    if interval.contains_value(found.parameter, -EFF_ZERO) {
        return Ok(found.distance <= length / 2.0);
    }
    Ok(false)
}

/// The tracing state of `SmSurfaceIntersector`: the curve being traced, the
/// start and through points, and the trace's flags.
#[derive(Clone, Debug, Default)]
pub struct Tracer {
    /// `m_vCurvePoints`.
    pub curve_points: Vec<TsectPoint>,
    /// `m_vStartPoints`, in list order.
    pub start_points: Vec<TsectPoint>,
    /// `m_vThroughPoints`, in list order.
    pub through_points: Vec<TsectPoint>,
    /// `m_pStartPoint`, which is the first curve point when set.
    pub has_start_point: bool,
    /// `m_pEndPoint`.
    pub end_point: Option<TsectPoint>,
    /// `m_dCurveTraceDirection`.
    pub trace_direction: f64,
    /// `m_dLastThroughParmeter`.
    pub last_through_parameter: f64,
    /// `m_dExtensionDistance`.
    pub extension_distance: f64,
    /// `m_bDoingExtensionClipping`.
    pub extension_clipping: bool,
    /// `m_bCurveIsClosed`.
    pub curve_is_closed: bool,
    /// `m_bClippedByThroughPoint`.
    pub clipped_by_through_point: bool,
}

/// Which list `sm_UpdateTsectPoint` snaps to.
#[derive(Clone, Copy)]
enum Target {
    Through(usize),
    Start(usize),
    FirstCurvePoint,
    End,
}

impl Tracer {
    fn target(&self, target: Target) -> TsectPoint {
        match target {
            Target::Through(i) => self.through_points[i],
            Target::Start(i) => self.start_points[i],
            Target::FirstCurvePoint => self.curve_points[0],
            Target::End => self.end_point.unwrap_or_default(),
        }
    }

    fn target_mut(&mut self, target: Target) -> Option<&mut TsectPoint> {
        match target {
            Target::Through(i) => self.through_points.get_mut(i),
            Target::Start(i) => self.start_points.get_mut(i),
            Target::FirstCurvePoint => self.curve_points.first_mut(),
            Target::End => self.end_point.as_mut(),
        }
    }

    /// `sm_UpdateTsectPoint`: moves `next` onto the target point, keeping the
    /// curve parameter proportional to the distance from the last point.
    fn update_point(
        &mut self,
        intersector: &Intersector<'_>,
        next: &mut TsectPoint,
        target: Target,
    ) -> Result<(), SolverRefusal> {
        let prev = *self.curve_points.last().ok_or(Fault::Degenerate)?;
        let distance = prev.curve_point.distance(next.curve_point);
        let step = next.curve_parameter - prev.curve_parameter;
        let snap = self.target(target);
        intersector.compute_point_values(
            snap.uv,
            next,
            Some(&prev),
            self.curve_points.last().map(|p| p.curve_tangent),
            None,
            self.trace_direction,
        )?;
        let new_distance = prev.curve_point.distance(next.curve_point);
        next.curve_parameter = prev.curve_parameter + new_distance / distance * step;
        if let Some(point) = self.target_mut(target)
            && point.curve_parameter < 0.0
        {
            point.curve_parameter = next.curve_parameter;
        }
        Ok(())
    }

    /// `TestPoint`: snaps `next` to the first through point, a start point,
    /// the trace's own start (closing the curve) or its end point when it
    /// lies on the span from the last point, and reports whether the trace
    /// is done there or because `next` leaves a domain across its side.
    #[allow(clippy::too_many_lines)]
    pub fn test_point(
        &mut self,
        intersector: &Intersector<'_>,
        next: &mut TsectPoint,
    ) -> Result<bool, SolverRefusal> {
        let tolerance = intersector.approx_tolerance;
        let prev_index = self
            .curve_points
            .len()
            .checked_sub(1)
            .ok_or(Fault::Degenerate)?;
        if !self.through_points.is_empty() {
            let prev = self.curve_points[prev_index];
            if is_point_on_span(&self.through_points[0], &prev, next, tolerance)? {
                self.update_point(intersector, next, Target::Through(0))?;
                let through = self.through_points.remove(0);
                if self.extension_clipping && self.through_points.is_empty() {
                    self.last_through_parameter = through.curve_parameter;
                }
                self.clipped_by_through_point = true;
                return Ok(false);
            }
        }
        for i in (0..self.start_points.len()).rev() {
            let prev = self.curve_points[prev_index];
            if self.has_start_point && prev_index == 0 {
                continue;
            }
            let start = self.start_points[i];
            if prev.curve_point.distance(start.curve_point) < tolerance {
                let size1 = (intersector.domains[0].max - intersector.domains[0].min).length();
                let size2 = (intersector.domains[1].max - intersector.domains[1].min).length();
                if prev.uv[0].distance(start.uv[0]) < size1 / 2.0
                    && prev.uv[1].distance(start.uv[1]) < size2 / 2.0
                {
                    continue;
                }
            }
            if is_point_on_span(&start, &prev, next, tolerance)? {
                self.update_point(intersector, next, Target::Start(i))?;
                let start = self.start_points[i];
                if next.kind == PointKind::Crossing
                    && start.tangent_plane_angle > intersector.stepping_parallel_tolerance
                {
                    self.start_points.remove(i);
                }
                if start.kind == PointKind::Crossing {
                    let through = |s: usize, a: usize| {
                        let along = start.param_types[s][a];
                        let across = start.param_types[s][1 - a];
                        let d = if a == 0 {
                            start.uv_deriv[s].x
                        } else {
                            start.uv_deriv[s].y
                        };
                        (along == PointType::Start || along == PointType::End)
                            && across == PointType::Inside
                            && crate::scalar::is_zero(d)
                    };
                    let any = through(0, 0) || through(1, 0) || through(0, 1) || through(1, 1);
                    return Ok(!any);
                }
                return Ok(true);
            }
        }
        if self.has_start_point && prev_index != 0 && prev_index != 1 {
            let prev = self.curve_points[prev_index];
            if is_point_on_span(&self.curve_points[0], &prev, next, tolerance)? {
                self.update_point(intersector, next, Target::FirstCurvePoint)?;
                self.curve_is_closed = true;
                return Ok(true);
            }
        }
        if let Some(end) = self.end_point {
            let prev = self.curve_points[prev_index];
            if is_point_on_span(&end, &prev, next, tolerance)? {
                let on_seam =
                    prev_index == 0 && prev.curve_point.distance(end.curve_point) < tolerance;
                if !on_seam {
                    self.update_point(intersector, next, Target::End)?;
                    return Ok(true);
                }
            }
        } else if self.extension_clipping
            && self.last_through_parameter != 0.0
            && next.curve_parameter > self.last_through_parameter + self.extension_distance
        {
            return Ok(true);
        }
        let mut on_boundary = false;
        for s in 0..2 {
            let domain = intersector.domains[s];
            let size = max(domain.min.max_dimension(), domain.max.max_dimension());
            let tol = EFF_ZERO * (1.0 + size) * 100.0;
            let angle_tol = EFF_ZERO_SQRT * (1.0 + size) * 10.0;
            let (on, binormal) = domain.point_on_boundary(next.uv[s], tol);
            if on
                && ((binormal.x != 0.0 && next.uv_deriv[s].x.abs() > angle_tol)
                    || (binormal.y != 0.0 && next.uv_deriv[s].y.abs() > angle_tol))
            {
                on_boundary = true;
                break;
            }
        }
        if !on_boundary {
            return Ok(false);
        }
        let near = |p: &TsectPoint| {
            (0..2).all(|k| {
                let domain = intersector.domains[k];
                let size = max(domain.min.max_dimension(), domain.max.max_dimension());
                let tol = EFF_ZERO * (1.0 + size) * 100.0;
                (p.uv[k] - next.uv[k]).length_squared() <= tol
            })
        };
        self.start_points.retain(|p| !near(p));
        Ok(true)
    }

    /// `TraceCurve`: steps from `start` along the meeting curve, shrinking a
    /// step that fails the span tests and growing one that passes easily,
    /// until a side, a start or end point, or the curve's own start. Points
    /// go to `curve_points`. A start already on one of `curves`, or a
    /// trace whose second point runs along one, gives no new points.
    #[allow(clippy::too_many_lines)]
    pub fn trace_curve(
        &mut self,
        intersector: &Intersector<'_>,
        start: &TsectPoint,
        curves: &[IntersectionCurve],
    ) -> Result<(), SolverRefusal> {
        let uv_tol = |s: usize| {
            let d = intersector.domains[s];
            EFF_ZERO * (1.0 + max(d.min.max_dimension(), d.max.max_dimension()))
        };
        if self.curve_points.is_empty() {
            let mut first = *start;
            if first.kind == PointKind::TangentPoint || first.kind == PointKind::Coincidence {
                if !intersector.is_point_on_curve(first.curve_point, curves, None)? {
                    first.curve_parameter = 0.0;
                    first.deviation = 0.0;
                    self.curve_points.push(first);
                    self.curve_is_closed = true;
                }
                return Ok(());
            }
            first.curve_parameter = 0.0;
            first.deviation = 0.0;
            self.curve_points.push(first);
            let probe = [0, 1].map(|s| first.uv[s] + (10.0 * EFF_ZERO_SQRT) * first.uv_deriv[s]);
            if !intersector.domains[0].contains_point(probe[0], uv_tol(0))
                || !intersector.domains[1].contains_point(probe[1], uv_tol(1))
            {
                if first.is_bounded()
                    && intersector.is_point_on_curve(first.curve_point, curves, None)?
                {
                    self.curve_points.clear();
                    return Ok(());
                }
                self.trace_direction = -self.trace_direction;
                self.reverse_curve_direction()?;
            }
        }
        self.curve_is_closed = false;
        let mut count = 0usize;
        let mut step = 0.0;
        // SMLib's `dDeviation` lives across the whole trace and keeps its
        // value when a step returns before the span test writes it.
        let mut deviation = 0.0;
        let check = |p: &TsectPoint| -> Result<(), SolverRefusal> {
            if !intersector.domains[0].contains_point(p.uv[0], uv_tol(0))
                || !intersector.domains[1].contains_point(p.uv[1], uv_tol(1))
            {
                return Err(Fault::Degenerate.into());
            }
            Ok(())
        };
        loop {
            let last = *self.curve_points.last().ok_or(Fault::Degenerate)?;
            step = intersector.compute_step_size(step);
            let mut next = TsectPoint::default();
            let mut outcome = match intersector.compute_next_point(
                &last,
                step,
                &mut next,
                self.curve_points.last().map(|p| p.curve_tangent),
                self.trace_direction,
            ) {
                Ok(outcome) => outcome,
                Err(SolverRefusal::Failed(_)) => break,
                Err(other) => return Err(other),
            };
            if outcome.boundary_hit && step < intersector.approx_tolerance * 100.0 {
                break;
            }
            if outcome.good {
                check(&next)?;
            }
            if outcome.deviation_written {
                deviation = outcome.deviation;
            }
            let mut decreased = false;
            while !outcome.good {
                step /= 1.5;
                if step < EFF_ZERO * 100.0 {
                    break;
                }
                decreased = true;
                outcome = match intersector.compute_next_point(
                    &last,
                    step,
                    &mut next,
                    self.curve_points.last().map(|p| p.curve_tangent),
                    self.trace_direction,
                ) {
                    Ok(outcome) => outcome,
                    Err(SolverRefusal::Failed(_)) => {
                        outcome.good = false;
                        break;
                    }
                    Err(other) => return Err(other),
                };
                if outcome.deviation_written {
                    deviation = outcome.deviation;
                }
                if outcome.good {
                    check(&next)?;
                }
            }
            if !outcome.good {
                break;
            }
            let last_parameter = next.curve_parameter;
            let done = self.test_point(intersector, &mut next)?;
            if done || next.curve_parameter != last_parameter {
                self.curve_points.push(next);
                if done {
                    break;
                }
                continue;
            }
            if !outcome.clipped
                && !decreased
                && deviation < intersector.approx_tolerance / 8.0
                && outcome.angle < intersector.angle_tolerance / 4.0
            {
                let mut last_parameter = next.curve_parameter;
                while outcome.good && !outcome.clipped {
                    step *= 1.5;
                    outcome = match intersector.compute_next_point(
                        &last,
                        step,
                        &mut next,
                        self.curve_points.last().map(|p| p.curve_tangent),
                        self.trace_direction,
                    ) {
                        Ok(outcome) => outcome,
                        Err(SolverRefusal::Failed(_)) => {
                            outcome.good = false;
                            break;
                        }
                        Err(other) => return Err(other),
                    };
                    if outcome.deviation_written {
                        deviation = outcome.deviation;
                    }
                    if !outcome.good {
                        continue;
                    }
                    check(&next)?;
                    if next.curve_parameter <= last_parameter + EFF_ZERO {
                        outcome.good = false;
                        break;
                    }
                    last_parameter = next.curve_parameter;
                    let done = self.test_point(intersector, &mut next)?;
                    if done || next.curve_parameter != last_parameter {
                        if done {
                            self.curve_points.push(next);
                            return Ok(());
                        }
                        break;
                    }
                }
                if !outcome.good {
                    step /= 1.5;
                    outcome = intersector.compute_next_point(
                        &last,
                        step,
                        &mut next,
                        self.curve_points.last().map(|p| p.curve_tangent),
                        self.trace_direction,
                    )?;
                    if outcome.deviation_written {
                        deviation = outcome.deviation;
                    }
                    if outcome.good {
                        check(&next)?;
                    }
                }
            }
            self.clipped_by_through_point = false;
            if step < EFF_ZERO * 100.0 {
                break;
            }
            let moved = next.curve_point.distance(last.curve_point);
            if moved < EFF_ZERO * 100.0 * (1.0 + last.curve_point.max_dimension()) {
                break;
            }
            // A second point on an earlier curve: this curve is traced
            // already, unless the step was clipped and its middle is off.
            if count == 0
                && next.kind != PointKind::Singularity
                && intersector.is_point_on_curve(next.curve_point, curves, None)?
            {
                if !outcome.clipped {
                    break;
                }
                let span = crate::hermite::Hermite::new(
                    last.curve_point,
                    last.curve_tangent * step,
                    next.curve_point,
                    next.curve_tangent * step,
                );
                if intersector.is_point_on_curve(span.evaluate_point(0.5), curves, None)? {
                    break;
                }
            }
            let done = self.test_point(intersector, &mut next)?;
            next.deviation = deviation;
            self.curve_points.push(next);
            if done {
                break;
            }
            count += 1;
            if count == 5000 {
                return Err(Fault::Degenerate.into());
            }
        }
        Ok(())
    }
}

/// Whether the point lies on a side of both domains with a positive
/// deviation, as `AddStartPoint` classifies start points.
fn bounded(point: &TsectPoint, scaled_zero: f64) -> bool {
    let on_side = |t: PointType| matches!(t, PointType::Start | PointType::End | PointType::Both);
    point.deviation > scaled_zero
        && (on_side(point.param_types[0][0]) || on_side(point.param_types[0][1]))
        && (on_side(point.param_types[1][0]) || on_side(point.param_types[1][1]))
}

/// `smsurf_EvaluateNormalSection` on the point's own evaluation matrix,
/// which SMLib updates in place, along a parameter direction read as a 3D
/// vector with zero `z` (as SMLib converts it). A failure leaves zero.
fn normal_section_curvature(point: &mut TsectPoint, s: usize) -> f64 {
    let mut jet: crate::curve_surface::SurfaceJet = [[Vector3::default(); 3]; 3];
    for k in 0..3 {
        for l in 0..3 {
            jet[k][l] = point.surface[s][k * 3 + l];
        }
    }
    let d = point.uv_deriv[s];
    let result = crate::curve_surface::normal_section_of(&mut jet, Vector3::new(d.x, d.y, 0.0));
    for k in 0..3 {
        for l in 0..3 {
            point.surface[s][k * 3 + l] = jet[k][l];
        }
    }
    result.map_or(0.0, |(_, curvature)| curvature)
}

impl Tracer {
    /// `AddStartPoint`: adds `start` (crossing points last, others first)
    /// unless an existing start point matches it: the same parameters, a
    /// bounded and an unbounded point close in space and parameters with
    /// tangents within a curvature-scaled angle, or nearby parameters whose
    /// surface points meet.
    pub fn add_start_point(
        &mut self,
        intersector: &Intersector<'_>,
        mut start: TsectPoint,
    ) -> Result<(), SolverRefusal> {
        let tolerance = intersector.approx_tolerance;
        let scaled_zero = EFF_ZERO * 100.0 * (1.0 + start.position(0).max_dimension());
        let bounded_start = bounded(&start, scaled_zero);
        let mut matched = false;
        for i in 0..self.start_points.len() {
            let target = self.start_points[i];
            let bounded_target = bounded(&target, scaled_zero);
            let d1 = (start.uv[0] - target.uv[0]).length_squared();
            let d2 = (start.uv[1] - target.uv[1]).length_squared();
            if d1 < EFF_ZERO_SQ && d2 < EFF_ZERO_SQ {
                matched = true;
                break;
            }
            if bounded_start != bounded_target
                && start.position(0).distance(target.position(0)) < 2.0 * tolerance
                && start.position(1).distance(target.position(1)) < 2.0 * tolerance
            {
                let size = |s: usize| {
                    let d = intersector.surfaces[s].natural_domain();
                    max(d.min.max_dimension(), d.max.max_dimension()) / 10.0
                };
                if start.uv[0].distance(target.uv[0]) < size(0)
                    && start.uv[1].distance(target.uv[1]) < size(1)
                {
                    let c0 = normal_section_curvature(&mut start, 0);
                    let c1 = normal_section_curvature(&mut start, 1);
                    let curvature = max(c0.abs(), c1.abs());
                    let radius = if crate::scalar::is_zero(curvature) {
                        100.0 * tolerance
                    } else {
                        min(1.0 / curvature, 100.0 * tolerance)
                    };
                    let angle_tolerance = 0.01 + tolerance.atan2(radius);
                    let between = start
                        .curve_tangent
                        .angle_between(target.curve_tangent)
                        .unwrap_or(0.0);
                    if between < angle_tolerance {
                        matched = true;
                        break;
                    }
                }
            }
            if d1 < EFF_ZERO && d2 < EFF_ZERO {
                let p1 =
                    crate::surface_eval::evaluate_point(intersector.surfaces[0], target.uv[0])?;
                let p2 =
                    crate::surface_eval::evaluate_point(intersector.surfaces[1], target.uv[1])?;
                if p1.distance(p2) < tolerance {
                    matched = true;
                    break;
                }
            }
        }
        if !matched {
            if start.kind == PointKind::Crossing {
                self.start_points.push(start);
            } else {
                self.start_points.insert(0, start);
            }
        }
        Ok(())
    }
}

impl Intersector<'_> {
    /// `ClassifyBoundaryStartPointIntersection`: whether a meeting found on
    /// a side curve is only a touch, not the start of an intersection curve.
    /// Tangent points touch; points whose curve leaves both domains both ways
    /// touch; tangent, singular and coincident points start curves; and a
    /// crossing on a side tangent to the other surface is probed along the
    /// curve both ways to see whether the gap opens monotonically away from
    /// the side (a touch) or crosses over.
    #[allow(clippy::too_many_lines)]
    pub fn classify_boundary_start(&self, point: &TsectPoint) -> Result<bool, SolverRefusal> {
        if point.kind == PointKind::TangentPoint {
            return Ok(true);
        }
        // SMLib adds and subtracts the scaled derivative; the same here.
        let p_inside = (0..2).all(|s| {
            self.domains[s].contains_point(
                point.uv[s] + EFF_ZERO_SQRT * point.uv_deriv[s],
                10.0 * EFF_ZERO,
            )
        });
        let n_inside = (0..2).all(|s| {
            self.domains[s].contains_point(
                point.uv[s] - EFF_ZERO_SQRT * point.uv_deriv[s],
                10.0 * EFF_ZERO,
            )
        });
        if !p_inside && !n_inside {
            return Ok(true);
        }
        if matches!(
            point.kind,
            PointKind::TangentCurve | PointKind::Singularity | PointKind::Coincidence
        ) {
            return Ok(false);
        }
        let angle_tolerance = 2.0 * PI / 180.0;
        let gap = point.position(1) - point.position(0);
        let mut tangent_u = [false; 2];
        let mut tangent_v = [false; 2];
        let mut boundary_u = [false; 2];
        let mut boundary_v = [false; 2];
        let mut binormal_u = [Vector3::default(); 2];
        let mut binormal_v = [Vector3::default(); 2];
        for s in 0..2 {
            let uv = point.uv[s];
            let domain = self.domains[s];
            let other_normal = point.normal(1 - s);
            let (su, sv) = (point.du(s), point.dv(s));
            if crate::scalar::are_same(domain.min.x, uv.x) {
                boundary_u[s] = true;
                tangent_u[s] = other_normal.is_perpendicular_to(sv, angle_tolerance);
                binormal_u[s] = su;
            } else if crate::scalar::are_same(domain.max.x, uv.x) {
                boundary_u[s] = true;
                tangent_u[s] = other_normal.is_perpendicular_to(sv, angle_tolerance);
                binormal_u[s] = -su;
            }
            if crate::scalar::are_same(domain.min.y, uv.y) {
                boundary_v[s] = true;
                tangent_v[s] = other_normal.is_perpendicular_to(su, angle_tolerance);
                binormal_v[s] = sv;
            } else if crate::scalar::are_same(domain.max.y, uv.y) {
                boundary_v[s] = true;
                tangent_v[s] = other_normal.is_perpendicular_to(su, angle_tolerance);
                binormal_v[s] = -sv;
            }
        }
        let on_tangent =
            (0..2).any(|s| (boundary_u[s] && tangent_u[s]) || (boundary_v[s] && tangent_v[s]));
        if !on_tangent {
            return Ok(false);
        }
        let scaled_zero = EFF_ZERO * (1.0 + point.curve_point.max_dimension());
        for s in 0..2 {
            if !tangent_u[s] && !tangent_v[s] {
                continue;
            }
            let (binormal, tangent) = if tangent_u[s] {
                (binormal_u[s], point.dv(s))
            } else {
                (binormal_v[s], point.du(s))
            };
            let stride = self.approx_tolerance / tangent.length() / 2.0 * 100.0;
            let mut side_degenerate = [false; 2];
            for side in 0..2 {
                let (mut zero, mut mono, mut div) = (true, true, true);
                if (side == 0 && !p_inside) || (side == 1 && !n_inside) {
                    side_degenerate[side] = true;
                    continue;
                }
                let mut curve_point = [Vector3::default(); 3];
                let mut gaps = [Vector3::default(); 3];
                let mut gap_length = [0.0; 3];
                let mut inside = [false; 3];
                let mut sample_distance = 0.0;
                curve_point[0] = point.position(s);
                gaps[0] = Vector3::new(0.0, 0.0, 0.0);
                gap_length[0] = 0.0;
                inside[0] = true;
                let mut delta = if side == 0 { stride } else { -stride };
                let mut ii = 1usize;
                loop {
                    let tgt = if ii > 1 { 2 } else { 1 };
                    inside[tgt] = false;
                    let uv0 = point.uv[s] + delta * point.uv_deriv[s];
                    let uv1 = point.uv[1 - s] + delta * point.uv_deriv[1 - s];
                    let mut other_pv = [Vector3::default(); 3];
                    let mut drop_uv = Vector2::new(0.0, 0.0);
                    if self.domains[s].contains_point(uv0, 10.0 * EFF_ZERO)
                        && self.domains[1 - s].contains_point(uv1, 10.0 * EFF_ZERO)
                    {
                        curve_point[tgt] =
                            crate::surface_eval::evaluate_point(self.surfaces[s], uv0)
                                .unwrap_or_default();
                        sample_distance = (point.position(1) - curve_point[tgt]).length();
                        let found = crate::surface_solve::local_point_solve(
                            self.surfaces[1 - s],
                            self.domains[1 - s],
                            SolverOperation::Normalize,
                            curve_point[tgt],
                            uv1,
                            None,
                        );
                        if let Ok(found) = found
                            && found.found
                        {
                            drop_uv = found.uv;
                            if let Ok(cells) = evaluate(
                                self.surfaces[1 - s],
                                drop_uv,
                                1,
                                1,
                                [true, true],
                                true,
                                true,
                                true,
                            ) {
                                // `[0][0]`, `[1][0]` (Su) and `[0][1]` (Sv).
                                other_pv = [cells[0], cells[2], cells[1]];
                            }
                            let shift = if s == 0 { -gap } else { gap };
                            gaps[tgt] = (other_pv[0] - curve_point[tgt]) + shift;
                            gap_length[tgt] = gaps[tgt].length();
                            inside[tgt] = true;
                        }
                    }
                    if !inside[tgt] {
                        let previous = (tgt + 3 - 1) % 3;
                        side_degenerate[side] = ii < 3
                            || (mono
                                && div
                                && binormal.dot(gaps[previous]) > 0.0
                                && sample_distance < self.approx_tolerance * 10.0);
                        break;
                    } else if ii >= 2 {
                        let tgt1 = (tgt + 3 - 1) % 3;
                        let tgt2 = (tgt + 3 - 2) % 3;
                        if gap_length[tgt] > scaled_zero {
                            zero = false;
                        }
                        let bxt = binormal.cross(tangent);
                        let sc = other_pv[0] - curve_point[tgt];
                        let (su, sv) = (other_pv[1], other_pv[2]);
                        let mut a = crate::matrix::Matrix::new(3, 3);
                        a[0][0] = 2.0 * su.dot(su);
                        a[0][1] = 2.0 * su.dot(sv);
                        a[0][2] = -bxt.dot(su);
                        a[1][0] = 2.0 * sv.dot(su);
                        a[1][1] = 2.0 * sv.dot(sv);
                        a[1][2] = -bxt.dot(sv);
                        a[2][0] = bxt.dot(su);
                        a[2][1] = bxt.dot(sv);
                        a[2][2] = 0.0;
                        let b = [-2.0 * sc.dot(su), -2.0 * sc.dot(sv), -bxt.dot(sc)];
                        let mut x = Vec::new();
                        let open = match a.solve_linear_system(&b, &mut x) {
                            Err(_) => true,
                            Ok(_) => {
                                let puv = Vector2::new(drop_uv.x + x[0], drop_uv.y + x[1]);
                                let pc = other_pv[0] + x[0] * su + x[1] * sv - curve_point[tgt];
                                !self.domains[1 - s].contains_point(puv, EFF_ZERO)
                                    || pc.dot(binormal) < -EFF_ZERO
                            }
                        };
                        if !zero {
                            if gaps[tgt].dot(gaps[tgt1]) <= 0.0 {
                                mono = false;
                            }
                            if gap_length[tgt] == 0.0
                                || gap_length[tgt]
                                    < gap_length[tgt1]
                                        + 1.0001 * (gap_length[tgt1] - gap_length[tgt2])
                            {
                                div = false;
                            }
                        }
                        if !zero && mono && div && open {
                            side_degenerate[side] = true;
                            break;
                        } else if zero && mono && div {
                            curve_point[1] = curve_point[2];
                            gaps[1] = gaps[2];
                            gap_length[1] = gap_length[2];
                            inside[1] = inside[2];
                        } else {
                            side_degenerate[side] = false;
                            break;
                        }
                    }
                    ii += 1;
                    delta *= 5.0;
                }
            }
            if side_degenerate[0] && side_degenerate[1] {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

/// A side curve of a surface (`SmSurfaceCache::GetIsoBoundaryCurves`): a
/// B-spline iso-curve trimmed to the domain, or a plane's line.
enum SideCurve {
    Spline(crate::spline::BSplineCurve),
    Line(crate::line::Line),
}

impl Intersector<'_> {
    /// `GetIsoBoundaryCurves` through `CreateIsoParametricCurve` with the
    /// domain: the sides at u min, v min, u max and v max. A plane's sides
    /// are `SmLine` segments between its corners (`EvaluatePointFast`).
    fn side_curves(&self, s: usize) -> Result<Vec<SideCurve>, SolverRefusal> {
        let surface = self.surfaces[s];
        let domain = self.domains[s];
        let mut curves = Vec::with_capacity(4);
        for (axis, value) in [
            (0, domain.min.x),
            (1, domain.min.y),
            (0, domain.max.x),
            (1, domain.max.y),
        ] {
            let interval = if axis == 0 {
                domain.v_interval()
            } else {
                domain.u_interval()
            };
            if let Some(frame) = surface.plane_frame() {
                let (mut start, mut end) = (domain.min, domain.max);
                if axis == 0 {
                    start.x = value;
                    end.x = value;
                } else {
                    start.y = value;
                    end.y = value;
                }
                let line = crate::line::Line::segment(
                    frame.point_fast(start),
                    frame.point_fast(end),
                    false,
                    Some(interval),
                )?;
                curves.push(SideCurve::Line(line));
                continue;
            }
            let iso = surface.iso_curve(axis, value)?;
            // SMLib ignores the trim's status.
            let curve = match iso.trim(interval, ZONE_TOL_3D) {
                Ok(trimmed) => trimmed.curve,
                Err(crate::curve_edit::TrimRefusal::Failed(_)) => iso,
                Err(crate::curve_edit::TrimRefusal::ShortSpan) => {
                    return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
                }
            };
            curves.push(SideCurve::Spline(curve));
        }
        Ok(curves)
    }
}

/// `SmSurface::EvaluateGeometric` at a point (from the left on both sides):
/// the principal curvatures and directions.
fn principal_at(
    surface: &crate::spline::BSplineSurface,
    uv: Vector2,
) -> Result<(f64, f64, Vector3, Vector3), SolverRefusal> {
    let mut jet = crate::curve_surface::surface_jet(surface, uv, [true, true])?;
    if jet[0][1].cross(jet[1][0]).length_squared() < EFF_ZERO {
        jet[2][2] = evaluate_normal(surface, uv, [true, true]);
    }
    let geometry = crate::curve_surface::surface_geometry(&mut jet)?;
    Ok((
        geometry.k1,
        geometry.k2,
        geometry.direction1,
        geometry.direction2,
    ))
}

/// `sm_ComputeCurveType`: crossing, or where the normals at the parameter
/// curves' middles are parallel, tangent or (matching principal curvatures)
/// coincident.
fn compute_curve_type(
    surfaces: [&crate::spline::BSplineSurface; 2],
    uv_curves: [&BSplineCurve; 2],
) -> Result<crate::surface_curves::CurveKind, SolverRefusal> {
    use crate::surface_curves::CurveKind;
    let mut uv = [Vector2::new(0.0, 0.0); 2];
    for i in 0..2 {
        let interval = uv_curves[i].natural_interval();
        let p = uv_curves[i].evaluate_point(interval.evaluate(0.5))?;
        uv[i] = Vector2::new(p.x, p.y);
    }
    let n1 = evaluate_normal(surfaces[0], uv[0], [true, true]);
    let n2 = evaluate_normal(surfaces[1], uv[1], [true, true]);
    let dot = n1.dot(n2);
    let sign = if dot < 0.0 { -1.0 } else { 1.0 };
    let limit = 0.0001;
    if dot.abs() <= 1.0 - limit {
        return Ok(CurveKind::Crossing);
    }
    // SMLib ignores a failed evaluation (`SE`), reading zeros.
    let zero = (0.0, 0.0, Vector3::default(), Vector3::default());
    let (k11, k12, v11, _) = principal_at(surfaces[0], uv[0]).unwrap_or(zero);
    let (k21, k22, v21, v22) = principal_at(surfaces[1], uv[1]).unwrap_or(zero);
    let mut kind = CurveKind::Tangent;
    let scale = 1.0 + k11.abs();
    if (k11 - k12).abs() < limit * scale {
        if (k21 - k22).abs() > limit * scale {
            return Ok(kind);
        }
        if (k11 - sign * k21).abs() > limit * scale {
            return Ok(kind);
        }
        return Ok(CurveKind::Coincident);
    }
    if (k21 - k22).abs() < limit * scale {
        return Ok(kind);
    }
    if v11.dot(v21).abs() > 1.0 - limit {
        if (k11 - sign * k21).abs() > limit * scale {
            return Ok(kind);
        }
        let scale2 = 1.0 + k12.abs();
        if (k12 - sign * k22).abs() > limit * scale2 {
            return Ok(kind);
        }
        kind = CurveKind::Coincident;
    }
    if v11.dot(v22).abs() > 1.0 - limit {
        if (k11 - sign * k22).abs() > limit * scale {
            return Ok(kind);
        }
        let scale2 = 1.0 + k12.abs();
        if (k12 - sign * k21).abs() > limit * scale2 {
            return Ok(kind);
        }
        kind = CurveKind::Coincident;
    }
    Ok(kind)
}

/// A parameter curve from `DropAndTrimCurve` as a B-spline (a plane's
/// `SmLine` through its NURBS form).
fn uv_spline(curve: crate::surface_drop_curve::UvCurve) -> Result<BSplineCurve, SolverRefusal> {
    match curve {
        crate::surface_drop_curve::UvCurve::Spline(spline) => Ok(spline),
        crate::surface_drop_curve::UvCurve::Line(line) => Ok(line.make_nurb()?),
    }
}

/// `SmBSplineCurve::Trim` as SMLib calls it here (status ignored).
fn trim_to(curve: &BSplineCurve, interval: Extent1) -> Result<BSplineCurve, SolverRefusal> {
    match curve.trim(interval, ZONE_TOL_3D) {
        Ok(trimmed) => Ok(trimmed.curve),
        Err(crate::curve_edit::TrimRefusal::Failed(_)) => Ok(curve.clone()),
        Err(crate::curve_edit::TrimRefusal::ShortSpan) => {
            Err(SolverRefusal::NotPorted(SolverOperation::Intersect))
        }
    }
}

impl Tracer {
    /// The long-range branch of `FindBoundaryStartPoints`: a side of
    /// surface `s` lying in the other surface over a range becomes an
    /// intersection curve (the side trimmed to the range, dropped onto both
    /// surfaces with `DropAndTrimCurve`), unless its two ends' gaps differ a
    /// hundredfold, it lies on an earlier curve, or it drops too far.
    #[allow(clippy::too_many_arguments)]
    fn side_range_curve(
        &mut self,
        intersector: &Intersector<'_>,
        side: &SideCurve,
        mut solution: crate::solution::Solution,
        curves: &mut Vec<IntersectionCurve>,
        edge_crossings: &mut Vec<TsectPoint>,
    ) -> Result<(), SolverRefusal> {
        let tolerance = intersector.approx_tolerance;
        if solution.end.value < EFF_ZERO {
            solution.end.value = EFF_ZERO;
        }
        if solution.start.value < EFF_ZERO {
            solution.start.value = EFF_ZERO;
        }
        let ratio = solution.start.value / solution.end.value;
        if ratio > 100.0 || ratio < 0.01 {
            return Ok(());
        }
        let middle = (solution.start.parameters[0] + solution.end.parameters[0]) / 2.0;
        let test = match side {
            SideCurve::Spline(curve) => curve.evaluate_point(middle)?,
            SideCurve::Line(line) => line.evaluate_point(middle),
        };
        if intersector.is_point_on_curve(test, curves, None)? {
            return Ok(());
        }
        let copy = match side {
            SideCurve::Spline(curve) => curve.clone(),
            SideCurve::Line(line) => line.make_nurb()?,
        };
        let mut trim = Extent1 {
            min: solution.start.parameters[0],
            max: solution.end.parameters[0],
        };
        let mut curve = trim_to(&copy, trim)?;
        let mut uv: Vec<BSplineCurve> = Vec::with_capacity(2);
        let mut distances = [0.0; 2];
        for target in 0..2 {
            let dropped = crate::surface_drop_curve::drop_and_trim_curve(
                intersector.surfaces[target],
                intersector.domains[target],
                &curve,
                curve.natural_interval(),
                tolerance,
                true,
                true,
            )?;
            if dropped.curves.len() != 1 {
                // Two curves (a seam) cannot occur without seam starts; none
                // drops the curve.
                return Ok(());
            }
            distances[target] = dropped.max_drop;
            let uv_curve = uv_spline(dropped.curves.into_iter().next().ok_or(Fault::Degenerate)?)?;
            let uv_interval = uv_curve.natural_interval();
            if !trim.is_contained_by(uv_interval, 0.0) {
                trim = uv_interval;
                curve = trim_to(&curve, trim)?;
                if target == 1 {
                    uv[0] = trim_to(&uv[0], trim)?;
                }
            }
            uv.push(uv_curve);
        }
        if distances[0] > tolerance || distances[1] > tolerance {
            return Ok(());
        }
        let kind = compute_curve_type(intersector.surfaces, [&uv[0], &uv[1]])?;
        let second = uv.pop().ok_or(Fault::Degenerate)?;
        let first = uv.pop().ok_or(Fault::Degenerate)?;
        curves.push(IntersectionCurve {
            curve,
            uv_curves: [first, second],
            kind,
            deviation: crate::scalar::max(distances[0], distances[1]),
        });
        let mut kept = Vec::with_capacity(self.start_points.len());
        for point in std::mem::take(&mut self.start_points) {
            if !intersector.is_point_on_curve(point.curve_point, curves, None)? {
                kept.push(point);
            }
        }
        self.start_points = kept;
        let mut kept = Vec::with_capacity(edge_crossings.len());
        for point in std::mem::take(edge_crossings) {
            if !intersector.is_point_on_curve(point.curve_point, curves, None)? {
                kept.push(point);
            }
        }
        *edge_crossings = kept;
        Ok(())
    }
}

impl Tracer {
    /// `FindBoundaryStartPoints` without trimming faces: meets each side
    /// curve of each surface with the other surface inside the overlap box;
    /// crossings away from `curves` become start points (those on the other
    /// surface's side last), and touches become touching curves appended to
    /// `curves`. Range answers (a side lying in the other surface) and a
    /// plane's side lines inside the box refuse with
    /// [`SolverRefusal::NotPorted`].
    #[allow(clippy::too_many_lines)]
    pub fn find_boundary_start_points(
        &mut self,
        intersector: &Intersector<'_>,
        overlap: crate::extent::Extent3,
        curves: &mut Vec<IntersectionCurve>,
    ) -> Result<(), SolverRefusal> {
        let tolerance = intersector.approx_tolerance;
        let mut edge_crossings: Vec<TsectPoint> = Vec::new();
        let mut degenerate: Vec<(TsectPoint, Vector3, usize)> = Vec::new();
        for s in 0..2 {
            if !intersector.use_surface_edges[s] {
                continue;
            }
            let sides = intersector.side_curves(s)?;
            for (index, side) in sides.iter().enumerate() {
                let along = index % 2;
                let norm = if index < 2 { 0.0 } else { 1.0 };
                let domain = intersector.domains[s];
                let interval = if along == 0 {
                    domain.v_interval()
                } else {
                    domain.u_interval()
                };
                let other = 1 - s;
                let solutions = match side {
                    SideCurve::Spline(curve) => {
                        if crate::curve_surface_global::degenerate_on(curve, interval) {
                            continue;
                        }
                        let mut bbox = curve.bounding_box(interval)?;
                        bbox.expand_absolute(tolerance);
                        if bbox.are_disjoint(overlap, 0.0) {
                            continue;
                        }
                        crate::curve_surface_global::global_curve_intersect(
                            intersector.surfaces[other],
                            intersector.domains[other],
                            curve,
                            interval,
                            tolerance,
                            None,
                        )?
                    }
                    SideCurve::Line(line) => {
                        if line.is_degenerate(EFF_ZERO, Some(interval)) {
                            continue;
                        }
                        let mut bbox = line.bounding_boxes(interval).0;
                        bbox.expand_absolute(tolerance);
                        if bbox.are_disjoint(overlap, 0.0) {
                            continue;
                        }
                        // `GlobalCurveIntersect` takes an `SmLine` straight
                        // to `GlobalLineIntersect`.
                        crate::curve_surface_global::line_intersect(
                            intersector.surfaces[other],
                            intersector.domains[other],
                            line.evaluate(interval.min, 1)?,
                            interval,
                            tolerance,
                        )?
                    }
                };
                for found in solutions.as_slice() {
                    let mut solution = *found;
                    if solution.kind == crate::solution::SolutionType::RangeOfValues {
                        let size = solution.end.parameters[0] - solution.start.parameters[0];
                        if size > interval.length() / 100.0 {
                            self.side_range_curve(
                                intersector,
                                side,
                                solution,
                                curves,
                                &mut edge_crossings,
                            )?;
                            continue;
                        }
                        // A short range: its tighter end, moved onto the
                        // nearest normal meeting nearby.
                        let SideCurve::Spline(curve) = side else {
                            return Err(SolverRefusal::NotPorted(SolverOperation::Intersect));
                        };
                        let t0 = solution.start.parameters[0];
                        let uv0 = Vector2::new(
                            solution.start.parameters[1],
                            solution.start.parameters[2],
                        );
                        let mut near =
                            Extent1::new(solution.start.parameters[0], solution.end.parameters[0]);
                        near.expand_absolute(interval.length() * EFF_ZERO_SQRT);
                        let near = near.intersect(curve.natural_interval()).unwrap_or_default();
                        if solution.start.value > solution.end.value
                            && !intersector.domains[other]
                                .point_on_boundary(uv0, EFF_ZERO)
                                .0
                        {
                            solution.start = solution.end;
                        }
                        if let Some(better) = crate::curve_surface::local_curve_solve(
                            intersector.surfaces[other],
                            intersector.domains[other],
                            curve,
                            near,
                            SolverOperation::Normalize,
                            (t0, uv0),
                        )? && better.distance < solution.start.value
                        {
                            solution.kind = crate::solution::SolutionType::SingleValue;
                            solution.start.parameters[0] = better.t;
                            solution.start.parameters[1] = better.uv.x;
                            solution.start.parameters[2] = better.uv.y;
                            solution.start.value = better.distance;
                        }
                    }
                    let t = solution.start.parameters[0];
                    let mut uvs = [Vector2::new(0.0, 0.0); 2];
                    uvs[other] =
                        Vector2::new(solution.start.parameters[1], solution.start.parameters[2]);
                    let side_value = domain.evaluate(norm, norm);
                    uvs[s] = if along == 0 {
                        Vector2::new(side_value.x, t)
                    } else {
                        Vector2::new(t, side_value.y)
                    };
                    if singular(intersector, s, uvs[s]) || singular(intersector, other, uvs[other])
                    {
                        continue;
                    }
                    let mut point = TsectPoint::default();
                    match intersector.compute_point_values(
                        uvs,
                        &mut point,
                        None,
                        self.curve_points.last().map(|p| p.curve_tangent),
                        None,
                        self.trace_direction,
                    ) {
                        Ok(()) => {}
                        Err(SolverRefusal::Failed(_)) => continue,
                        Err(other) => return Err(other),
                    }
                    let pv = match side {
                        SideCurve::Spline(curve) => curve.evaluate(t, 1, true, true)?,
                        SideCurve::Line(line) => line.evaluate(t, 1)?,
                    };
                    let degenerate_point = intersector.classify_boundary_start(&point)?;
                    if intersector.is_point_on_curve(point.curve_point, curves, None)? {
                        continue;
                    }
                    if degenerate_point {
                        degenerate.push((point, pv[0], other));
                        continue;
                    }
                    point.deviation = solution.start.value;
                    if intersector.domains[other]
                        .point_on_boundary(uvs[other], EFF_ZERO * 100.0)
                        .0
                    {
                        edge_crossings.push(point);
                    } else {
                        self.add_start_point(intersector, point)?;
                    }
                }
            }
        }
        for point in edge_crossings {
            self.add_start_point(intersector, point)?;
        }
        for (point, curve_point, other) in degenerate {
            let middle = (curve_point + point.position(other)).divided(2.0);
            if intersector.is_point_on_curve(middle, curves, None)? {
                continue;
            }
            let on_start = self.start_points.iter().any(|start| {
                start.position(0).distance_squared(point.position(0)) < intersector.approx_tolerance
                    && start.position(1).distance_squared(point.position(1))
                        < intersector.approx_tolerance
            });
            if on_start {
                continue;
            }
            curves.push(IntersectionCurve::touching(
                middle,
                point.uv,
                point.position(0).distance(point.position(1)),
            )?);
        }
        Ok(())
    }
}
