// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: NLib stacks and in-place edits become owned values,
// statuses become `Fault`, and every surface is checked by
// `BSplineSurface::new`. Translated directly from SMLib
// SmBSplineSurfaceNLib.cpp (`CreateRuledSurface`, `Reparameterize`),
// SmBSplineSurface.cpp (the NLib-surface constructor), NL_SrfCommon.cpp
// (`N_CreateRuledSrf`), NL_CrvTool.cpp (`N_CrvsMakeCompatible`,
// `N_CrvsMakeCompatibleKnotTol`), NL_BasisBasic.cpp
// (`N_GetCompatibleKnotArray`) and NL_SrfBasic.cpp
// (`N_SrfReparamToInterval`). rust-brep has no port of these as a unit.
//! Ruled surfaces, surfaces of revolution and parameter-range edits.
//!
//! Every surface is an unproved candidate. A surface built from NLib data
//! runs `sm_FixupKnotVector` on both knot vectors with `SM_EFF_ZERO_PARAM`,
//! and so does this port. Where SMLib returns a surface the checked storage
//! cannot hold (a parameter range reduced to a single value), the port
//! refuses with [`Fault::InvalidInput`]; the gate records it.

#![allow(clippy::needless_range_loop)]

use crate::{
    Fault,
    curve::polygon,
    curve_arc::CircleParameterization,
    curve_build::fixup_knot_vector,
    curve_compat::make_compatible,
    extent::Extent2,
    scalar::EFF_ZERO,
    spline::{BSplineCurve, BSplineSurface, ControlPoint, NO_WEIGHT, NO_Z, unique_knots},
    vector::Vector3,
};

/// `SM_EFF_ZERO_PARAM`: the knot tolerance of SMLib's surface constructor.
const KNOT_TOLERANCE: f64 = 1.2e-8;

/// `SmSurfParamType`: which parameter direction of a surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceParameter {
    /// `SM_SP_U`.
    U,
    /// `SM_SP_V`.
    V,
}

/// `NL_WMIN` and `NL_WMAX`: the weights `N_Iges128Srf` accepts.
const MIN_WEIGHT: f64 = 1.0e-3;
const MAX_WEIGHT: f64 = 100.0;

/// A knot vector as `CreateCanonical` rebuilds it from unique knots and
/// multiplicities.
fn expanded_knots(curve: &BSplineCurve) -> Vec<f64> {
    unique_knots(curve.knots(), None)
        .into_iter()
        .flat_map(|(knot, multiplicity)| std::iter::repeat_n(knot, multiplicity))
        .collect()
}

/// `SmBSplineCurve::GetControlPolygon`: Euclidean points (by the weight's
/// reciprocal) and the weights, empty for a polynomial curve.
fn control_polygon(curve: &BSplineCurve) -> (Vec<Vector3>, Vec<f64>) {
    let mut points = Vec::new();
    let mut weights = Vec::new();
    let mut rational = false;
    for point in curve.control_points() {
        let [x, y, z, w] = point.coordinates;
        let z = if z == NO_Z { 0.0 } else { z };
        if w == NO_WEIGHT {
            points.push(Vector3::new(x, y, z));
            weights.push(0.0);
        } else {
            rational = true;
            let reciprocal = 1.0 / w;
            points.push(Vector3::new(x * reciprocal, y * reciprocal, z * reciprocal));
            weights.push(w);
        }
    }
    if !rational {
        weights.clear();
    }
    (points, weights)
}

/// The circle of `radius` about the axis through `center`, `x_axis` outward,
/// as `CreateSurfOfRevolution` builds it.
fn sweep_circle(
    axis: Vector3,
    center: Vector3,
    mut x_axis: Vector3,
    radius: f64,
    angle: f64,
    unitize_y: bool,
) -> Result<BSplineCurve, Fault> {
    x_axis /= radius;
    let mut y_axis = axis.cross(x_axis);
    if unitize_y {
        y_axis.unitize()?;
    }
    BSplineCurve::circle_segment(
        3,
        center,
        x_axis,
        y_axis,
        radius,
        0.0,
        angle,
        CircleParameterization::Quadratic,
    )
}

impl BSplineSurface {
    /// `SmBSplineSurface::CreateSurfOfRevolution` for a B-spline generator
    /// with clamped knots, turned `angle_degrees` about the axis through
    /// `origin` along `axis`, with SMLib's quadratic circular sweep. With
    /// `swap_uv` the generator runs along U and the circle along V; otherwise
    /// the circle is U. A generator point on the axis becomes a pole.
    ///
    /// SMLib also takes a sweep circle to copy its parameterisation and
    /// clamps and corrects the knots of unclamped data through
    /// `N_IgesClampSrf` and `N_IgesCorrectSrfKnots`; the port carries none of
    /// that and refuses such inputs with [`Fault::Degenerate`].
    pub fn surface_of_revolution(
        generator: &BSplineCurve,
        origin: Vector3,
        axis: Vector3,
        angle_degrees: f64,
        swap_uv: bool,
    ) -> Result<Self, Fault> {
        let mut axis = axis;
        axis.unitize()?;
        let curve_interval = generator.natural_interval();
        let curve_degree = generator.degree();
        let (curve_points, curve_weights) = control_polygon(generator);
        let curve_rational = !curve_weights.is_empty();
        let curve_knots = expanded_knots(generator);
        let on_axis = |point: Vector3| {
            let scaled_zero = EFF_ZERO * (1.0 + point.max_dimension());
            let along = (point - origin).dot(axis);
            let center = origin + along * axis;
            let offset = point - center;
            (center, offset, offset.length(), scaled_zero)
        };
        // The first generator point off the axis gives the circle's shape.
        let mut sweep = None;
        for &point in &curve_points {
            let (center, offset, radius, scaled_zero) = on_axis(point);
            if radius < scaled_zero {
                continue;
            }
            sweep = Some(sweep_circle(
                axis,
                center,
                offset,
                radius,
                angle_degrees,
                false,
            )?);
            break;
        }
        let sweep = sweep.ok_or(Fault::Degenerate)?;
        let (circle_points, circle_weights) = control_polygon(&sweep);
        let circle_count = circle_points.len();
        let circle_degree = sweep.degree();
        let circle_knots = expanded_knots(&sweep);
        let circle_interval = sweep.natural_interval();
        let mut points = Vec::new();
        let mut weights = Vec::new();
        for (j, &point) in curve_points.iter().enumerate() {
            for weight in &circle_weights {
                let mut weight = *weight;
                if curve_rational {
                    weight *= curve_weights[j];
                }
                weights.push(weight);
            }
            let (center, offset, radius, scaled_zero) = on_axis(point);
            if radius < scaled_zero {
                points.extend(std::iter::repeat_n(point, circle_count));
                continue;
            }
            let circle = sweep_circle(axis, center, offset, radius, angle_degrees, true)?;
            let polygon = polygon(&circle, circle.natural_interval())?;
            points.extend(polygon.into_iter().take(circle_count));
        }
        // The net is curve-major; without `swap_uv` the circle comes first.
        let rows = curve_points.len();
        let (degrees, knot_lists) = if swap_uv {
            ([curve_degree, circle_degree], [curve_knots, circle_knots])
        } else {
            let mut transposed_points = Vec::with_capacity(points.len());
            let mut transposed_weights = Vec::with_capacity(weights.len());
            for k in 0..circle_count {
                for j in 0..rows {
                    transposed_points.push(points[j * circle_count + k]);
                    transposed_weights.push(weights[j * circle_count + k]);
                }
            }
            points = transposed_points;
            weights = transposed_weights;
            ([circle_degree, curve_degree], [circle_knots, curve_knots])
        };
        let domain = if swap_uv {
            [curve_interval, circle_interval]
        } else {
            [circle_interval, curve_interval]
        };
        canonical_surface(degrees, knot_lists, domain, points, weights)
    }
}

/// `CreateCanonical` and `N_Iges128Srf` for data that needs no clamping.
fn canonical_surface(
    degrees: [usize; 2],
    mut knots: [Vec<f64>; 2],
    domain: [crate::extent::Extent1; 2],
    points: Vec<Vector3>,
    weights: Vec<f64>,
) -> Result<BSplineSurface, Fault> {
    if degrees.iter().any(|&d| !(1..=32).contains(&d)) {
        return Err(Fault::InvalidInput);
    }
    if !weights.is_empty() && weights.len() != points.len() {
        return Err(Fault::InvalidInput);
    }
    let counts = [0, 1].map(|axis| knots[axis].len().wrapping_sub(degrees[axis] + 1));
    if knots.iter().zip(degrees).any(|(k, d)| k.len() < d + 2)
        || counts[0].checked_mul(counts[1]) != Some(points.len())
    {
        return Err(Fault::InvalidInput);
    }
    let non_unit = weights.iter().any(|w| (w - 1.0).abs() > EFF_ZERO);
    let weights = if non_unit {
        weights
    } else {
        vec![1.0; points.len()]
    };
    if weights
        .iter()
        .any(|&w| !(MIN_WEIGHT..=MAX_WEIGHT).contains(&w))
    {
        return Err(Fault::Degenerate);
    }
    let rational = weights.iter().any(|&w| w != weights[0]);
    for axis in 0..2 {
        let (u, p) = (&knots[axis], degrees[axis]);
        let m = u.len() - 1;
        let (low, high) = (domain[axis].min, domain[axis].max);
        // `N_IgesValidate`.
        if u[p] > low || low >= high || u[m - p] < high || u.windows(2).any(|w| w[0] > w[1]) {
            return Err(Fault::Degenerate);
        }
        // `N_IgesClampSrf` and `N_IgesCorrectSrfKnots` change nothing only
        // for knots clamped at exactly the domain with no over-repeated knot.
        let clamped = u[..=p].iter().all(|&k| k == low)
            && u[m - p..].iter().all(|&k| k == high)
            && u[p + 1] != u[p]
            && u[m - p - 1] != u[m - p]
            && unique_knots(u, None).iter().all(|&(_, mult)| mult <= p + 1);
        if !clamped {
            return Err(Fault::Degenerate);
        }
    }
    let controls = points
        .iter()
        .zip(&weights)
        .map(|(point, &w)| ControlPoint {
            coordinates: if rational {
                [w * point.x, w * point.y, w * point.z, w]
            } else {
                [point.x, point.y, point.z, NO_WEIGHT]
            },
        })
        .collect();
    for axis in &mut knots {
        fixup_knot_vector(axis, KNOT_TOLERANCE);
    }
    BSplineSurface::new(degrees, knots, counts, controls)
}

impl BSplineSurface {
    /// `SmBSplineSurface::CreateRuledSurface`: the surface of straight lines
    /// joining `first` and `second` after they are made compatible. With
    /// `linear` U the surface is linear in U and the curves vary in V; with V
    /// the reverse. The control points always carry z: a planar curve's z is
    /// zero.
    pub fn ruled(
        first: &BSplineCurve,
        second: &BSplineCurve,
        linear: SurfaceParameter,
    ) -> Result<Self, Fault> {
        let compatible = make_compatible(&[first.clone(), second.clone()], 0.0)?;
        let (a, b) = (&compatible[0], &compatible[1]);
        let planar = |curve: &BSplineCurve| curve.control_points()[0].coordinates[2] == NO_Z;
        let with_z = |point: &ControlPoint, planar: bool| {
            let mut point = *point;
            if planar {
                point.coordinates[2] = 0.0;
            }
            point
        };
        let (a_planar, b_planar) = (planar(a), planar(b));
        let count = a.control_points().len();
        let degree = a.degree();
        let mut curve_knots = a.knots().to_vec();
        let mut line_knots = vec![0.0, 0.0, 1.0, 1.0];
        let mut points = Vec::with_capacity(2 * count);
        let (degrees, counts) = match linear {
            SurfaceParameter::U => {
                points.extend(a.control_points().iter().map(|p| with_z(p, a_planar)));
                points.extend(b.control_points().iter().map(|p| with_z(p, b_planar)));
                ([1, degree], [2, count])
            }
            SurfaceParameter::V => {
                for i in 0..count {
                    points.push(with_z(&a.control_points()[i], a_planar));
                    points.push(with_z(&b.control_points()[i], b_planar));
                }
                ([degree, 1], [count, 2])
            }
        };
        fixup_knot_vector(&mut line_knots, KNOT_TOLERANCE);
        fixup_knot_vector(&mut curve_knots, KNOT_TOLERANCE);
        let knots = match linear {
            SurfaceParameter::U => [line_knots, curve_knots],
            SurfaceParameter::V => [curve_knots, line_knots],
        };
        BSplineSurface::new(degrees, knots, counts, points)
    }

    /// `SmBSplineSurface::Reparameterize` through `N_SrfReparamToInterval`:
    /// the same surface on the parameter range `domain`.
    pub fn reparameterized(&self, domain: Extent2) -> Result<Self, Fault> {
        if self.plane_frame().is_some() {
            return Err(Fault::InvalidInput);
        }
        let mut knots = [self.knots()[0].to_vec(), self.knots()[1].to_vec()];
        let ranges = [(domain.min.x, domain.max.x), (domain.min.y, domain.max.y)];
        for axis in 0..2 {
            let p = self.degrees()[axis];
            let u = &mut knots[axis];
            let m = u.len() - 1;
            let (a, b) = ranges[axis];
            if a != u[0] || b != u[m] {
                let u0 = u[0];
                let factor = (b - a) / (u[m] - u[0]);
                for i in 0..=p {
                    u[i] = a;
                }
                for i in p + 1..=m - p - 1 {
                    u[i] = factor * (u[i] - u0) + a;
                }
                for i in m - p..=m {
                    u[i] = b;
                }
            }
        }
        BSplineSurface::new(
            self.degrees(),
            knots,
            self.counts(),
            self.control_points().to_vec(),
        )
    }
}
