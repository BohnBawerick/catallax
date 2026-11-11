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
//! Ruled surfaces and parameter-range edits.
//!
//! Every surface is an unproved candidate. A surface built from NLib data
//! runs `sm_FixupKnotVector` on both knot vectors with `SM_EFF_ZERO_PARAM`,
//! and so does this port. Where SMLib returns a surface the checked storage
//! cannot hold (a parameter range reduced to a single value), the port
//! refuses with [`Fault::InvalidInput`]; the gate records it.

#![allow(clippy::needless_range_loop)]

use crate::{
    Fault,
    curve_build::fixup_knot_vector,
    curve_compat::make_compatible,
    extent::Extent2,
    spline::{BSplineCurve, BSplineSurface, ControlPoint, NO_Z},
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
