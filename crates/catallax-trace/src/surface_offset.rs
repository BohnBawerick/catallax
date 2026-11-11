// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: the offset surface borrows a checked B-spline
// surface instead of SMLib's polymorphic base, statuses become `Fault`, and
// nested offsets are not carried; the extended surface that SMLib keeps in a
// member and replaces from `const` methods is a shared cell. Translated
// directly from SMLib SmOffsetSurface.cpp (the constructor's flags,
// `EvaluatePoint`, `EvaluateNormal`, `Evaluate`, `EvaluateSimple`,
// `EvaluatePlainOrSimple`, `sm_ExtendSurface`, `sm_FindNonsingularStepoff`,
// `GetNaturalUVDomain`); rust-brep has no port of these.
//! The offset of a B-spline surface by a signed distance along its normal.
//!
//! `Offset(u, v) = Base(u, v) + distance * N(u, v)`. Derivatives come from
//! the base surface's derivatives one order higher through the derivatives of
//! the unit normal, as SMLib computes them, quirks included: its second
//! mixed and third pure derivatives in `V` store the base derivative of the
//! wrong index, and the port stores the same.
//!
//! With extension switched on ([`OffsetSurface::allow_extension`]) an
//! evaluation outside the base surface's domain extends the base surface
//! (`sm_ExtendSurface`, in `surface_extend.rs`), wraps the parameter where the
//! base is periodic and stops it at a singular side, and keeps the extension
//! for later evaluations, as SMLib does; a later evaluation outside that
//! extends it again. With it off the parameter is clamped to the domain. The
//! port carries no nested offsets.
//!
//! Everything returned is an unproved candidate.

#![allow(clippy::too_many_arguments)]

use crate::{
    Fault,
    curve_surface::get_singularities,
    extent::Extent2,
    scalar::{EFF_ZERO, EFF_ZERO_SQ, ZONE_TOL_3D},
    spline::BSplineSurface,
    surface_eval::{
        EvaluableSurface, SolverSurface, evaluate, evaluate_normal, evaluate_normal_of,
        evaluate_simple, first_derivatives,
    },
    surface_extend::{extend_for_evaluation, is_periodic},
    vector::{Vector2, Vector3},
};
use std::{cell::RefCell, rc::Rc};

/// SMLib's highest offset derivative.
const MAX_DERIVATIVES: usize = 3;

/// A surface offset by `distance` along the base surface's normal.
#[derive(Clone, Debug)]
pub struct OffsetSurface<'a> {
    base: &'a BSplineSurface,
    distance: f64,
    /// `m_bAllowExtension`.
    allow_extension: bool,
    /// `m_bIsPeriodicOnU` and `m_bIsPeriodicOnV`.
    periodic: [bool; 2],
    /// `m_lSingularities`.
    singularities: u32,
    /// `m_pExtendedSurface`, shared by the views of this surface at other
    /// distances ([`Self::with_distance`]), as SMLib's offset is one object.
    extended: Rc<RefCell<Option<Rc<BSplineSurface>>>>,
}

/// The surface an evaluation acts on: the base or its extension.
enum Acting<'s> {
    Base(&'s BSplineSurface),
    Extended(Rc<BSplineSurface>),
}

impl Acting<'_> {
    fn surface(&self) -> &BSplineSurface {
        match self {
            Acting::Base(surface) => surface,
            Acting::Extended(surface) => surface,
        }
    }
}

/// `sm_FindNonsingularStepoff`: a parameter step from `uv` toward the
/// domain's middle at which the base surface has a first derivative longer
/// than `tolerance` in the direction that degenerates.
fn find_nonsingular_stepoff(
    surface: &BSplineSurface,
    uv: Vector2,
    domain: Extent2,
    bu_length: f64,
    bv_length: f64,
    tolerance: f64,
) -> Result<Vector2, Fault> {
    let u_is_degenerate = bu_length < bv_length;
    let middle = Vector2::new(
        (domain.min.x + domain.max.x) / 2.0,
        (domain.min.y + domain.max.y) / 2.0,
    );
    let mut step_size = 0.001;
    while step_size < 0.2 {
        for along_v in [u_is_degenerate, !u_is_degenerate] {
            let step = if along_v {
                Vector2::new(0.0, (middle.y - uv.y) * step_size)
            } else {
                Vector2::new((middle.x - uv.x) * step_size, 0.0)
            };
            let cells = evaluate_simple(surface, uv + step, 1, 1, [true, true], true)?;
            // `sMat[2]` is the U derivative and `sMat[1]` the V derivative.
            let longer = if u_is_degenerate {
                cells[2].length() > tolerance
            } else {
                cells[1].length() > tolerance
            };
            if longer {
                return Ok(step);
            }
        }
        step_size *= 10.0;
    }
    Err(Fault::Degenerate)
}

impl<'a> OffsetSurface<'a> {
    /// `SmOffsetSurface::SmOffsetSurface`: `base` offset by `distance`, which
    /// is negative for an inset.
    pub fn new(base: &'a BSplineSurface, distance: f64) -> Self {
        let domain = base.natural_domain();
        Self {
            base,
            distance,
            allow_extension: false,
            periodic: [is_periodic(base, domain, 0), is_periodic(base, domain, 1)],
            singularities: get_singularities(base, ZONE_TOL_3D).0,
            extended: Rc::new(RefCell::new(None)),
        }
    }

    /// `SetAllowExtension`: whether an evaluation outside the base surface's
    /// domain extends it (true) or clamps the parameter (false, the default).
    #[must_use]
    pub fn allow_extension(mut self, allow: bool) -> Self {
        self.allow_extension = allow;
        self
    }

    /// The same surface at another offset distance, sharing the extension,
    /// as `SetOffsetDistance` leaves SMLib's one object.
    #[must_use]
    pub fn with_distance(&self, distance: f64) -> Self {
        Self {
            distance,
            ..self.clone()
        }
    }

    /// `m_bIsPeriodicOnU` (axis 0) or `m_bIsPeriodicOnV` (axis 1).
    pub fn is_periodic(&self, axis: usize) -> bool {
        self.periodic[axis]
    }

    /// `m_lSingularities`: the singular sides, as SMLib's `SM_SS_*` flags.
    pub fn singularities(&self) -> u32 {
        self.singularities
    }

    /// `m_pExtendedSurface`: the extension SMLib is holding, if any.
    pub fn extended_surface(&self) -> Option<Rc<BSplineSurface>> {
        self.extended.borrow().clone()
    }

    /// The surface an evaluation at `uv` acts on, with `uv` moved into its
    /// domain: the base, or the extension when `uv` is outside the base and
    /// one is held; when `uv` is still outside, extended (or wrapped, or
    /// stopped at a singular side) if allowed and clamped if not.
    fn acting(&self, uv: &mut Vector2) -> Result<Acting<'a>, Fault> {
        let mut surface = Acting::Base(self.base);
        let mut domain = self.base.natural_domain();
        if !domain.contains_point(*uv, 0.0)
            && let Some(extended) = self.extended.borrow().clone()
        {
            domain = extended.natural_domain();
            surface = Acting::Extended(extended);
        }
        if !domain.contains_point(*uv, 0.0) {
            if self.allow_extension {
                let (extension, wrapped) = extend_for_evaluation(
                    surface.surface(),
                    *uv,
                    self.periodic,
                    self.singularities,
                )?;
                *uv = wrapped;
                if let Some(extension) = extension {
                    let extension = Rc::new(extension);
                    *self.extended.borrow_mut() = Some(Rc::clone(&extension));
                    surface = Acting::Extended(extension);
                }
            } else {
                *uv = domain.clamp_point(*uv);
            }
        }
        Ok(surface)
    }

    /// `SmOffsetSurface::GetNaturalUVDomain`: the base surface's, or the
    /// extension's once one is held.
    pub fn natural_domain(&self) -> Extent2 {
        match self.extended.borrow().as_ref() {
            Some(extended) => extended.natural_domain(),
            None => self.base.natural_domain(),
        }
    }

    /// The signed offset distance.
    pub fn distance(&self) -> f64 {
        self.distance
    }

    /// `SmOffsetSurface::EvaluateNormal`: the base surface's unit normal
    /// inside the domain, the offset surface having the same normal. Outside
    /// it SMLib uses `SmSurface::EvaluateNormal` on the offset surface itself.
    pub fn evaluate_normal(&self, uv: Vector2, sides: [bool; 2]) -> Vector3 {
        if self.base.natural_domain().contains_point(uv, 0.0) {
            evaluate_normal(self.base, uv, sides)
        } else {
            evaluate_normal_of(self, uv, sides)
        }
    }

    /// `SmOffsetSurface::EvaluatePoint`.
    pub fn evaluate_point(&self, uv: Vector2) -> Result<Vector3, Fault> {
        let mut uv = uv;
        let acting = self.acting(&mut uv)?;
        let surface = acting.surface();
        if self.distance == 0.0 {
            return crate::surface_eval::evaluate_point(surface, uv);
        }
        let [point, su, sv] = evaluate_first(surface, uv)?;
        let mut normal = su.cross(sv);
        if normal.length_squared() < EFF_ZERO_SQ {
            normal = evaluate_normal(surface, uv, [true, true]);
        }
        normal.unitize()?;
        Ok(point + self.distance * normal)
    }

    /// `SmOffsetSurface::Evaluate`: the point and derivatives to
    /// `(hu, hv)` (at most three in total), with `non_zero_tangents`
    /// replacing vanishing first derivatives as SMLib does. The result is
    /// `(hu + 1) * (hv + 1)` cells, row `u` and column `v`; cells SMLib does
    /// not write keep SMLib's undefined marker.
    pub fn evaluate(
        &self,
        uv: Vector2,
        hu: usize,
        hv: usize,
        sides: [bool; 2],
        upper: bool,
        non_zero_tangents: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        let mut cells = self.plain_or_simple(!non_zero_tangents, uv, hu, hv, sides, upper)?;
        cells.truncate((hu + 1) * (hv + 1));
        Ok(cells)
    }

    /// `SmOffsetSurface::EvaluateSimple`: [`Self::evaluate`] without repair
    /// of vanishing first derivatives.
    pub fn evaluate_simple(
        &self,
        uv: Vector2,
        hu: usize,
        hv: usize,
        sides: [bool; 2],
        upper: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        let mut cells = self.plain_or_simple(true, uv, hu, hv, sides, upper)?;
        cells.truncate((hu + 1) * (hv + 1));
        Ok(cells)
    }

    /// `SmOffsetSurface::EvaluatePlainOrSimple`, returning SMLib's whole
    /// scratch matrix of 64 cells, which a non-square upper-half evaluation
    /// fills past the declared size.
    fn plain_or_simple(
        &self,
        simple: bool,
        uv_in: Vector2,
        hu: usize,
        hv: usize,
        sides: [bool; 2],
        upper: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        let offset = self.distance;
        // Each order is at most three whichever way the total is counted;
        // check before adding them.
        if hu > MAX_DERIVATIVES || hv > MAX_DERIVATIVES {
            return Err(Fault::Degenerate);
        }
        let highest_offset = if upper { hu.max(hv) } else { hu + hv };
        let highest_base = highest_offset + 1;
        if highest_offset > MAX_DERIVATIVES {
            return Err(Fault::Degenerate);
        }
        let mut uv = uv_in;
        let acting = self.acting(&mut uv)?;
        let surface = acting.surface();
        let domain = surface.natural_domain();
        if offset == 0.0 {
            let mut cells = if simple {
                evaluate_simple(surface, uv, hu, hv, sides, upper)?
            } else {
                evaluate(surface, uv, hu, hv, sides, upper, true, true)?
            };
            if upper && hu == hv {
                // SMLib's upper-half evaluation of a square matrix leaves the
                // other cells as the caller had them.
                for u in 0..=hu {
                    for v in 0..=hv {
                        if u + v > hu {
                            cells[u * (hv + 1) + v] = Vector3::default();
                        }
                    }
                }
            }
            return Ok(cells);
        }
        let stride = highest_base + 1;
        let base = evaluate_simple(surface, uv, highest_base, highest_base, [true, true], true)?;
        let cell = |u: usize, v: usize| base[u * stride + v];
        let b = cell(0, 0);
        let (mut bu, mut bv) = (cell(1, 0), cell(0, 1));
        let tolerance = 10.0 * EFF_ZERO * (1.0 + b.max_dimension());
        let (bu_length, bv_length) = (bu.length(), bv.length());
        // In upper-half mode SMLib writes cells past a non-square matrix;
        // the port keeps them in scratch and returns the declared cells.
        let mut out = vec![Vector3::default(); 64];
        let column = hv + 1;
        if !simple && (bu_length < tolerance || bv_length < tolerance) {
            // At a pole: extrapolate from two steps off it.
            let step =
                find_nonsingular_stepoff(surface, uv_in, domain, bu_length, bv_length, tolerance)?;
            let first = self.plain_or_simple(simple, uv_in + step, hu, hv, sides, upper)?;
            let second = self.plain_or_simple(simple, uv_in + 2.0 * step, hu, hv, sides, upper)?;
            let extrapolate = |index: usize| 2.0 * first[index] - second[index];
            let tangent = |index: usize| -> Result<Vector3, Fault> {
                let mut v1 = first[index];
                let mut v2 = second[index];
                let mut v0 = 2.0 * v1 - v2;
                if !simple && v0.length() < tolerance {
                    v1.unitize()?;
                    v2.unitize()?;
                    v0 = 2.0 * v1 - v2;
                    v0.unitize()?;
                    v0 *= tolerance * 1.1;
                }
                Ok(v0)
            };
            out[0] = extrapolate(0);
            if hu >= 1 {
                out[column] = tangent(column)?;
            }
            if hv >= 1 {
                out[1] = tangent(1)?;
            }
            if hu >= 2 {
                out[2 * column] = extrapolate(2 * column);
            }
            if (upper && highest_offset >= 2) || (!upper && hu >= 1 && hv >= 1) {
                out[column + 1] = extrapolate(column + 1);
            }
            if hv >= 2 {
                out[2] = extrapolate(2);
            }
            if hu >= 3 {
                out[3 * column] = extrapolate(3 * column);
            }
            if (upper && highest_offset >= 3) || (!upper && hu >= 2 && hv >= 1) {
                out[2 * column + 1] = extrapolate(2 * column + 1);
            }
            if (upper && highest_offset >= 3) || (!upper && hu >= 1 && hv >= 2) {
                out[column + 2] = extrapolate(column + 2);
            }
            if hv >= 3 {
                out[3] = extrapolate(3);
            }
            return Ok(out);
        }
        let mut n = bu.cross(bv);
        let mut length = n.length();
        if length < tolerance
            && let Ok([_, du, dv]) = first_derivatives(surface, uv, [true, true])
        {
            n = du.cross(dv);
            length = n.length();
            if !simple {
                bu = du;
                bv = dv;
            }
        }
        if !simple {
            let (bu_length, bv_length) = (bu.length(), bv.length());
            let smallest = length.min(bu_length).min(bv_length);
            if smallest < tolerance {
                let normal = evaluate_normal(surface, uv, sides);
                n = normal;
                length = n.length();
                if bu_length < bv_length {
                    let mut du = bv.cross(n);
                    let _ = du.unitize();
                    du *= tolerance * 1.1;
                    bu = du;
                } else {
                    let mut dv = n.cross(bu);
                    let _ = dv.unitize();
                    dv *= tolerance * 1.1;
                    bv = dv;
                }
            }
        }
        if length < EFF_ZERO {
            return Err(Fault::Degenerate);
        }
        let c = n;
        let m = c.length();
        let n = c / m;
        out[0] = b + offset * n;
        let high = highest_base;
        let c_at = |u: usize, v: usize| -> Vector3 {
            if u <= high && v <= high {
                cell(u, v)
            } else {
                Vector3::default()
            }
        };
        let (buu, buv, bvv) = (c_at(2, 0), c_at(1, 1), c_at(0, 2));
        let (buuu, buuv, buvv, bvvv) = (c_at(3, 0), c_at(2, 1), c_at(1, 2), c_at(0, 3));
        let (buuuu, buuuv, buuvv, buvvv, bvvvv) =
            (c_at(4, 0), c_at(3, 1), c_at(2, 2), c_at(1, 3), c_at(0, 4));
        let (mut cu, mut cv, mut cuu, mut cuv, mut cvv) = (
            Vector3::default(),
            Vector3::default(),
            Vector3::default(),
            Vector3::default(),
            Vector3::default(),
        );
        let (mut nu, mut nv, mut nuu, mut nuv, mut nvv) = (
            Vector3::default(),
            Vector3::default(),
            Vector3::default(),
            Vector3::default(),
            Vector3::default(),
        );
        let (mut mu, mut mv, mut muu, mut muv, mut mvv) = (0.0, 0.0, 0.0, 0.0, 0.0);
        if hu >= 1 {
            cu = buu.cross(bv) + bu.cross(buv);
            mu = n.dot(cu);
            nu = (cu - mu * n) / m;
            out[column] = bu + offset * nu;
        }
        if hv >= 1 {
            cv = buv.cross(bv) + bu.cross(bvv);
            mv = n.dot(cv);
            nv = (cv - mv * n) / m;
            out[1] = bv + offset * nv;
        }
        if hu >= 2 {
            cuu = buuu.cross(bv) + (2.0 * buu).cross(buv) + bu.cross(buuv);
            muu = (cu.dot(cu) + c.dot(cuu) - mu * mu) / m;
            nuu = (cuu - 2.0 * mu * nu - muu * n) / m;
            out[2 * column] = buu + offset * nuu;
        }
        if (upper && highest_offset >= 2) || (!upper && hu >= 1 && hv >= 1) {
            cuv = buuv.cross(bv) + buu.cross(bvv) + bu.cross(buvv);
            muv = (cu.dot(cv) + c.dot(cuv) - mu * mv) / m;
            nuv = (cuv - mu * nv - mv * nu - muv * n) / m;
            out[column + 1] = buv + offset * nuv;
        }
        if hv >= 2 {
            cvv = buvv.cross(bv) + (2.0 * buv).cross(bvv) + bu.cross(bvvv);
            mvv = (cv.dot(cv) + c.dot(cvv) - mv * mv) / m;
            nvv = (cvv - 2.0 * mv * nv - mvv * n) / m;
            out[2] = bvv + offset * nvv;
        }
        if hu >= 3 {
            let cuuu = buuuu.cross(bv)
                + (3.0 * buuu).cross(buv)
                + (3.0 * buu).cross(buuv)
                + bu.cross(buuuv);
            let muuu = (3.0 * cu.dot(cuu) + c.dot(cuuu) - 3.0 * mu * muu) / m;
            let nuuu = (cuuu - 3.0 * mu * nuu - 3.0 * muu * nu - muuu * n) / m;
            out[3 * column] = buuu + offset * nuuu;
        }
        if (upper && highest_offset >= 3) || (!upper && hu >= 2 && hv >= 1) {
            let cuuv = buuuv.cross(bv)
                + buuu.cross(bvv)
                + (3.0 * buuv).cross(buv)
                + (2.0 * buu).cross(buvv)
                + bu.cross(buuvv);
            let muuv =
                (cuu.dot(cv) + 2.0 * cu.dot(cuv) + c.dot(cuuv) - muu * mv - 2.0 * mu * muv) / m;
            let nuuv =
                (cuuv - mv * nuu - 2.0 * muv * nu - 2.0 * mu * nuv - muu * nv - muuv * n) / m;
            out[2 * column + 1] = buuv + offset * nuuv;
        }
        if (upper && highest_offset >= 3) || (!upper && hu >= 1 && hv >= 2) {
            let cuvv = buuvv.cross(bv)
                + (2.0 * buuv).cross(bvv)
                + (3.0 * buvv).cross(buv)
                + buu.cross(bvvv)
                + bu.cross(buvvv);
            let muvv =
                (cvv.dot(cu) + 2.0 * cv.dot(cuv) + c.dot(cuvv) - mvv * mu - 2.0 * mv * muv) / m;
            let nuvv =
                (cuvv - mu * nvv - 2.0 * muv * nv - 2.0 * mv * nuv - mvv * nu - muvv * n) / m;
            // SMLib stores `Buuv + offset * Nuvv` here, not `Buvv`.
            out[column + 2] = buuv + offset * nuvv;
        }
        if hv >= 3 {
            let cvvv = buvvv.cross(bv)
                + (3.0 * buvv).cross(bvv)
                + (3.0 * buv).cross(bvvv)
                + bu.cross(bvvvv);
            let mvvv = (3.0 * cv.dot(cvv) + c.dot(cvvv) - 3.0 * mv * mvv) / m;
            let nvvv = (cvvv - 3.0 * mv * nvv - 3.0 * mvv * nv - mvvv * n) / m;
            // SMLib stores `Buuv + offset * Nvvv` here, not `Bvvv`.
            out[3] = buuv + offset * nvvv;
        }
        if !simple {
            let tolerance = EFF_ZERO * (1.0 + out[0].max_dimension());
            if hu > 0
                && out[column].length_squared() < tolerance * tolerance
                && let Ok([_, mut du, _]) = first_derivatives(surface, uv, [true, true])
                && du.length_squared() > tolerance * tolerance
            {
                let _ = du.unitize();
                out[column] = du * tolerance * 1.1;
            }
            if hv > 0
                && out[1].length_squared() < tolerance * tolerance
                && let Ok([_, _, mut dv]) = first_derivatives(surface, uv, [true, true])
                && dv.length_squared() > tolerance * tolerance
            {
                let _ = dv.unitize();
                out[1] = dv * tolerance * 1.1;
            }
        }
        Ok(out)
    }
}

impl EvaluableSurface for OffsetSurface<'_> {
    fn normal_domain(&self) -> Extent2 {
        self.natural_domain()
    }

    fn plane_normal(&self) -> Option<Vector3> {
        None
    }

    fn evaluate_cells(
        &self,
        uv: Vector2,
        hu: usize,
        hv: usize,
        sides: [bool; 2],
        upper: bool,
        non_zero_tangents: bool,
        _zero_sampling: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        self.evaluate(uv, hu, hv, sides, upper, non_zero_tangents)
    }
}

impl SolverSurface for OffsetSurface<'_> {
    fn point_at(&self, uv: Vector2) -> Result<Vector3, Fault> {
        self.evaluate_point(uv)
    }

    fn normal_at(&self, uv: Vector2, sides: [bool; 2]) -> Vector3 {
        self.evaluate_normal(uv, sides)
    }

    fn simple_cells(
        &self,
        uv: Vector2,
        hu: usize,
        hv: usize,
        sides: [bool; 2],
        upper: bool,
    ) -> Result<Vec<Vector3>, Fault> {
        self.evaluate_simple(uv, hu, hv, sides, upper)
    }
}

/// `SmSurface::Evaluate1stDerivatives`: point, `Su`, `Sv`.
fn evaluate_first(surface: &BSplineSurface, uv: Vector2) -> Result<[Vector3; 3], Fault> {
    first_derivatives(surface, uv, [true, true])
}
