// SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES.
// SPDX-License-Identifier: Apache-2.0
// Modified for Catallax: written for the B-spline surface candidates of
// `spline.rs` without a trimming face; statuses become `Fault`, memory-block
// managers become owned vectors, tree nodes become an index arena in SMLib's
// allocation order, and the tessellation, trimming and non-B-spline branches
// are not carried.
// Adapted from rust-brep crates/kernel/src/surface/partial_net_decompose.rs
// (sm_decompose_srf), commit 4c83df59d777ccfcb14f441b1c8d984d1b5d93ab, a
// translation of SMLib SmNurbsSrf.cpp (sm_DecomposeSrf); the rest is
// translated directly from SMLib SmSurfaceCache.cpp (BuildTree,
// BuildTreeBase, SubdivideToTolerances, TestAgainstTolerances,
// SubdivideNode, BuildTreeTop), SmTrimSrfCache.cpp (BuildTree without a
// face), SmNurbsSrf.cpp (sm_ComputeNetConstants for a control net and the
// angle outputs of sm_ComputeSurfConstants),
// SmBSplineSurface.cpp (CalculateBoundingBox without the polar box) and
// SmSurface.cpp (CalculateContinuities, CacheMakeOrValidate for
// SM_OC_TRIMSRF).
//! The surface cache of step 5: a B-spline surface decomposed into Bezier
//! patches, subdivided to SMLib's chord-height and angle tolerances, under a
//! box tree for the global solvers.

use crate::Fault;
use crate::curve::calculate_continuities;
use crate::extent::{Extent2, Extent3};
use crate::geometry::{Continuity, line_closest_point};
use crate::global_solver::{AuxData, SearchTrees, TreeNode};
use crate::pseudo_box::PseudoBox;
use crate::scalar::{EFF_ZERO, EFF_ZERO_SQ, PI, ZONE_TOL_3D, acos, degrees, sqrt};
use crate::solution::NodeId;
use crate::spline::{BSplineSurface, combine, unique_knots};
use crate::surface_eval::{derivatives, eval_point, evaluate_normal, evaluate_point};
use crate::vector::{Vector2, Vector3};

/// Why a surface cache was not built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CacheRefusal {
    /// SMLib's `SM_ERR` or `SM_ERR_INVALID_INPUT`.
    Failed(Fault),
}

impl From<Fault> for CacheRefusal {
    fn from(fault: Fault) -> Self {
        Self::Failed(fault)
    }
}

/// `N_BasisGetSpanCount`: the knot spans of nonzero length.
fn span_count(knots: &[f64], degree: usize) -> usize {
    (degree..knots.len() - 1 - degree)
        .filter(|&i| knots[i] != knots[i + 1])
        .count()
}

/// `sm_DecomposeSrf`: the Bezier patches of `surface` by knot insertion,
/// strips in u first and patches in v within a strip, with the span counts.
pub fn decompose(surface: &BSplineSurface) -> Result<([usize; 2], Vec<BSplineSurface>), Fault> {
    let [p, q] = surface.degrees();
    let [u_knots, v_knots] = surface.knots();
    let m = surface.counts()[1] - 1;
    let (r, s) = (u_knots.len() - 1, v_knots.len() - 1);
    let pw = surface.control_points();
    let (nsu, nsv) = (span_count(u_knots, p), span_count(v_knots, q));
    if nsu < 1 || nsv < 1 {
        return Err(Fault::Degenerate);
    }
    let width = m + 1;
    let patch_width = q + 1;
    let mut bw = vec![pw[0]; (p + 1) * width];
    let mut nbw = bw.clone();
    let mut patches = vec![vec![pw[0]; (p + 1) * patch_width]; nsu * nsv];
    let mut patch_knots = vec![(vec![0.0; 2 * p + 2], vec![0.0; 2 * q + 2]); nsu * nsv];
    let (mut uals, mut omus) = (vec![0.0; p.max(1)], vec![0.0; p.max(1)]);
    let (mut vals, mut omvs) = (vec![0.0; q.max(1)], vec![0.0; q.max(1)]);
    let (mut isu, mut ieu) = (p, p + 1);
    let mut iq: isize = -1;
    bw[..(p + 1) * width].copy_from_slice(&pw[..(p + 1) * width]);
    while ieu < r {
        iq += 1;
        let start = ieu;
        while ieu < r && u_knots[ieu] == u_knots[ieu + 1] {
            ieu += 1;
        }
        let mlu = ieu - start + 1;
        let ru = p.saturating_sub(mlu);
        if mlu < p {
            let num = u_knots[ieu] - u_knots[isu];
            for i in (mlu + 1..=p).rev() {
                uals[i - mlu - 1] = num / (u_knots[isu + i] - u_knots[isu]);
                omus[i - mlu - 1] = 1.0 - uals[i - mlu - 1];
            }
            for i in 1..=ru {
                let su = mlu + i;
                let save = ru - i;
                for j in (su..=p).rev() {
                    for l in 0..=m {
                        bw[j * width + l] = combine(
                            uals[j - su],
                            bw[j * width + l],
                            omus[j - su],
                            bw[(j - 1) * width + l],
                        );
                    }
                }
                if ieu < r {
                    for l in 0..=m {
                        nbw[save * width + l] = bw[p * width + l];
                    }
                }
            }
        }
        let iq_index = iq as usize;
        if iq_index >= nsu {
            return Err(Fault::Degenerate);
        }
        let (mut isv, mut iev) = (q, q + 1);
        let mut jq: isize = -1;
        for j in 0..=q {
            for k in 0..=p {
                patches[iq_index * nsv][k * patch_width + j] = bw[k * width + j];
            }
        }
        while iev < s {
            jq += 1;
            let current = iq_index * nsv + jq as usize;
            let has_next = (jq as usize) < nsv - 1;
            let start = iev;
            while iev < s && v_knots[iev] == v_knots[iev + 1] {
                iev += 1;
            }
            let mlv = iev - start + 1;
            let rv = q.saturating_sub(mlv);
            if mlv < q {
                let num = v_knots[iev] - v_knots[isv];
                for i in (mlv + 1..=q).rev() {
                    vals[i - mlv - 1] = num / (v_knots[isv + i] - v_knots[isv]);
                    omvs[i - mlv - 1] = 1.0 - vals[i - mlv - 1];
                }
                for i in 1..=rv {
                    let sv = mlv + i;
                    let save = rv - i;
                    for j in (sv..=q).rev() {
                        for k in 0..=p {
                            let qw = &mut patches[current];
                            qw[k * patch_width + j] = combine(
                                vals[j - sv],
                                qw[k * patch_width + j],
                                omvs[j - sv],
                                qw[k * patch_width + j - 1],
                            );
                        }
                    }
                    if iev < s && has_next {
                        for k in 0..=p {
                            let value = patches[current][k * patch_width + q];
                            patches[current + 1][k * patch_width + save] = value;
                        }
                    }
                }
            }
            let (uq, vq) = &mut patch_knots[current];
            for i in 0..=p {
                uq[i] = u_knots[isu];
                uq[i + p + 1] = u_knots[ieu];
            }
            for j in 0..=q {
                vq[j] = v_knots[isv];
                vq[j + q + 1] = v_knots[iev];
            }
            if iev < s && has_next {
                for i in rv..=q {
                    for k in 0..=p {
                        patches[current + 1][k * patch_width + i] = bw[k * width + iev - q + i];
                    }
                }
            }
            isv = iev;
            iev += 1;
        }
        if ieu < r {
            for i in ru..=p {
                for l in 0..=m {
                    nbw[i * width + l] = pw[(ieu - p + i) * width + l];
                }
            }
        }
        isu = ieu;
        ieu += 1;
        std::mem::swap(&mut bw, &mut nbw);
    }
    let beziers = patches
        .into_iter()
        .zip(patch_knots)
        .map(|(points, (uq, vq))| BSplineSurface::new([p, q], [uq, vq], [p + 1, q + 1], points))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(([nsu, nsv], beziers))
}

/// `SmBSplineSurface::CalculateBoundingBox` without the polar box: the
/// axis-aligned and pseudo boxes of the whole control net, the pseudo basis
/// from the net corners or, when tighter, from `guess`, both grown by the
/// zone tolerance.
pub fn bounding_boxes(surface: &BSplineSurface, guess: Option<&PseudoBox>) -> (Extent3, PseudoBox) {
    let [nu, nv] = surface.counts();
    let point = |i: usize, j: usize| surface.control_points()[i * nv + j].euclidean();
    let (p00, p10, p01, p11) = (
        point(0, 0),
        point(nu - 1, 0),
        point(0, nv - 1),
        point(nu - 1, nv - 1),
    );
    let v1 = (p00 - p10) + (p01 - p11);
    let v2 = (p00 - p01) + (p10 - p11);
    let mut pseudo = PseudoBox::new();
    if v1.length_squared() > EFF_ZERO_SQ
        && let Ok([b1, b2, b3]) = v1.make_unit_ortho_vectors(Some(v2))
    {
        let _ = pseudo.set_basis(b1, b2, b3);
    }
    let mut trial = guess.map(|guess| {
        let mut trial = PseudoBox::new();
        let _ = trial.set_basis(guess.basis[0], guess.basis[1], guess.basis[2]);
        trial
    });
    let mut bbox = Extent3::singleton(point(0, 0));
    for i in 0..nu {
        for j in 0..nv {
            let p = point(i, j);
            if i != 0 || j != 0 {
                bbox.add_point(p);
            }
            pseudo.add_point(p);
            if let Some(trial) = &mut trial {
                trial.add_point(p);
            }
        }
    }
    if let Some(trial) = trial
        && trial.volume() < pseudo.volume() - EFF_ZERO
    {
        pseudo = trial;
    }
    bbox.expand_absolute(ZONE_TOL_3D);
    pseudo.expand_absolute(ZONE_TOL_3D);
    (bbox, pseudo)
}

/// `SmSurface::CalculateContinuities` along `axis`: the class at every
/// unique knot of `axis` from iso-curves at up to five parameters across,
/// with kinks whose normals agree within one degree raised to G1, and the
/// smallest interior class.
pub fn continuities(
    surface: &BSplineSurface,
    axis: usize,
) -> Result<(Continuity, Vec<Continuity>), Fault> {
    let other_knots: Vec<f64> = unique_knots(surface.knots()[axis], None)
        .iter()
        .map(|k| k.0)
        .collect();
    if other_knots.len() <= 1 {
        return Err(Fault::Degenerate);
    }
    let count = other_knots.len();
    let mut classes = vec![Continuity::G1; count];
    classes[0] = Continuity::Discontinuous;
    classes[count - 1] = Continuity::Discontinuous;
    let mut minimum = Continuity::G1;
    let across = 1 - axis;
    let mut knots: Vec<f64> = unique_knots(surface.knots()[across], None)
        .iter()
        .map(|k| k.0)
        .collect();
    if knots.len() == 2 {
        knots.push(0.5 * (knots[0] + knots[1]));
    }
    if knots.len() > 5 {
        let (first, last) = (knots[0], knots[knots.len() - 1]);
        knots.truncate(5);
        knots[4] = last;
        knots[3] = first + 0.845 * (last - first);
        knots[2] = 0.5 * (first + last);
        knots[1] = first + 0.345 * (last - first);
    }
    let mut first = true;
    for i in 0..knots.len() {
        let mut parameter = knots[i];
        let mut iso = surface.iso_curve(across, parameter)?;
        let iso_knots: Vec<f64> = unique_knots(iso.knots(), None)
            .iter()
            .map(|k| k.0)
            .collect();
        if iso.is_degenerate(EFF_ZERO) {
            if i > 0 && first {
                parameter = (knots[i - 1] + knots[i]) / 2.0;
                iso = surface.iso_curve(across, parameter)?;
            }
            if iso.is_degenerate(EFF_ZERO) {
                continue;
            }
        }
        let (_, mut iso_classes) = calculate_continuities(&iso);
        let mut iso_minimum = Continuity::CInfinity;
        for j in 1..iso_classes.len().saturating_sub(1) {
            if iso_classes[j] < Continuity::G1 {
                let uv = if axis == 0 {
                    Vector2::new(iso_knots[j], parameter)
                } else {
                    Vector2::new(parameter, iso_knots[j])
                };
                let n1 = evaluate_normal(surface, uv, [false, false]);
                let n2 = evaluate_normal(surface, uv, [true, true]);
                if n1.dot(n2) > 0.9998476951563913 {
                    iso_classes[j] = Continuity::G1;
                }
            }
            if iso_classes[j] < iso_minimum {
                iso_minimum = iso_classes[j];
            }
        }
        let same = iso_classes.len() == count;
        let find = |value: f64| other_knots.iter().position(|&k| k == value);
        if first {
            first = false;
            if same {
                classes = iso_classes;
            } else {
                classes = vec![Continuity::Undefined; count];
                for (j, &class) in iso_classes.iter().enumerate() {
                    if let Some(index) = find(iso_knots[j]) {
                        classes[index] = class;
                    }
                }
            }
            minimum = iso_minimum;
        } else {
            if iso_minimum < minimum {
                minimum = iso_minimum;
            }
            for (j, &class) in iso_classes.iter().enumerate() {
                if same {
                    if class < classes[j] {
                        classes[j] = class;
                    }
                } else if let Some(index) = find(iso_knots[j])
                    && class < classes[index]
                {
                    classes[index] = class;
                }
            }
        }
    }
    Ok((minimum, classes))
}

/// What `sm_ComputeNetConstants` measures on a control net: the largest
/// distances of control points from the chords of the net rows (u), net
/// columns (v) and net diagonals (uv), and the largest accumulated turning
/// angles of the rows and columns, in degrees.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NetConstants {
    /// Largest row chord height.
    pub u_chord_height: f64,
    /// Largest column chord height.
    pub v_chord_height: f64,
    /// Largest diagonal chord height.
    pub uv_chord_height: f64,
    /// Largest row turning angle.
    pub u_angle_deg: f64,
    /// Largest column turning angle.
    pub v_angle_deg: f64,
}

/// The chord heights and turning angles along one direction of the net.
/// `line(k, i)` is control point `i` of line `k`; `lines` lines of
/// `last + 1` points.
fn net_lines(
    point: &dyn Fn(usize, usize) -> Vector3,
    lines: usize,
    last: usize,
    chords: bool,
    angles: bool,
) -> (f64, f64) {
    const LIMIT_RAD: f64 = 10.0;
    let mut chord_squared = 0.0;
    let mut max_angle = 0.0;
    let mut over_limit = false;
    for k in 0..lines {
        let (mut line_point, mut line_vector) = (Vector3::default(), Vector3::default());
        if chords {
            line_point = point(k, 0);
            line_vector = point(k, last) - line_point;
            let mut try_index = last.wrapping_sub(1);
            while line_vector.length_squared() < EFF_ZERO_SQ {
                line_vector = point(k, try_index) - line_point;
                if try_index == 0 {
                    break;
                }
                try_index -= 1;
            }
            if line_vector.length_squared() < EFF_ZERO_SQ {
                continue;
            }
        }
        let mut total = 0.0;
        let mut last_vector = Vector3::new(0.0, 0.0, 0.0);
        let mut current = point(k, 0);
        for i in 1..=last {
            let next = point(k, i);
            if chords {
                let t = line_closest_point(line_point, line_vector, next);
                let height = next - (line_point + line_vector * t);
                let squared = height.length_squared();
                if squared > chord_squared {
                    chord_squared = squared;
                }
            }
            if angles {
                if total > LIMIT_RAD {
                    if !over_limit {
                        over_limit = true;
                        total *= (last / i) as f64;
                        if !chords {
                            break;
                        }
                    }
                    continue;
                }
                let step = next - current;
                let distance = step.length();
                if distance < ZONE_TOL_3D {
                    continue;
                }
                let step = step / distance;
                if last_vector.length_squared() > EFF_ZERO_SQ {
                    let cosine = step.dot(last_vector).clamp(-1.0, 1.0);
                    total += acos(cosine);
                }
                last_vector = step;
            }
            current = next;
        }
        if total > max_angle {
            max_angle = total;
        }
        if over_limit && !chords {
            break;
        }
    }
    (chord_squared, max_angle)
}

/// `sm_ComputeNetConstants` for a control net (`pSurface` null), with
/// SMLib's zone tolerance.
pub fn net_constants(
    surface: &BSplineSurface,
    chords: bool,
    angles: bool,
    diagonals: bool,
) -> NetConstants {
    let [nu, nv] = surface.counts();
    let (n, m) = (nu - 1, nv - 1);
    let net = |i: usize, j: usize| surface.control_points()[i * nv + j].euclidean();
    let (v_chord, v_angle) = net_lines(&|k, i| net(k, i), n + 1, m, chords, angles);
    let (u_chord, u_angle) = net_lines(&|k, i| net(i, k), m + 1, n, chords, angles);
    let mut uv_chord = 0.0;
    if diagonals {
        let domain = surface.natural_domain();
        let (ud, vd) = (domain.u_interval(), domain.v_interval());
        let at = |u: f64, v: f64| {
            eval_point(
                surface,
                Vector2::new(ud.evaluate(u), vd.evaluate(v)),
                [true, true],
            )
            .unwrap_or_default()
        };
        for jj in 0..2 {
            let ds = if jj == 0 { -0.25 } else { 0.25 };
            let mut dv = jj as f64;
            let line_point = at(0.0, dv);
            let last_point = at(1.0, 1.0 - dv);
            let mut line_vector = last_point - line_point;
            let mut tries = 0;
            while tries <= 3 && line_vector.length_squared() < EFF_ZERO_SQ {
                let du = 1.0 - 0.25 * tries as f64;
                dv = (1 - jj) as f64 + ds * tries as f64;
                line_vector = at(du, dv) - line_point;
                tries += 1;
            }
            if tries == 4 {
                continue;
            }
            for itr in 1..=3 {
                let du = 0.25 * itr as f64;
                let dv = jj as f64 - ds * itr as f64;
                let next = at(du, dv);
                let t = line_closest_point(line_point, line_vector, next);
                let squared = (next - (line_point + line_vector * t)).length_squared();
                if squared > uv_chord {
                    uv_chord = squared;
                }
            }
        }
        let count = m.min(n);
        for jj in 0..2 {
            let (start, end, step): (usize, usize, isize) = if jj == 0 {
                (0, count, 1)
            } else {
                (m, m - count, -1)
            };
            let line_point = net(0, start);
            let mut line_vector = net(count, end) - line_point;
            let mut itr = 1;
            while itr < count && line_vector.length_squared() < EFF_ZERO_SQ {
                let column = (end as isize - step * itr as isize) as usize;
                line_vector = net(count - itr, column) - line_point;
                itr += 1;
            }
            if itr == count {
                continue;
            }
            for itr in 1..=count {
                let column = (start as isize + step * itr as isize) as usize;
                let next = net(itr, column);
                let t = line_closest_point(line_point, line_vector, next);
                let squared = (next - (line_point + line_vector * t)).length_squared();
                if squared > uv_chord {
                    uv_chord = squared;
                }
            }
        }
    }
    NetConstants {
        u_chord_height: sqrt(u_chord),
        v_chord_height: sqrt(v_chord),
        uv_chord_height: sqrt(uv_chord),
        u_angle_deg: degrees(u_angle),
        v_angle_deg: degrees(v_angle),
    }
}

/// `SmBezierPatch`: one Bezier patch of the decomposition or of its
/// subdivision.
#[derive(Clone, Debug, PartialEq)]
pub struct PatchData {
    /// The patch (SMLib frees it after construction; the gate compares the
    /// rest).
    pub bezier: BSplineSurface,
    /// The parallelepiped bounding the patch.
    pub pseudo_box: PseudoBox,
    /// The parameter domain covered.
    pub uv_domain: Extent2,
    /// Continuity at the low and high u sides.
    pub u_continuities: [Continuity; 2],
    /// Continuity at the low and high v sides.
    pub v_continuities: [Continuity; 2],
    /// The direction its children split, 0 for u and 1 for v.
    pub split_axis: usize,
    /// The chord height tolerance held for this patch and its descendants.
    pub chord_height_satisfied: bool,
    /// The angle tolerance held for this patch and its descendants.
    pub angle_satisfied: bool,
}

/// `SmBezierAux2d` of an organising node above the first patches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeadData {
    /// The parameter domain covered.
    pub uv_domain: Extent2,
    /// The direction its children split, 0 for u and 1 for v.
    pub split_axis: usize,
}

/// What a surface-cache node's data refers to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceNodeData {
    /// A Bezier patch node: the index into the patches.
    Patch(usize),
    /// An organising node: the index into the heads.
    Head(usize),
}

/// The subdivision tolerances of an `SmSurfaceCache`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceCacheSettings {
    /// Largest pseudo-box thickness and control-point chord height; 0 ignores.
    pub chord_height_tolerance: f64,
    /// Largest control-net turning angle, in radians; 0 ignores.
    pub angle_tolerance_rad: f64,
    /// Smallest node side as a fraction of the domain; 0 means 0.001.
    pub min_side_length_ratio_uv: f64,
}

/// `SmSurfaceCache` (as `SmTrimSrfCache` builds it without a face) over a
/// B-spline surface.
#[derive(Clone, Debug)]
pub struct SurfaceCache {
    /// The parameter domain covered.
    pub domain: Extent2,
    /// The tolerances the tree was built to.
    pub settings: SurfaceCacheSettings,
    /// The continuity class at each unique u knot.
    pub u_continuities: Vec<Continuity>,
    /// The continuity class at each unique v knot.
    pub v_continuities: Vec<Continuity>,
    /// The first-generation patch counts in u and v.
    pub first_generation: [usize; 2],
    /// Every patch, in SMLib's allocation order.
    pub patches: Vec<PatchData>,
    /// Every organising node's data, in allocation order.
    pub heads: Vec<HeadData>,
    /// The tree, in SMLib's node allocation order.
    pub trees: SearchTrees,
    /// The data of each tree node, by node identity.
    pub node_data: Vec<SurfaceNodeData>,
    /// The parent of each tree node, by node identity.
    pub parents: Vec<Option<NodeId>>,
    /// The root of the tree.
    pub top: NodeId,
}

/// The angle outputs of `sm_ComputeSurfConstants`: the largest turning angle,
/// in degrees, of the polyline through five samples of the surface over
/// `domain` along each of five rows (V) and five columns (U). A pair of samples
/// closer than the scaled zero ends the walk along that line, and SMLib's two
/// loops skip the first angle differently: the V walk by whether the last
/// direction is nonzero, the U walk by the sample index.
fn sampled_angles(surface: &BSplineSurface, domain: Extent2) -> Result<(f64, f64), Fault> {
    let at = |u: f64, v: f64| evaluate_point(surface, domain.evaluate(u, v));
    let mut max_v = 0.0;
    let mut max_u = 0.0;
    for jj in 0..=4 {
        let u = f64::from(jj) / 4.0;
        let mut total = 0.0;
        let mut last = Vector3::new(0.0, 0.0, 0.0);
        let mut current = at(u, 0.0)?;
        for ii in 1..=4 {
            let next = at(u, f64::from(ii) / 4.0)?;
            let mut direction = next - current;
            let distance = direction.length();
            if distance < EFF_ZERO {
                break;
            }
            direction /= distance;
            if last.length_squared() > EFF_ZERO_SQ {
                total += acos(direction.dot(last).clamp(-1.0, 1.0));
            }
            last = direction;
            current = next;
        }
        if total > max_v {
            max_v = total;
        }
    }
    for jj in 0..=4 {
        let v = f64::from(jj) / 4.0;
        let mut total = 0.0;
        let mut last = Vector3::new(0.0, 0.0, 0.0);
        let mut current = at(0.0, v)?;
        for ii in 1..=4 {
            let next = at(f64::from(ii) / 4.0, v)?;
            let mut direction = next - current;
            let distance = direction.length();
            if distance < EFF_ZERO {
                break;
            }
            direction /= distance;
            if ii > 1 {
                total += acos(direction.dot(last).clamp(-1.0, 1.0));
            }
            last = direction;
            current = next;
        }
        if total > max_u {
            max_u = total;
        }
    }
    Ok((degrees(max_u), degrees(max_v)))
}

/// A node while the tree is built.
struct BuildNode {
    bbox: Extent3,
    children: Option<[usize; 2]>,
    parent: Option<usize>,
    data: SurfaceNodeData,
}

/// The state of `BuildTree`.
struct Builder<'a> {
    /// `m_cpSurface`: the whole surface, which the sampled angles use.
    surface: &'a BSplineSurface,
    domain: Extent2,
    settings: SurfaceCacheSettings,
    nodes: Vec<BuildNode>,
    patches: Vec<PatchData>,
    heads: Vec<HeadData>,
}

impl Builder<'_> {
    /// `TestAgainstTolerances` for a patch node: whether to subdivide and
    /// along which axis.
    fn test_against_tolerances(&mut self, node: usize) -> Result<Option<usize>, CacheRefusal> {
        let SurfaceNodeData::Patch(index) = self.nodes[node].data else {
            return Err(Fault::Degenerate.into());
        };
        if self.nodes[node].children.is_some() {
            return Err(Fault::Degenerate.into());
        }
        let settings = self.settings;
        let chord_squared = settings.chord_height_tolerance * settings.chord_height_tolerance;
        let angle_tolerance_deg = settings.angle_tolerance_rad * 180.0 / PI;
        let patch = &mut self.patches[index];
        if chord_squared < EFF_ZERO {
            patch.chord_height_satisfied = true;
        }
        if angle_tolerance_deg < EFF_ZERO {
            patch.angle_satisfied = true;
        }
        let mut subdivide = false;
        let mut axis = 0;
        if chord_squared >= EFF_ZERO {
            let height = patch.pseudo_box.intervals[2].length();
            if height * height <= chord_squared {
                patch.chord_height_satisfied = true;
            } else {
                subdivide = true;
                let constants = net_constants(&patch.bezier, true, false, true);
                let tolerance = settings.chord_height_tolerance;
                if constants.u_chord_height <= tolerance
                    && constants.v_chord_height <= tolerance
                    && constants.uv_chord_height > tolerance
                {
                    let domain = patch.bezier.natural_domain();
                    let (u, v) = (domain.u_interval(), domain.v_interval());
                    let middle = Vector2::new((u.min + u.max) / 2.0, (v.min + v.max) / 2.0);
                    let d = derivatives(&patch.bezier, middle, [true, true], true, 1, 1)?;
                    if d[2].length_squared() * (u.max - u.min)
                        < d[1].length_squared() * (v.max - v.min)
                    {
                        axis = 1;
                    }
                } else if constants.v_chord_height > constants.u_chord_height {
                    axis = 1;
                }
            }
        }
        if !subdivide && angle_tolerance_deg >= EFF_ZERO {
            let constants = net_constants(&patch.bezier, false, true, false);
            let (mut u_angle, mut v_angle) = (constants.u_angle_deg, constants.v_angle_deg);
            // More than one and a half turns says the control net is no
            // guide to the angle: use surface evaluations. [B99]
            if u_angle > 540.0 || v_angle > 540.0 {
                (u_angle, v_angle) = sampled_angles(self.surface, patch.uv_domain)?;
            }
            if angle_tolerance_deg < EFF_ZERO
                || (u_angle <= angle_tolerance_deg && v_angle <= angle_tolerance_deg)
            {
                patch.angle_satisfied = true;
            } else {
                if v_angle > u_angle {
                    axis = 1;
                }
                subdivide = true;
            }
        }
        if !subdivide {
            return Ok(None);
        }
        let whole = self.domain.size();
        let size = patch.uv_domain.size();
        let ratio = if settings.min_side_length_ratio_uv != 0.0 {
            settings.min_side_length_ratio_uv
        } else {
            0.001
        };
        if (axis == 0 && size.x < whole.x * ratio) || (axis == 1 && size.y < whole.y * ratio) {
            if axis == 0 && size.y > whole.y / 100.0 {
                axis = 1;
            } else if axis == 1 && size.x > whole.x / 100.0 {
                axis = 0;
            } else {
                return Ok(None);
            }
        }
        Ok(Some(axis))
    }

    /// `SubdivideNode` for a patch node: halves its domain along `axis` and
    /// splits its patch into two child patch nodes.
    fn subdivide_node(&mut self, node: usize, axis: usize) -> Result<[usize; 2], CacheRefusal> {
        let SurfaceNodeData::Patch(index) = self.nodes[node].data else {
            return Err(Fault::Degenerate.into());
        };
        let parent = &mut self.patches[index];
        parent.split_axis = axis;
        let (low, high) = (parent.uv_domain.min, parent.uv_domain.max);
        let (mut domain1, mut domain2) = (parent.uv_domain, parent.uv_domain);
        let split = if axis == 0 {
            let split = (high.x + low.x) / 2.0;
            domain1.set_min_max(low, Vector2::new(split, high.y))?;
            domain2.set_min_max(Vector2::new(split, low.y), high)?;
            split
        } else {
            let split = (high.y + low.y) / 2.0;
            domain1.set_min_max(low, Vector2::new(high.x, split))?;
            domain2.set_min_max(Vector2::new(low.x, split), high)?;
            split
        };
        let [bezier1, bezier2] = parent.bezier.split(axis, split)?;
        let guess = parent.pseudo_box;
        let (bbox1, pseudo1) = bounding_boxes(&bezier1, Some(&guess));
        let (bbox2, pseudo2) = bounding_boxes(&bezier2, Some(&guess));
        let mut child1 = PatchData {
            bezier: bezier1,
            pseudo_box: pseudo1,
            uv_domain: domain1,
            u_continuities: parent.u_continuities,
            v_continuities: parent.v_continuities,
            split_axis: 0,
            chord_height_satisfied: parent.chord_height_satisfied,
            angle_satisfied: parent.angle_satisfied,
        };
        let mut child2 = PatchData {
            bezier: bezier2,
            pseudo_box: pseudo2,
            uv_domain: domain2,
            ..child1.clone()
        };
        if axis == 0 {
            child1.u_continuities[1] = Continuity::C1C2;
            child2.u_continuities[0] = Continuity::C1C2;
        } else {
            child1.v_continuities[1] = Continuity::C1C2;
            child2.v_continuities[0] = Continuity::C1C2;
        }
        let mut ids = [0; 2];
        for (k, (child, bbox)) in [(child1, bbox1), (child2, bbox2)].into_iter().enumerate() {
            self.patches.push(child);
            self.nodes.push(BuildNode {
                bbox,
                children: None,
                parent: Some(node),
                data: SurfaceNodeData::Patch(self.patches.len() - 1),
            });
            ids[k] = self.nodes.len() - 1;
        }
        self.nodes[node].children = Some(ids);
        Ok(ids)
    }

    /// `SubdivideToTolerances`: subdivides `node` until every leaf meets
    /// the tolerances, then bounds it by its children.
    fn subdivide_to_tolerances(&mut self, node: usize) -> Result<(), CacheRefusal> {
        let Some(axis) = self.test_against_tolerances(node)? else {
            return Ok(());
        };
        let [child1, child2] = self.subdivide_node(node, axis)?;
        self.subdivide_to_tolerances(child1)?;
        self.subdivide_to_tolerances(child2)?;
        self.nodes[node].bbox = self.nodes[child1].bbox.union(self.nodes[child2].bbox);
        Ok(())
    }

    /// `BuildTreeTop`: pairs neighbouring first-generation nodes along the
    /// direction whose pairing makes the smaller box, until one root is left.
    fn build_tree_top(&mut self, first: Vec<usize>, [mut nu, mut nv]: [usize; 2]) -> usize {
        let mut current = first;
        while nu > 1 || nv > 1 {
            let mut axis = 0;
            if nu == 1 {
                axis = 1;
            } else if nv != 1 {
                let b00 = self.nodes[current[0]].bbox;
                let ext_v = b00.union(self.nodes[current[1]].bbox);
                let ext_u = b00.union(self.nodes[current[nv]].bbox);
                let diff_u = ext_u.max - ext_u.min;
                let diff_v = ext_v.max - ext_v.min;
                if diff_v.length_squared() < diff_u.length_squared() {
                    axis = 1;
                }
            }
            let mut next = Vec::new();
            let pair = |builder: &mut Self, a: usize, b: usize, axis: usize| {
                let domain = builder.uv_domain(a).union(builder.uv_domain(b));
                builder.heads.push(HeadData {
                    uv_domain: domain,
                    split_axis: axis,
                });
                let bbox = builder.nodes[a].bbox.union(builder.nodes[b].bbox);
                builder.nodes.push(BuildNode {
                    bbox,
                    children: Some([a, b]),
                    parent: None,
                    data: SurfaceNodeData::Head(builder.heads.len() - 1),
                });
                let id = builder.nodes.len() - 1;
                builder.nodes[a].parent = Some(id);
                builder.nodes[b].parent = Some(id);
                id
            };
            if axis == 1 {
                for i in 0..nu {
                    for k in (1..nv).step_by(2) {
                        let id = pair(self, current[i * nv + k - 1], current[i * nv + k], 1);
                        next.push(id);
                    }
                    if nv % 2 == 1 {
                        next.push(current[i * nv + nv - 1]);
                    }
                }
                nv = nv.div_ceil(2);
            } else {
                for i in (1..nu).step_by(2) {
                    for k in 0..nv {
                        let id = pair(self, current[(i - 1) * nv + k], current[i * nv + k], 0);
                        next.push(id);
                    }
                }
                if nu % 2 == 1 {
                    for k in 0..nv {
                        next.push(current[(nu - 1) * nv + k]);
                    }
                }
                nu = nu.div_ceil(2);
            }
            current = next;
        }
        current[0]
    }

    fn uv_domain(&self, node: usize) -> Extent2 {
        match self.nodes[node].data {
            SurfaceNodeData::Patch(index) => self.patches[index].uv_domain,
            SurfaceNodeData::Head(index) => self.heads[index].uv_domain,
        }
    }
}

impl SurfaceCache {
    /// The parameter domain of a node of the cache's own tree.
    pub fn uv_domain(&self, id: NodeId) -> Extent2 {
        match self.node_data[id.0 as usize] {
            SurfaceNodeData::Patch(index) => self.patches[index].uv_domain,
            SurfaceNodeData::Head(index) => self.heads[index].uv_domain,
        }
    }

    /// `SmSurfaceCache::BuildTree` over the natural domain: the Bezier
    /// decomposition, each patch subdivided to `settings`, and the tree top.
    pub fn build(
        surface: &BSplineSurface,
        settings: SurfaceCacheSettings,
    ) -> Result<Self, CacheRefusal> {
        crate::surface_eval::check_supported(surface)?;
        let domain = surface.natural_domain();
        let (counts, beziers) = decompose(surface)?;
        let (_, u_continuities) = continuities(surface, 0)?;
        let (_, v_continuities) = continuities(surface, 1)?;
        let mut builder = Builder {
            surface,
            domain,
            settings,
            nodes: Vec::new(),
            patches: Vec::new(),
            heads: Vec::new(),
        };
        let mut first = Vec::new();
        for (index, bezier) in beziers.into_iter().enumerate() {
            let (i, j) = (index / counts[1], index % counts[1]);
            let (bbox, pseudo_box) = bounding_boxes(&bezier, None);
            let uv_domain = bezier.natural_domain();
            builder.patches.push(PatchData {
                bezier,
                pseudo_box,
                uv_domain,
                u_continuities: [u_continuities[i], u_continuities[i + 1]],
                v_continuities: [v_continuities[j], v_continuities[j + 1]],
                split_axis: 0,
                chord_height_satisfied: false,
                angle_satisfied: false,
            });
            builder.nodes.push(BuildNode {
                bbox,
                children: None,
                parent: None,
                data: SurfaceNodeData::Patch(index),
            });
            first.push(index);
        }
        for node in 0..first.len() {
            builder.subdivide_to_tolerances(node)?;
        }
        let top = builder.build_tree_top(first, counts);
        let mut trees = SearchTrees::new();
        let mut node_data = Vec::new();
        let mut parents = Vec::new();
        for node in &builder.nodes {
            trees.add(TreeNode {
                bbox: node.bbox,
                children: node
                    .children
                    .map(|[a, b]| [NodeId(a as u32), NodeId(b as u32)]),
                aux: match node.data {
                    SurfaceNodeData::Patch(_) => AuxData::BezierSurface,
                    SurfaceNodeData::Head(_) => AuxData::BezierSurfaceHead,
                },
            });
            node_data.push(node.data);
            parents.push(node.parent.map(|p| NodeId(p as u32)));
        }
        Ok(Self {
            domain,
            settings,
            u_continuities,
            v_continuities,
            first_generation: counts,
            patches: builder.patches,
            heads: builder.heads,
            trees,
            node_data,
            parents,
            top: NodeId(top as u32),
        })
    }

    /// `SmSurface::CacheMakeOrValidate` for `SM_OC_TRIMSRF` without a face:
    /// the chord height is a fortieth of a size estimate from five surface
    /// points and the angle 20 degrees, or a fifteenth and 40 degrees when
    /// both degrees are at most two.
    pub fn for_trimmed_surface(surface: &BSplineSurface) -> Result<Self, CacheRefusal> {
        let domain = surface.natural_domain();
        let at = |uv: Vector2| evaluate_point(surface, uv);
        let middle = at(domain.evaluate(0.5, 0.5))?;
        let p00 = at(domain.min)?;
        let p10 = at(domain.evaluate(0.9, 0.1))?;
        let p01 = at(domain.evaluate(0.1, 0.9))?;
        let p11 = at(domain.max)?;
        let size = (middle.distance(p00)
            + middle.distance(p11)
            + middle.distance(p10)
            + middle.distance(p01))
            / 2.0;
        let mut chord = size / 40.0;
        let mut angle = 20.0 * PI / 180.0;
        let [p, q] = surface.degrees();
        if p <= 2 && q <= 2 {
            chord = size / 15.0;
            angle *= 2.0;
        }
        Self::build(
            surface,
            SurfaceCacheSettings {
                chord_height_tolerance: chord,
                angle_tolerance_rad: angle,
                min_side_length_ratio_uv: 0.001,
            },
        )
    }
}

/// One side's plane from `SmPatchBoundaryBoundingPlanes`: the side curve
/// lies in it and the whole control net on one side of it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryPlane {
    /// A point of the plane.
    pub base: Vector3,
    /// The plane normal, pointing away from the control net.
    pub normal: Vector3,
    /// The largest distance of the side curve's samples and polygon from
    /// the plane.
    pub tolerance: f64,
    /// The side curve at parameter fractions 0, 1/4, 1/2, 3/4 and 1.
    pub samples: [Vector3; 5],
}

/// `SmPatchBoundaryBoundingPlanes` after `ComputeBoundingPlanes`, for the
/// sides in SMLib's order: v low, u high, v high, u low.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryPlanes {
    /// The plane of each side, where one exists.
    pub planes: [Option<BoundaryPlane>; 4],
    /// Whether each side curve collapses to a point.
    pub degenerate: [bool; 4],
}

/// `sm_ComputePlaneFrom3Points`: the plane through three side samples
/// unless they are too close or nearly collinear (`Ok(None)`), with
/// whether they are collinear.
fn plane_from_3_points(
    p1: Vector3,
    p2: Vector3,
    p3: Vector3,
    tolerance: f64,
) -> Result<(Option<(Vector3, Vector3)>, bool), Fault> {
    let (mut a, mut b, mut c) = (p1 - p2, p2 - p3, p3 - p1);
    if a.length() < tolerance || b.length() < tolerance || c.length() < tolerance {
        return Ok((None, false));
    }
    let _ = a.unitize();
    let _ = b.unitize();
    let _ = c.unitize();
    let angle_tolerance = 0.01;
    let (n1, n2, n3) = (a.cross(b), b.cross(c), c.cross(a));
    let (s1, s2, s3) = (n1.length(), n2.length(), n3.length());
    if s1 < angle_tolerance || s2 < angle_tolerance || s3 < angle_tolerance {
        let collinear = s1 < crate::scalar::EFF_ZERO_SQRT
            && s2 < crate::scalar::EFF_ZERO_SQRT
            && s3 < crate::scalar::EFF_ZERO_SQRT;
        return Ok((None, collinear));
    }
    if n1.dot(n2) < 0.0 || n2.dot(n3) < 0.0 || n3.dot(n1) < 0.0 {
        return Err(Fault::Degenerate);
    }
    let base = (p1 + p2 + p3) / 3.0;
    let mut normal = (n1 + n2 + n3) / 3.0;
    let _ = normal.unitize();
    Ok((Some((base, normal)), false))
}

/// `sm_CheckBoundingPlaneSeparatingHalfspace`: whether every control point
/// of `surface` off the plane lies on one side, and whether that side is
/// the normal's.
fn separating_halfspace(surface: &BSplineSurface, base: Vector3, normal: Vector3) -> (bool, bool) {
    let mut sign = 0;
    for point in surface.control_points() {
        let [x, y, z, w] = point.coordinates;
        let mut p = Vector3::new(x, y, z);
        if w != crate::spline::NO_WEIGHT {
            p = Vector3::new(x / w, y / w, z / w);
        }
        let distance = normal.dot(p - base);
        if distance.abs() < crate::scalar::EFF_ZERO_SQRT {
            continue;
        }
        if (distance > 0.0 && sign == -1) || (distance < 0.0 && sign == 1) {
            return (false, false);
        }
        sign = if distance > 0.0 { 1 } else { -1 };
    }
    (true, sign == 1)
}

/// `SmBSplineCurve::IsOnPlane` over the natural interval: whether the ends,
/// middle and control polygon lie within `tolerance` of the plane, and the
/// largest of those distances.
fn curve_on_plane(
    curve: &crate::spline::BSplineCurve,
    base: Vector3,
    normal: Vector3,
    tolerance: f64,
) -> Result<(bool, f64), Fault> {
    use crate::geometry::{plane_point_distance, point_project_to_plane};
    let interval = curve.natural_interval();
    let start = curve.evaluate_point(interval.min)?;
    let end = curve.evaluate_point(interval.max)?;
    let middle = curve.evaluate_point(interval.evaluate(0.5))?;
    let gap = |p: Vector3| (p - point_project_to_plane(p, base, normal)).length();
    let (d_start, d_end, d_middle) = (gap(start), gap(end), gap(middle));
    if d_start > tolerance || d_end > tolerance || d_middle > tolerance {
        return Ok((false, 0.0));
    }
    // smos_3Max(start, middle, end).
    let mut largest = if d_start > d_middle {
        if d_start > d_end { d_start } else { d_end }
    } else if d_middle > d_end {
        d_middle
    } else {
        d_end
    };
    for point in crate::curve::polygon(curve, interval)? {
        let distance = plane_point_distance(base, normal, point)?;
        if distance > largest {
            largest = distance;
        }
    }
    Ok((largest < tolerance, largest))
}

/// `SmPatchBoundaryBoundingPlanes::ComputeBoundingPlanes` for `surface`
/// over `domain`: for each side curve, a plane through three samples (or,
/// for a straight side, the plane through it across the surface normal at
/// its corner), kept when the side lies in it within 1e-4 and the whole
/// control net on one side.
pub fn boundary_planes(surface: &BSplineSurface, domain: Extent2) -> Result<BoundaryPlanes, Fault> {
    let tolerance = 0.0001;
    let mut result = BoundaryPlanes {
        planes: [None; 4],
        degenerate: [false; 4],
    };
    for side in 0..4 {
        let low_corner = side == 0 || side == 3;
        let corner = if low_corner { domain.min } else { domain.max };
        let (axis, parameter) = if side == 0 || side == 2 {
            (1, corner.y)
        } else {
            (0, corner.x)
        };
        let curve = surface.iso_curve(axis, parameter)?;
        if curve.is_degenerate(EFF_ZERO) {
            result.degenerate[side] = true;
            continue;
        }
        let interval = curve.natural_interval();
        let at = |t: f64| curve.evaluate_point(interval.evaluate(t));
        let samples = [at(0.0)?, at(0.25)?, at(0.5)?, at(0.75)?, at(1.0)?];
        let (plane, collinear) =
            plane_from_3_points(samples[0], samples[2], samples[4], tolerance)?;
        let (base, mut normal) = if collinear {
            let base = (samples[0] + samples[2] + samples[4]) / 3.0;
            let from_left = !low_corner;
            let surface_normal = evaluate_normal(surface, corner, [from_left, from_left]);
            let t = interval.evaluate(if low_corner { 0.0 } else { 1.0 });
            let d = curve.evaluate(t, 1, from_left, true)?;
            let mut direction = d[1];
            direction.unitize()?;
            (base, direction.cross(surface_normal))
        } else if let Some(plane) = plane {
            plane
        } else {
            continue;
        };
        let (planar, planarity) = curve_on_plane(&curve, base, normal, tolerance)?;
        if !planar {
            continue;
        }
        let (separating, invert) = separating_halfspace(surface, base, normal);
        if !separating {
            continue;
        }
        if invert {
            normal = -normal;
        }
        result.planes[side] = Some(BoundaryPlane {
            base,
            normal,
            tolerance: planarity,
            samples,
        });
    }
    Ok(result)
}
