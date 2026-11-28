use std::cmp::Ordering;

use catallax_numbers::Rational;
use catallax_shapes::PatchSectionBranch;

use super::{DisplaySurface, Xyz, edges::sample_span, evaluate::*, spline};

pub(super) fn branch(
    branch: &PatchSectionBranch,
    ends: [Option<Xyz>; 2],
    origin: Xyz,
    tolerance: f64,
    max: usize,
) -> Option<Vec<Xyz>> {
    let source = branch.surface().surface();
    let cutting = describe(&source, origin)?;
    let nets = branch
        .charts()
        .iter()
        .map(|c| spline::net(c.carrier(), origin))
        .collect::<Option<Vec<_>>>()?;
    let mut patches = nets
        .iter()
        .map(spline::Patch::new)
        .collect::<Option<Vec<_>>>()?;
    let mut points = Vec::new();
    for (chart, solved, bounds, parameters) in branch.display_cells() {
        let patch = &mut patches[chart];
        let ends = parameters.each_ref().map(|uv| uv.each_ref().map(scalar));
        if ends.iter().flatten().any(|t| !t.is_finite()) {
            return None;
        }
        let free = 1 - solved;
        let interval = ends.map(|uv| uv[free]);
        let range = bounds[solved].each_ref().map(Rational::to_f64_nearest);
        let start = points.last().copied().unwrap_or_else(|| patch.at(ends[0]));
        let end = patch.at(ends[1]);
        if points.is_empty() {
            points.push(start);
        }
        let mut at = |t| {
            root_point(patch, &cutting, solved, t, range, tolerance * 0.001)
                .unwrap_or([f64::NAN; 3])
        };
        sample_span(&mut at, interval, [start, end], tolerance, &mut |_, p| {
            (points.len() < max).then(|| points.push(p))
        })?;
    }
    if branch.is_closed() {
        let first = *points.first()?;
        *points.last_mut()? = first;
    }
    clip(points, ends, branch.is_closed(), tolerance)
}

/// A retained fillet cap cut is a single graph over its across-fillet parameter.
pub(super) fn fillet_cap(
    spline: &super::DisplaySpline,
    plane: &DisplaySurface,
    ends: [Xyz; 2],
    tolerance: f64,
    max: usize,
) -> Option<Vec<Xyz>> {
    let mut patch = spline::Patch::new(spline)?;
    let interval = ends.map(|p| patch.project(p, None, [None; 2])[1]);
    if interval[0].partial_cmp(&interval[1]) == Some(Ordering::Equal) {
        return None;
    }
    let range = patch.bounds[0];
    let mut at =
        |v| root_point(&mut patch, plane, 0, v, range, tolerance * 0.001).unwrap_or([f64::NAN; 3]);
    let mut points = vec![ends[0]];
    sample_span(&mut at, interval, ends, tolerance, &mut |_, p| {
        (points.len() < max).then(|| points.push(p))
    })?;
    Some(points)
}

fn root_point(
    patch: &mut spline::Patch<'_>,
    cutting: &DisplaySurface,
    solved: usize,
    t: f64,
    [mut lo, mut hi]: [f64; 2],
    accuracy: f64,
) -> Option<Xyz> {
    let mut at = |s| {
        let uv = if solved == 0 { [s, t] } else { [t, s] };
        let p = patch.at(uv);
        Some((p, signed_distance(cutting, p)?))
    };
    let (a, fa) = at(lo)?;
    let (b, fb) = at(hi)?;
    if fa.abs() <= accuracy {
        return Some(a);
    }
    if fb.abs() <= accuracy {
        return Some(b);
    }
    if fa.is_sign_negative() == fb.is_sign_negative() {
        return None;
    }
    for _ in 0..52 {
        let mid = lo.midpoint(hi);
        let (p, f) = at(mid)?;
        if f.abs() <= accuracy
            || [lo, hi]
                .iter()
                .any(|end| mid.partial_cmp(end) == Some(Ordering::Equal))
        {
            return Some(p);
        }
        if f.is_sign_negative() == fa.is_sign_negative() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(at(lo.midpoint(hi))?.0)
}

fn signed_distance(surface: &DisplaySurface, p: Xyz) -> Option<f64> {
    use DisplaySurface as D;
    let delta = |f: &super::DisplayFrame| {
        let d = sub(p, f.origin);
        let z = dot(d, f.z);
        [norm(sub(d, scale(f.z, z))), z]
    };
    finite(match surface {
        D::Plane(f) => dot(sub(p, f.origin), f.z),
        D::Cylinder { frame, radius } => delta(frame)[0] - radius,
        D::Cone { frame, slope } => {
            let [radial, z] = delta(frame);
            radial - slope * z
        }
        D::Sphere { frame, radius } => norm(sub(p, frame.origin)) - radius,
        D::Torus {
            frame,
            major,
            minor,
        } => {
            let [radial, z] = delta(frame);
            libm::hypot(radial - major, z) - minor
        }
        _ => return None,
    })
}

fn clip(
    mut points: Vec<Xyz>,
    ends: [Option<Xyz>; 2],
    closed: bool,
    tolerance: f64,
) -> Option<Vec<Xyz>> {
    let [start, end] = match ends {
        [None, None] if closed => return Some(points),
        [Some(a), Some(b)] => [a, b],
        _ => return None,
    };
    let nearest = |p| {
        points
            .windows(2)
            .enumerate()
            .filter_map(|(i, pair)| {
                let d = sub(pair[1], pair[0]);
                let length = positive(dot(d, d))?;
                let t = (dot(sub(p, pair[0]), d) / length).clamp(0.0, 1.0);
                let error = norm(sub(p, add(pair[0], scale(d, t))));
                Some((i as f64 + t, error))
            })
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .filter(|(_, e)| *e <= tolerance * 2.0)
    };
    let (a, _) = nearest(start)?;
    let (mut b, _) = nearest(end)?;
    let count = points.len() - 1;
    if closed && b <= a {
        b += count as f64;
    }
    if b < a {
        return None;
    }
    let lo = libm::floor(a) as usize + 1;
    let hi = libm::ceil(b) as usize;
    if hi <= count {
        points.drain(hi..);
        points.drain(..lo.min(points.len()));
        points.insert(0, start);
        points.push(end);
        return Some(points);
    }
    // Wrapped traversal omits the duplicate end of the complete branch.
    let mut arc = Vec::with_capacity(hi - lo + 2);
    arc.push(start);
    arc.extend((lo..hi).map(|i| points[i % count]));
    arc.push(end);
    Some(arc)
}
