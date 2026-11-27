use catallax_mesh::float::{DisplayFace, DisplayFrame, DisplaySpline, DisplaySurface};

use super::{FaceRefusal, SurfaceParams, position32, scalar32};

type Bounds = [[f64; 2]; 2];

pub(super) const GRID: u32 = 2;
const SPLINE: u32 = 3;
const MAX_SPLINE_DEGREE: usize = 12;

pub(super) fn pack<F>(
    face: &DisplayFace<F>,
    bounds: Bounds,
    detail: [u32; 2],
    positions: &mut Vec<[f32; 4]>,
    max_positions: usize,
    knots: &mut Vec<f32>,
    max_knots: usize,
) -> Result<SurfaceParams, FaceRefusal> {
    let span = bounds.map(|[lower, upper]| upper - lower);
    let facing = if face.reversed { -1.0 } else { 1.0 };
    let mut params = SurfaceParams::default();
    match &face.surface {
        DisplaySurface::Plane(frame) => {
            params.origin = position32(point(frame, bounds[0][0], bounds[1][0], 0.0))?;
            params.x = vector(frame.x.map(|x| x * span[0]), 0.0)?;
            params.y = vector(frame.y.map(|y| y * span[1]), 0.0)?;
            params.z = vector(frame.z, 0.0)?;
        }
        DisplaySurface::Cylinder { frame, radius } => {
            if !radius.is_finite() || *radius <= 0.0 {
                return Err(FaceRefusal::InvalidSurface);
            }
            params.mode[0] = 1;
            params.origin = position32(point(frame, 0.0, 0.0, bounds[1][0]))?;
            params.x = vector(frame.x, bounds[0][0].rem_euclid(1.0))?;
            params.y = vector(frame.y, span[0])?;
            params.z = vector(frame.z.map(|z| z * span[1]), *radius)?;
        }
        DisplaySurface::Grid(grid) => {
            let source_span = face.bounds.map(|[lower, upper]| upper - lower);
            let mut base = [0.0; 2];
            let mut scale = [0.0; 2];
            let mut wrap = [0.0; 2];
            for axis in 0..2 {
                base[axis] = (bounds[axis][0] - face.bounds[axis][0]) / source_span[axis];
                scale[axis] = span[axis] / source_span[axis];
                if base[axis] < -1e-10 || base[axis] + scale[axis] > 1.0 + 1e-10 {
                    if face.periods[axis]
                        .is_none_or(|p| (source_span[axis] / p - 1.0).abs() > 1e-10)
                    {
                        return Err(FaceRefusal::UnsupportedTrim);
                    }
                    wrap[axis] = 1.0;
                }
            }
            params.x = vector([base[0], base[1], 0.0], 0.0)?;
            params.y = vector([scale[0], scale[1], 0.0], 0.0)?;
            params.z = vector([wrap[0], wrap[1], 0.0], 0.0)?;
            let count = check_grid(grid.size, positions, max_positions)?;
            if grid.points.len() != count {
                return Err(FaceRefusal::InvalidSurface);
            }
            params.mode = grid_mode(positions.len(), grid.size);
            positions.reserve(count);
            for &point in &grid.points {
                positions.push(position32(point)?);
            }
        }
        DisplaySurface::BSpline(spline) => {
            params.mode[0] = SPLINE;
            params.spline =
                pack_spline(spline, bounds, positions, max_positions, knots, max_knots)?;
        }
        surface => {
            // Transitional sampler for surfaces that do not yet have a WGSL
            // evaluator. The GPU still evaluates preview vertices and normals.
            let size = detail.map(|n| (n as usize).saturating_add(1).max(2));
            let count = check_grid(size, positions, max_positions)?;
            params.mode = grid_mode(positions.len(), size);
            params.y = [1.0, 1.0, 0.0, 0.0];
            positions.reserve(count);
            sample_analytic(surface, bounds, size, positions)?;
        }
    }
    params.origin[3] = facing;
    Ok(params)
}

fn grid_mode(first: usize, size: [usize; 2]) -> [u32; 4] {
    [GRID, first as u32, size[0] as u32, size[1] as u32]
}

fn check_grid(size: [usize; 2], points: &[[f32; 4]], limit: usize) -> Result<usize, FaceRefusal> {
    if size[0] < 2 || size[1] < 2 {
        return Err(FaceRefusal::InvalidSurface);
    }
    size[0]
        .checked_mul(size[1])
        .filter(|&n| n <= limit - points.len())
        .ok_or(FaceRefusal::MeshSize)
}

fn vector(v: [f64; 3], w: f64) -> Result<[f32; 4], FaceRefusal> {
    Ok([
        scalar32(v[0])?,
        scalar32(v[1])?,
        scalar32(v[2])?,
        scalar32(w)?,
    ])
}

fn point(frame: &DisplayFrame, x: f64, y: f64, z: f64) -> [f64; 3] {
    std::array::from_fn(|i| frame.origin[i] + x * frame.x[i] + y * frame.y[i] + z * frame.z[i])
}

fn tick(bounds: Bounds, axis: usize, i: usize, count: usize) -> f64 {
    if i + 1 == count {
        bounds[axis][1]
    } else {
        bounds[axis][0] + (bounds[axis][1] - bounds[axis][0]) * (i as f64 / (count - 1) as f64)
    }
}

fn sample_analytic(
    surface: &DisplaySurface,
    bounds: Bounds,
    size: [usize; 2],
    positions: &mut Vec<[f32; 4]>,
) -> Result<(), FaceRefusal> {
    let frame = match surface {
        DisplaySurface::Cone { frame, .. }
        | DisplaySurface::Sphere { frame, .. }
        | DisplaySurface::Torus { frame, .. }
        | DisplaySurface::EllipticCylinder { frame, .. } => frame,
        _ => unreachable!("only unsupported analytic carriers are sampled here"),
    };
    let azimuths: Vec<_> = (0..size[0])
        .map(|u| (tick(bounds, 0, u, size[0]).rem_euclid(1.0) * std::f64::consts::TAU).sin_cos())
        .collect();
    for v in 0..size[1] {
        let v = tick(bounds, 1, v, size[1]);
        let [rx, ry, z] = match surface {
            DisplaySurface::Cone { slope, .. } => [v * slope, v * slope, v],
            DisplaySurface::Sphere { radius, .. } => {
                let (sin, cos) = (v.rem_euclid(1.0) * std::f64::consts::TAU).sin_cos();
                [radius * cos, radius * cos, radius * sin]
            }
            DisplaySurface::Torus { major, minor, .. } => {
                let (sin, cos) = (v.rem_euclid(1.0) * std::f64::consts::TAU).sin_cos();
                let r = major + minor * cos;
                [r, r, minor * sin]
            }
            DisplaySurface::EllipticCylinder { radii: [x, y], .. } => [*x, *y, v],
            _ => unreachable!(),
        };
        for &(sin, cos) in &azimuths {
            positions.push(position32(point(frame, rx * cos, ry * sin, z))?);
        }
    }
    Ok(())
}

fn pack_spline(
    spline: &DisplaySpline,
    bounds: Bounds,
    controls: &mut Vec<[f32; 4]>,
    max_controls: usize,
    knots: &mut Vec<f32>,
    max_knots: usize,
) -> Result<[u32; 8], FaceRefusal> {
    let [nu, nv] = spline.size;
    let [du, dv] = spline.degree;
    let count = nu.checked_mul(nv).ok_or(FaceRefusal::InvalidSurface)?;
    if spline.points.len() != count || spline.weights.as_ref().is_some_and(|w| w.len() != count) {
        return Err(FaceRefusal::InvalidSurface);
    }
    if count > max_controls - controls.len() {
        return Err(FaceRefusal::MeshSize);
    }
    let control_offset = controls.len() as u32;
    let [ku, kv] = pack_knots(spline, bounds, knots, max_knots)?;
    controls.reserve(count);
    // The description is U-fastest, while spline_sample reads V-fastest.
    for u in 0..nu {
        for v in 0..nv {
            let i = v * nu + u;
            let w = spline.weights.as_ref().map_or(1.0, |weights| weights[i]);
            if !w.is_finite() || w <= 0.0 {
                return Err(FaceRefusal::InvalidSurface);
            }
            let weight = scalar32(w)?;
            if weight == 0.0 {
                return Err(FaceRefusal::GpuRange);
            }
            let [x, y, z] = spline.points[i].map(|p| scalar32(p * w));
            controls.push([x?, y?, z?, weight]);
        }
    }
    Ok([
        control_offset,
        ku,
        kv,
        du as u32,
        dv as u32,
        nu as u32,
        nv as u32,
        0,
    ])
}

fn pack_knots(
    spline: &DisplaySpline,
    bounds: Bounds,
    packed: &mut Vec<f32>,
    limit: usize,
) -> Result<[u32; 2], FaceRefusal> {
    let mut offsets = [0; 2];
    for (axis, &[lower, upper]) in bounds.iter().enumerate() {
        let knots = &spline.knots[axis];
        let degree = spline.degree[axis];
        let n = spline.size[axis];
        if degree > MAX_SPLINE_DEGREE
            || degree >= n
            || n.checked_add(degree + 1) != Some(knots.len())
            || knots[degree] >= knots[n]
        {
            return Err(FaceRefusal::InvalidSurface);
        }
        let lo = knots[degree];
        let hi = knots[n];
        let tolerance = (hi - lo) * 1e-10;
        if lower < lo - tolerance || upper > hi + tolerance {
            return Err(FaceRefusal::UnsupportedTrim);
        }
        if knots.len() > limit - packed.len() {
            return Err(FaceRefusal::MeshSize);
        }
        offsets[axis] = packed.len() as u32;
        packed.reserve(knots.len());
        // Match the trim's normalized UVs before converting to binary32.
        let mut previous = f64::NEG_INFINITY;
        let mut previous32 = f32::NEG_INFINITY;
        let span = upper - lower;
        for &knot in knots {
            if !knot.is_finite() || knot < previous {
                return Err(FaceRefusal::InvalidSurface);
            }
            let converted = scalar32((knot - lower) / span)?;
            if knot > previous && converted <= previous32 {
                return Err(FaceRefusal::GpuRange);
            }
            packed.push(converted);
            previous = knot;
            previous32 = converted;
        }
    }
    Ok(offsets)
}
