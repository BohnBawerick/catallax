use catallax_mesh::float::{DisplayFace, DisplayFrame, DisplaySpline, DisplaySurface};

use crate::eval::analytic::Kind;

use super::{FaceRefusal, SurfaceParams, position32, scalar32};

type Bounds = [[f64; 2]; 2];

pub(super) const GRID: u32 = 2;
const SPLINE: u32 = 3;
const ANALYTIC: u32 = 4;
const MAX_SPLINE_DEGREE: usize = 12;

pub(super) fn pack<F>(
    face: &DisplayFace<F>,
    bounds: Bounds,
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
            params = cylinder(frame, *radius, None, bounds, span)?;
        }
        DisplaySurface::EllipticCylinder { frame, radii } => {
            params = cylinder(frame, radii[0], Some(radii[1]), bounds, span)?;
        }
        DisplaySurface::Cone { frame, slope } => {
            params = analytic(frame, Kind::Cone, *slope, None, bounds)?;
        }
        DisplaySurface::Sphere { frame, radius } => {
            params = analytic(frame, Kind::Sphere, *radius, None, bounds)?;
        }
        DisplaySurface::Torus {
            frame,
            major,
            minor,
        } => {
            params = analytic(frame, Kind::Torus, *major, Some(*minor), bounds)?;
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
            let offsets = pack_knots(spline, bounds, face.periods, knots, max_knots, &mut params)?;
            params.spline = pack_spline(spline, positions, max_positions, offsets)?;
        }
    }
    params.origin[3] = facing;
    Ok(params)
}

fn cylinder(
    frame: &DisplayFrame,
    radius: f64,
    second_radius: Option<f64>,
    bounds: Bounds,
    span: [f64; 2],
) -> Result<SurfaceParams, FaceRefusal> {
    let rx = radius32(radius)?;
    let ry = second_radius.map(radius32).transpose()?.unwrap_or(0.0);
    Ok(SurfaceParams {
        origin: position32(point(frame, 0.0, 0.0, bounds[1][0]))?,
        x: vector(frame.x, bounds[0][0].rem_euclid(1.0))?,
        y: vector(frame.y, span[0])?,
        z: vector(frame.z.map(|z| z * span[1]), f64::from(rx))?,
        mode: [1, ry.to_bits(), 0, 0],
        ..Default::default()
    })
}

fn analytic(
    frame: &DisplayFrame,
    kind: Kind,
    first: f64,
    second: Option<f64>,
    [[u0, u1], [v0, v1]]: Bounds,
) -> Result<SurfaceParams, FaceRefusal> {
    let k0 = radius32(first)?;
    let k1 = second.map(radius32).transpose()?.unwrap_or(0.0);
    let v_base = if kind == Kind::Torus {
        v0.rem_euclid(1.0)
    } else {
        v0
    };
    Ok(SurfaceParams {
        origin: position32(frame.origin)?,
        x: vector(frame.x, f64::from(k0))?,
        y: vector(frame.y, f64::from(k1))?,
        z: vector(frame.z, v1 - v0)?,
        mode: [
            ANALYTIC + kind as u32,
            scalar32(u0.rem_euclid(1.0))?.to_bits(),
            scalar32(v_base)?.to_bits(),
            scalar32(u1 - u0)?.to_bits(),
        ],
        ..Default::default()
    })
}

fn radius32(radius: f64) -> Result<f32, FaceRefusal> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(FaceRefusal::InvalidSurface);
    }
    let radius = scalar32(radius)?;
    if radius == 0.0 {
        return Err(FaceRefusal::GpuRange);
    }
    Ok(radius)
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

fn pack_spline(
    spline: &DisplaySpline,
    controls: &mut Vec<[f32; 4]>,
    max_controls: usize,
    [ku, kv]: [u32; 2],
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
    periods: [Option<f64>; 2],
    packed: &mut Vec<f32>,
    limit: usize,
    params: &mut SurfaceParams,
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
        let span = hi - lo;
        let tolerance = span * 1e-10;
        let (knot_base, knot_span) = if lower < lo - tolerance || upper > hi + tolerance {
            if periods[axis].is_none_or(|p| (span - p).abs() > tolerance) {
                return Err(FaceRefusal::UnsupportedTrim);
            }
            params.x[axis] = scalar32(((lower - lo) / span).rem_euclid(1.0))?;
            params.y[axis] = scalar32((upper - lower) / span)?;
            params.z[axis] = 1.0;
            (lo, span)
        } else {
            params.y[axis] = 1.0;
            (lower, upper - lower)
        };
        if knots.len() > limit - packed.len() {
            return Err(FaceRefusal::MeshSize);
        }
        offsets[axis] = packed.len() as u32;
        packed.reserve(knots.len());
        // Keep knots in trim coordinates unless the carrier must wrap.
        let mut previous = f64::NEG_INFINITY;
        let mut previous32 = f32::NEG_INFINITY;
        for &knot in knots {
            if !knot.is_finite() || knot < previous {
                return Err(FaceRefusal::InvalidSurface);
            }
            let converted = scalar32((knot - knot_base) / knot_span)?;
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
