use catallax_mesh::float::{DisplayFace, DisplayFrame, DisplaySpline, DisplaySurface};

use super::{FaceRefusal, SurfaceParams, position32, scalar32};

type Bounds = [[f64; 2]; 2];

pub(super) fn pack<F>(
    face: &DisplayFace<F>,
    bounds: Bounds,
    detail: [u32; 2],
    positions: &mut Vec<[f32; 4]>,
    max_positions: usize,
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
        surface => {
            // Transitional sampler for surfaces that do not yet have a WGSL
            // evaluator. The GPU still evaluates preview vertices and normals.
            let size = detail.map(|n| (n as usize).saturating_add(1).max(2));
            let count = check_grid(size, positions, max_positions)?;
            params.mode = grid_mode(positions.len(), size);
            params.y = [1.0, 1.0, 0.0, 0.0];
            positions.reserve(count);
            match surface {
                DisplaySurface::BSpline(spline) => sample_spline(spline, bounds, size, positions)?,
                _ => sample_analytic(surface, bounds, size, positions)?,
            }
        }
    }
    params.origin[3] = facing;
    Ok(params)
}

fn grid_mode(first: usize, size: [usize; 2]) -> [u32; 4] {
    [2, first as u32, size[0] as u32, size[1] as u32]
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
    let (frame, radii, axial): (_, [f64; 2], fn([f64; 2], f64) -> [f64; 3]) = match surface {
        DisplaySurface::Cone { frame, slope } => (frame, [*slope, 0.0], |[slope, _], v| {
            [v * slope, v * slope, v]
        }),
        DisplaySurface::Sphere { frame, radius } => (frame, [*radius, 0.0], |[radius, _], v| {
            let (s, c) = (v.rem_euclid(1.0) * std::f64::consts::TAU).sin_cos();
            [radius * c, radius * c, radius * s]
        }),
        DisplaySurface::Torus {
            frame,
            major,
            minor,
        } => (frame, [*major, *minor], |[major, minor], v| {
            let (s, c) = (v.rem_euclid(1.0) * std::f64::consts::TAU).sin_cos();
            [major + minor * c, major + minor * c, minor * s]
        }),
        DisplaySurface::EllipticCylinder { frame, radii } => (frame, *radii, |[x, y], v| [x, y, v]),
        _ => unreachable!("only unsupported analytic carriers are sampled here"),
    };
    let azimuths: Vec<_> = (0..size[0])
        .map(|u| (tick(bounds, 0, u, size[0]).rem_euclid(1.0) * std::f64::consts::TAU).sin_cos())
        .collect();
    for v in 0..size[1] {
        let [rx, ry, z] = axial(radii, tick(bounds, 1, v, size[1]));
        for &(sin, cos) in &azimuths {
            positions.push(position32(point(frame, rx * cos, ry * sin, z))?);
        }
    }
    Ok(())
}

fn sample_spline(
    spline: &DisplaySpline,
    bounds: Bounds,
    size: [usize; 2],
    positions: &mut Vec<[f32; 4]>,
) -> Result<(), FaceRefusal> {
    let count = spline.size[0]
        .checked_mul(spline.size[1])
        .ok_or(FaceRefusal::InvalidSurface)?;
    if spline.points.len() != count
        || spline
            .weights
            .as_ref()
            .is_some_and(|w| w.len() != count || w.iter().any(|w| !w.is_finite() || *w <= 0.0))
    {
        return Err(FaceRefusal::InvalidSurface);
    }
    let u = Basis::sample(spline, 0, bounds, size[0])?;
    let v = Basis::sample(spline, 1, bounds, size[1])?;
    for (first_v, bv) in v.rows() {
        for (first_u, bu) in u.rows() {
            let mut p = [0.0; 3];
            let mut sum = 0.0;
            for (j, &b) in bv.iter().enumerate() {
                let first = (first_v + j) * spline.size[0] + first_u;
                for (i, &a) in bu.iter().enumerate() {
                    let weight = a
                        * b
                        * spline
                            .weights
                            .as_ref()
                            .map_or(1.0, |weights| weights[first + i]);
                    for k in 0..3 {
                        p[k] += weight * spline.points[first + i][k];
                    }
                    sum += weight;
                }
            }
            positions.push(position32(p.map(|p| p / sum))?);
        }
    }
    Ok(())
}

struct Basis {
    starts: Vec<usize>,
    values: Vec<f64>,
    stride: usize,
}

impl Basis {
    fn sample(
        spline: &DisplaySpline,
        axis: usize,
        bounds: Bounds,
        count: usize,
    ) -> Result<Self, FaceRefusal> {
        let knots = &spline.knots[axis];
        let degree = spline.degree[axis];
        let n = spline.size[axis];
        if degree >= n
            || n.checked_add(degree + 1) != Some(knots.len())
            || knots.iter().any(|k| !k.is_finite())
            || knots.windows(2).any(|k| k[0] > k[1])
            || knots[degree] >= knots[n]
        {
            return Err(FaceRefusal::InvalidSurface);
        }
        let (lo, hi) = (knots[degree], knots[n]);
        let tolerance = (hi - lo) * 1e-10;
        let period = hi - lo;
        let periodic = bounds[0][axis] >= lo - tolerance && bounds[1][axis] <= hi + tolerance;
        if !periodic {
            return Err(FaceRefusal::UnsupportedTrim);
        }
        let stride = degree + 1;
        let len = count.checked_mul(stride).ok_or(FaceRefusal::MeshSize)?;
        let mut basis = Self {
            starts: Vec::with_capacity(count),
            values: vec![0.0; len],
            stride,
        };
        let mut left = vec![0.0; stride];
        let mut right = vec![0.0; stride];
        for (i, values) in basis.values.chunks_exact_mut(stride).enumerate() {
            let t = tick(bounds, axis, i, count).clamp(lo, lo + period);
            let span = knots
                .partition_point(|&k| k <= t)
                .saturating_sub(1)
                .clamp(degree, n - 1);
            basis.starts.push(span - degree);
            values[0] = 1.0;
            for j in 1..=degree {
                left[j] = t - knots[span + 1 - j];
                right[j] = knots[span + j] - t;
                let mut saved = 0.0;
                for r in 0..j {
                    let denom = right[r + 1] + left[j - r];
                    let term = if denom > 0.0 { values[r] / denom } else { 0.0 };
                    values[r] = saved + right[r + 1] * term;
                    saved = left[j - r] * term;
                }
                values[j] = saved;
            }
        }
        Ok(basis)
    }

    fn rows(&self) -> impl Iterator<Item = (usize, &[f64])> {
        self.starts
            .iter()
            .copied()
            .zip(self.values.chunks_exact(self.stride))
    }
}
