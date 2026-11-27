//! Standalone GPU surface evaluation; mesher integration supplies the body origin and trims.
#![allow(dead_code)] // The mesher consumes this module in a separate lane.

#[cfg(test)]
use catallax_shapes::{Point3, SplineSupport, Surface, Vector3};

pub const WGSL: &str = include_str!("analytic.wgsl");

#[cfg(test)]
type D3 = [f64; 3];

#[repr(u32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Cone,
    Sphere,
    Torus,
    Sweep,
    Pipe,
    Ruled,
    LinearFillet,
}

/// Binary64, body-local inputs before upload. The chart orientation is supplied
/// by the caller for analytic supports and differential charts; for pipes it
/// is relative to the increasing-u/v chart, not the ball's outward normal.
pub struct Definition {
    pub kind: Kind,
    pub orientation: f64,
    pub data: [[f64; 4]; 7],
}

/// Storage layout shared with `AnalyticSurface` in WGSL.
#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub struct AnalyticSurface {
    kind: u32,
    orientation: f32,
    reserved: [u32; 2],
    data: [[f32; 4]; 7],
}

#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct AnalyticVertex {
    pub position: [f32; 4],
    pub normal: [f32; 4],
}

impl Definition {
    /// Reads stored support data, without sampling the surface or calling a proof route.
    /// `same` selects the natural or reversed face sense. Fitted or general spline
    /// supports are left to the separate control-net evaluator.
    #[cfg(test)]
    pub fn from_surface(surface: &Surface, body_origin: D3, same: bool) -> Option<Self> {
        let orientation = if same { 1.0 } else { -1.0 };
        let mut data = [[0.0; 4]; 7];
        let kind = match surface {
            Surface::Cone(s) => {
                frame(&mut data, &s.apex, &s.axis, &s.ref_dir, body_origin);
                data[5][0] = s.slope.to_f64();
                Kind::Cone
            }
            Surface::Sphere(s) => {
                frame(&mut data, &s.center, &s.axis, &s.ref_dir, body_origin);
                data[5][0] = s.radius.to_f64();
                Kind::Sphere
            }
            Surface::Torus(s) => {
                frame(&mut data, &s.center, &s.axis, &s.ref_dir, body_origin);
                data[5][..2].copy_from_slice(&[s.major.to_f64(), s.minor.to_f64()]);
                Kind::Torus
            }
            Surface::TorusSheet(s) => {
                frame(&mut data, &s.center, &s.axis, &s.ref_dir, body_origin);
                data[5][..2].copy_from_slice(&[s.major.to_f64(), s.minor.to_f64()]);
                Kind::Torus
            }
            Surface::BSpline(SplineSupport::BezierExtrusion(s)) => {
                for (row, control) in data.iter_mut().zip(s.controls()) {
                    *row = lane(local(control, body_origin));
                }
                data[4] = lane(vector(s.direction()));
                data[5][0] = (s.controls().len() - 1) as f64;
                Kind::Sweep
            }
            Surface::BSpline(SplineSupport::RollingPipe(s)) => {
                let law = s.law();
                let axis = vector(&law.axis);
                let across = vector(&law.across);
                let side = vector(&law.side);
                let reach = law.reach.to_f64();
                let slope = law.slope.to_f64();
                let wall = law.wall.to_f64() / law.scale.to_f64();
                let first = vector(&law.normal).map(|x| x * reach);
                let last = std::array::from_fn(|i| (across[i] - slope * axis[i]) * wall);
                let chart = dot(cross(side, last), first).signum();
                data[0] = lane(local(&law.base, body_origin));
                data[1] = lane(axis);
                data[2] = lane(across);
                data[3] = lane(side);
                data[4] = lane(first);
                data[5] = [
                    law.offset.to_f64(),
                    law.lift.to_f64(),
                    slope,
                    slope * law.tilt.to_f64() / law.along.to_f64(),
                ];
                data[6] = [wall, law.radius.to_f64(), 0.0, 0.0];
                return Some(Self {
                    kind: Kind::Pipe,
                    orientation: orientation * chart,
                    data,
                });
            }
            Surface::BSpline(SplineSupport::ObliqueDraft(s)) => {
                data[0] = lane(local(s.origin(), body_origin));
                data[1] = lane(vector(s.axis()));
                data[2] = lane(vector(s.across()));
                data[3] = lane(vector(s.side()));
                data[4][..2].copy_from_slice(&s.lean().each_ref().map(|x| x.to_f64()));
                let lo = s.heights()[0].to_f64();
                let hi = s.heights()[1].to_f64();
                let slope = s.tangent().to_f64();
                data[5] = [
                    lo,
                    hi,
                    s.amplitude(&catallax_numbers::Real::zero()).to_f64(),
                    slope,
                ];
                Kind::Ruled
            }
            Surface::BSpline(SplineSupport::LinearFillet(s)) => {
                let [first, last] = s.rays().each_ref().map(vector);
                data[0] = lane(local(&s.edge()[0], body_origin));
                data[1] = lane(vector(s.along()));
                let shift = vector(s.shift());
                data[2] = lane(std::array::from_fn(|i| shift[i] + first[i]));
                data[3] = lane(std::array::from_fn(|i| shift[i] + last[i]));
                let [lo, hi] = s.radii().each_ref().map(|r| r.to_f64());
                data[5] = [lo, hi, ((1.0 + dot(first, last)) * 0.5).sqrt(), 0.0];
                Kind::LinearFillet
            }
            _ => return None,
        };
        Some(Self {
            kind,
            orientation,
            data,
        })
    }

    pub fn upload(&self) -> AnalyticSurface {
        AnalyticSurface {
            kind: self.kind as u32,
            orientation: self.orientation as f32,
            reserved: [0; 2],
            data: self.data.map(|row| row.map(|x| x as f32)),
        }
    }
}

#[cfg(test)]
fn frame(data: &mut [[f64; 4]; 7], p: &Point3, a: &Vector3, r: &Vector3, origin: D3) {
    let axis = unit(vector(a));
    let across = unit(vector(r));
    data[0] = lane(local(p, origin));
    data[1] = lane(axis);
    data[2] = lane(across);
    data[3] = lane(unit(cross(axis, across)));
}

#[cfg(test)]
fn local(p: &Point3, origin: D3) -> D3 {
    let point = [&p.x, &p.y, &p.z].map(|v| v.to_f64());
    std::array::from_fn(|i| point[i] - origin[i])
}

#[cfg(test)]
fn vector(v: &Vector3) -> D3 {
    [&v.x, &v.y, &v.z].map(|x| x.to_f64())
}

#[cfg(test)]
fn lane([x, y, z]: D3) -> [f64; 4] {
    [x, y, z, 0.0]
}

#[cfg(test)]
fn dot(a: D3, b: D3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[cfg(test)]
fn cross([x, y, z]: D3, [a, b, c]: D3) -> D3 {
    [y * c - z * b, z * a - x * c, x * b - y * a]
}

#[cfg(test)]
fn unit(v: D3) -> D3 {
    let length = v[0].hypot(v[1]).hypot(v[2]);
    v.map(|x| x / length)
}
