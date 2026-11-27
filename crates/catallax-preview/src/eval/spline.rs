use bytemuck::{Pod, Zeroable};

pub(crate) const WGSL: &str = include_str!("spline.wgsl");

/// Layout of `SplineNet` in WGSL. Offsets count elements, not bytes. The
/// control array contains `[w*x, w*y, w*z, w]` in u-major, v-minor order;
/// the knot array contains both expanded knot vectors. Degrees are at most 8
/// per axis. The converted knot vectors must preserve every nonzero span.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub(crate) struct Net {
    pub control_offset: u32,
    pub knot_u_offset: u32,
    pub knot_v_offset: u32,
    pub degree_u: u32,
    pub degree_v: u32,
    pub count_u: u32,
    pub count_v: u32,
}
