// Unproved preview batch kernels. Each invocation owns one output value;
// there is no reduction, atomic or cross-invocation ordering.
struct Params {
    origin: vec4<i32>,
    counts: vec4<u32>,
    start: u32,
    len: u32,
    pitch: f32,
    pad: u32,
}
@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read_write> values: array<f32>;
@group(0) @binding(2) var<storage, read> points: array<vec4<f32>>;

@compute @workgroup_size(64)
fn evaluate_points(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= params.len) { return; }
    values[id.x] = field(points[id.x].xyz);
}

// Grid point (i, j, k) is (origin + (i, j, k) + 1/2) * pitch. The host keeps
// every index below 2^23 in magnitude and the pitch a power of two, so the
// point is exact in binary32 and equal to the processor's binary64 point.
@compute @workgroup_size(64)
fn evaluate_grid(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= params.len) { return; }
    let index = params.start + id.x;
    let nx = params.counts.x;
    let ny = params.counts.y;
    let n = vec3<i32>(i32(index % nx), i32((index / nx) % ny), i32(index / (nx * ny)));
    values[id.x] = field((vec3<f32>(params.origin.xyz + n) + vec3<f32>(0.5)) * params.pitch);
}
