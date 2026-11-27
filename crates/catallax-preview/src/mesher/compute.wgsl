struct InputVertex {
    uv: vec2<f32>,
    boundary: u32,
    surface: u32,
}

struct Vertex {
    position: vec4<f32>,
    normal: vec4<f32>,
}

struct Surface {
    origin: vec4<f32>,
    x: vec4<f32>,
    y: vec4<f32>,
    z: vec4<f32>,
    mode: vec4<u32>,
}

@group(0) @binding(0) var<storage, read> inputs: array<InputVertex>;
@group(0) @binding(1) var<storage, read> positions: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> surfaces: array<Surface>;
@group(0) @binding(3) var<storage, read_write> vertices: array<Vertex>;

@compute @workgroup_size(64)
fn mesh(@builtin(global_invocation_id) id: vec3<u32>, @builtin(num_workgroups) groups: vec3<u32>) {
    let i = id.x + id.y * (groups.x * 64u);
    if i >= arrayLength(&inputs) {
        return;
    }
    let input = inputs[i];
    let surface = surfaces[input.surface];
    var vertex: Vertex;
    switch surface.mode.x {
        case 0u: { vertex = plane(surface, input.uv); }
        case 1u: { vertex = cylinder(surface, input.uv); }
        default: { vertex = grid(surface, input.uv); }
    }
    if input.boundary != 0xffffffffu {
        vertex.position = positions[input.boundary];
    }
    vertex.normal = vec4(vertex.normal.xyz * surface.origin.w, 0.0);
    vertices[i] = vertex;
}
