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
    spline: SplineNet,
}

@group(0) @binding(0) var<storage, read> inputs: array<InputVertex>;
// Shared by boundary/grid positions and homogeneous spline controls.
@group(0) @binding(1) var<storage, read> spline_controls: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read> surfaces: array<Surface>;
@group(0) @binding(3) var<storage, read_write> vertices: array<Vertex>;
@group(0) @binding(4) var<storage, read> spline_knots: array<f32>;

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
        case 3u: {
            let uv = surface.x.xy + input.uv * surface.y.xy;
            let wrapped = select(uv, fract(uv), surface.z.xy != vec2(0.0));
            let sample = spline_sample(surface.spline, wrapped);
            var normal = sample.normal;
            if all(normal == vec3(0.0)) {
                // A collapsed boundary tangent uses its limit from the chart interior.
                let near = mix(input.uv, vec2(0.5), 0.0001) * surface.y.xy + surface.x.xy;
                let inside = select(near, fract(near), surface.z.xy != vec2(0.0));
                normal = spline_sample(surface.spline, inside).normal;
            }
            vertex = Vertex(vec4(sample.position, 1.0), vec4(normal, 0.0));
        }
        case 4u, 5u, 6u: { vertex = analytic_face(surface, input.uv); }
        case 7u: { vertex = stereographic_sphere(surface, input.uv); }
        default: { vertex = grid(surface, input.uv); }
    }
    if input.boundary != 0xffffffffu {
        vertex.position = spline_controls[input.boundary];
    }
    vertex.normal = vec4(vertex.normal.xyz * surface.origin.w, 0.0);
    vertices[i] = vertex;
}
