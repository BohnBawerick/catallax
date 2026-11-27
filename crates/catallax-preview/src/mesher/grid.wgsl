fn grid(surface: Surface, uv: vec2<f32>) -> Vertex {
    let size = surface.mode.zw;
    let span = size - vec2(1u);
    let source_uv = surface.x.xy + uv * surface.y.xy;
    let grid_uv = select(source_uv, fract(source_uv), surface.z.xy != vec2(0.0));
    let xy = clamp(grid_uv, vec2(0.0), vec2(1.0)) * vec2<f32>(span);
    let cell = min(vec2<u32>(xy), span - vec2(1u));
    let f = xy - vec2<f32>(cell);
    let stride = size.x;
    let i = surface.mode.y + cell.y * stride + cell.x;
    let p00 = positions[i].xyz;
    let p10 = positions[i + 1u].xyz;
    let p01 = positions[i + stride].xyz;
    let p11 = positions[i + stride + 1u].xyz;
    let p = mix(mix(p00, p10, f.x), mix(p01, p11, f.x), f.y);
    let du = mix(p10 - p00, p11 - p01, f.y);
    let dv = mix(p01 - p00, p11 - p10, f.x);
    let normal = normalize(cross(du, dv));
    return Vertex(vec4(p, 1.0), vec4(normal, 0.0));
}
