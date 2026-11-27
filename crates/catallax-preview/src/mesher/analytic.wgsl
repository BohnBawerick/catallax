fn plane(surface: Surface, uv: vec2<f32>) -> Vertex {
    let p = surface.origin.xyz + uv.x * surface.x.xyz + uv.y * surface.y.xyz;
    return Vertex(vec4(p, 1.0), vec4(surface.z.xyz, 0.0));
}

fn cylinder(surface: Surface, uv: vec2<f32>) -> Vertex {
    let angle = fract(surface.x.w + uv.x * surface.y.w) * 6.283185307179586;
    let x = cos(angle) * surface.x.xyz;
    let y = sin(angle) * surface.y.xyz;
    let rx = surface.z.w;
    let ry = bitcast<f32>(surface.mode.y);
    let base = surface.origin.xyz + uv.y * surface.z.xyz;
    if ry == 0.0 {
        let radial = x + y;
        return Vertex(vec4(base + rx * radial, 1.0), vec4(radial, 0.0));
    }
    let scale = max(rx, ry);
    let normal = normalize((ry / scale) * x + (rx / scale) * y);
    return Vertex(vec4(base + rx * x + ry * y, 1.0), vec4(normal, 0.0));
}
