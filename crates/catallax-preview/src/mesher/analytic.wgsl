fn plane(surface: Surface, uv: vec2<f32>) -> Vertex {
    let p = surface.origin.xyz + uv.x * surface.x.xyz + uv.y * surface.y.xyz;
    return Vertex(vec4(p, 1.0), vec4(surface.z.xyz, 0.0));
}

fn cylinder(surface: Surface, uv: vec2<f32>) -> Vertex {
    let angle = fract(surface.x.w + uv.x * surface.y.w) * 6.283185307179586;
    let radial = cos(angle) * surface.x.xyz + sin(angle) * surface.y.xyz;
    let p = surface.origin.xyz + surface.z.w * radial + uv.y * surface.z.xyz;
    return Vertex(vec4(p, 1.0), vec4(radial, 0.0));
}
