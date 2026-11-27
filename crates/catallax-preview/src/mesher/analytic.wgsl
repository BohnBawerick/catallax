fn plane(surface: Surface, uv: vec2<f32>) -> Vertex {
    let p = surface.origin.xyz + uv.x * surface.x.xyz + uv.y * surface.y.xyz;
    return Vertex(vec4(p, 1.0), vec4(surface.z.xyz, 0.0));
}

fn cylinder(surface: Surface, uv: vec2<f32>) -> Vertex {
    let ray = analytic_circle(surface.x.w + uv.x * surface.y.w);
    let x = ray.x * surface.x.xyz;
    let y = ray.y * surface.y.xyz;
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

fn analytic_face(surface: Surface, uv: vec2<f32>) -> Vertex {
    let base = bitcast<vec2<f32>>(surface.mode.yz);
    let span = vec2(bitcast<f32>(surface.mode.w), surface.z.w);
    let data = array(
        surface.origin,
        surface.z,
        surface.x,
        surface.y,
        vec4(0.0),
        vec4(surface.x.w, surface.y.w, 0.0, 0.0),
        vec4(0.0),
    );
    let result = eval_analytic(
        AnalyticSurface(surface.mode.x - 4u, 1.0, vec2(0u), data),
        base + uv * span,
    );
    return Vertex(result.position, result.normal);
}

fn stereographic_sphere(surface: Surface, uv: vec2<f32>) -> Vertex {
    let base = bitcast<vec2<f32>>(surface.mode.yz);
    let st = base + uv * vec2(surface.z.w, surface.y.w);
    let r2 = dot(st, st);
    let n = (2.0 * (st.x * surface.x.xyz + st.y * surface.y.xyz)
        + (1.0 - r2) * surface.z.xyz) / (1.0 + r2);
    return Vertex(vec4(surface.origin.xyz + surface.x.w * n, 1.0), vec4(n, 0.0));
}
