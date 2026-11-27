// Display-only evaluators. Positions are relative to the body's upload origin.
// Circular charts use turns; pipe and draft use the kernel's rational quarters.
struct AnalyticSurface {
    kind: u32,
    orientation: f32,
    reserved: vec2<u32>,
    data: array<vec4<f32>, 7>,
}

struct AnalyticVertex {
    position: vec4<f32>,
    normal: vec4<f32>,
}

fn analytic_vertex(p: vec3<f32>, n: vec3<f32>, orientation: f32) -> AnalyticVertex {
    return AnalyticVertex(vec4(p, 1.0), vec4(orientation * n, 0.0));
}

fn analytic_circle(u: f32) -> vec2<f32> {
    let turn = fract(u);
    if turn == 0.0 { return vec2(1.0, 0.0); }
    if turn == 0.25 { return vec2(0.0, 1.0); }
    if turn == 0.5 { return vec2(-1.0, 0.0); }
    if turn == 0.75 { return vec2(0.0, -1.0); }
    let angle = 6.283185307179586 * turn;
    return vec2(cos(angle), sin(angle));
}

fn analytic_quarters(u: f32) -> vec2<f32> {
    let quarter = u - 4.0 * floor(u * 0.25);
    let t = fract(quarter);
    let ray = vec2(1.0 - t * t, 2.0 * t) / (1.0 + t * t);
    switch u32(quarter) {
        case 0u: { return ray; }
        case 1u: { return vec2(-ray.y, ray.x); }
        case 2u: { return -ray; }
        default: { return vec2(ray.y, -ray.x); }
    }
}

fn eval_analytic(surface: AnalyticSurface, uv: vec2<f32>) -> AnalyticVertex {
    let d = surface.data;
    let origin = d[0].xyz;
    let axis = d[1].xyz;
    let across = d[2].xyz;
    let side = d[3].xyz;
    let k = d[5];
    var p: vec3<f32>;
    var n: vec3<f32>;
    switch surface.kind {
        case 0u: { // Cone: u is turns, v is axial distance.
            // The apex chooses the seam's limiting normal, independent of u.
            let ray = analytic_circle(select(uv.x, 0.0, uv.y == 0.0));
            let radial = ray.x * across + ray.y * side;
            p = origin + uv.y * (axis + k.x * radial);
            n = normalize(radial - k.x * axis);
        }
        case 1u: { // Sphere: v is latitude in turns.
            let ray = analytic_circle(uv.x);
            let lat = analytic_circle(uv.y);
            n = lat.x * (ray.x * across + ray.y * side) + lat.y * axis;
            p = origin + k.x * n;
        }
        case 2u: { // Torus or its outer sheet on the retained domain.
            let ray = analytic_circle(uv.x);
            let section = analytic_circle(uv.y);
            let radial = ray.x * across + ray.y * side;
            n = section.x * radial + section.y * axis;
            p = origin + k.x * radial + k.y * n;
        }
        case 3u: { // Quadratic or cubic Bezier translation sweep.
            let a = mix(d[0].xyz, d[1].xyz, uv.x);
            let b = mix(d[1].xyz, d[2].xyz, uv.x);
            let ab = mix(a, b, uv.x);
            var profile = ab;
            if k.x == 3.0 {
                let c = mix(d[2].xyz, d[3].xyz, uv.x);
                let bc = mix(b, c, uv.x);
                profile = mix(ab, bc, uv.x);
            }
            p = profile + uv.y * d[4].xyz;
            // Form the tangent from control differences, not translated points.
            let edge0 = d[1].xyz - d[0].xyz;
            let edge1 = d[2].xyz - d[1].xyz;
            var tangent = mix(edge0, edge1, uv.x);
            if k.x == 3.0 {
                let edge2 = d[3].xyz - d[2].xyz;
                tangent = mix(tangent, mix(edge1, edge2, uv.x), uv.x);
            }
            n = normalize(cross(tangent, d[4].xyz));
        }
        case 4u: { // Rolling pipe, on either a cylinder or a cone wall.
            let ray = analytic_quarters(uv.x);
            let radial = ray.x * across + ray.y * side;
            let centre = origin + (k.x * radial + k.y * ray.x * axis) / (1.0 + k.w * ray.x);
            let first = d[4].xyz;
            let last = (radial - k.z * axis) * d[6].x;
            let chord = (2.0 - uv.y) * first + uv.y * last;
            n = (2.0 * dot(first, chord) / dot(chord, chord)) * chord - first;
            p = centre + d[6].y * n;
        }
        default: { // Oblique draft: a ruled polar section, affine in v.
            let ray = analytic_quarters(uv.x);
            let radial = ray.x * across + ray.y * side;
            let azimuth = -ray.y * across + ray.x * side;
            let lean = d[4].xy;
            let z = mix(k.x, k.y, uv.y);
            let radius = k.z - k.w * z + dot(lean, ray);
            let du = radius * azimuth + dot(lean, vec2(-ray.y, ray.x)) * radial;
            let dv = axis - k.w * radial;
            p = origin + z * axis + radius * radial;
            n = normalize(cross(du, dv));
        }
    }
    return analytic_vertex(p, n, surface.orientation);
}
