// Tensor-product carrier evaluation. The enclosing shader binds
// spline_controls: array<vec4f> and spline_knots: array<f32> as read-only storage.
// Controls are (w * local_xyz, w), u-major with v contiguous. Unit weights give
// the polynomial case. Offsets index elements, not bytes.
struct SplineNet {
    control_offset: u32,
    knot_u_offset: u32,
    knot_v_offset: u32,
    degree_u: u32,
    degree_v: u32,
    count_u: u32,
    count_v: u32,
}

struct SplineSample {
    position: vec3f,
    // Parametric orientation (du × dv); zero when the tangents do not span a plane.
    normal: vec3f,
}

struct SplineBasis {
    first: u32,
    value: array<f32, 9>,
    // Scaled degree-(p-1) weights for adjacent control differences.
    difference: array<f32, 8>,
}

// Find the right-hand span at interior knots, and the left-hand span at the
// upper endpoint. The active interval is closed; repeated knots are allowed.
fn spline_basis(offset: u32, degree: u32, count: u32, parameter: f32) -> SplineBasis {
    let lo = spline_knots[offset + degree];
    let hi = spline_knots[offset + count];
    let t = clamp(parameter, lo, hi);
    var left = degree;
    var right = count;
    loop {
        if right - left <= 1u { break; }
        let mid = (left + right) / 2u;
        if t < spline_knots[offset + mid] { right = mid; }
        else { left = mid; }
    }
    // A non-clamped net can have repeated knots at its active upper endpoint.
    var span = left;
    while spline_knots[offset + span] == hi { span -= 1u; }

    var basis: SplineBasis;
    basis.value[0] = 1.0;
    // Cox–de Boor on the nonzero span. At the last level, use the degree-(p-1)
    // basis to get the derivative without differentiating a sampled position.
    for (var level = 1u; level <= degree; level += 1u) {
        var carry = 0.0;
        for (var i = 0u; i < level; i += 1u) {
            let a = spline_knots[offset + span + i + 1u - level];
            let b = spline_knots[offset + span + i + 1u];
            let term = basis.value[i] / (b - a);
            if level == degree {
                basis.difference[i] = f32(degree) * term;
            }
            basis.value[i] = carry + (b - t) * term;
            carry = (t - a) * term;
        }
        basis.value[level] = carry;
    }
    return basis;
}

fn spline_scaled(v: vec3f) -> vec3f {
    let m = max(max(abs(v.x), abs(v.y)), abs(v.z));
    if !(m > 0.0) { return vec3f(0.0); }
    return v / m;
}

// Direct evaluation, with at most 9×3 control reads for the 8×2 case.
// Both axes can use the same basis routine (including when a net is transposed).
fn spline_sample(net: SplineNet, uv: vec2f) -> SplineSample {
    let u = spline_basis(net.knot_u_offset, net.degree_u, net.count_u, uv.x);
    let v = spline_basis(net.knot_v_offset, net.degree_v, net.count_v, uv.y);
    var h = vec4f(0.0);
    var hu = vec4f(0.0);
    var hv = vec4f(0.0);
    var previous_row = vec4f(0.0);
    for (var i = 0u; i <= net.degree_u; i += 1u) {
        let row = net.control_offset + (u.first + i) * net.count_v + v.first;
        var value = vec4f(0.0);
        var dv = vec4f(0.0);
        var previous_control = vec4f(0.0);
        for (var j = 0u; j <= net.degree_v; j += 1u) {
            let control = spline_controls[row + j];
            value += v.value[j] * control;
            if j > 0u {
                dv += v.difference[j - 1u] * (control - previous_control);
            }
            previous_control = control;
        }
        h += u.value[i] * value;
        if i > 0u {
            hu += u.difference[i - 1u] * (value - previous_row);
        }
        previous_row = value;
        hv += u.value[i] * dv;
    }
    let p = h.xyz / h.w;
    // The common positive 1/h.w factor has no effect on the normal.
    let du = spline_scaled(hu.xyz - p * hu.w);
    let dv = spline_scaled(hv.xyz - p * hv.w);
    let n = spline_scaled(cross(du, dv));
    var normal = vec3f(0.0);
    if any(n != vec3f(0.0)) { normal = normalize(n); }
    return SplineSample(p, normal);
}
