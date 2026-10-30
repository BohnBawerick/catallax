/* field.js: the hero gear as a signed distance field, ray marched in WebGL2.
   The same parameters as Wire.gen.gear, so the field and the edges line up.
   Illustration only. This is not kernel output. */
(function (root) {
  'use strict';

  const TAU = Math.PI * 2, TEETH = 18, A = TAU / TEETH;
  const OUTER = 1, RT = 0.84, BORE = 0.2, HUB = 0.4, H = 0.13, HUB_H = 0.182, HOLE_R = 0.11, HOLE_C = 0.62;
  const v2 = (r, t) => `vec2(${(r * Math.cos(t)).toFixed(6)},${(r * Math.sin(t)).toFixed(6)})`;

  const VERT = `#version 300 es
  in vec2 p; void main() { gl_Position = vec4(p, 0., 1.); }`;

  const FRAG = `#version 300 es
  precision highp float;
  uniform vec2 uRes; uniform mat4 uInv; uniform float uFov, uT;
  uniform vec3 uInk, uBone, uHot, uSlate, uGraphite;
  out vec4 outColor;
  const float PI2 = 6.28318530718;
  float seg(vec2 p, vec2 a, vec2 b) { vec2 pa = p - a, ba = b - a; return length(pa - ba * clamp(dot(pa, ba) / dot(ba, ba), 0., 1.)); }

  // One tooth sector of the profile. The chain runs counterclockwise, so the left of a segment is inside.
  float profile(vec2 p) {
    float a = PI2 / ${TEETH}.;
    float th = atan(p.y, p.x);
    float t = th - floor(th / a) * a;
    vec2 q = length(p) * vec2(cos(t), sin(t));
    float u = t / a;
    const vec2 P0 = ${v2(RT, 0)}, P1 = ${v2(OUTER, 0.27 * A)}, P2 = ${v2(OUTER, 0.48 * A)}, P3 = ${v2(RT, 0.75 * A)}, P4 = ${v2(RT, A)};
    float d = min(seg(q, ${v2(OUTER, -0.52 * A)}, ${v2(RT, -0.25 * A)}), seg(q, ${v2(RT, -0.25 * A)}, P0));
    d = min(d, min(seg(q, P0, P1), seg(q, P1, P2)));
    d = min(d, min(seg(q, P2, P3), seg(q, P3, P4)));
    d = min(d, seg(q, P4, ${v2(OUTER, 1.27 * A)}));
    vec2 s0 = u < .27 ? P0 : u < .48 ? P1 : u < .75 ? P2 : P3;
    vec2 s1 = u < .27 ? P1 : u < .48 ? P2 : u < .75 ? P3 : P4;
    float side = (s1.x - s0.x) * (q.y - s0.y) - (s1.y - s0.y) * (q.x - s0.x);
    return side > 0. ? -d : d;
  }
  float extrude(float d, float y, float h) { vec2 w = vec2(d, abs(y) - h); return min(max(w.x, w.y), 0.) + length(max(w, 0.)); }
  float map(vec3 p) {
    vec2 xz = vec2(p.x, -p.z);
    float r = length(xz);
    float d = min(extrude(profile(xz), p.y, ${H}), extrude(r - ${HUB}, p.y, ${HUB_H}));
    d = max(d, -(r - ${BORE}));
    float a = PI2 / 6., th = atan(xz.y, xz.x), k = th - a * floor(th / a + .5);
    return max(d, -(length(r * vec2(cos(k), sin(k)) - vec2(${HOLE_C}, 0.)) - ${HOLE_R}));
  }
  vec3 normal(vec3 p) {
    vec2 e = vec2(.0015, 0.);
    return normalize(vec3(map(p + e.xyy) - map(p - e.xyy), map(p + e.yxy) - map(p - e.yxy), map(p + e.yyx) - map(p - e.yyx)));
  }
  vec3 contours(vec3 q) {
    float d = map(q), s = 9.;
    float flow = d * s - uT * .35;
    float w = fwidth(flow);
    float fine = 1. - smoothstep(0., 1.3, abs(fract(flow - .5) - .5) / w);
    float major = 1. - smoothstep(0., 1.6, abs(fract(flow / 5. - .5) - .5) * 5. / w);
    float fade = exp(-max(d, 0.) * .72);
    vec3 line = mix(uSlate, uHot, exp(-max(d, 0.) * 3.2));
    vec3 c = uInk + line * fine * .55 * fade + uBone * major * .4 * fade;
    return c + uHot * .10 * exp(-max(d, 0.) * 14.);
  }
  vec3 shade(vec3 p, vec3 rd) {
    vec3 n = normal(p), l = normalize(vec3(.5, .9, .45));
    float dif = max(dot(n, l), 0.), fres = pow(1. - max(dot(n, -rd), 0.), 3.);
    float inBore = 1. - smoothstep(${BORE + 0.01}, ${BORE + 0.05}, length(p.xz));
    vec3 base = mix(uGraphite * 1.5, uHot * .75, inBore);
    return base * (.35 + .9 * dif) + uBone * fres * .5 + uBone * pow(max(dot(reflect(rd, n), l), 0.), 40.) * .35;
  }
  void main() {
    vec2 ndc = gl_FragCoord.xy / uRes * 2. - 1.;
    float tf = tan(uFov * .5);
    vec3 rv = normalize(vec3(ndc.x * uRes.x / uRes.y * tf, ndc.y * tf, -1.));
    vec3 ro = (uInv * vec4(0., 0., 0., 1.)).xyz;
    vec3 rd = normalize((uInv * vec4(rv, 0.)).xyz);
    float tp = abs(rd.y) > 1e-4 ? -ro.y / rd.y : -1.;
    vec3 col = uInk;
    if (tp > 0.) col = contours(ro + rd * tp);

    // Bounding sphere keeps the march short. The closest miss antialiases the silhouette.
    float b = dot(ro, rd), c = dot(ro, ro) - 1.25 * 1.25, h = b * b - c;
    if (h > 0.) {
      float px = 2. * tf / uRes.y, t = max(-b - sqrt(h), .001), tEnd = -b + sqrt(h), hit = -1., cover = 1e3, near = t;
      for (int i = 0; i < 72 && t < tEnd; i++) {
        float d = map(ro + rd * t), k = d / (t * px);
        if (k < cover) { cover = k; near = t; }
        if (d < .0012) { hit = t; break; }
        t += d * .92;
      }
      float amt = hit > 0. ? 1. : 1. - clamp(cover, 0., 1.);
      float at = hit > 0. ? hit : near;
      if (amt > 0. && (tp < 0. || at < tp)) col = mix(col, shade(ro + rd * at, rd), amt);
    }
    outColor = vec4(pow(col, vec3(.95)), 1.);
  }`;

  const rgb = hex => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16) / 255);

  function FieldGear(canvas, palette, { fov = 0.72 } = {}) {
    const gl = canvas.getContext('webgl2', { antialias: false, alpha: false, powerPreference: 'high-performance' });
    if (!gl) return null;
    const compile = (type, src) => {
      const s = gl.createShader(type);
      gl.shaderSource(s, src);
      gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      return s;
    };
    const prog = gl.createProgram();
    gl.attachShader(prog, compile(gl.VERTEX_SHADER, VERT));
    gl.attachShader(prog, compile(gl.FRAGMENT_SHADER, FRAG));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
    gl.useProgram(prog);
    gl.bindBuffer(gl.ARRAY_BUFFER, gl.createBuffer());
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(prog, 'p');
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
    const U = n => gl.getUniformLocation(prog, n);
    gl.uniform1f(U('uFov'), fov);
    gl.uniform3fv(U('uInk'), rgb(palette.ink));
    gl.uniform3fv(U('uBone'), rgb(palette.bone));
    gl.uniform3fv(U('uHot'), rgb(palette.oxide));
    gl.uniform3fv(U('uSlate'), rgb(palette.slate));
    gl.uniform3fv(U('uGraphite'), rgb('#131619'));
    const uRes = U('uRes'), uInv = U('uInv'), uT = U('uT');
    return {
      resize(w, h, scale) {
        canvas.width = Math.max(1, Math.round(w * scale));
        canvas.height = Math.max(1, Math.round(h * scale));
        gl.viewport(0, 0, canvas.width, canvas.height);
        gl.uniform2f(uRes, canvas.width, canvas.height);
      },
      render(inverse, time) {
        gl.uniformMatrix4fv(uInv, false, Float32Array.from(inverse));
        gl.uniform1f(uT, time);
        gl.drawArrays(gl.TRIANGLES, 0, 3);
      },
    };
  }

  root.FieldGear = FieldGear;
})(window);
