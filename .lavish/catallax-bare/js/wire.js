/* wire.js: a tiny wireframe engine. Column-major 4x4 matrices, perspective
   projection, Canvas 2D strokes batched by edge tag and depth band.
   No dependencies. Shared by the site and the video render pages. */
(function (root) {
  'use strict';

  const F = Float64Array;
  const mat = {
    identity: () => new F([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]),
    multiply(a, b) {
      const o = new F(16);
      for (let c = 0; c < 4; c++)
        for (let r = 0; r < 4; r++)
          o[c * 4 + r] = a[r] * b[c * 4] + a[4 + r] * b[c * 4 + 1] + a[8 + r] * b[c * 4 + 2] + a[12 + r] * b[c * 4 + 3];
      return o;
    },
    chain: (...ms) => ms.reduce((a, b) => mat.multiply(a, b)),
    rotateX(t) { const c = Math.cos(t), s = Math.sin(t); return new F([1, 0, 0, 0, 0, c, s, 0, 0, -s, c, 0, 0, 0, 0, 1]); },
    rotateY(t) { const c = Math.cos(t), s = Math.sin(t); return new F([c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1]); },
    rotateZ(t) { const c = Math.cos(t), s = Math.sin(t); return new F([c, s, 0, 0, -s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]); },
    translate: (x, y, z) => new F([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, x, y, z, 1]),
    scale: (x, y = x, z = x) => new F([x, 0, 0, 0, 0, y, 0, 0, 0, 0, z, 0, 0, 0, 0, 1]),
    perspective(fovy, aspect, near, far) {
      const f = 1 / Math.tan(fovy / 2), nf = 1 / (near - far);
      return new F([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0]);
    },
  };

  // Mesh: flat vertex array, edge index pairs, one tag per edge.
  // Tags pick a stroke style: 0 primary, 1 secondary, 2 accent, 3 hidden/dashed.
  class Mesh {
    constructor() { this.v = []; this.e = []; this.tag = []; this.dots = []; }
    get vertexCount() { return this.v.length / 3; }
    get edgeCount() { return this.tag.length; }
    vert(x, y, z) { this.v.push(x, y, z); return this.v.length / 3 - 1; }
    edge(a, b, tag = 0) { this.e.push(a, b); this.tag.push(tag); }
    loop(pts, tag = 0, closed = true) {
      const ids = pts.map(p => this.vert(p[0], p[1], p[2]));
      for (let i = 0; i < ids.length - (closed ? 0 : 1); i++) this.edge(ids[i], ids[(i + 1) % ids.length], tag);
      return ids;
    }
  }

  const TAU = Math.PI * 2;
  const circle = (r, n, at) => Array.from({ length: n }, (_, i) => at(Math.cos(i / n * TAU) * r, Math.sin(i / n * TAU) * r));

  const gen = {
    // Spur gear in the XY plane, axis Z. Trapezoid teeth, hub, bore, lightening holes.
    gear({ teeth = 18, outer = 1, root = 0.84, bore = 0.2, hub = 0.4, thick = 0.26, holes = 6 } = {}) {
      const m = new Mesh(), h = thick / 2, a = TAU / teeth, prof = [];
      for (let i = 0; i < teeth; i++) {
        const t = i * a;
        prof.push([root, t], [outer, t + 0.27 * a], [outer, t + 0.48 * a], [root, t + 0.75 * a]);
      }
      const ring = z => prof.map(([r, t]) => [Math.cos(t) * r, Math.sin(t) * r, z]);
      const f = m.loop(ring(h)), b = m.loop(ring(-h));
      f.forEach((id, i) => m.edge(id, b[i], 1));
      const hubF = m.loop(circle(hub, 48, (x, y) => [x, y, h * 1.4]), 1);
      const hubB = m.loop(circle(hub, 48, (x, y) => [x, y, -h * 1.4]), 1);
      for (let i = 0; i < 48; i += 6) m.edge(hubF[i], hubB[i], 1);
      const boreF = m.loop(circle(bore, 40, (x, y) => [x, y, h * 1.4]), 2);
      const boreB = m.loop(circle(bore, 40, (x, y) => [x, y, -h * 1.4]), 2);
      for (let i = 0; i < 40; i += 5) m.edge(boreF[i], boreB[i], 2);
      for (let k = 0; k < holes; k++) {
        const c = k / holes * TAU, cx = Math.cos(c) * 0.62, cy = Math.sin(c) * 0.62;
        m.loop(circle(0.11, 20, (x, y) => [cx + x, cy + y, h]), 0);
        m.loop(circle(0.11, 20, (x, y) => [cx + x, cy + y, -h]), 3);
      }
      return m;
    },

    torus({ R = 0.74, r = 0.3, nu = 56, nv = 16 } = {}) {
      const m = new Mesh();
      const p = (u, v) => [(R + r * Math.cos(v)) * Math.cos(u), r * Math.sin(v), (R + r * Math.cos(v)) * Math.sin(u)];
      for (let j = 0; j < nv; j++) m.loop(Array.from({ length: nu }, (_, i) => p(i / nu * TAU, j / nv * TAU)), j % 4 === 0 ? 0 : 1);
      for (let i = 0; i < nu; i += 2) m.loop(Array.from({ length: nv }, (_, j) => p(i / nu * TAU, j / nv * TAU)), i % 8 === 0 ? 2 : 1);
      return m;
    },

    // Cube with a construction grid on every face.
    lattice({ s = 0.8, n = 6 } = {}) {
      const m = new Mesh();
      for (let axis = 0; axis < 3; axis++)
        for (const side of [-s, s])
          for (let k = 0; k <= n; k++) {
            const t = -s + (2 * s * k) / n, tag = k === 0 || k === n ? 0 : 1;
            const pt = (u, w) => { const q = [0, 0, 0]; q[axis] = side; q[(axis + 1) % 3] = u; q[(axis + 2) % 3] = w; return q; };
            m.loop([pt(t, -s), pt(t, s)], tag, false);
            m.loop([pt(-s, t), pt(s, t)], tag, false);
          }
      const d = s * 0.42;
      const inner = [[-d, -d, -d], [d, -d, -d], [d, d, -d], [-d, d, -d], [-d, -d, d], [d, -d, d], [d, d, d], [-d, d, d]].map(p => m.vert(...p));
      [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]].forEach(([a, b]) => m.edge(inner[a], inner[b], 2));
      return m;
    },

    box(sx, sy, sz, tag = 0, at = [0, 0, 0]) {
      const m = new Mesh();
      const c = [[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]]
        .map(([x, y, z]) => m.vert(at[0] + x * sx, at[1] + y * sy, at[2] + z * sz));
      [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [0, 4], [1, 5], [2, 6], [3, 7]].forEach(([a, b]) => m.edge(c[a], c[b], tag));
      m.dots = c;
      return m;
    },

    // Cylinder along Y.
    cylinder({ r = 0.4, y0 = -0.5, y1 = 0.5, seg = 48, rails = 12, tag = 0, railTag = tag } = {}) {
      const m = new Mesh();
      const a = m.loop(circle(r, seg, (x, z) => [x, y0, z]), tag);
      const b = m.loop(circle(r, seg, (x, z) => [x, y1, z]), tag);
      for (let i = 0; i < seg; i += seg / rails) m.edge(a[i], b[i], railTag);
      return m;
    },

    // Triangulated UV sphere: the "mesh kernel" view.
    sphere({ r = 0.9, lat = 10, lon = 22 } = {}) {
      const m = new Mesh(), id = [];
      const top = m.vert(0, r, 0), bot = m.vert(0, -r, 0);
      for (let i = 1; i < lat; i++) {
        const phi = i / lat * Math.PI, row = [];
        for (let j = 0; j < lon; j++) {
          const th = j / lon * TAU + (i % 2) * (Math.PI / lon);
          row.push(m.vert(Math.sin(phi) * Math.cos(th) * r, Math.cos(phi) * r, Math.sin(phi) * Math.sin(th) * r));
        }
        id.push(row);
      }
      for (let i = 0; i < id.length; i++)
        for (let j = 0; j < lon; j++) {
          const a = id[i][j], b = id[i][(j + 1) % lon];
          m.edge(a, b, 0);
          if (i === 0) m.edge(top, a, 1);
          if (i === id.length - 1) m.edge(bot, a, 1);
          if (i < id.length - 1) {
            const sh = i % 2 ? 1 : 0;
            m.edge(a, id[i + 1][j], 1);
            m.edge(a, id[i + 1][(j + lon - 1 + sh) % lon], 1);
          }
        }
      return m;
    },

    // Bicubic Bezier patch (a NURBS with uniform weights). Pass `into` to update vertices in place.
    patch(ctrl, n = 14, into) {
      const B = t => { const s = 1 - t; return [s * s * s, 3 * s * s * t, 3 * s * t * t, t * t * t]; };
      const P = (u, v) => {
        const bu = B(u), bv = B(v), o = [0, 0, 0];
        for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) for (let k = 0; k < 3; k++) o[k] += bu[i] * bv[j] * ctrl[i][j][k];
        return o;
      };
      const pts = [];
      for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) pts.push(ctrl[i][j]);
      for (let i = 0; i <= n; i++) for (let j = 0; j <= n; j++) pts.push(P(i / n, j / n));
      if (into) { pts.forEach((p, k) => { into.v[k * 3] = p[0]; into.v[k * 3 + 1] = p[1]; into.v[k * 3 + 2] = p[2]; }); return into; }
      const m = new Mesh();
      pts.forEach(p => m.vert(...p));
      for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
        if (j < 3) m.edge(i * 4 + j, i * 4 + j + 1, 3);
        if (i < 3) m.edge(i * 4 + j, (i + 1) * 4 + j, 3);
      }
      const s = k => 16 + k;
      for (let i = 0; i <= n; i++) for (let j = 0; j < n; j++) {
        const edgeRow = i === 0 || i === n;
        m.edge(s(i * (n + 1) + j), s(i * (n + 1) + j + 1), edgeRow ? 2 : 0);
        m.edge(s(j * (n + 1) + i), s((j + 1) * (n + 1) + i), edgeRow ? 2 : 1);
      }
      m.dots = Array.from({ length: 16 }, (_, k) => k);
      return m;
    },

    merge(...meshes) {
      const m = new Mesh();
      for (const src of meshes) {
        const off = m.vertexCount;
        m.v.push(...src.v);
        src.e.forEach(i => m.e.push(i + off));
        m.tag.push(...src.tag);
        m.dots.push(...src.dots.map(i => i + off));
      }
      return m;
    },
  };

  // Reads brand colors from CSS custom properties, so nothing is hardcoded here.
  function palette(el = document.documentElement) {
    const cs = getComputedStyle(el), v = n => cs.getPropertyValue(n).trim();
    return { ink: v('--ink'), steel: v('--steel'), slate: v('--slate'), bone: v('--bone'), oxide: v('--oxide'), signal: v('--signal'), datum: v('--datum') };
  }

  function defaultStyles(p) {
    return [
      { color: p.bone, width: 1, alpha: 0.95 },
      { color: p.slate, width: 1, alpha: 0.7 },
      { color: p.oxide, width: 1.6, alpha: 1 },
      { color: p.slate, width: 1, alpha: 0.55, dash: [3, 4] },
    ];
  }

  const BANDS = 4;

  class Renderer {
    constructor(canvas, { fov = 0.72, maxDpr = 2 } = {}) {
      this.canvas = canvas;
      this.ctx = canvas.getContext('2d');
      this.fov = fov;
      this.maxDpr = maxDpr;
      this.view = mat.translate(0, 0, -4);
      this.resize();
    }
    resize(w = this.canvas.clientWidth, h = this.canvas.clientHeight) {
      this.dpr = Math.min(window.devicePixelRatio || 1, this.maxDpr);
      this.w = w; this.h = h;
      this.canvas.width = Math.round(w * this.dpr);
      this.canvas.height = Math.round(h * this.dpr);
      this.proj = mat.perspective(this.fov, w / Math.max(h, 1), 0.1, 100);
    }
    clear() {
      const c = this.ctx;
      c.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
      c.clearRect(0, 0, this.w, this.h);
    }
    // Projects every vertex: clip = P * V * M * v, then perspective divide to pixels.
    project(mesh, model) {
      const m = mat.chain(this.proj, this.view, model), n = mesh.vertexCount;
      if (!mesh._s || mesh._s.length !== n * 3) mesh._s = new Float32Array(n * 3);
      const s = mesh._s, v = mesh.v, W = this.w, H = this.h;
      let zMin = Infinity, zMax = -Infinity;
      for (let i = 0; i < n; i++) {
        const x = v[i * 3], y = v[i * 3 + 1], z = v[i * 3 + 2];
        const cx = m[0] * x + m[4] * y + m[8] * z + m[12];
        const cy = m[1] * x + m[5] * y + m[9] * z + m[13];
        const cw = m[3] * x + m[7] * y + m[11] * z + m[15];
        s[i * 3] = (cx / cw * 0.5 + 0.5) * W;
        s[i * 3 + 1] = (0.5 - cy / cw * 0.5) * H;
        s[i * 3 + 2] = cw;
        if (cw < zMin) zMin = cw;
        if (cw > zMax) zMax = cw;
      }
      return { s, zMin, span: Math.max(zMax - zMin, 1e-6), mvp: m };
    }
    draw(mesh, model, { progress = 1, alpha = 1, styles, head = null, dots = null, fog = 0.78 } = {}) {
      if (alpha <= 0.002 || progress <= 0) return null;
      const ctx = this.ctx, { s, zMin, span, mvp } = this.project(mesh, model);
      const E = mesh.edgeCount, shown = Math.min(E, progress * E), full = Math.floor(shown), frac = shown - full;
      const e = mesh.e, tags = mesh.tag;
      const band = i => {
        const z = (s[e[i * 2] * 3 + 2] + s[e[i * 2 + 1] * 3 + 2]) * 0.5;
        return Math.min(BANDS - 1, ((z - zMin) / span * BANDS) | 0);
      };
      ctx.lineCap = 'round';
      for (let t = 0; t < styles.length; t++) {
        const st = styles[t];
        ctx.strokeStyle = st.color;
        ctx.lineWidth = st.width;
        ctx.setLineDash(st.dash || []);
        for (let b = 0; b < BANDS; b++) {
          ctx.globalAlpha = alpha * st.alpha * (1 - fog * b / (BANDS - 1));
          ctx.beginPath();
          let any = false;
          for (let i = 0; i < full; i++) {
            if (tags[i] !== t || band(i) !== b) continue;
            const a = e[i * 2] * 3, c = e[i * 2 + 1] * 3;
            ctx.moveTo(s[a], s[a + 1]);
            ctx.lineTo(s[c], s[c + 1]);
            any = true;
          }
          if (any) ctx.stroke();
        }
      }
      ctx.setLineDash([]);
      let tip = null;
      if (full < E && frac > 0) {
        const a = e[full * 2] * 3, c = e[full * 2 + 1] * 3, st = styles[tags[full]];
        tip = [s[a] + (s[c] - s[a]) * frac, s[a + 1] + (s[c + 1] - s[a + 1]) * frac];
        ctx.globalAlpha = alpha * st.alpha;
        ctx.strokeStyle = st.color;
        ctx.lineWidth = st.width;
        ctx.beginPath(); ctx.moveTo(s[a], s[a + 1]); ctx.lineTo(tip[0], tip[1]); ctx.stroke();
        if (head) {
          ctx.globalAlpha = alpha;
          ctx.fillStyle = head;
          ctx.beginPath(); ctx.arc(tip[0], tip[1], 2.6, 0, TAU); ctx.fill();
        }
      }
      if (dots && mesh.dots.length) {
        ctx.fillStyle = dots;
        for (const i of mesh.dots) {
          ctx.globalAlpha = alpha * Math.min(1, progress * 1.5) * (1 - 0.6 * (s[i * 3 + 2] - zMin) / span);
          ctx.fillRect(s[i * 3] - 2.5, s[i * 3 + 1] - 2.5, 5, 5);
        }
      }
      ctx.globalAlpha = 1;
      return { tip, mvp };
    }
  }

  root.Wire = { mat, Mesh, gen, palette, defaultStyles, Renderer };
})(window);
