/* Site motion: rigid metal, expo.out, no overshoot.
   The page reads complete without this file. */
(() => {
  'use strict';
  const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
  const { mat, gen, Renderer, palette, defaultStyles } = Wire;
  const P = palette();
  const styles = defaultStyles(P);
  const $ = (s, el = document) => el.querySelector(s);
  const $$ = (s, el = document) => [...el.querySelectorAll(s)];

  gsap.registerPlugin(ScrollTrigger);
  gsap.defaults({ ease: 'expo.out', duration: 0.9 });

  let paused = reduce;

  function loopWhileVisible(el, frame) {
    let raf = 0, last = 0, on = false;
    const tick = now => {
      frame(Math.min(0.05, (now - last) / 1000 || 0), now);
      last = now;
      raf = requestAnimationFrame(tick);
    };
    new IntersectionObserver(([e]) => {
      if (e.isIntersecting && !on) { on = true; last = performance.now(); raf = requestAnimationFrame(tick); }
      else if (!e.isIntersecting && on) { on = false; cancelAnimationFrame(raf); }
    }).observe(el);
  }

  function changes() {
    let last = '';
    return (...values) => {
      const key = values.join();
      if (key === last) return false;
      last = key;
      return true;
    };
  }

  const flatGear = () => {
    const g = gen.gear();
    for (let i = 0; i < g.v.length; i += 3) { const y = g.v[i + 1]; g.v[i + 1] = g.v[i + 2]; g.v[i + 2] = -y; }
    return g;
  };

  /* Pause: one icon, everything that moves on its own ---------------------- */
  const pauseBtn = $('.motion'), videos = $$('video'), near = new Set();
  pauseBtn.hidden = reduce;
  pauseBtn.addEventListener('click', () => {
    paused = !paused;
    pauseBtn.setAttribute('aria-label', paused ? 'Resume motion' : 'Pause motion');
    $('use', pauseBtn).setAttribute('href', paused ? '#s-play' : '#s-pause');
    videos.forEach(v => { if (paused) v.pause(); else if (near.has(v)) v.play().catch(() => {}); });
  });
  videos.forEach(v => {
    if (reduce) { v.controls = true; return; }
    new IntersectionObserver(([e]) => {
      if (e.isIntersecting) { near.add(v); if (!paused) v.play().catch(() => {}); }
      else { near.delete(v); v.pause(); }
    }, { rootMargin: '25% 0px' }).observe(v);
  });

  /* Hero: wireframe with live leader lines --------------------------------- */
  const heroCanvas = $('.hero__canvas');
  const hero = new Renderer(heroCanvas, { fov: 0.72 });
  const shapes = { gear: flatGear(), torus: gen.torus(), lattice: gen.lattice() };
  const S = { mesh: shapes.gear, yaw: 0.6, pitch: -0.5, vyaw: 0, vpitch: 0, progress: reduce ? 1 : 0, x: 1.3, y: 0.05, k: 1.25 };
  const AUTO = 0.14;
  let drag = null, heroStale = true, heroDrawn = 0;
  const leaders = $('.leaders'), lines = $$('line', leaders), dots = $$('circle', leaders), balloons = $$('.balloon');
  let balloonAt = [];

  function layoutHero() {
    hero.resize();
    heroStale = true;
    const wide = innerWidth > 900;
    Object.assign(S, wide ? { x: hero.w / hero.h * 0.68, y: 0.38, k: 0.9 } : { x: 0, y: 0.92, k: 0.8 });
    const box = heroCanvas.getBoundingClientRect();
    balloonAt = balloons.map(b => { const r = b.getBoundingClientRect(); return [r.left - box.left, r.top - box.top + r.height / 2]; });
  }
  layoutHero();
  new ResizeObserver(layoutHero).observe(heroCanvas);

  heroCanvas.addEventListener('pointerdown', e => {
    if (drag) return;
    drag = { id: e.pointerId, x: e.clientX, y: e.clientY };
    heroCanvas.setPointerCapture(e.pointerId);
  });
  heroCanvas.addEventListener('pointermove', e => {
    if (!drag || e.pointerId !== drag.id) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    drag.x = e.clientX; drag.y = e.clientY;
    S.vyaw = dx * 0.008; S.vpitch = dy * 0.006;
    S.yaw += S.vyaw;
    S.pitch = Math.max(-1.4, Math.min(1.4, S.pitch + S.vpitch));
  });
  const endDrag = e => {
    if (!drag || (e.pointerId !== undefined && e.pointerId !== drag.id)) return;
    drag = null;
    if (reduce) S.vyaw = S.vpitch = 0;
  };
  ['pointerup', 'pointercancel', 'lostpointercapture'].forEach(t => heroCanvas.addEventListener(t, endDrag));
  addEventListener('blur', endDrag);

  // Callouts land on real features of the gear: the tooth ring, the bore, the nearest hole.
  const ringAt = (r, y, n) => Array.from({ length: n }, (_, i) => {
    const t = i / n * Math.PI * 2;
    return [r * Math.cos(t), y, -r * Math.sin(t)];
  });
  const FEATURES = [ringAt(1, 0.13, 90), ringAt(0.2, 0.182, 48), ringAt(0.62, 0.13, 6)];
  const leaderAt = [];
  function drawLeaders(m, dt) {
    FEATURES.forEach((points, i) => {
      const [bx, by] = balloonAt[i] || [0, 0];
      let best = null, bestD = Infinity;
      for (const [x, y, z] of points) {
        const w = m[3] * x + m[7] * y + m[11] * z + m[15];
        const sx = ((m[0] * x + m[4] * y + m[8] * z + m[12]) / w * 0.5 + 0.5) * hero.w;
        const sy = (0.5 - (m[1] * x + m[5] * y + m[9] * z + m[13]) / w * 0.5) * hero.h;
        const d = (sx - bx) ** 2 + (sy - by) ** 2;
        if (d < bestD) { bestD = d; best = [sx, sy]; }
      }
      const glide = 1 - Math.exp(-dt * 12);
      const at = leaderAt[i] || (leaderAt[i] = best);
      at[0] += (best[0] - at[0]) * glide; at[1] += (best[1] - at[1]) * glide;
      lines[i].setAttribute('x1', at[0]); lines[i].setAttribute('y1', at[1]);
      lines[i].setAttribute('x2', bx); lines[i].setAttribute('y2', by);
      dots[i].setAttribute('cx', at[0]); dots[i].setAttribute('cy', at[1]);
    });
  }

  const heroChanged = changes();
  loopWhileVisible(heroCanvas, (dt, now) => {
    if (!drag) {
      S.yaw += S.vyaw + (paused ? 0 : AUTO * dt);
      S.pitch = Math.max(-1.4, Math.min(1.4, S.pitch + S.vpitch));
      S.vyaw *= 0.94; S.vpitch *= 0.9;
    }
    const busy = drag || S.progress < 1 || heroStale;
    if (!busy && now - heroDrawn < 32) return;
    if (!heroChanged(S.yaw.toFixed(4), S.pitch.toFixed(4), S.progress, S.x, S.y, S.k, S.mesh.edgeCount) && !heroStale) return;
    heroStale = false; heroDrawn = now;
    const model = mat.chain(mat.translate(S.x, S.y, 0), mat.rotateX(S.pitch), mat.rotateY(S.yaw), mat.scale(S.k));
    hero.clear();
    hero.draw(S.mesh, model, { progress: S.progress, styles, head: P.oxide });
    if (S.progress > 0.8) document.body.classList.add('is-drawn');
    if (S.mesh === shapes.gear && innerWidth > 900) drawLeaders(hero.project(S.mesh, model).mvp, dt);
  });
  if (!reduce) gsap.to(S, { progress: 1, duration: 2.4, ease: 'power2.inOut', delay: 0.2 });

  const shapeBtns = $$('.shape');
  shapeBtns.forEach(btn => btn.addEventListener('click', () => {
    const next = shapes[btn.dataset.shape];
    if (next === S.mesh) return;
    shapeBtns.forEach(b => { b.classList.toggle('is-on', b === btn); b.setAttribute('aria-pressed', b === btn); });
    document.body.classList.toggle('is-gear', btn.dataset.shape === 'gear');
    gsap.killTweensOf(S, 'progress');
    if (reduce) { S.mesh = next; S.progress = 1; heroStale = true; return; }
    gsap.timeline()
      .to(S, { progress: 0, duration: 0.45, ease: 'power2.in' })
      .add(() => { S.mesh = next; })
      .to(S, { progress: 1, duration: 1.6, ease: 'power2.inOut' });
  }));
  document.body.classList.add('is-gear');

  /* Twin: the same gear as edges and as a field, wiped by a handle ---------- */
  const twin = $('.twin'), twinWire = $('.twin__wire'), twinField = $('.twin__field'), handle = $('.twin__handle');
  const gearEdges = flatGear();
  const wire = new Renderer(twinWire, { fov: 0.72 });
  let field = null;
  try { field = FieldGear(twinField, P, { fov: 0.72 }); } catch (e) { console.warn(e); }
  $('.recipe').classList.toggle('has-field', !!field);
  const T = { yaw: 0.5, pitch: 0.62, split: 0.5, t: 0 };
  const setSplit = v => {
    T.split = Math.max(0.04, Math.min(0.96, v));
    twin.style.setProperty('--split', `${T.split * 100}%`);
    handle.setAttribute('aria-valuenow', Math.round(T.split * 100));
  };
  setSplit(0.5);
  let twinStale = true, fieldScale = 1, slowFrames = 0;
  const FIELD_PIXELS = 1.6e6;
  function layoutTwin() {
    wire.resize();
    const w = twin.clientWidth, h = twin.clientHeight;
    if (field) field.resize(w, h, Math.min(devicePixelRatio || 1, Math.sqrt(FIELD_PIXELS / (w * h))) * fieldScale);
    twinStale = true;
  }
  layoutTwin();
  new ResizeObserver(layoutTwin).observe(twin);

  let dragging = null;
  const fromPointer = e => setSplit((e.clientX - twin.getBoundingClientRect().left) / twin.clientWidth);
  handle.addEventListener('pointerdown', e => { dragging = e.pointerId; handle.setPointerCapture(e.pointerId); });
  handle.addEventListener('pointermove', e => { if (e.pointerId === dragging) fromPointer(e); });
  ['pointerup', 'pointercancel', 'lostpointercapture'].forEach(t => handle.addEventListener(t, () => { dragging = null; }));
  handle.addEventListener('keydown', e => {
    const step = { ArrowLeft: -0.05, ArrowRight: 0.05, ArrowDown: -0.05, ArrowUp: 0.05 }[e.key];
    if (step) { e.preventDefault(); setSplit(T.split + step); }
    if (e.key === 'Home') { e.preventDefault(); setSplit(0.04); }
    if (e.key === 'End') { e.preventDefault(); setSplit(0.96); }
  });

  const twinChanged = changes();
  let twinDrawn = 0;
  loopWhileVisible(twin, (dt, now) => {
    if (!paused) { T.yaw += 0.1 * dt; T.t += dt; }
    // A GPU that cannot hold the rate gets fewer pixels, not a stutter.
    slowFrames = dt > 0.045 ? slowFrames + 1 : 0;
    if (field && slowFrames > 24 && fieldScale > 0.4) { fieldScale *= 0.7; slowFrames = 0; layoutTwin(); }
    if (!twinStale && now - twinDrawn < 32) return;
    if (!twinChanged(T.yaw.toFixed(4), T.t.toFixed(3)) && !twinStale) return;
    twinStale = false; twinDrawn = now;
    const wide = innerWidth > 900, k = wide ? 1.3 : 1.02, y = wide ? 0.34 : 0.22;
    const model = mat.chain(mat.translate(0, y, 0), mat.rotateX(T.pitch), mat.rotateY(T.yaw), mat.scale(k));
    wire.clear();
    wire.draw(gearEdges, model, { styles });
    if (field) field.render(mat.chain(mat.scale(1 / k), mat.rotateY(-T.yaw), mat.rotateX(-T.pitch), mat.translate(0, -y, 4)), T.t);
  });

  /* Kernel anatomy --------------------------------------------------------- */
  const anatomy = $('.anatomy'), aCanvas = $('.anatomy__canvas');
  const ana = new Renderer(aCanvas, { fov: 0.7 });
  ana.view = mat.translate(0, 0, -4.1);
  let anaStale = true;
  new ResizeObserver(() => { ana.resize(); anaStale = true; }).observe(aCanvas);
  const PIN_ROOM = '(min-height: 640px)';
  anatomy.classList.toggle('is-pinned', matchMedia(PIN_ROOM).matches);

  const brep = gen.box(0.6, 0.6, 0.6);
  [0, 1, 2, 3].forEach(i => { brep.tag[i] = 2; });
  const ctrl = Array.from({ length: 4 }, (_, i) => Array.from({ length: 4 }, (_, j) =>
    [(-1 + 2 * i / 3) * 1.05, 0.32 * Math.sin(i * 1.7) * Math.cos(j * 1.3) - 0.1, (-1 + 2 * j / 3) * 1.05]));
  const nurbs = gen.patch(ctrl);
  const bool = gen.merge(
    gen.box(0.9, 0.42, 0.9),
    gen.cylinder({ r: 0.44, y0: -1, y1: 1, seg: 48, rails: 8, tag: 3 }),
    gen.cylinder({ r: 0.44, y0: -0.42, y1: 0.42, seg: 48, rails: 16, tag: 2, railTag: 1 }));
  const layers = [brep, nurbs, bool, gen.sphere({ r: 1 })].map((mesh, i) => ({ mesh, alpha: i ? 0 : 1, progress: reduce ? 1 : 0 }));
  let cur = 0, rot = 0.6, aT = 0, anaStarted = reduce, pinned = false;
  const anaChanged = changes();
  const steps = $$('.step'), stepNum = $('[data-step]');

  loopWhileVisible(aCanvas, dt => {
    if (!anaStarted) showStep(0);
    if (!paused) { rot += dt * 0.3; aT += dt; }
    if (!paused && layers[1].alpha > 0) {
      ctrl[1][2][1] = 0.55 * Math.sin(aT * 1.1) + 0.1;
      ctrl[2][1][1] = -0.45 * Math.sin(aT * 0.8 + 1);
      gen.patch(ctrl, 14, nurbs);
    }
    if (!anaChanged(rot, aT, ...layers.flatMap(l => [l.alpha, l.progress])) && !anaStale) return;
    anaStale = false;
    const wide = innerWidth > 900;
    const model = mat.chain(mat.translate(wide ? ana.w / ana.h * 0.5 : 0, wide ? 0.06 : 0, 0), mat.rotateX(-0.5), mat.rotateY(rot), mat.scale(wide ? 0.92 : Math.min(1, ana.w / ana.h * 0.98)));
    ana.clear();
    for (const l of layers) ana.draw(l.mesh, model, { progress: l.progress, alpha: l.alpha, styles, head: P.oxide, dots: P.oxide });
  });

  function showStep(i) {
    if (i === cur && anaStarted) return;
    const from = cur;
    cur = i;
    anaStarted = true;
    stepNum.textContent = i + 1;
    const out = steps[from], inn = steps[i];
    if (reduce) {
      layers.forEach((l, k) => { l.alpha = k === i ? 1 : 0; l.progress = 1; });
      if (pinned) { gsap.set(out, { opacity: 0 }); gsap.set(inn, { opacity: 1, y: 0 }); }
      return;
    }
    if (from !== i) gsap.to(layers[from], { alpha: 0, duration: 0.45, ease: 'power2.in', overwrite: true });
    gsap.fromTo(layers[i], { alpha: 1, progress: 0 }, { progress: 1, duration: 1.4, ease: 'power2.inOut', overwrite: true });
    if (!pinned) return;
    if (out !== inn) gsap.to(out, { opacity: 0, y: -24, duration: 0.3, ease: 'power2.in', overwrite: true });
    gsap.fromTo(inn, { opacity: 0, y: 30 }, { opacity: 1, y: 0, duration: 0.7, delay: from !== i ? 0.15 : 0, overwrite: true });
  }

  /* Support: the cube in ink on oxide ------------------------------------- */
  const sCanvas = $('.support__canvas');
  const sup = new Renderer(sCanvas, { fov: 0.72 });
  const ink = [
    { color: P.ink, width: 1.3, alpha: 1 },
    { color: P.ink, width: 1, alpha: 0.5 },
    { color: P.ink, width: 2, alpha: 1 },
    { color: P.ink, width: 1, alpha: 0.5, dash: [3, 4] },
  ];
  const cube = gen.lattice();
  let supStale = true, supYaw = 0.7;
  new ResizeObserver(() => { sup.resize(); supStale = true; }).observe(sCanvas);
  const supChanged = changes();
  loopWhileVisible(sCanvas, dt => {
    if (!paused) supYaw += dt * 0.12;
    if (!supChanged(supYaw.toFixed(4)) && !supStale) return;
    supStale = false;
    const wide = innerWidth > 900;
    sup.clear();
    sup.draw(cube, mat.chain(mat.translate(wide ? sup.w / sup.h * 0.92 : 0, wide ? 0 : -0.8, 0), mat.rotateX(0.55), mat.rotateY(supYaw), mat.scale(wide ? 1 : 0.7)), { styles: ink, fog: 0.5 });
  });

  /* Headings: dimensioned by the page itself ------------------------------- */
  const dimmed = $$('main h1, main h2').filter(h => !h.closest('.step'));
  const dims = dimmed.map(h => {
    const d = document.createElement('span');
    d.className = 'dim js-only';
    d.setAttribute('aria-hidden', 'true');
    d.innerHTML = '<i></i><b></b><i></i>';
    h.after(d);
    return d;
  });
  const measure = () => dimmed.forEach((h, i) => {
    const w = Math.max(...$$('.ln > span', h).map(s => s.offsetWidth));
    dims[i].style.setProperty('--w', `${Math.min(w, h.clientWidth)}px`);
    $('b', dims[i]).textContent = w.toFixed(2);
  });
  new ResizeObserver(measure).observe(document.body);

  /* Reveal: transforms inside a clip, so text stays in the reading order ---- */
  const seen = new IntersectionObserver(entries => entries.forEach(e => {
    if (!e.isIntersecting) return;
    e.target.classList.add('in');
    seen.unobserve(e.target);
  }), { threshold: 0.2 });
  $$('.sheet, .numbers').forEach(el => seen.observe(el));

  /* Section rail ----------------------------------------------------------- */
  const rail = $$('.rail a');
  const byId = Object.fromEntries(rail.map(a => [a.hash.slice(1), a]));
  const spy = new IntersectionObserver(entries => entries.forEach(e => {
    if (!e.isIntersecting) return;
    rail.forEach(a => a.classList.remove('is-on'));
    byId[e.target.id].classList.add('is-on');
    document.body.classList.toggle('on-hot', e.target.id === 'support');
  }), { rootMargin: '-50% 0px -50% 0px' });
  Object.keys(byId).forEach(id => spy.observe(document.getElementById(id)));

  /* Clock ------------------------------------------------------------------ */
  const clock = $('[data-clock]');
  const fmt = new Intl.DateTimeFormat('en-GB', { timeZone: 'Australia/Perth', hour: '2-digit', minute: '2-digit' });
  const tick = () => { clock.textContent = fmt.format(new Date()); };
  tick();
  setInterval(tick, 30000);

  document.fonts.ready.then(() => {
    measure();
    gsap.matchMedia().add(PIN_ROOM, () => {
      pinned = true;
      anatomy.classList.add('is-pinned');
      ScrollTrigger.create({
        trigger: anatomy, start: 'top top', end: '+=240%', pin: true,
        onEnter: () => showStep(cur),
        onUpdate: self => showStep(Math.min(3, Math.floor(self.progress * 4))),
      });
      return () => {
        pinned = false;
        anatomy.classList.remove('is-pinned');
        gsap.set(steps, { clearProps: 'opacity,transform' });
        showStep(0);
      };
    });
    gsap.to('.progress', { scaleX: 1, ease: 'none', scrollTrigger: { start: 0, end: 'max', scrub: true } });
    if (!reduce) {
      gsap.fromTo('.film video', { scale: 1.16 }, { scale: 1, ease: 'none', scrollTrigger: { trigger: '.film', start: 'top bottom', end: 'bottom top', scrub: true } });
      gsap.from('.chart i', { scaleX: 0, duration: 1.3, stagger: 0.04, scrollTrigger: { trigger: '.chart', start: 'top 75%', once: true } });
    }
    ScrollTrigger.refresh();
  });
})();
