/* Site motion. Personality: premium, rigid metal. Signature ease expo.out, no overshoot. */
(() => {
  'use strict';
  const reduce = matchMedia('(prefers-reduced-motion: reduce)').matches;
  const fine = matchMedia('(hover: hover) and (pointer: fine)').matches;
  const { mat, gen, Renderer, palette, defaultStyles } = Wire;
  const P = palette();
  const styles = defaultStyles(P);
  const $ = (s, el = document) => el.querySelector(s);
  const $$ = (s, el = document) => [...el.querySelectorAll(s)];

  gsap.registerPlugin(ScrollTrigger, SplitText);
  gsap.defaults({ ease: 'expo.out', duration: 0.9 });

  // Run a canvas loop only while its element is on screen.
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

  /* Hero: draggable wireframe ------------------------------------------- */
  const heroCanvas = $('.hero__canvas');
  const hero = new Renderer(heroCanvas, { fov: 0.72 });
  // Lay the gear flat (axis Y) so the orbit never shows it edge-on.
  const gear = gen.gear();
  for (let i = 0; i < gear.v.length; i += 3) { const y = gear.v[i + 1]; gear.v[i + 1] = gear.v[i + 2]; gear.v[i + 2] = -y; }
  const shapes = { gear, torus: gen.torus(), lattice: gen.lattice() };
  const S = { mesh: shapes.gear, yaw: 0.6, pitch: -0.5, vyaw: 0, vpitch: 0, progress: reduce ? 1 : 0, mx: 0, my: 0, x: 1.3, y: 0.05, k: 1.25 };
  const AUTO = reduce ? 0 : 0.16; // rad/s
  let drag = null;

  function layoutHero() {
    hero.resize();
    const wide = innerWidth > 900;
    Object.assign(S, wide ? { x: hero.w / hero.h * 0.82, y: 0.1, k: 1.15 } : { x: 0, y: 0.8, k: 0.72 });
  }
  layoutHero();
  new ResizeObserver(layoutHero).observe(heroCanvas);

  heroCanvas.addEventListener('pointerdown', e => {
    drag = { x: e.clientX, y: e.clientY };
    heroCanvas.setPointerCapture(e.pointerId);
  });
  heroCanvas.addEventListener('pointermove', e => {
    const r = heroCanvas.getBoundingClientRect();
    S.mx = (e.clientX - r.left) / r.width - 0.5;
    S.my = (e.clientY - r.top) / r.height - 0.5;
    if (!drag) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    drag = { x: e.clientX, y: e.clientY };
    S.vyaw = dx * 0.008;
    S.vpitch = dy * 0.006;
    S.yaw += S.vyaw;
    S.pitch = Math.max(-1.4, Math.min(1.4, S.pitch + S.vpitch));
  });
  const endDrag = () => { drag = null; if (reduce) S.vyaw = S.vpitch = 0; };
  heroCanvas.addEventListener('pointerup', endDrag);
  heroCanvas.addEventListener('pointercancel', endDrag);

  const hud = Object.fromEntries($$('[data-hud]').map(el => [el.dataset.hud, el]));
  let fps = 60, hudT = 0;
  const deg = r => ((r * 180 / Math.PI) % 360 + 360) % 360;
  const pmx = { x: 0, y: 0 };

  loopWhileVisible(heroCanvas, (dt, now) => {
    if (!drag) {
      S.yaw += S.vyaw + AUTO * dt;
      S.pitch = Math.max(-1.4, Math.min(1.4, S.pitch + S.vpitch));
      S.vyaw *= 0.94; S.vpitch *= 0.9;
    }
    // Hover parallax eases toward the pointer; skipped under reduced motion.
    if (!reduce) { pmx.x += (S.mx - pmx.x) * 0.06; pmx.y += (S.my - pmx.y) * 0.06; }
    const model = mat.chain(mat.translate(S.x, S.y, 0), mat.rotateX(S.pitch + pmx.y * 0.25), mat.rotateY(S.yaw + pmx.x * 0.35), mat.scale(S.k));
    hero.clear();
    const out = hero.draw(S.mesh, model, { progress: S.progress, styles, head: P.oxide });
    if (dt) fps += (1 / dt - fps) * 0.05;
    if (out && now - hudT > 120) {
      hudT = now;
      hud.yaw.textContent = deg(S.yaw).toFixed(1).padStart(5, '0');
      hud.pitch.textContent = (S.pitch * 180 / Math.PI).toFixed(1);
      hud.v.textContent = S.mesh.vertexCount.toLocaleString('en-US');
      hud.e.textContent = S.mesh.edgeCount.toLocaleString('en-US');
      hud.fps.textContent = Math.round(Math.min(fps, 999));
      const m = mat.chain(hero.view, model);
      let txt = 'MODELVIEW\n';
      for (let r = 0; r < 4; r++) txt += [0, 4, 8, 12].map(c => m[c + r].toFixed(2).padStart(6)).join(' ') + '\n';
      hud.matrix.textContent = txt;
    }
  });

  const shapeBtns = $$('.shape');
  shapeBtns.forEach(btn => btn.addEventListener('click', () => {
    const next = shapes[btn.dataset.shape];
    if (next === S.mesh) return;
    shapeBtns.forEach(b => { b.classList.toggle('is-active', b === btn); b.setAttribute('aria-pressed', b === btn); });
    gsap.killTweensOf(S, 'progress');
    if (reduce) { S.mesh = next; S.progress = 1; return; }
    gsap.timeline()
      .to(S, { progress: 0, duration: 0.45, ease: 'power2.in' })
      .add(() => { S.mesh = next; })
      .to(S, { progress: 1, duration: 1.6, ease: 'power2.inOut' });
  }));

  /* Kernel anatomy: one canvas, four layers ------------------------------ */
  const aCanvas = $('.anatomy__canvas');
  const ana = new Renderer(aCanvas, { fov: 0.7 });
  ana.view = mat.translate(0, 0, -4.1);
  new ResizeObserver(() => ana.resize()).observe(aCanvas);

  const brep = gen.box(0.6, 0.6, 0.6);
  [0, 1, 2, 3].forEach(i => { brep.tag[i] = 2; });     // one face loop, in oxide
  const ctrl = Array.from({ length: 4 }, (_, i) => Array.from({ length: 4 }, (_, j) =>
    [(-1 + 2 * i / 3) * 1.05, 0.32 * Math.sin(i * 1.7) * Math.cos(j * 1.3) - 0.1, (-1 + 2 * j / 3) * 1.05]));
  const nurbs = gen.patch(ctrl);
  const bool = gen.merge(
    gen.box(0.9, 0.42, 0.9),
    gen.cylinder({ r: 0.44, y0: -1, y1: 1, seg: 48, rails: 8, tag: 3 }),
    gen.cylinder({ r: 0.44, y0: -0.42, y1: 0.42, seg: 48, rails: 16, tag: 2, railTag: 1 }));
  const sphere = gen.sphere({ r: 1 });
  const layers = [brep, nurbs, bool, sphere].map((mesh, i) => ({ mesh, alpha: i ? 0 : 1, progress: reduce ? 1 : 0 }));
  let cur = 0, rot = 0.6, aT = 0, anaStarted = reduce;

  loopWhileVisible(aCanvas, dt => {
    if (!anaStarted) showStep(0);   // build layer 01 as the section scrolls in
    if (!reduce) { rot += dt * 0.3; aT += dt; }
    if (!reduce && layers[1].alpha > 0) {   // let two control points breathe
      ctrl[1][2][1] = 0.55 * Math.sin(aT * 1.1) + 0.1;
      ctrl[2][1][1] = -0.45 * Math.sin(aT * 0.8 + 1);
      gen.patch(ctrl, 14, nurbs);
    }
    const model = mat.chain(mat.rotateX(-0.5), mat.rotateY(rot), mat.scale(1.05));
    ana.clear();
    for (const l of layers) ana.draw(l.mesh, model, { progress: l.progress, alpha: l.alpha, styles, head: P.oxide, dots: P.oxide });
  });

  const steps = $$('.step'), idx = $$('.anatomy__index li'), stepNum = $('[data-step-num]');
  idx[0].classList.add('is-active');
  function showStep(i) {
    if (i === cur && anaStarted) return;
    const from = cur;
    cur = i;
    anaStarted = true;
    idx.forEach((li, k) => li.classList.toggle('is-active', k === i));
    stepNum.textContent = String(i + 1).padStart(2, '0');
    const out = steps[from], inn = steps[i];
    if (reduce) {
      layers.forEach((layer, k) => { layer.alpha = k === i ? 1 : 0; layer.progress = 1; });
      gsap.set(out, { autoAlpha: 0 });
      gsap.set(inn, { autoAlpha: 1, y: 0 });
      gsap.set($$('.step__specs li', inn), { autoAlpha: 1, y: 0 });
      return;
    }
    if (from !== i) gsap.to(layers[from], { alpha: 0, duration: 0.45, ease: 'power2.in', overwrite: true });
    gsap.fromTo(layers[i], { alpha: 1, progress: 0 }, { progress: 1, duration: 1.4, ease: 'power2.inOut', overwrite: true });
    if (out !== inn) gsap.to(out, { autoAlpha: 0, y: -28, duration: 0.35, ease: 'power2.in', overwrite: true });
    gsap.fromTo(inn, { autoAlpha: 0, y: 36 }, { autoAlpha: 1, y: 0, duration: 0.8, delay: from !== i ? 0.18 : 0, overwrite: true });
    gsap.fromTo($$('.step__specs li', inn), { autoAlpha: 0, y: 10 }, { autoAlpha: 1, y: 0, duration: 0.5, stagger: 0.06, delay: 0.35 });
  }

  /* Code highlighting ------------------------------------------------------ */
  const RX = /(\/\/[^\n]*)|("(?:\\.|[^"\\])*")|(#!?\[[^\]]*\])|(\b\d[\d_]*(?:\.\d+)?\b)|(\b[a-z_][a-z0-9_]*!)|(\b(?:pub|struct|impl|fn|let|mut|const|use|for|in|as|type|if|else|return|where|self|Self|crate|mod|enum|match|loop|while)\b)|(\b(?:[A-Z][A-Za-z0-9_]*|f32|f64|u8|u32|u64|usize|i32|bool|f64x4)\b)|(\b[a-z_][a-z0-9_]*(?=\())/g;
  const CLS = [null, 'c-com', 'c-str', 'c-attr', 'c-num', 'c-mac', 'c-kw', 'c-ty', 'c-fn'];
  const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  const highlight = line => {
    let html = '', at = 0;
    for (const m of line.matchAll(RX)) {
      const g = m.findIndex((x, k) => k > 0 && x !== undefined);
      html += esc(line.slice(at, m.index)) + `<span class="${CLS[g]}">${esc(m[0])}</span>`;
      at = m.index + m[0].length;
    }
    return html + esc(line.slice(at));
  };
  $$('code[data-rust]').forEach(code => {
    code.innerHTML = code.textContent.split('\n').map(l => `<span class="ln">${highlight(l)}</span>`).join('');
  });

  /* Counters + benchmark bars (DOM prep, safe before fonts) --------------- */
  $$('.bench__row:not(.bench__row--head)').forEach(row => {
    const cells = $$('[data-v]', row), max = Math.max(...cells.map(c => +c.dataset.v));
    cells.forEach(c => {
      const bar = document.createElement('i');
      bar.className = 'bar';
      bar.style.setProperty('--w', (+c.dataset.v / max).toFixed(3));
      c.prepend(bar);
    });
  });

  /* Footer clock, Perth time --------------------------------------------- */
  const clock = $('[data-clock]');
  const fmt = new Intl.DateTimeFormat('en-GB', { timeZone: 'Australia/Perth', hour: '2-digit', minute: '2-digit' });
  const tickClock = () => { clock.textContent = fmt.format(new Date()); };
  tickClock();
  setInterval(tickClock, 30000);

  /* Videos: pause offscreen; reduced motion gets poster + controls -------- */
  $$('video').forEach(v => {
    if (reduce) { v.removeAttribute('autoplay'); v.pause(); v.controls = true; return; }
    new IntersectionObserver(([e]) => { if (e.isIntersecting) v.play().catch(() => {}); else v.pause(); }).observe(v);
  });

  /* Crosshair cursor ------------------------------------------------------- */
  if (fine) {
    const c = $('.cursor'), cx = $('.cursor__x'), cy = $('.cursor__y'), lab = $('.cursor__label');
    addEventListener('pointermove', e => {
      c.classList.add('is-on');
      cx.style.transform = `translateX(${e.clientX}px)`;
      cy.style.transform = `translateY(${e.clientY}px)`;
      lab.style.transform = `translate(${e.clientX + 12}px, ${e.clientY + 10}px)`;
      lab.textContent = `${String(e.clientX).padStart(4, '0')} · ${String(e.clientY).padStart(4, '0')}`;
    }, { passive: true });
    document.documentElement.addEventListener('pointerleave', () => c.classList.remove('is-on'));
  }

  /* Magnetic buttons + tilt cards (pointer only, motion allowed) ---------- */
  if (fine && !reduce) {
    $$('.magnetic').forEach(el => {
      const inner = el.firstElementChild;
      const q = (t, p) => gsap.quickTo(t, p, { duration: 0.5, ease: 'power3.out' });
      const x = q(el, 'x'), y = q(el, 'y'), ix = q(inner, 'x'), iy = q(inner, 'y');
      el.addEventListener('pointermove', e => {
        const r = el.getBoundingClientRect(), dx = e.clientX - r.left - r.width / 2, dy = e.clientY - r.top - r.height / 2;
        x(dx * 0.28); y(dy * 0.38); ix(dx * 0.12); iy(dy * 0.16);
      });
      el.addEventListener('pointerleave', () => { x(0); y(0); ix(0); iy(0); });
    });
    $$('[data-tilt]').forEach(card => {
      gsap.set(card, { transformPerspective: 1000 });
      const rx = gsap.quickTo(card, 'rotationX', { duration: 0.7, ease: 'power3.out' });
      const ry = gsap.quickTo(card, 'rotationY', { duration: 0.7, ease: 'power3.out' });
      card.addEventListener('pointermove', e => {
        const r = card.getBoundingClientRect(), px = (e.clientX - r.left) / r.width, py = (e.clientY - r.top) / r.height;
        ry((px - 0.5) * 10); rx((0.5 - py) * 8);
        card.style.setProperty('--mx', `${px * 100}%`);
        card.style.setProperty('--my', `${py * 100}%`);
      });
      card.addEventListener('pointerleave', () => { rx(0); ry(0); });
    });
  }

  /* Everything that measures text waits for fonts ------------------------ */
  document.fonts.ready.then(() => {
    // Preserve the four-layer scroller; reduced motion switches layers without scrubbing.
    gsap.to('.anatomy__rail i', {
      scaleY: 1, ease: 'none',
      scrollTrigger: {
        trigger: '.anatomy', start: 'top top', end: '+=300%', pin: true, scrub: !reduce,
        onEnter: () => showStep(cur),
        onUpdate: self => showStep(Math.min(3, Math.floor(self.progress * 4))),
      },
    });

    if (reduce) {
      $$('[data-count]').forEach(el => { el.textContent = (+el.dataset.count).toLocaleString('en-US', { minimumFractionDigits: +(el.dataset.dec || 0) }); });
      ScrollTrigger.refresh();
      return;
    }

    // Intro: headline lines rise from a mask, wireframe builds edge by edge.
    const title = SplitText.create('.hero__title .line', { type: 'lines', mask: 'lines' });
    gsap.timeline({ delay: 0.1 })
      .from('.nav > *', { y: -16, autoAlpha: 0, stagger: 0.06, duration: 0.7 }, 0)
      .from(title.lines, { yPercent: 108, duration: 1.2, stagger: 0.09 }, 0.15)
      .from('.hero__kicker span', { y: 12, autoAlpha: 0, stagger: 0.06, duration: 0.7 }, 0.25)
      .from('.hero__lede', { y: 24, autoAlpha: 0 }, 0.6)
      .from('.hero__ctas', { y: 20, autoAlpha: 0 }, 0.72)
      .from('.hud, .hero__frame i', { autoAlpha: 0, duration: 0.8, stagger: 0.05, ease: 'power2.out' }, 0.8)
      .to(S, { progress: 1, duration: 2.4, ease: 'power2.inOut' }, 0.2);

    // Hero parallax out: copy lifts faster than the solid.
    gsap.to('.hero__copy', { yPercent: -22, autoAlpha: 0.15, ease: 'none', scrollTrigger: { trigger: '.hero', start: 'top top', end: 'bottom top', scrub: true } });
    gsap.to(heroCanvas, { yPercent: 18, ease: 'none', scrollTrigger: { trigger: '.hero', start: 'top top', end: 'bottom top', scrub: true } });

    // Progress bar and hide-on-scroll nav.
    const nav = $('.nav');
    gsap.to('.progress', { scaleX: 1, ease: 'none', scrollTrigger: { start: 0, end: 'max', scrub: true } });
    ScrollTrigger.create({ start: 0, end: 'max', onUpdate: s => nav.classList.toggle('is-hidden', s.direction === 1 && s.scroll() > 240) });

    // Manifesto: words brighten as they scroll through the reading line.
    const words = SplitText.create('[data-words]', { type: 'words' });
    gsap.fromTo(words.words, { opacity: 0.14 }, {
      opacity: 1, stagger: 0.1, ease: 'none',
      scrollTrigger: { trigger: '[data-words]', start: 'top 80%', end: 'bottom 50%', scrub: true },
    });

    // Section headings: line mask reveal, re-split on resize.
    $$('[data-lines]').forEach(el => SplitText.create(el, {
      type: 'lines', mask: 'lines', autoSplit: true,
      onSplit: self => gsap.from(self.lines, { yPercent: 110, duration: 1.1, stagger: 0.08, scrollTrigger: { trigger: el, start: 'top 85%', once: true } }),
    }));

    // Assembly video: slow push-in parallax, copy drifts against it.
    gsap.fromTo('.assembly__media video', { scale: 1.2, yPercent: -5 }, { scale: 1, yPercent: 5, ease: 'none', scrollTrigger: { trigger: '.assembly', start: 'top bottom', end: 'bottom top', scrub: true } });
    gsap.fromTo('.assembly__copy', { y: 80 }, { y: -40, ease: 'none', scrollTrigger: { trigger: '.assembly', start: 'top bottom', end: 'bottom top', scrub: true } });

    // Card art parallax inside the tilt cards.
    $$('.card__art').forEach(svg => gsap.fromTo(svg, { y: 24 }, { y: -24, ease: 'none', scrollTrigger: { trigger: svg, start: 'top bottom', end: 'bottom top', scrub: true } }));

    // Generic entrance: rise 28px, stagger in batches.
    const rise = '.lineage__copy > p, .lineage__copy .btn, .lineage__lessons .label, .lineage__list li, .works__list li, .works__not, .support__works > .label, .manifesto__sig, .stat, .code__lede, .code__notes li, .card, .bench__lede, .bench__row, .bench__foot, .ci__copy > p:not(.label), .ci__list li, .ci__term, .contact__lede, .contact__ctas, .assembly__note';
    gsap.set(rise, { y: 28, autoAlpha: 0 });
    ScrollTrigger.batch(rise, { start: 'top 90%', once: true, onEnter: b => gsap.to(b, { y: 0, autoAlpha: 1, stagger: 0.07, duration: 0.9, overwrite: 'auto' }) });

    // Counters.
    $$('[data-count]').forEach(el => {
      const end = +el.dataset.count, dec = +(el.dataset.dec || 0), o = { v: 0 };
      const show = () => { el.textContent = o.v.toLocaleString('en-US', { minimumFractionDigits: dec, maximumFractionDigits: dec }); };
      gsap.to(o, { v: end, duration: 2.2, ease: 'expo.out', onUpdate: show, scrollTrigger: { trigger: el, start: 'top 88%', once: true } });
    });

    // Benchmark bars grow from zero.
    gsap.from('.bench .bar', { scaleX: 0, duration: 1.4, stagger: 0.03, scrollTrigger: { trigger: '.bench__table', start: 'top 75%', once: true } });

    // Code windows type themselves out, line by line, then park a caret.
    $$('.window').forEach(win => {
      const lines = $$('.ln', win);
      gsap.set(lines, { clipPath: 'inset(0 100% 0 0)' });
      const tl = gsap.timeline({ paused: true, onComplete: () => lines.at(-1).insertAdjacentHTML('beforeend', '<span class="caret"></span>') });
      let t = 0;
      lines.forEach(ln => {
        const n = ln.textContent.length, d = Math.min(0.32, 0.004 * n + 0.03);
        tl.to(ln, { clipPath: 'inset(0 0% 0 0)', duration: d, ease: `steps(${Math.max(1, Math.min(n, 48))})` }, t);
        t += d * 0.55;
      });
      ScrollTrigger.create({ trigger: win, start: 'top 78%', once: true, onEnter: () => tl.play() });
    });

    ScrollTrigger.refresh();
  });
})();
