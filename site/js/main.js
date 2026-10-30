/* Site motion. Personality: premium, rigid metal. Signature ease expo.out, no overshoot.
   The page is complete without this file: motion decorates content that is already readable and reachable. */
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

  // Set up a below-the-fold element only as it nears the viewport, which keeps load light.
  function whenNear(el, setup) {
    const io = new IntersectionObserver(([e]) => {
      if (e.isIntersecting) { io.disconnect(); setup(); }
    }, { rootMargin: '100% 0px' });
    io.observe(el);
  }

  // True when any value differs from the last call, so a still canvas can skip its draw.
  function changes() {
    let last = '';
    return (...values) => {
      const key = values.join();
      if (key === last) return false;
      last = key;
      return true;
    };
  }

  /* Pause control: one state, shown by a button beside each moving element ---- */
  const PAUSE_KEY = 'catallax-motion';
  const stored = () => { try { return localStorage.getItem(PAUSE_KEY); } catch { return null; } };
  const remember = value => { try { localStorage.setItem(PAUSE_KEY, value); } catch { /* the page still works */ } };
  let paused = reduce || stored() === 'paused';
  const motionBtns = $$('.motion'), videos = $$('video'), nearVideos = new Set();
  const showPaused = () => motionBtns.forEach(b => { $('[data-motion-verb]', b).textContent = paused ? 'Resume' : 'Pause'; });
  motionBtns.forEach(btn => {
    btn.hidden = reduce;
    btn.addEventListener('click', () => {
      paused = !paused;
      remember(paused ? 'paused' : 'playing');
      showPaused();
      videos.forEach(v => { if (paused) v.pause(); else if (nearVideos.has(v)) v.play().catch(() => {}); });
    });
  });
  showPaused();

  /* Videos: start near the viewport, never offscreen. Reduced motion gets poster + controls. */
  videos.forEach(v => {
    if (reduce) return;
    v.controls = false;
    new IntersectionObserver(([e]) => {
      if (e.isIntersecting) { nearVideos.add(v); if (!paused) v.play().catch(() => {}); }
      else { nearVideos.delete(v); v.pause(); }
    }, { rootMargin: '25% 0px' }).observe(v);
  });

  /* Navigation menu under the desktop breakpoint ------------------------------ */
  const nav = $('.nav'), menu = $('.nav__toggle'), menuLinks = $('.nav__links');
  const setMenu = open => {
    nav.classList.toggle('is-open', open);
    menu.setAttribute('aria-expanded', open);
    menu.textContent = open ? 'Close' : 'Menu';
  };
  const menuOpen = () => nav.classList.contains('is-open');
  menu.addEventListener('click', () => setMenu(!menuOpen()));
  menuLinks.addEventListener('click', e => { if (e.target.closest('a')) setMenu(false); });
  document.addEventListener('pointerdown', e => { if (menuOpen() && !nav.contains(e.target)) setMenu(false); });
  addEventListener('keydown', e => { if (e.key === 'Escape' && menuOpen()) { setMenu(false); menu.focus(); } });

  /* Hero: draggable wireframe ------------------------------------------- */
  const heroCanvas = $('.hero__canvas');
  const hero = new Renderer(heroCanvas, { fov: 0.72 });
  // Lay the gear flat (axis Y) so the orbit never shows it edge-on.
  const gear = gen.gear();
  for (let i = 0; i < gear.v.length; i += 3) { const y = gear.v[i + 1]; gear.v[i + 1] = gear.v[i + 2]; gear.v[i + 2] = -y; }
  const shapes = { gear, torus: gen.torus(), lattice: gen.lattice() };
  const shapeNames = { gear: 'a spur gear', torus: 'a torus', lattice: 'a cube with a construction grid' };
  const S = { mesh: shapes.gear, yaw: 0.6, pitch: -0.5, vyaw: 0, vpitch: 0, progress: reduce ? 1 : 0, mx: 0, my: 0, x: 1.3, y: 0.05, k: 1.25 };
  const AUTO = 0.16; // rad/s
  let drag = null, heroStale = true, heroDrawn = 0;

  function layoutHero() {
    hero.resize();
    heroStale = true;
    const wide = innerWidth > 900;
    Object.assign(S, wide ? { x: hero.w / hero.h * 0.82, y: 0.1, k: 1.15 } : { x: 0, y: 0.8, k: 0.72 });
  }
  layoutHero();
  new ResizeObserver(layoutHero).observe(heroCanvas);

  heroCanvas.addEventListener('pointerdown', e => {
    if (drag) return; // a second finger never steals the orbit
    drag = { id: e.pointerId, x: e.clientX, y: e.clientY };
    heroCanvas.setPointerCapture(e.pointerId);
  });
  heroCanvas.addEventListener('pointermove', e => {
    const r = heroCanvas.getBoundingClientRect();
    S.mx = (e.clientX - r.left) / r.width - 0.5;
    S.my = (e.clientY - r.top) / r.height - 0.5;
    if (!drag || e.pointerId !== drag.id) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    drag.x = e.clientX; drag.y = e.clientY;
    S.vyaw = dx * 0.008;
    S.vpitch = dy * 0.006;
    S.yaw += S.vyaw;
    S.pitch = Math.max(-1.4, Math.min(1.4, S.pitch + S.vpitch));
  });
  const endDrag = e => {
    if (!drag || (e.pointerId !== undefined && e.pointerId !== drag.id)) return;
    drag = null;
    if (reduce) S.vyaw = S.vpitch = 0;
  };
  ['pointerup', 'pointercancel', 'lostpointercapture'].forEach(type => heroCanvas.addEventListener(type, endDrag));
  addEventListener('blur', endDrag);

  const hud = Object.fromEntries($$('[data-hud]').map(el => [el.dataset.hud, el]));
  let fps = 60, hudT = 0;
  const deg = r => ((r * 180 / Math.PI) % 360 + 360) % 360;
  const pmx = { x: 0, y: 0 };
  const heroChanged = changes();

  loopWhileVisible(heroCanvas, (dt, now) => {
    if (!drag) {
      S.yaw += S.vyaw + (paused ? 0 : AUTO * dt);
      S.pitch = Math.max(-1.4, Math.min(1.4, S.pitch + S.vpitch));
      S.vyaw *= 0.94; S.vpitch *= 0.9;
    }
    // Hover parallax eases toward the pointer; skipped under reduced motion.
    if (!reduce) { pmx.x += (S.mx - pmx.x) * 0.06; pmx.y += (S.my - pmx.y) * 0.06; }
    if (dt) fps += (1 / dt - fps) * 0.05;
    // Slow ambient rotation needs no more than 30 fps; a still solid needs no draw at all.
    const busy = drag || S.progress < 1 || heroStale;
    if (!busy && now - heroDrawn < 32) return;
    const moved = heroChanged(S.yaw.toFixed(4), S.pitch.toFixed(4), S.progress, pmx.x.toFixed(3), pmx.y.toFixed(3), S.x, S.y, S.k, S.mesh.edgeCount);
    if (!moved && !heroStale) return;
    heroStale = false;
    heroDrawn = now;
    const model = mat.chain(mat.translate(S.x, S.y, 0), mat.rotateX(S.pitch + pmx.y * 0.25), mat.rotateY(S.yaw + pmx.x * 0.35), mat.scale(S.k));
    hero.clear();
    const out = hero.draw(S.mesh, model, { progress: S.progress, styles, head: P.oxide });
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

  // The solid builds edge by edge. The headline's entrance lives in CSS, so it runs from first paint.
  if (!reduce) gsap.to(S, { progress: 1, duration: 2.4, ease: 'power2.inOut', delay: 0.2 });

  const shapeBtns = $$('.shape');
  shapeBtns.forEach(btn => btn.addEventListener('click', () => {
    const next = shapes[btn.dataset.shape];
    if (next === S.mesh) return;
    shapeBtns.forEach(b => { b.classList.toggle('is-active', b === btn); b.setAttribute('aria-pressed', b === btn); });
    heroCanvas.setAttribute('aria-label', `Rotating wireframe of ${shapeNames[btn.dataset.shape]}. Drag to orbit.`);
    gsap.killTweensOf(S, 'progress');
    if (reduce) { S.mesh = next; S.progress = 1; heroStale = true; return; }
    gsap.timeline()
      .to(S, { progress: 0, duration: 0.45, ease: 'power2.in' })
      .add(() => { S.mesh = next; })
      .to(S, { progress: 1, duration: 1.6, ease: 'power2.inOut' });
  }));

  /* Kernel anatomy: one canvas, four layers ------------------------------ */
  const anatomy = $('.anatomy');
  const aCanvas = $('.anatomy__canvas');
  const ana = new Renderer(aCanvas, { fov: 0.7 });
  ana.view = mat.translate(0, 0, -4.1);
  let anaStale = true;
  new ResizeObserver(() => { ana.resize(); anaStale = true; }).observe(aCanvas);

  // The section pins only when its text fits. Shorter viewports and high zoom read it as a plain list.
  const PIN_ROOM = '(min-height: 640px)';
  anatomy.classList.toggle('is-pinned', matchMedia(PIN_ROOM).matches);

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
  let cur = 0, rot = 0.6, aT = 0, anaStarted = reduce, pinned = false;
  const anaChanged = changes();

  loopWhileVisible(aCanvas, dt => {
    if (!anaStarted) showStep(0);   // build layer 01 as the section scrolls in
    if (!paused) { rot += dt * 0.3; aT += dt; }
    if (!paused && layers[1].alpha > 0) {   // let two control points breathe
      ctrl[1][2][1] = 0.55 * Math.sin(aT * 1.1) + 0.1;
      ctrl[2][1][1] = -0.45 * Math.sin(aT * 0.8 + 1);
      gen.patch(ctrl, 14, nurbs);
    }
    const moved = anaChanged(rot, aT, ...layers.flatMap(l => [l.alpha, l.progress]));
    if (!moved && !anaStale) return;
    anaStale = false;
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
      if (pinned) { gsap.set(out, { opacity: 0 }); gsap.set(inn, { opacity: 1, y: 0 }); }
      return;
    }
    if (from !== i) gsap.to(layers[from], { alpha: 0, duration: 0.45, ease: 'power2.in', overwrite: true });
    gsap.fromTo(layers[i], { alpha: 1, progress: 0 }, { progress: 1, duration: 1.4, ease: 'power2.inOut', overwrite: true });
    if (!pinned) return;
    // Opacity, not visibility: every layer stays in the reading order for screen readers.
    if (out !== inn) gsap.to(out, { opacity: 0, y: -28, duration: 0.35, ease: 'power2.in', overwrite: true });
    gsap.fromTo(inn, { opacity: 0, y: 36 }, { opacity: 1, y: 0, duration: 0.8, delay: from !== i ? 0.18 : 0, overwrite: true });
    gsap.fromTo($$('.step__specs li', inn), { opacity: 0, y: 10 }, { opacity: 1, y: 0, duration: 0.5, stagger: 0.06, delay: 0.35 });
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

  // Highlight at load: it changes the block height, and a shift later would move anchor targets mid-scroll.
  $$('code[data-rust]').forEach(code => {
    code.innerHTML = code.textContent.split('\n').map(l => `<span class="ln">${highlight(l)}</span>`).join('');
  });

  // Code windows type themselves out line by line as they arrive, then park a caret.
  if (!reduce) $$('.window').forEach(win => whenNear(win, () => {
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
  }));

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
  const nextTask = () => new Promise(resolve => { setTimeout(resolve); });

  document.fonts.ready.then(async () => {
    // Preserve the four-layer scroller. matchMedia reverts it if the viewport loses the room for it.
    gsap.matchMedia().add(PIN_ROOM, () => {
      pinned = true;
      anatomy.classList.add('is-pinned');
      gsap.to('.anatomy__rail i', {
        scaleY: 1, ease: 'none',
        scrollTrigger: {
          trigger: '.anatomy', start: 'top top', end: '+=300%', pin: true, scrub: true,
          onEnter: () => showStep(cur),
          onUpdate: self => showStep(Math.min(3, Math.floor(self.progress * 4))),
        },
      });
      return () => {
        pinned = false;
        anatomy.classList.remove('is-pinned');
        gsap.set([...steps, ...$$('.step__specs li')], { clearProps: 'opacity,transform' });
        showStep(0);
      };
    });

    if (reduce) {
      $$('[data-count]').forEach(el => { el.textContent = (+el.dataset.count).toLocaleString('en-US', { minimumFractionDigits: +(el.dataset.dec || 0) }); });
      ScrollTrigger.refresh();
      return;
    }
    await nextTask();

    // Hero parallax out: copy lifts faster than the solid.
    gsap.to('.hero__copy', { yPercent: -22, opacity: 0.15, ease: 'none', scrollTrigger: { trigger: '.hero', start: 'top top', end: 'bottom top', scrub: true } });
    gsap.to(heroCanvas, { yPercent: 18, ease: 'none', scrollTrigger: { trigger: '.hero', start: 'top top', end: 'bottom top', scrub: true } });

    // Progress bar and hide-on-scroll nav.
    gsap.to('.progress', { scaleX: 1, ease: 'none', scrollTrigger: { start: 0, end: 'max', scrub: true } });
    ScrollTrigger.create({ start: 0, end: 'max', onUpdate: s => { if (!menuOpen()) nav.classList.toggle('is-hidden', s.direction === 1 && s.scroll() > 240); } });

    // Manifesto: text brightens as it crosses the reading line. A mask leaves the DOM, link and reading order alone.
    const manifesto = $('.manifesto__text');
    manifesto.classList.add('is-scrubbed');
    gsap.fromTo(manifesto, { '--reveal': '0%' }, { '--reveal': '110%', ease: 'none', scrollTrigger: { trigger: manifesto, start: 'top 80%', end: 'bottom 50%', scrub: true } });
    await nextTask();

    // Section headings: line mask reveal, split only as each nears the viewport, re-split on resize.
    $$('[data-lines]').forEach(el => whenNear(el, () => SplitText.create(el, {
      type: 'lines', mask: 'lines', linesClass: 'split-line', autoSplit: true,
      onSplit: self => gsap.from(self.lines, { yPercent: 135, duration: 1.1, stagger: 0.08, scrollTrigger: { trigger: el, start: 'top 85%', once: true } }),
    })));

    // Assembly video: slow push-in parallax, copy drifts against it.
    gsap.fromTo('.assembly__media video', { scale: 1.2, yPercent: -5 }, { scale: 1, yPercent: 5, ease: 'none', scrollTrigger: { trigger: '.assembly', start: 'top bottom', end: 'bottom top', scrub: true } });
    gsap.fromTo('.assembly__copy', { y: 80 }, { y: -40, ease: 'none', scrollTrigger: { trigger: '.assembly', start: 'top bottom', end: 'bottom top', scrub: true } });

    // Card art parallax inside the tilt cards.
    $$('.card__art').forEach(svg => gsap.fromTo(svg, { y: 24 }, { y: -24, ease: 'none', scrollTrigger: { trigger: svg, start: 'top bottom', end: 'bottom top', scrub: true } }));
    await nextTask();

    // Generic entrance: rise 28px in batches. Opacity, not visibility, so hidden copy stays focusable and readable by assistive tech.
    const RISE = '.why__copy > p, .why__copy .btn, .why__terms .label, .why__list li, .thanks__copy > p, .thanks__sources .label, .thanks__list li, .works__list li, .works__not, .support__works > .label, .manifesto__sig, .stat, .code__lede, .code__notes li, .card, .bench__lede, .bench__row, .bench__foot, .ci__copy > p:not(.label), .ci__list li, .ci__term, .contact__lede, .contact__ctas, .assembly__note';
    const below = $$(RISE).filter(el => el.getBoundingClientRect().top > innerHeight);
    gsap.set(below, { y: 28, opacity: 0 });
    const riser = new IntersectionObserver(entries => {
      const batch = entries.filter(e => e.isIntersecting).map(e => e.target);
      if (!batch.length) return;
      batch.forEach(el => riser.unobserve(el));
      gsap.to(batch, { y: 0, opacity: 1, stagger: 0.07, duration: 0.9, overwrite: 'auto' });
    }, { rootMargin: '0px 0px -10% 0px' });
    below.forEach(el => riser.observe(el));
    // A keyboard user can reach a link before its block has risen.
    addEventListener('focusin', e => {
      const el = e.target.closest(RISE);
      if (!el || !below.includes(el)) return;
      riser.unobserve(el);
      gsap.to(el, { y: 0, opacity: 1, duration: 0.3, overwrite: 'auto' });
    });

    // Counters.
    $$('[data-count]').forEach(el => {
      const end = +el.dataset.count, dec = +(el.dataset.dec || 0), o = { v: 0 };
      const show = () => { el.textContent = o.v.toLocaleString('en-US', { minimumFractionDigits: dec, maximumFractionDigits: dec }); };
      gsap.to(o, { v: end, duration: 2.2, ease: 'expo.out', onUpdate: show, scrollTrigger: { trigger: el, start: 'top 88%', once: true } });
    });

    // Benchmark bars grow from zero.
    gsap.from('.bench .bar', { scaleX: 0, duration: 1.4, stagger: 0.03, scrollTrigger: { trigger: '.bench__table', start: 'top 75%', once: true } });

    ScrollTrigger.refresh();
  });
})();
