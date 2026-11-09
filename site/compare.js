/* Definition-of-done demo. Synthetic placeholder data only, never measured results. */
(() => {
  'use strict';
  const CASES = 24;
  // One distribution for every row, so differences between rows mean nothing.
  const KERNELS = [['Catallax', 0x9e3779b9], ['OpenCASCADE', 0x85ebca6b], ['NVIDIA SMLib', 0xc2b2ae35]];

  const random = seed => () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  // Each demo case: no valid solid, over the time limit, or a geometric error in mm.
  const demo = seed => {
    const next = random(seed);
    return Array.from({ length: CASES }, () => {
      const r = next(), e = next();
      if (r < 0.08) return { kind: 'none' };
      if (r < 0.16) return { kind: 'slow' };
      return { kind: 'solid', error: 10 ** (-13 + 12 * e) };
    });
  };

  const root = document.querySelector('[data-demo]');
  const slider = root.querySelector('input[type="range"]');
  const shown = root.querySelector('output');
  const rows = root.querySelector('.demo__rows');
  const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
  const power = k => `10${String(k).replace(/./g, c => SUP[c])} mm`;

  const series = KERNELS.map(([name, seed]) => {
    const row = document.createElement('div');
    row.className = 'demo__row';
    row.innerHTML = `<b>${name}</b><div class="demo__cells" aria-hidden="true">${'<i></i>'.repeat(CASES)}</div><p class="demo__count"><span></span><small>Demo</small></p>`;
    rows.append(row);
    return { cases: demo(seed), cells: [...row.querySelectorAll('i')], count: row.querySelector('span') };
  });

  const draw = () => {
    const k = Number(slider.value), tolerance = 10 ** k;
    shown.textContent = power(k);
    slider.setAttribute('aria-valuetext', `${power(k)}, demo tolerance`);
    for (const s of series) {
      let done = 0;
      s.cases.forEach((c, i) => {
        const ok = c.kind === 'solid' && c.error <= tolerance;
        done += ok;
        s.cells[i].className = ok ? 'ok' : c.kind === 'solid' ? '' : c.kind;
      });
      s.count.textContent = `${done} / ${CASES}`;
    }
  };
  slider.addEventListener('input', draw);
  draw();
})();
