/* The page measures its own heading. Everything else on these sheets is static. */
(() => {
  'use strict';
  const dim = document.querySelector('.dim');
  const lines = [...document.querySelectorAll('h1 .ln > span')];
  const measure = () => {
    const w = Math.max(...lines.map(s => s.offsetWidth));
    dim.style.setProperty('--w', `${Math.min(w, dim.parentElement.clientWidth)}px`);
    dim.querySelector('b').textContent = w.toFixed(2);
  };
  document.fonts.ready.then(measure);
  new ResizeObserver(measure).observe(document.body);

})();
