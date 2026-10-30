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

  // The logo keeps the visitor's place on the main page when they came straight from it.
  // A sandboxed review frame can only reach the main page through history, so sheets
  // replace each other there and the main page stays one step back.
  const ref = document.referrer;
  const sandboxed = self.origin === 'null';
  if (sandboxed || (ref && /\/(index\.html)?$/.test(new URL(ref).pathname))) {
    document.querySelector('[data-home]').addEventListener('click', e => {
      e.preventDefault();
      history.back();
    });
  }
  if (sandboxed) {
    document.querySelectorAll('a[href="why.html"], a[href="thanks.html"]').forEach(a => a.addEventListener('click', e => {
      e.preventDefault();
      location.replace(a.getAttribute('href'));
    }));
  }
})();
