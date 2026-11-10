/* Tolerance slider over the recorded comparison. The page is complete at the default tolerance without this file. */
(() => {
  'use strict';
  const { tolerances, outcomes } = JSON.parse(document.getElementById('sweep').textContent);
  const slider = document.getElementById('tolerance');
  const shown = document.querySelector('[data-sweep] output');
  const rows = [...document.querySelectorAll('[data-sweep] .demo__row')].map(row => ({
    codes: outcomes[row.dataset.engine],
    cells: [...row.querySelectorAll('.demo__cells i')],
    count: row.querySelector('.demo__count span'),
  }));

  const draw = () => {
    const i = Number(slider.value), label = `${tolerances[i]} mm`;
    shown.textContent = label;
    slider.setAttribute('aria-valuetext', label);
    for (const row of rows) {
      let done = 0;
      row.codes.forEach((codes, n) => {
        done += codes[i] === 'D';
        row.cells[n].className = `o-${codes[i]}`;
      });
      row.count.textContent = `${done} / ${row.codes.length}`;
    }
  };
  slider.addEventListener('input', draw);
  draw();
})();
