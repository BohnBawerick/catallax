#!/usr/bin/env python3
"""Import a finished three-kernel comparison, then generate the compare page.

`--import-run` projects a `catallax-kernel-comparison/1` run from
data/kernel-compare-runs into site/kernel-compare.json. Rendering reads only
that projection. Never hand-edit the projection or the page.
"""
import argparse
from collections import Counter
import hashlib
from html import escape
import json
from pathlib import Path

from pages import REPO, chrome, heading

ROOT = Path(__file__).resolve().parents[2]
SITE = ROOT/'site'
PROJECTION = SITE/'kernel-compare.json'
ENGINES = (('catallax', 'Catallax'), ('occt', 'OpenCASCADE'), ('smlib', 'SMLib'))
OUTCOMES = {'done-sampled': 'D', 'not-done': 'N', 'refused': 'R', 'unsolved-time-limit': 'T', 'unresolved': 'U'}
LABELS = {'D': 'Done', 'N': 'Not done', 'R': 'Refused', 'T': 'Over time limit', 'U': 'Unresolved'}
CORPORA = {'fillets': 'Fillet corpus', 'five-tool': 'Five-tool targets', 'comparison-fixtures': 'Comparison fixtures',
           'comparison-measure': 'Representative recipes'}
DEFAULT_TOLERANCE = 0.001


def reason(engine, solve, check, code):
    """A short, factual reason for any outcome other than done."""
    if code == 'R':
        refusal = solve.get('refusal') or solve.get('native') or {}
        return refusal.get('code') or solve.get('detail', '').split(': ', 1)[-1]
    if code == 'U':
        oracle = check.get('oracle')
        if oracle:
            return 'reference not available: ' + oracle['detail'].split(':', 1)[0]
        if check.get('distance_upper_mm') is None:
            return 'no sampled distance bound'
        return 'sampled bound not settled'
    if code == 'N':
        if check.get('valid') is False:
            return 'invalid raw STEP'
        failures = len(check.get('witness_failures') or [])
        if failures:
            return f'{failures} material witness{"es" if failures != 1 else ""} failed'
        return 'sampled distance outside tolerance'
    return None


def import_run(source, destination, ref):
    data = json.loads(source.read_text())
    if data.get('format') != 'catallax-kernel-comparison/1' or data.get('status') != 'completed':
        raise SystemExit('expected a completed catallax-kernel-comparison/1 run')
    tolerances = data['policy']['tolerances_mm']
    if DEFAULT_TOLERANCE not in tolerances:
        raise SystemExit(f'run does not sweep the site tolerance {DEFAULT_TOLERANCE} mm')
    cases, volumes = [], []
    for case in data['cases']:
        row = {'id': case['id'], 'corpus': case['corpus'], 'engines': {}}
        volume = {}
        for engine, _ in ENGINES:
            result = case['engines'][engine]
            sweep = {s['tolerance_mm']: s for s in result['sweep']}
            if sorted(sweep) != sorted(tolerances):
                raise SystemExit(f'{case["id"]} {engine}: sweep does not cover every tolerance')
            unknown = {s['outcome'] for s in sweep.values()} - set(OUTCOMES)
            if unknown:
                raise SystemExit(f'{case["id"]} {engine}: unknown outcome {unknown}')
            codes = ''.join(OUTCOMES[sweep[t]['outcome']] for t in tolerances)
            check = result['check']
            row['engines'][engine] = {
                'outcomes': codes,
                'upper_mm': check.get('distance_upper_mm'),
                'reasons': {c: reason(engine, result['solve'], check, c) for c in sorted(set(codes) - {'D'})},
            }
            if 'analytic_volume_check' in check:
                volume[engine] = check['analytic_volume_check']
        cases.append(row)
        if volume:
            first = next(iter(volume.values()))
            volumes.append({'id': case['id'], 'expected_mm3': first['expected_mm3'], 'method': first['method'],
                            'error_mm3': {e: volume[e]['error_mm3'] if e in volume else None for e, _ in ENGINES}})
    projection = {
        'format': 'catallax-site-kernel-compare/1',
        'source': {'path': 'data/kernel-compare-runs/' + source.name, 'ref': ref,
                   'sha256': hashlib.sha256(source.read_bytes()).hexdigest()},
        'source_tag': data['source_tag'], 'source_revision': data['source_revision'],
        'harness_revision': data['harness_revision'], 'raw_result_sha256': data['raw_result_sha256'],
        'started_utc': data['started_utc'], 'completed_utc': data['completed_utc'],
        'machine_class': data['machine']['class'], 'attempts': data['executed_producer_attempts'],
        'policy': data['policy'],
        'tools': {name: info['version'] for name, info in data['observation_runs'][-1]['tools']['packages'].items()},
        'checker': next(c['engines']['occt']['check']['checker'] for c in data['cases'] if c['engines']['occt']['check'].get('checker')),
        'tolerances_mm': tolerances, 'cases': cases, 'volume_checks': volumes,
        'summary': {str(t): data['summary'][str(t)]['outcomes'] for t in tolerances},
    }
    verify(projection)
    destination.write_text(json.dumps(projection, indent=1, sort_keys=True) + '\n')


def counts(projection, index):
    tally = {e: Counter() for e, _ in ENGINES}
    for case in projection['cases']:
        for engine, _ in ENGINES:
            tally[engine][case['engines'][engine]['outcomes'][index]] += 1
    return tally


def verify(projection):
    """The projection's recount must equal the run's own published summary."""
    names = {v: k for k, v in OUTCOMES.items()}
    for i, tolerance in enumerate(projection['tolerances_mm']):
        key = str(tolerance)
        tally = counts(projection, i)
        for engine, _ in ENGINES:
            ours = {names[c]: n for c, n in tally[engine].items()}
            if ours != projection['summary'][key][engine]:
                raise SystemExit(f'recount differs from the run summary at {key} mm for {engine}')


def tol(value):
    """A tolerance in plain decimals, as the slider prints it."""
    return format(value, 'f').rstrip('0')


def mm(value):
    if value is None:
        return '—'
    return f'{value:.3g}' if value else '0'


def render(projection_path):
    p = json.loads(projection_path.read_text())
    verify(p)
    tolerances, cases = p['tolerances_mm'], p['cases']
    at = tolerances.index(DEFAULT_TOLERANCE)
    total, limit = len(cases), p['policy']['solve_limit_seconds']
    tally = counts(p, at)
    data_link = f'{REPO}/blob/{p["source"]["ref"]}/{p["source"]["path"]}'

    summary_rows = ''.join(
        f'<tr><th scope="row">{name}</th>' + ''.join(f'<td>{tally[engine][c]}</td>' for c in 'DNRTU') + '</tr>'
        for engine, name in ENGINES)
    def tallied(corpus, engine, code):
        return sum(c['corpus'] == corpus and c['engines'][engine]['outcomes'][at] == code for c in cases)
    corpus_rows = ''.join(
        f'<tr><th scope="row">{CORPORA[corpus]}</th><td>{sum(c["corpus"] == corpus for c in cases)}</td>' +
        ''.join(f'<td>{tallied(corpus, e, "D")}<small>{tallied(corpus, e, "U")} unresolved</small></td>' for e, _ in ENGINES) + '</tr>'
        for corpus in CORPORA)
    options = ''.join(f'<option value="{i}" label="{tol(t)} mm"></option>' for i, t in enumerate(tolerances))
    cells = ''.join(f'''
        <div class="demo__row" data-engine="{engine}"><b>{name}</b><div class="demo__cells" aria-hidden="true">{''.join(f'<i class="o-{c["engines"][engine]["outcomes"][at]}" title="{escape(c["id"])}"></i>' for c in cases)}</div><p class="demo__count"><span>{tally[engine]["D"]} / {total}</span><small>Done</small></p></div>''' for engine, name in ENGINES)
    outcomes = {engine: [c['engines'][engine]['outcomes'] for c in cases] for engine, _ in ENGINES}

    def short(case_id):
        return case_id.replace('comparison/crates/catallax-features/tests/fixtures/', 'fixtures/')

    def cell(case, engine):
        e = case['engines'][engine]
        code = e['outcomes'][at]
        note = e['reasons'].get(code)
        if code == 'D':
            note = f'≤ {mm(e["upper_mm"])} mm sampled'
        return f'<td class="o o-{code}"><span>{LABELS[code]}</span>{f"<small>{escape(note)}</small>" if note else ""}</td>'

    def link(case_id):
        if case_id.startswith('fillets/'):
            name = case_id.split('/', 1)[1]
            return f'<a href="gallery/{name}.html">{escape(short(case_id))}</a>'
        return escape(short(case_id))
    case_rows = ''.join(f'<tr><th scope="row">{link(c["id"])}</th>' + ''.join(cell(c, e) for e, _ in ENGINES) + '</tr>' for c in cases)

    def error(value):
        return 'not answered' if value is None else f'{value:+.3g}'
    volume_rows = ''.join(
        f'<tr><th scope="row">{link(v["id"])}</th><td>{v["expected_mm3"]:.10g}</td>' +
        ''.join(f'<td>{error(v["error_mm3"][e])}</td>' for e, _ in ENGINES) + '</tr>'
        for v in p['volume_checks'])

    body = f'''<main id="main" class="read ev">
  <header class="read__head">
{heading(['One definition', 'of <em>done</em>.'])}
  </header>
  <p class="ev__lede">Catallax {escape(p["source_tag"])}, OpenCASCADE and SMLib, run on the same {total} cases and scored by the same rule.</p>

  <section class="done" aria-labelledby="rule">
    <h2 class="ev__h" id="rule">The rule</h2>
    <p class="done__rule">A case is done when the kernel returns one valid solid within <var>T</var> = {limit:g} s, its authored material witnesses pass, and every sampled boundary distance to the reference is at most <var>δ</var> = {tol(DEFAULT_TOLERANCE)} mm. Anything else is not done.</p>
    <dl class="done__terms">
      <div><dt class="mono">Tolerance <var>δ</var></dt><dd>{tol(DEFAULT_TOLERANCE)} mm. The run also records {", ".join(tol(t) for t in tolerances if t != DEFAULT_TOLERANCE)} mm; slide below.</dd></div>
      <div><dt class="mono">Time limit <var>T</var></dt><dd>{limit:g} s for the whole process, including start-up and STEP export. Over the limit is unsolved.</dd></div>
      <div><dt class="mono">Unresolved</dt><dd>The check could not settle a bound, usually because the reference does not define that case yet. Not a loss for the kernel.</dd></div>
      <div><dt class="mono">Applied uniformly</dt><dd>Same cases, same <var>δ</var>, same <var>T</var> for all three kernels. One run per kernel, no retries.</dd></div>
    </dl>
  </section>

  <section aria-labelledby="counts">
    <h2 class="ev__h" id="counts">Results at <var>δ</var> = {tol(DEFAULT_TOLERANCE)} mm</h2>
    <div class="table-wrap"><table class="ev__table ev__table--num">
      <thead><tr><th scope="col">Kernel</th><th scope="col">Done</th><th scope="col">Not done</th><th scope="col">Refused</th><th scope="col">Over <var>T</var></th><th scope="col">Unresolved</th></tr></thead>
      <tbody>{summary_rows}</tbody>
    </table></div>
    <p>Out of {total} cases each. Done means sampled done, not proven correct everywhere.</p>
  </section>

  <section aria-labelledby="slide">
    <h2 class="ev__h" id="slide">Slide the tolerance</h2>
    <figure class="demo demo--real" data-sweep>
      <p class="demo__flag"><strong>Recorded data</strong><span>Each kernel's unchanged result, checked at every tolerance in the run. One square per case.</span></p>
      <div class="demo__control js-only">
        <label for="tolerance"><span class="mono">Tolerance <var>δ</var></span><output for="tolerance">{tol(DEFAULT_TOLERANCE)} mm</output></label>
        <input id="tolerance" type="range" min="0" max="{len(tolerances) - 1}" step="1" value="{at}" list="tolerances">
        <datalist id="tolerances">{options}</datalist>
        <div class="demo__scale" aria-hidden="true"><span>{tol(tolerances[0])} mm · loose</span><span>{tol(tolerances[-1])} mm · strict</span></div>
      </div>
      <div class="demo__rows">{cells}
      </div>
      <figcaption class="demo__key mono"><span><i class="o-D"></i>Done</span><span><i class="o-N"></i>Not done</span><span><i class="o-R"></i>Refused</span><span><i class="o-T"></i>Over <var>T</var></span><span><i class="o-U"></i>Unresolved</span></figcaption>
    </figure>
    <script type="application/json" id="sweep">{json.dumps({'tolerances': tolerances, 'outcomes': outcomes}, separators=(',', ':'))}</script>
  </section>

  <section aria-labelledby="corpora">
    <h2 class="ev__h" id="corpora">Done, by corpus</h2>
    <p>Most unresolved fillet cases wait on the reference: the field does not yet define those fillets, so no kernel's answer can be scored there.</p>
    <div class="table-wrap"><table class="ev__table ev__table--num">
      <thead><tr><th scope="col">Corpus</th><th scope="col">Cases</th>{"".join(f'<th scope="col">{n}</th>' for _, n in ENGINES)}</tr></thead>
      <tbody>{corpus_rows}</tbody>
    </table></div>
  </section>

  <section aria-labelledby="volumes">
    <h2 class="ev__h" id="volumes">Closed-form volume checks</h2>
    <p>Fillet cases with an independently derived exact volume. Volume error is a separate diagnostic, not part of the rule above.</p>
    <div class="table-wrap"><table class="ev__table ev__table--num">
      <thead><tr><th scope="col">Case</th><th scope="col">Exact volume, mm³</th>{"".join(f'<th scope="col">{n} error, mm³</th>' for _, n in ENGINES)}</tr></thead>
      <tbody>{volume_rows}</tbody>
    </table></div>
  </section>

  <section aria-labelledby="cases">
    <h2 class="ev__h" id="cases">Every case at <var>δ</var> = {tol(DEFAULT_TOLERANCE)} mm</h2>
    <details class="code"><summary class="mono">Show all {total} cases</summary>
    <div class="table-wrap"><table class="ev__table ev__cases">
      <thead><tr><th scope="col">Case</th>{"".join(f'<th scope="col">{n}</th>' for _, n in ENGINES)}</tr></thead>
      <tbody>{case_rows}</tbody>
    </table></div>
    </details>
  </section>

  <section aria-labelledby="scope">
    <h2 class="ev__h" id="scope">What these numbers are</h2>
    <div class="scope">
      <p>Sampled checks, not a certificate over the whole part. The reference is Catallax's field definition, compiled independently of its B-rep but written by the same authors. OpenCASCADE both competes and imports every STEP file for checking.</p>
      <p>Run on a shared development host. Solve times are only checked against <var>T</var>; they are not a speed ranking.</p>
      <p>Catallax <a href="{REPO}/releases/tag/{escape(p["source_tag"])}" target="_blank" rel="noopener">{escape(p["source_tag"])}</a> at <code>{p["source_revision"][:10]}</code>. OpenCASCADE via cadquery-ocp {escape(p["tools"]["cadquery-ocp"])}. SMLib via smlibpy {escape(p["tools"]["smlibpy"])}. {p["attempts"]} attempts, {escape(p["started_utc"][:10])}. <a href="{data_link}" target="_blank" rel="noopener">Run data</a> · <a href="kernel-compare.json">site projection</a> · <a href="{REPO}/tree/{p["harness_revision"]}/tools/kernel-compare" target="_blank" rel="noopener">harness</a>.</p>
    </div>
  </section>
</main>'''
    script = '\n<script src="compare.js"></script>'
    return chrome('Compare', f'Catallax {p["source_tag"]}, OpenCASCADE and SMLib on {total} cases, scored by one definition of done.',
                  'compare.html', 0, 'compare.html', body,
                  ('Data', data_link, f'Run {p["started_utc"][:10]}'), script)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--import-run', type=Path, help='completed run JSON from data/kernel-compare-runs')
    parser.add_argument('--source-ref', default='main', help='git ref that holds the imported run file')
    parser.add_argument('--projection', type=Path, default=PROJECTION)
    parser.add_argument('--output', type=Path, default=SITE/'compare.html')
    parser.add_argument('--check', action='store_true', help='fail if the page differs from the committed projection')
    args = parser.parse_args()
    if args.check and args.import_run:
        parser.error('--check cannot import a run')
    if args.import_run:
        import_run(args.import_run, args.projection, args.source_ref)
    page = render(args.projection)
    if args.check:
        if args.output.read_text() != page:
            raise SystemExit('compare.html is stale; run site/render/kernel_compare.py')
        print('compare.html matches the committed comparison projection')
    else:
        args.output.write_text(page)
