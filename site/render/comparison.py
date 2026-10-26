#!/usr/bin/env python3
"""Import a finished benchmark, then generate the static comparison page."""
import argparse
import hashlib
from html import escape
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / '.claude/skills/verify-catallax/scripts'
sys.path.insert(0, str(SCRIPTS))
from head_to_head import DATA, sha, summarize

ENGINES = ('catallax', 'solidcraft', 'cadquery', 'smlibpy', 'pali_ijo')


def render_index(source, index, check=False):
    data = json.loads(source.read_text())
    rows = summarize({'engines': ENGINES, 'pairs': data['observations']})
    parts = {row['part'] for row in rows}
    summary_rows = []
    for engine in ENGINES:
        subset = [row for row in rows if row['engine'] == engine]
        construction, edit = (
            sum(row['phase_outcomes']['build'].get('answered', 0)
                for row in subset if row['phase'] == phase)
            for phase in ('construction', 'edit'))
        total = sum(row['samples'] for row in subset if row['phase'] == 'construction')
        highlight = ' class="is-us"' if engine == 'catallax' else ''
        summary_rows.append(
            '    <div class="bench__row" role="row">\n'
            f'      <span role="cell">{escape(data["tools"][engine]["name"])}</span>\n'
            f'      <span role="cell"{highlight} data-v="{construction}"><b>{construction}</b></span>'
            f'<span role="cell"{highlight} data-v="{edit}"><b>{edit}</b></span>'
            f'<span role="cell" data-v="{total}"><b>{total}</b></span>\n'
            '    </div>')
    table = (
        '  <div class="bench__table" role="table" aria-label="Five-tool completed builds">\n'
        '    <div class="bench__row bench__row--head label" role="row">\n'
        '      <span role="columnheader">Tool</span><span role="columnheader">Construction</span>'
        '<span role="columnheader">Parameter edit</span><span role="columnheader">Samples per phase</span>\n'
        '    </div>\n' + '\n'.join(summary_rows) + '\n  </div>')
    paragraph = (
        '<p class="bench__foot" id="comparison-summary">'
        f'Completed builds across {len(parts)} targets and {data["timing"]["measured_pairs"]} measured samples per tool and target. '
        'An edit skipped after a failed construction remains incomplete. Completion is not correctness. '
        '<a href="comparison.html">The five-tool comparison</a> reports timings, correctness checks, versions and method.</p>')
    original = index.read_text()
    generated, count = re.subn(
        r'  <div class="bench__table"[^>]*>.*?</div>\n  <p class="bench__foot" id="comparison-summary">[^\n]*</p>',
        lambda _: table + '\n  ' + paragraph, original, flags=re.DOTALL)
    if count != 1:
        raise ValueError('index.html must have exactly one benchmark table and comparison-summary paragraph')
    if check:
        if generated != original:
            raise SystemExit('index.html comparison summary is stale; run site/render/comparison.py')
    else:
        index.write_text(generated)


def import_run(path, destination):
    report = json.loads(path.read_text())
    if report['status'] != 'completed' or set(report['engines']) != set(ENGINES):
        raise ValueError('publication requires a completed five-tool run')
    parts = {json.loads(p.read_text())['recipe'] for p in (DATA / 'cases').rglob('*.json')}
    if {pair['part'] for pair in report['pairs']} != parts:
        raise ValueError('publication requires every committed P0-P5 target')
    for part in parts:
        samples = [p for p in report['pairs'] if p['part'] == part and not p['warmup']]
        if len(samples) != report['timing']['measured_pairs'] or not samples:
            raise ValueError('missing measured samples')
        if any(set(p['engines']) != set(ENGINES) for p in samples):
            raise ValueError('missing engine observations')
    if report['host'].get('machine_class') != 'quiet benchmark machine':
        raise ValueError('publication timing must come from the designated quiet benchmark machine')
    revision = report['source_revision']
    def committed_bytes(path):
        relative = path.relative_to(ROOT).as_posix()
        return subprocess.check_output(['git', 'show', f'{revision}:{relative}'], cwd=ROOT)
    def committed_sha(path):
        return hashlib.sha256(committed_bytes(path)).hexdigest()
    config = json.loads(committed_bytes(DATA/'tools.json'))
    tools = {'catallax': {'name': 'Catallax', 'version': report['source_revision'],
                         'sha256': report['binary_sha256'], 'source': 'https://github.com/BohnBawerick/catallax',
                         'kernel': 'Catallax B-rep', 'license': 'MIT'},
             'solidcraft': {'name': 'SolidCraft', 'version': report['reference']['version'],
                           'sha256': report['reference']['binary_sha256'], 'source': report['reference']['source'],
                           'kernel': 'released SolidCraft', 'license': report['reference']['license']},
             **report['tools']}
    phase_keys = ('status', 'wall_ms', 'work', 'work_limit', 'operations_rebuilt',
                  'refusal_ms', 'refusal_upper_bound_ms', 'claimed', 'expected_faces')
    check_keys = ('status', 'oracle', 'oracle_source', 'numeric_comparison', 'property_coverage',
                  'claimed_properties', 'step', 'membership', 'step_sha256', 'binding', 'version')
    observations = []
    for pair in report['pairs']:
        observation = {key: pair[key] for key in ('part', 'sample', 'warmup', 'order', 'load_average_before', 'load_average_after')}
        observation['engines'] = {}
        for engine, result in pair['engines'].items():
            public = {'status': result['status'], 'variants': {}}
            public['edit_survival'] = {key: value for key, value in result.get('edit_survival', {}).items()
                                       if key in ('status', 'missing', 'unnamed', 'coverage', 'reason')}
            for label, variant in result.get('variants', {}).items():
                public_variant = {}
                for phase in ('build', 'display', 'measurement', 'step'):
                    public_variant[phase] = {key: value for key, value in variant.get(phase, {}).items() if key in phase_keys}
                public_variant['correctness'] = {key: value for key, value in variant.get('correctness', {}).items() if key in check_keys}
                public['variants'][label] = public_variant
            observation['engines'][engine] = public
        observations.append(observation)
    # No local paths, subprocess commands or private input names are published.
    result = {'format': 'catallax-site-comparison/1', 'raw_sha256': sha(path),
              'source_revision': report['source_revision'], 'tools': tools,
              'dart_version': config['dart']['version'], 'checker_version': config['checker']['version'],
              'host': report['host'], 'timing': report['timing'],
              'observations': observations,
              'installation_hashes': {key: value for key, value in report['installation'].items() if key.endswith('_sha256')},
              'load_average_range': [min(p[k][0] for p in report['pairs'] for k in ('load_average_before','load_average_after')),
                                     max(p[k][0] for p in report['pairs'] for k in ('load_average_before','load_average_after'))],
              'case_sha256': {p.relative_to(DATA).as_posix(): committed_sha(p) for p in sorted((DATA/'cases').rglob('*.json'))},
              'script_sha256': {p.name: committed_sha(p) for p in [SCRIPTS/'head_to_head.py', SCRIPTS/'external_tools.py',
                                  *[DATA/f for f in ('targets.py','python_tools.py','pali_parts.dart','analytic.py','check.py','tools.json','requirements.lock')]]}}
    destination.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')


def render(source, output):
    data = json.loads(source.read_text())
    rows = summarize({'engines': ENGINES, 'pairs': data['observations']})
    tools = data['tools']
    def text(value):
        return escape(str(value), quote=True)
    def tool(engine):
        value = tools[engine]
        return f'<a href="{text(value["source"])}">{text(value["name"])}</a>'
    def cell(value):
        return f'<td>{value}</td>'
    def ratio(row, yes, total):
        return 'not checked' if not row[total] else f'{row[yes]} / {row[total]}'
    parts = sorted({r['part'] for r in rows})
    values = [r['build_median_ms'] for r in rows if r['build_median_ms'] is not None and r['build_median_ms'] > 0]
    low = math.floor(math.log10(min(values))) if values else -2
    high = max(low+1, math.ceil(math.log10(max(values)))) if values else 4
    def bar(row):
        ms = row['build_median_ms']
        if ms is None:
            if row.get('not_expressible'):
                reason = 'not expressible'
            elif row['failed_rebuilds']:
                reason = 'build refused or failed'
            else:
                reason = 'incomplete'
            return f'<span class="unmeasured">{reason}; no successful median</span>'
        if ms <= 0 or not math.isfinite(ms):
            raise ValueError('successful timing must be finite and positive')
        width = 100*(math.log10(ms)-low)/(high-low)
        check = 'property / STEP checks passed' if (row['properties_checked'] == row['samples'] and row['membership_checked'] == row['samples'] and row['step_passed'] == row['samples'] and not row['wrong_properties'] and not row['wrong_membership']) else 'property / STEP checks failed or incomplete'
        return f'<span class="bar-track"><span class="bar-fill {text(row["engine"])}" style="width:{width:.6f}%"></span><span class="bar-number">{ms:.3f} ms</span></span><small>{check}</small>'
    options = ''.join(f'<option value="{text(p)}">{text(p)}</option>' for p in parts)
    charts = []
    for part in parts:
        group = [r for r in rows if r['part']==part]
        lines = []
        for engine in ENGINES:
            a = next(r for r in group if r['engine']==engine and r['phase']=='construction')
            b = next(r for r in group if r['engine']==engine and r['phase']=='edit')
            lines.append(f'<tr><th scope="row">{tool(engine)}</th>{cell(bar(a))}{cell(bar(b))}</tr>')
        ticks = ' · '.join(f'{10**power:g}' for power in range(low,high+1))
        charts.append(f'<section class="part-chart" data-part="{text(part)}"><h3>{text(part)}</h3><p class="axis">Log scale in milliseconds. Powers of ten: {ticks}. Shorter means less elapsed build time, not greater correctness.</p><div class="table-scroll" tabindex="0" role="region" aria-label="Build and edit timings for {text(part)}"><table class="timings"><thead><tr><th>Tool</th><th>Construction</th><th>Parameter edit</th></tr></thead><tbody>{"".join(lines)}</tbody></table></div></section>')
    checks = []
    for row in rows:
        check_cells = [text(row['part']), tool(row['engine']), text(row['phase']),
            f'{row["wrong_properties"]} / {row["properties_checked"]}',
            f'{row["wrong_membership"]} / {row.get("membership_points_checked",0)}',
            f'{row["step_opened"]} / {row["step_checked"]} / {row["step_passed"]}',
            ratio(row,'names_passed','names_checked') if row['phase']=='edit' else 'not applicable',
            str(row['failed_rebuilds']), str(row['incomplete']), str(row.get('not_expressible',0))]
        checks.append(f'<tr data-part="{text(row["part"])}">'+''.join(cell(c) for c in check_cells)+'</tr>')
    versions = []
    for engine in ENGINES:
        pin = tools[engine]
        version = pin['version'] + (' @ ' + pin['revision'][:12] if pin.get('revision') else '')
        versions.append('<tr>'+cell(tool(engine))+cell(text(version))+cell(text(pin['kernel']))+cell(f'<code>{text(pin.get("sha256",pin.get("archive_sha256","")))}</code>')+'</tr>')
    features = []
    for engine in ENGINES:
        cells = [tool(engine)]
        for family in ('P0/','P1-','P2-','P3-','P4-','P5-'):
            subset = [r for r in rows if r['engine']==engine and r['part'].startswith(family)]
            answered = sum(r['phase_outcomes']['build'].get('answered',0) for r in subset)
            count = sum(r['samples'] for r in subset)
            cells.append(f'{answered} / {count}')
        features.append('<tr>'+''.join(cell(c) for c in cells)+'</tr>')
    host = data['host']
    timing = data['timing']
    output.write_text(f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Target-part comparison | Catallax</title><meta name="description" content="Measured construction, parameter rebuilds and correctness on thirteen fixed CAD targets across five tools.">
<link rel="icon" href="brand/mark.svg" type="image/svg+xml"><link rel="stylesheet" href="brand/tokens.css"><link rel="stylesheet" href="css/site.css"><link rel="stylesheet" href="css/comparison.css"></head>
<body><a class="skip" href="#main">Skip to content</a><header class="comparison-nav"><a href="/">Catallax</a><nav aria-label="Comparison sections"><a href="#speed">Timings</a><a href="#correctness">Correctness</a><a href="#range">Target range</a><a href="#method">Method</a></nav></header>
<main id="main" class="comparison"><h1>Same targets.<br>Five tools.</h1><p class="intro">Construction and one parameter edit on thirteen fixed P0–P5 parts. Timings and correctness are reported separately. A quick refusal is not a completed part.</p>
<p>Measured on {text(host['machine_class'])}, {text(host.get('cpu_model',host['cpu']))}. {text(host['system'])}. One core, {timing['warmups']} warmup and {timing['measured_pairs']} measured samples per tool and part. Catallax work allowance: {timing['max_work_per_catallax_phase']:,} per phase. <a href="comparison-results.json">Committed measurements and checksums</a>.</p>
<section id="speed"><h2>Build and edit time</h2><label for="part">Target part</label><select id="part"><option value="all">All parts</option>{options}</select><p>Median successful build time. Edit means a complete script rebuild for Catallax, CadQuery, smlibpy and Pali Ijo; SolidCraft uses its native parameter-change command. Startup is outside internal build timers and retained in the raw run. Construction and edit can include tool-specific validation or display work. These are application operations, not isolated kernel instruction timings.</p>{''.join(charts)}</section>
<section id="correctness"><h2>What the answers passed</h2><p>The target selector also filters this table. Property cells are discrepancies / checked sets of volume, area and three centroid coordinates. Membership cells are wrong / checked classifications of imported STEP, not native point queries in each tool. STEP cells are opened / checked / full pass. Names are passed / checked edits. Zero checked is not a pass.</p><div class="table-scroll" tabindex="0" role="region" aria-label="Per-part correctness"><table><thead><tr><th>Part</th><th>Tool</th><th>Variant</th><th>Property errors / sets</th><th>STEP membership errors / points</th><th>STEP opened / checked / pass</th><th>Names</th><th>Failed builds</th><th>Incomplete</th><th>Not expressible</th></tr></thead><tbody>{''.join(checks)}</tbody></table></div></section>
<section id="range"><h2>Target range by tool</h2><p>Completed construction and edit builds / scheduled builds across measured samples. An edit skipped after a failed construction remains incomplete. These counts cover only this workload, not each tool's complete feature set. A completed build can still fail the correctness checks above.</p><div class="table-scroll" tabindex="0" role="region" aria-label="Measured feature range"><table><thead><tr><th>Tool</th><th>P0 curved booleans</th><th>P1 bracket</th><th>P2 chamfers</th><th>P3 fillets</th><th>P4 shell</th><th>P5 housing</th></tr></thead><tbody>{''.join(features)}</tbody></table></div></section>
<section id="method"><h2>Versions and method</h2><div class="table-scroll" tabindex="0" role="region" aria-label="Tool versions"><table class="versions"><thead><tr><th>Tool</th><th>Version or revision</th><th>Kernel</th><th>Artifact SHA-256</th></tr></thead><tbody>{''.join(versions)}</tbody></table></div>
<p>CadQuery was selected for its released headless Python API, solid primitives, blends, shells and STEP export. smlibpy is a public MIT Python tool bundling SMLib under Apache-2.0. Pali Ijo is MIT; upstream publishes no release archive or pub.dev package, so its official 0.1.2 source snapshot is pinned by commit and checksum and compiled with Dart {text(data['dart_version'])}. No upstream code or assets are copied into this site or repository.</p>
<p>OpenCASCADE is both CadQuery's modeling kernel and the STEP checker. The checker uses cadquery-ocp {text(data['checker_version'])}; CadQuery uses 7.8.1.1.post1. Properties and membership are checked against independently authored analytic definitions, not another kernel's constructed part. Four P0 curved-intersection families have analytic membership but no implemented property formula here: block-sphere, crossing-pipes, general-cone and general-torus. Their property cells remain unmeasured. P5's junction blend meets the countersink openings, so a complete trimmed-blend oracle is unavailable here; its properties and membership remain unmeasured rather than checking a full-torus approximation. The other eight families use elementary property formulas and converged quadrature of smooth analytic fillet sections. Numerical diagnostics use the unchanged absolute 1e-6 or relative 1e-8 allowance, not a modeling tolerance or exactness certificate.</p>
<p>Membership uses OpenCASCADE to classify imported STEP on a fixed 216-point grid, with expected values from analytic solid predicates. It measures that export/import/classification path, not native point queries in the producing tool. Boundary points remain unresolved. Passing samples do not prove set equality. STEP is opened without repair and checked for topology validity and the producing tool's face count. Names check string survival, not semantic correspondence. CadQuery and smlibpy expose no persistent per-face/edge rebuild-name contract here; Pali Ijo's ancestry IDs are process-local, so those name checks are not run. Missing operations are not approximated and failed prefixes are not exported as answers.</p>
<p>Tool order: {text(timing['order'])}. One-minute load averages ranged from {data['load_average_range'][0]:.2f} to {data['load_average_range'][1]:.2f}; the machine label alone does not prove an idle CPU. Process timeouts are incomplete, never geometric refusals. No ratios or general speed ranking are asserted. <a href="https://github.com/BohnBawerick/catallax/tree/main/.claude/skills/verify-catallax/scripts/head-to-head">Scripts, release pins and rerun commands</a>.</p></section>
<footer><a href="/">Back to Catallax</a><p>Generated from <a href="comparison-results.json">comparison-results.json</a>. Raw run SHA-256: <code>{text(data['raw_sha256'])}</code>.</p></footer></main>
<script>const select=document.getElementById('part');function filter(){{for(const item of document.querySelectorAll('[data-part]'))item.hidden=select.value!=='all'&&item.dataset.part!==select.value;}}select.value='P1-bracket.json';select.addEventListener('change',filter);filter();</script></body></html>\n''')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--import-run', type=Path)
    parser.add_argument('--results', type=Path, default=ROOT/'site/comparison-results.json')
    parser.add_argument('--output', type=Path, default=ROOT/'site/comparison.html')
    parser.add_argument('--index', type=Path, default=ROOT/'site/index.html')
    parser.add_argument('--check', action='store_true', help='fail if the generated page or homepage summary differs from results')
    args = parser.parse_args()
    if args.check and args.import_run:
        parser.error('--check cannot import a run')
    if args.import_run:
        import_run(args.import_run, args.results)
    if args.check:
        with tempfile.TemporaryDirectory() as directory:
            generated = Path(directory)/'comparison.html'
            render(args.results, generated)
            if generated.read_bytes() != args.output.read_bytes():
                raise SystemExit('comparison.html is stale; run site/render/comparison.py')
        render_index(args.results, args.index, check=True)
        print('comparison.html and index summary match committed results')
    else:
        render(args.results, args.output)
        render_index(args.results, args.index)
