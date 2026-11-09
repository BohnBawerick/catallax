#!/usr/bin/env python3
"""Generate the static fillet gallery from the corpus recipes.

The gallery shows geometry only. Kernel outcomes stay off the page until the
full runs are published, so no expected or measured status is read here.
"""
import argparse
from html import escape
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
CASES = ROOT/'parts/fillets/cases.jsonl'
SITE = ROOT/'site'
REPO = 'https://github.com/BohnBawerick/catallax'

FAMILIES = {
    'box': 'Plane edges', 'variable': 'Variable radius', 'corner': 'Corners and mitres',
    'L': 'Concave', 'cylinder': 'Cylinder', 'cone': 'Cone', 'dome': 'Sphere',
    'torus': 'Torus', 'pipes': 'Crossing cylinders', 'pyramid': 'Four-face vertex',
}


def family(tags):
    tags = set(tags)
    if tags & {'crossing', 'ring'} and 'cylinder' in tags and 'plane' not in tags:
        return 'pipes'
    for name, wanted in (('variable', {'variable-radius'}), ('torus', {'torus'}), ('dome', {'sphere'}),
                         ('cone', {'cone'}), ('cylinder', {'cylinder'}), ('pyramid', {'vertex-4'}),
                         ('corner', {'vertex-3', 'mitre', 'setback'}), ('L', {'concave'})):
        if tags & wanted and not (name == 'torus' and 'cylinder' in tags):
            return name
    return 'box'


# Illustration art: one isometric wireframe per family, not kernel output.
X, Y, Z = (1, 0, 0), (0, 1, 0), (0, 0, 1)


def arc(c, u, v, r, t0, t1, n=18):
    return [tuple(c[k] + math.cos(t)*r*u[k] + math.sin(t)*r*v[k] for k in range(3))
            for t in (t0 + (t1 - t0)*i/n for i in range(n + 1))]


def front_back(points):
    """Split a closed loop into the halves facing and hidden from the viewer."""
    seen, hidden, run, facing = [], [], None, None
    for p in points:
        f = p[0] + p[1] >= -1e-9
        if run is None or f != facing:
            if run:
                run.append(p)
                (seen if facing else hidden).append(run)
            run, facing = [p], f
        else:
            run.append(p)
    if run:
        (seen if facing else hidden).append(run)
    return seen, hidden


def ring(c, r, n=48):
    return arc(c, X, Y, r, 0, 2*math.pi, n)


def box(r0=3, r1=None):
    w, d, h = 20, 20, 10
    r1 = r0 if r1 is None else r1
    edges = [[(0, 0, h), (w - r0, 0, h)], [(w - r1, d, h), (0, d, h), (0, 0, h)],
             [(w, 0, h - r0), (w, 0, 0), (w, d, 0), (w, d, h - r1)], [(0, d, h), (0, d, 0), (w, d, 0)]]
    accent = [[(w - r0, 0, h), (w - r1, d, h)], [(w, 0, h - r0), (w, d, h - r1)],
              arc((w - r0, 0, h - r0), X, Z, r0, 0, math.pi/2), arc((w - r1, d, h - r1), X, Z, r1, 0, math.pi/2)]
    hidden = [[(0, 0, 0), (w, 0, 0)], [(0, 0, 0), (0, d, 0)], [(0, 0, 0), (0, 0, h)]]
    return edges, accent, hidden


def corner():
    w, d, h, r = 18, 18, 14, 4
    c = (w - r, d - r, h - r)
    edges = [[(0, 0, h), (w - r, 0, h)], [(0, d - r, h), (0, 0, h)], [(0, d, h - r), (0, d, 0), (w - r, d, 0)],
             [(w, 0, h - r), (w, 0, 0), (w, d - r, 0)]]
    accent = [[(0, d - r, h), (w - r, d - r, h)], [(0, d, h - r), (w - r, d, h - r)],
              [(w - r, 0, h), (w - r, d - r, h)], [(w, 0, h - r), (w, d - r, h - r)],
              [(w - r, d, 0), (w - r, d, h - r)], [(w, d - r, 0), (w, d - r, h - r)],
              arc((0, d - r, h - r), Y, Z, r, 0, math.pi/2), arc((w - r, 0, h - r), X, Z, r, 0, math.pi/2),
              arc((w - r, d - r, 0), X, Y, r, 0, math.pi/2), arc(c, X, Y, r, 0, math.pi/2),
              arc(c, X, Z, r, 0, math.pi/2), arc(c, Y, Z, r, 0, math.pi/2)]
    hidden = [[(0, 0, 0), (w, 0, 0)], [(0, 0, 0), (0, d, 0)], [(0, 0, 0), (0, 0, h)]]
    return edges, accent, hidden


def ell():
    w, t, h, d, r = 22, 5, 18, 14, 4
    profile = [(0, 0), (w, 0), (w, t)]
    profile += [(p[0], p[2]) for p in arc((t + r, 0, t + r), X, Z, r, -math.pi/2, -math.pi, 10)]
    profile += [(t, h), (0, h), (0, 0)]
    front = [(a, d, b) for a, b in profile]
    back = [(a, 0, b) for a, b in profile]
    fa, ba = front[3:14], back[3:14]
    edges = [front, back[1:-1], [(w, 0, 0), (w, d, 0)], [(w, 0, t), (w, d, t)], [(t, 0, h), (t, d, h)], [(0, 0, h), (0, d, h)]]
    accent = [fa, ba, [fa[0], ba[0]], [fa[-1], ba[-1]]]
    hidden = [[(0, 0, h), (0, 0, 0), (w, 0, 0)], [(0, 0, 0), (0, d, 0)]]
    return edges, accent, hidden


def cylinder():
    R, H, r = 9, 16, 3
    edges, hidden = front_back(ring((0, 0, 0), R))
    mid_seen, mid_hidden = front_back(ring((0, 0, H - r), R))
    s = R*math.sqrt(0.5)
    k = (R - r)/R
    edges += [[(s, -s, 0), (s, -s, H - r)], [(-s, s, 0), (-s, s, H - r)]]
    accent = mid_seen + [ring((0, 0, H), R - r),
                         arc((s*k, -s*k, H - r), (math.sqrt(0.5), -math.sqrt(0.5), 0), Z, r, 0, math.pi/2, 8),
                         arc((-s*k, s*k, H - r), (-math.sqrt(0.5), math.sqrt(0.5), 0), Z, r, 0, math.pi/2, 8)]
    return edges, accent, hidden + mid_hidden


def cone():
    R, R2, H, r = 10, 5, 14, 1.6
    edges, hidden = front_back(ring((0, 0, 0), R))
    s = math.sqrt(0.5)
    top = H - r*0.6
    edges += [[(R*s, -R*s, 0), (R2*s, -R2*s, top)], [(-R*s, R*s, 0), (-R2*s, R2*s, top)]]
    return edges, [ring((0, 0, top), R2 + 0.4), ring((0, 0, H), R2 - r)], hidden


def dome():
    edges, _, hidden = box(0.001, 0.001)
    c, R = (10, 10, 10), 7
    u, v = (math.sqrt(0.5), -math.sqrt(0.5), 0), (0.4082, 0.4082, 0.8165)
    seen, behind = front_back(ring(c, R))
    return edges + [arc(c, u, v, R, 0, math.pi, 36)] + seen, [ring(c, R + 1.4)], hidden + behind


def torus():
    R, r = 9, 3
    edges, hidden = front_back(ring((0, 0, 0), R + r, 56))
    inner_seen, inner_hidden = front_back(ring((0, 0, 0), R - r, 40))
    return edges + inner_hidden, [ring((0, 0, r), R, 56)], hidden + inner_seen


def pipes():
    R, L, k = 5, 14, math.sqrt(0.5)
    edges, accent, hidden = [], [], []
    for a, u in ((X, Y), (Y, X)):
        for t in (-L, L):
            loop = arc(tuple(a[i]*t for i in range(3)), u, Z, R, 0, 2*math.pi, 40)
            (edges if t > 0 else hidden).append(loop)
        for cs in ((-k, k), (k, -k)):
            off = (u[0]*cs[0]*R, u[1]*cs[0]*R, cs[1]*R)
            edges.append([tuple(a[i]*-L + off[i] for i in range(3)), tuple(a[i]*L + off[i] for i in range(3))])
    for sign in (1, -1):
        accent.append([(R*math.cos(t), sign*R*math.cos(t), R*math.sin(t)) for t in (math.pi*i/24 for i in range(25))])
    return edges, accent, hidden


def pyramid():
    w, h, r = 20, 14, 2.6
    c = (10, 10, h)
    edges = [[(w, 0, 0), (w, w, 0), (0, w, 0)]]
    hidden = [[(0, w, 0), (0, 0, 0), (w, 0, 0)], [(0, 0, 0), (c[0] - r*0.5, c[1] - r*0.5, h - r)]]
    for p in ((w, 0, 0), (w, w, 0), (0, w, 0)):
        edges.append([p, (c[0] + (p[0] - c[0])*r/14, c[1] + (p[1] - c[1])*r/14, h - r)])
    return edges, [ring((c[0], c[1], h - r), r*0.95, 30), ring((c[0], c[1], h - r*0.35), r*0.55, 24)], hidden


SHAPES = {'box': box, 'variable': lambda: box(4.5, 0.6), 'corner': corner, 'L': ell, 'cylinder': cylinder,
          'cone': cone, 'dome': dome, 'torus': torus, 'pipes': pipes, 'pyramid': pyramid}


def art(name, w=200, h=130, pad=14):
    edges, accent, hidden = SHAPES[name]()
    project = lambda p: ((p[0] - p[1])*0.866, (p[0] + p[1])*0.5 - p[2])
    flat = [project(p) for line in edges + accent + hidden for p in line]
    x0, x1 = min(q[0] for q in flat), max(q[0] for q in flat)
    y0, y1 = min(q[1] for q in flat), max(q[1] for q in flat)
    s = min((w - 2*pad)/(x1 - x0), (h - 2*pad)/(y1 - y0))
    ox, oy = (w - (x1 - x0)*s)/2 - x0*s, (h - (y1 - y0)*s)/2 - y0*s

    def path(line):
        return ''.join(f'{"L" if i else "M"}{project(p)[0]*s + ox:.1f} {project(p)[1]*s + oy:.1f}' for i, p in enumerate(line))
    return ''.join(f'<path class="{kind}" d="{path(line)}"/>'
                   for kind, lines in (('h', hidden), ('e', edges), ('a', accent)) for line in lines if len(line) > 1)


# Recipe text: what each step builds, with the corpus values as written.
def size(values):
    return ' × '.join(values)


def describe(step):
    op = step['op']
    if op == 'block':
        return f'{size(step["size"])} mm'
    if op == 'cylinder':
        return f'r {step["radius"]}, h {step["height"]} mm'
    if op == 'cone':
        return f'r {step["base_radius"]} to {step["top_radius"]}, h {step["height"]} mm'
    if op == 'sphere':
        return f'r {step["radius"]} mm'
    if op == 'torus':
        return f'R {step["major"]}, r {step["minor"]} mm'
    if op == 'sketch':
        segments = sum(len(loop) for loop in step['loops'])
        return f'{segments}-segment closed profile'
    if op == 'extrude':
        return f'profile {step["profile"]}, {step["distance"]} mm'
    if op == 'boolean':
        return f'{step["kind"]} {step["target"]} and {step["tool"]}'
    if op == 'fillet':
        edges = len(step['edges']) if isinstance(step['edges'], list) else step['edges']
        count = f'{edges} edge{"s" if edges != 1 else ""}' if isinstance(edges, int) else 'all edges'
        return f'r {step["radius"]} mm on {count} of {step["body"]}'
    if op == 'fillet_variable':
        law = step['edges'][0]['radius']
        return f'r {law["start"]} to {law["end"]} mm, {law["kind"]} law, on {step["body"]}'
    raise ValueError(f'no gallery description for operation {op!r}')


def headline(recipe):
    """The first stock size and every blend radius, for the card."""
    steps = recipe['steps']
    stock = next(s for s in steps if s['op'] not in ('sketch', 'boolean'))
    blends = [s for s in steps if s['op'].startswith('fillet')]
    radii = []
    for blend in blends:
        if blend['op'] == 'fillet':
            radii.append(f'r {blend["radius"]}')
        else:
            law = blend['edges'][0]['radius']
            radii.append(f'r {law["start"]}–{law["end"]}')
    first = f'extruded profile, {stock["distance"]} mm' if stock['op'] == 'extrude' else f'{stock["op"]} {describe(stock)}'
    return f'{first} · fillet {", then ".join(radii)} mm'


def chrome(title, description, canonical, depth, current, body, scripts=''):
    up = '../'*depth
    pages = [('two-kernels.html', 'Two kernels'), ('gallery.html', 'Gallery'), ('compare.html', 'Compare'), ('kernel.html', 'Kernel')]
    links = ''.join(f'\n  <a class="bar__link" href="{up}{href}"{" aria-current=\"page\"" if href == current else ""}>{name}</a>'
                    for href, name in pages)
    return f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)} | Catallax</title>
<meta name="description" content="{escape(description)}">
<link rel="canonical" href="https://kernel.zigrawbohnbawerick.org/{canonical}">
<meta name="color-scheme" content="dark">
<meta name="theme-color" content="#0b0c0d">
<link rel="icon" href="{up}brand/mark.svg" type="image/svg+xml">
<script>document.documentElement.classList.add('js')</script>
<link rel="stylesheet" href="{up}brand/fonts.css">
<link rel="stylesheet" href="{up}brand/tokens.css">
<link rel="stylesheet" href="{up}bare.css">
<link rel="stylesheet" href="{up}read.css">
<link rel="stylesheet" href="{up}evidence.css">
</head>
<body>
<svg width="0" height="0" style="position:absolute" aria-hidden="true">
<symbol id="s-projection" viewBox="0 0 48 24"><path d="M3 8.5 25 3v18L3 15.5Z"/><circle cx="37" cy="12" r="9"/><circle cx="37" cy="12" r="4"/></symbol>
</svg>
<a class="skip" href="#main">Skip to content</a>
<div class="frame" aria-hidden="true"></div>
<header class="bar">
  <a class="bar__brand" href="{up}index.html" aria-label="Catallax, home"><img src="{up}brand/mark.svg" alt="" width="28" height="28"><span>catallax</span></a>{links}
  <a class="bar__link" href="{REPO}" target="_blank" rel="noopener">Code</a>
</header>
{body}
<footer class="block mono">
  <div class="block__cell block__cell--sym"><svg class="sym sym--wide"><use href="#s-projection"/></svg></div>
  <div class="block__cell"><small>Title</small><span>{escape(title)}</span></div>
  <div class="block__cell"><small>Source</small><a href="{REPO}/blob/main/parts/fillets/cases.jsonl" target="_blank" rel="noopener">Fillet corpus</a></div>
  <div class="block__cell"><small>Drawn</small><a href="https://zigrawbohnbawerick.org/" target="_blank" rel="noopener">Zigrawbohnbawerick Pty Ltd</a></div>
  <div class="block__cell"><small>Licence</small><a href="{REPO}/blob/main/LICENSE" target="_blank" rel="noopener">MIT</a></div>
  <div class="block__cell"><small>Owner</small><a href="https://zigrawbohnbawerick.org/" target="_blank" rel="noopener">© 2026 Zigrawbohnbawerick Pty Ltd</a></div>
  <div class="block__cell"><small>Notes</small><span class="block__links"><a href="{up}why.html">Why</a><a href="{up}thanks.html">Thanks</a><a href="{up}THIRD-PARTY.txt">Third-party licences</a></span></div>
</footer>
<script src="{up}read.js"></script>{scripts}
</body>
</html>
'''


def heading(lines):
    spans = '\n'.join(f'      <span class="ln"><span>{line}</span></span>' for line in lines)
    return f'''    <h1>
{spans}
    </h1>
    <span class="dim js-only" aria-hidden="true"><i></i><b></b><i></i></span>'''


def render_grid(cases):
    used = [name for name in FAMILIES if any(c['family'] == name for c in cases)]
    symbols = ''.join(f'\n<symbol id="f-{name}" viewBox="0 0 200 130">{art(name)}</symbol>' for name in used)
    options = ''.join(f'<option value="{name}">{FAMILIES[name]}</option>' for name in used)
    cards = ''.join(f'''
      <li class="card" data-family="{c["family"]}" data-text="{escape(c["id"] + " " + " ".join(c["tags"]))}">
        <a href="gallery/{c["id"]}.html">
          <svg class="wire" aria-hidden="true"><use href="#f-{c["family"]}"/></svg>
          <span class="card__id">{escape(c["id"])}</span>
          <span class="card__geo">{escape(headline(c["recipe"]))}</span>
          <span class="card__tags mono">{escape(" · ".join(c["tags"]))}</span>
        </a>
      </li>''' for c in cases)
    body = f'''<svg width="0" height="0" style="position:absolute" aria-hidden="true">{symbols}
</svg>
<main id="main" class="read ev">
  <header class="read__head">
{heading(['Every case,', '<em>browsable</em>.'])}
  </header>
  <p class="ev__lede">{len(cases)} fillet cases from the Catallax test corpus. Geometry only. Kernel results are published after the full runs.</p>
  <section aria-labelledby="cases">
    <h2 class="sr" id="cases">Cases</h2>
    <div class="filters js-only" role="search">
      <label><span class="mono">Search</span><input id="q" type="search" autocomplete="off" spellcheck="false"></label>
      <label><span class="mono">Family</span><select id="family"><option value="">All</option>{options}</select></label>
      <output class="mono" id="shown" for="q family" aria-live="polite">{len(cases)} of {len(cases)}</output>
    </div>
    <ul class="cards">{cards}
    </ul>
    <p class="illus-note mono">Thumbnails: illustration by family</p>
  </section>
</main>'''
    script = '''
<script>
(() => {
  const q = document.getElementById('q'), family = document.getElementById('family'), shown = document.getElementById('shown');
  const cards = [...document.querySelectorAll('.card')];
  const filter = () => {
    const text = q.value.trim().toLowerCase();
    let n = 0;
    for (const card of cards) {
      const on = (!family.value || card.dataset.family === family.value) && (!text || card.dataset.text.includes(text));
      card.hidden = !on;
      n += on;
    }
    shown.textContent = `${n} of ${cards.length}`;
  };
  q.addEventListener('input', filter);
  family.addEventListener('change', filter);
})();
</script>'''
    return chrome('Gallery', f'{len(cases)} fillet cases from the Catallax test corpus, browsable by geometry.',
                  'gallery.html', 0, 'gallery.html', body, script)


def render_case(cases, i):
    c = cases[i]
    steps = ''.join(f'<tr><td><code>{escape(s["id"])}</code></td><td>{escape(s["op"].replace("_", " "))}</td><td>{escape(describe(s))}</td></tr>'
                    for s in c['recipe']['steps'])
    recipe = escape(json.dumps(c['recipe'], indent=2))
    prev, nxt = cases[i - 1], cases[(i + 1) % len(cases)]
    kernels = ''.join(f'''
      <li class="kernel{" kernel--empty" if empty else ""}">
        <h3>{name}</h3>
        <div class="kernel__view">{f'<svg class="wire" viewBox="0 0 200 130" aria-hidden="true">{art(c["family"])}</svg><span class="kernel__illus mono">Illustration</span>' if not empty else ''}</div>
        <p class="kernel__state mono">{state}</p>
      </li>''' for name, state, empty in (('Catallax', 'Result pending', False), ('OpenCASCADE', 'Result pending', False),
                                           ('NVIDIA SMLib', 'Placeholder · not run yet', True)))
    body = f'''<main id="main" class="read ev">
  <header class="read__head">
    <p class="mono ev__crumb"><a href="../gallery.html">Gallery</a> / {FAMILIES[c["family"]]}</p>
{heading([escape(c["id"])])}
  </header>
  <p class="ev__lede">{escape(c["target"][0].upper() + c["target"][1:])}.</p>
  <section aria-labelledby="side">
    <h2 class="ev__h" id="side">Three kernels, one recipe</h2>
    <ul class="kernels">{kernels}
    </ul>
  </section>
  <section aria-labelledby="recipe">
    <h2 class="ev__h" id="recipe">Recipe</h2>
    <div class="table-wrap"><table class="ev__table">
      <thead><tr><th scope="col">Step</th><th scope="col">Operation</th><th scope="col">Geometry</th></tr></thead>
      <tbody>{steps}</tbody>
    </table></div>
    <p class="mono ev__tags">{escape(" · ".join(c["tags"]))}</p>
    <details class="code"><summary class="mono">Recipe JSON</summary><pre><code>{recipe}</code></pre></details>
  </section>
  <nav class="ev__pager mono" aria-label="Cases">
    <a href="{prev["id"]}.html" rel="prev">← {escape(prev["id"])}</a>
    <a href="{nxt["id"]}.html" rel="next">{escape(nxt["id"])} →</a>
  </nav>
</main>'''
    return chrome(c['id'], f'Fillet case {c["id"]} from the Catallax corpus: recipe and geometry.',
                  f'gallery/{c["id"]}.html', 1, None, body)


def pages(source):
    cases = [json.loads(line) for line in source.read_text().splitlines() if line.strip()]
    for c in cases:
        c['family'] = family(c['tags'])
    out = {'gallery.html': render_grid(cases)}
    for i, c in enumerate(cases):
        if not re.fullmatch(r'[a-z0-9-]+', c['id']):
            raise ValueError(f'case id {c["id"]!r} is not a safe file name')
        out[f'gallery/{c["id"]}.html'] = render_case(cases, i)
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', type=Path, default=CASES)
    parser.add_argument('--site', type=Path, default=SITE)
    parser.add_argument('--check', action='store_true', help='fail if the committed gallery differs from the corpus')
    args = parser.parse_args()
    generated = pages(args.cases)
    if args.check:
        committed = {str(p.relative_to(args.site)) for p in (args.site/'gallery').glob('*.html')} | {'gallery.html'}
        stale = sorted(name for name in generated if not (args.site/name).is_file()
                       or (args.site/name).read_text() != generated[name])
        extra = sorted(committed - set(generated))
        if stale or extra:
            raise SystemExit(f'gallery is stale ({", ".join(stale + extra)}); run site/render/gallery.py')
        print(f'gallery matches {len(generated) - 1} corpus cases')
    else:
        (args.site/'gallery').mkdir(exist_ok=True)
        for old in (args.site/'gallery').glob('*.html'):
            if f'gallery/{old.name}' not in generated:
                old.unlink()
        for name, text in generated.items():
            (args.site/name).write_text(text)
