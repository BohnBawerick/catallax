#!/usr/bin/env python3
"""Static census behind docs/R10-BODY-CHECKER-AUDIT.md. Read-only.

Counts, in non-test source, the code a finished B-rep body's validity rests
on: lines per crate, the topology check and the fact book, every site that
forms a recorded fact (with its origin), every site that reads one, and
every site that sets a face or coedge sense. Test files (`tests/`,
`*tests.rs`), benches, examples and a file's trailing `#[cfg(test)] mod`
are skipped, as in tools/tolerance-inventory/not_exact_sites.py.

A code line is a line that is not blank and not only a comment. The scan is
textual: a site is a line, not a parsed call.
"""
import argparse
import collections
import json
from pathlib import Path
import re

SKIPPED = ('/tests/', '/benches/', '/examples/')
FUNCTION = re.compile(r'^\s*(pub(\([^)]*\))?\s+)?(const\s+)?(async\s+)?fn\s+(\w+)')
ORIGIN = re.compile(r'FactOrigin::(Construction|Proof)\(\s*([^)]*)\)')
READS = re.compile(r'\.holds\(|\bfacts_about\(|\bhas_conic_incidence\(|\.facts\(\)\s*\.iter\(')
FACES = re.compile(r'\badd_(named_)?face\(')
COEDGES = re.compile(r'\badd_coedge\(')
TOPOLOGY_CHECK = re.compile(r'\bmodel\s*\.\s*check\(|\.model\.check\(')

# The modules named in the audit's trusted-base tables.
CHECK = 'crates/catallax-topology/src/check.rs'
FACT_BOOK = 'crates/catallax-numbers/src/facts.rs'
STORE = 'crates/catallax-shapes/src/store.rs'
CHECKED_VIEWS = [
    'crates/catallax-topology/src/spline_prism.rs',
    'crates/catallax-topology/src/oblique_draft.rs',
    'crates/catallax-topology/src/rational_patch.rs',
    'crates/catallax-topology/src/bezier_extrusion.rs',
    'crates/catallax-topology/src/linear_fillet.rs',
    'crates/catallax-topology/src/rolling_mitre.rs',
    'crates/catallax-topology/src/rational_atlas.rs',
]


def test_modules(root):
    """Files that `#[cfg(test)] mod name;` declares, and their subdirectories."""
    found = set()
    for path in (root / 'crates').rglob('*.rs'):
        lines = path.read_text().split('\n')
        base = path.parent if path.stem in ('lib', 'main', 'mod') else path.parent / path.stem
        for index, line in enumerate(lines[:-1]):
            declared = re.match(r'\s*(pub(\([^)]*\))?\s+)?mod (\w+);', lines[index + 1])
            if line.strip() == '#[cfg(test)]' and declared:
                name = declared.group(3)
                found.update({(base / f'{name}.rs').relative_to(root).as_posix(),
                              (base / name).relative_to(root).as_posix() + '/'})
    return found


def sources(root):
    tests = test_modules(root)
    for path in sorted((root / 'crates').rglob('*.rs')):
        rel = path.relative_to(root).as_posix()
        if (any(part in rel for part in SKIPPED) or rel.endswith('tests.rs')
                or rel in tests or any(rel.startswith(t) for t in tests if t.endswith('/'))):
            continue
        lines = path.read_text().split('\n')
        for index, line in enumerate(lines):
            if (line.strip() == '#[cfg(test)]' and index + 1 < len(lines)
                    and re.match(r'\s*mod \w+\s*\{', lines[index + 1])):
                lines = lines[:index]
                break
        yield rel, lines


def code_lines(lines):
    count, block = 0, False
    for line in lines:
        text = line.strip()
        if block:
            if '*/' in text:
                block = False
                text = text.split('*/', 1)[1].strip()
            else:
                continue
        if text.startswith('/*'):
            block = '*/' not in text
            continue
        if text and not text.startswith('//'):
            count += 1
    return count


def crate_of(rel):
    return rel.split('/')[1]


def census(root):
    files = {rel: lines for rel, lines in sources(root)}
    per_crate = collections.Counter()
    for rel, lines in files.items():
        per_crate[crate_of(rel)] += code_lines(lines)

    writes, reads, faces, coedges, checks = [], [], [], [], []
    for rel, lines in files.items():
        function = None
        for number, line in enumerate(lines, 1):
            found = FUNCTION.match(line)
            if found:
                function = found.group(5)
            site = {'file': rel, 'line': number, 'function': function}
            for kind, value in ORIGIN.findall(line):
                origin = value.strip()
                writes.append({**site, 'origin': kind,
                               'name': origin.strip('"') if kind == 'Construction' else origin})
            if READS.search(line):
                reads.append(site)
            if FACES.search(line) and 'fn add_' not in line:
                faces.append(site)
            if COEDGES.search(line) and 'fn add_' not in line:
                coedges.append(site)
            if TOPOLOGY_CHECK.search(line):
                checks.append(site)

    def by_crate(sites):
        counted = collections.Counter(crate_of(site['file']) for site in sites)
        return dict(sorted(counted.items()))

    construction = [w for w in writes if w['origin'] == 'Construction']
    names = collections.Counter(w['name'] for w in construction)
    module = {path: code_lines(files[path]) for path in
              [CHECK, FACT_BOOK, STORE, *CHECKED_VIEWS] if path in files}
    return {
        'format': 'catallax-r10-census/1',
        'code_lines_by_crate': dict(sorted(per_crate.items())),
        'code_lines_total': sum(per_crate.values()),
        'modules': module,
        'fact_writes': {
            'sites': len(writes),
            'construction_sites': len(construction),
            'proof_sites': len(writes) - len(construction),
            'by_crate': by_crate(writes),
            'construction_names': dict(sorted(names.items())),
            'list': writes,
        },
        'fact_reads': {'sites': len(reads), 'by_crate': by_crate(reads), 'list': reads},
        'face_sense_writes': {'sites': len(faces), 'by_crate': by_crate(faces),
                              'files': len({s['file'] for s in faces}), 'list': faces},
        'coedge_sense_writes': {'sites': len(coedges), 'by_crate': by_crate(coedges),
                                'files': len({s['file'] for s in coedges})},
        'topology_check_calls': {'sites': len(checks), 'by_crate': by_crate(checks),
                                 'list': checks},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    print(json.dumps(census(args.root), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
