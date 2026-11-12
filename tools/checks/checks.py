#!/usr/bin/env python3
"""Conservative merge selection and exact-run evidence. Never prints replay inputs."""
import argparse
from collections import Counter
from fractions import Fraction
import fnmatch
import hashlib
import importlib.util
import io
import json
import os
import re
from pathlib import Path
import subprocess
import shutil
import sys
import tarfile
import tempfile
import time

ROOT = Path(__file__).resolve().parents[2]
RULES = 'tools/checks/rules.json'
INVENTORY = 'tools/checks/earlier-answers.json'
DRIVER = '.claude/skills/verify-catallax/scripts/control-catallax.py'
RAW_POLICY = 'tools/checks/known-raw-invalid.json'
CORRECTIONS = 'tools/checks/replay-corrections.json'
UPGRADES = 'tools/checks/replay-upgrades.json'
UPGRADE_FORMAT = 'catallax-merge-upgrades/1'
UPGRADE_CLASS = 'budget-refusal-to-answer'
# A former budget refusal that now answers with an empty body. It writes no file
# for CI to audit, so its row names an emptiness check that the replay job
# recomputes exactly from the question, without any kernel Boolean.
EMPTY_CLASS = 'budget-refusal-to-empty'
EMPTINESS_METHOD = 'contained-in-subtracted-block'
BUDGET_REFUSAL = 'numbers.operation-budget'
STEP_CLASS = 'invalid-to-valid-step'
STEP_CORRECTION_FIELDS = {
    'step.added_edges', 'step.added_vertices',
    'step.max_deviation', 'step.approximated_curves'}
# A written decimal that changed only in its last stated digit because the
# enclosure now settles the text. Only the answer's own stated STEP fields may
# move, and only inside the stated bound; the replay job checks the actual bytes.
# A certified property interval (volume, area, centroid) is exact, so no length
# bound can license moving it: `properties.*` stays byte-identical.
LAST_DIGIT_CLASS = 'stable-last-digit'
LAST_DIGIT_FIELDS = frozenset(STEP_CORRECTION_FIELDS)
# Captain's 29 October 2026 option-A ruling. Data cannot grant itself authority.
RAW_STEP_RULING_DIGESTS = {
    '7755b0bd5ef52fd26b40197d73a8e9d767683b78644b24d496795dc0b4217e9e'}


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, sort_keys=True, indent=2) + '\n')


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def optional_sha(path):
    return sha(path) if Path(path).exists() else None


def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args])


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def matches(path, patterns):
    return any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns)


def cargo_graph(directory):
    # No target filter: metadata retains optional, build, dev and target-specific
    # dependencies. All features are requested; only workspace edges select jobs.
    metadata = json.loads(subprocess.check_output([
        'cargo', 'metadata', '--format-version', '1', '--no-deps', '--all-features'],
        cwd=directory))
    names = {package['name'] for package in metadata['packages']}
    all_edges, build_edges = {}, {}
    for package in metadata['packages']:
        deps = package['dependencies']
        all_edges[package['name']] = {dep['name'] for dep in deps if dep['name'] in names}
        build_edges[package['name']] = {dep['name'] for dep in deps
                                        if dep['name'] in names and dep['kind'] != 'dev'}
    return all_edges, build_edges


def closure(seeds, edges):
    result = set(seeds)
    pending = list(seeds)
    while pending:
        for target in edges.get(pending.pop(), ()):
            if target not in result:
                result.add(target)
                pending.append(target)
    return result


def replay_catalog(metadata):
    """Use only the committed public inventory, never infer ownership."""
    if (not metadata or any(type(row['index']) is not int for row in metadata)
            or [row['index'] for row in metadata] != list(range(len(metadata)))):
        raise ValueError('replay inventory must be nonempty, ordered and complete')
    families = {}
    for row in metadata:
        family = row.get('family')
        if not isinstance(family, str) or not family or not all(c.isalnum() or c == '-' for c in family):
            raise ValueError('invalid replay family')
        families.setdefault(family, []).append(row['index'])
    return [dict(family=family, partition=part, indices=families[family])
            for part, family in enumerate(sorted(families), 1)]


def replay_partition(plan, partition):
    rows = plan['selection']['replay_partitions']
    if type(partition) is not int or not 1 <= partition <= len(rows):
        raise ValueError('missing replay partition')
    row = rows[partition - 1]
    if row['partition'] != partition or not row['indices']:
        raise ValueError('empty or malformed replay partition')
    return row


def indexed_replay_cases(report):
    cases = report['cases']
    indices = report.get('replay_indices', list(range(len(cases))))
    if (not cases or len(indices) != len(cases) or len(set(indices)) != len(indices)
            or any(type(index) is not int or index < 0 for index in indices)):
        raise ValueError('empty or malformed replay case identities')
    return list(zip(indices, cases))


def selection_errors(selection, metadata=None):
    jobs = selection.get('jobs', {})
    errors = []
    if jobs.get('fmt') != 1 or any(type(n) is not int or n < 1 for n in jobs.values()):
        errors.append('missing mandatory fmt or invalid job partitions')
    if jobs.get('tests') and not jobs.get('test-build'):
        errors.append('tests require test-build')
    if any(name in jobs for name in ('replay', 'step-regressions', 'scoreboard')) and not jobs.get('build'):
        errors.append('binary consumers require build')
    if set(jobs) & {'clippy', 'tests', 'test-build'} and not selection.get('packages'):
        errors.append('selected Cargo checks require nonempty packages')
    rows = selection.get('replay_partitions', [])
    indices = [index for row in rows for index in row.get('indices', [])]
    if 'replay' in jobs:
        if (len(rows) != jobs['replay'] or not indices
                or len(indices) != len(set(indices))
                or len(indices) != selection.get('selected_cases')
                or any(row.get('partition') != part or not row.get('indices')
                       for part, row in enumerate(rows, 1))):
            errors.append('replay requires nonempty, disjoint, complete partitions')
        if metadata is not None and rows != replay_catalog(metadata):
            errors.append('unknown replay owners require the full public catalog')
    elif rows or selection.get('selected_cases') != 0:
        errors.append('unselected replay has cases')
    return errors


def changed_check_jobs(before, after):
    """Compare effective job definitions and retain their runnable prerequisites."""
    old_jobs, new_jobs = before.get('jobs', {}), after.get('jobs', {})
    inherited = any(before.get(key) != after.get(key) for key in ('env', 'defaults', 'permissions'))
    changed = {name for name, job in new_jobs.items()
               if inherited or old_jobs.get(name) != job}
    runnable = set(new_jobs) - {'select', 'check', 'timings'}
    dependencies = {}
    for name, job in new_jobs.items():
        needs = job.get('needs', [])
        dependencies[name] = set([needs] if isinstance(needs, str) else needs) & runnable
    return closure(changed & runnable, dependencies)


def select(paths, rules, graph, build_graph, event='pull_request', ref='', changed_jobs=()):
    names = set(graph)
    reverse = {name: {consumer for consumer, deps in graph.items() if name in deps}
               for name in names}
    binary = closure({'catallax-cli'}, build_graph)
    packages, reasons = set(), []
    replay = scoreboard = build = False
    full_event = event in {'workflow_dispatch', 'schedule'} or ref.startswith('refs/heads/release/')
    substantive = [path for path in paths if not matches(path, rules['documentation_paths'])]
    ci_only = (not full_event and not (set(changed_jobs) & {'tests', 'replay'}) and bool(substantive)
               and all(matches(path, rules.get('ci_paths', [])) for path in substantive))
    if ci_only:
        selected = set(changed_jobs) & rules['jobs'].keys() - {'tests', 'replay'}
        selected.add('fmt')
        packages = names if selected & {'clippy', 'test-build'} else set()
        return {'packages': sorted(packages),
                'jobs': {name: rules['jobs'][name] for name in sorted(selected)},
                'selected_cases': 0, 'ci_only': True,
                'reasons': ['CI-only change: tooling checks and changed job definitions, without replay or test shards']}
    for path in paths:
        if full_event or matches(path, rules['full_paths']):
            packages |= names
            replay = scoreboard = build = True
            reasons.append(f'{path}: check policy or shared build input')
        elif path.startswith('crates/') and len(path.split('/')) > 2 and path.split('/')[1] in names:
            owner = path.split('/')[1]
            packages |= closure({owner}, reverse)
            build = build or owner in binary
            scoreboard = scoreboard or owner in binary
            replay = replay or (owner in binary and matches(path, rules['replay_paths']))
            reasons.append(f'{path}: crate and reverse consumers, unknown corpus owners')
        elif path.startswith('parts/'):
            packages |= closure({'catallax-app', 'catallax-cli', 'catallax-compare',
                                 'catallax-features', 'catallax-field', 'catallax-iges'}, reverse)
            scoreboard = build = True
            reasons.append(f'{path}: target-part fixtures and consumers')
        elif path.startswith('.claude/skills/verify-catallax/'):
            packages |= closure({'catallax-app', 'catallax-cli', 'catallax-measure'}, reverse)
            build = True
            replay = replay or matches(path, rules['replay_paths'])
            scoreboard = True
            reasons.append(f'{path}: verification driver or shared example')
        elif matches(path, rules['documentation_paths']):
            reasons.append(f'{path}: documentation only')
        else:
            packages |= names
            replay = scoreboard = build = True
            reasons.append(f'{path}: unknown mapping, full single-configuration checks')
    if full_event:
        packages = names
        replay = scoreboard = build = True
        reasons.append('nightly, release push or manual dispatch: every check')
    jobs = {'fmt': 1}
    if packages:
        jobs.update(clippy=1, tests=rules['jobs']['tests'])
        if 'test-build' in rules['jobs']:
            jobs['test-build'] = 1
    if build:
        jobs['build'] = 1
        if 'step-regressions' in rules['jobs']:
            jobs['step-regressions'] = 1
    if replay:
        jobs['replay'] = rules['jobs']['replay']
    if scoreboard:
        jobs['scoreboard'] = 1
    # A newly declared job is mandatory until its own reviewed routing exists.
    for name, partitions in rules['jobs'].items():
        if name not in {'fmt', 'clippy', 'tests', 'test-build', 'build', 'step-regressions',
                        'replay', 'scoreboard'}:
            jobs[name] = partitions
    return {'packages': sorted(packages), 'jobs': jobs,
            'selected_cases': len(read(ROOT / INVENTORY)) if replay else 0,
            'ci_only': False, 'reasons': reasons}


def union_selection(base, candidate):
    jobs = {name: max(base['jobs'].get(name, 0), candidate['jobs'].get(name, 0))
            for name in base['jobs'].keys() | candidate['jobs'].keys()}
    ci_only = base.get('ci_only', False) and candidate.get('ci_only', False)
    # Base rules can predate the split. Preserve their checks after unioning.
    if not ci_only:
        if 'tests' in jobs:
            jobs['test-build'] = 1
        if 'build' in jobs:
            jobs['step-regressions'] = 1
    return {'packages': sorted(set(base['packages']) | set(candidate['packages'])),
            'jobs': jobs,
            'ci_only': ci_only,
            'selected_cases': max(base['selected_cases'], candidate['selected_cases']),
            'base_selected': base['selected_cases'], 'candidate_selected': candidate['selected_cases'],
            'reasons': sorted(set(base['reasons'] + candidate['reasons']))}


def metadata_errors(base, candidate, corrections=None, base_policy=None, upgrades=None):
    if base == candidate:
        return []
    failure = ['accepted-answer metadata changed without base approval']
    rows = upgrade_rows(upgrades)
    if rows is None or len(base) != len(candidate):
        return failure
    # A listed upgrade settles its own row; every other change stays on the
    # correction path, unchanged.
    remaining = []
    for before, after in zip(base, candidate):
        if before != after and before['index'] in rows:
            if upgrade_metadata_errors(before, after, rows[before['index']]):
                return failure
        else:
            remaining.append((before, after))
    return correction_metadata_errors(
        [before for before, _ in remaining], [after for _, after in remaining],
        corrections, base_policy)


def upgrade_rows(upgrades):
    if upgrades is None:
        return {}
    if upgrades.get('format') != UPGRADE_FORMAT or not isinstance(upgrades.get('upgrades'), list):
        return None
    rows = upgrades['upgrades']
    if len({row.get('index') for row in rows}) != len(rows):
        return None
    return {row['index']: row for row in rows}


def upgrade_metadata_errors(before, after, row):
    """A former budget refusal may become an answer, and nothing else may change.

    The row carries no verification of its own: the new answer must write a
    STEP file, which the replay job audits directly (upgrade_audit_errors).
    """
    failure = [f"index {before['index']}: upgrade does not match the base and candidate rows"]
    if (row.get('class') not in (UPGRADE_CLASS, EMPTY_CLASS)
            or row.get('request_sha256') != before['request_sha256']
            or row.get('before_sha256') != before['accepted_answer_sha256']
            or row.get('before_manifest_sha256') != before['accepted_manifest_sha256']
            or any(value is not None for value in before['accepted_files'].values())
            or not isinstance(row.get('after_files'), dict)
            or set(row['after_files']) != set(before['accepted_files'])):
        return failure
    if row['class'] == UPGRADE_CLASS and not valid_digest(row['after_files'].get('answer.step')):
        return failure
    if row['class'] == EMPTY_CLASS and (
            any(value is not None for value in row['after_files'].values())
            or emptiness_claim(row) is None):
        return failure
    expected = json.loads(json.dumps(before))
    expected['accepted_answer_sha256'] = row.get('after_sha256')
    expected['accepted_manifest_sha256'] = row.get('after_manifest_sha256')
    expected['accepted_files'] = dict(row['after_files'])
    if (expected != after or not valid_digest(row.get('after_sha256'))
            or not valid_digest(row.get('after_manifest_sha256'))
            or any(value is not None and not valid_digest(value)
                   for value in row['after_files'].values())):
        return failure
    return []


def upgrade_manifest_errors(manifest, base_manifest, row=None, problem=None, files=()):
    """The base verifier must refuse for budget and the candidate must answer.

    An empty-class row also needs an exactly empty answer that writes no file,
    and its emptiness check must hold when recomputed from `problem`.
    """
    refusal = base_manifest.get('refusal') or {}
    if refusal.get('code') != BUDGET_REFUSAL:
        return ['upgrade base answer is not a work-budget refusal']
    if manifest.get('status') != 'answered' or manifest.get('refusal'):
        return ['upgrade candidate does not answer']
    if row is None or row.get('class') != EMPTY_CLASS:
        return []
    volume = (manifest.get('properties') or {}).get('volume') or {}
    try:
        zero = Fraction(volume['lo']) == 0 and Fraction(volume['hi']) == 0
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        zero = False
    if not zero or manifest.get('step') is not None or files:
        return ['empty upgrade candidate is not an exactly empty answer without files']
    return emptiness_errors(row, problem)


def emptiness_claim(row):
    """The row's declared emptiness check, or None when it is missing or malformed."""
    claim = row.get('emptiness')
    if (not isinstance(claim, dict) or set(claim) != {'method', 'block'}
            or claim['method'] != EMPTINESS_METHOD or not isinstance(claim['block'], str)):
        return None
    return claim


def emptiness_errors(row, problem):
    """Recompute the row's emptiness check exactly from the question.

    The result step must subtract the named block operand from a step whose
    every piece lies in that block. A piece lies in it when it is a primitive
    sharing the block's exact rotation and contained in it, a union of such
    pieces, an intersection with one such piece, or a difference whose first
    argument is one. Closed containment suffices: what is left has no interior.
    """
    claim = emptiness_claim(row)
    failure = ['empty upgrade emptiness check does not hold for the question']
    if claim is None or not isinstance(problem, dict):
        return failure
    try:
        operands, steps = problem['operands'], {step['id']: step for step in problem['steps']}
        result = steps[problem['result']]
        block = operands[claim['block']]
        if (result['op'] != 'subtract' or result['args'][1] != claim['block']
                or block['shape'] != 'block'):
            return failure
        frame = pose(block)
        size = [Fraction(value) for value in block['size']]

        def inside(name, depth=0):
            if depth > len(steps):
                return False
            if name in operands:
                return contained(operands[name], frame, size)
            step = steps[name]
            first, second = step['args']
            if step['op'] == 'unite':
                return inside(first, depth + 1) and inside(second, depth + 1)
            if step['op'] == 'intersect':
                return inside(first, depth + 1) or inside(second, depth + 1)
            if step['op'] == 'subtract':
                return inside(first, depth + 1)
            return False
        return [] if inside(result['args'][0]) else failure
    except (KeyError, TypeError, ValueError, ZeroDivisionError, IndexError):
        return failure


def pose(operand):
    """The operand's exact rotation quaternion and translation."""
    placed = operand.get('pose') or {}
    rotation = tuple(Fraction(value) for value in placed.get('rotation', ['1', '0', '0', '0']))
    translation = [Fraction(value) for value in placed.get('translation', ['0', '0', '0'])]
    if len(rotation) != 4 or len(translation) != 3 or not any(rotation):
        raise ValueError('bad pose')
    return rotation, translation


def contained(operand, frame, size):
    """Whether a block, cylinder or sphere lies in the block at `frame` of `size`."""
    rotation, translation = pose(operand)
    if rotation != frame[0]:
        return False
    w, x, y, z = rotation
    norm = w * w + x * x + y * y + z * z
    matrix = ((w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y)),
              (2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x)),
              (2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z))
    offset = [a - b for a, b in zip(translation, frame[1])]
    # The block's own coordinates of the operand's origin: the transposed rotation.
    local = [sum(matrix[row][column] * offset[row] for row in range(3)) / norm
             for column in range(3)]
    if operand['shape'] == 'block':
        low, high = local, [a + Fraction(b) for a, b in zip(local, operand['size'])]
    elif operand['shape'] == 'sphere':
        radius = Fraction(operand['radius'])
        low, high = [a - radius for a in local], [a + radius for a in local]
    elif operand['shape'] == 'cylinder':
        radius, height = Fraction(operand['radius']), Fraction(operand['height'])
        low = [local[0] - radius, local[1] - radius, local[2]]
        high = [local[0] + radius, local[1] + radius, local[2] + height]
    else:
        return False
    return all(0 <= a and b <= c for a, b, c in zip(low, high, size))


def upgrade_audit_errors(report, upgrades, audit):
    """Independently audit each newly upgraded answer in this replay.

    `report` is the replay's run.json; its `upgrade_indices` are the listed
    upgrades whose base verifier refused for budget. Each must have its own
    replayed STEP file, byte-identical to the listed after STEP, and `audit`
    (the direct OpenCASCADE property audit) must accept it with that case's
    manifest. Absent, stale, mismatched or failed evidence is an error.
    """
    errors = []
    cases = dict(indexed_replay_cases(report))
    for index in report.get('upgrade_indices', []):
        row, case = upgrades.get(index), cases.get(index)
        if row is None or case is None:
            errors.append(f'index {index}: upgrade audit has no listed row or replayed case')
            continue
        if row.get('class') == EMPTY_CLASS:
            # No file to audit; the replay job recomputed its emptiness check.
            continue
        path = Path(case['path']) / 'answer.step'
        if optional_sha(path) != row['after_files'].get('answer.step') or not path.exists():
            errors.append(f'index {index}: upgraded STEP is absent or not the listed answer')
            continue
        try:
            passed = audit(path, case['response']['manifest'])
        except Exception:
            passed = False
        if not passed:
            errors.append(f'index {index}: upgraded STEP failed the direct property audit')
    return errors


def correction_metadata_errors(base, candidate, corrections, base_policy):
    if base == candidate:
        return []
    failure = ['accepted-answer metadata changed without base approval']
    if not corrections or corrections.get('format') != 'catallax-replay-corrections/1':
        return failure
    if not base_policy or raw_step_errors([], base_policy):
        return failure
    rows = corrections.get('corrections', [])
    if len({row['index'] for row in rows}) != len(rows) or len(base) != len(candidate):
        return failure
    approved = {row['index']: row['step_sha256']
                for ruling in base_policy['rulings'] for row in ruling['approved']}
    by_index = {row['index']: row for row in rows}
    for before, after in zip(base, candidate):
        if before == after:
            continue
        index = before['index']
        correction = by_index.get(index)
        if not correction:
            return failure
        common = (
            correction.get('request_sha256') == before['request_sha256']
            and correction.get('before_sha256') == before['accepted_answer_sha256']
            and correction.get('before_manifest_sha256') == before['accepted_manifest_sha256']
            and correction.get('before_step_sha256') == before['accepted_files']['answer.step']
            and valid_digest(correction.get('unchanged_manifest_sha256'))
            and valid_digest(correction.get('after_sha256'))
            and valid_digest(correction.get('after_manifest_sha256'))
            and valid_digest(correction.get('after_step_sha256')))
        fields = correction.get('manifest_fields')
        if correction.get('class') == LAST_DIGIT_CLASS:
            if (not common or not isinstance(fields, list) or len(fields) != len(set(fields))
                    or any(not isinstance(field, str) or field not in LAST_DIGIT_FIELDS
                           for field in fields)):
                return failure
        elif correction.get('class') == STEP_CLASS:
            if (not common or index not in base_policy['known_raw_invalid']
                    or correction.get('before_step_sha256') != approved.get(index)
                    or not set(fields or []) <= STEP_CORRECTION_FIELDS):
                return failure
        else:
            return failure
        expected = json.loads(json.dumps(before))
        expected['accepted_answer_sha256'] = correction.get('after_sha256')
        expected['accepted_manifest_sha256'] = correction.get('after_manifest_sha256')
        expected['accepted_files']['answer.step'] = correction.get('after_step_sha256')
        if expected != after:
            return failure
    return []


def valid_digest(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def field_parent(value, path):
    node = value
    for part in path.split('.')[:-1]:
        if not isinstance(node, dict) or part not in node:
            raise KeyError(path)
        node = node[part]
    return node, path.rsplit('.', 1)[1]


def protected_manifest_digest(manifest, driver, fields=None):
    comparable = json.loads(json.dumps(driver.comparable(manifest)))
    for field in (STEP_CORRECTION_FIELDS if fields is None else fields):
        parent, leaf = field_parent(comparable, field)
        if isinstance(parent, dict):
            parent.pop(leaf, None)
    return driver.digest(comparable)


def correction_manifest_errors(manifest, correction, driver, base_manifest):
    if driver.digest(driver.comparable(base_manifest)) != correction['before_manifest_sha256']:
        return ['base manifest does not match the accepted before manifest']
    fields = correction.get('manifest_fields', [])
    if correction.get('class') == STEP_CLASS:
        fields = STEP_CORRECTION_FIELDS
    try:
        protected = protected_manifest_digest(base_manifest, driver, fields)
    except KeyError:
        return ['correction names a manifest field that does not exist']
    if protected != correction['unchanged_manifest_sha256']:
        return ['correction commitment does not match the authenticated base manifest']
    try:
        if protected_manifest_digest(manifest, driver, fields) != protected:
            return ['correction changes a manifest field outside its listed fields']
    except KeyError:
        return ['correction names a manifest field that does not exist']
    return []


DECIMAL_TOKEN = re.compile(rb'(?<!#)[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?')


def decimal_value(text):
    try:
        return Fraction(text)
    except (ValueError, ZeroDivisionError):
        return None


def stated_tolerance(manifest):
    value = (manifest.get('step') or {}).get('max_deviation')
    bound = decimal_value(value) if isinstance(value, str) else None
    # A deviation bound is a distance: a negative or unparseable value states
    # no accuracy at all, so it cannot license any change.
    return bound if bound is not None and bound >= 0 else None


def within_tolerance(before, after, tolerance, field=''):
    """Two decimal texts agree when their values differ inside the stated bound."""
    first, second = decimal_value(before), decimal_value(after)
    if first is None or second is None:
        return False
    if field == 'step.max_deviation':
        # A corrected enclosure may tighten the file's own bound, never widen
        # it, and a tightened bound is still never negative.
        return second >= 0 and second <= first
    return abs(second - first) <= tolerance


def last_digit_leaf_errors(before, after, tolerance, path=''):
    if before == after:
        return []
    if type(before) is not type(after) or not isinstance(before, (dict, list)):
        if (isinstance(before, str) and isinstance(after, str)
                and path in LAST_DIGIT_FIELDS
                and within_tolerance(before, after, tolerance, path)):
            return []
        return [f'{path}: changed outside its stated accuracy']
    errors = []
    if isinstance(before, dict):
        for key in sorted(set(before) | set(after)):
            child = f'{path}.{key}' if path else key
            if key not in before or key not in after:
                errors.append(f'{child}: field appeared or disappeared')
            else:
                errors += last_digit_leaf_errors(before[key], after[key], tolerance, child)
    else:
        if len(before) != len(after):
            errors.append(f'{path}: list length changed')
        else:
            for index, (a, b) in enumerate(zip(before, after)):
                errors += last_digit_leaf_errors(a, b, tolerance, f'{path}[{index}]')
    return errors


def last_digit_errors(before_dir, after_dir, driver):
    """The saved answer may move only inside its own stated accuracy."""
    errors = []
    before_manifest = read(Path(before_dir) / 'manifest.json')
    after_manifest = read(Path(after_dir) / 'manifest.json')
    tolerance = stated_tolerance(after_manifest)
    if tolerance is None:
        return ['corrected answer states no valid nonnegative accuracy bound']
    errors += last_digit_leaf_errors(driver.comparable(before_manifest),
                                     driver.comparable(after_manifest), tolerance)
    before_step, after_step = Path(before_dir) / 'answer.step', Path(after_dir) / 'answer.step'
    if before_step.exists() != after_step.exists():
        return errors + ['STEP file appeared or disappeared']
    if before_step.exists():
        before_bytes, after_bytes = before_step.read_bytes(), after_step.read_bytes()
        if before_bytes != after_bytes:
            before_tokens = DECIMAL_TOKEN.findall(before_bytes)
            after_tokens = DECIMAL_TOKEN.findall(after_bytes)
            if (DECIMAL_TOKEN.sub(b'#', before_bytes) != DECIMAL_TOKEN.sub(b'#', after_bytes)
                    or len(before_tokens) != len(after_tokens)):
                errors.append('STEP text changed outside its decimal tokens')
            else:
                for token, (a, b) in enumerate(zip(before_tokens, after_tokens)):
                    if a != b and not within_tolerance(
                            a.decode(), b.decode(), tolerance):
                        errors.append(f'STEP token {token}: changed outside the stated accuracy')
                        break
    if not errors and driver.comparable(before_manifest) == driver.comparable(after_manifest) and (
            not before_step.exists() or before_step.read_bytes() == after_step.read_bytes()):
        errors.append('saved answer changed no byte')
    return errors


def raw_step_errors(observations, policy):
    if policy.get('format') != 'catallax-known-raw-invalid/1':
        return ['unsupported known-raw-invalid policy']
    rulings = policy.get('rulings')
    ruling_digest = hashlib.sha256(json.dumps(
        rulings, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()
    if ruling_digest not in RAW_STEP_RULING_DIGESTS:
        return ['known-raw-invalid expansion requires an explicit trusted captain ruling']
    approved = {row['index']: row['step_sha256']
                for ruling in rulings for row in ruling['approved']}
    known = policy.get('known_raw_invalid')
    if (not isinstance(known, list) or any(type(index) is not int for index in known)
            or len(known) != len(set(known))):
        return ['invalid or duplicate known-raw-invalid indices']
    errors = [f'index {index}: no captain ruling for new raw-invalid entry'
              for index in set(known) - approved.keys()]
    for row in observations:
        index, fingerprint = row['index'], row['step_sha256']
        if fingerprint is None:
            if index in known:
                errors.append(f'index {index}: listed raw-invalid STEP is missing')
        elif row['valid']:
            if index in known:
                errors.append(f'index {index}: raw-valid upgrade must leave the known-raw-invalid list')
        elif index not in known:
            errors.append(f'index {index}: unlisted raw-invalid STEP')
        elif fingerprint != approved.get(index):
            errors.append(f'index {index}: STEP bytes changed while still raw-invalid')
    return errors


def decide(plan, reported, records, binary_sha256=None, actual_jobs=()):
    errors = selection_errors(plan['selection'])
    if errors:
        return errors
    required, declared = plan['selection']['jobs'], plan['declared']
    allowed_names = set(declared) | {'select', 'check'}
    # GitHub leaves the matrix expression unexpanded when replay is skipped.
    if ('replay' in declared and 'replay' not in required
            and reported.get('replay', {}).get('result') == 'skipped'):
        allowed_names.add('replay (${{ matrix.partition }})')
    for name, count in required.items():
        if count > 1:
            allowed_names.discard(name)
            allowed_names.update(f'{name} ({part})' for part in range(1, count + 1))
    for actual in actual_jobs:
        if actual.split(' / ', 1)[0] not in allowed_names:
            errors.append(f'{actual}: unexpected job')
    for name in reported.keys() - set(declared):
        errors.append(f'{name}: unexpected job')
    for name in set(declared) | set(required):
        result = reported.get(name, {}).get('result')
        if name not in required:
            if result not in (None, 'skipped'):
                errors.append(f'{name}: unexpected selected job ({result})')
            continue
        if result != 'success':
            errors.append(f'{name}: {result or "missing"}')
    seen = set()
    for record in records:
        name, part = record.get('job'), record.get('partition')
        key = (name, part)
        if name not in required:
            errors.append(f'{name}: unexpected evidence')
            continue
        if type(part) is not int or not 1 <= part <= required[name]:
            errors.append(f'{name}: unexpected partition {part}')
            continue
        if key in seen:
            errors.append(f'{name}: duplicate partition {part}')
        seen.add(key)
        expected = dict(plan['identity'])
        expected['result'] = 'success'
        if name in {'build', 'step-regressions', 'replay', 'scoreboard'}:
            expected['binary_sha256'] = binary_sha256
        else:
            expected['binary_sha256'] = None
        if name == 'replay':
            selected = replay_partition(plan, part)
            expected.update(workers=4, fresh=True, cache_hits=0,
                            family=selected['family'], indices=selected['indices'],
                            selected_cases=len(selected['indices']),
                            verified_corrections=[row['index'] for row in plan['correction_base_answers']
                                                  if row['index'] in selected['indices']])
        for field, value in expected.items():
            if type(record.get(field)) is not type(value) or record.get(field) != value:
                errors.append(f'{name}: stale {field}: evidence={record.get(field)!r}, candidate={value!r}')
        if name in {'build', 'step-regressions', 'replay', 'scoreboard'} and not binary_sha256:
            errors.append(f'{name}: missing binary_sha256 from current build')
    for name, count in required.items():
        if name in {'release-evidence', 'release-dry-run'}:
            continue  # Reusable release workflows retain their own artifact gates.
        for part in range(1, count + 1):
            if (name, part) not in seen:
                errors.append(f'{name}: partition {part}/{count} missing')
    return errors


def prepare(args):
    base_rules = json.loads(git('show', f'{args.base}:{RULES}')) if git_has(args.base, RULES) else read(ROOT / RULES)
    base_metadata = json.loads(git('show', f'{args.base}:{INVENTORY}')) if git_has(args.base, INVENTORY) else read(ROOT / INVENTORY)
    candidate_rules, candidate_metadata = read(ROOT / RULES), read(ROOT / INVENTORY)
    base_policy = (json.loads(git('show', f'{args.base}:{RAW_POLICY}'))
                   if git_has(args.base, RAW_POLICY) else read(ROOT / RAW_POLICY))
    errors = metadata_errors(base_metadata, candidate_metadata, read(ROOT / CORRECTIONS), base_policy,
                             read(ROOT / UPGRADES))
    errors += raw_step_errors([], read(ROOT / RAW_POLICY))
    errors += raw_step_errors([], base_policy)
    if errors:
        raise ValueError('; '.join(errors))
    if base_rules['workers'] != 4 or candidate_rules['workers'] != 4:
        raise ValueError('fixed workers must be 4 under both rules')
    paths = [p.decode() for p in git('diff', '--name-only', '--no-renames', '-z', args.base, 'HEAD').split(b'\0') if p]
    with tempfile.TemporaryDirectory(prefix='catallax-base-metadata-') as temporary:
        with tarfile.open(fileobj=io.BytesIO(git('archive', args.base))) as archive:
            archive.extractall(temporary, filter='data')
        base_graph, base_build = cargo_graph(temporary)
    graph, build = cargo_graph(ROOT)
    import yaml
    workflow = '.github/workflows/check.yml'
    old_workflow = yaml.safe_load(git('show', f'{args.base}:{workflow}')) or {}
    new_workflow = yaml.safe_load((ROOT / workflow).read_text()) or {}
    changed_jobs = changed_check_jobs(old_workflow, new_workflow)
    before = select(paths, base_rules, base_graph, base_build, args.event, args.ref, changed_jobs)
    after = select(paths, candidate_rules, graph, build, args.event, args.ref, changed_jobs)
    selection = union_selection(before, after)
    # A removed crate contributes its old consumers, not an impossible -p flag.
    selection['packages'] = sorted(set(selection['packages']) & set(graph))
    catalog = replay_catalog(candidate_metadata)
    selection['replay_partitions'] = catalog if 'replay' in selection['jobs'] else []
    if 'replay' in selection['jobs']:
        # Historical policy had a single full replay. It still contributes every
        # case; only execution is partitioned, with no ownership assumed.
        selection['jobs']['replay'] = len(catalog)
        selection['selected_cases'] = len(candidate_metadata)
    errors = selection_errors(selection, candidate_metadata)
    if errors:
        raise ValueError('; '.join(errors))
    if args.release_evidence:
        selection['jobs']['release-evidence'] = 1
    if args.release_dry_run:
        selection['jobs']['release-dry-run'] = 1
    plan = {'format': 'catallax-check-plan/1', 'base': args.base, 'head': args.head,
            'correction_base_answers': [before for before, after in zip(base_metadata, candidate_metadata)
                                        if before != after],
            'declared': sorted(set(base_rules['jobs']) | set(candidate_rules['jobs']) |
                               {'select', 'release-evidence', 'release-dry-run', 'timings'}),
            'identity': {'commit': git('rev-parse', 'HEAD').decode().strip(),
                         'head': args.head,
                         'run_id': str(args.run_id), 'run_attempt': str(args.run_attempt),
                         'corpus_sha256': base_rules['corpus_sha256'],
                         'raw_policy_sha256': sha(ROOT / RAW_POLICY),
                         'corrections_sha256': sha(ROOT / CORRECTIONS),
                         'upgrades_sha256': sha(ROOT / UPGRADES),
                         'metadata_sha256': sha(ROOT / INVENTORY), 'driver_sha256': sha(ROOT / DRIVER)},
            'selection': selection}
    plan['decisions'] = {
        name: {'selected': name in selection['jobs'],
               'partitions': selection['jobs'].get(name, 0),
               'reasons': ['always: formatting and check policy'] if name == 'fmt'
               else selection['reasons'] if name in selection['jobs']
               else ['no affected owner or consumer under either rule set']}
        for name in plan['declared'] if name != 'select'}
    plan['decisions']['timings'] = {
        'selected': True, 'partitions': 1, 'required_by_gate': False,
        'reasons': ['observational timing artifact after the final gate']}
    write(args.out, plan)
    if args.github_output:
        with Path(args.github_output).open('a') as out:
            for name in candidate_rules['jobs']:
                out.write(f'{name}={str(name in selection["jobs"]).lower()}\n')
            out.write('packages=' + ' '.join('-p ' + name for name in selection['packages']) + '\n')
            out.write('base_corrections=' + str(bool(plan['correction_base_answers'])).lower() + '\n')
            out.write('base_commit=' + args.base + '\n')
            matrix = [{key: row[key] for key in ('family', 'partition')}
                      for row in selection['replay_partitions']]
            out.write('replay_matrix=' + json.dumps(matrix, separators=(',', ':')) + '\n')
    print(json.dumps(selection, sort_keys=True))


def git_has(commit, path):
    return subprocess.run(['git', '-C', str(ROOT), 'cat-file', '-e', f'{commit}:{path}'],
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0


def evidence(args):
    plan = read(args.plan)
    record = {**plan['identity'], 'job': args.job, 'partition': args.partition,
              'binary_sha256': sha(args.binary) if args.binary else None, 'result': args.result}
    if args.job == 'replay':
        summary = read(args.summary)
        record.update(summary)
    write(args.out, record)


def gate(args):
    plan = read(args.plan)
    errors = selection_errors(plan['selection'], read(ROOT / INVENTORY))
    if errors:
        raise ValueError('; '.join(errors))
    current = {'commit': git('rev-parse', 'HEAD').decode().strip(),
               'run_id': os.environ.get('GITHUB_RUN_ID', 'local'),
               'run_attempt': os.environ.get('GITHUB_RUN_ATTEMPT', '1'),
               'raw_policy_sha256': sha(ROOT / RAW_POLICY),
               'corrections_sha256': sha(ROOT / CORRECTIONS),
               'upgrades_sha256': sha(ROOT / UPGRADES),
               'metadata_sha256': sha(ROOT / INVENTORY), 'driver_sha256': sha(ROOT / DRIVER)}
    for field, value in current.items():
        if plan['identity'].get(field) != value:
            raise ValueError(f'select: stale {field}')
    # The selector is mandatory even if a candidate deletes its own dependency.
    reported = json.loads(os.environ['RESULTS'])
    if reported.get('select', {}).get('result') != 'success':
        raise ValueError('select: missing or failed')
    reported.pop('select', None)
    records = [read(path) for path in Path(args.evidence).glob('*.json')]
    builds = [record for record in records if record.get('job') == 'build']
    binary_hash = builds[0].get('binary_sha256') if len(builds) == 1 else None
    errors = decide(plan, reported, records, binary_hash, Path(args.jobs).read_text().splitlines())
    if errors:
        raise ValueError('\n'.join(errors))
    print('Every selected check and partition passed with exact-candidate evidence')


def build_base(args):
    plan = read(args.plan)
    if not plan['correction_base_answers']:
        raise ValueError('no new corrections or upgrades require a base verifier')
    source = Path(args.source).resolve()
    head = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD']).decode().strip()
    if head != plan['base']:
        raise ValueError('base verifier checkout does not match the selected base')
    if subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain']).strip():
        raise ValueError('base verifier checkout is not clean')
    target = ROOT / 'target' / 'check-base'
    subprocess.run([
        'cargo', 'build', '--release', '--locked', '-p', 'catallax-cli',
        '--bin', 'catallax-verify', '--manifest-path', str(source / 'Cargo.toml'),
        '--target-dir', str(target)], cwd=source, check=True)
    destination = Path(args.out)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(target / 'release' / 'catallax-verify', destination)
    print(json.dumps({'base_commit': plan['base'], 'binary_sha256': sha(destination)}))


def run_fresh(driver, binary, requests, out, workers):
    run_args = argparse.Namespace(binary=binary, inputs=[str(requests)], source='auto',
                                 adapter=None, out=str(out), cache=None, cores=None,
                                 max_work=100_000_000, workers=workers, fresh=True, dry_run=False)
    try:
        driver.run(run_args)
    except Exception:
        raise ValueError('replay execution failed; private diagnostics retained on runner') from None
    return read(Path(out) / 'run.json')


def correction_bases(args, plan, driver, release, inputs, indices):
    accepted = [row for row in plan['correction_base_answers'] if row['index'] in indices]
    if not accepted:
        return {}, {}
    if not args.base_binary:
        raise ValueError('new corrections or upgrades require the base verifier')
    source = Path(args.out) / 'base-requests.jsonl'
    source.parent.mkdir(parents=True, exist_ok=True)
    with source.open('w') as stream:
        for row in accepted:
            stream.write(driver.canonical({'max_work': 100_000_000, **inputs[row['index']],
                                           'format': driver.CONTRACT}) + '\n')
    report = run_fresh(driver, args.base_binary, source, Path(args.out) / 'base', 1)
    if (len(report['cases']) != len(accepted)
            or report['execution']['worker_exits'] != {'0': 0}
            or report['config']['binary_sha256'] != sha(args.base_binary)):
        raise ValueError('incomplete base correction replay')
    manifests, directories = {}, {}
    for row, case in zip(accepted, report['cases']):
        directory = Path(case['path'])
        request = read(directory / 'request.json')
        request.pop('out', None)
        manifest = case['response']['manifest']
        if (driver.digest(request) != row['request_sha256']
                or driver.digest(driver.comparable(manifest)) != row['accepted_manifest_sha256']
                or release.answer_digest(driver, directory) != row['accepted_answer_sha256']
                or optional_sha(directory / 'answer.step') != row['accepted_files']['answer.step']
                or case['cache_hit'] or case['response']['work']['used'] > 100_000_000):
            raise ValueError(f"index {row['index']}: base replay does not match the accepted answer")
        manifests[row['index']] = manifest
        directories[row['index']] = directory
    return manifests, directories


def replay(args):
    started = time.monotonic()
    rules, metadata = read(ROOT / RULES), read(ROOT / INVENTORY)
    corrections = {row['index']: row for row in read(ROOT / CORRECTIONS)['corrections']}
    upgrades = upgrade_rows(read(ROOT / UPGRADES))
    if upgrades is None:
        raise ValueError('malformed replay upgrade list')
    plan = read(args.plan)
    errors = selection_errors(plan['selection'], metadata)
    if errors:
        raise ValueError('; '.join(errors))
    selected = replay_partition(plan, args.partition)
    indices = selected['indices']
    if (plan['identity']['metadata_sha256'] != sha(ROOT / INVENTORY)
            or plan['identity']['corrections_sha256'] != sha(ROOT / CORRECTIONS)
            or plan['identity'].get('upgrades_sha256') != sha(ROOT / UPGRADES)
            or plan['identity']['commit'] != git('rev-parse', 'HEAD').decode().strip()):
        raise ValueError('stale correction replay plan')
    if sha(args.inputs) != rules['corpus_sha256']:
        raise ValueError('corpus_sha256 mismatch')
    if len(os.sched_getaffinity(0)) < 4:
        raise ValueError('fixed replay requires at least four available cores')
    driver = load_module('catallax_merge_driver', ROOT / DRIVER)
    release = load_module('catallax_merge_release', ROOT / 'tools/release/release.py')
    inputs = read(args.inputs)
    if len(inputs) != len(metadata):
        raise ValueError('replay inputs must match the complete public inventory')
    base_manifests, base_directories = correction_bases(args, plan, driver, release, inputs, indices)
    upgrade_indices = sorted(index for index in base_manifests if index in upgrades)
    # Do not invoke the CLI: its ordinary JSON summary contains private ids.
    # The pinned file is an array; the public driver accepts request JSONL.
    normalized = Path(args.out) / 'replay-requests.jsonl'
    normalized.parent.mkdir(parents=True, exist_ok=True)
    with normalized.open('w') as stream:
        for index in indices:
            row = inputs[index]
            stream.write(driver.canonical({'max_work': 100_000_000, **row,
                                           'format': driver.CONTRACT}) + '\n')
    report = run_fresh(driver, args.binary, normalized, args.out, 4)
    execution = report['execution']
    if execution['workers'] != 4 or execution['fresh'] is not True:
        raise ValueError('replay requires four fresh workers')
    if execution['worker_exits'] != {str(i): 0 for i in range(4)}:
        raise ValueError('replay has missing or failed workers')
    if report['config']['binary_sha256'] != sha(args.binary):
        raise ValueError('binary_sha256 changed during replay')
    failures, outcomes, work = [], Counter(), 0
    work_cases = []
    verified_corrections = []
    cases = report['cases']
    if len(cases) != len(indices):
        raise ValueError('replay needs every selected case')
    # Preserve original corpus positions for raw STEP and correction checks.
    report['replay_indices'] = indices
    # The direct STEP audit independently checks exactly these new answers.
    report['upgrade_indices'] = upgrade_indices
    write(Path(args.out) / 'run.json', report)
    for index, case in zip(indices, cases):
        accepted = metadata[index]
        directory = Path(case['path'])
        request = read(directory / 'request.json')
        request.pop('out', None)
        fingerprint = driver.digest(request)
        response = case['response']
        manifest = response['manifest']
        cause = None
        if fingerprint != accepted['request_sha256']:
            cause = 'request identity'
        answer_sha256 = release.answer_digest(driver, directory)
        if cause is None and answer_sha256 != accepted['accepted_answer_sha256']:
            manifest_sha256 = driver.digest(driver.comparable(manifest))
            if manifest_sha256 != accepted['accepted_manifest_sha256']:
                cause = 'manifest'
            else:
                changed = [name for name in release.OUTPUTS
                           if (sha(directory / name) if (directory / name).exists() else None)
                           != accepted['accepted_files'][name]]
                cause = 'STEP bytes' if changed == ['answer.step'] else 'output bytes'
        if index in base_manifests and index in upgrades:
            written = [name for name in release.OUTPUTS if (directory / name).exists()]
            if upgrade_manifest_errors(manifest, base_manifests[index], upgrades[index],
                                       inputs[index].get('problem'), written):
                cause = 'unapproved upgrade'
        elif index in base_manifests and correction_manifest_errors(
                manifest, corrections[index], driver, base_manifests[index]):
            cause = 'unapproved correction manifest field'
        if (cause is None and index in base_manifests
                and corrections.get(index, {}).get('class') == LAST_DIGIT_CLASS):
            if index not in base_directories or last_digit_errors(
                    base_directories[index], directory, driver):
                cause = 'answer changed beyond its last stated digit'
        if response['work']['used'] > 100_000_000:
            cause = 'above 100M work'
        if case['cache_hit'] or manifest.get('status') == 'internal-error':
            cause = 'cache or internal error'
        if cause:
            failures.append({'index': index, 'request_sha256': fingerprint, 'cause': cause,
                             'answer_sha256': answer_sha256})
        elif index in base_manifests:
            verified_corrections.append(index)
        refusal = manifest.get('refusal')
        outcomes[refusal['code'] if refusal else manifest['status']] += 1
        work += response['work']['used']
        work_cases.append({'index': index, 'used': response['work']['used']})
    summary = {'workers': execution['workers'], 'fresh': execution['fresh'],
               'commit': plan['identity']['commit'],
               'family': selected['family'], 'indices': indices,
               'cache_hits': sum(c['cache_hit'] for c in cases),
               'selected_cases': len(cases), 'elapsed_ms': (time.monotonic() - started) * 1000,
               'outcomes': dict(outcomes), 'work_total': work, 'work_cases': work_cases,
               'verified_corrections': verified_corrections,
               'equal_cases': len(cases) - len(failures), 'failures': failures}
    write(args.summary, summary)
    print(json.dumps(summary, sort_keys=True))
    if failures:
        raise ValueError('replay: changed earlier answers')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('select')
    p.add_argument('--base', required=True)
    p.add_argument('--head', default='')
    p.add_argument('--event', default='pull_request')
    p.add_argument('--ref', default='')
    p.add_argument('--run-id', default=os.environ.get('GITHUB_RUN_ID', 'local'))
    p.add_argument('--run-attempt', default=os.environ.get('GITHUB_RUN_ATTEMPT', '1'))
    p.add_argument('--release-evidence', action='store_true')
    p.add_argument('--release-dry-run', action='store_true')
    p.add_argument('--github-output')
    p.add_argument('--out', required=True)
    p.set_defaults(function=prepare)
    p = commands.add_parser('build-base')
    p.add_argument('--plan', required=True)
    p.add_argument('--source', required=True)
    p.add_argument('--out', required=True)
    p.set_defaults(function=build_base)
    p = commands.add_parser('evidence')
    p.add_argument('--plan', required=True)
    p.add_argument('--job', required=True)
    p.add_argument('--partition', type=int, default=1)
    p.add_argument('--result', default='success')
    p.add_argument('--binary')
    p.add_argument('--summary')
    p.add_argument('--out', required=True)
    p.set_defaults(function=evidence)
    p = commands.add_parser('gate')
    p.add_argument('--plan', required=True)
    p.add_argument('--evidence', required=True)
    p.add_argument('--jobs', required=True)
    p.set_defaults(function=gate)
    p = commands.add_parser('replay')
    p.add_argument('--plan', required=True)
    p.add_argument('--partition', type=int, required=True)
    p.add_argument('--base-binary')
    p.add_argument('--inputs', required=True)
    p.add_argument('--binary', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--summary', required=True)
    p.set_defaults(function=replay)
    args = parser.parse_args()
    try:
        args.function(args)
    except (ValueError, OSError, KeyError) as error:
        print(str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
