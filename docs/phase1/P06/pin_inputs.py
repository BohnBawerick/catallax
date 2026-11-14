"""Static P06 manifest generator; reads tracked public text only; executes no kernel."""
import copy
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASELINE = 'f4db345f2e7c51835d1ebeba9f6446d117144117'
OUT = ROOT / 'docs/phase1/P06/inputs.json'
H2H = '.claude/skills/verify-catallax/scripts/head-to-head/'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def read(path):
    return json.loads((ROOT / path).read_text())

def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2) + '\n').encode()

names = ['P1-bracket', 'P2-collar', 'P3-filleted-boss', 'P4-shelled-tray',
         'P5-flanged-housing', 'P6-lattice-bracket']
later = ['field-driven bracket with gyroid infill', 'helical spring and threaded boss',
         'lofted duct', 'topology-optimised bracket']
sources = ['LICENSE', 'NOTICE', 'docs/TARGET-PARTS.md', 'parts/README.md',
           'parts/edit-checks.json', 'tools/research/m5_scoreboard.py',
           'tools/research/r8_step_census.py', 'docs/milestones/M5-scoreboard.md',
           'docs/milestones/M5-scoreboard.json', 'docs/milestones/R8-public-step-census.md',
           'docs/milestones/R8-public-step-census.json', 'docs/milestones/R8-public-step-census.csv',
           H2H+'README.md', H2H+'reference.json', H2H+'tools.json', H2H+'targets.py',
           H2H+'python_tools.py', H2H+'pali_parts.dart']
ordered = []
for i in range(1, 11):
    row = {'order': i, 'id': f'P{i}', 'license': 'MIT (repository LICENSE)',
           'chain_policy': 'Construct full original, then full edited recipe; no successful prefix substitution.',
           'additional_edit_chain': None}
    if i <= 6:
        path = 'parts/'+names[i-1]+'.json'
        sources.append(path)
        recipe = read(path)
        row.update(recipe=path, construction_operations=[{'id':s['id'], 'op':s['op']} for s in recipe['steps']], result=recipe['result'])
        row['chain'] = [{'ordinal': 0, 'intent': 'construct', 'source': path}]
        if i <= 5:
            case = H2H+'cases/'+names[i-1]+'.json'
            sources.append(case)
            d = read(case)
            assert sha((ROOT/path).read_bytes()) == d['recipe_sha256']
            edit = d['edit']
            changed = copy.deepcopy(recipe)
            step = next(s for s in changed['steps'] if s['id'] == edit['step'])
            target = step
            for key in edit['path'][:-1]:
                target = target[key]
            old = target[edit['path'][-1]]
            target[edit['path'][-1]] = edit['value']
            row['chain'].append({'ordinal': 1, 'intent': 'rebuild complete edited recipe', 'edit':edit,
                'old_value':old, 'edited_recipe_canonical_sha256':sha(canonical(changed)), 'definition_source':case,
                'reference_edit':d['reference']['edit']})
            row['reference_intent'] = 'Same scalar value, placement, complete target and result as committed reference case; operation graphs may differ; equivalence unexecuted.'
            row['status'] = 'one committed edit pinned; longer edit sequence undefined'
        else:
            row.update(status='construction pinned; edit and reference implementation missing',
                       missing=['P6 parameter edit and ordered continuation', 'reference lattice/Boolean full-target definition and edit'])
    else:
        row.update(name=later[i-7], recipe=None, chain=[], status='named intent only',
                   missing=['fixed recipe, dimensions, placements and operation order', 'parameter edit and ordered continuation', 'reference full-target operation and edit definitions'])
    ordered.append(row)
score = read('docs/milestones/M5-scoreboard.json')
census = {r['model_id']:r for r in csv.DictReader((ROOT/'docs/milestones/R8-public-step-census.csv').open())}
for case in score['groups']['real']['cases']:
    if case['status'] != 'answered':
        continue
    model = case['name'].split()[-1]
    ordered.append({'order':len(ordered)+1, 'id':'ABC '+model, 'model_id':model,
        'selection':'ten answering M5 controls in published order; excludes refused candidate 00000549',
        'source':'ABC chunk 0000 STEP v00; see R8 provenance',
        'source_geometry_sha256':None, 'source_geometry_availability':'not accessed; individual raw hashes absent from inspected public scoreboard/census',
        'license':'dataset MIT grant reported by R8; individual source exceptions unresolved; no raw redistribution',
        'historical_m5_status':case['status'], 'historical_export_step_sha256':case['step_sha256'],
        'census_observation':census[model],
        'chain':[{'ordinal':0,'intent':'import complete unchanged STEP; definition only, not executed'}],
        'reference_intent':'same immutable source bytes, placement and edit in both kernels; currently not executable',
        'missing':['public individual source-byte hash/authorized input', 'specific ordered post-import edit definition', 'reference import/edit definition'],
        'status':'identity pinned; edit comparison blocked'})
manifest = {'format':'catallax-phase1-p06-input-pin/1', 'baseline':BASELINE,
 'repository':'BohnBawerick/catallax', 'evidence_kind':'static source inspection; no new geometry execution',
 'canonical_edited_recipe_encoding':'UTF-8 json.dumps(sort_keys=True, indent=2) plus LF; informational hash, not a copied recipe',
 'reference_release_pin':read(H2H+'reference.json'),
 'census_archive_pin':{k:v for k,v in read('docs/milestones/R8-public-step-census.json').items() if k.startswith('archive_') or k in ('seed','licence_scope')},
 'ordered_inputs':ordered,
 'source_files':[{'path':p, 'sha256':sha((ROOT/p).read_bytes()), 'git_revision':BASELINE} for p in sorted(set(sources))]}
assert len(ordered) == 20
assert len({r['id'] for r in ordered}) == 20
OUT.write_bytes(canonical(manifest))
print('Pinned 20 ordered identities, 5 edits and',len(manifest['source_files']),'source hashes; static inspection only.')
