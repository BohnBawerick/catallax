"""Recompute the P02 comparison from retained driver output; no kernel runs."""
from pathlib import Path
import copy, hashlib, json
E=Path(__file__).resolve().parent
inventory=json.loads((E/'inventory.json').read_text())
def comparable(response):
    r=copy.deepcopy(response)
    r.pop('wall_ms',None)
    r['manifest'].pop('stats',None)
    for operation in r['operation_log']['operations']:
        operation.pop('wall_ms',None)
    return r
rows=[]
for suite,offset,count in [('m10',0,14),('later-cuts',14,13)]:
    runs=[json.loads((E/f'{suite}-workers-{n}/run.json').read_text()) for n in [1,4]]
    assert all(len(r['cases'])==count for r in runs)
    for i,(a,b) in enumerate(zip(*(r['cases'] for r in runs),strict=True)):
        same=comparable(a['response'])==comparable(b['response'])
        assert not a['cache_hit'] and not b['cache_hit']
        artifacts=[]
        for case in [a,b]:
            # Resolve the retained case by key, without relying on its original absolute path.
            n=1 if case is a else 4
            directory=E/f'{suite}-workers-{n}'/'cases'/case['key']
            files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in directory.iterdir() if p.is_file()}
            artifacts.append({'worker_count':n,'files':files,'step_present':(directory/'answer.step').exists()})
            assert case['response']['work']['limit']==20_000_000
        step_equal=artifacts[0]['files'].get('answer.step')==artifacts[1]['files'].get('answer.step')
        r=a['response']; rows.append({'id':inventory[offset+i]['id'],'name':inventory[offset+i].get('name'), 'status':r['manifest']['status'],'refusal':r['manifest'].get('refusal'),'work':r['work'],'response_equal_excluding_diagnostic_times_and_manifest_stats':same,'step_bytes_equal':step_equal,'artifacts':artifacts})
assert all(r['response_equal_excluding_diagnostic_times_and_manifest_stats'] and r['step_bytes_equal'] for r in rows)
(E/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps({'cases':len(rows),'comparison_equal':True,'m10_work':sum(r['work']['used'] for r in rows[:14]),'later_cuts_work':sum(r['work']['used'] for r in rows[14:])},indent=2))
