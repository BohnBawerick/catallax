"""Extract public baseline recipes; changes mirror cited tests or explicit selections."""
import json, copy
from pathlib import Path
root=Path(__file__).resolve().parents[3]
out=Path(__file__).resolve().parent/'inputs'
def load(p): return json.loads((root/p).read_text())
def save(n,v): (out/(n+'.json')).write_text(json.dumps(v,indent=2)+'\n')
def side(step,s): return {'step':step,'role':{'side':{'segment':s}}}
def edge(a,b): return {'between':[a,b]}
rows=[json.loads(l) for l in (root/'parts/fillets/cases.jsonl').read_text().splitlines()]
def case(n): return copy.deepcopy(next(r['recipe'] for r in rows if r['id']==n))
v=case('l-mixed-all'); save('mixed-all-control',v)
v['steps'][3]['edges']=[edge({'step':'base','role':'end'},side('wall','y_max')),edge({'step':'base','role':'end'},side('base','x_min'))]; save('mixed-unselected',v)
save('partial-oblique',case('partial-oblique'))
v=load('crates/catallax-features/tests/fixtures/bored-shell.json'); save('bored-cap-control',v)
v['steps'][2]['open']=[{'face':side('stock','x_min')}]; save('bored-side-opening',v)
v=load('crates/catallax-features/tests/fixtures/draft-oblique.json'); save('oblique-solid-control',v)
v['steps'][0]['loops'].append([{'kind':'circle','id':'bore','centre':['0','0'],'radius':'2'}]); save('hollow-oblique-draft',v)
v={'format':'catallax-recipe/3','name':'plane bevel','steps':[{'id':'base','op':'block','size':['10','12','8'],'pose':{'rotation':['1','0','0','0'],'translation':['0','0','0']}},{'id':'bevel','op':'chamfer','body':'base','distance':'2','edges':[edge(side('base','x_min'),side('base','y_min'))]}],'result':'bevel'}
save('bevel-single-control',v)
v['steps'][1]['edges'] += [edge(side('base','x_min'),{'step':'base','role':'start'}),edge(side('base','y_min'),{'step':'base','role':'start'})]; save('bevel-three-edge-corner',v)
plane={'rotation':['1','0','0','0'],'translation':['0','0','0']}
v={'format':'catallax-recipe/3','name':'collar plain hole','steps':[{'id':'section','op':'sketch','plane':plane,'loops':[[{'kind':'line','id':i,'to':p} for i,p in [('bottom',['12','0']),('right',['12','12']),('top',['6','12']),('left',['6','0'])]]]},{'id':'ring','op':'revolve','profile':'section','axis':{'point':['0','0','0'],'direction':['0','1','0']},'angle':'360'},{'id':'collar','op':'chamfer','body':'ring','distance':'1','edges':[edge(side('ring','right'),side('ring','bottom')),edge(side('ring','right'),side('ring','top'))]},{'id':'drilled','op':'hole','target':'collar','at':[{'id':'radial-0','point':['21/2','12','0'],'direction':['1','-1','1']}],'diameter':'2','depth':'through','tip':'flat'}],'result':'drilled'}
save('skew-hole',v)
v['steps'][3]['at'][0].update(point=['9','12','0'],direction=['0','-1','0']); save('radial-hole-control',v)
