"""Read the native pcurve endpoints and their physical plane/patch mapping.

This oracle handles only the authored linear cube controls. It does not ask an
importer to replace a planar pcurve with a projection of its 3D edge.
"""
from pathlib import Path
from fractions import Fraction as Q
import re,json,sys

def fields(value):
 out=[];depth=0;start=0;quoted=False
 for i,c in enumerate(value):
  if c=="'":quoted=not quoted
  elif not quoted:
   if c=='(':depth+=1
   elif c==')':depth-=1
   elif c==',' and depth==0:out.append(value[start:i]);start=i+1
 return out+[value[start:]]
def refs(s):return [int(v) for v in re.findall(r'#(\d+)',s)]
def parse(expr):return fields(expr[expr.index('(')+1:-1])
def point(es,i):return [Q(v) for v in fields(parse(es[i])[1][1:-1])]
def add(a,b):return [x+y for x,y in zip(a,b)]
def scale(a,s):return [x*s for x in a]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def check(root,name):
 text=(root/name/'answer.step').read_text();es={int(i):v for i,v in re.findall(r'^#(\d+)=(.*);$',text,re.M)};meta=json.loads((root/name/'control.json').read_text())
 face=next(parse(e) for e in es.values() if e.startswith('ADVANCED_FACE(') and parse(e)[0]=="'"+meta['face']+"'");surface=refs(face[2])[0]
 def xyz(uv):
  s=es[surface]
  if s.startswith('PLANE('):
   axis=parse(es[refs(parse(s)[1])[0]]);o=point(es,refs(axis[1])[0]);n=point(es,refs(axis[2])[0]);x=point(es,refs(axis[3])[0]);assert sum(v*v for v in n)==sum(v*v for v in x)==1;return add(o,add(scale(x,uv[0]),scale(cross(n,x),uv[1])))
  m=re.search(r'B_SPLINE_SURFACE\(1,1,\(\((#[\d,#]+)\),\((#[\d,#]+)\)\)',s);assert m,s
  ids=refs(m[1])+refs(m[2]);ps=[point(es,i) for i in ids];u=(uv[0]+2)/4;v=(uv[1]+2)/4
  return [ps[0][a]*(1-u)*(1-v)+ps[1][a]*(1-u)*v+ps[2][a]*u*(1-v)+ps[3][a]*u*v for a in range(3)]
 bound=parse(es[refs(face[1])[0]]);loop=parse(es[refs(bound[1])[0]]);uses=refs(loop[1]);ends=[]
 for i in meta['uses']:
  use=parse(es[uses[i]]);edge=parse(es[refs(use[3])[0]]);sc=parse(es[refs(edge[3])[0]]);pc=next(parse(es[p]) for p in refs(sc[2]) if refs(parse(es[p])[1])[0]==surface)
  rep=parse(es[refs(pc[2])[0]]);trim=parse(es[refs(rep[1])[0]]);basis=es[refs(trim[1])[0]];controls=refs(re.search(r'B_SPLINE_CURVE\(1,(\(#[\d,#]+\))',basis)[1]);uv=[point(es,j) for j in controls];assert len(uv)==2
  if edge[4]=='.F.':uv.reverse()
  if use[4]=='.F.':uv.reverse()
  ends.append([xyz(p) for p in uv])
 delta=[abs(a-b) for a,b in zip(ends[0][1],ends[1][0])];assert sum(v!=0 for v in delta)<=1;gap=sum(delta);tol=Q(re.search(r'source_tolerance_mm=([^;]+);',text)[1]);return {'control':name,'actual_step_join_gap_mm':str(gap),'tolerance_mm':str(tol),'exceeds_tolerance':gap>tol}

def main():
 import argparse
 parser=argparse.ArgumentParser(description="Check native join controls from the actual STEP numerals.")
 parser.add_argument("root",type=Path);parser.add_argument("--baseline",action="store_true");parser.add_argument("--out",type=Path);args=parser.parse_args()
 rows=[check(args.root,"valid")];assert rows[0]["actual_step_join_gap_mm"]=="0"
 if args.baseline:
  rows.extend(check(args.root,n) for n in ["source-gap","rounded-gap"])
  assert [r["actual_step_join_gap_mm"] for r in rows]==["0","9/500","3/250"]
  assert all(r["exceeds_tolerance"] for r in rows[1:])
 else:
  for name,stage in [("source-gap","source admission:"),("rounded-gap","decimal admission:")]:
   response=json.loads((args.root/name/"response.json").read_text());answer=response["result"]
   assert answer["status"]=="refused" and stage in answer["reason"] and "trim.native-join-position" in answer["reason"]
   assert not response["work"]["exhausted"] and response["work"]["limit"]==20_000_000
   assert not (args.root/name/"answer.step").exists()
   rows.append({"control":name,"refused":True,"stage":stage,"no_stale_step":True})
 result={"method":"actual STEP numerals in the authored cube controls","cases":rows}
 if args.out:args.out.write_text(json.dumps(result,indent=2)+"\n")
 print(json.dumps(result,indent=2))
if __name__=="__main__":main()
