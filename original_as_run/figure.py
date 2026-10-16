import sys; sys.path.insert(0,'W:/temp/claude/yasen_viz')
from common import *
import itertools, collections, statistics as st
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUT='W:/temp/claude/yasen_viz/'
g18=read('18Salignment_curated.fas'); g28=read('28Salignment_curated.fas')
def by_sp(s):
    d=collections.defaultdict(list)
    for n in s: d[species(n)].append(n)
    return d
sp18,sp28=by_sp(g18),by_sp(g28)
shared=sorted(set(sp18)&set(sp28))
print('shared species',len(shared), 'only18',sorted(set(sp18)-set(sp28)),'only28',sorted(set(sp28)-set(sp18)))

def mean_between(s,A,B): return st.mean(pdist(s[a],s[b]) for a in A for b in B)
def mean_within(s,A): return st.mean(pdist(s[a],s[b]) for a,b in itertools.combinations(A,2))

pts=[]  # (x18,y28,level,label)
for a,b in itertools.combinations(shared,2):
    lvl='Same genus' if genus(a)==genus(b) else 'Different genera'
    pts.append((mean_between(g18,sp18[a],sp18[b]),mean_between(g28,sp28[a],sp28[b]),lvl,(a,b)))
for a in shared:
    if len(sp18[a])>1 and len(sp28[a])>1:
        pts.append((mean_within(g18,sp18[a]),mean_within(g28,sp28[a]),'Same species',(a,)))

above=sum(1 for x,y,*_ in pts if y>x); below=sum(1 for x,y,*_ in pts if y<x)
print('pairs',len(pts),'28S higher',above,'18S higher',below,'equal',len(pts)-above-below)
for lvl in ['Same species','Same genus','Different genera']:
    xs=[p for p in pts if p[2]==lvl]
    print(lvl,len(xs),'median 18S %.3f 28S %.3f'%(st.median(p[0] for p in xs),st.median(p[1] for p in xs)),
          '28S higher in',sum(p[1]>p[0] for p in xs))
failed=[p for p in pts if p[2]!='Same species' and p[0]<0.0025 and p[1]>p[0]]
print('18S ~0 but 28S >0:',[(p[3],round(p[0],4),round(p[1],4)) for p in failed])
for p in pts:
    if p[2]=='Same species': print('  intra',p[3],round(p[0],4),round(p[1],4))
for p in pts:
    if p[2]!='Same species' and p[1]<p[0]: print('  18S>28S',p[3],round(p[0],4),round(p[1],4))

# --- site variability, one sequence per shared species (equal n for both genes) ---
def profile(s,spd):
    reps=[s[spd[x][0]] for x in shared]
    L=len(reps[0]); var=[];inf=0;nv=0
    for i in range(L):
        col=[r[i] for r in reps if r[i] in 'ACGT']
        c=collections.Counter(col)
        if not col: var.append(np.nan);continue
        var.append(1-c.most_common(1)[0][1]/len(col))
        if len(c)>1: nv+=1
        if sum(v>=2 for v in c.values())>=2: inf+=1
    return np.array(var),nv,inf,L
prof={}
for g,s,spd in [('18S',g18,sp18),('28S',g28,sp28)]:
    prof[g]=profile(s,spd); print(g,'variable',prof[g][1],'informative',prof[g][2],'length',prof[g][3])
np.save(OUT+'pts.npy',np.array([(p[0],p[1]) for p in pts]))
import pickle; pickle.dump(dict(pts=pts,prof=prof,shared=shared,failed=failed),open(OUT+'data.pkl','wb'))
