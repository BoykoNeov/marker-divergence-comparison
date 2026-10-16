import re
D='M:/claud_projects/article/yasen/'
def read(f):
    seqs={};k=None
    for l in open(D+f):
        l=l.rstrip('\n').strip()
        if l.startswith('>'): k=l[1:].strip(); seqs[k]=''
        elif k: seqs[k]+=l.upper()
    return seqs
def species(n):
    p=re.split(r'[_ ]',n)
    return '_'.join(p[:3]) if p[0]=='Plagiorhynchus' else '_'.join(p[:2])
def genus(sp):  # subgenus counts as the genus-level unit
    p=sp.split('_'); return '_'.join(p[:2]) if p[0]=='Plagiorhynchus' else p[0]
def pdist(a,b,cols=None):
    d=n=0
    idx=cols if cols is not None else range(len(a))
    for i in idx:
        x,y=a[i],b[i]
        if x in 'ACGT' and y in 'ACGT':
            n+=1; d+= x!=y
    return d/n
def complete_cols(seqs):
    L=len(next(iter(seqs.values())))
    return [i for i in range(L) if all(s[i] in 'ACGT' for s in seqs.values())]
