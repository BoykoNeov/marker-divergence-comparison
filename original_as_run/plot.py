import pickle, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT='W:/temp/claude/yasen_viz/'
d=pickle.load(open(OUT+'data.pkl','rb')); pts=d['pts']; prof=d['prof']
pts=[(x,y,('Different genera/subgenera' if l=='Different genera' else l),k) for x,y,l,k in pts]
plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.family':'Arial','font.size':8,'axes.edgecolor':'#898781','axes.linewidth':0.6,
    'xtick.color':'#52514e','ytick.color':'#52514e','axes.labelcolor':'#0b0b0b','xtick.major.width':0.6,'ytick.major.width':0.6})
INK='#0b0b0b'; INK2='#52514e'; GRID='#e1e0d9'
COL={'Same species':'#2a78d6','Same genus':'#eb6834','Different genera/subgenera':'#1baf7a'}
def it(n): return n.replace('Plagiorhynchus_','P. (').replace('Prosthorhynchus_','Prosthorhynchus) ') if 0 else n
def short(sp):
    p=sp.split('_'); return f'{p[0][0]}. {p[-1]}'

fig=plt.figure(figsize=(6.85,6.2))
gs=fig.add_gridspec(3,1,height_ratios=[1,1,3.3],hspace=0.55,left=0.09,right=0.97,top=0.95,bottom=0.11)
win=15
for k,(g,ax) in enumerate(zip(['18S','28S'],[fig.add_subplot(gs[0]),fig.add_subplot(gs[1])])):
    v,nv,inf,L=prof[g]
    vv=np.nan_to_num(v); sm=np.convolve(vv,np.ones(win)/win,mode='same')
    x=np.arange(1,L+1)
    ax.fill_between(x,0,sm,color=INK2,alpha=0.35,lw=0)
    ax.plot(x,sm,color=INK2,lw=1)
    ax.vlines(x[vv>0],-0.06,-0.015,color=INK2,lw=0.4)  # rug: every variable position
    ax.set_xlim(0,600); ax.set_ylim(-0.07,0.45)
    ax.set_yticks([0,0.2,0.4]); ax.grid(axis='y',color=GRID,lw=0.5); ax.set_axisbelow(True)
    for s in ['top','right']: ax.spines[s].set_visible(False)
    ax.text(0.0,1.08,f'{g} rRNA',transform=ax.transAxes,fontweight='bold',color=INK,va='bottom')
    ax.text(1.0,1.08,f'{nv} variable sites, {inf} parsimony-informative, of {L} aligned positions',
            transform=ax.transAxes,ha='right',va='bottom',color=INK2)
    if k==0: ax.text(-0.075,1.08,'A',transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom')
    ax.set_ylabel('Site\nvariability',fontsize=7.5)
ax.set_xlabel('Alignment position')

sub=gs[2].subgridspec(1,2,wspace=0.28)
ax=fig.add_subplot(sub[0]); ins=fig.add_subplot(sub[1])
ax.text(-0.2,1.02,'B',transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom')
ins.text(-0.2,1.02,'C',transform=ins.transAxes,fontsize=11,fontweight='bold',va='bottom')
order=['Different genera/subgenera','Same genus','Same species']
def panel(a,lim,big):
    a.plot([0,lim],[0,lim],color='#898781',lw=0.8,ls=(0,(4,3)),zorder=1)
    for lvl in order:
        P=[p for p in pts if p[2]==lvl and p[0]<=lim and p[1]<=lim]
        a.scatter([p[0] for p in P],[p[1] for p in P],s=(12 if lvl=='Different genera/subgenera' else 20) if big else 26,
                  color=COL[lvl],edgecolor='white',linewidth=0.6,zorder=3)
    pad=lim*0.03
    a.set_xlim(-pad,lim); a.set_ylim(-pad,lim); a.set_aspect('equal')
    a.set_xlabel('18S uncorrected p-distance'); a.set_ylabel('28S uncorrected p-distance')
    a.grid(color=GRID,lw=0.5); a.set_axisbelow(True)
    for sp in ['top','right']: a.spines[sp].set_visible(False)
    if big: a.text(0.04,0.93,'28S more divergent',transform=a.transAxes,color=INK2,fontsize=7.5,style='italic')
    if big: a.text(0.96,0.04,'18S more divergent',transform=a.transAxes,ha='right',color=INK2,fontsize=7.5,style='italic')
panel(ax,0.33,True); ax._quad=True
z=0.05
panel(ins,z,False)
ax.add_patch(plt.Rectangle((-0.001,-0.001),z+0.001,z+0.001,fill=False,ec=INK2,lw=0.7,zorder=4))
ax.text(z+0.006,0.004,'C',color=INK2,fontsize=7.5,fontweight='bold')
ax.set_title('All species pairs',fontsize=8,color=INK,loc='left')
ins.set_title('Close-up of the low-divergence corner',fontsize=8,color=INK,loc='left')
n={l:sum(p[2]==l for p in pts) for l in order}
hi={l:sum(p[2]==l and p[1]>p[0] for p in pts) for l in order}
handles=[Line2D([],[],marker='o',ls='',color=COL[l],markeredgecolor='white',markersize=6,
         label=f'{l} (n = {n[l]})') for l in ['Same species','Same genus','Different genera/subgenera']]
fig.legend(handles=handles,loc='lower center',ncol=3,frameon=False,fontsize=7.5,handletextpad=0.3,
           bbox_to_anchor=(0.5,0.0),columnspacing=1.5)
ann={('Centrorhynchus_buteonis','Centrorhynchus_microcephalus'):(-0.0006,0.0205),
     ('Polymorphus_minutus','Polymorphus_obtusus'):(0.025,0.0050),
     ('Southwellina_hispida',):(0.025,0.0013),
     ('Pseudocorynosoma_anatarium','Pseudocorynosoma_tepehuanesi'):(0.025,0.0088),
     ('Profilicollis_altmani','Profilicollis_bullocki'):(0.025,0.0126)}
def short(sp):
    p=sp.split('_'); return f'{p[0][0]}. {p[-1]}'
for p in pts:
    if p[3] in ann:
        tx,ty=ann[p[3]]
        lab=' vs '.join(short(s) for s in p[3]) if len(p[3])==2 else f'{short(p[3][0])} (within species)'
        ins.annotate(lab,(p[0],p[1]),(tx,ty),fontsize=6.5,color=INK,style='italic',va='center',annotation_clip=False,
                     arrowprops=dict(arrowstyle='-',color='#898781',lw=0.6,shrinkA=0,shrinkB=3,
                                     relpos=(0.03,0) if 'buteonis' in lab else (0.5,0.5)))
fig.savefig(OUT+'Fig_18S_vs_28S_divergence.png',dpi=600)
fig.savefig(OUT+'Fig_18S_vs_28S_divergence.pdf')
fig.savefig(OUT+'Fig_18S_vs_28S_divergence.svg')
print('ok')
