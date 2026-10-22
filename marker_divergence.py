"""Compare sequence divergence at two DNA markers sampled from the same species.

Reads two aligned FASTA files, calculates uncorrected p-distances (pairwise
deletion of gaps and ambiguous positions) for the species present in both
alignments, prints a short summary and draws a three-panel figure (PNG, PDF and
SVG): site variability along each alignment, and marker 1 against marker 2
distances for every comparison, with a close-up of the low-divergence corner.
See README.md for details.
"""
import argparse
import collections
import itertools
import json
import math
import re
import statistics as st
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

SPECIES, GENUS, BETWEEN = 'species', 'genus', 'between'
DRAW_ORDER = [BETWEEN, GENUS, SPECIES]

# Settings that can be given on the command line or in a --settings JSON file.
DEFAULTS = {
    'names': ['Marker 1', 'Marker 2'],  # short marker names, used in axis labels and the summary
    'titles': None,                     # titles of the two site-variability panels (default: names)
    'species_parts': 2,                 # header parts that make up the species name
    'subgenus_genera': [],              # genera whose headers read Genus_Subgenus_species_...
    'window': 15,                       # moving-average window for the site-variability profile
    'position_max': None,               # x-axis end of the site-variability panels (default: automatic)
    'scatter_max': None,                # axis end of the all-comparisons panel (default: automatic)
    'zoom': 0.05,                       # axis end of the close-up panel
    'out_name': None,                   # base name of the figure files (default: <name1>_vs_<name2>_divergence)
    'labels': [],                       # labelled points in the close-up panel (settings file only)
}


# ---------------------------------------------------------------- sequences --

def read_fasta(path):
    """Return {header: sequence} in file order; sequences upper-cased."""
    seqs, name = {}, None
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].strip(); seqs[name] = ''
            elif name:
                seqs[name] += line.upper()
    if not seqs:
        sys.exit(f'No sequences found in {path}')
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) > 1:
        sys.exit(f'{path}: sequences differ in length ({min(lengths)}-{max(lengths)}); is the file aligned?')
    return seqs


class Naming:
    """Turns sequence headers into species and genus-level units."""

    def __init__(self, parts, subgenus_genera):
        self.parts, self.subgenus = parts, set(subgenus_genera)

    def species(self, header):
        """Species from a header: the first `parts` parts, plus one for a genus with subgenera."""
        p = re.split(r'[_ ]', header)
        return '_'.join(p[:self.parts + (p[0] in self.subgenus)])

    def genus(self, sp):
        """Genus-level unit of a species; for a genus with subgenera, the subgenus."""
        p = sp.split('_')
        return '_'.join(p[:2]) if p[0] in self.subgenus else p[0]


def pdist(a, b):
    """Uncorrected p-distance, ignoring positions with a gap or ambiguity in either sequence.
    NaN when the two sequences share no comparable position."""
    d = n = 0
    for x, y in zip(a, b):
        if x in 'ACGT' and y in 'ACGT':
            n += 1; d += x != y
    return d / n if n else math.nan


def by_species(seqs, naming):
    groups = collections.defaultdict(list)
    for name in seqs:
        groups[naming.species(name)].append(name)
    return groups


# ---------------------------------------------------------------- analysis --

def comparisons(s1, s2, sp1, sp2, shared, naming):
    """One point per species pair (mean between-species distance) and per species
    with >= 2 sequences at both markers (mean within-species distance)."""
    def between(s, A, B): return st.mean(pdist(s[a], s[b]) for a in A for b in B)
    def within(s, A): return st.mean(pdist(s[a], s[b]) for a, b in itertools.combinations(A, 2))

    pts = []  # (distance at marker 1, distance at marker 2, level, species involved)
    for a, b in itertools.combinations(shared, 2):
        lvl = GENUS if naming.genus(a) == naming.genus(b) else BETWEEN
        pts.append((between(s1, sp1[a], sp1[b]), between(s2, sp2[a], sp2[b]), lvl, (a, b)))
    for a in shared:
        if len(sp1[a]) > 1 and len(sp2[a]) > 1:
            pts.append((within(s1, sp1[a]), within(s2, sp2[a]), SPECIES, (a,)))
    kept = []
    for p in pts:
        if math.isnan(p[0]) or math.isnan(p[1]):
            print(f'Warning: no comparable positions for {" vs ".join(p[3])}; comparison left out')
        else:
            kept.append(p)
    return kept


def site_profile(seqs, groups, shared):
    """Variability per aligned position, using the first sequence of each shared species.
    Returns (variability array, n variable sites, n parsimony-informative sites, length)."""
    reps = [seqs[groups[sp][0]] for sp in shared]
    length = len(reps[0])
    var, n_var, n_inf = [], 0, 0
    for i in range(length):
        col = [r[i] for r in reps if r[i] in 'ACGT']
        if not col:
            var.append(np.nan); continue
        c = collections.Counter(col)
        var.append(1 - c.most_common(1)[0][1] / len(col))
        if len(c) > 1: n_var += 1
        if sum(v >= 2 for v in c.values()) >= 2: n_inf += 1
    return np.array(var), n_var, n_inf, length


def summary(pts, prof, shared, only, names, level_names):
    n1, n2 = names
    print(f'Species in both alignments: {len(shared)}')
    for g, sps in zip(names, only):
        if sps: print(f'  only in {g}:', ', '.join(sps))
    above = sum(y > x for x, y, *_ in pts); below = sum(y < x for x, y, *_ in pts)
    print(f'Comparisons: {len(pts)}; {n2} higher {above}, {n1} higher {below}, equal {len(pts) - above - below}')
    for lvl in reversed(DRAW_ORDER):
        xs = [p for p in pts if p[2] == lvl]
        if not xs:
            print(f'  {level_names[lvl]}: n = 0'); continue
        print(f'  {level_names[lvl]}: n = {len(xs)}, median {n1} {st.median(p[0] for p in xs):.3f}, '
              f'{n2} {st.median(p[1] for p in xs):.3f}, {n2} higher in {sum(p[1] > p[0] for p in xs)}')
    for g, (_, nv, ni, L) in zip(names, prof):
        print(f'{g}: {nv} variable and {ni} parsimony-informative sites of {L} aligned positions')


# ---------------------------------------------------------------- figure --

INK, INK2, GRID, RULE = '#0b0b0b', '#52514e', '#e1e0d9', '#898781'
COLOURS = {SPECIES: '#2a78d6', GENUS: '#eb6834', BETWEEN: '#1baf7a'}


def short(sp):
    p = sp.split('_'); return f'{p[0][0]}. {p[-1]}'


def scatter_panel(a, pts, lim, big, names):
    a.plot([0, lim], [0, lim], color=RULE, lw=0.8, ls=(0, (4, 3)), zorder=1)
    for lvl in DRAW_ORDER:
        P = [p for p in pts if p[2] == lvl and p[0] <= lim and p[1] <= lim]
        size = (12 if lvl == BETWEEN else 20) if big else 26
        a.scatter([p[0] for p in P], [p[1] for p in P], s=size,
                  color=COLOURS[lvl], edgecolor='white', linewidth=0.6, zorder=3)
    pad = lim * 0.03
    a.set_xlim(-pad, lim); a.set_ylim(-pad, lim); a.set_aspect('equal')
    a.set_xlabel(f'{names[0]} uncorrected p-distance'); a.set_ylabel(f'{names[1]} uncorrected p-distance')
    a.grid(color=GRID, lw=0.5); a.set_axisbelow(True)
    for s in ['top', 'right']: a.spines[s].set_visible(False)
    if big:
        a.text(0.04, 0.93, f'{names[1]} more divergent', transform=a.transAxes, color=INK2, fontsize=7.5,
               style='italic')
        a.text(0.96, 0.04, f'{names[0]} more divergent', transform=a.transAxes, ha='right', color=INK2,
               fontsize=7.5, style='italic')


def draw(pts, prof, cfg, level_names, out_dir):
    plt.rcParams.update({'pdf.fonttype': 42, 'svg.fonttype': 'none', 'font.family': 'Arial', 'font.size': 8,
                         'axes.edgecolor': RULE, 'axes.linewidth': 0.6, 'xtick.color': INK2,
                         'ytick.color': INK2, 'axes.labelcolor': INK, 'xtick.major.width': 0.6,
                         'ytick.major.width': 0.6})
    fig = plt.figure(figsize=(6.85, 6.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 3.3], hspace=0.55, left=0.09, right=0.97, top=0.95, bottom=0.11)

    # A: site-variability profiles
    window = cfg['window']
    smooth = []
    for v, *_ in prof:
        w = min(window, len(v))   # a window longer than the alignment would change the curve's length
        smooth.append(np.convolve(np.nan_to_num(v), np.ones(w) / w, mode='same'))
    pos_max = cfg['position_max'] or math.ceil(max(L for *_, L in prof) / 100) * 100
    var_max = max(0.45, max(s.max() for s in smooth) + 0.05)
    for k, (title, (v, nv, ni, L), sm) in enumerate(zip(cfg['titles'], prof, smooth)):
        ax = fig.add_subplot(gs[k])
        vv = np.nan_to_num(v)
        x = np.arange(1, L + 1)
        ax.fill_between(x, 0, sm, color=INK2, alpha=0.35, lw=0)
        ax.plot(x, sm, color=INK2, lw=1)
        ax.vlines(x[vv > 0], -0.06, -0.015, color=INK2, lw=0.4)   # tick at every variable position
        ax.set_xlim(0, pos_max); ax.set_ylim(-0.07, var_max)
        ax.set_yticks(np.arange(0, var_max, 0.2)); ax.grid(axis='y', color=GRID, lw=0.5); ax.set_axisbelow(True)
        for s in ['top', 'right']: ax.spines[s].set_visible(False)
        ax.text(0.0, 1.08, title, transform=ax.transAxes, fontweight='bold', color=INK, va='bottom')
        ax.text(1.0, 1.08, f'{nv} variable sites, {ni} parsimony-informative, of {L} aligned positions',
                transform=ax.transAxes, ha='right', va='bottom', color=INK2)
        if k == 0:
            ax.text(-0.075, 1.08, 'A', transform=ax.transAxes, fontsize=11, fontweight='bold', va='bottom')
        ax.set_ylabel('Site\nvariability', fontsize=7.5)
    ax.set_xlabel('Alignment position')

    # B: all comparisons; C: close-up of the low-divergence corner
    sub = gs[2].subgridspec(1, 2, wspace=0.28)
    axB, axC = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
    axB.text(-0.2, 1.02, 'B', transform=axB.transAxes, fontsize=11, fontweight='bold', va='bottom')
    axC.text(-0.2, 1.02, 'C', transform=axC.transAxes, fontsize=11, fontweight='bold', va='bottom')
    lim = cfg['scatter_max'] or math.ceil(max(max(p[0], p[1]) for p in pts) * 100) / 100 + 0.01
    zoom, u = cfg['zoom'], lim / 0.33   # u scales the marker of the close-up area to the axis length
    scatter_panel(axB, pts, lim, True, cfg['names'])
    scatter_panel(axC, pts, zoom, False, cfg['names'])
    axB.add_patch(plt.Rectangle((-0.001 * u, -0.001 * u), zoom + 0.001 * u, zoom + 0.001 * u, fill=False, ec=INK2,
                                lw=0.7, zorder=4))
    axB.text(zoom + 0.006 * u, 0.004 * u, 'C', color=INK2, fontsize=7.5, fontweight='bold')
    axB.set_title('All species pairs', fontsize=8, color=INK, loc='left')
    axC.set_title('Close-up of the low-divergence corner', fontsize=8, color=INK, loc='left')

    n = {l: sum(p[2] == l for p in pts) for l in DRAW_ORDER}
    handles = [Line2D([], [], marker='o', ls='', color=COLOURS[l], markeredgecolor='white', markersize=6,
                      label=f'{level_names[l]} (n = {n[l]})') for l in reversed(DRAW_ORDER) if n[l]]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=7.5, handletextpad=0.3,
               bbox_to_anchor=(0.5, 0.0), columnspacing=1.0)

    labels = {tuple(sorted(lab['species'])): lab for lab in cfg['labels']}
    found = set()
    for p in pts:
        lab = labels.get(tuple(sorted(p[3])))
        if lab is None: continue
        found.add(tuple(sorted(p[3])))
        text = (' vs '.join(short(s) for s in p[3]) if len(p[3]) == 2
                else f'{short(p[3][0])} (within species)')
        axC.annotate(text, (p[0], p[1]), tuple(lab['text_at']), fontsize=6.5, color=INK, style='italic',
                     va='center', annotation_clip=False,
                     arrowprops=dict(arrowstyle='-', color=RULE, lw=0.6, shrinkA=0, shrinkB=3,
                                     relpos=tuple(lab.get('line_from', (0.5, 0.5)))))
    for key in labels.keys() - found:
        print(f'Warning: label for {" vs ".join(key)} matches no comparison; not drawn')

    out_dir.mkdir(parents=True, exist_ok=True)
    for ext, kw in [('png', {'dpi': 600}), ('pdf', {}), ('svg', {})]:
        fig.savefig(out_dir / f'{cfg["out_name"]}.{ext}', **kw)
    print(f'Figure saved to {out_dir / cfg["out_name"]}.png/.pdf/.svg')


# ---------------------------------------------------------------- settings --

def check_settings(cfg, source):
    """Stop with a message if a setting is unknown or has the wrong shape."""
    unknown = set(cfg) - set(DEFAULTS)
    if unknown:
        sys.exit(f'{source}: unknown setting(s): {", ".join(sorted(unknown))}. '
                 f'Known settings: {", ".join(DEFAULTS)}')
    for key in ('names', 'titles'):
        if cfg.get(key) is not None and (not isinstance(cfg[key], list) or len(cfg[key]) != 2):
            sys.exit(f'{source}: "{key}" must be a list of two texts')
    for lab in cfg.get('labels', []):
        if (not isinstance(lab, dict) or set(lab) - {'species', 'text_at', 'line_from'}
                or not 1 <= len(lab.get('species', [])) <= 2 or len(lab.get('text_at', [])) != 2
                or len(lab.get('line_from', [0, 0])) != 2):
            sys.exit(f'{source}: each label needs "species" (one or two names), "text_at" ([x, y]) '
                     f'and optionally "line_from" ([x, y]); problem with: {lab}')


def settings(args):
    """Defaults, overridden by the settings file, overridden by command-line options."""
    cfg = dict(DEFAULTS)
    if args.settings:
        with open(args.settings, encoding='utf-8') as fh:
            from_file = json.load(fh)
        check_settings(from_file, args.settings)
        cfg.update(from_file)
    cfg.update({k: v for k, v in vars(args).items() if k in DEFAULTS and v is not None})
    check_settings({k: v for k, v in cfg.items() if k != 'labels'}, 'command line')
    cfg['titles'] = cfg['titles'] or cfg['names']
    cfg['out_name'] = cfg['out_name'] or '{}_vs_{}_divergence'.format(*cfg['names']).replace(' ', '_')
    return cfg


# ---------------------------------------------------------------- main --

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('aln1', type=Path, help='aligned FASTA file of marker 1 (x-axis)')
    ap.add_argument('aln2', type=Path, help='aligned FASTA file of marker 2 (y-axis)')
    ap.add_argument('--settings', type=Path, help='JSON file with any of the settings below, and labels')
    ap.add_argument('--names', nargs=2, metavar=('NAME1', 'NAME2'), help='short marker names (default: Marker 1, Marker 2)')
    ap.add_argument('--titles', nargs=2, metavar=('TITLE1', 'TITLE2'),
                    help='titles of the site-variability panels (default: the names)')
    ap.add_argument('--species-parts', type=int, help='header parts that make up the species name (default: 2)')
    ap.add_argument('--subgenus-genera', nargs='+', metavar='GENUS',
                    help='genera whose headers read Genus_Subgenus_species_...')
    ap.add_argument('--window', type=int, help='moving-average window, in positions (default: 15)')
    ap.add_argument('--position-max', type=int, help='end of the alignment-position axis (default: automatic)')
    ap.add_argument('--scatter-max', type=float, help='axis end of panel B (default: automatic)')
    ap.add_argument('--zoom', type=float, help='axis end of the close-up panel C (default: 0.05)')
    ap.add_argument('--out', type=Path, default=Path('output'), help='folder for the figure files (default: output)')
    ap.add_argument('--out-name', help='base name of the figure files (default: <name1>_vs_<name2>_divergence)')
    args = ap.parse_args(argv)
    cfg = settings(args)

    naming = Naming(cfg['species_parts'], cfg['subgenus_genera'])
    level_names = {SPECIES: 'Intraspecific, within species', GENUS: 'Interspecific within genus',
                   BETWEEN: 'Interspecific between genera/subgenera' if cfg['subgenus_genera']
                   else 'Interspecific between genera'}
    s1, s2 = read_fasta(args.aln1), read_fasta(args.aln2)
    sp1, sp2 = by_species(s1, naming), by_species(s2, naming)
    shared = sorted(set(sp1) & set(sp2))
    if len(shared) < 2:
        sys.exit(f'Only {len(shared)} species found in both alignments; at least 2 are needed. '
                 f'Check the header names (see README).')
    pts = comparisons(s1, s2, sp1, sp2, shared, naming)
    prof = [site_profile(s1, sp1, shared), site_profile(s2, sp2, shared)]
    summary(pts, prof, shared, [sorted(set(sp1) - set(sp2)), sorted(set(sp2) - set(sp1))], cfg['names'],
            level_names)
    draw(pts, prof, cfg, level_names, args.out)
    return pts, prof


if __name__ == '__main__':
    main()
