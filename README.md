# marker-divergence-comparison

Compare how fast two DNA markers diverge, using the same species for both.

`marker_divergence.py` reads two aligned FASTA files, for example 18S and 28S
rDNA, or COI and ITS2. It calculates uncorrected p-distances between and within
species at each marker. It plots one marker's distances against the other's and
shows how variability changes along each alignment. It also prints a short
summary.

It was written for, and used to draw, a figure in a study of acanthocephalan
parasites of birds (citation below). It works for any group and any pair of
DNA markers whose sequence names follow the rule described below.

## What it does

1. **Reads the two aligned FASTA files.** All sequences in a file must have the
   same length, gaps included.
2. **Groups the sequences into species**, using the name in each header (see
   [Sequence names](#sequence-names)).
3. **Keeps only the species present in both alignments.** It lists the species
   that occur in only one of them.
4. **Calculates uncorrected p-distances.** These are the proportion of compared
   positions at which two sequences differ. A position is skipped for a given
   pair if either sequence has a gap or an ambiguous base there (pairwise
   deletion). Only A, C, G and T are compared.
   - **Between species:** the mean over all combinations of the two species'
     sequences.
   - **Within a species:** the mean over all pairs of that species' sequences.
     This is calculated only for species with at least two sequences at both
     markers.
   - A pair of sequences with no position in common that can be compared is
     reported and left out.
5. **Assigns each comparison to a level:** same species, same genus, or
   different genera.
6. **Measures variability along each alignment**, using one sequence per
   species: the first one listed in the file. At each position, variability is
   one minus the frequency of the most common base. It is smoothed with a
   moving average (15 positions by default).
   - A position is *variable* if more than one base occurs there.
   - It is *parsimony-informative* if at least two bases each occur in at least
     two sequences.
7. **Prints a summary and draws the figure.**

## Sequence names

The species is read from the start of each FASTA header. The parts of the name
are separated by `_` or spaces.

- By default the species is the first two parts: `>Genus_species_anything_else`
  (for example `>Polymorphus_minutus_isolate1`). Every sequence of a species
  must start with exactly the same two parts, in both files.
- `--species-parts N` uses the first N parts instead, for example 3 for
  subspecies.
- `--subgenus-genera GENUS ...` handles genera whose headers include a subgenus:
  `>Genus_Subgenus_species_...`. For these genera the species is read from one
  more part, and the subgenus, not the genus, counts as the genus-level unit.
  Comparisons between subgenera are then grouped with "different
  genera/subgenera".

## The figure

- **Panel A:** variability along each alignment. Tick marks below each curve
  mark every variable position.
- **Panel B:** for every comparison, the marker 1 distance (x-axis) against the
  marker 2 distance (y-axis). The dashed line marks equal divergence at both
  markers.
- **Panel C:** a close-up of the low-divergence corner of panel B (distances up
  to 0.05 by default). Selected points can be labelled (see
  [Settings file](#settings-file)).

## Requirements

- Python 3 with NumPy and Matplotlib. Tested with Python 3.14.3, NumPy 2.4.6
  and Matplotlib 3.10.9.
- The figure uses the Arial font. On a computer without Arial, Matplotlib uses
  another font and prints a warning; the content is unchanged.

Install the libraries, if needed:

    python -m pip install -r requirements.txt

## How to run it

    python marker_divergence.py marker1.fas marker2.fas --names 18S 28S

The first file goes on the x-axis and the second on the y-axis. The figure is
saved in a folder called `output` in the current folder.

Options:

| Option | Meaning | Default |
|---|---|---|
| `--names NAME1 NAME2` | short marker names, used in the axis labels and the summary | `Marker 1` `Marker 2` |
| `--titles TITLE1 TITLE2` | titles of the two panels in A | the names |
| `--species-parts N` | header parts that make up the species name | 2 |
| `--subgenus-genera GENUS ...` | genera whose headers include a subgenus | none |
| `--window N` | moving-average window for panel A, in positions | 15 |
| `--position-max N` | end of the position axis in panel A | alignment length, rounded up to 100 |
| `--scatter-max X` | axis end of panel B | largest distance, rounded up, plus 0.01 |
| `--zoom X` | axis end of panel C | 0.05 |
| `--out FOLDER` | folder for the figure files | `output` |
| `--out-name NAME` | base name of the figure files | `<name1>_vs_<name2>_divergence` |
| `--settings FILE` | JSON file with any of the settings above, and labels | none |

## Settings file

Any option above, except `--out` and `--settings`, can also be given in a JSON
file. Write the name with underscores instead of dashes. Options given on the
command line override the file. The file is also the only place to define
labels for panel C:

```json
{
  "names": ["18S", "28S"],
  "subgenus_genera": ["Plagiorhynchus"],
  "labels": [
    {"species": ["Polymorphus_minutus", "Polymorphus_obtusus"], "text_at": [0.025, 0.005]},
    {"species": ["Southwellina_hispida"], "text_at": [0.025, 0.0013], "line_from": [0, 0.5]}
  ]
}
```

- `species`: two species for a between-species point, or one for a
  within-species point. They are written as the script names them, with
  underscores, and in either order.
- `text_at`: where the label text goes, in panel C's axis units.
- `line_from`: optional. The point on the text box where the connecting line
  starts: `[0, 0]` is its lower left corner and `[1, 1]` its upper right. The
  default is its centre, `[0.5, 0.5]`.

An unknown setting stops the run with a message. A label that matches no
comparison is reported and skipped.

## Output

- `<name>.png`: 600 dpi.
- `<name>.pdf` and `<name>.svg`: vector versions, with editable text.
- On screen, a summary:
  - the number of species in both alignments;
  - how many comparisons are more divergent at marker 2, at marker 1, or equal;
  - for each level, the number of comparisons and the median distances at both
    markers;
  - the numbers of variable and parsimony-informative sites.

## Reproducing the article figure

The alignments used in the article are not included in this repository. With
the two curated alignments, the article figure is made by:

    python marker_divergence.py 18Salignment_curated.fas 28Salignment_curated.fas --settings examples/article_settings.json

[`examples/article_settings.json`](examples/article_settings.json) holds the
article's marker names, panel titles, axis limits, subgenus rule and labels.

Run this way, the script gives exactly the same distances and site values as
the scripts originally used for the article, in the same order. Its PNG is
identical, pixel for pixel, to the published figure. Those original scripts are
kept unchanged in [`original_as_run/`](original_as_run/) for the record. They
will not run as they are, because they point to folders on the authors'
computer.

The distance matrices in the article's supplementary tables were calculated
separately, not with this script, so a few values there may differ slightly
from the script's.

## Citation

If you use this script, please cite:

> [Article citation to be added on publication.]

## Licence

MIT, see [LICENSE](LICENSE).
