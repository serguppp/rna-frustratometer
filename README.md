# RNA Frustratometer

Mutational frustration of RNA contacts with cgRNASP and SimRNA. The package scores
all 16 A/C/G/U combinations at each selected pair, including the native identity.
It computes F using the mutual interaction energy (`F_pair`) and the full RNA
energy (`F_total`), with sequence composition weights and a uniform control.

## Setup and run

Requires Python 3.12, uv, and the local energy engines. For cgRNASP, build the
[cgRNASP fork](https://github.com/serguppp/cgRNASP/tree/0fc74f8) using
`python build.py`. For SimRNA, obtain the
[licensed SimRNA 3.20 dynamic Linux distribution](https://iimcb.genesilico.pl/index.php/software/simrna).
Its local engine folder must contain `SimRNA`, `config_n0.dat` and `data/`.

```powershell
uv sync --locked
docker build -f Dockerfile.simrna -t rna-frustratometer-simrna:local .
uv run scripts/analysis.py
```

Edit the paths and chosen structure in `scripts/analysis.py`. The example uses
the prepared, single-model `tests/data/sim/2KOC.pdb` for both input and base
templates. Use a new output directory each time. The SimRNA container reads
the package source directly and mounts the local engine and data read-only.
Its input needs canonical RNA, the five scoring sites and the atoms defining
the glycosidic frame. Check missing atoms and alternate conformations during
input preparation.

- `pipeline.py` prepares the reference and calls the selected scorer.
- `cgrnasp.py` invokes the C++ scan, which keeps bead coordinates fixed.
- `simrna_mutation.py` places each new base template in the native glycosidic
  frame. Default protocol 3 preserves the sugar and the N1/N9 anchor. Protocol
  4 also keeps ring coordinates in same-class substitutions.
- `simrna.py` scores cached variants in one container. A separate diagnostic
  score isolates each pair and removes local backbone strain from its energy.
- `frustration.py` computes `F = (mean(E_variant) − E_native) / SD(E_variant)`
  using population moments. Undefined scores are empty. The SimRNA pair SD
  cutoff accounts for printed energy precision and is unrelated to F classes.

Each potential writes its raw `contact_variants.csv`, calculated
`frustration_scores.csv` and weight metadata. `shared_scores.csv` joins the
same mutation targets after scoring. The pair and total contact maps plot the
intersection with defined F using a common color scale within each figure.
The two engines use different criteria to choose mutation targets and their
statistical scores have model-specific units.

## Compare the homolog panel

`uv run scripts/homologs.py` reads the local panel manifest, the existing Rfam
and US-align mapping, and the saved energies from the 40-structure experiment
under `results/`. These local inputs are outside Git. Edit the paths at the
top of the script for another panel.

The comparison uses the same aligned pairs for both potentials and both members
of each homolog comparison. It saves the matched F and native pair energies in
`matched_contacts.csv`, and coverage, Spearman correlation and median absolute
F difference in `homolog_statistics.csv`. Results are stratified by all mapped
positions, Rfam/3D agreement, geometric core and conserved identities. Pooled
family composition and uniform weights are computed from the same saved
energies. PNG profiles show mean incident F and contact counts for selected
families. These are descriptive comparisons of homologs.

Run `uv run pytest`, `uv run ruff check .` and `uv run mypy` for local checks.
The SimRNA Python bridge was adapted from
[`simrna-adapter` at `31be7d6`](https://github.com/serguppp/simrna-adapter/tree/31be7d6)
under [its MIT license](LICENSE-SIMRNA-ADAPTER).
