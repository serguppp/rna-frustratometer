# Panel 16

Prepared single-model, single-chain RNA PDBs for the panel listed in panel16.json.
The manifest records author chain IDs, selected model index, sources and SHA256,
terminal approximations, removed non-RNA residues and sequence continuity.
Original downloads remain in _sources/.

Model index 0 was selected for 2K95 and 2K96 (NMR ensembles). Six newly prepared structures use
terminal O5-prime coordinates as a P proxy: 4KQY, 6UBU, 6UC7, 6UC8, 6UC9 and 2K95.
Existing 2K96 uses the same approximation.
The terminal GDP and CCC residues of 1Y27 were explicitly reduced to G and C:
the alpha phosphate PA becomes P in GDP and extra terminal phosphate groups
are removed. These prepared ends are model input conventions.

Ligands, ions and water are removed. Comparing ligand-bound and apo geometries
with these RNA-only scorers measures RNA geometry effects without direct ligand
interaction energies. Homologous positions require alignment before cross-RNA
profile comparisons. The per-structure scan preserves the provided residue IDs.

Edit INPUTS in scripts/analysis.py to choose prepared PDBs. The default example is 2KOC. Both potentials write independent score tables and their common pairs are joined in shared_scores.csv. SimRNA defaults to geometry protocol 3 and uses full-atom templates from 2KOC.
