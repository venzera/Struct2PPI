# Struct2PPI - Protein-Protein Interaction Network Visualizer from protein complexes

Generate interactive protein-protein interaction graphs from PDB/CIF files of Cryo-EM complexes (and not just Cryo-EM), or from **structure predictions** (AlphaFold3 / Boltz / Chai-1) scored with the Local Interaction Score framework (LIS / cLIS / iLIS).

![PPI Graph Example - 8XKS](Scheme.jpg)

## Features

- Parse PDB and mmCIF structure files
- Extract chain labels from COMPND (PDB) or entity information (CIF)
- Detect protein-protein interactions using distance-based cutoff (default 5.0 Å)
- Filter non-protein molecules (ligands, ions, lipids) automatically
- Generate interactive HTML visualizations with:
  - Draggable nodes
  - Colored edges with legend
  - Hover information (chain names, contact counts)
- Export chain descriptions and interacting residue pairs to text files

## Scripts

### `ppi_graph.py` - Standard 2D Network Graph
Interactive Plotly-based visualization with draggable nodes and colored edges.

```bash
python ppi_graph.py structure.cif
python ppi_graph.py structure.pdb --cutoff 5.0
```

### `ppi_graph_3d.py` - STRING DB Style with 3D Structures
Network graph with 3D protein structures rendered inside nodes using 3Dmol.js.

```bash
python ppi_graph_3d.py structure.cif
python ppi_graph_3d.py structure.pdb --cutoff 4.0
```

### `ppi_graph_3d_dg.py` - 3D Structures with Binding Energy (PRODIGY + FoldX)
Network graph with 3D protein structures and binding energy (ΔG) calculated using PRODIGY for each interacting chain pair. Supports FoldX integration for force-field-based scoring and mutation analysis. Accepts both PDB and mmCIF input — CIF files (e.g. AlphaFold3 output) are transparently converted to PDB via [gemmi](https://gemmi.readthedocs.io/) before running PRODIGY / FoldX.

Features:
- Calculates binding affinity (ΔG in kcal/mol) for each chain pair using PRODIGY
- Displays ΔG values as edge labels on the graph
- Legend sorted by binding strength (strongest interactions first)
- Creates separate PDB files for each interacting chain pair
- **Clickable edges** — click any edge (or legend item) to open a modal with the 3D pair complex, with the two chains colored distinctly (blue/red). Hover an edge to preview contacts + ΔG before clicking
- **Double-click a node** — opens the same modal showing just that single chain's 3D structure
- Chain labels come from CIF entity descriptions (or PDB `COMPND`); when unavailable (e.g. raw AlphaFold3 CIF) the chain ID itself is used as the label instead of erroring
- Outputs binding strength ranking file
- FoldX RepairPDB + AnalyseComplex scoring (`--fx_score`)
- FoldX BuildModel mutation ΔΔG analysis (`--fx_mut`)
- One-vs-all PRODIGY mode: each chain vs all others combined (`--one_vs_all`)

```bash
# PRODIGY binding energy (PDB or CIF)
python ppi_graph_3d_dg.py structure.pdb
python ppi_graph_3d_dg.py structure.cif
python ppi_graph_3d_dg.py structure.pdb --cutoff 4.0
python ppi_graph_3d_dg.py structure.pdb --skip-prodigy  # Skip PRODIGY (for testing)

# FoldX scoring (RepairPDB + AnalyseComplex)
python ppi_graph_3d_dg.py structure.pdb --fx_score --fx_path /path/to/foldx

# FoldX mutation ddG (BuildModel, skips repair)
python ppi_graph_3d_dg.py structure.pdb --fx_mut individual_list.txt --fx_path /path/to/foldx

# One-vs-all: each chain vs all others combined (PRODIGY)
python ppi_graph_3d_dg.py structure.pdb --one_vs_all
```

#### Flags

| Flag | Description |
|------|-------------|
| `--fx_score` | Run FoldX RepairPDB on the input structure, then AnalyseComplex to score interaction energies for each chain pair |
| `--fx_path` | Path to the FoldX executable (required with `--fx_score` or `--fx_mut`) |
| `--fx_mut` | Path to `individual_list.txt` with mutations for FoldX BuildModel. Skips repair (copies input with `_Repair.pdb` suffix), runs BuildModel, then AnalyseComplex on both WT and mutant structures to calculate binding ΔΔG per chain pair |
| `--one_vs_all` | For each chain, combine all other chains as a single partner and calculate binding energy via PRODIGY (`--selection A B,C,D,...`) |

#### Why one-vs-all mode?

According to Kastritis & Bonvin (*J Mol Biol* 2014), PRODIGY's scoring function explicitly depends on the non-interacting surface (NIS) — polar and charged residues on the solvent-exposed surface outside the interface contribute to ΔG through hydration shell stability and long-range electrostatics. I assumed that in a multi-chain complex, the NIS composition changes when additional partners are bound (some surfaces get buried), so simply summing pairwise ΔG values would use wrong NIS contexts for each pair. The one-vs-all mode lets PRODIGY see the actual NIS of the full assembly, which should give a more realistic estimate.

### `ppi_graph_predicted.py` - Predicted complexes (LIVIA / AFM-LIS integration)

Build a PPI graph from a **predicted** complex (AlphaFold3, Boltz, Chai-1) and score every chain pair with the **Local Interaction Score** framework instead of PRODIGY/FoldX binding energy. This is the structure-prediction counterpart to `ppi_graph_3d_dg.py`, integrating the LIS / cLIS / iLIS metrics from [LIVIA / AFM-LIS](https://github.com/flyark/AFM-LIS) (Kim & Perrimon, 2026).

The HTML is labelled **PREDICTED MODE** and shows the platform + model name. **Clicking an edge** opens a panel that, for that chain–chain interaction, shows:

- the **3D pair complex** (two chains colored)
- the inter-chain **PAE interaction plot**
- the **contact-filtered LIA (cLIA) map** (PAE ≤ 12 Å **and** Cβ–Cβ ≤ 8 Å)
- **LIS / cLIS / iLIS**, LIA / cLIA / iLIA, ipTM, and per-chain + interface **pLDDT** scores

A pair is flagged a likely interaction when **iLIS ≥ 0.223** (Kim et al., 2026). **Double-click a node** for that single chain's 3D structure and average pLDDT. Edge thickness scales with iLIS; positive pairs are drawn solid, sub-threshold pairs dashed.

```bash
# Single prediction (folder of AF3 / Boltz / Chai output) → one graph
python ppi_graph_predicted.py /path/to/af3_output_folder/
python ppi_graph_predicted.py /path/to/boltz_folder/ --model 0
python ppi_graph_predicted.py /path/to/predictions/ --output-dir ./results

# Batch over many 2-chain predictions: pick the #1 model by iLIS per prediction,
# emit a ranking + a graph for each dimer
python ppi_graph_predicted.py /path/to/dimers_parent_folder/ --batch-dimer-predicted
```

> **Input must include the predictor's PAE / confidence file — not just the structure.**
> LIS is computed from the Predicted Aligned Error, which is **not** stored in the `.cif`/`.pdb`. Point the tool at the *whole output folder*, which must contain:
> - **AlphaFold3:** `*_model_N.cif` **+** `*_full_data_N.json` (PAE) **+** `*_summary_confidences_N.json`
> - **Boltz:** `*_model_N.cif/.pdb` **+** `confidence_*.json` / `pae_*.npz`
> - **Chai-1:** `pred.*.cif` **+** `scores.*.json` **+** `pae.*.npy/.npz`
>
> A structure file on its own (e.g. an AlphaFold Server `fold_*_model_N.cif` downloaded without its `full_data` JSON) cannot be scored and the tool exits with a message saying which companion file is missing.

#### Flags

| Flag | Description |
|------|-------------|
| `--predicted` | Predicted mode (default): one prediction → one graph; uses the best model by mean iLIS |
| `--batch-dimer-predicted` | Scan a parent folder of 2-chain predictions; for each, select the #1 model by iLIS and write `batch_dimer_predicted_ranking.txt` plus a per-dimer graph |
| `--model` | Model rank to visualize (default: best mean iLIS) |
| `--edge-min-ilis` | Only draw edges with iLIS ≥ this value (default: 0, draw all) |
| `--pae-cutoff` | PAE cutoff in Å for the confident interface (default: 12) |
| `--cb-cutoff` | Cβ–Cβ contact cutoff in Å for cLIS/cLIA (default: 8) |

#### The LIS metrics

- **LIS** — average interface confidence over all residue pairs with PAE ≤ 12 Å.
- **cLIS** — same, but restricted to pairs that are also in direct contact (Cβ–Cβ ≤ 8 Å).
- **iLIS = √(LIS × cLIS)** — the recommended single metric; the geometric mean forces iLIS → 0 when no confident contacts exist, removing the "confident but physically distant" false positives.
- **LIA / cLIA** — counts of confident / confident-and-contacting residue pairs (the *area* of the interface).

LIS scoring is provided by the bundled `lis_lib.py`, a trimmed adaptation of AFM-LIS's `lis.py` for Struct2PPI: parallel processing, zip/tar/gz archive handling, and ChimeraX `.cxc` export are removed; the core scoring matches upstream exactly.

```bash
pip install biopython networkx scipy numpy gemmi   # predicted mode needs no extra deps
```

## Output Files

For input `structure.pdb` (or `.cif`, supported by all three scripts):
- `structure_ppi_graph.html` - Interactive 2D network (ppi_graph.py)
- `structure_ppi_3d.html` - Network with 3D structure nodes (ppi_graph_3d.py)
- `structure_ppi_3d_dg.html` - Network with 3D nodes + binding energy (ppi_graph_3d_dg.py)
- `structure_chain_info.txt` - Chain ID to protein name mapping
- `structure_residue_contacts.txt` - All interacting residue pairs (with ΔG for ppi_graph_3d_dg.py)
- `structure_binding_strength.txt` - Chain pairs sorted by binding strength (ppi_graph_3d_dg.py)
- `structure_complexes/` - PDB files for each chain pair (ppi_graph_3d_dg.py)
- `structure_foldx_scores.txt` - FoldX interaction energy ranking (`--fx_score`)
- `structure_foldx_ddg.txt` - FoldX binding ΔΔG per chain pair (`--fx_mut`)
- `structure_one_vs_all.txt` - One-vs-all PRODIGY ranking (`--one_vs_all`)
- `structure_foldx/` - FoldX working directory with intermediate files (`--fx_score` or `--fx_mut`)

For predicted mode (`ppi_graph_predicted.py`, input = prediction folder named `<name>`):
- `<name>_predicted_ppi.html` - Interactive predicted PPI graph with PAE / cLIA maps + LIS scores
- `<name>_lis_scores.txt` - Per-chain-pair LIS / cLIS / iLIS / LIA / cLIA / ipTM / pLDDT, sorted by iLIS, with a POSITIVE/negative call
- `batch_dimer_predicted_ranking.txt` - Dimers ranked by best-model iLIS (`--batch-dimer-predicted`)

## Dependencies

```bash
pip install biopython networkx plotly scipy numpy prodigy-prot gemmi
```

| Package | Purpose |
|---------|---------|
| biopython | PDB/CIF parsing, structure manipulation |
| networkx | Graph building and layout algorithms |
| plotly | Interactive 2D visualization |
| scipy | Distance calculations (cdist) |
| numpy | Numerical operations |
| prodigy-prot | Binding energy calculation (ppi_graph_3d_dg.py) |
| gemmi | CIF → PDB conversion for ppi_graph_3d_dg.py (PRODIGY / FoldX require PDB input) |

For FoldX features (`--fx_score`, `--fx_mut`), [FoldX](https://foldxsuite.crg.eu/) must be installed separately (academic license required).

For 3D visualization (`ppi_graph_3d.py`, `ppi_graph_3d_dg.py`), 3Dmol.js is loaded from CDN (no installation required).

## Usage Examples

```bash
# Basic usage with CIF file
python ppi_graph.py 8xks.cif

# Custom distance cutoff
python ppi_graph.py structure.pdb --cutoff 4.0

# 3D structure nodes (STRING DB style)
python ppi_graph_3d.py 8xks.cif

# 3D structure nodes with binding energy (PDB or CIF — CIF is auto-converted via gemmi)
python ppi_graph_3d_dg.py 8xks.pdb
python ppi_graph_3d_dg.py af3_prediction.cif

# Specify output directory
python ppi_graph.py structure.cif --output-dir ./results
python ppi_graph_3d_dg.py structure.pdb --output-dir ./results

# FoldX: repair structure and score all chain pair interactions
python ppi_graph_3d_dg.py 8xks.pdb --fx_score --fx_path /path/to/foldx

# FoldX: calculate mutation effect on binding between chains
python ppi_graph_3d_dg.py 5L08.pdb --fx_mut individual_list.txt --fx_path /path/to/foldx --skip-prodigy

# One-vs-all: rank chains by total binding contribution
python ppi_graph_3d_dg.py 8xks.pdb --one_vs_all --skip-prodigy
```

## Test Data

The `test/` folder contains example outputs generated from PDB entry 8XKS (20 protein chains, 81 interactions).

### Predicted mode — CTCF / RAD21 (AlphaFold3)

`test/fold_ctcf_rad21_branchio/` is a complete AlphaFold3 server output (3 chains: CTCF + two RAD21 fragments) used to exercise predicted mode end-to-end:

```bash
python ppi_graph_predicted.py test/fold_ctcf_rad21_branchio --output-dir test/fold_ctcf_rad21_branchio
```

All three chain pairs score as positive interactions (iLIS ≥ 0.223):

| Pair | iLIS | LIS | cLIS | LIA | cLIA | ipTM |
|------|------|-----|------|-----|------|------|
| A–B | 0.766 | 0.667 | 0.879 | 18507 | 134 | 0.940 |
| B–C | 0.651 | 0.572 | 0.740 | 539 | 10 | 0.420 |
| A–C | 0.630 | 0.538 | 0.736 | 5626 | 40 | 0.830 |

`test/fold_p53_mdm2_model_2.cif` is included as the *counter-example*: a structure-only AlphaFold Server download with **no** `full_data` PAE JSON. Running predicted mode on it exits with a message stating which companion file is required — illustrating that the scores/PAE file (not just the model) is mandatory for LIS.

### 5L08 — Pairwise vs One-vs-All (Caspase-8 tDED filament)

PDB ID [5L08](https://www.rcsb.org/structure/5L08) is a Caspase-8 complex with 9 identical chains (A–I, 19 interacting pairs).


Sum of one-vs-all ΔG = **-120.80 kcal/mol** (9 chains, average -13.42 kcal/mol per chain).

The one-vs-all ranking reveals which chains are most embedded in the complex: chain B (ΔG = -18.10) has the strongest total binding contribution, while chain D (ΔG = -9.20) is most peripheral. This information is not directly extractable from pairwise scoring, where individual pair values lack the NIS context of the full assembly.

**Σ pairwise ΔG per chain vs one-vs-all ΔG:**

| Rank | Chain | Σ pairwise ΔG | # pairs | One-vs-all ΔG | Δ |
|------|-------|---------------|---------|---------------|-------|
| 1 | B | -39.7 | 6 | -18.1 | +21.6 |
| 2 | H | -33.5 | 5 | -16.6 | +16.9 |
| 3 | E | -31.3 | 5 | -14.3 | +17.0 |
| 4 | A | -26.9 | 4 | -14.7 | +12.2 |
| 5 | F | -26.0 | 4 | -13.6 | +12.4 |
| 6 | C | -24.1 | 4 | -11.6 | +12.5 |
| 7 | G | -23.3 | 4 | -10.7 | +12.6 |
| 8 | I | -20.3 | 3 | -12.0 | +8.3 |
| 9 | D | -17.3 | 3 | -9.2 | +8.1 |

The relative ranking is preserved — chains B, H, E remain the most connected, D remains the most peripheral — but the absolute values diverge substantially. Summing pairwise ΔG overestimates binding by +8 to +22 kcal/mol per chain. The discrepancy scales with the number of interaction partners: chain B (6 pairs, Δ = +21.6) vs chain D (3 pairs, Δ = +8.1). This is consistent with pairwise scoring using incorrect NIS contexts — each isolated pair calculation assumes the remaining surface is fully solvent-exposed, when in reality other partners bury parts of it.

## ΔG Evaluation
### Binding Energy Prediction
PRODIGY-based binding energy prediction enables rapid assessment of protein-protein interaction strengths directly from experimental structures. While PRODIGY excels at providing reliable relative rankings between complexes, users seeking accurate absolute binding free energy values should consider more rigorous methods such as FoldX, Rosetta, or MM/GBSA.
### Validation
Benchmarking on PDB ID 8XKS demonstrated strong correlation between PRODIGY and FoldX 5.1 ΔG values (r = 0.966). Notably, structure preprocessing with FoldX RepairPDB (repairing missing atoms, energy minimization) did not substantially improve this correlation (r = 0.964), suggesting that PRODIGY performs robustly even on unprocessed experimental structures. And there is no need to make FoldX's PDB repairing which can take a lot of time for Cryo-EM complexes.
![Binding Energy Correlation](evaluation/evaluation_8xks.png)

> Delgado J., Reche R., Cianferoni D., Orlando G., van der Kant R., Rousseau F., Schymkowitz J., Serrano L. "FoldX force field revisited, an improved version." *Bioinformatics*, Volume 41, Issue 2, btaf064 (2025). DOI: [10.1093/bioinformatics/btaf064](https://doi.org/10.1093/bioinformatics/btaf064)

## References

If you use the binding energy prediction feature (`ppi_graph_3d_dg.py`), please cite PRODIGY:

> Xue L., Rodrigues J., Kastritis P., Bonvin A.M.J.J., Vangone A. "PRODIGY: a web server for predicting the binding affinity of protein-protein complexes." *Bioinformatics* (2016). DOI: [10.1093/bioinformatics/btw514](https://doi.org/10.1093/bioinformatics/btw514)

> Vangone A. and Bonvin A.M.J.J. "Contacts-based prediction of binding affinity in protein-protein complexes." *eLife*, 4:e07454 (2015). DOI: [10.7554/eLife.07454](https://doi.org/10.7554/eLife.07454)

> Kastritis P.L., Rodrigues J.P.G.L.M., Folkers G.E., Boelens R., Bonvin A.M.J.J. "Proteins Feel More Than They See: Fine-Tuning of Binding Affinity by Properties of the Non-Interacting Surface." *Journal of Molecular Biology*, 14, 2632–2652 (2014). DOI: [10.1016/j.jmb.2014.04.017](https://doi.org/10.1016/j.jmb.2014.04.017)

If you use predicted mode (`ppi_graph_predicted.py`) and the LIS / cLIS / iLIS metrics, please cite LIVIA / AFM-LIS:

> Kim, Ah-Ram, and Norbert Perrimon. "LIVIA: a browser-based tool for assessing and visualizing predicted protein interactions." *bioRxiv* (2026): 2026-05.

> Kim A.-R., et al. "Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer." *bioRxiv* (2024). DOI: [10.1101/2024.02.19.580970](https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1)

The LIS implementation in `lis_lib.py` is adapted from the AFM-LIS reference code: [github.com/flyark/AFM-LIS](https://github.com/flyark/AFM-LIS).

CIF → PDB conversion in `ppi_graph_3d_dg.py` and `ppi_graph_predicted.py` uses GEMMI:

> Wojdyr, Marcin. "GEMMI: A library for structural biology." *Journal of Open Source Software* 7.73 (2022): 4200. DOI: [10.21105/joss.04200](https://doi.org/10.21105/joss.04200)
