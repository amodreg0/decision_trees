# Somatic SV Filtering and Validation Pipeline

This repository contains scripts for filtering somatic structural variants (SVs) called by [Severus](https://github.com/KolmogorovLab/Severus) and evaluating their accuracy against ground-truth call sets using [Truvari](https://github.com/ACEnglish/truvari). Filtering thresholds are derived from decision trees trained on labelled SV data across eight cancer cell lines.

---

## Repository structure

```
.
├── annotate_nd.py                   # Annotates Severus VCFs with ND/ND2 scores
├── annotate_lcr.py                  # Annotates Severus VCFs with LCR (low-complexity region) overlap
├── truvari_filter.sh                # SLURM array job: filters VCFs and runs Truvari benchmarks
├── decission-tree-per-cell-line.R   # Fits one decision tree per cell line + pooled ALL model
├── LOCO-cross-validation.R          # Leave-one-cell-line-out (LOCO) cross-validation
└── README.md
```

---

## Overview

### 1. VCF annotation

Before any filtering or modelling, Severus VCF files are annotated with two additional INFO fields:

**`annotate_nd.py`**
Computes and appends `ND` and `ND2` (neighbourhood depth) scores to each variant record. These scores capture the local read-depth context around each breakpoint and are among the most informative features in the downstream decision trees.

**`annotate_lcr.py`**
Flags variants whose breakpoints overlap low-complexity regions (LCRs) by appending an `in_bed` INFO/FORMAT field. LCR overlap is used as a feature in the decision trees; it can optionally also be used as a hard filter.

Both scripts write bgzipped, tabix-indexed VCFs and are designed to be run on the T2T-CHM13 (hs1) reference.

---

### 2. VCF filtering and Truvari benchmarking (`run_truvari_slurm_DT.sh`)

A SLURM array job (one task per cell line, 8 total) that:

1. **Filters** each Severus VCF using three hard thresholds derived from the pruned decision tree trained on all cell lines combined:

   | Feature | Threshold | Behaviour when missing |
   |---|---|---|
   | `SOMATICSCORE` (SS) | ≥ 50 | Missing → **keep** |
   | `ND2` | ≥ 2.9 | Missing → **keep** |
   | `PR_SR_ratio` = PR_alt / (PR_alt + SR_alt) | ≥ 0.44 | Missing → **remove** |

   Filters are applied identically regardless of SV type (BND, INS, DEL, DUP, INV). Ground-truth VCFs are passed through without filtering.

2. **Runs Truvari bench** for four comparison pairs per cell line:
   - T2T calls vs. hg38 calls (lifted)
   - Ground truth vs. hg38 calls
   - Ground truth vs. T2T calls
   - Ground truth vs. the intersection of T2T and hg38 calls (ALL3)

**Dependencies:** `conda` environment `truvari` with `truvari`, `bgzip`, `tabix`, `python3`.

---

### 3. Decision tree per cell line (`decission-tree-per-cell-line.R`)

Fits an `rpart` classification tree for each cell line independently, and a pooled tree on all cell lines combined (`ALL`), to identify which SV features best discriminate true somatic calls (`is_correct = TRUE`) from false positives.

**Input:** per-cell-line TSV files (`{cell_line}_SV_variants_table_annotated.tsv`) and a pooled file (`ALL_cell_lines_SV_variants_table_annotated.tsv`), each with columns:

| Column | Type | Description |
|---|---|---|
| `SV_type` | factor | BND, INS, DEL, DUP, INV |
| `SS` | numeric | Severus somatic score |
| `SR` | numeric | Split-read alt count |
| `PR` | numeric | Paired-read alt count |
| `PR/SR` | numeric | PR alt / SR alt ratio (renamed `PR_SR_ratio`) |
| `ND` | numeric | Neighbourhood depth score |
| `ND2` | numeric | Neighbourhood depth score (variant 2) |
| `in_bed` | logical | Breakpoint overlaps an LCR region |
| `is_correct` | logical | Ground-truth label |

**Key design choices:**
- Equal class priors (`prior = c(0.5, 0.5)`) to counteract the heavy class imbalance (~87 % FALSE).
- Two pruning levels are saved for every tree: `pruned_1se` (1-SE rule, simpler/more conservative) and `pruned_min` (minimum CV error).
- Performance metrics (Precision, Recall, F1) are computed on the **training set** — these are optimistic estimates; use LOCO results for unbiased estimates.

**Outputs** (written to `Trees/`):
- `Trees/Plots/{cell_line}_{pruning}_tree.png` — tree diagram
- `Trees/Rules/{cell_line}_{pruning}_rules.txt` — human-readable decision rules, variable importance, CP table
- `Trees/Plots/_variable_importance_across_cell_lines.png` — feature importance comparison across all cell lines
- `Trees/Rules/_top_splits_summary.tsv` — root split feature and threshold per cell line
- `Trees/Rules/_per_cell_line_summary.tsv` — Precision / Recall / F1 per cell line + mean ± SD

---

### 4. LOCO cross-validation (`LOCO-cross-validation.R`)

Performs **Leave-One-Cell-line-Out (LOCO)** cross-validation to obtain honest, generalisation-oriented performance estimates: each of the 8 cell lines is held out in turn, a tree is trained on the remaining 7, and evaluated on the held-out cell line.

This answers two questions:
- **Performance:** what Precision / Recall / F1 can we expect when applying a tree-derived filter to a new, unseen cell line?
- **Stability:** does the tree consistently select the same root split feature and threshold across all 8 folds? Consistent root splits across folds support the use of a single universal filtering rule.

**Outputs** (written to `LOCO/`):
- `LOCO/Plots/LOCO_holdout-{cell_line}_{pruning}_tree.png` — tree trained on 7 cell lines
- `LOCO/Rules/LOCO_holdout-{cell_line}_{pruning}_rules.txt` — rules + held-out cell line performance
- `LOCO/Plots/_loco_performance.png` — Precision / Recall / F1 per fold with mean dashed lines
- `LOCO/Plots/_loco_root_splits.png` — root split feature per fold (stability check)
- `LOCO/Rules/_loco_stability_summary.tsv` — per-fold metrics + mean ± SD

---

## Dependencies

### R packages
```r
install.packages(c("data.table", "rpart", "rpart.plot", "dplyr", "tidyr", "ggplot2"))
```

### Python / conda
```bash
conda activate truvari
# requires: truvari, bgzip, tabix, python >= 3.8
```

---

## Cell lines

All scripts are parameterised for the following eight cancer cell lines:

`H1437` · `H2009` · `HCC1395` · `HCC1937` · `HCC1954` · `Hs578T` · `HG008` · `COLO829`

---

## Reference genome

Filtering and Truvari benchmarking use the **T2T-CHM13v2 (hs1)** reference genome. hg38 Severus calls are liftover-converted to hs1 coordinates before benchmarking.
