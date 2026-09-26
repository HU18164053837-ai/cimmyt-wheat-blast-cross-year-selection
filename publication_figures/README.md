# Publication figure package

This directory contains the journal-ready figure set for the manuscript **Historical wheat blast nurseries identify resistance candidates for prospective field validation**.

## Contents

- `rebuild_phytoparasitica_figures.R`: R-only figure rebuild script.
- `outputs/`: SVG, PDF, 600-dpi TIFF, and PNG versions of five main figures and two supplementary figures.
- `outputs/Source_Data/`: source-data CSV files for every quantitative figure.
- `outputs/Figure_QA_Notes.txt`: scope and visual-quality checks.
- `outputs/R_sessionInfo.txt`: R runtime information from the verified export.

## Rebuild

From the repository root, run:

```powershell
Rscript publication_figures/rebuild_phytoparasitica_figures.R
```

The script reads the tracked tables under `history/` and overwrites the files under `publication_figures/outputs/`. It requires R 4.4.1 or later and the `svglite` package.

The plotting script does not modify analysis data, screening rules, candidate membership, or uncertainty estimates.
