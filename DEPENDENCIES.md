# Software and English file map

## Runtime requirements

| Component | Tested version | Required packages |
|---|---:|---|
| PowerShell | 5.1 or later | none |
| Python | 3.11 or later | pandas, openpyxl |
| R | 4.4.1 | base R; svglite and ragg for preferred vector and 600-dpi TIFF export |

All random-number analyses set an explicit seed. `RUN_ALL.ps1` accepts explicit paths to the Python and R executables and records the executed script hashes.

## English map of numbered scripts

| Step | Script function |
|---:|---|
| 04 | Historical field audit, harmonization, duplicate control, and cross-year nursery mapping |
| 05 | Environmental pressure drift, phenotypic spread, and recurring-GID robustness |
| 06 | Evidence-duration stratification and pedigree-source analysis |
| 07 | Nonredundant portfolio construction and historical validation-environment profile |
| 08 | LOCAL CHECK-anchored candidate comparisons |
| 09 | Rolling past-only temporal evaluation |
| 10 | Background-adjusted candidate effects and observed failure boundaries |
| 11 | Threshold multiverse, full-funnel portfolio stability, and same-size simple baselines |
| 12 | Clustered temporal uncertainty and transition-matched permutation control |
| 13 | Research workflow figure |
| 14 | Editor-requested method-transparency and submission-consistency checks |

## Analysis-unit glossary

| Unit | Locked count | Meaning |
|---|---:|---|
| Location-year-sowing environment | 38 | Broader site, environmental year, and sowing-period combination used for equal-nursery environmental summaries |
| Nursery-environment stratum | 148 | Nursery-specific stratum nested within a location-year-sowing environment |
| Candidate by nursery-environment event | 275 | One focal candidate summarized within one observed nursery-environment stratum after within-cell median aggregation |
| Cross-year-eligible GID | 104 | GID observed in at least two years, two nurseries, and four phenotype records |
| Selected GID-transition event | 462 | A training-selected GID that was reobserved in the subsequent year with sufficient evaluation records |
