# Cancer gene-expression module similarity across cancer types

**Research portfolio / public code reconstruction.** Project under the guidance of Dr. Lujia Chen at the University of Pittsburgh. This repository illustrates a method used to investigate gene-expression module (GEM) similarities between colorectal and gastric cancers and between colorectal and liver cancers. **The original project inputs, collaborator-derived gene-module tables, full notebooks, and actual output results are not released here.**

## Repository structure
```text
Cancer-GEM-Comparison/
  README.md
  scripts/compare_gems.R
  examples/expression_SYNTHETIC.csv
  examples/modules_SYNTHETIC.csv
  docs/PRIVATE_VS_PUBLIC.md
  .gitignore
```

## Biological question
Do gene-expression modules derived from different tumor contexts exhibit similar transcriptomic patterns in corresponding cell-type categories?

## Comparison design
The supervised exploration considered **colorectal versus gastric** and **colorectal versus liver** cancer GEMs. Comparisons were organized by B-cell, myeloid, stromal, and T/NK-cell categories, using TCGA expression as a common scoring reference rather than treating TCGA samples as original single cells. These analyses explored similarity, not proof of a shared mechanism or independent biological validation.

## Work performed in the supervised project
The exploratory project compared GEMs across four cellular compartments: **B cells, myeloid cells, stromal cells and T/NK cells**. Work included using TCGA COAD/READ expression as a reference; Ensembl-to-gene-symbol harmonization; scoring GEMs with **gene set variation analysis (GSVA)**; comparing cross-cancer module score profiles using **cosine similarity and Pearson correlation**; examining overlapping genes; and visualizing comparisons as heatmaps. Exploratory notebooks also cover pathway enrichment. Some source notebooks and runs are exploratory and may contain incomplete intermediate results.

## What is publicly runnable
`scripts/compare_gems.R` accepts two **local, approved** CSV inputs: a prepared gene-by-sample normalized expression matrix and a long-format GEM membership table. It creates cross-cancer similarity tables and heatmaps separated by comparison and cell type. This is a **clean reconstruction**, not the original exact R Markdown pipeline or a demonstration that the historical results were independently reproduced.

```text
Cancer-GEM-Comparison/
  README.md
  scripts/compare_gems.R
  examples/expression_SYNTHETIC.csv
  examples/modules_SYNTHETIC.csv
  docs/PRIVATE_VS_PUBLIC.md
  .gitignore
```

## Input formats
1. `expression.csv` — first column `gene` (gene symbol), followed by normalized expression values for each sample. Input harmonization and data approval take place privately.
2. `modules.csv` — columns `comparison,cell_type,origin,gem,gene`; one row per member gene. Each comparison/cell-type group must have two origins (for example CRC and Gastric). Minimum five intersecting genes per module.

Packages: R, Bioconductor **GSVA**, and CRAN **ggplot2**. Install GSVA with `BiocManager::install('GSVA')` if necessary. Example run on clearly **synthetic** matrices:

```bash
Rscript scripts/compare_gems.R \
  --expression examples/expression_SYNTHETIC.csv \
  --modules examples/modules_SYNTHETIC.csv \
  --outdir private_output
```

For real analysis, use authorized, prepared inputs kept outside the repository. Never upload TCGA-derived inputs or GEM sets unless their specific reuse/sharing rights have been established.

## What is and is not included
The public R script implements module scoring, cross-origin cosine similarity, Pearson correlation, shared-gene counts, and heatmaps on local inputs. It does **not** re-create the original expression downloads, Ensembl annotation steps, notebook-specific exploratory filtering, or the pathway-enrichment notebooks. The notebooks are documented as project work but are not distributed. Synthetic matrices are example input only and do not represent cancer biology.

## Interpretation and limitations
Cosine and correlation quantify similarities of **GSVA module-score profiles across the chosen reference samples**; they do not prove equivalent pathways, cell states, causal mechanisms, or clinical relevance. Shared CRC reference matrices and the different source cancers can influence apparent similarity. The example data are **random synthetic data, not biological evidence**. The repository omits the original datasets, collaborator notebooks, and historical results pending explicit permission. It cannot reproduce the original project end-to-end without authorized input data and version-matched processing.

## Disclosure and ownership
No public license is assigned pending clarification of research code ownership and intellectual-property permissions. Confirm release of even newly reconstructed code with the supervisor/institution if needed. This repository provides portfolio-level methodology and independently written demonstration code, not a release of unpublished collaborator work.
