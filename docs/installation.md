# Installation

## Python

Use Python 3.9 or newer. The core package is installed with `python -m pip install -e .`. The regulatory-network plotting module requires Graphviz only for the best layered layout; without it, the code falls back to a deterministic bipartite layout.

```bash
conda create -n sccrs python=3.10 -y
conda activate sccrs
pip install -e .
# Optional Graphviz executable for dot layout:
conda install -c conda-forge graphviz
```

## R / Monocle3

```r
install.packages(c("optparse", "igraph", "mgcv", "dplyr", "tidyr", "ggplot2", "patchwork", "pheatmap", "Cairo"))
install.packages("BiocManager")
BiocManager::install(c("monocle3", "UCell", "SingleCellExperiment"))
install.packages("Seurat")
```

Use an R version compatible with your installed Monocle3. The trajectory script must be UTF-8 without BOM; this repository version is saved accordingly.
