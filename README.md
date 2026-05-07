# ModGenePlexus

**A module-based, network-ML approach for gene classification of long, noisy, heterogeneous gene lists from GWAS and transcriptomics studies**

[![Preprint](https://img.shields.io/badge/Preprint-bioRxiv-blue)](https://doi.org/10.1101/2025.08.11.669721)
[![Zenodo](https://img.shields.io/badge/Archive-Zenodo-blue)](https://doi.org/10.5281/zenodo.19857910)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-green)](LICENSE)

> **Paper:** McKim A, Mancuso CA, Krishnan A. *A module-based approach for post-omics, post-GWAS network-based gene classification.* bioRxiv (2025). https://doi.org/10.1101/2025.08.11.669721
>
> **Complete archive (code + data + figures + results):** https://doi.org/10.5281/zenodo.19857910


## What is ModGenePlexus?

Complex traits and diseases involve hundreds to thousands of genes that span multiple biological
processes, making them poor candidates for single-model network-based gene classification methods
like [GenePlexus](https://github.com/krishnanlab/PyGenePlexus). The core problem: GWAS and
transcriptomic gene lists are large, noisy (false positives and negatives from low power and biological
heterogeneity), and functionally heterogeneous (genes spread across multiple disconnected
neighborhoods in the genome-scale network).

**ModGenePlexus** addresses this with a divide-and-conquer strategy:

1. **Module discovery** — Input genes (e.g., DEGs from RNA-seq, MAGMA-prioritized GWAS genes)
   are clustered into topologically coherent network modules using DOMINO, which simultaneously
   denoises the list by dropping weakly connected genes and expanding each module via semi-supervised
   label propagation.
2. **Module-specific classification** — A supervised GenePlexus classifier (logistic regression on
   STRING network features) is trained independently for each module, producing genome-wide gene
   rankings per module.
3. **Score aggregation** — Module-level predictions are combined using the tau score to produce a
   final ranked gene list across all modules.

Benchmarked across simulated traits (combined GOBP gene sets), 1,517 transcriptomic gene lists
from CREEDS (diseases, gene perturbations, drug treatments), and 691 GWAS-derived gene lists from
GWAS Atlas, ModGenePlexus consistently and significantly outperforms GenePlexus — with the
performance advantage growing as gene sets become larger and more heterogeneous. Beyond
improved classification, ModGenePlexus reveals more granular and interpretable biological processes:
in a Type 2 Diabetes case study, it recovered 366 enriched GO Biological Process terms compared
to 52 from GenePlexus, including copper ion transport and iron homeostasis pathways that the
single-model approach missed entirely.


## Repository contents

This repository is a **snapshot of the analysis code** used to produce the results in the associated
paper, with emphasis on the study-bias holdout validation framework. For the full archive including
large result files, see Zenodo.

```
ModGenePlexus/
├── pygeneplexus/        # GenePlexus Python code used within the ModGenePlexus workflow
├── src/                 # Main ModGenePlexus analysis scripts
│   └── diabetes/        # Type 2 diabetes case-study scripts
├── figures/             # Figure assembly scripts and selected figure outputs
└── tsne/                # Network embedding visualization scripts and selected outputs
```

The complete Zenodo archive additionally contains all intermediate and final result files,
processed data, and supplementary figure outputs that are too large for GitHub.


## Reproducing the paper results

All analyses in the paper can be reproduced using the scripts in `src/`. The core evaluation pipeline follows these steps:

1. **Network and gene set processing** — Build the STRING v10 network (threshold: edge weight > 0.7;
   16,624 nodes, 400,729 edges) and compile gene set collections (GOBP for simulations; CREEDS
   and GWAS Atlas for real-world validation). See `src/` and File S1 in the paper for processing details.
2. **Module discovery** — Run DOMINO on each input gene list to generate network modules.
3. **Study-bias holdout evaluation** — Genes in the top two-thirds of PubMed mention frequency serve
   as training positives; understudied genes (bottom third) are held out as the test set. See Methods
   for full details on negative gene selection via PyGenePlexus.
4. **Model training and aggregation** — Train GenePlexus classifiers per module; aggregate with the
   tau score.
5. **Enrichment analysis** — Run GOBP enrichment (clusterProfiler) on top-ranked predictions from
   each module and from GenePlexus; compare information content and term specificity.

> **Note:** Some scripts assume the original file paths and compute environment from the study.
> The Zenodo archive contains all input data and result files needed to run the scripts without
> re-downloading or re-processing primary data.


## Dependencies

ModGenePlexus builds on:

- [PyGenePlexus](https://github.com/krishnanlab/PyGenePlexus) (v1.0.1) — network-based gene
  classification
- [DOMINO](https://github.com/Shamir-Lab/DOMINO) — active module identification
- [scikit-learn](https://scikit-learn.org/) — logistic regression with L2 regularization
- [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) (R) —
  GOBP enrichment analysis
- [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html) +
  [circlize](https://cran.r-project.org/web/packages/circlize/index.html) (R) — coefficient visualization


## Citation

If you use this code or data, please cite both the paper and the archive:

```bibtex
@article{mckim2025modgeneplexus,
  title   = {A module-based approach for post-omics, post-{GWAS} network-based gene classification},
  author  = {McKim, Alexander and Mancuso, Christopher A. and Krishnan, Arjun},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.08.11.669721}
}
```

**Zenodo archive:** McKim A, Mancuso CA, Krishnan A. ModGenePlexus (code + data + results).
Zenodo. https://doi.org/10.5281/zenodo.19857910


## License

BSD 3-Clause License. See [`LICENSE`](LICENSE).


## Contact

Questions about the code or method: open a [GitHub Issue](../../issues).
For broader questions about the Krishnan Lab's work on network-based gene classification and
data reuse, visit [thekrishnanlab.org](https://www.thekrishnanlab.org).
