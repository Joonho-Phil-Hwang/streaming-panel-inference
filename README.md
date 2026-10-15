# Streaming panel inference

Replication files for **Fixed-Effects Estimation and Cluster-Robust Inference
for Streaming Panel Data**, by Joonho Hwang and Seojeong Lee.

The package reproduces Tables 1–3 and Figures 1–3: the baseline and
few-treated-cluster simulations, storage calculations, and the Community
Health Center (CHC) application. All input data and saved results are included.

## Requirements

Tested with R 4.5.2 on macOS (Apple Silicon). Full estimation requires a C++
compiler for Rcpp. Install the R packages once:

```r
install.packages(c(
  "Rcpp", "RcppArmadillo", "plm", "clubSandwich", "ggplot2", "patchwork",
  "data.table", "haven", "fixest", "etwfe", "digest", "tidyselect"
))
```

Report-only execution needs `ggplot2`, `patchwork`, and `data.table`.
[Tested package versions](results/software_versions.csv) are included.

## Run

From the repository folder:

```sh
Rscript replicate.R                              # Full replication
Rscript replicate.R --report-only                # Tables and figures from saved results
Rscript replicate.R --quick --stage=simulation    # Eight-replication installation check
```

The full run uses 2,000 replications per simulation design and reconstructs
the CHC analysis from the original data. It took about seven minutes on the
tested machine, excluding initial compilation. Paths are detected automatically;
no code edits or data downloads are needed.

Generated files are overwritten in `results/`; quick-run output goes to
`results_quick/`. The default is four Monte Carlo workers on macOS/Linux and
one on Windows, which has not been tested. Use `--workers=1` for serial
execution, `--stage=simulation` or `--stage=application` for one analysis,
and `--help` for all options.

## Files and outputs

| Program | Purpose |
| --- | --- |
| `replicate.R` | Run all analyses and check manuscript benchmarks |
| `code/simulation.R` | Simulations and storage calculations; Tables 1–2 and Figures 1–2 |
| `code/application.R` | CHC analysis; Table 3 and Figure 3 |
| `code/numerics.cpp` | Shared numerical routines called from R |

Tables are in [results/tables/](results/tables/) and figures in
[results/figures/](results/figures/). Figure 1 shows storage requirements;
Figure 2, baseline covariance accuracy; and Figure 3, CHC event estimates.
Draws, estimates, and numerical checks are retained in the corresponding
subfolders of `results/`. LaTeX tables use the manuscript's preamble and commands.

## Data and verification

The CHC data come from the replication archive for
[Bailey and Goodman-Bacon (2015)](https://doi.org/10.1257/aer.20120070).
The unchanged data and original documentation are in `data/`; see
[source and integrity notes](data/SOURCE_AND_INTEGRITY.md) for provenance
and the data checksum. These materials retain their original rights and terms.

The full run passed 234 numerical checks and 11 manuscript benchmark comparisons.
Files in `reference/` contain only comparison benchmarks, read after estimation.
See [VERIFICATION.txt](VERIFICATION.txt) for details.
[FILE_MANIFEST_SHA256.txt](FILE_MANIFEST_SHA256.txt) records the delivered
files' checksums; regenerated outputs may differ after a rerun.
