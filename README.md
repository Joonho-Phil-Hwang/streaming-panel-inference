# Replication package

**Fixed-Effects Estimation and Cluster-Robust Inference for Streaming Panel Data**  
Joonho Hwang and Seojeong Lee  
Prepared for the October 15, 2026 manuscript revision (`main8`).

This package reproduces the baseline simulations, the few-treated-cluster
experiment, the storage calculation, and the Community Health Center (CHC)
application. It produces Tables 1–3 and Figures 1–3. All input data and the
validated production results are included; running the analysis downloads
nothing. The manuscript itself is maintained separately.

## Start here

Install the R packages below once. Then, from this folder, run:

```sh
Rscript replicate.R
```

This reruns the simulations with 2,000 replications per design, reconstructs
the application from the original data, runs the numerical checks, and
overwrites the generated files in `results/`. No source or input data are
changed. The run stops if a numerical check or manuscript benchmark fails.

To regenerate the tables and figures from the included saved results:

```sh
Rscript replicate.R --report-only
```

This second command checks the reporting pipeline, but does not rerun estimation
or Monte Carlo draws. To check installation with eight replications instead:

```sh
Rscript replicate.R --quick --stage=simulation
```

The short run writes to `results_quick/`; its numbers are not the paper's
results. Full production results remain in `results/`.

The entry point locates its own folder. It can be called by absolute path from
another working directory, including a path containing spaces. No path in the
code needs editing.

## Software

The package was tested on macOS with Apple Silicon and R 4.5.2. Full estimation
requires a working C++ compiler for Rcpp; the numerical kernel is compiled
automatically when the run starts. Python, Stata, MATLAB, and LaTeX are not
needed to run the analysis. The generated table fragments are intended for
inclusion in the manuscript's LaTeX preamble.

Install dependencies in R:

```r
install.packages(c(
  "Rcpp", "RcppArmadillo", "plm", "clubSandwich", "ggplot2", "patchwork",
  "data.table", "haven", "fixest", "etwfe", "digest", "tidyselect"
))
```

The versions used for the verified full run are recorded in
`results/software_versions.csv` and `results/sessionInfo.txt`. Report-only
execution needs just `ggplot2`, `patchwork`, and `data.table`.

The default is four Monte Carlo workers on macOS/Linux and one on Windows.
Use `--workers=1` for serial execution. Windows parallel execution is disabled;
the serial Windows path has not been tested. BLAS and OpenMP threads are set
to one to avoid nested parallelism. For libraries that read thread settings
before R starts, the equivalent macOS/Linux command is:

```sh
VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript replicate.R
```

See `VERIFICATION.txt` for the measured run time and the checks performed on
this delivery. Runtime depends on the machine, compiler, and R package versions.

## Files

There are four code files:

| File | Purpose |
| --- | --- |
| `replicate.R` | Options, dependencies, execution order, and manuscript benchmark checks |
| `code/simulation.R` | Both Monte Carlo experiments, numerical validation, storage counts, Tables 1–2, Figures 1–2 |
| `code/application.R` | Data preparation, CHC estimation, independent checks, Table 3, Figure 3 |
| `code/numerics.cpp` | Shared deterministic numerical routines called from R |

Other folders contain inputs and outputs:

| Folder | Contents |
| --- | --- |
| `data/` | Unchanged CHC source data, original documentation, and source/integrity notes |
| `reference/` | Small frozen summary tables used only to check the newly computed results |
| `results/simulation/` | Baseline draws, summaries, and validation results |
| `results/few_treated/` | Few-treated draws, summaries, and validation results |
| `results/application/` | CHC estimates, covariance matrices, sample counts, and checks |
| `results/storage/` | Storage counts underlying Figure 1 |
| `results/tables/`, `results/figures/` | The three manuscript tables and three figures |

`FILE_MANIFEST_SHA256.txt` records the delivered files' checksums. Generated
files may change after a rerun, so the manifest describes the delivery rather
than every subsequent state of the folder. No local R library, build cache,
development script, previous manuscript, or duplicate input archive is included.

## Outputs and manuscript correspondence

| Manuscript item | Generated file | Underlying numerical output |
| --- | --- | --- |
| Table 1: baseline approximation and coverage | `results/tables/table_1_simulation.tex` | `results/simulation/accuracy_summary.csv`, `coverage_summary.csv` |
| Table 2: few treated clusters | `results/tables/table_2_few_treated.tex` | `results/few_treated/summary.csv` |
| Table 3: CHC estimates | `results/tables/table_3_chc.tex` | `results/application/manuscript_effect_table.csv` |
| Figure 1: symmetric storage count | `results/figures/figure_1_storage.pdf` | `results/storage/storage_counts.csv` |
| Figure 2: baseline covariance accuracy | `results/figures/figure_2_accuracy.pdf` | `results/simulation/accuracy_summary.csv` |
| Figure 3: CHC event estimates | `results/figures/figure_3_chc.pdf` | `results/application/event_plot_data.csv` |

Figure numbers in these filenames follow their order in the current manuscript.
The manuscript's older figure asset filenames use a different order.

The table fragments use the manuscript's commands and packages, including
`\b`, `\citet`, `booktabs`, `array`, `float`, and `setspace`. The package
regenerates them but does not compile or modify the manuscript.

## Designs and reproducibility

**Baseline.** Section 5 uses three regressors and unit and time effects, with
2,000 replications at each of N = 100, 200, 400, 800, 1,600, and 3,200.
Each path starts with 50 units and five observed periods per unit. The three
endpoints add units, restore the held-out observations, and add period 7.
The same simulated panel and arrival priorities are used across the N grid.
The coefficient is shared by all covariance estimators at an endpoint.

The master seed is 20260903; replication r uses seed 20260903 + 100000 + r
for the panel. Separate fixed offsets set the unit permutation, held-out period,
observation priority, and new-period priority. The R RNG is Mersenne-Twister
with Inversion and Rejection. Each replication sets its seeds explicitly, so
the draws do not depend on the worker count. Wald coverage uses `qnorm(.975)`.

**Few treated clusters.** Section 5.4 keeps N = 100 and six periods, adding
treated × post to the baseline DGP. Post starts in period 4 and the true
treatment effect is zero. G = 4, 8, or 16 units are treated. The initial
20-unit batch contains two treated and 18 untreated units; the remaining
80 units arrive in random order. There are 2,000 replications per G, with
seed 4999 + r and the same explicit R RNG. Intervals use critical value 1.96.
CR3 results are retained as auxiliary output; Table 2 reports the four methods
discussed in the text. All methods use identical data and arrivals within a draw.

Panel A of Table 2 divides the mean estimated variance by the sample variance
of the 2,000 treatment-effect estimates. The benchmark is one. Panel B reports
empirical coverage in percent. Paired differences and Monte Carlo standard
errors are available in the auxiliary CSVs, although they are not displayed
in the manuscript table.

**Storage.** Figure 1 counts free entries of symmetric matrices, rather than
measuring peak process memory. With q = p(p+1)/2, exact aggregate storage is
8[q(q+1)/2 + pq + 2p² + p] bytes and approximate aggregate storage is
8[2p² + p] bytes. One GiB is 2³⁰ bytes. Per-unit summaries and scores needed
when units can return are outside these aggregate counts.

**Application.** Section 6 starts from the original Bailey–Goodman-Bacon
replication dataset. `data/SOURCE_AND_INTEGRITY.md` identifies the source and
checksum. The code reproduces the original sample restrictions, then builds
the early-cohort CHC sample: 2,541 counties and 76,230 county-year observations
from 1959–1988. The adjusted ETWFE comparison uses 2,539 counties; the pooled
comparison uses 3,043. Sample counts and exclusions are written to CSV files.

Counties enter in ascending and descending FIPS order within each year. All
counties are revisited when a calendar period opens. Covariance estimates use
county clustering, fixed population weights, working variance σ²/wᵢ I,
and no scalar degrees-of-freedom multiplier. The application uses the
manuscript's support convention and eigenvalue cutoff 10⁻¹². The output retains
full-sample CR0/CR2 references and the corrected stored-score estimates.

## Verification and saved results

A full run checks recursive estimates against direct batch calculations,
compares R and C++ implementations, and checks covariance calculations against
`plm`, `clubSandwich`, and `fixest`/`etwfe` where appropriate. It also checks
the data hash, sample construction, arrival-order invariance of exact
estimates, score identities, and correction support conditions.

The three RDS files under `results/` retain the raw Monte Carlo draws and the
application estimates needed for reporting. `reference/paper_results.rds`
contains only compact manuscript benchmarks. It is read after estimation
finishes and never supplies data or coefficients to the estimators.
`results/replication_checks.csv` records the 11 manuscript benchmark comparisons;
the analysis subfolders contain the more detailed numerical checks.

The baseline simulation writes resumable ten-replication chunks while running.
On interruption, rerun the same command: compatible chunks are reused after
the numerical validation stage. Settings and source fingerprints must match;
otherwise choose a new output folder. Completed runs remove their checkpoints.
The shorter few-treated experiment and CHC application restart from the
beginning after interruption.

Additional options, which can be combined where applicable:

```sh
Rscript replicate.R --stage=simulation --workers=1
Rscript replicate.R --stage=application
Rscript replicate.R --reps=100 --output=trial_results
Rscript replicate.R --help
```

`--output` changes the location of all generated files. With `--report-only`,
the chosen output folder must already contain the three saved RDS files
for the selected stages. Production benchmark checks for simulations apply
only to 2,000-replication results; application benchmarks always apply.

## Attribution and organization

The CHC data and original documentation are credited to Martha J. Bailey and
Andrew Goodman-Bacon (2015). They retain their original rights and terms;
this package does not assign a new license to third-party materials. The R
packages are separately licensed by their maintainers.

The organization follows useful practices in two author-maintained JBES
replication packages: a single main entry point, explicit dependencies and
data information, and a clear map from outputs to manuscript items. Examples:
[Jun and Lee](https://github.com/sokbae/replication-JunLee-JBES) and
[Angelini, Fanelli, and Neri](https://github.com/lcnnri/proxy_svar_stability_restrictions).
These were organizational references; their analysis code was not copied.
This is not a claim that JBES mandates this particular directory layout.
