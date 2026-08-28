# Revision tutorial

`tutorial.qmd` is the canonical source for the revision tutorial and renders to
`index.html`. Do not use the root-level tutorial files as the revision target.

## Reproduce the rendered tutorial

From this directory, first install the packages named in `tutorial.qmd` and the
additional packages required by `R/study_level_heterogeneity.R` (`dplyr`,
`metafor`, `brms`, `posterior`, `cmdstanr`, `drmTMB`, and `here`). `brms` requires a working Stan
toolchain. The `blsmeta` package is an optional direct Bayesian study-scale
sensitivity route, not a dependency of the rendered primary results. Its
official installation route is `remotes::install_github("donaldRwilliams/blsmeta")`;
it also requires `coda`, `rjags`, and a compatible JAGS 4.x installation.
Current `rjags` releases do not load against Homebrew JAGS 5, so record the
package and JAGS versions before attempting that sensitivity fit.

The fitted-object cache is intentionally ignored by Git. Regenerate it locally:

```sh
Rscript R/regenerate_tutorial_artifacts.R
Rscript R/study_level_heterogeneity.R --fit-bayesian
quarto render tutorial.qmd
```

The first command recreates the baseline tutorial objects. The second reads
`data/ponisio2014dataset.csv`, recalculates lnCVR with `correct = TRUE`, checks
the 318/36 full and 232/30 matched-subset counts, and writes
`Rdata/study_level_artifacts.rds`. Rendering reads these local artifacts and
does not download, fabricate, or refit them implicitly.

The Bayesian study-level script uses the explicitly recovered primary
specification that generated the historical study-SD result: four CmdStan
chains, 2,000 iterations per chain (1,000 warmup), `adapt_delta = 0.99`,
diagonal category-specific study-level covariance, `normal(0, 1)` location
coefficients, `normal(-1, 1)` residual-scale coefficients, and
`exponential(2)` study-level SDs. It also fits a prespecified `exponential(1)`
study-SD-prior sensitivity while holding data, likelihood, sampling-error
treatment, and all other priors fixed. `--mc-audit` reruns five prespecified
seeds under the recovered primary specification to quantify Monte Carlo
variation. The script exponentiates `brms` `b_sigma_*` draws before making SD
ratios or prediction intervals and uses response-by-fertiliser sampling-variance
medians for the displayed prediction scenarios. Inspect the diagnostics table
in the rendered tutorial before interpreting its estimates.

To run the optional direct Bayesian study-scale sensitivity after installing a
compatible `blsmeta`/JAGS toolchain, use:

```sh
Rscript R/study_level_heterogeneity.R --fit-blsmeta
```

This writes `Rdata/blsmeta_sensitivity.rds` and reports ratios only if the
post-warmup split-R-hat and effective-sample-size criteria pass. It is not read
by the primary rendered tutorial.
