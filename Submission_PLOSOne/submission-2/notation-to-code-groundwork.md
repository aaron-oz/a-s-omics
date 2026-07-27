# Groundwork for R2.2 (notation-to-code table) and R2.1 (model specification)

Verified 2026-07-27 by reading the canonical pipeline. Every path, object name, and line
number below was checked against the files, not recalled. This is drafting input, not
finished prose: the supplementary table itself still has to be written.

Root for all paths: `Supplements for Publication/code/spatial-modeling/supp-4-lik-comparison/`,
abbreviated `supp-4/` below. This is the pipeline that produces Figure 4, Table 1, and the
per-feature fitted fields that everything downstream consumes.

## The fitted model, as implemented

For feature $i$ at location $s$, the Poisson branch (`supp-4/02-array-run-model.R`, lines
60 to 118) fits

$$y_i(s) \sim \text{Poisson}\big(N(s)\,\lambda_i(s)\big), \qquad
  \log \lambda_i(s) = \alpha_i + F_i(s)$$

with $N(s)$ the observed nUMI entering as a known offset (`E = total.count`), $\alpha_i$ an
intercept, and $F_i(s)$ a zero-mean Matérn Gaussian process represented by the SPDE
approach.

## Symbol to code mapping

| Symbol | Meaning | Code object | Where |
|---|---|---|---|
| $y_i(s)$ | observed count of feature $i$ at bin $s$ | `feat.count` | `run.dat`, built in `supp-4/00-launch-experiments.R` |
| $N(s)$ | observed total counts (nUMI) at bin $s$, used as a known offset | `total.count`, passed as `E =` | `02-array-run-model.R:77` |
| $\lambda_i(s)$ | modeled intensity per unit exposure | `lambda <- exp(feat.count.int + feat.count.field)` | `02-array-run-model.R:105` |
| $N(s)\lambda_i(s)$ | posterior predictive mean count | `expect <- lambda * total.count` | `02-array-run-model.R:106` |
| $\alpha_i$ | feature-specific intercept | `feat.count.int(1)` | `02-array-run-model.R:70` |
| $F_i(s)$ | spatial Gaussian process for feature $i$ | `feat.count.field(cbind(x, y), model = matern.feat.count)` | `02-array-run-model.R:71` |
| $\Sigma(\theta_i)$ | Matérn covariance of $F_i$, via SPDE | `matern.feat.count <- inla.spde2.pcmatern(...)` | `02-array-run-model.R:64` |
| $\theta_i = (\rho_i, \sigma_i)$ | range and marginal SD of the Matérn field | `matern.pri.feat.count <- c(500, .95, .1, .05)` | `02-array-run-model.R:63` |
| mesh | SPDE triangulation of the FOV | `mesh <- fm_mesh_2d(...)` | `00-launch-experiments.R:45-54` |
| $\theta_N$, $\Sigma(\theta_N)$ | latent total-count model, **not used in the reported analyses** | `matern.pri.total`, `matern.remain` | commented out at `02-array-run-model.R:31,35-38`; live version in `supp-3-.../02-array-run-model.R:30` |

Note for R2.3: the strongest possible evidence that the final analyses used the observed
offset rather than a latent $N(s)$ is that the latent-total code is **commented out** in the
canonical script. The joint latent model lives in the `supp-3-raw-v-norm-joint-feat-total-model/`
pipeline, which is what Supplement 3 reports.

## Priors, exactly as set

Penalized-complexity priors via `inla.spde2.pcmatern`, with `alpha = 2` (Matérn smoothness
$\nu = 1$ in two dimensions) and `constr = TRUE` (sum-to-zero constraint for identifiability
against the intercept):

- `prior.range = c(500, 0.95)`, that is $P(\rho_i < 500) = 0.95$, with range in raw chip
  coordinate units
- `prior.sigma = c(0.1, 0.05)`, that is $P(\sigma_i > 0.1) = 0.05$

**Watch out:** the file-level default at `02-array-run-model.R:32` is `c(500, .95, 1, .05)`,
with marginal SD 1 rather than 0.1. It is overridden inside the Poisson branch at line 63.
Quote the line-63 values. The `feat.present` prior `c(1000, .95, 10, .05)` at line 33 belongs
to the ZAP model's binomial component, not to the Poisson model.

Fitting uses `bru(..., bru_max_iter = 1)` and posterior summaries use `n.samp <- 1000` draws
(`02-array-run-model.R:55`).

## Mesh and field of view

- FOV polygon `sub.bound`: x from 15000 to 25000, y from 11500 to 18000, raw chip
  coordinates (`01-load-pkg-set-io-set-params.R:43`)
- bin size `agg.size <- 50` (`01-load-pkg-set-io-set-params.R:39`)
- `max.edge = diff(range(x))/(3*5)`, inner multiplier 0.15 at bin 50, `cutoff = max.edge/10`,
  `offset = c(max.edge, bound.outer/1.5)` (`00-launch-experiments.R:44-54`)
- resulting mesh is about 10,231 nodes over 26,000 observation bins

## Metrics (R2.4)

Code names differ from the manuscript's names, which the table should reconcile:

| Manuscript | Code | Notes |
|---|---|---|
| MAE | `MAE` | mean absolute error |
| RMSE | `RMSE` | root mean squared error |
| DS | `MDS` | mean Dawid-Sebastiani score |
| LS | `MLG` | mean log score |

Per-feature scores are written to `{feature}-scores.csv` and collated in
`03-array-post-process.R:49-64`. Two facts the reviewer's R2.4 asks about:

1. **The parenthetical percentages** are the mean over features of (that model's score
   divided by the Poisson score *for the same feature*), as a percentage. Built at
   `03-array-post-process.R:57-59` (`rel[2, 3:6] <- rel[2, 3:6] / rel[1, 3:6]`). They are not
   the ratio of the two averages printed in the table, which is why ZIP MAE shows 111% when
   0.29 / 0.21 is 138%. The caption's phrase "average percentage difference" is a misnomer
   for a ratio.
2. **The averages exclude nonconverged features** via `MAE < 100 & abs(MDS) < 100`
   (`03-array-post-process.R:70-71`), which the code comments say drops 4 and 3 features
   respectively. So the Table 1 caption's "all 918 modeled ligand/receptor features" is wrong
   twice over: the fitted set is 912, and the averaged set is smaller still.

## Verified counts (for the "Additional corrections" section)

| Quantity | Value | How established |
|---|---|---|
| genes cross-referenced against FANTOM5 | 974 | manuscript line 213, not independently re-derived |
| curated feature list | 918 | `nrow(to.use)` in `features.to.use.2025-04-04.Robj` |
| dropped before fitting | 3 (*a*, *Calcb*, *Cckar*) | `f.to.drop`, `01-load-pkg-set-io-set-params.R` |
| failed to converge | 3 (*Fgf1*, *Fzd7*, *Trhr*) | in the 915 list, no `-pred-poi.rds` produced |
| features with fitted fields | 912 | count of `*-pred-poi.rds` in the 2025-04-10 run |
| L-R pairs computable from 918 genes | 1493 | FANTOM5 pairs with both partners in the list |
| L-R pairs computable from 912 fitted | **1477** | the correct figure for the manuscript |
| rows in the clustered assay | 1481 | `dim(m.est.min)` = 1481 x 26000 |
| of which exact duplicates | 4 | PDGFC/PDGFD x PDGFRA/PDGFRB, `identical() == TRUE` against their clean twins |

The four duplicates entered as feature names carrying a cloud-sync conflicted-copy suffix,
alongside their correctly named originals. `1481 - 4 = 1477`, which reconciles the assay with
the fitted feature set.
