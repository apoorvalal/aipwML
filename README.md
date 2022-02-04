# `aipwML` -- Regression adjustment, IPW, and AIPW estimation

The [`aipwML`](https://github.com/apoorvalal/aipwML) package computes
causal effects using nuisance functions estimated using linear /
logistic regression, regularised regression (fit with `glmnet`), random
forests (fit with `ranger`), and gradient boosted trees (fit with
`xgboost`). It is written to be as modular and possible so that users
can specify different choices for the outcome and propensity score
models in `mhatter` and `ehatter`.

## installation

``` r
library(remotes)
remotes::install_github("apoorvalal/aipwML")
```

This writeup demonstrates the estimation functions using the Lalonde
observational dataset where experimental controls were replaced with
control units from the PSID, and standard estimators are badly biased
for the experimental effect of ≈ $1700.

Data Prep
=========

``` r
data(lalonde.psid); df = lalonde.psid
y = 're78'; w = 'treat'
x = setdiff(colnames(df), c(y, w))

# outcome model formula
fo = re78 ~ (age + education + black + hispanic + married + nodegree +
    re74 + re75 + u74 + u75)
# pscore formula
fp = treat ~ (age + education + black + hispanic + married + nodegree +
    re74 + re75 + u74 + u75)
```

We have data
{*y*<sub>*i*</sub>, *w*<sub>*i*</sub>, *x*<sub>*i*</sub>}<sub>*i* = 1</sub><sup>*N*</sup> ∈ ℝ × {0, 1} × ℝ<sup>*k*</sup>.
Under selection on observables assumptions, we can compute the ATE by
imputing the missing potential outcome.

Regression Adjustment
=====================

$$
\\hat{\\tau}\_{\\text{reg}}^{\\text{ATE}}  = \\frac{1}{N} \\sum\_{i=1}^N (
    \\hat{\\mu}\_1 (x\_i) - \\hat{\\mu}\_0 (x\_i)
)
$$

``` r
regadjusts = c(
  ate_reg('ols',      w = w, y = y, df = df, fml = fo),
  ate_reg('lasso',    w = w, y = y, df = df, fml = fo),
  ate_reg('ridge',    w = w, y = y, df = df, fml = fo),
  ate_reg('rforest',  w = w, y = y, df = df, fml = fo),
  ate_reg('xgboost',  w = w, y = y, df = df, fml = fo)
)
regadjusts |> round(3)
```

    ## [1]  -8746 -13824 -12260  -9677 -11623

pretty bad.

Inverse Propensity Weighting (IPW)
==================================

$$
\\hat{\\tau}\_{\\text{ipw}}^{\\text{ATE}} = \\frac{1}{N} \\sum\_{i=1}^N
\\frac{y\_i (w\_i - \\hat{e}(x\_i)) }{\\hat{e}(x\_i) (1 - \\hat{e}(x\_i)) }
$$

``` r
ipws = c(
  ate_ipw('logit',   w = w, y = y, df = df, fml = fp),
  ate_ipw('lasso',   w = w, y = y, df = df, fml = fp),
  ate_ipw('ridge',   w = w, y = y, df = df, fml = fp),
  ate_ipw('rforest', w = w, y = y, df = df, fml = fp),
  ate_ipw('xgboost', w = w, y = y, df = df, fml = fp)
)

ipws |> round(3)
```

    ## [1] -10454 -15260 -17703 -19568 -19442

Still pretty bad. Now trim extreme pscores.

``` r
psr = c(0.05, 0.95)
ipws2 = c(
  ate_ipw('logit',   w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('lasso',   w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('ridge',   w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('rforest', w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('xgboost', w = w, y = y, df = df, fml = fp, psrange = psr)
)

ipws2 |> round(3)
```

    ## [1] -1356.4 -1144.8 -1623.8   361.6  2971.2

Better.

Augmented IPW
=============

$$
\\hat{\\tau}\_{\\mathrm{AIPW}}^{\\text{ATE}} =
  \\frac{1}{N} \\sum\_{i=1}^{N}
  \\left\[\\left(
    \\hat{m}\_{1}\\left(x\_{i}\\right)+\\frac{w\_{i}}{\\hat{e}\\left(x\_{i}\\right)}
    \\left(y\_{i}-\\hat{m}\_{1}\\left(x\_{i}\\right)\\right)\\right) -
    \\left(\\hat{m}\_{0}\\left(x\_{i}\\right)+\\frac{1-w\_{i}}{1-\\hat{e}\\left(x\_{i}\\right)}
    \\left(y\_{i}-\\hat{m}\_{0}\\left(x\_{i}\\right)\\right)\\right)\\right\]
$$

Need to chose how to estimate *e* and *m*, so we perform an exhaustive
search. For each choice of outcome model, I try every other fitter for
the pscore.

OLS outcome model
-----------------

``` r
ols_mean = c(
  ate_aipw(fit_me(meanfn = 'ols', pscorefn = 'logit',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ols', pscorefn = 'lasso',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ols', pscorefn = 'ridge',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ols', pscorefn = 'rforest',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ols', pscorefn = 'xgboost',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr)
)
```

Ridge outcome model
-------------------

``` r
ridge_mean = c(
  ate_aipw(fit_me(meanfn = 'ridge', pscorefn = 'logit',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ridge', pscorefn = 'lasso',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ridge', pscorefn = 'ridge',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ridge', pscorefn = 'rforest',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'ridge', pscorefn = 'xgboost',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr)
)
```

LASSO outcome model
-------------------

``` r
lasso_mean = c(
  ate_aipw(fit_me(meanfn = 'lasso', pscorefn = 'logit',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'lasso', pscorefn = 'lasso',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'lasso', pscorefn = 'ridge',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'lasso', pscorefn = 'rforest',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'lasso', pscorefn = 'xgboost',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr)
)
```

Random Forest outcome model
---------------------------

``` r
rforest_mean = c(
  ate_aipw(fit_me(meanfn = 'rforest', pscorefn = 'logit',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'rforest', pscorefn = 'lasso',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'rforest', pscorefn = 'ridge',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'rforest', pscorefn = 'rforest',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'rforest', pscorefn = 'xgboost',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr)
)
```

GBM outcome model
-----------------

``` r
xgboost_mean = c(
  ate_aipw(fit_me(meanfn = 'xgboost', pscorefn = 'logit',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'xgboost', pscorefn = 'lasso',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'xgboost', pscorefn = 'ridge',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'xgboost', pscorefn = 'rforest',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr),
  ate_aipw(fit_me(meanfn = 'xgboost', pscorefn = 'xgboost',
    mean_fml = fo, psc_fml = fp, y = y, w = w, df = df), psrange = psr)
)
```

AIPW table
----------

``` r
# stack estimates
aipw_estimates = rbind(ols_mean, lasso_mean, ridge_mean, rforest_mean, xgboost_mean)
colnames(aipw_estimates) = c('PS: logit', 'PS: lasso', 'PS: ridge',
  'PS: rforest', 'PS: xgboost')
rownames(aipw_estimates)= c('Outcome: ols', 'Outcome: lasso', 'Outcome: ridge',
  'Outcome: rforest', 'Outcome: xgboost')
aipw_estimates |> kbl() %>%
  kable_styling()
```

<table class="table" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
PS: logit
</th>
<th style="text-align:right;">
PS: lasso
</th>
<th style="text-align:right;">
PS: ridge
</th>
<th style="text-align:right;">
PS: rforest
</th>
<th style="text-align:right;">
PS: xgboost
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Outcome: ols
</td>
<td style="text-align:right;">
56.49
</td>
<td style="text-align:right;">
1.199
</td>
<td style="text-align:right;">
-519.4
</td>
<td style="text-align:right;">
439.2
</td>
<td style="text-align:right;">
3227
</td>
</tr>
<tr>
<td style="text-align:left;">
Outcome: lasso
</td>
<td style="text-align:right;">
-236.96
</td>
<td style="text-align:right;">
-74.295
</td>
<td style="text-align:right;">
-909.8
</td>
<td style="text-align:right;">
558.7
</td>
<td style="text-align:right;">
3213
</td>
</tr>
<tr>
<td style="text-align:left;">
Outcome: ridge
</td>
<td style="text-align:right;">
-245.79
</td>
<td style="text-align:right;">
-478.324
</td>
<td style="text-align:right;">
-620.9
</td>
<td style="text-align:right;">
189.8
</td>
<td style="text-align:right;">
3247
</td>
</tr>
<tr>
<td style="text-align:left;">
Outcome: rforest
</td>
<td style="text-align:right;">
-197.46
</td>
<td style="text-align:right;">
110.549
</td>
<td style="text-align:right;">
-355.1
</td>
<td style="text-align:right;">
764.2
</td>
<td style="text-align:right;">
3202
</td>
</tr>
<tr>
<td style="text-align:left;">
Outcome: xgboost
</td>
<td style="text-align:right;">
-320.45
</td>
<td style="text-align:right;">
6.479
</td>
<td style="text-align:right;">
-879.7
</td>
<td style="text-align:right;">
307.1
</td>
<td style="text-align:right;">
3249
</td>
</tr>
</tbody>
</table>

A relatively stable *row or column* in the above table suggests that we
got one of the two nuisance functions ‘right’. In this case, it looks
like the GBM pscore function yields stable estimates across all choices
of outcome models.

Manual use for inference, other estimands
=========================================

the `fit_me` functions fits the `m` functions and `e` function for each
observation and returns a dataset that can then be used for manual
calculations.

``` r
library(data.table)
fit_mod = fit_me(meanfn = 'xgboost', pscorefn = 'xgboost',
    mean_fml = fo, psc_fml  = fp, y = y, w = w, df = df)
setDT(fit_mod)
fit_mod |> head()
```

    ##          y     w    m0      m1     eh
    ##      <num> <num> <num>   <num>  <num>
    ## 1:  9930.0     1  9017  9929.8 0.9657
    ## 2:  3595.9     1 11787  3593.9 1.0012
    ## 3: 24909.5     1  3748 24906.3 0.9887
    ## 4:  7506.1     1  8709  2685.8 1.0012
    ## 5:   289.8     1  2021   291.1 0.9965
    ## 6:  4056.5     1  6751  4060.6 0.9889

trim extreme pscores before AIPW
--------------------------------

``` r
fit_mod |> ate_aipw(c(0.1, 0.9)) |> round(3)
```

    ## [1] 2893

bootstrap
---------

``` r
library(boot); library(MASS)
boot.fn <- function(data, ind){
  d = data[ind, ]
  fit_mod = fit_me(meanfn = 'lasso', pscorefn = 'lasso',
    mean_fml = fo, psc_fml  = fp, y = y, w = w, df = d) |>
    ate_aipw(c(0.1, 0.9))
}
out = boot(df, boot.fn, R = 100)
out |> print()
```

    ## 
    ## ORDINARY NONPARAMETRIC BOOTSTRAP
    ## 
    ## 
    ## Call:
    ## boot(data = df, statistic = boot.fn, R = 100)
    ## 
    ## 
    ## Bootstrap Statistics :
    ##     original  bias    std. error
    ## t1*      342    -101        1023

ATT
---

the ATT subsets to treated units and computes the average between the
realised *Y* and imputed *Y*(0), which can be done easily with our
estimates.

``` r
fit_mod[w == 1, mean(y - m0)] |> round(3)
```

    ## [1] 1796
