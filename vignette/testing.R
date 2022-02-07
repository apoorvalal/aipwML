# %% ####################################################
rm(list = ls())
library(causalsens); library(knitr); library(kableExtra)
# ml libraries
library(glmnet); library(ranger); library(xgboost)
# this library
source("../R/fns.R")
# library(aipwML)

# %%
data(lalonde.psid); df = lalonde.psid
y = 're78'; w = 'treat'
x = setdiff(colnames(df), c(y, w))

# outcome model formula
fo = re78 ~ (age + education + black + hispanic + married + nodegree +
    re74 + re75 + u74 + u75)
# pscore formula
fp = treat ~ (age + education + black + hispanic + married + nodegree +
    re74 + re75 + u74 + u75)

# %%

regadjusts = c(
  ate_reg('ols',      w = w, y = y, df = df, fml = fo),
  ate_reg('lasso',    w = w, y = y, df = df, fml = fo),
  ate_reg('ridge',    w = w, y = y, df = df, fml = fo),
  ate_reg('rforest',  w = w, y = y, df = df, fml = fo),
  ate_reg('xgboost',  w = w, y = y, df = df, fml = fo)
)
regadjusts |> round(3)


# %%
ipws = c(
  ate_ipw('logit',   w = w, y = y, df = df, fml = fp),
  ate_ipw('lasso',   w = w, y = y, df = df, fml = fp),
  ate_ipw('ridge',   w = w, y = y, df = df, fml = fp),
  ate_ipw('rforest', w = w, y = y, df = df, fml = fp),
  ate_ipw('xgboost', w = w, y = y, df = df, fml = fp)
)

ipws |> round(3)


# %%

ipws2 = c(
  ate_ipw('logit',   w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('lasso',   w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('ridge',   w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('rforest', w = w, y = y, df = df, fml = fp, psrange = psr),
  ate_ipw('xgboost', w = w, y = y, df = df, fml = fp, psrange = psr)
)

ipws2 |> round(3)


# %%
psr = c(0.05, 0.95)

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


# %%
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


# %%
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


# %%
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

# %%
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


# %% # stack estimates
aipw_estimates = rbind(ols_mean, lasso_mean, ridge_mean, rforest_mean, xgboost_mean)
colnames(aipw_estimates) = c('PS: logit', 'PS: lasso', 'PS: ridge',
  'PS: rforest', 'PS: xgboost')
rownames(aipw_estimates)= c('Outcome: ols', 'Outcome: lasso', 'Outcome: ridge',
  'Outcome: rforest', 'Outcome: xgboost')
aipw_estimates |> kbl() %>%
  kable_styling() |> chr_nb()


# %%
library(data.table)
fit_mod = fit_me(meanfn = 'xgboost', pscorefn = 'xgboost',
    mean_fml = fo, psc_fml  = fp, y = y, w = w, df = df)
setDT(fit_mod)
fit_mod |> head()


# %%
fit_mod |> ate_aipw(c(0.1, 0.9)) |> round(3)


# %%
library(boot); library(MASS)
boot.fn <- function(data, ind){
  d = data[ind, ]
  fit_mod = fit_me(meanfn = 'lasso', pscorefn = 'lasso',
    mean_fml = fo, psc_fml  = fp, y = y, w = w, df = d) |>
    ate_aipw(c(0.1, 0.9))
}
out = boot(df, boot.fn, R = 100)
out |> print()


# %%
fit_mod[w == 1, mean(y - m0)] |> round(3)
