## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
if(!requireNamespace("OpenMx", quietly = TRUE)){
  knitr::opts_chunk$set(eval = FALSE)
}

## ----data---------------------------------------------------------------------
library(tidySEM)
library(OpenMx)
set.seed(1)
df <- data.frame(rbind(
  matrix(rnorm(300, mean = 2), ncol = 3),
  matrix(rnorm(600), ncol = 3)
))
df$X1 <- ordered(cut(df$X1, 2, labels = FALSE))
df$X2 <- ordered(cut(df$X2, 2, labels = FALSE))
df$X3 <- ordered(cut(df$X3, 3, labels = FALSE))

## ----threshold-model----------------------------------------------------------
mod <- as_ram("
X1 | t1
X2 | t1
X3 | t1
X3 | t2", data = df)

## ----threshold-constrain------------------------------------------------------
mod_constrain <- as_ram("
X1 | c*t1
X2 | c*t1
X3 | t1
X3 | t2", data = df)

## ----fit-default--------------------------------------------------------------
res <- run_mx(mod)

## ----table-default------------------------------------------------------------
table_results(res)

## ----wald-default-x1----------------------------------------------------------
wald_test(res, "model.Thresholds[1,1] = 0")

## ----wald-default-x3----------------------------------------------------------
wald_test(res, "model.Thresholds[2,3] = 1")

## ----threshold-model-deviances------------------------------------------------
mod2 <- as_ram("
X1 | t1
X2 | t1
X3 | t1
X3 | t2", data = df, threshold_method = "mx_deviances")

## ----fit-deviances------------------------------------------------------------
res2 <- run_mx(mod2)
table_results(res2)

## ----wald-deviances-x1--------------------------------------------------------
wald_test(res2, "model.mat_dev[1,1] = 0")

## -----------------------------------------------------------------------------
wald_test(res, "model.Thresholds[1,1] = 0")

## ----wald-deviances-x3--------------------------------------------------------
wald_test(res2, "model.mat_dev[1,3]+model.mat_dev[2,3] = 1")

## -----------------------------------------------------------------------------
wald_test(res, "model.Thresholds[2,3] = 1")

## -----------------------------------------------------------------------------
coef(res2)
res2_thresholds <- deviances_to_thresholds(res2)
coef(res2_thresholds)
wald_test(res2_thresholds, "model.Thresholds[1,1] = 0")

