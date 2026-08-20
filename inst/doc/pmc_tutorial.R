## ----include = FALSE----------------------------------------------------------

run_everything = suppressWarnings(tryCatch({isTRUE(as.logical(readLines("run_everything.txt")))}, error = function(e){FALSE}))

knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  eval = run_everything,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# library(tidySEM)
# library(OpenMx)

## -----------------------------------------------------------------------------
# set.seed(10)
# n <- 200
# # Set class-specific means
# class_means <- c(rep(0, floor(0.3 * n)), rep(2, ceiling(0.7 * n)))
# 
# # Simulate continuous indicators
# df <- rnorm(2 * n, mean = rep(class_means, 2))
# df <- data.frame(matrix(df, nrow = n))
# names(df) <- paste0("X", 1:2)

## -----------------------------------------------------------------------------
# res <- mx_profiles(data = df, classes = 1:3)

## ----eval = FALSE-------------------------------------------------------------
# library(future)
# plan(multisession)
# res_blrt <- BLRT(res, replications = 100)
# res_blrt

## ----eval = run_everything, echo = FALSE--------------------------------------
# library(future)
# plan(multisession)
# res_blrt <- BLRT(res, replications = 100)
# write.csv(res_blrt, "pmc_res_blrt.csv", row.names = FALSE)

## ----eval = TRUE, echo=FALSE--------------------------------------------------
res_blrt <- read.csv("pmc_res_blrt.csv", stringsAsFactors = FALSE)
class(res_blrt) <- c("LRT", "data.frame")
attr(res_blrt, "type") <- "Bootstrapped"
res_blrt

## ----eval = FALSE-------------------------------------------------------------
# set.seed(1)
# res_pmc <- pmc(res)

## ----eval = run_everything, echo = FALSE--------------------------------------
# set.seed(1)
# res_pmc <- pmc(res)
# write.csv(res_pmc, "pmc_res_pmc.csv", row.names = FALSE)

## ----eval = TRUE, echo=FALSE--------------------------------------------------
res_pmc <- read.csv("pmc_res_pmc.csv", stringsAsFactors = FALSE)
class(res_pmc) <- c("pmc_df", "data.frame")
attr(res_pmc, "stat") <- "SRMR"
res_pmc

## ----eval = FALSE-------------------------------------------------------------
# # Convert the indicators to ordinal
# df[] <- lapply(df, cut, breaks = 3, labels = FALSE)
# df[] <- lapply(df, mxFactor, levels = 1:3)
# res_cat <- mx_lca(df, classes = 1:3)
# pmc(res_cat, reps = 20)

## ----eval = run_everything, echo = FALSE--------------------------------------
# df[] <- lapply(df, cut, breaks = 3, labels = FALSE)
# df[] <- lapply(df, mxFactor, levels = 1:3)
# res_cat <- mx_lca(df, classes = 1:3)
# chisq_function <- function(x, y){
#   tab_obs <- table(x)
#   tab_sim <- table(y)
#   tab_sim[tab_sim == 0] <- NA
#   sum((tab_obs - tab_sim)^2 / tab_sim, na.rm = TRUE)
# }
# res_pmc <- pmc(res_cat, reps = 20)
# write.csv(res_pmc, "pmc_res_chi2.csv", row.names = FALSE)

## ----eval = TRUE, echo=FALSE--------------------------------------------------
res_pmc <- read.csv("pmc_res_chi2.csv", stringsAsFactors = FALSE)
class(res_pmc) <- c("pmc_df", "data.frame")
attr(res_pmc, "stat") <- "chi squared"
res_pmc

