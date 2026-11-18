## ----setup, include = FALSE---------------------------------------------------
run_everything <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false"))))
show_results <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_SHOW_RESULTS", "false"))))
run_openmx <- requireNamespace("OpenMx", quietly = TRUE)
if (!run_openmx) run_everything <- show_results <- FALSE

knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>",
  out.width = "100%"
)


## ----eval = run_openmx--------------------------------------------------------
library(tidySEM)
library(OpenMx)

## ----eval = run_openmx--------------------------------------------------------
set.seed(10)
n <- 200
# Set class-specific means
class_means <- c(rep(0, floor(0.3 * n)), rep(2, ceiling(0.7 * n)))

# Simulate continuous indicators
df <- rnorm(2 * n, mean = rep(class_means, 2))
df <- data.frame(matrix(df, nrow = n))
names(df) <- paste0("X", 1:2)

## ----profiles-example, eval = FALSE-------------------------------------------
# res <- mx_profiles(data = df, classes = 1:3)

## ----profiles-cache, echo = FALSE, eval = run_openmx && (run_everything || show_results), results = "hide"----
# res <- tidySEM:::run_cached(mx_profiles(data = df, classes = 1:3))

## ----eval = FALSE-------------------------------------------------------------
# library(future)
# plan(multisession)
# res_blrt <- BLRT(res, replications = 100)
# res_blrt

## ----blrt-cache, eval = run_openmx && (run_everything || show_results), echo = FALSE, results = "hide"----
# res_blrt <- tidySEM:::run_cached({
#   library(future)
#   plan(multisession)
#   BLRT(res, replications = 100)
# })

## ----blrt-results, eval = run_openmx && show_results, echo = FALSE, results = "markup"----
# res_blrt

## ----eval = FALSE-------------------------------------------------------------
# set.seed(1)
# res_pmc <- pmc(res)

## ----pmc-continuous-cache, eval = run_openmx && (run_everything || show_results), echo = FALSE, results = "hide"----
# res_pmc <- tidySEM:::run_cached({
#   set.seed(1)
#   pmc(res)
# })

## ----pmc-continuous-results, eval = run_openmx && show_results, echo = FALSE, results = "markup"----
# res_pmc

## ----eval = FALSE-------------------------------------------------------------
# # Convert the indicators to ordinal
# df[] <- lapply(df, cut, breaks = 3, labels = FALSE)
# df[] <- lapply(df, mxFactor, levels = 1:3)
# res_cat <- mx_lca(df, classes = 1:3)
# pmc(res_cat, reps = 20)

## ----pmc-ordinal-cache, eval = run_openmx && (run_everything || show_results), echo = FALSE, results = "hide"----
# res_pmc <- tidySEM:::run_cached({
#   df[] <- lapply(df, cut, breaks = 3, labels = FALSE)
#   df[] <- lapply(df, mxFactor, levels = 1:3)
#   res_cat <- mx_lca(df, classes = 1:3)
#   chisq_function <- function(x, y){
#     tab_obs <- table(x)
#     tab_sim <- table(y)
#     tab_sim[tab_sim == 0] <- NA
#     sum((tab_obs - tab_sim)^2 / tab_sim, na.rm = TRUE)
#   }
#   pmc(res_cat, reps = 20)
# })

## ----pmc-ordinal-results, eval = run_openmx && show_results, echo = FALSE, results = "markup"----
# res_pmc

