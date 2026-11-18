## ----setup, include = FALSE---------------------------------------------------
run_openmx <- requireNamespace("OpenMx", quietly = TRUE)
run_everything <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false")))) && run_openmx
show_results <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_SHOW_RESULTS", "false")))) && run_openmx
run_chunks <- run_openmx && (run_everything || show_results)
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  warning = FALSE,
  comment = "#>",
  out.width = "100%"
)

## ----eval = run_chunks--------------------------------------------------------
# library(tidySEM)
# library(OpenMx)

## ----eval = run_chunks--------------------------------------------------------
# set.seed(10)
# n <- 200
# 
# # Set class-specific means
# class_means <- c(rep(0, floor(0.3 * n)),
# rep(2, ceiling(0.7 * n)))
# 
# # Simulate continuous indicators
# df <- rnorm(4 * n, mean = rep(class_means, 4))
# df <- matrix(df, nrow = n)
# df <- t(t(df) * c(1, 2, 0.5, 1))
# df <- data.frame(df)
# names(df) <- paste0("X", 1:4)
# 
# # Convert one indicator to ordinal
# df$X4 <- cut(df$X4, breaks = 3, labels = FALSE)
# df$X4 <- mxFactor(df$X4, levels = 1:3)

## ----eval = FALSE-------------------------------------------------------------
# res_2 <- mx_mixed_lca(
# data = df,
# classes = 2
# )

## ----mixed-fit-two, echo = FALSE, eval = run_chunks---------------------------
# res_2 <- tidySEM:::run_cached(mx_mixed_lca(data = df, classes = 2))

## ----eval = show_results------------------------------------------------------
# class(res_2)

## ----eval = FALSE-------------------------------------------------------------
# res_1_3 <- mx_mixed_lca(
#   data = df,
#   classes = 1:3
# )

## ----mixed-fit-one-three, eval = run_chunks, echo = FALSE---------------------
# res_1_3 <- tidySEM:::run_cached(mx_mixed_lca(data = df, classes = 1:3))

## ----eval = FALSE-------------------------------------------------------------
# table_fit(res_1_3)

## ----mixed-fit-table, eval = run_chunks, echo = FALSE-------------------------
# tab <- tidySEM:::run_cached(table_fit(res_1_3))
# if (show_results) tab

## ----eval = FALSE, echo = TRUE------------------------------------------------
# library(future)
# library(progressr)
# plan(multisession) # Parallel processing for Windows
# handlers("progress") # Progress bar
# set.seed(1)
# res_blrt <- BLRT(res_1_3, replications = 100)
# res_blrt

## ----mixed-blrt, eval = run_chunks, echo = FALSE------------------------------
# res_blrt <- tidySEM:::run_cached({
#   library(future)
#   library(progressr)
#   plan(multisession) # Parallel processing for Windows
#   handlers("progress") # Progress bar
#   set.seed(1)
#   BLRT(res_1_3, replications = 100)
# })
# if (show_results) res_blrt

## ----mixed-pmc, eval = run_chunks, echo = FALSE-------------------------------
# res_pmc <- tidySEM:::run_cached({
#   set.seed(1)
#   pmc(res_1_3)
# })
# if (show_results) res_pmc

## ----echo = TRUE, eval = FALSE------------------------------------------------
# set.seed(1)
# res_pmc <- pmc(res_1_3)
# res_pmc

## ----mixed-class-prob, eval = run_chunks, echo = FALSE------------------------
# class_probs <- tidySEM:::run_cached(class_prob(res_1_3[[2]], c("sum.posterior", "sum.mostlikely")))
# if (show_results) class_probs

## ----eval = FALSE, echo = TRUE------------------------------------------------
# class_prob(res_1_3[[2]], c("sum.posterior", "sum.mostlikely"))

## ----eval = FALSE, echo = TRUE------------------------------------------------
# table_fit(res_1_3[[2]])

## ----mixed-fit-two-table, eval = run_chunks, echo = FALSE---------------------
# fit_two <- tidySEM:::run_cached(table_fit(res_1_3[[2]]))
# if (show_results) fit_two

## ----eval = FALSE, echo = TRUE------------------------------------------------
# table_results(res_1_3[[2]])

## ----mixed-results-two, eval = run_chunks, echo = FALSE-----------------------
# results_two <- tidySEM:::run_cached(table_results(res_1_3[[2]]))
# if (show_results) results_two

## ----eval = FALSE, echo = TRUE------------------------------------------------
# table_prob(res_1_3[[2]])

## ----mixed-prob-two, eval = run_chunks, echo = FALSE--------------------------
# prob_two <- tidySEM:::run_cached(table_prob(res_1_3[[2]]))
# if (show_results) prob_two

## ----eval = FALSE-------------------------------------------------------------
# res_2_free <- mx_mixed_lca(
#   data = df,
#   classes = 2,
#   variances = "varying"
# )

## ----mixed-free-variances, eval = run_chunks, echo = FALSE--------------------
# res_2_free <- tidySEM:::run_cached(mx_mixed_lca(data = df, classes = 2, variances = "varying"))

## ----echo = TRUE, eval = FALSE------------------------------------------------
# compare <- list(
#   fixed_covs = res_1_3[[2]],
#   free_covs = res_2_free)
# table_fit(compare)

## ----tabfitcomp, echo = FALSE, eval = run_chunks------------------------------
# fit_compare <- tidySEM:::run_cached({
#   compare <- list(fixed_covs = res_1_3[[2]], free_covs = res_2_free)
#   table_fit(compare)
# })
# if (show_results) fit_compare

## ----eval = FALSE, echo = TRUE------------------------------------------------
# plot_profiles(res_1_3[[2]], variables = c("X1", "X2", "X3"))

## ----mixed-plot-profiles, eval = run_chunks, echo = FALSE, out.width = "100%"----
# tidySEM:::run_cached(
#   plot_profiles(res_1_3[[2]], variables = c("X1", "X2", "X3")),
#   type = "png", width = 4, height = 3, dpi = 150
# )

## ----eval = FALSE, echo = TRUE------------------------------------------------
# plot_bivariate(res_1_3[[2]], variables = c("X1", "X2", "X3"))

## ----mixed-plot-bivariate, eval = run_chunks, echo = FALSE, out.width = "100%"----
# tidySEM:::run_cached(
#   plot_bivariate(res_1_3[[2]], variables = c("X1", "X2", "X3")),
#   type = "png", width = 4, height = 4, dpi = 150
# )

## ----eval = FALSE, echo = TRUE------------------------------------------------
# plot_prob(res_1_3[[2]])

## ----mixed-plot-prob, eval = run_chunks, echo = FALSE, out.width = "100%"-----
# tidySEM:::run_cached(
#   plot_prob(res_1_3[[2]]),
#   type = "png", width = 3, height = 3, dpi = 150
# )

