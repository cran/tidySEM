## ----setup, include = FALSE---------------------------------------------------
library(yaml)
library(scales)

run_openmx <- requireNamespace("OpenMx", quietly = TRUE)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  results = "hide",
  tidy.opts = list(width.cutoff = 60),
  tidy = requireNamespace("formatR", quietly = TRUE)
)

options(scipen = 1, digits = 2)

run_everything <- isTRUE(as.logical(toupper(
  Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false")
)))

show_results <- isTRUE(as.logical(toupper(
  Sys.getenv("TIDYSEM_SHOW_RESULTS", "false")
)))

if (show_results) {
  knitr::opts_chunk$set(results = "markup")
}

if (!run_openmx) {
  run_everything <- show_results <- FALSE
}

## -----------------------------------------------------------------------------
library(tidySEM)

set.seed(123)
x <- data.frame(
  depression = c(
    rnorm(50, mean = 6, sd = 4),
    rbeta(100, 2, 3)
  )
)
hist(x$depression)

## ----res-depression, eval = run_everything || show_results--------------------
# res <- tidySEM:::run_cached(mx_profiles(
#   x,
#   classes = 1:3
# ))

## ----eval = show_results------------------------------------------------------
# table_fit(res)

## ----blrt-depression, echo = FALSE, eval = run_everything || show_results-----
# res_blrt <- tidySEM:::run_cached(BLRT(res))

## ----echo = TRUE, eval = FALSE------------------------------------------------
# res_blrt <- BLRT(res)

## ----eval = show_results------------------------------------------------------
# res_blrt

## ----eval = FALSE, echo = TRUE------------------------------------------------
# res_theory <- mx_profiles(
#   x,
#   classes = 2,
#   run = FALSE
# )
# library(OpenMx)
# omxGetParameters(res_theory)

## ----runrestheory, echo = FALSE, eval = run_everything || show_results--------
# res_theory <- tidySEM:::run_cached(mx_profiles(
#   x,
#   classes = 2,
#   run = FALSE
# ))
# library(OpenMx)
# omxGetParameters(res_theory)

## ----eval = run_everything || show_results------------------------------------
# res_theory <- omxSetParameters(
#   res_theory,
#   labels = c("m11", "m21"),
#   values = c(0, 10),
#   free = FALSE
# )

## ----eval = run_everything || show_results------------------------------------
# res_theory_bounds <- omxSetParameters(
#   res_theory,
#   labels = c("m11", "m21"),
#   values = c(0, 10),
#   free = TRUE,
#   lbound = c(NA, 10),
#   ubound = c(1, NA)
# )

## ----theory-fits, echo = FALSE, eval = run_everything || show_results---------
# theory_fits <- tidySEM:::run_cached({
#   res_theory <- run_mx(res_theory)
#   res_theory_bounds <- run_mx(res_theory_bounds)
#   list(theory = res_theory, theory_bounds = res_theory_bounds)
# })
# res_theory <- theory_fits$theory
# res_theory_bounds <- theory_fits$theory_bounds

## ----echo = TRUE, eval = FALSE------------------------------------------------
# res_theory <- run_mx(res_theory)
# res_theory_bounds <- run_mx(res_theory_bounds)

## ----eval = run_everything || show_results------------------------------------
# compare_these <- list(
#   baseline = res[[2]],
#   theory = res_theory,
#   theory_bounds = res_theory_bounds
# )
# 
# table_fit(compare_these)

## ----blrt-theory, echo = FALSE, eval = run_everything || show_results---------
# blrt_theory <- tidySEM:::run_cached(BLRT(list(res[[2]], res_theory)))

## ----echo = TRUE, eval = FALSE------------------------------------------------
# blrt_theory <- BLRT(list(res[[2]], res_theory))

## ----eval = show_results------------------------------------------------------
# blrt_theory

## ----theory-starts, echo = FALSE, eval = run_everything || show_results-------
# res_theory_starts <- omxSetParameters(
#   res_theory,
#   labels = c("m11", "m21"),
#   values = c(0, 10),
#   free = TRUE
# )
# res_theory_starts <- tidySEM:::run_cached(run_mx(res_theory_starts))

## ----eval = FALSE, echo = TRUE------------------------------------------------
# res_theory_starts <- omxSetParameters(
#   res_theory,
#   labels = c("m11", "m21"),
#   values = c(0, 10),
#   free = TRUE
# )
# res_theory_starts <- run_mx(res_theory_starts)

## ----eval = run_everything || show_results------------------------------------
# omxGetParameters(res_theory_starts)
# wald_test(res_theory_starts, "m11 = 0; m21 = 10")

## ----echo = TRUE, eval=run_openmx---------------------------------------------
# Load required packages
library(tidySEM) 
library(ggplot2)
# Load data
df <- zegwaard_carecompass[, c("burdened", "trapped", "negaffect", "loneliness")]

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# desc <- tidySEM::descriptives(df)
# desc <- desc[, c("name", "n", "missing", "unique",
#                  "mean", "median", "sd", "min", "max",
#                  "skew_2se", "kurt_2se")]
# desc

## ----tabdesc, echo = FALSE, eval=run_openmx-----------------------------------
desc <- tidySEM::descriptives(df)
desc <- desc[, c("name", "n", "missing", "unique",
                 "mean", "median", "sd", "min", "max",
                 "skew_2se", "kurt_2se")]
knitr::kable(desc, caption = "Descriptive statistics")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# df_plot <- df
# names(df_plot) <- paste0("Value.", names(df_plot))
# df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                    timevar = "Variable")
# ggplot(df_plot, aes(x = Value)) +
#   geom_density() +
#   facet_wrap(~Variable)+
#   theme_bw()

## ----figdesc, echo = FALSE, eval = run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
#   # mcartest <- mice::mcar(df)
#   df_plot <- df
#   names(df_plot) <- paste0("Value.", names(df_plot))
#   df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                      timevar = "Variable")
#   ggplot(df_plot, aes(x = Value)) +
#     geom_density() +
#     facet_wrap(~Variable) +
#     theme_bw()
# }, type = "png", width = 100 / 25.4, height = 100 / 25.4)

## ----fitlca, eval = run_everything || show_results, echo = FALSE--------------
# res <- tidySEM:::run_cached({
#   set.seed(123) # setting seed
#   mx_profiles(data = df,
#               classes = 1:5)
# })
# 

## ----eval = FALSE, echo = TRUE------------------------------------------------
# set.seed(123)
# res <- mx_profiles(data = df,
#                    classes = 1:5)

## ----eval = run_everything || show_results, echo = TRUE-----------------------
# fit <- table_fit(res)

## ----fit_table, include = TRUE, eval=F----------------------------------------
# fit <- table_fit(res) # model fit table
# fit[ , c("Name", "LL", "Parameters", "n",
#          "BIC", "Entropy",
#          "prob_min", "prob_max",
#          "n_min", "n_max",
#          "np_ratio", "np_local")]

## ----tabfit, echo = FALSE, eval = show_results--------------------------------
# tbl <- fit[ , c("Name", "LL", "Parameters", "n",
#          "BIC", "Entropy",
#          "prob_min", "prob_max",
#          "n_min", "n_max")]
# names(tbl) <- c("Name", "LL", "p", "n",
#          "BIC", "Entropy",
#          "p_min", "p_max",
#          "n_min", "n_max")
# knitr::kable(tbl, caption = "Model fit table")

## ----tablmr, echo = FALSE, eval = run_everything || show_results--------------
# tbl <- tidySEM:::run_cached(lr_lmr(res))

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# lr_lmr(res)

## ----echo = FALSE, eval = show_results----------------------------------------
# knitr::kable(tbl, caption = "LMR test table", digits = 2)

## ----caregiver-blrt, echo = FALSE, eval = run_everything || show_results------
# res_blrt <- tidySEM:::run_cached({
#   library(future)
#   library(progressr)
#   plan(multisession) # Parallel processing for Windows
#   handlers("progress") # Progress bar
#   set.seed(1)
#   BLRT(res, replications = 100)
#   })

## ----eval = FALSE, echo = TRUE------------------------------------------------
# library(future)
# library(progressr)
# plan(multisession) # Parallel processing for Windows
# handlers("progress") # Progress bar
# set.seed(1)
# res_blrt <- BLRT(res, replications = 100)
# res_blrt

## ----eval = show_results, echo = FALSE----------------------------------------
# knitr::kable(res_blrt, caption = "BLRT test table", digits = 2)

## ----alternative-fit, eval = run_everything || show_results, echo = FALSE-----
# fit_compare <- tidySEM:::run_cached({
#   res_alt <- mx_profiles(df, classes = 4, variances = "varying")
#   compare <- list(res[[4]], res_alt)
#   table_fit(compare)
# })

## ----echo = TRUE, eval = FALSE------------------------------------------------
# res_alt <- mx_profiles(df, classes = 4, variances = "varying")
# compare <- list(res[[4]], res_alt)
# table_fit(compare)

## ----tabfitcomp, echo = FALSE, eval = show_results----------------------------
# knitr::kable(fit_compare[ , c("Name", "LL", "Parameters",
#          "BIC", "Entropy",
#          "prob_min", "prob_max",
#          "n_min", "n_max")], caption = "Comparing competing theoretical models")

## ----echo = TRUE, eval =FALSE-------------------------------------------------
# res_final <- mx_switch_labels(res[[4]])

## ----final-class, echo = FALSE, eval = run_everything || show_results---------
# res_final <- mx_switch_labels(res[[4]])
# cp <- class_prob(res_final)
# out <- list(counts = cp$sum.posterior$proportion)

## ----final-results, echo = FALSE, eval = run_everything || show_results-------
# tab <- tidySEM:::run_cached(table_results(res_final, columns = c("label", "est", "se", "confint", "class")))

## ----echo = TRUE, eval = FALSE------------------------------------------------
# table_results(res_final, columns = c("label", "est", "se", "confint", "class"))

## ----eval = show_results, echo=FALSE------------------------------------------
# knitr::kable(tab, caption = "Four-class model results")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# plot_bivariate(res_final)

## ----lpa-bivariate-plot, echo = FALSE, eval = run_everything || show_results, fig.cap="Bivariate profile plot", out.width="100%"----
# tidySEM:::run_cached({
#   p <- plot_bivariate(res_final, return_list = TRUE)
#   p[[1]] <- p[[1]] + scale_y_continuous(breaks = c(0, .1, .2, .3), labels = c(".0", ".1", ".2", ".3"))
#   tidySEM:::merge_corplots(p)
# }, type = "png", width = 315 / 25.4, height = 150 / 25.4)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# plot_profiles(res_final)

## ----lpa-profiles-plot, echo = FALSE, eval = run_everything || show_results, fig.cap="Profile plot", out.width="100%"----
# tidySEM:::run_cached(
#   plot_profiles(res_final),
#   type = "png",
#   width = 100 / 25.4,
#   height = 100 / 25.4
# )

## ----echo = TRUE, eval = run_everything || show_results-----------------------
# df_new <- data.frame(
#   burdened = 2,
#   trapped = 0.5,
#   negaffect = 1.5,
#   loneliness = 4
# )
# predict_class(res_final, newdata = df_new)

