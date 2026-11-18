## ----setup, include = FALSE---------------------------------------------------
library(yaml)
library(scales)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  results = "hide",
  tidy.opts = list(width.cutoff = 60), tidy = requireNamespace("formatR", quietly = TRUE)
)
options(scipen = 1, digits = 2)
run_openmx <- requireNamespace("OpenMx", quietly = TRUE)
run_everything <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false"))))
show_results <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_SHOW_RESULTS", "false"))))
if (!run_openmx) run_everything <- show_results <- FALSE
if (show_results) knitr::opts_chunk$set(results = "markup")


## ----echo = TRUE, eval=TRUE---------------------------------------------------
# Load required packages
library(tidySEM) 
library(ggplot2)
# Load data
# Keep one indicator for each identity domain and collapse adjacent categories.
# This deliberately smaller example makes ordinal mixture estimation faster.
ordinal3 <- function(x) {
  z <- as.integer(x)
  observed <- sort(unique(z[!is.na(z)]))
  if (length(observed) < 3L) stop("Indicator has fewer than three observed categories")
  # Collapse ordered response categories into three approximately equal ranges.
  groups <- cut(seq_along(observed), breaks = 3L, labels = FALSE)
  ordered(groups[match(z, observed)], levels = 1:3)
}
df <- maene_identity[c("Ethnic_1", "Belgian", "Flemish")]
df[] <- lapply(df, ordinal3)
# Auxiliary variables must correspond to the same rows as the indicators.
aux_data <- maene_identity

## ----echo = TRUE, eval=FALSE, results='asis'----------------------------------
# desc <- tidySEM::descriptives(df)
# desc <- desc[, c("name", "type", "n", "missing", "unique",
# "mode", "mode_value", "v")
# ]
# desc

## ----cache-desc, echo = FALSE, eval = run_everything || show_results----------
# desc <- tidySEM:::run_cached({
#   desc <- tidySEM::descriptives(df)
#   desc[, c("name", "type", "n", "missing", "unique", "mode", "mode_value", "v")]
# })

## ----tabdesc, echo = FALSE, eval = show_results-------------------------------
# knitr::kable(desc, caption = "Descriptive statistics for recoded ordinal items")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# df_plot <- df
# names(df_plot) <- paste0("Value.", names(df_plot))
# df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long")
# ggplot(df_plot, aes(x = Value)) +
#   geom_bar() +
#   facet_wrap(~time, scales = "free")+
#   theme_bw()

## ----figdesc, echo = FALSE, eval = run_everything || show_results, fig.cap="Bar charts for ordinal indicators", out.width="100%"----
# tidySEM:::run_cached({
#   df_plot <- df
#   names(df_plot) <- paste0("Value.", names(df_plot))
#   df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long")
#   ggplot(df_plot, aes(x = Value)) +
#     geom_bar() + facet_wrap(~time, scales = "free") + theme_bw()
# }, type = "png", width = 100 / 25.4, height = 100 / 25.4)

## ----fitlca, eval = run_everything || show_results, echo = FALSE--------------
# res <- tidySEM:::run_cached({
#   set.seed(123)
#   mx_lca(data = df, classes = 1:3)
# })

## ----eval = FALSE, echo = TRUE------------------------------------------------
# set.seed(123)
# res <- mx_lca(data = df, classes = 1:3)

## ----echo = TRUE, eval=F------------------------------------------------------
# fit <- table_fit(res)
# fit[ , c("Name", "LL", "n", "Parameters", "BIC",
# "Entropy", "prob_min", "n_min", "np_ratio",
# "np_local")]

## ----cache-fit, echo = FALSE, eval = run_everything || show_results-----------
# fit <- tidySEM:::run_cached(table_fit(res))

## ----tabfit, echo = FALSE, eval = show_results--------------------------------
# knitr::kable(fit[ , c("Name", "LL", "n", "Parameters", "BIC",
# "Entropy", "prob_min", "n_min", "np_ratio", "np_local")],
# caption = "Model fit table", digits = 2)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# res_final <- res[[3]]

## ----select-final, echo = FALSE, eval = run_everything || show_results--------
# res_final <- res[[3]]

## ----echo = TRUE, eval = FALSE------------------------------------------------
# tab <- table_prob(res_final)
# reshape(tab, direction = "wide", v.names = "Probability", timevar = "group", idvar = c("Variable", "Category"))

## ----cache-prob-table, echo = FALSE, eval = run_everything || show_results----
# tab <- tidySEM:::run_cached({
#   tab <- table_prob(res_final)
#   reshape(tab, direction = "wide", v.names = "Probability",
#           timevar = "group", idvar = c("Variable", "Category"))
# })

## ----tabprob, eval = show_results, echo = FALSE-------------------------------
# knitr::kable(tab, caption = "Three-class model results in probability scale")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# plot_prob(res_final, bw = TRUE)

## ----figprob, echo = FALSE, eval = run_everything || show_results, fig.cap="Probability plot", out.width="100%"----
# tidySEM:::run_cached(
#   plot_prob(res_final, bw = TRUE),
#   type = "png",
#   width = 210 / 25.4 * 1.5,
#   height = 100 / 25.4 * 1.5
# )

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# aux_dep <- BCH(res_final, data = aux_data$depression)

## ----cache-aux-dep, echo = FALSE, eval = run_everything || show_results-------
# aux_dep <- tidySEM:::run_cached(BCH(res_final, data = aux_data$depression))

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# df_aux <- aux_data[, c("vict_teacher", "depression")]
# # Dummy-code vict_teacher
# df_aux$vict_teacher <- (as.integer(df_aux$vict_teacher)-1)
# aux_model <- BCH(res_final,
#                  model = "depression ~ vict_teacher",
#                  data = df_aux)

## ----cache-aux-model, echo = FALSE, eval = run_everything || show_results-----
# aux_model <- tidySEM:::run_cached({
#   df_aux <- aux_data[, c("vict_teacher", "depression")]
#   df_aux$vict_teacher <- (as.integer(df_aux$vict_teacher)-1)
#   BCH(res_final, model = "depression ~ vict_teacher", data = df_aux)
# })

## ----echo=TRUE, eval=FALSE----------------------------------------------------
# table_results(aux_model, columns = c("label", "est_sig", "se", "pval", "confint", "group", "matrix"))

## ----aux-results, echo = FALSE, eval = show_results, results = "markup"-------
# table_results(aux_model, c("label", "est_sig", "se", "pval", "confint", "group", "matrix"))

## ----echo=TRUE, eval=FALSE----------------------------------------------------
# lr_test(aux_model, compare = "A")

## ----aux-lr, echo = FALSE, eval = show_results, results = "markup"------------
# lr_test(aux_model, compare = "A")

## ----echo=TRUE, eval=FALSE----------------------------------------------------
# wald_test(
#   aux_model,
#   "aux1.A[1,2]=aux2.A[1,2]&aux1.A[1,2]=aux3.A[1,2]"
# )

## ----aux-wald, echo = FALSE, eval = show_results, results = "markup"----------
# wald_test(aux_model, "aux1.A[1,2]=aux2.A[1,2]&aux1.A[1,2]=aux3.A[1,2]")

