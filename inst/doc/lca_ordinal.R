## ----setup, include = FALSE---------------------------------------------------
library(yaml)
library(scales)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE,
  results = "hide",
  tidy.opts = list(width.cutoff = 60), tidy = TRUE
)
options(scipen = 1, digits = 2)
eval_results <- FALSE
if(suppressWarnings(tryCatch({isTRUE(as.logical(readLines("pkgdown.txt")))}, error = function(e){FALSE}))){
  eval_results <- TRUE
  knitr::opts_chunk$set(
  results = "markup"
)
}
run_everything = suppressWarnings(tryCatch({isTRUE(as.logical(readLines("run_everything.txt")))}, error = function(e){FALSE}))

## ----echo = TRUE, eval=TRUE---------------------------------------------------
# Load required packages
library(tidySEM) 
library(ggplot2)
# Load data
df <- maene_identity[1:5]

## ----echo = TRUE, eval=FALSE, results='asis'----------------------------------
#  desc <- tidySEM::descriptives(df)
#  desc <- desc[, c("name", "type", "n", "missing", "unique",
#  "mode", "mode_value", "v")
#  ]
#  desc

## ----echo = FALSE, eval=run_everything----------------------------------------
#  desc <- tidySEM::descriptives(df)
#  desc <- desc[, c("name", "type", "n", "missing", "unique",
#  "mode", "mode_value", "v")
#  ]
#  write.csv(desc, "lca_desc.csv", row.names = FALSE)

## ----tabdesc, echo = FALSE, eval=eval_results---------------------------------
#  desc <- read.csv("lca_desc.csv", stringsAsFactors = FALSE)
#  knitr::kable(desc, caption = "Descriptive statistics for ordinal items")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df_plot <- df
#  names(df_plot) <- paste0("Value.", names(df_plot))
#  df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long")
#  ggplot(df_plot, aes(x = Value)) +
#    geom_bar() +
#    facet_wrap(~time, scales = "free")+
#    theme_bw()

## ----echo = FALSE, eval = run_everything--------------------------------------
#  df_plot <- df
#  names(df_plot) <- paste0("Value.", names(df_plot))
#  df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long")
#  p = ggplot(df_plot, aes(x = Value)) +
#    geom_bar() +
#    facet_wrap(~time, scales = "free")+
#    theme_bw()
#  ggsave("lca_plot_desc.png", p, device = "png", width = 100, height = 100, units = "mm")

## ----echo = FALSE, eval=eval_results, fig.cap="Bar charts for ordinal indicators", out.width="80%"----
#  knitr::include_graphics("lca_plot_desc.png")

## ----fitlca, eval = run_everything, echo = FALSE------------------------------
#  set.seed(123)
#  res <- mx_lca(data = df, classes = 1:6)
#  saveRDS(res, "lca_res.RData")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  set.seed(123)
#  res <- mx_lca(data = df, classes = 1:6)

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  res <- readRDS("lca_res.RData")

## ----echo = TRUE, eval=F------------------------------------------------------
#  fit <- table_fit(res)
#  fit[ , c("Name", "LL", "n", "Parameters", "BIC",
#  "Entropy", "prob_min", "n_min", "np_ratio",
#  "np_local")]

## ----echo = FALSE, eval = run_everything--------------------------------------
#  fit <- table_fit(res)
#  write.csv(fit, "lca_fit.csv", row.names = FALSE)

## ----tabfit, echo = FALSE, eval = eval_results--------------------------------
#  fit <- read.csv("lca_fit.csv", stringsAsFactors = FALSE)
#  class(fit) <- c("tidy_fit", "data.frame")
#  knitr::kable(fit[ , c("Name", "LL", "n", "Parameters", "BIC",
#  "Entropy", "prob_min", "n_min", "np_ratio",
#  "np_local")], caption = "Model fit table", digits = 2)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  res_final <- res[[3]]

## ----echo = FALSE, eval = run_everything--------------------------------------
#  res_final <- res[[3]]

## ----echo = FALSE, eval = run_everything--------------------------------------
#  cp <- class_prob(res_final)
#  out <- list(counts = cp$sum.posterior$proportion)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  tab <- table_prob(res_final)
#  reshape(tab, direction = "wide", v.names = "Probability", timevar = "group", idvar = c("Variable", "Category"))

## ----echo = FALSE, eval = run_everything--------------------------------------
#  tab <- table_prob(res_final)
#  tab <- reshape(tab, direction = "wide", v.names = "Probability", timevar = "group", idvar = c("Variable", "Category"))
#  write.csv(tab, "lca_tab_prob.csv", row.names = FALSE)

## ----eval = eval_results, echo=FALSE------------------------------------------
#  tab <- read.csv("lca_tab_prob.csv", stringsAsFactors = FALSE)
#  knitr::kable(tab, caption = "Three-class model results in probability scale")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  plot_prob(res_final, bw = TRUE)

## ----echo = FALSE, eval = run_everything--------------------------------------
#  p <- plot_prob(res_final, bw = TRUE)
#  ggsave("lca_prob.png", p, device = "png", width = 210, height = 100, units = "mm", scale = 1.5)

## ----echo = FALSE, eval = eval_results, fig.cap="Probability plot", out.width="100%"----
#  knitr::include_graphics("lca_prob.png")

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  aux_dep <- BCH(res_final, data = maene_identity$depression)

## ----echo = FALSE, eval=run_everything----------------------------------------
#  aux_dep <- BCH(res_final, data = maene_identity$depression)

## ----echo = FALSE, eval=FALSE-------------------------------------------------
#  aux_dep <- readRDS("lca_aux_dep.RData")

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  df_aux <- maene_identity[, c("vict_teacher", "depression")]
#  # Dummy-code vict_teacher
#  df_aux$vict_teacher <- (as.integer(df_aux$vict_teacher)-1)
#  aux_model <- BCH(res_final,
#                   model = "depression ~ vict_teacher",
#                   data = df_aux)

## ----echo = FALSE, eval=run_everything----------------------------------------
#  df_aux <- maene_identity[, c("vict_teacher", "depression")]
#  # Dummy-code vict_teacher
#  df_aux$vict_teacher <- (as.integer(df_aux$vict_teacher)-1)
#  aux_model <- BCH(res_final,
#                   model = "depression ~ vict_teacher",
#                   data = df_aux)

## ----echo = FALSE, eval=FALSE-------------------------------------------------
#  aux_model <- readRDS("lca_aux_model.RData")

