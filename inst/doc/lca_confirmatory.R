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
df <- zegwaard_carecompass[, c("burdened", "trapped", "negaffect", "loneliness")]

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  desc <- tidySEM::descriptives(df)
#  desc <- desc[, c("name", "n", "missing", "unique",
#                   "mean", "median", "sd", "min", "max",
#                   "skew_2se", "kurt_2se")]
#  desc

## ----tabdesc, echo = FALSE, eval=TRUE-----------------------------------------
desc <- tidySEM::descriptives(df)
desc <- desc[, c("name", "n", "missing", "unique",
                 "mean", "median", "sd", "min", "max",
                 "skew_2se", "kurt_2se")]
knitr::kable(desc, caption = "Descriptive statistics")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df_plot <- df
#  names(df_plot) <- paste0("Value.", names(df_plot))
#  df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                     timevar = "Variable")
#  ggplot(df_plot, aes(x = Value)) +
#    geom_density() +
#    facet_wrap(~Variable)+
#    theme_bw()

## ----echo = FALSE, eval = run_everything--------------------------------------
#  # mcartest <- mice::mcar(df)
#  df_plot <- df
#  names(df_plot) <- paste0("Value.", names(df_plot))
#  df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                     timevar = "Variable")
#  p <- ggplot(df_plot, aes(x = Value)) +
#    geom_density() +
#    facet_wrap(~Variable)+
#    theme_bw()
#  ggsave("plot_lpa_desc.png", p, device = "png", width = 100, height = 100, units = "mm")

## ----figdesc, echo = FALSE, eval = eval_results, out.width="80%"--------------
#  knitr::include_graphics("plot_lpa_desc.png")

## ----fitlca, eval = run_everything, echo = FALSE------------------------------
#  set.seed(123) # setting seed
#  res <- mx_profiles(data = df,
#                     classes = 1:5)
#  saveRDS(res, "res_lpa.RData")
#  fit <- table_fit(res)
#  write.csv(fit, "lpatabfit.csv", row.names = FALSE)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  set.seed(123)
#  res <- mx_profiles(data = df,
#                     classes = 1:5)

## ----eval = eval_results, echo = FALSE----------------------------------------
#  fit <- read.csv("lpatabfit.csv", stringsAsFactors = FALSE)
#  class(fit) <- c("tidy_fit", "data.frame")

## ----fit_table, include = TRUE, eval=F----------------------------------------
#  fit <- table_fit(res) # model fit table
#  fit[ , c("Name", "LL", "Parameters", "n",
#           "BIC", "Entropy",
#           "prob_min", "prob_max",
#           "n_min", "n_max",
#           "np_ratio", "np_local")]

## ----tabfit, echo = FALSE, eval = eval_results--------------------------------
#  tbl <- fit[ , c("Name", "LL", "Parameters", "n",
#           "BIC", "Entropy",
#           "prob_min", "prob_max",
#           "n_min", "n_max")]
#  names(tbl) <- c("Name", "LL", "p", "n",
#           "BIC", "Entropy",
#           "p_min", "p_max",
#           "n_min", "n_max")
#  knitr::kable(tbl, caption = "Model fit table")

## ----lmr_table, echo = TRUE, eval=FALSE---------------------------------------
#  lr_lmr(res)

## ----tablmr, echo = FALSE, eval = run_everything------------------------------
#  res <- readRDS("res_lpa.RData")
#  tbl <- lr_lmr(res)
#  write.csv(tbl, "lpatablmr.csv", row.names = FALSE)

## ----echo = FALSE, eval = TRUE------------------------------------------------
tbl <- read.csv("lpatablmr.csv", stringsAsFactors = FALSE)
knitr::kable(tbl, caption = "LMR test table", digits = 2)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  library(future)
#  library(progressr)
#  plan(multisession) # Parallel processing for Windows
#  handlers("progress") # Progress bar
#  set.seed(1)
#  res_blrt <- BLRT(res, replications = 20)

## ----eval = run_everything, echo = FALSE--------------------------------------
#  library(future)
#  library(progressr)
#  plan(multisession) # Parallel processing for Windows
#  handlers("progress") # Progress bar
#  set.seed(1)
#  res_blrt <- BLRT(res, replications = 20)
#  write.csv(res_blrt, "appendixbresblrt.csv", row.names = FALSE)

## ----eval = eval_results, echo = FALSE----------------------------------------
#  res_blrt <- read.csv("appendixbresblrt.csv", stringsAsFactors = FALSE)
#  knitr::kable(res_blrt, caption = "BLRT test table", digits = 2)

## ----eval = run_everything, echo = FALSE--------------------------------------
#  res_alt <- mx_profiles(df, classes = 4, variances = "varying")
#  compare <- list(res[[4]], res_alt)
#  fit_compare <- table_fit(compare)
#  write.csv(fit_compare, "lpa_fit_compare.csv", row.names = FALSE)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  res_alt <- mx_profiles(df, classes = 4, variances = "varying")
#  compare <- list(res[[4]], res_alt)
#  table_fit(compare)

## ----tabfitcomp, echo = FALSE, eval = eval_results----------------------------
#  fit_compare <- read.csv("lpa_fit_compare.csv", stringsAsFactors = FALSE)
#  class(fit_compare) <- c("tidy_fit", "data.frame")
#  knitr::kable(fit_compare[ , c("Name", "LL", "Parameters",
#           "BIC", "Entropy",
#           "prob_min", "prob_max",
#           "n_min", "n_max")], caption = "Comparing competing theoretical models")

## ----echo = TRUE, eval =FALSE-------------------------------------------------
#  res_final <- mx_switch_labels(res[[4]])

## ----echo = FALSE, eval = run_everything--------------------------------------
#  res_final <- mx_switch_labels(res[[4]])
#  cp <- class_prob(res_final)
#  out <- list(counts = cp$sum.posterior$proportion)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  table_results(res_final, columns = c("label", "est", "se", "confint", "class"))

## ----echo = FALSE, eval = run_everything--------------------------------------
#  tab <- table_results(res_final, columns = c("label", "est", "se", "confint", "class"))
#  write.csv(tab, "lpa_tab_res.csv", row.names = FALSE)

## ----eval = eval_results, echo=FALSE------------------------------------------
#  tab <- read.csv("lpa_tab_res.csv", stringsAsFactors = FALSE)
#  knitr::kable(tab, caption = "Four-class model results")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  plot_bivariate(res_final)

## ----echo = FALSE, eval = run_everything--------------------------------------
#  p <- plot_bivariate(res_final, return_list = TRUE)
#  p[[1]] <- p[[1]] + scale_y_continuous(breaks= c(0, .1, .2, .3), labels = c(".0", ".1", ".2", ".3"))
#  p <- tidySEM:::merge_corplots(p)
#  ggsave("lpa_bivariate.png", p, device = "png", width = 210, height = 100, units = "mm", scale = 1.5)

## ----echo = FALSE, eval = eval_results, fig.cap="Bivariate profile plot", out.width="80%"----
#  knitr::include_graphics("lpa_bivariate.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  plot_profiles(res_final)

## ----echo = FALSE, eval = run_everything--------------------------------------
#  p1 <- plot_profiles(res_final, ci = NULL, sd= FALSE, add_line = TRUE, rawdata = FALSE)
#  p1 <- p1 + theme(legend.position = c(.85, .2))
#  p2 <- plot_profiles(res_final)
#  p2 <- p2 + theme(legend.position = "none")
#  p <- ggpubr::ggarrange(p1, p2)
#  ggsave("lpa_profiles.png", p, device = "png", width = 210, height = 100, units = "mm", scale = 1)

## ----echo = FALSE, eval = eval_results, fig.cap="Bivariate profile plot", out.width="80%"----
#  knitr::include_graphics("lpa_profiles.png")

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  aux_sex <- BCH(res_final, data = zegwaard_carecompass$sexpatient)

## ----echo = FALSE, eval=run_everything----------------------------------------
#  aux_sex <- BCH(res_final, data = zegwaard_carecompass$sexpatient)
#  saveRDS(aux_sex, "lpa_aux_sex.RData")

## ----echo = FALSE, eval=FALSE-------------------------------------------------
#  aux_sex <- readRDS("lpa_aux_sex.RData")

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  df_aux <- zegwaard_carecompass[, c("freqvisit", "distance")]
#  df_aux$freqvisit <- as.numeric(df_aux$freqvisit)
#  aux_model <- BCH(res_final, model = "freqvisit ~ distance",
#                   data = df_aux)

## ----echo = FALSE, eval=run_everything----------------------------------------
#  df_aux <- zegwaard_carecompass[, c("freqvisit", "distance")]
#  df_aux$freqvisit <- as.numeric(df_aux$freqvisit)
#  aux_model <- BCH(res_final, model = "freqvisit ~ distance",
#                 data = df_aux)
#  saveRDS(aux_model, "lpa_aux_model.RData")

## ----echo = FALSE, eval=FALSE-------------------------------------------------
#  aux_model <- readRDS("lpa_aux_model.RData")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df_new <- data.frame(
#    burdened = 2,
#    trapped = 0.5,
#    negaffect = 1.5,
#    loneliness = 4
#  )
#  predict(res_final, newdata = df_new)

## ----echo = FALSE, eval = TRUE------------------------------------------------
structure(c(0.000808074819286418, 1.35412723878036e-15, 0.999191879285034, 
4.58956785816116e-08, 3), dim = c(1L, 5L), dimnames = list(NULL, 
    c("class1", "class2", "class3", "class4", "predicted")))

