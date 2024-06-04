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
df_analyze <- 
  alkema_microplastics[alkema_microplastics$category == "Fragment", ]
df <- df_analyze[ ,c("length", "width")]

## ----tabdesc, echo = TRUE, eval=TRUE------------------------------------------
desc <- tidySEM::descriptives(df)
desc <- desc[, c("name", "type", "n", "unique", 
"mean", "median", "sd", "min", "max", "skew_2se", "kurt_2se")]
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
#  df_plot <- df
#  names(df_plot) <- paste0("Value.", names(df_plot))
#  df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                     timevar = "Variable")
#  p <- ggplot(df_plot, aes(x = Value)) +
#    geom_density() +
#    facet_wrap(~Variable)+
#    theme_bw()
#  ggsave("plot_gmm_desc.png", p, device = "png", width = 100, height = 100, units = "mm")

## ----figdesc, echo = FALSE, eval = eval_results-------------------------------
#  df_plot <- df
#  names(df_plot) <- paste0("Value.", names(df_plot))
#  df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                     timevar = "Variable")
#  knitr::include_graphics("plot_gmm_desc.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df_plot$Value <- log(df_plot$Value)
#  ggplot(df_plot, aes(x = Value)) +
#    geom_density() +
#    facet_wrap(~Variable)+
#    theme_bw()

## ----echo = FALSE, eval = run_everything--------------------------------------
#  df_plot$Value <- log(df_plot$Value)
#  p <- ggplot(df_plot, aes(x = Value)) +
#    geom_density() +
#    facet_wrap(~Variable)+
#    theme_bw()
#  ggsave("plot_gmm_desc_log.png", p, device = "png", width = 100, height = 100, units = "mm")

## ----figdesclog, echo = FALSE, eval = eval_results----------------------------
#  knitr::include_graphics("plot_gmm_desc_log.png")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  df <- reshape(df_plot, direction = "wide", v.names = "Value")[, -1]
#  names(df) <- gsub("Value.", "", names(df), fixed = TRUE)
#  ggplot(df, aes(x = length, y = width)) +
#    geom_point(alpha = .1) +
#    theme_bw()

## ----eval = run_everything, echo = FALSE--------------------------------------
#  df <- reshape(df_plot, direction = "wide", v.names = "Value")[, -1]
#  names(df) <- gsub("Value.", "", names(df), fixed = TRUE)
#  p <- ggplot(df, aes(x = length, y = width)) +
#    geom_point(alpha = .1) +
#    theme_bw()
#  ggsave("plot_gmm_scatter.png", p, device = "png", width = 100, height = 100, units = "mm")

## ----figscatter, eval = eval_results, echo = FALSE----------------------------
#  knitr::include_graphics("plot_gmm_scatter.png")

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  pca <- princomp(df)
#  df <- data.frame(pca$scores)
#  names(df) <- c("pc1", "pc2")

## ----fitlca, eval = FALSE, echo = TRUE----------------------------------------
#  set.seed(123)
#  res <- mx_profiles(data = df,
#                     classes = 1:3,
#                     variances = c("equal", "varying"),
#                     covariances = c("equal", "varying"),
#                     expand_grid = TRUE)
#  saveRDS(res, "res_gmm.RData")

## ----eval = run_everything, echo = FALSE--------------------------------------
#  set.seed(123)
#  res <- mx_profiles(data = df,
#                     classes = 1:3,
#                     variances = c("equal", "varying"),
#                     covariances = c("equal", "varying"),
#                     expand_grid = TRUE)
#  saveRDS(res, "res_gmm.RData")
#  fit <- table_fit(res)
#  write.csv(fit, "gmm_tabfit.csv", row.names = FALSE)
#  #"Warning: In model 'mix4' Optimizer returned a non-zero status code 6. The model does not satisfy the first-order optimality conditions to the required accuracy, and no improved point for the merit function could be found during the final linesearch (Mx status RED)"

## ----eval = eval_results, echo = FALSE----------------------------------------
#  fit <- read.csv("gmm_tabfit.csv", stringsAsFactors = FALSE)
#  class(fit) <- c("tidy_fit", "data.frame")

## ----echo = TRUE, eval=F------------------------------------------------------
#  fit <- table_fit(res)

## ----tabfit, echo = TRUE, eval = eval_results---------------------------------
#  tbl <- fit[ , c("Name", "LL", "Parameters",
#           "BIC", "Entropy",
#           "prob_min",
#           "n_min",
#           "np_ratio", "np_local"
#           )]
#  names(tbl) <- c("Name", "LL", "p", "BIC", "Ent.",
#           "p_min",
#           "n_min",
#           "np_ratio", "np_local")
#  knitr::kable(tbl,
#                    caption = "Model fit table.")

## ----echo = TRUE, eval = run_everything---------------------------------------
#  fit <- fit[!fit$n_min < .1, ]

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  plot(fit) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## ----echo = FALSE, eval = run_everything--------------------------------------
#  p <- plot(fit) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  ggsave("gmm_plotfit.png", p, device = "png", width = 100, height = 100, units = "mm", scale = 1)

## ----echo = FALSE, eval = eval_results, fig.cap="Bivariate profile plot", out.width="100%"----
#  knitr::include_graphics("gmm_plotfit.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  res_bic <- res[["free var, free cov 2"]]
#  cp <- class_prob(res_bic)
#  results <- table_results(res_bic, columns = c("label", "est", "std_est"))
#  results

## ----echo = FALSE, eval = run_everything--------------------------------------
#  res_bic <- res[["free var, free cov 2"]]
#  cp <- class_prob(res_bic)
#  results <- table_results(res_bic, columns = c("label", "est", "est_std"))
#  write.csv(results, "gmm_results.csv", row.names = FALSE)

## ----eval = eval_results, echo = FALSE----------------------------------------
#  results <- read.csv("gmm_results.csv", stringsAsFactors = FALSE)
#  knitr::kable(results, digits = 2, caption = "Results of a 2-class model with free (co)variances")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  plot_bivariate(res_bic)

## ----echo = FALSE, eval = run_everything--------------------------------------
#  p <- plot_bivariate(res_bic)
#  ggsave("gmm_bivariate_bic.png", p, device = "png", width = 210, height = 100, units = "mm", scale = 1.5)

## ----echo = FALSE, eval = eval_results, fig.cap="Bivariate profile plot", out.width="100%"----
#  knitr::include_graphics("gmm_bivariate_bic.png")

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  df_pt <- mx_dummies(df_analyze$poly_type)
#  aux_pt <- BCH(res_bic, model = "poly_typeOther | t1
#                                  poly_typePE | t1
#                                  poly_typePP | t1", data = df_pt)
#  aux_pt <- mxTryHardOrdinal(aux_pt)

## ----echo = FALSE, eval=run_everything----------------------------------------
#  df_pt <- mx_dummies(df_analyze$poly_type)
#  aux_pt <- BCH(res_bic, model = "poly_typeOther | t1
#                                    poly_typePE | t1
#                                    poly_typePP | t1", data = df_pt)
#  aux_pt <- mxTryHardOrdinal(aux_pt)
#  saveRDS(aux_pt, "gmm_aux_pt.RData")

## ----echo = FALSE, eval=FALSE-------------------------------------------------
#  aux_pt <- readRDS("gmm_aux_pt.RData")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  wald_test(aux_pt, "class1.Thresholds[1,1] = class2.Thresholds[1,1];
#            class1.Thresholds[1,2] = class2.Thresholds[1,2];
#            class1.Thresholds[1,3] = class2.Thresholds[1,3]")

