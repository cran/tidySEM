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
  tidy.opts = list(width.cutoff = 60), tidy = requireNamespace("formatR", quietly = TRUE)
)
options(scipen = 1, digits = 2)
run_everything <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false"))))
show_results <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_SHOW_RESULTS", "false"))))
if (!run_openmx) run_everything <- show_results <- FALSE
if (show_results) knitr::opts_chunk$set(results = "markup")

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

## ----prepare-plot-data, include = FALSE---------------------------------------
df_plot <- df
names(df_plot) <- paste0("Value.", names(df_plot))
df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
                   timevar = "Variable")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# df_plot <- df
# names(df_plot) <- paste0("Value.", names(df_plot))
# df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                    timevar = "Variable")
# ggplot(df_plot, aes(x = Value)) +
#   geom_density() +
#   facet_wrap(~Variable)+
#   theme_bw()

## ----image-desc, echo = FALSE, eval = run_everything || show_results, out.width = "100%"----
# tidySEM:::run_cached({
# df_plot <- df
# names(df_plot) <- paste0("Value.", names(df_plot))
# df_plot <- reshape(df_plot, varying = names(df_plot), direction = "long",
#                    timevar = "Variable")
# ggplot(df_plot, aes(x = Value)) +
#   geom_density() +
#   facet_wrap(~Variable)+
#   theme_bw()
# }, type = "png", width = 100/25.4, height = 100/25.4, dpi = 150)

## ----prepare-log-data, include = FALSE----------------------------------------
df_plot$Value <- log(df_plot$Value)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# df_plot$Value <- log(df_plot$Value)
# ggplot(df_plot, aes(x = Value)) +
#   geom_density() +
#   facet_wrap(~Variable)+
#   theme_bw()

## ----image-desc-log, echo = FALSE, eval = run_everything || show_results, out.width = "100%"----
# tidySEM:::run_cached({
# ggplot(df_plot, aes(x = Value)) +
#   geom_density() +
#   facet_wrap(~Variable)+
#   theme_bw()
# }, type = "png", width = 100/25.4, height = 100/25.4, dpi = 150)

## ----prepare-analysis-data, include = FALSE-----------------------------------
df <- reshape(df_plot, direction = "wide", v.names = "Value")[, -1]
names(df) <- gsub("Value.", "", names(df), fixed = TRUE)

## ----eval = FALSE, echo = TRUE------------------------------------------------
# df <- reshape(df_plot, direction = "wide", v.names = "Value")[, -1]
# names(df) <- gsub("Value.", "", names(df), fixed = TRUE)
# ggplot(df, aes(x = length, y = width)) +
#   geom_point(alpha = .1) +
#   theme_bw()

## ----image-scatter, echo = FALSE, eval = run_everything || show_results, out.width = "100%"----
# tidySEM:::run_cached({
# ggplot(df, aes(x = length, y = width)) +
#   geom_point(alpha = .1) +
#   theme_bw()
# }, type = "png", width = 100/25.4, height = 100/25.4, dpi = 150)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
# pca <- princomp(df)
# df <- data.frame(pca$scores)
# names(df) <- c("pc1", "pc2")

## ----fitlca, eval = FALSE, echo = TRUE----------------------------------------
# set.seed(123)
# res <- mx_profiles(data = df,
#                    classes = 1:3,
#                    variances = c("equal", "varying"),
#                    covariances = c("equal", "varying"),
#                    expand_grid = TRUE)

## ----fitlca-cached, eval = run_everything || show_results, echo = FALSE-------
# res <- tidySEM:::run_cached({
#   set.seed(123)
#   mx_profiles(data = df,
#               classes = 1:3,
#               variances = c("equal", "varying"),
#               covariances = c("equal", "varying"),
#               expand_grid = TRUE)
# })

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# fit <- table_fit(res)

## ----fit-table-compute, eval = run_everything || show_results, echo = FALSE----
# fit <- table_fit(res)

## ----tabfit, echo = TRUE, eval = show_results---------------------------------
# tbl <- fit[ , c("Name", "LL", "Parameters",
#          "BIC", "Entropy",
#          "prob_min",
#          "n_min",
#          "np_ratio", "np_local"
#          )]
# names(tbl) <- c("Name", "LL", "p", "BIC", "Ent.",
#          "p_min",
#          "n_min",
#          "np_ratio", "np_local")
# knitr::kable(tbl,
#                   caption = "Model fit table.")

## ----echo = TRUE, eval = run_everything || show_results-----------------------
# fit <- fit[!fit$n_min < .1, ]

## ----echo = TRUE, eval = FALSE------------------------------------------------
# plot(fit) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## ----image-fit, echo = FALSE, eval = run_everything || show_results, out.width = "100%", fig.cap="Bivariate profile plot"----
# tidySEM:::run_cached({
# plot(fit) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# }, type = "png", width = 100/25.4, height = 100/25.4, dpi = 150)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# res_bic <- res[["free var, free cov 2"]]
# cp <- class_prob(res_bic)
# results <- table_results(res_bic, columns = c("label", "est", "std_est"))
# results

## ----results-compute, echo = FALSE, eval = run_everything || show_results-----
# res_bic <- res[["free var, free cov 2"]]
# cp <- class_prob(res_bic)
# results <- table_results(res_bic, columns = c("label", "est", "est_std"))

## ----results-table, eval = show_results, echo = FALSE-------------------------
# knitr::kable(results, digits = 2, caption = "Results of a 2-class model with free (co)variances")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# plot_bivariate(res_bic)

## ----image-bivariate, echo = FALSE, eval = run_everything || show_results, out.width = "100%", fig.cap="Bivariate profile plot"----
# tidySEM:::run_cached({
# plot_bivariate(res_bic)
# }, type = "png", width = 315/25.4, height = 150/25.4, dpi = 150)

