## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", warning = FALSE,
                      message = FALSE, fig.width = 7, fig.height = 4.5,
                      out.width = "100%")
run_openmx <- requireNamespace("OpenMx", quietly = TRUE)
run_everything <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false"))))
show_results <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_SHOW_RESULTS", "false"))))
if (!run_openmx) run_everything <- show_results <- FALSE
# run_cached() manages RDS and PNG files; this path remains for CSV and SVG exports.
input <- knitr::current_input(dir = TRUE)
cache_dir <- if (is.null(input) || !nzchar(input)) {
  file.path(tempdir(), "cache_interactive")
} else {
  file.path(
    dirname(input),
    paste0("cache_", tools::file_path_sans_ext(basename(input)))
  )
}
if (run_everything) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

## ----intro-1------------------------------------------------------------------
library(tidySEM)
library(ggplot2)
df <- zhao_2025

## ----eval = FALSE-------------------------------------------------------------
# # For CSV data:
# df <- read.csv("your_file.csv")
# # For Excel spreadsheets:
# df <- readxl::read_excel("your_file.xlsx", sheet = 1)
# # For SPSS data:
# df <- foreign::read.spss("your_file.sav", to.data.frame = TRUE)
# # For Stata data:
# df <- foreign::read.dta("your_file.dta")

## ----intro-3------------------------------------------------------------------
table_desc <- descriptives(df)
table_desc

## ----intro-4, eval = run_openmx-----------------------------------------------
# This is a character variable, because someone misstyped a t:
v1 <- c("5.1", "4.9", "4.7", "4.6t", "5", "5.4")
class(v1)
# You can turn it into numeric; the mistyped value is lost:
as.numeric(v1)
# This is an integer variable representing a Likert scale:
v2 <- c(2L, 3L, 2L, 5L, 3L, 4L)
class(v2)
# You can turn it into ordered:
ordered(v2)
# This is a factor variable:
v3 <- factor(c("Medicine", "Humanities", "Humanities", "Medicine", "Social science", 
"Medicine"))
class(v3)
# You can encode it as multiple binary dummmies:
tidySEM::mx_dummies(v3)

## ----intro-5------------------------------------------------------------------
df$Gender <- factor(df$Gender,
                    levels = c(1,2),
                    labels = c("M", "F"))

## ----intro-6------------------------------------------------------------------
# First, insert a separator character between scale names and item numbers:
names(df)
names(df) <- gsub("(?=\\d)", "_", names(df), perl = TRUE)
names(df)

## ----intro-7------------------------------------------------------------------
# Make a tidy_sem object
df_tidysem <- tidy_sem(df)
# Construct mean scale scores, removing missing values
df_scales <- score_items(df_tidysem, na.rm = TRUE)
# Add them to df:
df <- data.frame(df, df_scales)

## ----intro-8------------------------------------------------------------------
ord_vars <- grep("\\d$", names(df), value = TRUE)
ord_vars

## ----intro-9------------------------------------------------------------------
df[ord_vars] <- lapply(df[ord_vars], ordered)

## ----eval = FALSE-------------------------------------------------------------
# library(mice)
# df <- ampute(df)$amp

## ----intro-11-----------------------------------------------------------------
df_lpa <- df[, c("AM", "CM")]

## ----intro-code-12, eval=FALSE------------------------------------------------
# res_lpa <- mx_profiles(
#   df_lpa,
#   classes = 1:5,
#   variances = c("equal", "varying"),
#   covariances = c("zero", "equal")
# )

## ----intro-fit-12, echo=FALSE, eval=run_everything || show_results------------
# res_lpa_different <- tidySEM:::run_cached({
# mx_profiles(
#   df_lpa,
#   classes = 1:5,
#   variances = c("equal", "varying"),
#   covariances = c("zero", "equal")
# )
# })

## ----intro-result-13, eval=show_results---------------------------------------
# format_numeric(table_fit(res_lpa_different))

## ----intro-14-----------------------------------------------------------------
k <- 1:10
nrow(df) / (1+3*k)

## ----intro-code-15, eval=FALSE------------------------------------------------
# res_lpa <- mx_profiles(
#   df_lpa,
#   classes = 1:5,
#   variances = "equal",
#   covariances = "zero"
# )

## ----intro-fit-15, echo=FALSE, eval=run_everything || show_results------------
# res_lpa <- tidySEM:::run_cached({
# mx_profiles(
#   df_lpa,
#   classes = 1:5,
#   variances = "equal",
#   covariances = "zero"
# )
# })

## ----eval = show_results------------------------------------------------------
# tab_lpa_fit <- table_fit(res_lpa)
# tab_lpa_fit

## ----eval = show_results------------------------------------------------------
# which.min(tab_lpa_fit$BIC)

## ----intro-code-19, eval=FALSE------------------------------------------------
# res_lpa_blrt <- res_lpa
# res_lpa_blrt[c(4, 5)] <- NULL
# res_blrt <- BLRT(res_lpa_blrt)

## ----intro-fit-19, echo=FALSE, eval=run_everything || show_results------------
# res_blrt <- tidySEM:::run_cached({
# res_lpa_blrt <- res_lpa
# res_lpa_blrt[c(4, 5)] <- NULL
# BLRT(res_lpa_blrt)
# })

## ----intro-result-20, eval=show_results---------------------------------------
# res_blrt

## ----eval = show_results------------------------------------------------------
# sapply(res_lpa, function(i) mx_diagnostics(i)$checks_pass)

## ----eval = show_results------------------------------------------------------
# tab_lpa_fit[, c("Entropy", "prob_min", "prob_max")]

## ----intro-result-23-show, eval=FALSE-----------------------------------------
# # Standard R formatting, some columns get 3 digits, some get 1, but they remain
# # numeric type:
# tab_lpa_fit_formatted <- format(tab_lpa_fit, digits = 2)
# tab_lpa_fit_formatted
# # Using format_numeric; reports 2 digits, but all columns become
# # character type:
# tab_lpa_fit_formatted <- format_numeric(tab_lpa_fit, digits = 2)
# write.csv(tab_lpa_fit_formatted, "tab_lpa_fit.csv", row.names = FALSE)

## ----intro-result-23, eval=show_results, echo = FALSE-------------------------
# # Standard R formatting, some columns get 3 digits, some get 1, but they remain
# # numeric type:
# tab_lpa_fit_formatted <- format(tab_lpa_fit, digits = 2)
# tab_lpa_fit_formatted
# # Using format_numeric; reports 2 digits, but all columns become
# # character type:
# tab_lpa_fit_formatted <- format_numeric(tab_lpa_fit, digits = 2)

## ----eval = show_results------------------------------------------------------
# class_prob(res_lpa[[3]], "sum.posterior")

## ----eval = show_results------------------------------------------------------
# tab_lpa_res <- table_results(res_lpa[[3]])
# tab_lpa_res

## ----intro-plot-code-26, eval=FALSE-------------------------------------------
# plot_profiles(res_lpa[[3]])
# plot_density(res_lpa[[3]])

## ----intro-plot-26-1, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# plot_profiles(res_lpa[[3]])
# }, type = "png", width = 7, height = 4.5, dpi = 150)

## ----intro-plot-26-2, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# plot_density(res_lpa[[3]])
# }, type = "png", width = 7, height = 4.5, dpi = 150)

## ----intro-27-----------------------------------------------------------------
df_lca <- df[, c("AM_1", "AM_2", "AM_3")]

## ----intro-28-----------------------------------------------------------------
k <- 1:10
L <- length(unique(df_lca$AM_1))
M <- ncol(df_lca)
nrow(df) / ((k-1)+k*((L-1)*M))

## ----intro-code-29, eval=FALSE------------------------------------------------
# res_lca <- mx_lca(
#   df_lca,
#   classes = 1:3
# )

## ----intro-fit-29, echo=FALSE, eval=run_everything || show_results------------
# res_lca <- tidySEM:::run_cached({
# mx_lca(
#   df_lca,
#   classes = 1:3
# )
# })
# 

## ----intro-result-30, eval=show_results---------------------------------------
# tab_lca_fit <- table_fit(res_lca)
# tab_lca_fit

## ----eval = show_results------------------------------------------------------
# tab_lca_res <- table_results(res_lca[[2]])
# tab_lca_res

## ----intro-plot-code-32, eval=FALSE-------------------------------------------
# library(ggplot2)
# thresholds <- res_lca[[2]]$class1$Thresholds$result[, "AM_1"]
# ggplot(data.frame(x = c(-4, 4)), aes(x = x)) +
#   stat_function(fun = dnorm, color = "red") +
#   geom_vline(xintercept = thresholds) +
#   scale_x_continuous(breaks = thresholds, labels = 1:6) +
#   theme_bw()

## ----intro-plot-32, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# library(ggplot2)
# thresholds <- res_lca[[2]]$class1$Thresholds$result[, "AM_1"]
# ggplot(data.frame(x = c(-4, 4)), aes(x = x)) +
#   stat_function(fun = dnorm, color = "red") +
#   geom_vline(xintercept = thresholds) +
#   scale_x_continuous(breaks = thresholds, labels = 1:6) +
#   theme_bw()
# }, type = "png", width = 7, height = 4.5, dpi = 150)

## ----eval = show_results------------------------------------------------------
# individual_class <- predict_class(res_lca[[2]])
# table(df$AM_1[individual_class == 1])

## ----eval = show_results------------------------------------------------------
# tab_lca_prob <- table_prob(res_lca[[2]])
# # Let's format it for easier visual interpretation:
# format_numeric(tab_lca_prob, digits = 2)

## ----intro-plot-code-35, eval=FALSE-------------------------------------------
# plot_prob(res_lca[[2]])

## ----intro-plot-35, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# plot_prob(res_lca[[2]])
# }, type = "png", width = 7, height = 4.5, dpi = 150)

## ----intro-36-----------------------------------------------------------------
df_bin <- df_lca
df_bin[] <- lapply(df_bin, function(x) ordered(as.integer(x <= 4)))

## ----intro-37-----------------------------------------------------------------
k <- 1:10
L <- 2
M <- ncol(df_lca)
nrow(df) / ((k-1)+k*((L-1)*M))

## ----intro-code-38, eval=FALSE------------------------------------------------
# res_bin <- mx_lca(
#   df_bin,
#   classes = 1:4
# )

## ----intro-fit-38, echo=FALSE, eval=run_everything || show_results------------
# res_bin <- tidySEM:::run_cached({
# mx_lca(
#   df_bin,
#   classes = 1:4
# )
# })

## ----intro-result-39, eval=show_results---------------------------------------
# tab_bin_fit <- table_fit(res_bin)
# tab_bin_fit

## ----eval = show_results------------------------------------------------------
# tab_bin_prob <- table_prob(res_bin[[2]])
# tab_bin_prob <- tab_bin_prob[!tab_bin_prob$Category == 1, -2]
# # Let's format it for easier visual interpretation:
# format_numeric(tab_bin_prob, digits = 2)

## ----intro-plot-code-41, eval=FALSE-------------------------------------------
# plot_prob(res_bin[[2]])

## ----intro-plot-41, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# plot_prob(res_bin[[2]])
# }, type = "png", width = 7, height = 4.5, dpi = 150)

## ----intro-42-----------------------------------------------------------------
df_mixed <- data.frame(df[, "CM", drop = FALSE], df_bin)

## ----intro-43-----------------------------------------------------------------
k <- 1:10
L <- 2
M <- ncol(df_lca)
nrow(df) / ((k-1)+k*((L-1)*M) + 1 + k)

## ----intro-code-44, eval=FALSE------------------------------------------------
# res_mixed <- mx_mixed_lca(
#   df_mixed,
#   classes = 1:3
# )

## ----intro-fit-44, echo=FALSE, eval=run_everything || show_results------------
# res_mixed <- tidySEM:::run_cached({
# mx_mixed_lca(
#   df_mixed,
#   classes = 1:3
# )
# })

## ----intro-result-45, eval=show_results---------------------------------------
# tab_mixed_fit <- table_fit(res_mixed)
# tab_mixed_fit

## ----eval = show_results------------------------------------------------------
# tab_mixed_res <- table_results(res_mixed[[2]])
# tab_mixed_res

## ----intro-plot-code-47, eval=FALSE-------------------------------------------
# plot_density(res_mixed[[2]], variables = "CM")

## ----intro-plot-47, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# plot_density(res_mixed[[2]], variables = "CM")
# }, type = "png", width = 7, height = 4.5, dpi = 150)

## ----intro-plot-code-48, eval=FALSE-------------------------------------------
# plot_prob(res_mixed[[2]])

## ----intro-plot-48, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# plot_prob(res_mixed[[2]])
# }, type = "png", width = 7, height = 4.5, dpi = 150)

## ----eval = show_results------------------------------------------------------
# cp <- class_prob(res_lpa[[3]])
# format_numeric(head(cp$individual))

## ----eval = show_results------------------------------------------------------
# format_numeric(head(cp$sum.posterior))

## ----eval = show_results------------------------------------------------------
# format_numeric(head(cp$sum.mostlikely))

## ----eval = show_results------------------------------------------------------
# format_numeric(cp$mostlikely.class)

## ----intro-result-53-show, eval=FALSE-----------------------------------------
# tab <- table_results(res_lpa[[2]], columns = NULL)
# tab <- tab[tab$Category == "Means", c("lhs", "class", "est", "se")]
# tab$Mean <- paste0(tab$est, "(", tab$se, ")")
# tab <- reshape(tab[, c("class", "lhs", "Mean")], direction = "wide", idvar = "class", timevar = "lhs", sep = " ")
# write.csv(tab, "tab_res.csv", row.names = FALSE)
# tab

## ----intro-result-53, eval=show_results, echo = FALSE-------------------------
# tab <- table_results(res_lpa[[2]], columns = NULL)
# tab <- tab[tab$Category == "Means", c("lhs", "class", "est", "se")]
# tab$Mean <- paste0(tab$est, "(", tab$se, ")")
# tab <- reshape(tab[, c("class", "lhs", "Mean")], direction = "wide", idvar = "class", timevar = "lhs", sep = " ")
# tab

## ----intro-plot-code-54, eval=FALSE-------------------------------------------
# p <- plot_density(res_lpa[[2]])
# ggsave(filename = "figure1.svg", plot = p, width = 148, height = 105, units = "mm", dpi = 300)
# p

## ----intro-plot-54, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# p <- plot_density(res_lpa[[2]])
# p
# }, type = "png", width = 7, height = 4.5, dpi = 150)

