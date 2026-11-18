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

## ----echo = TRUE--------------------------------------------------------------
library(tidySEM)
library(ggplot2)
library(MASS)

## -----------------------------------------------------------------------------
# Get descriptives
df <- plas_depression
desc <- descriptives(df)
desc <- desc[, c("name", "mean", "median", "sd", "min", "max", 
"skew_2se", "kurt_2se")
]
knitr::kable(desc, caption = "Item descriptives")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# df_plot <- reshape(df, direction = "long", varying = names(df))
# ggplot(df_plot, aes(x = scl)) +
#   geom_density() +
#   facet_wrap(~time) + theme_bw()

## ----lcga-plot-dist, echo = FALSE, eval = run_everything || show_results, out.width = "100%"----
# tidySEM:::run_cached({
# df_plot <- reshape(df, direction = "long", varying = names(df))
# p = ggplot(df_plot, aes(x = scl)) +
#   geom_density() +
#   facet_wrap(~time) + theme_bw()
# p
# }, type = "png", width = 150 / 25.4, height = 120 / 25.4, dpi = 150)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# df_scores <- df_plot
# # Store original range of SCL
# rng_scl <- range(df_scores$scl)
# # Log-transform
# df_scores$log <- scales::rescale(log(df_scores$scl), to = c(0, 1))
# # Square root transform
# df_scores$sqrt <- scales::rescale(sqrt(df_scores$scl), to = c(0, 1))
# # Cube root transform
# df_scores$qrt <- scales::rescale(df_scores$scl^.33, to = c(0, 1))
# # Reciprocal transform
# df_scores$reciprocal <- scales::rescale(1/df_scores$scl, to = c(0, 1))
# # Define function for Box-Cox transformation
# bc <- function(x, lambda){
#   (((x ^ lambda) - 1) / lambda)
# }
# # Inverse Box-Cox transformation
# invbc <- function(x, lambda){
#   ((x*lambda)+1)^(1/lambda)
# }
# # Box-Cox transform
# fit_bc <- lm(scl ~ 1, data = df_scores, y = TRUE)
# b <- MASS::boxcox(fit_bc, plotit = FALSE)
# lambda <- b$x[which.max(b$y)]
# df_scores$boxcox <- bc(df_scores$scl, lambda)
# # Store range of Box-Cox transformed data
# rng_bc <- range(df_scores$boxcox)
# df_scores$boxcox <- scales::rescale(df_scores$boxcox, to = c(0, 1))
# # Rescale SCL
# df_scores$scl <- scales::rescale(df_scores$scl, to = c(0, 1))

## ----eval = FALSE, echo = TRUE------------------------------------------------
# # Make plot data
# df_plot <- do.call(rbind, lapply(c("scl", "log", "sqrt", "qrt", "boxcox"), function(n){
#   data.frame(df_scores[c("time", "id")],
#              Value = df_scores[[n]],
#              Transformation = n)
# }))
# # Plot
# ggplot(df_plot, aes(x = Value, colour = Transformation)) +
#   geom_density() +
#   facet_wrap(~time) +
#   scale_y_sqrt() +
#   xlab("scl (rescaled to 0-1)") +
#   theme_bw()

## ----lcga-transform-data, echo=FALSE, include=FALSE---------------------------
df_plot <- reshape(df, direction = "long", varying = names(df))
df_scores <- df_plot

# Store original range of SCL
rng_scl <- range(df_scores$scl)

# Log-transform
df_scores$log <- scales::rescale(log(df_scores$scl), to = c(0, 1))

# Square root transform
df_scores$sqrt <- scales::rescale(sqrt(df_scores$scl), to = c(0, 1))

# Cube root transform
df_scores$qrt <- scales::rescale(df_scores$scl^.33, to = c(0, 1))

# Reciprocal transform
df_scores$reciprocal <- scales::rescale(1 / df_scores$scl, to = c(0, 1))

# Define Box-Cox functions
bc <- function(x, lambda) {
  (((x ^ lambda) - 1) / lambda)
}

invbc <- function(x, lambda) {
  ((x * lambda) + 1)^(1 / lambda)
}

# Box-Cox transformation
fit_bc <- lm(scl ~ 1, data = df_scores, y = TRUE)
b <- MASS::boxcox(fit_bc, plotit = FALSE)

lambda <- b$x[which.max(b$y)]

df_scores$boxcox <- bc(df_scores$scl, lambda)

# Store original Box-Cox range
rng_bc <- range(df_scores$boxcox)

df_scores$boxcox <- scales::rescale(
  df_scores$boxcox,
  to = c(0, 1)
)

# Rescale SCL
df_scores$scl <- scales::rescale(
  df_scores$scl,
  to = c(0, 1)
)

## ----lcga-plot-trans, echo = FALSE, eval = run_everything || show_results, out.width = "100%"----
# tidySEM:::run_cached({
# df_plot <- do.call(rbind, lapply(c("scl", "log", "sqrt", "qrt", "boxcox"), function(n){
#   data.frame(df_scores[c("time", "id")],
#              Value = df_scores[[n]],
#              Transformation = n)
# }))
# # Plot
# ggplot(df_plot, aes(x = Value, colour = Transformation)) +
#   geom_density() +
#   facet_wrap(~time) +
#   scale_y_sqrt() +
#   xlab("scl (rescaled to 0-1)") +
#   theme_bw()
# }, type = "png", width = 150 / 25.4, height = 120 / 25.4, dpi = 150)

## ----lcga-data-wide, echo = FALSE, include = FALSE----------------------------
dat <- df_scores[, c("id", "time", "boxcox")]
dat <- reshape(dat, direction = "wide", v.names = "boxcox", timevar = "time", idvar = "id")
names(dat) <- gsub("boxcox.", "scl", names(dat))

## ----lcga, echo = TRUE, eval = FALSE------------------------------------------
# set.seed(27796)
# dat[["id"]] <- NULL
# res_step <- mx_growth_mixture(
#   model =
#   "
#   i =~ 1*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
#   step =~ 0*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
#   s =~ 0*scl1 + 0*scl2 + 1*scl3 +2*scl4 +3*scl5 +4*scl6
#   scl1 ~~ vscl1*scl1
#   scl2 ~~ vscl2*scl2
#   scl3 ~~ vscl3*scl3
#   scl4 ~~ vscl4*scl4
#   scl5 ~~ vscl5*scl5
#   scl6 ~~ vscl6*scl6
#   i ~~ 0*i
#   step ~~ 0*step
#   s ~~ 0*s
#   i ~~ 0*s
#   i ~~ 0*step
#   s ~~ 0*step
#   scl1~0*1
#   scl2~0*1
#   scl3~0*1
#   scl4~0*1
#   scl5~0*1
#   scl6~0*1
#   i~NA*1
#   s~NA*1
#   step~NA*1
#   ", classes = 1:4,
#   data = dat)

## ----lcga_cached, echo = FALSE, eval = (run_everything || show_results)-------
# res_step <- tidySEM:::run_cached({
# set.seed(27796)
# dat[["id"]] <- NULL
# res_step <- mx_growth_mixture(
#   model =
# "i =~ 1*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
# step =~ 0*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
# s =~ 0*scl1 + 0*scl2 + 1*scl3 +2*scl4 +3*scl5 +4*scl6
# scl1 ~~ vscl1*scl1
# scl2 ~~ vscl2*scl2
# scl3 ~~ vscl3*scl3
# scl4 ~~ vscl4*scl4
# scl5 ~~ vscl5*scl5
# scl6 ~~ vscl6*scl6
# i ~~ 0*i
# step ~~ 0*step
# s ~~ 0*s
# i ~~ 0*s
# i ~~ 0*step
# s ~~ 0*step
# scl1~0*1
#   scl2~0*1
#   scl3~0*1
#   scl4~0*1
#   scl5~0*1
#   scl6~0*1
#   i~NA*1
#   s~NA*1
#   step~NA*1", classes = 1:4,
#   data = dat)
# res_step
# })

## ----eval = FALSE, echo = TRUE------------------------------------------------
# # Get fit table fit
# tab_fit <- table_fit(res_step)
# # Select columns
# tab_fit[, c("Name", "Classes", "LL", "Parameters", "BIC", "Entropy", "prob_min", "n_min", "warning", "lmr_p")]

## ----tabfit, eval = (run_everything || show_results), echo = FALSE------------
# tab_fit <- tidySEM:::run_cached({
#   table_fit(res_step)
# })
# knitr::kable(tab_fit[, c("Name", "Classes", "LL", "Parameters", "BIC", "Entropy", "prob_min", "n_min")], digits = 2, caption = "Fit of LCGA models")

## ----plotscree, echo = TRUE, eval = FALSE-------------------------------------
# plot(tab_fit, statistics = c("AIC", "BIC", "saBIC"))

## ----lcga-plot-fit, echo = FALSE, eval = run_everything || show_results, out.width = "100%"----
# tidySEM:::run_cached({
# p <- plot(tab_fit, statistics = c("AIC", "BIC", "saBIC"))
# p
# }, type = "png", width = 150 / 25.4, height = 120 / 25.4, dpi = 150)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# res_final <- mx_switch_labels(res_step[[3]], param = "M[1,7]", decreasing = FALSE)
# tab_res <- table_results(res_final, columns = NULL)
# # Select rows and columns
# tab_res <- tab_res[tab_res$Category %in% c("Means", "Variances"), c("Category", "lhs", "est", "se", "pval", "confint", "name")]
# tab_res

## ----tabres_cached, echo = FALSE, eval = (run_everything || show_results)-----
# tab_res <- tidySEM:::run_cached({
# res_final <- mx_switch_labels(res_step[[3]], param = "M[1,7]", decreasing = FALSE)
# table_results(res_final, columns = NULL)
# })

## ----tabres, echo = FALSE, eval = show_results--------------------------------
# tab_res <- tab_res[tab_res$Category %in% c("Means", "Variances"), c("Category", "lhs", "est", "se", "pval", "confint", "name")]
# knitr::kable(tab_res, digits = 2, caption = "Results from 3-class LCGA model", longtable = TRUE)

## ----params, echo = TRUE, eval = FALSE----------------------------------------
# names(coef(res_final))

## ----echo = FALSE, eval = TRUE------------------------------------------------
c("mix3.weights[1,2]", "mix3.weights[1,3]", "vscl1", "vscl2", 
"vscl3", "vscl4", "vscl5", "vscl6", "class1.M[1,7]", "class1.M[1,8]", 
"class1.M[1,9]", "class2.M[1,7]", "class2.M[1,8]", "class2.M[1,9]", 
"class3.M[1,7]", "class3.M[1,8]", "class3.M[1,9]")

## ----echo = TRUE, eval = FALSE------------------------------------------------
# wald_tests <- wald_test(res_final,
#                    "
#                    class1.M[1,7] = class2.M[1,7]&
#                    class1.M[1,7] = class3.M[1,7];
#                    class1.M[1,8] = class2.M[1,8]&
#                    class1.M[1,8] = class3.M[1,8];
#                    class1.M[1,9] = class2.M[1,9]&
#                    class1.M[1,9] = class3.M[1,9]")
# # Rename the hypothesis
# wald_tests$Hypothesis <- c("Mean i", "Mean step", "Mean slope")
# knitr::kable(wald_tests, digits = 2, caption = "Wald tests")

## ----wald_cached, echo = FALSE, eval = (run_everything || show_results)-------
# wald_tests <- tidySEM:::run_cached({
# res_final <- mx_switch_labels(res_step[[3]], param = "M[1,7]", decreasing = FALSE)
# wald_tests <- wald_test(res_final,
#                    "class1.M[1,7] = class2.M[1,7]&class1.M[1,7] = class3.M[1,7];class1.M[1,8] = class2.M[1,8]&class1.M[1,8] = class3.M[1,8];class1.M[1,9] = class2.M[1,9]&class1.M[1,9] = class3.M[1,9]")
# wald_tests$Hypothesis <- c("Mean i", "Mean step", "Mean slope")
# wald_tests
# })

## ----waldtests, echo = FALSE, eval = show_results-----------------------------
# knitr::kable(wald_tests, digits = 2, caption = "Wald tests")

## ----makelcgaplot, echo = TRUE, eval = FALSE----------------------------------
# p <- plot_growth(res_step[[3]], rawdata = TRUE, alpha_range = c(0, .05))
# # Add Y-axis breaks in original scale
# brks <- seq(0, 1, length.out = 5)
# labs <- round(invbc(scales::rescale(brks, from = c(0, 1), to = rng_bc), lambda))
# p <- p + scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = labs) + ylab("SCL (rescaled from Box-Cox)")
# p

## ----lcga-plot-traj, echo = FALSE, eval = run_everything || show_results, out.width = "100%"----
# tidySEM:::run_cached({
# p <- plot_growth(res_step[[3]], rawdata = TRUE, alpha_range = c(0, .05))
# # Add Y-axis breaks in original scale
# brks <- seq(0, 1, length.out = 5)
# labs <- round(invbc(scales::rescale(brks, from = c(0, 1), to = rng_bc), lambda))
# p <- p + scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = labs) + ylab("SCL (rescaled from Box-Cox)")
# p
# }, type = "png", width = 150 / 25.4, height = 120 / 25.4, dpi = 150)

