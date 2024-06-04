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

## ----echo = TRUE, eval = TRUE-------------------------------------------------
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
#  df_plot <- reshape(df, direction = "long", varying = names(df))
#  ggplot(df_plot, aes(x = scl)) +
#    geom_density() +
#    facet_wrap(~time) + theme_bw()

## ----eval = run_everything, echo = FALSE--------------------------------------
#  df_plot <- reshape(df, direction = "long", varying = names(df))
#  p = ggplot(df_plot, aes(x = scl)) +
#    geom_density() +
#    facet_wrap(~time) + theme_bw()
#  ggsave("plot_dist.png", p, width = 150, height = 120, units = "mm")

## ----echo = FALSE, out.width="80%", eval = eval_results-----------------------
#  knitr::include_graphics("plot_dist.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df_scores <- df_plot
#  # Store original range of SCL
#  rng_scl <- range(df_scores$scl)
#  # Log-transform
#  df_scores$log <- scales::rescale(log(df_scores$scl), to = c(0, 1))
#  # Square root transform
#  df_scores$sqrt <- scales::rescale(sqrt(df_scores$scl), to = c(0, 1))
#  # Cube root transform
#  df_scores$qrt <- scales::rescale(df_scores$scl^.33, to = c(0, 1))
#  # Reciprocal transform
#  df_scores$reciprocal <- scales::rescale(1/df_scores$scl, to = c(0, 1))
#  # Define function for Box-Cox transformation
#  bc <- function(x, lambda){
#    (((x ^ lambda) - 1) / lambda)
#  }
#  # Inverse Box-Cox transformation
#  invbc <- function(x, lambda){
#    ((x*lambda)+1)^(1/lambda)
#  }
#  # Box-Cox transform
#  b <- MASS::boxcox(lm(df_scores$scl ~ 1), plotit = FALSE)
#  lambda <- b$x[which.max(b$y)]
#  df_scores$boxcox <- bc(df_scores$scl, lambda)
#  # Store range of Box-Cox transformed data
#  rng_bc <- range(df_scores$boxcox)
#  df_scores$boxcox <- scales::rescale(df_scores$boxcox, to = c(0, 1))
#  # Rescale SCL
#  df_scores$scl <- scales::rescale(df_scores$scl, to = c(0, 1))

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  # Make plot data
#  df_plot <- do.call(rbind, lapply(c("scl", "log", "sqrt", "qrt", "boxcox"), function(n){
#    data.frame(df_scores[c("time", "id")],
#               Value = df_scores[[n]],
#               Transformation = n)
#  }))
#  # Plot
#  ggplot(df_plot, aes(x = Value, colour = Transformation)) +
#    geom_density() +
#    facet_wrap(~time) +
#    scale_y_sqrt() +
#    xlab("scl (rescaled to 0-1)") +
#    theme_bw()

## ----eval = run_everything, echo = FALSE--------------------------------------
#  df_plot <- reshape(df, direction = "long", varying = names(df))
#  df_scores <- df_plot
#  # Store original range of SCL
#  rng_scl <- range(df_scores$scl)
#  # Log-transform
#  df_scores$log <- scales::rescale(log(df_scores$scl), to = c(0, 1))
#  # Square root transform
#  df_scores$sqrt <- scales::rescale(sqrt(df_scores$scl), to = c(0, 1))
#  # Cube root transform
#  df_scores$qrt <- scales::rescale(df_scores$scl^.33, to = c(0, 1))
#  # Reciprocal transform
#  df_scores$reciprocal <- scales::rescale(1/df_scores$scl, to = c(0, 1))
#  # Define function for Box-Cox transformation
#  bc <- function(x, lambda){
#    (((x ^ lambda) - 1) / lambda)
#  }
#  # Inverse Box-Cox transformation
#  invbc <- function(x, lambda){
#    ((x*lambda)+1)^(1/lambda)
#  }
#  # Box-Cox transform
#  b <- MASS::boxcox(lm(df_scores$scl ~ 1), plotit = FALSE)
#  lambda <- b$x[which.max(b$y)]
#  df_scores$boxcox <- bc(df_scores$scl, lambda)
#  # Store range of Box-Cox transformed data
#  rng_bc <- range(df_scores$boxcox)
#  df_scores$boxcox <- scales::rescale(df_scores$boxcox, to = c(0, 1))
#  # Rescale SCL
#  df_scores$scl <- scales::rescale(df_scores$scl, to = c(0, 1))
#  
#  df_plot <- do.call(rbind, lapply(c("scl", "log", "sqrt", "qrt", "boxcox"), function(n){
#    data.frame(df_scores[c("time", "id")],
#               Value = df_scores[[n]],
#               Transformation = n)
#  }))
#  # Plot
#  p <- ggplot(df_plot, aes(x = Value, colour = Transformation)) +
#    geom_density() +
#    facet_wrap(~time) +
#    scale_y_sqrt() +
#    xlab("scl (rescaled to 0-1)") +
#    theme_bw()
#  ggsave("plot_trans.png", p, width = 150, height = 120, units = "mm")

## ----echo = FALSE, eval = eval_results, out.width="80%"-----------------------
#  knitr::include_graphics("plot_trans.png")

## ----echo = TRUE, eval = run_everything---------------------------------------
#  dat <- df_scores[, c("id", "time", "boxcox")]
#  dat <- reshape(dat, direction = "wide", v.names = "boxcox", timevar = "time", idvar = "id")
#  names(dat) <- gsub("boxcox.", "scl", names(dat))

## ----lcga, echo = TRUE, eval = FALSE------------------------------------------
#  set.seed(27796)
#  dat[["id"]] <- NULL
#  res_step <- mx_growth_mixture(
#    model =
#    "
#    i =~ 1*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
#    step =~ 0*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
#    s =~ 0*scl1 + 0*scl2 + 1*scl3 +2*scl4 +3*scl5 +4*scl6
#    scl1 ~~ vscl1*scl1
#    scl2 ~~ vscl2*scl2
#    scl3 ~~ vscl3*scl3
#    scl4 ~~ vscl4*scl4
#    scl5 ~~ vscl5*scl5
#    scl6 ~~ vscl6*scl6
#    i ~~ 0*i
#    step ~~ 0*step
#    s ~~ 0*s
#    i ~~ 0*s
#    i ~~ 0*step
#    s ~~ 0*step", classes = 1:5,
#    data = dat)
#  # Additional iterations because of
#  # convergence problems for model 1:
#  res_step[[1]] <- mxTryHardWideSearch(res_step[[1]], extraTries = 50)

## ----echo = FALSE, eval = run_everything--------------------------------------
#  set.seed(27796)
#  dat[["id"]] <- NULL
#  res_step <- mx_growth_mixture(
#    model =
#  "i =~ 1*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
#  step =~ 0*scl1 + 1*scl2 + 1*scl3 +1*scl4 +1*scl5 +1*scl6
#  s =~ 0*scl1 + 0*scl2 + 1*scl3 +2*scl4 +3*scl5 +4*scl6
#  scl1 ~~ vscl1*scl1
#  scl2 ~~ vscl2*scl2
#  scl3 ~~ vscl3*scl3
#  scl4 ~~ vscl4*scl4
#  scl5 ~~ vscl5*scl5
#  scl6 ~~ vscl6*scl6
#  i ~~ 0*i
#  step ~~ 0*step
#  s ~~ 0*s
#  i ~~ 0*s
#  i ~~ 0*step
#  s ~~ 0*step", classes = 1:5,
#    data = dat)
#  # Additional iterations because of
#  # convergence problems for model 1:
#  res_step[[1]] <- mxTryHardWideSearch(res_step[[1]], extraTries = 50)
#  saveRDS(res_step, "res_step.RData")
#  fit <- table_fit(res_step)
#  write.csv(fit, "lcga_fit.csv", row.names = FALSE)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  res_step <- readRDS("res_step.RData")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  # Get fit table fit
#  tab_fit <- table_fit(res_step)
#  # Select columns
#  tab_fit[, c("Name", "Classes", "LL", "Parameters", "BIC", "Entropy", "prob_min", "n_min", "warning", "lmr_p")]

## ----tabfit, eval = eval_results, echo = FALSE--------------------------------
#  # Get fit table fit
#  tab_fit <- read.csv("lcga_fit.csv", stringsAsFactors = FALSE)
#  class(tab_fit) <- c("tidy_fit", "data.frame")
#  knitr::kable(tab_fit[, c("Name", "Classes", "LL", "Parameters", "BIC", "Entropy", "prob_min", "n_min")], digits = 2, caption = "Fit of LCGA models")

## ----plotscree, echo = TRUE, eval = FALSE-------------------------------------
#  plot(tab_fit, statistics = c("AIC", "BIC", "saBIC"))

## ----echo = FALSE, eval = run_everything--------------------------------------
#  p <- plot(tab_fit, statistics = c("AIC", "BIC", "saBIC"))
#  ggsave("lcga_plot_fit.png", p, width = 150, height = 120, units = "mm")

## ----echo = FALSE, eval = eval_results, out.width="80%"-----------------------
#  knitr::include_graphics("lcga_plot_fit.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  res_final <- mx_switch_labels(res_step[[3]], param = "M[1,7]", decreasing = FALSE)
#  tab_res <- table_results(res_final, columns = NULL)
#  # Select rows and columns
#  tab_res <- tab_res[tab_res$Category %in% c("Means", "Variances"), c("Category", "lhs", "est", "se", "pval", "confint", "name")]
#  tab_res

## ----echo = FALSE, eval = run_everything--------------------------------------
#  res_final <- mx_switch_labels(res_step[[3]], param = "M[1,7]", decreasing = FALSE)
#  tab_res <- table_results(res_final, columns = NULL)
#  write.csv(tab_res, "lcga_tab_res.csv", row.names = FALSE)

## ----tabres, echo=FALSE, eval = eval_results----------------------------------
#  tab_res <- read.csv("lcga_tab_res.csv", stringsAsFactors = FALSE)
#  class(tab_res) <- c("tidy_results", "data.frame")
#  tab_res <- tab_res[tab_res$Category %in% c("Means", "Variances"), c("Category", "lhs", "est", "se", "pval", "confint", "name")]
#  knitr::kable(tab_res, digits = 2, caption = "Results from 3-class LCGA model", longtable = TRUE)

## ----params, echo = TRUE, eval = FALSE----------------------------------------
#  names(coef(res_final))

## ----echo = FALSE, eval = TRUE------------------------------------------------
c("mix3.weights[1,2]", "mix3.weights[1,3]", "vscl1", "vscl2", 
"vscl3", "vscl4", "vscl5", "vscl6", "class1.M[1,7]", "class1.M[1,8]", 
"class1.M[1,9]", "class2.M[1,7]", "class2.M[1,8]", "class2.M[1,9]", 
"class3.M[1,7]", "class3.M[1,8]", "class3.M[1,9]")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  wald_tests <- wald_test(res_final,
#                     "
#                     class1.M[1,7] = class2.M[1,7]&
#                     class1.M[1,7] = class3.M[1,7];
#                     class1.M[1,8] = class2.M[1,8]&
#                     class1.M[1,8] = class3.M[1,8];
#                     class1.M[1,9] = class2.M[1,9]&
#                     class1.M[1,9] = class3.M[1,9]")
#  # Rename the hypothesis
#  wald_tests$Hypothesis <- c("Mean i", "Mean step", "Mean slope")
#  knitr::kable(wald_tests, digits = 2, caption = "Wald tests")

## ----echo = FALSE, eval = run_everything--------------------------------------
#  wald_tests <- wald_test(res_final,
#                     "class1.M[1,7] = class2.M[1,7]&class1.M[1,7] = class3.M[1,7];class1.M[1,8] = class2.M[1,8]&class1.M[1,8] = class3.M[1,8];class1.M[1,9] = class2.M[1,9]&class1.M[1,9] = class3.M[1,9]")
#  # Rename the hypothesis
#  wald_tests$Hypothesis <- c("Mean i", "Mean step", "Mean slope")
#  write.csv(wald_tests, "lcga_wald_tests.csv", row.names = FALSE)

## ----waldtests, echo = FALSE, eval = eval_results-----------------------------
#  wald_tests <- read.csv("lcga_wald_tests.csv", stringsAsFactors = FALSE)
#  class(wald_tests) <- c("wald_test", "data.frame")
#  knitr::kable(wald_tests, digits = 2, caption = "Wald tests")

## ----makelcgaplot, echo = TRUE, eval = FALSE----------------------------------
#  p <- plot_growth(res_step[[3]], rawdata = TRUE, alpha_range = c(0, .05))
#  # Add Y-axis breaks in original scale
#  brks <- seq(0, 1, length.out = 5)
#  labs <- round(invbc(scales::rescale(brks, from = c(0, 1), to = rng_bc), lambda))
#  p <- p + scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = labs) + ylab("SCL (rescaled from Box-Cox)")
#  p

## ----echo = FALSE, eval = run_everything--------------------------------------
#  p <- plot_growth(res_step[[3]], rawdata = TRUE, alpha_range = c(0, .05))
#  # Add Y-axis breaks in original scale
#  brks <- seq(0, 1, length.out = 5)
#  labs <- round(invbc(scales::rescale(brks, from = c(0, 1), to = rng_bc), lambda))
#  p <- p + scale_y_continuous(breaks = seq(0, 1, length.out = 5), labels = labs) + ylab("SCL (rescaled from Box-Cox)")
#  ggsave("plot_traj.png", p, width = 150, height = 120, units = "mm")

## ----echo = FALSE, eval = eval_results, out.width="80%"-----------------------
#  knitr::include_graphics("plot_traj.png")

