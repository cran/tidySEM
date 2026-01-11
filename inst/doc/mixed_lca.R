## ----include = FALSE----------------------------------------------------------
run_chunks <- requireNamespace("OpenMx", quietly = TRUE)
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  eval = run_chunks,
  comment = "#>"
)

show_plots <- FALSE
if(suppressWarnings(tryCatch({isTRUE(as.logical(readLines("pkgdown.txt")))}, error = function(e){FALSE}))){
  show_plots <- TRUE
}
run_everything = suppressWarnings(tryCatch({isTRUE(as.logical(readLines("run_everything.txt")))}, error = function(e){FALSE}))

## ----eval = run_chunks--------------------------------------------------------
library(tidySEM)
library(OpenMx)

## ----eval = run_chunks--------------------------------------------------------
set.seed(10)
n <- 200

# Set class-specific means
class_means <- c(rep(0, floor(0.3 * n)),
rep(2, ceiling(0.7 * n)))

# Simulate continuous indicators
df <- rnorm(4 * n, mean = rep(class_means, 4))
df <- matrix(df, nrow = n)
df <- t(t(df) * c(1, 2, 0.5, 1))
df <- data.frame(df)
names(df) <- paste0("X", 1:4)

# Convert one indicator to ordinal
df$X4 <- cut(df$X4, breaks = 3, labels = FALSE)
df$X4 <- mxFactor(df$X4, levels = 1:3)

## ----eval = run_chunks--------------------------------------------------------
res_2 <- mx_mixed_lca(
data = df,
classes = 2
)

## ----eval = run_chunks--------------------------------------------------------
class(res_2)

## ----eval = FALSE-------------------------------------------------------------
# res_1_3 <- mx_mixed_lca(
#   data = df,
#   classes = 1:3
# )

## ----eval = run_everything, echo = FALSE--------------------------------------
# res_1_3 <- mx_mixed_lca(
#   data = df,
#   classes = 1:3
# )

## ----eval = FALSE-------------------------------------------------------------
# table_fit(res_1_3)

## ----eval = run_everything, echo = FALSE--------------------------------------
# tab <- table_fit(res_1_3)
# write.csv(tab, "mixed_lca_table_fit.csv", row.names = FALSE)

## ----eval = TRUE, echo=FALSE--------------------------------------------------
tab <- read.csv("mixed_lca_table_fit.csv", stringsAsFactors = FALSE)
class(tab) <- c("tidy_fit"   , "data.frame")
tab

## ----eval = FALSE, echo = TRUE------------------------------------------------
# library(future)
# library(progressr)
# plan(multisession) # Parallel processing for Windows
# handlers("progress") # Progress bar
# set.seed(1)
# res_blrt <- BLRT(res_1_3, replications = 100)
# res_blrt

## ----eval = run_everything, echo = FALSE, warning=FALSE-----------------------
# library(future)
# library(progressr)
# plan(multisession) # Parallel processing for Windows
# handlers("progress") # Progress bar
# set.seed(1)
# res_blrt <- BLRT(res_1_3, replications = 100)
# write.csv(res_blrt, "mixed_lca_res_blrt.csv", row.names = FALSE)

## ----eval = TRUE, echo = FALSE------------------------------------------------
res_blrt <- read.csv("mixed_lca_res_blrt.csv", stringsAsFactors = FALSE)
class(res_blrt) <- c("LRT", "data.frame")
attr(res_blrt, "type") <- "Bootstrapped"
res_blrt

## ----eval = run_everything, echo = FALSE--------------------------------------
# set.seed(1)
# res_pmc_srmr <- pmc_srmr(res_1_3)
# write.csv(res_pmc_srmr, "mixed_lca_res_pmc_srmr.csv", row.names = F)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# set.seed(1)
# res_pmc_srmr <- pmc_srmr(res_1_3)
# res_pmc_srmr

## ----echo = FALSE-------------------------------------------------------------
res_pmc_srmr <- read.csv("mixed_lca_res_pmc_srmr.csv", stringsAsFactors = F)
res_pmc_srmr

## ----eval = run_everything, echo = FALSE--------------------------------------
# tmp <- class_prob(res_1_3[[2]], c("sum.posterior", "sum.mostlikely"))
# write.csv(tmp$sum.posterior, "mixed_lca_posterior.csv", row.names = F)
# write.csv(tmp$sum.mostlikely, "mixed_lca_mostlikely.csv", row.names = F)

## ----eval = FALSE, echo = TRUE------------------------------------------------
# class_prob(res_1_3[[2]], c("sum.posterior", "sum.mostlikely"))

## ----eval = TRUE, echo = FALSE------------------------------------------------
list(sum.posterior = read.csv("mixed_lca_posterior.csv", stringsAsFactors = F), sum.mostlikely = read.csv("mixed_lca_mostlikely.csv", stringsAsFactors = F))

## ----eval = FALSE, echo = TRUE------------------------------------------------
# table_fit(res_1_3[[2]])

## ----eval = run_everything, echo = FALSE--------------------------------------
# tmp <- table_fit(res_1_3[[2]])
# write.csv(tmp, "mixed_lca_fit2.csv", row.names = F)

## ----eval = TRUE, echo = FALSE------------------------------------------------
tmp <- read.csv("mixed_lca_fit2.csv", stringsAsFactors = F)
class(tmp) <- c("tidy_fit", "data.frame")
tmp

## ----eval = FALSE, echo = TRUE------------------------------------------------
# table_results(res_1_3[[2]])

## ----eval = run_everything, echo = FALSE--------------------------------------
# tmp <- table_results(res_1_3[[2]])
# write.csv(tmp, "mixed_lca_res2.csv", row.names = F)

## ----eval = TRUE, echo = FALSE------------------------------------------------
tmp <- read.csv("mixed_lca_res2.csv", stringsAsFactors = F)
class(tmp) <- c("tidy_results","data.frame"  )
tmp

## ----eval = FALSE, echo = TRUE------------------------------------------------
# table_prob(res_1_3[[2]])

## ----eval = run_everything, echo = FALSE--------------------------------------
# tmp <- table_prob(res_1_3[[2]])
# write.csv(tmp, "mixed_lca_prob2.csv", row.names = F)

## ----eval = TRUE, echo = FALSE------------------------------------------------
read.csv("mixed_lca_prob2.csv", stringsAsFactors = F)

## ----eval = run_everything, echo = TRUE---------------------------------------
# res_2_free <- mx_mixed_lca(
#   data = df,
#   classes = 2,
#   variances = "varying"
# )

## ----echo = TRUE, eval = FALSE------------------------------------------------
# compare <- list(
#   fixed_covs = res_1_3[[2]],
#   free_covs = res_2_free)
# table_fit(compare)

## ----eval = run_everything, echo = FALSE--------------------------------------
# compare <- list(
#   fixed_covs = res_1_3[[2]],
#   free_covs = res_2_free)
# fit_compare <- table_fit(compare)
# write.csv(fit_compare, "mixed_lca_compare.csv", row.names = FALSE)

## ----tabfitcomp, echo = FALSE, eval = TRUE------------------------------------
fit_compare <- read.csv("mixed_lca_compare.csv", stringsAsFactors = FALSE)
class(fit_compare) <- c("tidy_fit", "data.frame")
fit_compare

## ----eval = run_everything, echo = FALSE--------------------------------------
# p <- plot_profiles(res_1_3[[2]], variables = c("X1", "X2", "X3"))
# ggplot2::ggsave("mixed_lca_profiles.png", p, device = "png", width = 4, height = 3, dpi = 150)

## ----eval = FALSE, echo = TRUE------------------------------------------------
# plot_profiles(res_1_3[[2]], variables = c("X1", "X2", "X3"))

## ----eval = TRUE, echo = FALSE------------------------------------------------
knitr::include_graphics("mixed_lca_profiles.png")

## ----eval = run_everything, echo = FALSE--------------------------------------
# p <- plot_bivariate(res_1_3[[2]], variables = c("X1", "X2", "X3"))
# ggplot2::ggsave("mixed_lca_bivariate.png", p, device = "png", width = 4, height = 4, dpi = 150)

## ----eval = FALSE, echo = TRUE------------------------------------------------
# plot_bivariate(res_1_3[[2]], variables = c("X1", "X2", "X3"))

## ----eval = TRUE, echo = FALSE------------------------------------------------
knitr::include_graphics("mixed_lca_bivariate.png")

## ----eval = run_everything, echo = FALSE--------------------------------------
# p <- plot_prob(res_1_3[[2]])
# ggplot2::ggsave("mixed_lca_prob.png", p, device = "png", width = 3, height = 3, dpi = 150)

## ----eval = FALSE, echo = TRUE------------------------------------------------
# plot_prob(res_1_3[[2]])

## ----eval = TRUE, echo = FALSE------------------------------------------------
knitr::include_graphics("mixed_lca_prob.png")

