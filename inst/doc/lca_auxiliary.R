## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, comment = "#>",
  eval = requireNamespace("OpenMx", quietly = TRUE),
  out.width = "100%"
)
run_everything <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false"))))
show_results <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_SHOW_RESULTS", "false"))))
run_openmx <- requireNamespace("OpenMx", quietly = TRUE)
if (!run_openmx) run_everything <- show_results <- FALSE
options(scipen = 1, digits = 2)
if (show_results) knitr::opts_chunk$set(results = "markup")

## ----data---------------------------------------------------------------------
library(tidySEM)
df <- zegwaard_carecompass
indicators <- c("burdened", "trapped", "negaffect", "loneliness")
measurement_data <- df[, indicators]

## ----fit-model-code, eval=FALSE-----------------------------------------------
# set.seed(297)
# lpa_model <- mx_profiles(measurement_data, classes = 2)

## ----fit-model, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# set.seed(297)
# lpa_model <- tidySEM:::run_cached(mx_profiles(measurement_data, classes = 2))

## ----fit-results, eval=show_results && run_openmx-----------------------------
# table_fit(lpa_model)

## ----classification-code, eval=FALSE------------------------------------------
# cp <- class_prob(lpa_model)
# head(cp$individual)
# cp$mostlikely.class

## ----classification-results, echo=FALSE, eval=show_results && run_openmx------
# cp <- class_prob(lpa_model)
# format_numeric(head(cp$individual))
# format_numeric(cp$mostlikely.class)

## ----bch-mean-code, eval=FALSE------------------------------------------------
# bch_mean <- BCH(lpa_model, data = df$distance)
# table_results(bch_mean, columns = NULL)

## ----bch-mean-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# bch_mean <- tidySEM:::run_cached(BCH(lpa_model, data = df$distance))

## ----bch-mean-results, echo=FALSE, eval=show_results && run_openmx------------
# table_results(bch_mean, columns = NULL)

## ----bch-lr-code, eval=FALSE--------------------------------------------------
# lr_test(bch_mean, compare = "M")

## ----bch-lr-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# bch_lr <- tidySEM:::run_cached(lr_test(bch_mean, compare = "M"))

## ----bch-lr-results, echo=FALSE, eval=show_results && run_openmx--------------
# bch_lr

## ----bch-ordinal-code, eval=FALSE---------------------------------------------
# sex <- ordered(df$sexpatient, levels = c("male", "female"))
# bch_sex <- BCH(lpa_model, data = sex)
# table_prob(bch_sex)

## ----bch-ordinal-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# sex <- ordered(df$sexpatient, levels = c("male", "female"))
# bch_sex <- tidySEM:::run_cached(BCH(lpa_model, data = sex))

## ----bch-ordinal-results, echo=FALSE, eval=show_results && run_openmx---------
# table_prob(bch_sex)

## ----bch-regression-code, eval=FALSE------------------------------------------
# df_aux <- df[, c("freqvisit", "distance"), drop = FALSE]
# df_aux$freqvisit <- as.numeric(df_aux$freqvisit)
# bch_reg <- BCH(lpa_model,
#                model = "freqvisit ~ distance",
#                data = df_aux)
# table_results(bch_reg, columns = c("label", "est_sig", "se", "pval", "confint", "group", "matrix"))
# lr_test(bch_reg, compare = "A")

## ----bch-regression-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# df_aux <- df[, c("freqvisit", "distance"), drop = FALSE]
# df_aux$freqvisit <- as.numeric(df_aux$freqvisit)
# bch_reg <- tidySEM:::run_cached(BCH(lpa_model, model = "freqvisit ~ distance", data = df_aux))

## ----bch-regression-lr-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# bch_reg_lr <- tidySEM:::run_cached(lr_test(bch_reg, compare = "A"))

## ----bch-regression-results, echo=FALSE, eval=show_results && run_openmx------
# table_results(bch_reg, columns = NULL)
# bch_reg_lr

## ----two-step-code, eval=FALSE------------------------------------------------
# two_mean <- two_step(lpa_model, data = df$distance)
# table_results(two_mean, columns = NULL)

## ----two-step-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# two_mean <- tidySEM:::run_cached(two_step(lpa_model, data = df$distance))

## ----two-step-results, echo=FALSE, eval=show_results && run_openmx------------
# table_results(two_mean, columns = NULL)

## ----two-step-se-code, eval=FALSE---------------------------------------------
# attr(two_mean, "se_robust")

## ----two-step-se-results, echo=FALSE, eval=show_results && run_openmx---------
# attr(two_mean, "se_robust")

## ----two-reg-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# two_reg <- tidySEM:::run_cached(two_step(
#   lpa_model,
#   model = "freqvisit ~ distance",
#   data = df_aux
# ))

## ----echo = TRUE, eval = FALSE------------------------------------------------
# two_reg <- two_step(
#   lpa_model,
#   model = "freqvisit ~ distance",
#   data = df_aux
# )
# table_results(two_reg)

## ----pseudo-lm-code, eval=FALSE-----------------------------------------------
# set.seed(2027)
# pseudo_mean <- pseudo_class(
#   x = lpa_model,
#   model = lm(distance ~ I(factor(class)) - 1, data = data),
#   data = df[, "distance", drop = FALSE],
#   m = 20
# )
# pseudo_mean

## ----pseudo-lm-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# set.seed(2027)
# pseudo_mean <- tidySEM:::run_cached(pseudo_class(
#   x = lpa_model,
#   model = lm(distance ~ I(factor(class)) - 1, data = data),
#   data = df[, "distance", drop = FALSE], m = 20
# ))

## ----pseudo-lm-results, echo=FALSE, eval=show_results && run_openmx-----------
# pseudo_mean

## ----pseudo-ram-code, eval=FALSE----------------------------------------------
# set.seed(2028)
# pseudo_ram <- pseudo_class(
#   x = lpa_model,
#   model = "distance ~ 1",
#   data = df[, "distance", drop = FALSE],
#   m = 20
# )
# pseudo_ram

## ----pseudo-ram-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# set.seed(2028)
# pseudo_ram <- tidySEM:::run_cached(pseudo_class(
#   x = lpa_model, model = "distance ~ 1",
#   data = df[, "distance", drop = FALSE], m = 20
# ))

## ----pseudo-ram-results, echo=FALSE, eval=show_results && run_openmx----------
# pseudo_ram

## ----pseudo-predict-code, eval=FALSE------------------------------------------
# set.seed(2030)
# pseudo_pred <- pseudo_class(
#   x = lpa_model,
#   model = function(data) {
#     nnet::multinom(class ~ distance + freqvisit,
#                    data = data, trace = FALSE)
#   },
#   data = df_aux,
#   m = 20
# )
# pseudo_pred

## ----pseudo-predict-fit, include=FALSE, eval=run_openmx && (run_everything || show_results)----
# if (requireNamespace("nnet", quietly = TRUE)) {
#   set.seed(2030)
#   pseudo_pred <- tidySEM:::run_cached(pseudo_class(
#     x = lpa_model,
#     model = function(data) nnet::multinom(
#       class ~ distance + freqvisit, data = data, trace = FALSE),
#     data = df_aux, m = 20
#   ))
# }

## ----pseudo-predict-results, echo=FALSE, eval=show_results && run_openmx------
# if (exists("pseudo_pred", inherits = FALSE)) pseudo_pred

