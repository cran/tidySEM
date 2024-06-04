## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
run_mplus <- FALSE

## ----setup, message=FALSE-----------------------------------------------------
library(tidySEM)
library(lavaan)
library(MplusAutomation)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df <- HolzingerSwineford1939
#  names(df)[7:15] <- paste0(rep(c("vis", "tex", "spe"), each = 3), "_", rep(1:3, 3))
#  df |>
#    subset(select = c("school", "vis_1", "vis_2", "vis_3", "tex_1", "tex_2", "tex_3", "spe_1",
#  "spe_2", "spe_3")) -> df

## ----echo = FALSE, eval = TRUE------------------------------------------------
df <- HolzingerSwineford1939
names(df)[7:15] <- paste0(rep(c("vis", "tex", "spe"), each = 3), "_", rep(1:3, 3))
subset(df, select = c("school", "vis_1", "vis_2", "vis_3", "tex_1", "tex_2", "tex_3", "spe_1", 
"spe_2", "spe_3")) -> df

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df |>
#    tidy_sem() |>
#    measurement() -> model

## ----eval = TRUE, echo = FALSE------------------------------------------------
model <- measurement(tidy_sem(df))

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  model |>
#    estimate_lavaan() -> fit_lav

## ----echo = FALSE, eval = TRUE------------------------------------------------
estimate_lavaan(model) -> fit_lav

## -----------------------------------------------------------------------------
table_results(fit_lav)

## -----------------------------------------------------------------------------
table_fit(fit_lav)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  model |>
#    estimate_mx() -> fit_mx
#  table_results(fit_mx)
#  table_fit(fit_mx)

## ----eval = TRUE, echo = FALSE------------------------------------------------
estimate_mx(model) -> fit_mx
table_results(fit_mx)
table_fit(fit_mx)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  fit_mplus <- mplusModeler(mplusObject(VARIABLE = "grouping IS school (1 = GW 2 = Pas);",
#                                  MODEL = c("visual BY vis_1 vis_2 vis_3;",
#                                            "textual BY tex_1 tex_2 tex_3;",
#                                            "speed BY spe_1 spe_2 spe_3;"),
#                                  usevariables = names(df),
#                                  rdata = df),
#                      modelout = "example.inp",
#                      run = 1L)
#  table_results(fit_mplus)

## ----eval = run_mplus, echo = FALSE-------------------------------------------
#  # fit <- mplusModeler(mplusObject(VARIABLE = "grouping IS school (1 = GW 2 = Pas);",
#  #                                 MODEL = c("visual BY x1 x2 x3;",
#  #                                           "textual BY x4 x5 x6;",
#  #                                           "speed BY x7 x8 x9;"),
#  #                                 usevariables = c(paste0("x", 1:9), "school"),
#  #                                 rdata = HolzingerSwineford1939),
#  #                     modelout = "example.inp",
#  #                     run = 1L)
#  # file.remove(list.files(pattern = "^example.+(inp|out|dat)$"))
#  #dput(fit$results$parameters)
#  #dput(fit, file = "mplusfit.R")

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Read the results
fit_mplus <- source("mplusfit.R")
fit_mplus <- fit_mplus$value

## -----------------------------------------------------------------------------
table_results(fit_mplus)
table_fit(fit_mplus)

