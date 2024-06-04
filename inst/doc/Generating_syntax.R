## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
run_mplus <- FALSE

## ----setup, echo = FALSE, message = FALSE-------------------------------------
library(tidySEM)
library(lavaan)
library(dplyr)

## -----------------------------------------------------------------------------
df <- HolzingerSwineford1939
names(df)

## -----------------------------------------------------------------------------
names(df)[grepl("^x", names(df))] <- c("vis_1", "vis_2", "vis_3", "tex_1", "tex_2", "tex_3", "spe_1", "spe_2", "spe_3")

## -----------------------------------------------------------------------------
model <- tidy_sem(df)
model

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  model |>
#    measurement() -> model
#  model

## ----echo = FALSE, eval = TRUE------------------------------------------------
model <- measurement(model)
model

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  model |>
#    estimate_lavaan()

## ----echo = FALSE, eval = TRUE------------------------------------------------
estimate_lavaan(model)

## ----include=FALSE------------------------------------------------------------
estimate_mx(model) -> res_mx

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  model |>
#    estimate_mx()

## ----echo = FALSE, eval = run_mplus, message=FALSE, warning=FALSE-------------
#  library(MplusAutomation)
#  model |>
#    estimate_mplus() -> res
#  #dput(capture.output(summary(res)))

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  library(MplusAutomation)
#  model |>
#    estimate_mplus()

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  # Display the results
#  cat(c("Estimated using ML ", "Number of obs: 301, number of (free) parameters: 30 ",
#  "", "Model: Chi2(df = 24) = 85.306, p = 0 ", "Baseline model: Chi2(df = 36) = 918.852, p = 0 ",
#  "", "Fit Indices: ", "", "CFI = 0.931, TLI = 0.896, SRMR = 0.06 ",
#  "RMSEA = 0.092, 90% CI [0.071, 0.114], p < .05 = 0.001 ", "AIC = 7535.49, BIC = 7646.703 "
#  ), sep = "\n")

## -----------------------------------------------------------------------------
dictionary(model)

## -----------------------------------------------------------------------------
syntax(model)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  dictionary(model) |>
#    mutate(label = ifelse(label == "vis", "Visual", label))

## ----echo = FALSE, eval = TRUE------------------------------------------------
tmp <- dictionary(model)
mutate(tmp, label = ifelse(label == "vis", "Visual", label))

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  syntax(model) |>
#    mutate(lhs = ifelse(lhs == "spe" & op == "=~", "tex", lhs)) |>
#    filter(!(lhs == "spe" | rhs == "spe")) -> syntax(model)

## ----eval = TRUE, echo = FALSE------------------------------------------------
tmp <- syntax(model)
tmp <- mutate(tmp, lhs = ifelse(lhs == "spe" & op == "=~", "tex", lhs))
tmp <- filter(tmp, !(lhs == "spe" | rhs == "spe"))
syntax(model) <- tmp

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  syntax(model) |>
#    mutate(free = ifelse(rhs == "spe_1", 1, free),
#    ustart = ifelse(rhs == "spe_1", NA, ustart)) -> syntax(model)

## ----eval = TRUE, echo = FALSE------------------------------------------------
syntax(model) |>
  mutate(free = ifelse(rhs == "spe_1", 1, free),
  ustart = ifelse(rhs == "spe_1", NA, ustart)) -> syntax(model)

## -----------------------------------------------------------------------------
estimate_lavaan(model)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  model |>
#    add_paths("vis ~ tex") |>
#    estimate_lavaan() |>
#    summary(estimates = TRUE)

## ----eval = TRUE, echo = FALSE------------------------------------------------
tmp <- add_paths(model, "vis ~ tex")
tmp <- estimate_lavaan(tmp)
summary(tmp, estimates = TRUE)

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  model |>
#    add_paths("vis ~ tex", vis =~ spe_1) |>
#    estimate_lavaan()

## ----eval = TRUE, echo = FALSE------------------------------------------------
tmp <- add_paths(model, "vis ~ tex", vis =~ spe_1)
estimate_lavaan(tmp)

