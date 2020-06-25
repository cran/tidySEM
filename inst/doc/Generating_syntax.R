## ---- include = FALSE---------------------------------------------------------
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

## ---- eval = TRUE-------------------------------------------------------------
model %>%
  measurement() -> model
model

## -----------------------------------------------------------------------------
res <- lavaan(as_lavaan(model), data = df)
summary(res, estimates = FALSE)

## -----------------------------------------------------------------------------
model %>%
  estimate_lavaan()

## ----echo = FALSE, eval = run_mplus, message=FALSE, warning=FALSE-------------
#  library(MplusAutomation)
#  model %>%
#    estimate_mplus() -> res
#  #dput(capture.output(summary(res)))

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  library(MplusAutomation)
#  model %>%
#    estimate_mplus()

## ----eval = TRUE, echo = FALSE------------------------------------------------
# Display the results
cat(c("Estimated using ML ", "Number of obs: 301, number of (free) parameters: 30 ", 
"", "Model: Chi2(df = 24) = 85.306, p = 0 ", "Baseline model: Chi2(df = 36) = 918.852, p = 0 ", 
"", "Fit Indices: ", "", "CFI = 0.931, TLI = 0.896, SRMR = 0.06 ", 
"RMSEA = 0.092, 90% CI [0.071, 0.114], p < .05 = 0.001 ", "AIC = 7535.49, BIC = 7646.703 "
), sep = "\n")

## -----------------------------------------------------------------------------
dictionary(model)

## -----------------------------------------------------------------------------
syntax(model)

## -----------------------------------------------------------------------------
dictionary(model) %>%
  mutate(label = ifelse(label == "vis", "Visual", label))

## -----------------------------------------------------------------------------
syntax(model) %>%
  mutate(lhs = ifelse(lhs == "spe" & op == "=~", "tex", lhs)) %>%
  filter(!(lhs == "spe" | rhs == "spe")) -> syntax(model)

## -----------------------------------------------------------------------------
estimate_lavaan(model)

## -----------------------------------------------------------------------------
model %>%
  add_paths("vis ~ tex") %>%
  estimate_lavaan() %>%
  summary(estimates = TRUE)

## -----------------------------------------------------------------------------
model %>%
  add_paths("vis ~ tex", vis =~ spe_1) %>%
  estimate_lavaan()

