## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(tidySEM)
library(lavaan)
library(dplyr)
library(ggplot2)
include_code <- function(expr){
  txt <- deparse(substitute(expr))
  eval.parent(expr)
  cat(c("``` r", txt, "```"), sep = "\n")
}
run_everything <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_RUN_EVERYTHING", "false"))))
show_results <- isTRUE(as.logical(toupper(Sys.getenv("TIDYSEM_SHOW_RESULTS", "false"))))

## ----echo = TRUE, eval = FALSE------------------------------------------------
# graph_sem(layout = matrix("x")) + coord_fixed()

## ----intro-graph-01, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# p <- graph_sem(layout = matrix("x")) + coord_fixed()
# p
# }, type = "png", width = 1, height = .7, dpi = 150)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# graph_sem(nodes = data.frame(name = "x", shape = "oval"), layout = matrix("x"), fix_coord = TRUE)

## ----intro-graph-02, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# p <- graph_sem(nodes = data.frame(name = "x", shape = "oval"), layout = matrix("x"), fix_coord = TRUE)
# p
# }, type = "png", width = 1, height = .7, dpi = 150)

## ----intro-graph-03, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# p <- ggplot(data.frame(x = 1, y = 1, xend = 4, yend = 1), aes(x=x, y=y,xend=xend, yend=yend)) +geom_segment(arrow = arrow(type = "closed"))+theme_void()
# p
# }, type = "png", width = 3, height = .7, dpi = 150)

## ----intro-graph-04, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# p <- ggplot(data.frame(x = 1, y = 1, xend = 2, yend = 1), aes(x=x, y=y,xend=xend, yend=yend)) + geom_curve(linetype = 2)+theme_void()+scale_y_continuous(limits = c(0,2))
# p
# }, type = "png", width = 4, height = 2, dpi = 150)

## ----intro-graph-05, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
#   p <- tidySEM:::.plot_variances(p = ggplot(NULL), df = data.frame(edge_xmin = 1, edge_ymin = 1,
#                                                    edge_xmax = 1, edge_ymax = 1,
#                                                    connect_from = "top",
#                                                    connect_to = "top", label = "e",
#                                                    arrow = "both"), diameter = 1,
#                           text_size = 8)+theme_void()
# p
# }, type = "png", width = 2, height = 2, dpi = 150)

## ----echo = TRUE, eval = FALSE------------------------------------------------
# df <- iris[, 1:2]
# names(df) <- c("x", "y")
# sem("y ~ x", df) %>%
#   graph_sem(spacing_x = 2.5, fix_coord = TRUE)

## ----intro-graph-06, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# df <- iris[, 1:2]
# names(df) <- c("x", "y")
# sem("y ~ x", df) %>%
#   graph_sem(spacing_x = 2.5, fix_coord = TRUE) -> p
# p
# }, type = "png", width = 4, height = 1, dpi = 150)

## ----echo = TRUE, eval=FALSE--------------------------------------------------
# df <- iris[ , c(1,3:4)]
# names(df) <- paste0("y_", 1:3)
# 
# tidy_sem(df) %>%
#   measurement() %>%
#   estimate_lavaan() %>%
#   graph_sem()

## ----intro-graph-07, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# df <- iris[ , c(1,3:4)]
# names(df) <- paste0("y_", 1:3)
# 
# tidy_sem(df) %>%
#   measurement() %>%
#   estimate_lavaan() %>%
#   graph_sem(fix_coord = TRUE) -> p
# p
# }, type = "png", width = 4, height = 4, dpi = 150)

## ----eval= FALSE--------------------------------------------------------------
# df <- iris[ , 1:4]
# names(df) <- c("y_1", "x", "y_2", "y_3")
# 
# tidy_sem(df) %>%
#   measurement() %>%
#   add_paths(y ~ x, x ~~ x, x ~1) %>%
#   estimate_lavaan() %>%
#   graph_sem()

## ----intro-graph-08, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# df <- iris[ , 1:4]
# names(df) <- c("y_1", "x", "y_2", "y_3")
# set.seed(58)
# tidy_sem(df) %>%
#   measurement() %>%
#   add_paths(y ~ x, x ~~ x, x ~1) %>%
#   estimate_lavaan() %>%
#   graph_sem(fix_coord = TRUE) -> p
# p
# }, type = "png", width = 6, height = 6, dpi = 150)

## ----eval= FALSE--------------------------------------------------------------
# tidy_sem(df) %>%
#   measurement() %>%
#   add_paths(y ~ x, x ~~ x) %>%
#   estimate_lavaan() %>%
#   graph_sem(layout =
#               get_layout("",     "x",    "",
#                          "",     "y",    "",
#                          "y_1", "y_2", "y_3", rows = 3))

## ----intro-graph-09, echo=FALSE, eval=run_everything || show_results, out.width="100%"----
# tidySEM:::run_cached({
# tidy_sem(df) %>%
#   measurement() %>%
#   add_paths(y ~ x, x ~~ x) %>%
#   estimate_lavaan() %>%
#   graph_sem(layout =
#               get_layout("",     "x",    "",
#                          "",     "y",    "",
#                          "y_1", "y_2", "y_3", rows = 3), fix_coord = TRUE) -> p
# p
# }, type = "png", width = 6, height = 6, dpi = 150)

