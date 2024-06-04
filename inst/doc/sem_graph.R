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
generate_pngs <- TRUE

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  graph_sem(layout = matrix("x")) + coord_fixed()

## ----echo = FALSE, out.width='20%'--------------------------------------------
p <- graph_sem(layout = matrix("x")) + coord_fixed()
if(generate_pngs) ggsave("var_obs.png", p, device = "png", width= 1, height = .7)
knitr::include_graphics("var_obs.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  graph_sem(nodes = data.frame(name = "x", shape = "oval"), layout = matrix("x"), fix_coord = TRUE)

## ----echo = FALSE, out.width='20%'--------------------------------------------
p <- graph_sem(nodes = data.frame(name = "x", shape = "oval"), layout = matrix("x"), fix_coord = TRUE)
if(generate_pngs) ggsave("var_lat.png", p, device = "png", width= 1, height = .7)
knitr::include_graphics("var_lat.png")

## ----echo = FALSE, out.width='20%'--------------------------------------------
p <- ggplot(data.frame(x = 1, y = 1, xend = 4, yend = 1), aes(x=x, y=y,xend=xend, yend=yend)) +geom_segment(arrow = arrow(type = "closed"))+theme_void()
if(generate_pngs) ggsave("arrow.png", p, device = "png", width= 3, height = .7)
knitr::include_graphics("arrow.png")

## ----echo = FALSE, out.width='20%'--------------------------------------------
p <- ggplot(data.frame(x = 1, y = 1, xend = 2, yend = 1), aes(x=x, y=y,xend=xend, yend=yend)) + geom_curve(linetype = 2)+theme_void()+scale_y_continuous(limits = c(0,2))
if(generate_pngs) ggsave("curve.png", p, device = "png", width= 4, height = 2)
knitr::include_graphics("curve.png")

## ----echo = FALSE, out.width='10%'--------------------------------------------
if(generate_pngs){
  p <- tidySEM:::.plot_variances(p = ggplot(NULL), df = data.frame(edge_xmin = 1, edge_ymin = 1,
                                                   edge_xmax = 1, edge_ymax = 1,
                                                   connect_from = "top",
                                                   connect_to = "top", label = "e",
                                                   arrow = "both"), diameter = 1,
                          text_size = 8)+theme_void()
  ggsave("error.png", p, device = "png", width= 2, height = 2)
} 
knitr::include_graphics("error.png")

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  df <- iris[, 1:2]
#  names(df) <- c("x", "y")
#  sem("y ~ x", df) %>%
#    graph_sem(spacing_x = 2.5, fix_coord = TRUE)

## ----echo = FALSE, out.width='50%'--------------------------------------------
df <- iris[, 1:2]
names(df) <- c("x", "y")
sem("y ~ x", df) %>%
  graph_sem(spacing_x = 2.5, fix_coord = TRUE) -> p
if(generate_pngs) ggsave("mod_reg.png", p, device = "png", width= 4, height = 1)
knitr::include_graphics("mod_reg.png")

## ----echo = TRUE, eval=FALSE--------------------------------------------------
#  df <- iris[ , c(1,3:4)]
#  names(df) <- paste0("y_", 1:3)
#  
#  tidy_sem(df) %>%
#    measurement() %>%
#    estimate_lavaan() %>%
#    graph_sem()

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='50%'----------
df <- iris[ , c(1,3:4)]
names(df) <- paste0("y_", 1:3)

tidy_sem(df) %>%
  measurement() %>%
  estimate_lavaan() %>%
  graph_sem(fix_coord = TRUE) -> p
if(generate_pngs) ggsave("mod_meas.png", p, device = "png", width= 4, height = 4)
knitr::include_graphics("mod_meas.png")

## ----eval= FALSE--------------------------------------------------------------
#  df <- iris[ , 1:4]
#  names(df) <- c("y_1", "x", "y_2", "y_3")
#  
#  tidy_sem(df) %>%
#    measurement() %>%
#    add_paths(y ~ x, x ~~ x, x ~1) %>%
#    estimate_lavaan() %>%
#    graph_sem()

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='70%'----------
df <- iris[ , 1:4]
names(df) <- c("y_1", "x", "y_2", "y_3")
set.seed(58)
tidy_sem(df) %>%
  measurement() %>%
  add_paths(y ~ x, x ~~ x, x ~1) %>%
  estimate_lavaan() %>%
  graph_sem(fix_coord = TRUE) -> p
if(generate_pngs) ggsave("mod_sem1.png", p, device = "png", width= 6, height = 6)
knitr::include_graphics("mod_sem1.png")

## ----eval= FALSE--------------------------------------------------------------
#  tidy_sem(df) %>%
#    measurement() %>%
#    add_paths(y ~ x, x ~~ x) %>%
#    estimate_lavaan() %>%
#    graph_sem(layout =
#                get_layout("",     "x",    "",
#                           "",     "y",    "",
#                           "y_1", "y_2", "y_3", rows = 3))

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='70%'----------
tidy_sem(df) %>%
  measurement() %>%
  add_paths(y ~ x, x ~~ x) %>%
  estimate_lavaan() %>%
  graph_sem(layout =
              get_layout("",     "x",    "",
                         "",     "y",    "",
                         "y_1", "y_2", "y_3", rows = 3), fix_coord = TRUE) -> p
if(generate_pngs) ggsave("mod_sem2.png", p, device = "png", width= 6, height = 6)
knitr::include_graphics("mod_sem2.png")

