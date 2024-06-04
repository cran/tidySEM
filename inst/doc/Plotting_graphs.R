## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
generate_pngs <- TRUE

## ----setup, message=FALSE, warning = FALSE------------------------------------
library(tidySEM)
library(lavaan)
library(ggplot2)
library(dplyr)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  library(lavaan)
#  HS.model <- ' visual  =~ x1 + x2 + x3
#                textual =~ x4 + x5 + x6
#                speed   =~ x7 + x8 + x9 '
#  fit <- cfa(HS.model, data=HolzingerSwineford1939)

## ----eval = TRUE, echo = FALSE, message=FALSE---------------------------------
library(lavaan)
suppressWarnings({
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
fit <- cfa(HS.model, data=HolzingerSwineford1939)
})

## ----echo = TRUE, eval = FALSE------------------------------------------------
#  graph_sem(model = fit)

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='100%'---------
p <- graph_sem(model = fit, text_size = 2, fix_coord = TRUE)
if(generate_pngs) ggsave("pgfig1.png", p, device = "png", width= 9.5, height = 3)
knitr::include_graphics("pgfig1.png")

## ----eval = TRUE, echo = TRUE-------------------------------------------------
get_layout(fit)

## ----eval = FALSE, echo = FALSE, message=FALSE, warning = FALSE---------------
#  library(lavaan)
#  suppressWarnings({
#  fit <- cfa(' visual  =~ x1 + x2 + x3 ',
#             data = HolzingerSwineford1939[1:50, ])
#  get_layout(fit)
#  })

## ----message=FALSE, warning = FALSE-------------------------------------------
get_layout(fit, layout_algorithm = "layout_in_circle")
get_layout(fit, layout_algorithm = "layout_on_grid")

## -----------------------------------------------------------------------------
get_layout("c", NA,  "d",
           NA,  "e", NA, rows = 2)

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  read.csv("example.csv")

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  write.csv(matrix(c("x1", "x2",  "x3", "",  "visual", ""), nrow = 2, byrow = TRUE), file = "example.csv", row.names = FALSE)
#  read.csv("example.csv")
#  file.remove("example.csv")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  read.table("clipboard", sep = "\t")

## ----eval = FALSE, echo = TRUE------------------------------------------------
#  read.table(pipe("pbpaste"), sep="\t")

## ----echo = FALSE, eval = TRUE------------------------------------------------
structure(list(V1 = structure(2:1, .Label = c("", "x1"), class = "factor"), 
    V2 = structure(2:1, .Label = c("visual", "x2"), class = "factor"), 
    V3 = structure(2:1, .Label = c("", "x3"), class = "factor")), class = "data.frame", row.names = c(NA, 
-2L))
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
fit <- cfa(HS.model, data=HolzingerSwineford1939)

## -----------------------------------------------------------------------------
get_layout("x", "y", rows = 1)

## -----------------------------------------------------------------------------
get_layout("", "m", "",
           "x", "", "y", rows = 2)

## -----------------------------------------------------------------------------
get_layout("", "F", "",
           "y1", "y2", "y3", rows = 2)

## -----------------------------------------------------------------------------
lay <- get_layout("", "", "visual","","textual","","speed","", "",
                  "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", rows = 2)

## ----eval = FALSE-------------------------------------------------------------
#  graph_sem(fit, layout = lay)

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='100%'---------
p <- graph_sem(fit, layout = lay) + coord_fixed()
if(generate_pngs) ggsave("pgfig2_1.png", p, device = "png", width= 9.5, height = 3)
knitr::include_graphics("pgfig2_1.png")

## ----eval = TRUE--------------------------------------------------------------
get_nodes(fit)

## ----eval = FALSE, echo = FALSE, results= "asis"------------------------------
#  knitr::kable(get_nodes(fit))

## -----------------------------------------------------------------------------
get_edges(fit)

## -----------------------------------------------------------------------------
get_edges(fit, label = paste(est, confint))

## -----------------------------------------------------------------------------
fit <- cfa(HS.model, data=HolzingerSwineford1939, meanstructure = TRUE)

## -----------------------------------------------------------------------------
get_nodes(fit)

## -----------------------------------------------------------------------------
get_nodes(fit, label = paste0(name, "\n", est, " ", confint))

## -----------------------------------------------------------------------------
p <- prepare_graph(fit)
edges(p)

## -----------------------------------------------------------------------------
prepare_graph(fit) %>%
  edit_graph({ label = paste(est_sig_std, "\n", confint_std) }) %>%
  plot()

## -----------------------------------------------------------------------------
prepare_graph(fit) %>%
  edit_graph({ label = paste(est_sig_std, "\n", confint_std) }) %>%
  edit_graph({ label = paste(name, "\n", est_sig_std, "\n", confint_std) }, element = "nodes") %>%
  plot()

## -----------------------------------------------------------------------------
prepare_graph(fit) %>%
  edit_graph({ label_color = "blue" }) %>%
  plot()

## -----------------------------------------------------------------------------
graph_data <- prepare_graph(model = fit, layout = lay)

## -----------------------------------------------------------------------------
nodes(graph_data)
edges(graph_data)

## ----message=FALSE------------------------------------------------------------
library(dplyr)
library(stringr)
nodes(graph_data) <- nodes(graph_data) %>%
  mutate(label = str_to_title(label))

## ----echo = FALSE-------------------------------------------------------------
# $label[1:3] <- str_to_title(nodes(graph_data)$label[1:3])
graph_data <- prepare_graph(model = fit, layout = lay)

## ----eval = FALSE-------------------------------------------------------------
#  edges(graph_data) %>%
#    mutate(connect_from = replace(connect_from, is.na(curvature), "bottom")) %>%
#    mutate(connect_to = replace(connect_to, is.na(curvature), "top")) -> edges(graph_data)

## ----eval = FALSE-------------------------------------------------------------
#  plot(graph_data)

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='100%'---------
p <- plot(graph_data) + coord_fixed()
if(generate_pngs) ggsave("pgfig2.png", p, device = "png", width= 9.5, height = 3)
knitr::include_graphics("pgfig2.png")

## ----eval = FALSE-------------------------------------------------------------
#  graph_sem(model = fit, layout = lay, angle = 170)

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='100%'---------
p <- graph_sem(model = fit, layout = lay, angle = 170)
if(generate_pngs) ggsave("pgfig3.png", p, device = "png", width= 9.5, height = 3)
knitr::include_graphics("pgfig3.png")

## ----eval = FALSE-------------------------------------------------------------
#  edg <- data.frame(from = "x",
#                    to = "y",
#                    linetype = 2,
#                    colour = "red",
#                    size = 2,
#                    alpha = .5)
#  
#  graph_sem(edges = edg, layout = get_layout("x", "y", rows = 1))

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='30%'----------
edg <- data.frame(from = "x",
                  to = "y",
                  linetype = 2,
                  colour = "red",
                  size = 2,
                  alpha = .5)

p <- graph_sem(edges = edg, layout = get_layout("x", "y", rows = 1))

if(generate_pngs) ggsave("pgfig4.png", p, device = "png", width= 4, height = 1)
knitr::include_graphics("pgfig4.png")

## ----eval = FALSE-------------------------------------------------------------
#  edg <- data.frame(from = "x",
#                    to = "y")
#  nod <- data.frame(name = c("x", "y"),
#                      shape = c("rect", "oval"),
#                      linetype = c(2, 2),
#                      colour = c("blue", "blue"),
#                      fill = c("blue", "blue"),
#                      size = c(2, 2),
#                      alpha = .5)
#  graph_sem(edges = edg, nodes = nod, layout = get_layout("x", "y", rows = 1))

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='30%'----------
edg <- data.frame(from = "x",
                  to = "y")
nod <- data.frame(name = c("x", "y"),
                    shape = c("rect", "oval"),
                    linetype = c(2, 2),
                    colour = c("blue", "blue"),
                    fill = c("blue", "blue"),
                    size = c(2, 2),
                    alpha = .5)
p <- graph_sem(edges = edg, nodes = nod, layout = get_layout("x", "y", rows = 1))

if(generate_pngs) ggsave("pgfig5.png", p, device = "png", width= 4, height = 1)
knitr::include_graphics("pgfig5.png")

## ----eval = FALSE-------------------------------------------------------------
#  edges(graph_data) %>%
#    mutate(colour = "black") %>%
#    mutate(colour = replace(colour, from == "visual" & to == "x2", "red")) %>%
#    mutate(linetype = 1) %>%
#    mutate(linetype = replace(linetype, from == "visual" & to == "x2", 2)) %>%
#    mutate(alpha = 1) %>%
#    mutate(alpha = replace(alpha, from == "visual" & to == "x2", .5)) -> edges(graph_data)
#  plot(graph_data)

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='100%'---------
edges(graph_data) %>%
  mutate(colour = "black") %>%
  mutate(colour = replace(colour, from == "visual" & to == "x2", "red")) %>%
  mutate(linetype = 1) %>%
  mutate(linetype = replace(linetype, from == "visual" & to == "x2", 2)) %>%
  mutate(alpha = 1) %>%
  mutate(alpha = replace(alpha, from == "visual" & to == "x2", .5)) -> edges(graph_data)
p <- plot(graph_data)

if(generate_pngs) ggsave("pgfig6.png", p, device = "png", width= 9.5, height = 3)
knitr::include_graphics("pgfig6.png")

## ----eval = FALSE-------------------------------------------------------------
#  edg <- data.frame(from = "x",
#                    to = "y",
#                    label = "text",
#                    label_colour = "blue",
#                    label_fill = "red",
#                    label_size = 6,
#                    label_alpha = .5,
#                    label_family = "mono",
#                    label_fontface = "bold",
#                    label_hjust = "left",
#                    label_vjust = "top",
#                    label_lineheight = 1.5,
#                    label_location = .2
#                    )
#  
#  graph_sem(edges = edg, layout = get_layout("x", "y", rows = 1))

## ----echo = FALSE, warning = FALSE, message = FALSE, out.width='30%'----------
edg <- data.frame(from = "x",
                  to = "y",
                  label = "text",
                  label_colour = "blue",
                  label_fill = "red",
                  label_size = 6,
                  label_alpha = .5,
                  label_family = "mono",
                  label_fontface = "bold",
                  label_hjust = "left",
                  label_vjust = "top",
                  label_lineheight = 1.5,
                  label_location = .2
                  )

p = graph_sem(edges = edg, layout = get_layout("x", "y", rows = 1))

if(generate_pngs) ggsave("pgfig7.png", p, device = "png", width= 4, height = 1)
knitr::include_graphics("pgfig7.png")

## ----eval = TRUE, echo = FALSE, warning = FALSE, message = FALSE, out.width='30%'----
fit <- sem("mpg ~ cyl
           mpg ~ am", data = mtcars, meanstructure = TRUE)
get_edges(fit)

## ----eval = FALSE, out.width="300px"------------------------------------------
#  set.seed(6)
#  prepare_graph(fit) %>%
#    color_pos_edges("green") %>%
#    plot()

## ----echo = FALSE, out.width="300px"------------------------------------------
set.seed(6)
prepare_graph(fit) %>%
  color_pos_edges("green") %>%
  plot() -> p
ggsave("pgfig8.png", p, device = "png", width= 4, height = 4)
knitr::include_graphics("pgfig8.png")

## ----eval = FALSE-------------------------------------------------------------
#  prepare_graph(fit) %>%
#    color_pos_edges("green") %>%
#    color_neg_edges("red") %>%
#    color_var("black") %>%
#    alpha_var(.2) %>%
#    plot()

## ----echo = FALSE, out.width="300px"------------------------------------------
set.seed(6)
prepare_graph(fit) %>%
  color_pos_edges("green") %>%
  color_neg_edges("red") %>%
  color_var("black") %>%
  alpha_var(.2) %>%
  plot() -> p
ggsave("pgfig9.png", p, device = "png", width= 4, height = 4)
knitr::include_graphics("pgfig9.png")

## ----eval = FALSE-------------------------------------------------------------
#  prepare_graph(fit) %>%
#    # Add color column to the graph elements
#    edit_graph({ color = "black" }) %>%
#    # Conditionally change color to blue when the column 'est' contains the number 4
#    if_edit(grepl("4", est), {color = "blue"}) %>%
#    plot()

## ----echo = FALSE, out.width="300px"------------------------------------------
set.seed(6)
prepare_graph(fit) %>%
  # Add color column to the graph elements
  edit_graph({ color = "black" }) %>% 
  # Conditionally change color to blue when the column 'est' contains the number 4
  if_edit(grepl("4", est), {color = "blue"}) %>%
  plot()->p
ggsave("pgfig10.png", p, device = "png", width= 4, height = 4)
knitr::include_graphics("pgfig10.png")

## ----eval = FALSE-------------------------------------------------------------
#  model <- "
#    Sepal.Length ~ Sepal.Width + Petal.Length
#    Sepal.Width ~ Petal.Length
#  "
#  # fit model
#  fit <- sem(model, data = iris)
#  # specify layout for consistency
#  layout <- get_layout("", "Sepal.Width", "",
#                       "Petal.Length", "", "Sepal.Length", rows = 2)
#  # get data from prepare_graph
#  p <- prepare_graph(fit, layout = layout, angle = 180)
#  
#  # standard graph
#  plot(p)
#  
#  # Duplicate node data.frame
#  df_nodes <- p$nodes
#  # Add mathematical notation to node label
#  df_nodes$label <- paste("atop(", p$nodes$label, ", ", c("alpha-div", # Add a Greek letter
#                                                   paste0("R^2 ==", formatC(inspect(fit, "r2"), digits = 2, format = "f"))), ")")  # Add R2 to node labels
#  # Set original labels to blank
#  p$nodes$label <- ""
#  
#  # Plot, treat as ggplot object and add parsed node labels
#  plot(p) + geom_text(data = df_nodes, aes(x=x, y=y, label=label), parse = TRUE)
#  

## ----echo = FALSE, out.width="300px"------------------------------------------
set.seed(6)
model <- "
  Sepal.Length ~ Sepal.Width + Petal.Length
  Sepal.Width ~ Petal.Length
"
# fit model
fit <- sem(model, data = iris)
# specify layout for consistency
layout <- get_layout("", "Sepal.Width", "",
                     "Petal.Length", "", "Sepal.Length", rows = 2)
# get data from prepare_graph
p <- prepare_graph(fit, layout = layout, angle = 180)

# standard graph
# plot(p)

# Duplicate node data.frame
df_nodes <- p$nodes
# Add mathematical notation to node label
df_nodes$label <- paste("atop(", p$nodes$label, ", ", c("alpha-div", # Add a Greek letter
                                                 paste0("R^2 ==", formatC(inspect(fit, "r2"), digits = 2, format = "f"))), ")")  # Add R2 to node labels
# Set original labels to blank
p$nodes$label <- "" 

# Plot, treat as ggplot object and add parsed node labels
p <- plot(p) + geom_text(data = df_nodes, aes(x=x, y=y, label=label), parse = TRUE)
ggsave("pgfig11.png", p, device = "png", width= 5, height = 5)
knitr::include_graphics("pgfig11.png")

