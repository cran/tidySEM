layout <- get_layout("ne", "phys", "",
                     "plea", "", "",
                     "dist", "", "dep",
                     "saf", "", "",
                     "coh", "stress", "", rows = 5)

df_nodes <- tidySEM:::long_layout(layout)
df_nodes <-  data.frame(node_id = 1:length(df_nodes$name), name = df_nodes$name, stringsAsFactors = FALSE)
df_nodes$shape <- "oval"
df_nodes$shape[grepl("(phys)", df_nodes$name)] <- "rect"
labels <- list("ne" = "Natural environment",
               "plea" = "Pleasantness",
               "dist" = "Disturbance",
               "saf"  = "Safety",
               "coh" = "Cohesion",
               "phys" = "Physical activity",
               "stress" = "Stress",
               "dep" ="Depression"
)
df_nodes$label <- unlist(labels[match(df_nodes$name, names(labels))])

df_edges <- data.frame(matrix(c(
  1, 6, "last", "+",
  2, 6, "last", "+",
  3, 6, "last", "-",
  4, 6, "last", "+",
  5, 6, "last", "+",

  1, 7, "last", "-",
  2, 7, "last", "-",
  3, 7, "last", "+",
  4, 7, "last", "-",
  5, 7, "last", "-",
  6, 8, "last", "-",
  7, 8, "last", "+",
  6, 7, "none", "a",
  1, 2, "none", "b",
  1, 3, "none", "c"), ncol = 4, byrow = TRUE), stringsAsFactors = FALSE)
names(df_edges) <- c("from", "to", "arrow", "label")

df_edges$from <- df_nodes$name[as.numeric(df_edges$from)]
df_edges$to <- df_nodes$name[as.numeric(df_edges$to)]
df_edges$curvature <- c(rep(NA, 12), rep(60, 3))
#df_edges$label[grepl("^ri", df_edges$from)] <- "1"
prep <- prepare_graph(nodes = df_nodes, layout = layout, edges = df_edges, angle = 15, fix_coord = FALSE)

prep$edges[13, c("connect_from", "connect_to")] <- c("bottom", "top")
p <- plot(prep)

test_that("Plot works", {
  expect_s3_class(p, "ggplot")
})
