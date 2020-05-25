#' Create graphViz layout from graph
#' 
#' Creates a matrix containing the graphViz layout for an igraph 
#' object.
#'
#' @param g The graph to process.
#' @param mode The rank direction passed to the \link[Rgraphviz]{layoutGraph} function.
#'
#' @return A matrix with two columns containing the graph layout.
#' 
#' @export
#' 
#' @importFrom graph graphNEL
#' @importFrom Rgraphviz layoutGraph
#' @importFrom igraph get.data.frame
#' 
graphVizLayout <- function(g, mode = "LR") {
  
  # get simple graph format
  simple <- igraph::get.data.frame(g)[,1:2]
  # create graphNEL
  rEG <- new("graphNEL", nodes = V(g)$name, edgemode = "directed")
  for(i in 1:nrow(simple)){
    rEG <- graph::addEdge(simple[i,1], simple[i,2], rEG, 1)
  }
  # set graph atts
  att <- list(graph = list(rankdir = mode))
  # render layout
  rEG <- Rgraphviz::layoutGraph(rEG, attrs = att)
  # get matrix with X and Y
  l <- cbind(rEG@renderInfo@nodes$nodeX, y = rEG@renderInfo@nodes$nodeY)
  colnames(l) <- c("x","y")
  return(l)
  
}

#' Beautiful hipathia graph plot
#' 
#' Creates a plot from the hipathia graph object using \link[Rgraphviz]{layoutGraph} and \link[ggraph]{ggraph}.
#'
#' @param g The graph to plot.
#'
#' @return A ggplot2 object containing the plot.
#' 
#' @export
#' 
#' @importFrom igraph V E remove.vertex.attribute induced_subgraph
#' @importFrom dplyr %>% case_when
#' @importFrom grid arrow
#' @importFrom ggraph ggraph geom_edge_link scale_edge_color_manual geom_node_label
#' @importFrom ggplot2 theme
#' 
#'
beautyHipathiaGraph <- function(g) {
  
  if(!is.null(igraph::V(g)$x)) {
    g <- igraph::remove.vertex.attribute(g, "x") %>% 
      igraph::remove.vertex.attribute(graph = ., name = "y")
  }
  # remove func nodes
  notFuncNodes <- V(g)[!grepl("_func", V(g)$name)]
  g <- igraph::induced_subgraph(g, vids = notFuncNodes)
  # prepare colors
  direction <- dplyr::case_when(igraph::E(g)$relation == 1 ~ "Activation", TRUE ~ "Inhibition")
  # create labels replacing complexs by ...
  labs <- igraph::V(g)$label %>%
    gsub(pattern = " .*", replacement = " (C)", x = .)
  # plot and return graph
  p <- ggraph::ggraph(g, layout = graphVizLayout(g, mode = "LR")) +
    ggraph::geom_edge_diagonal(aes(color = direction), 
                               arrow = grid::arrow(length = unit(2, 'mm'), type = "closed"), 
                               start_cap = rectangle(8, 5, 'mm'),
                               end_cap = rectangle(10, 5, 'mm')) +
    ggraph::geom_node_label(aes(label = labs), size = 2) +
    ggraph::scale_edge_color_manual(values = c("Activation"="red", "Inhibition"="blue")) +
    ggplot2::theme(panel.background = element_blank(), legend.position = "bottom")
  return(p)
  
}
