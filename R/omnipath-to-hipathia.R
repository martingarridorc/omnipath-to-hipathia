#' Create hiPathia metaginfo from Omnipath interactions
#'
#' Creates an hipathia metaginfo object from a set of omnipath-formatted interactions. 
#' Uses a temporary directory to write the sif and att files required by the \link[hipathia]{mgi_from_sif}
#' function.
#'
#' @param omnipathInteractions A set of interactions generated with the OmnipathR package.
#' @param db The database used to transform gene symbols into entrez ids, which are required by hipathia.
#'
#' @return A metaginfo object ready to use with the \link[hipathia]{hipathia} function.
#' 
#' @export
#' 
#' @importFrom dplyr mutate select group_by summarise ungroup transmute case_when
#' @importFrom OmnipathR interaction_graph
#' @importFrom igraph V E induced_subgraph
#'
omnipathToHipathia <- function(omnipathInteractions, db = org.Hs.eg.db, addGraphVizLayout = TRUE) {
  
  # create igraph object from interactions (for future graph-level modifications)
  g <- OmnipathR::interaction_graph(omnipathInteractions)
  nodes <- igraph::V(g)$name
  # create attributes data frame from graph in hipathia-ready format
  att <- data.frame(symbol = nodes) %>%
    dplyr::mutate(symbolSplitted = strsplit(symbol, "_")) %>%
    tidyr::unnest(cols = symbolSplitted) %>%
    dplyr::mutate(entrez = AnnotationDbi::mapIds(x = db, keys = symbolSplitted, keytype = "SYMBOL", column = "ENTREZID")) %>%
    subset(!is.na(entrez)) %>%
    dplyr::select(symbol, entrez) %>%
    dplyr::group_by(symbol) %>% 
    dplyr::summarise(entrez = paste(entrez, collapse = ",/,")) %>% 
    dplyr::ungroup() %>%
    dplyr::mutate(id = paste0("N-hsa00-", 1:nrow(.))) %>%
    dplyr::transmute(ID = id, label = symbol, X = 0, Y = 0, 
                     color = "white", shape = "rectangle", type = "gene", 
                     label.cex = 0.5, label.color = "black", 
                     width = 46, height = 17,	genesList = entrez,	
                     tooltip = "A")
  # create label to node mapping vector
  labelToNode <- att$ID
  names(labelToNode) <- att$label
  # subset graph to nodes in att
  g <- igraph::induced_subgraph(g, vids = igraph::V(g)$name[igraph::V(g)$name %in% att$label])
  message("Parsing a graph with ", length(igraph::V(g)), " nodes and ", length(igraph::E(g)), " edges.")
  # create sif file in hipathia format
  ints <- igraph::as_data_frame(x = g)
  sif <- dplyr::transmute(ints, 
                          from = from, 
                          type = dplyr::case_when(consensus_stimulation == 1 ~ "activation",
                                                  consensus_inhibition == 1 ~ "inhibition"),
                          to = to) %>%
    dplyr::mutate(from = labelToNode[from], type = type, to = labelToNode[to] ) %>%
    subset(complete.cases(.))
  # write outputs to temporal directory and use mgi_from_sif to reconstruct new metaginfo
  tempDir <- tempdir()
  message("Writing on ", tempDir)
  # write sif and att
  write.table(x = att, file = file.path(tempDir, "hsa00.att"), sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(x = sif, file = file.path(tempDir, "hsa00.sif"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  # write pathway index
  ind <- data.frame(id = "\"00\"", path = "Omnipath")
  write.table(x = ind, file = file.path(tempDir, "name.pathways_hsa.txt"), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  # read as metaginfo
  message("Starting hipathia mgi creation from sif...")
  metaginfo <- hipathia::mgi_from_sif(tempDir, spe = "hsa")
  message("Done!")
  return(metaginfo)
  
}

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
graphVizLayout <- function(g, mode = "TB") {
  
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




