#' Title
#'
#' @param adjlist
#' @param directed
#' @param selfloops
#'
#' @return
#' @export
#'
#' @examples
CreateIgGraphs <- function(adjlist, directed, selfloops){
  if(directed)
    mode <- 'directed'
  if(!directed)
    mode <- 'undirected'

  lapply(X = adjlist, FUN = graph_from_adjacency_matrix, mode=mode, diag=selfloops)
}
