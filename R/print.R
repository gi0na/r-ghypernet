### print
#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.ghype <- function(x, ...){
  directed <- 'undirected'
  if(x$directed)
    directed <- 'directed'
  selfloops <- 'no selfloops'
  if(x$selfloops)
    selfloops <- 'selfloops'
  out <- paste('ghypne',directed,',',selfloops,'\n')
  cat(out)
  out <- paste(x$n[1], 'vertices,', x$m, 'edges','\n')
  cat(out)
  cat('Loglikelihood:\n')
  cat(x$loglikelihood)
}
