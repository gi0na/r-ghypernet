#' Check graph input type (for whether it's a graph or a edgelist).
#'
#' Returns TRUE if the supplied object \code{graph} is an adjacency matrix. Returns FALSE if the provided object is an edgelist. The function checks whether the edgelist conforms to our standards (sender, target, edgecount).
#'
#' @param graph A graph adjacency matrix or an edgelist.
#' @return TRUE or FALSE. Returns TRUE if the provided object \code{graph} is an adjacency matrix.
checkGraphtype <- function(graph) {
  #TODO: add igraph compatibility
  # is matrix?
  isMatrix <-  FALSE
  if (dim(graph)[1] == dim(graph)[2]) {
    isMatrix <- TRUE
  }
  # if edgelist = check if 3 columns are there, 2 characters/factors, 1 numeric
  if (!isTRUE(isMatrix)) {
    if (ncol(graph) != 3) {
      stop(
        "graph needs to be an adjacecy matrix or an edgelist. The edgelist needs to have exactly 3 columns: sender, target, edgecounts."
      )
    }
    if (!is.character(graph[, 1]) & !is.factor(graph[, 1])) {
      stop(
        "The first row in your edgelist needs to be a factor/character vector with sender node IDs."
      )
    }
    if (!is.character(graph[, 2]) & !is.factor(graph[, 2])) {
      stop(
        "The second row in your edgelist needs to be a factor/character vector with target node IDs."
      )
    }
    if (!is.numeric(graph[, 3]) & !is.integer(graph[, 3])) {
      stop(
        "The third row in your edgelist needs to be a numeric/integer vector with edge counts."
      )
    }
  }
  # return:
  isMatrix
}

################################################################################
#' Calculate weighted reciprocity change statistics for multi-edge graphs.
#'
#' The function takes either an edgelist or an adjacency matrix and returns an
#' adjacency matrix with the reciprocity change statistic. This reciprocity
#' matrix can then be used as a predictor in the gHypEG regression.
#'
#' @param graph A graph adjacency matrix or an edgelist.  The edgelist needs to
#'   have 3 columns: a sender vector, a target vector and an edgecount vector.
#' @param nodes optional character/factor vector. If an edgelist is provied,
#'   you have to provide a list of unique identifiers of your nodes in the graph.
#'   This is because in the edgelist, isolates are usually not recorded.
#'   If you do not specify isolates in your nodes object, they are excluded
#'   from the analysis (falsifies data).
#' @param zero_values optional numeric value. Use this to substitute zero-values
#'   in your reciprocity change statistic matrix. Zero values in the predictors
#'   are recognized in the gHypEG regression as structural zeros. To ensure this
#'   does not happen, please recode your zero-values in all your predictors.
#'   If \code{zero_values} is not specified, the minmal value divided by 10 is
#'   used instead.
#' @return Reciprocity change statistic matrix.
#' @author LB, GC
#' @seealso \code{\link{sharedPartner_stat}} or \code{\link{homophily_stat}}
#' @export
#' @example 
#' 
reciprocity_stat <- function(graph, nodes = NULL, zero_values = NULL){
  ## preprocess:
  # is graph object a matrix or edgelist?
  isMatrix <- checkGraphtype(graph) # returns isMatrix==TRUE/FALSE
  # if not isMatrix transform edgelist into matrix
  if(!isTRUE(isMatrix)){
    
    el <- graph #bc graph will be overwritten
    
    # check if they provided a list of nodes for the graph (matching attribute files)
    if(is.null(nodes)){
      nodes <- unique(c(el[,1], el[,2]))
    }
    #
    
    #transform edgelist to matrix
    # save it to object "graph"
    graph <- el2adj(el, nodes)
    
  }else{
    nodes = rownames(graph)
  }
  
  ## calculate reciprocity
  recip_mat <- t(graph)
  
  ## transform zero values
  if(is.null(zero_values)){
    recip_mat[recip_mat == 0] <- min(recip_mat[recip_mat > 0])/10
  }else{
    recip_mat[recip_mat == 0] <- zero_values
  }
  
  return(recip_mat)
}

################################################################################
#' Calculate (un-)weighted shared partner change statistics for multi-edge graphs.
#'
#' The function calculates the change statistic for shared partners for each
#' dyad in the graph. Shared partner statistics count for each dyad involving
#' nodes i and j in the graph, how many nodes k these two nodes have in common
#' (or share). The shared partner $k$ counts are weighted by their
#' interactions with the focal nodes $i$ and $j$. This is neccessary in
#' dense multi-edge graphs to ensure that meaningful triadic closure is
#' detected. The statistic can be calculated in 3 different forms: undirected,
#' incoming shared partners (where shared partner k: k->i and k->j) and outgoing
#' shared partners (where shared partner k: k<-i and k<-j).
#'
#' @param graph A graph adjacency matrix or an edgelist. The edgelist needs to
#'   have 3 columns: a sender vector, a target vector and an edgecount vector.
#' @param weighted set to TRUE.
#' @param triad.type set to \code{undirected}. Can be set to \code{incoming}
#'   or \code{outgoing} instead. This then corresponds to directd triadic closure
#'   in the multi-edge graph.
#' @param nodes optional character/factor vector. If an edgelist is provied,
#'   you have to provide a list of unique identifiers of your nodes in the graph.
#'   This is because in the edgelist, isolates are usually not recorded.
#'   If you do not specify isolates in your nodes object, they are excluded
#'   from the analysis (falsifies data).
#' @param zero_values optional numeric value. Use this to substitute zero-values
#'   in your shared partner change statistic matrix. Zero values in the predictors
#'   are recognized in the gHypEG regression as structural zeros. To ensure this
#'   does not happen, please recode your zero-values in all your predictors.
#' @param directed boolean. Is the graph directed?
#'   If \code{zero_values} is not specified, the 0.1 is used instead.
#' @return Shared partner change statistic matrix.
#' @author LB, GC, GV
#' @seealso \code{\link{reciprocity_stat}} or \code{\link{homophily_stat}}
#' @export
#' @import dplyr
#' @importFrom utils setTxtProgressBar txtProgressBar
sharedPartner_stat <- function(graph,
                               directed,
                               weighted = TRUE,
                               triad.type = 'undirected',
                               nodes = NULL,
                               zero_values = NULL) {
  ## preprocess:
  # is graph object a matrix or edgelist?
  isMatrix <- checkGraphtype(graph) # returns isMatrix==TRUE/FALSE
  # transform edgelist into matrix
  if (isTRUE(isMatrix)) {
    adj <- graph
    el <- adj2el(adj, directed = directed)
  } else{
    adj <- el2adj(graph, nodes)
  }
  if (is.null(nodes))
    nodes <- rownames(adj)
  
  # check triad.type
  if (!triad.type == 'undirected' &
      !triad.type == 'directed.incoming' &
      !triad.type == 'directed.outgoing') {
    stop(
      'triad.type needs to be one of three: "undirected", "directed.incoming" or "directed.outgoing".
      Check help files for additional infos.'
    )
  }
  
  if (!directed & !isSymmetric(adj))
    adj <- adj + t(adj)
  diag(adj) <- 0
  
  # progressbar
  pb = txtProgressBar(
    min = 1,
    max = nrow(adj) - 1,
    initial = 0,
    style = 3
  )
  
  ## now count triangles
  if (triad.type == 'directed.outgoing')
    partners <- sapply(nodes, findPartners_target, el = el)
  if (triad.type == 'directed.incoming')
    partners <- sapply(nodes, findPartners_sender, el = el)
  if (triad.type == 'undirected' | isFALSE(directed))
    partners <- sapply(nodes, findPartners_all, el = el)
  
  tri <-
    matrix(0, nrow(adj), ncol(adj), dimnames = list(rownames(adj), rownames(adj)))
  
  for (i in nodes[-length(nodes)]) {
    for (k in nodes[(which(nodes == i) + 1):length(nodes)]) {
      js <- intersect(partners[[i]], partners[[k]])
      if(isFALSE(weighted)){
        tri[i, k] <- length(js)
        next
      }
      if (length(js) > 0) {
        if (triad.type == 'directed.outgoing' |
            triad.type == 'undirected') {
          tri[i, k] <-
            tri[i, k] + sum(sapply(js, function(j)
              adj[j, i] * adj[j, k] / max(adj[j, i], adj[j, k], 1)))
        }
        if (directed & (triad.type == 'directed.incoming' |
                        triad.type == 'undirected')) {
          tri[i, k] <-
            tri[i, k] + sum(sapply(js, function(j)
              adj[i, j] * adj[k, j] / max(adj[i, j], adj[k, j], 1)))
        }
      }
    }
    setTxtProgressBar(pb, which(nodes==i))
  }
  tri <- tri + t(tri)
  
  ## transform zero values
  if (!is.null(zero_values)) {
    tri[tri == 0] <- zero_values
  }
  
  return(tri)
  }

findPartners_sender <- function(node, el){
  colnames(el) <- c('sender','target')
  return(el$target[el$sender == node])
}

findPartners_target <- function(node, el){
  colnames(el) <- c('sender','target')
  return(el$sender[el$target == node])
}

findPartners_all <- function(node, el){
  unique(c(findPartners_sender(node,el),findPartners_target(node,el)))
}

################################################################################
## graph statistic: homophily
#Function: Creates a matrix for categorical attribute matches between two nodes
# as well as absolute difference effects.

#' Calculate homophily in multi-edge graphs.
#'
#' The function calculates homophily matrices. If you supply a categorical
#' variable (factor, character), the function returns attribute matches for dyads
#' from the same group. If you supply a continuous variable
#' (numeric, integers), the function returns absolute difference effects for
#' each dyad in the graph.
#'
#' @param variable A attribute variable. Can be categorical (attribute matches) or
#'   continuous (absolute difference effects).
#' @param type set to \code{categorical}. Can be set to \code{absdiff} instead.
#'   If set to \code{categorical}, the homophily statistic calculates matches
#'   between dyads from the same group (analogous to dummy variables measuring
#'   attribute match between two nodes (=10) and attribute mismatch (=1)). If
#'   set to \code{absdiff} it calculates the difference in values from variable
#'   for each dyad in the graph.
#' @param nodes optional character/factor vector. If an edgelist is provied,
#'   you have to provide a list of unique identifiers of your nodes in the graph.
#'   This is because in the edgelist, isolates are usually not recorded.
#'   If you do not specify isolates in your nodes object, they are excluded
#'   from the analysis (falsifies data).
#' @param zero_values optional numeric value. Use this to substitute zero-values
#'   in your reciprocity change statistic matrix. Zero values in the predictors
#'   are recognized in the gHypEG regression as structural zeros. To ensure this
#'   does not happen, please recode your zero-values in all your predictors.
#'   If \code{zero_values} is not specified, the minmal value divided by 10 is
#'   used instead.
#' @author LB, GC
#' @seealso \code{\link{reciprocity_stat}} or \code{\link{sharedPartner_stat}}
#' @export
homophily_stat <- function(variable = variable,
                           type = 'categorical',
                           #type = categorical, absdiff
                           nodes = nodes,
                           these.categories.only = NULL,
                           zero_values = NULL) {
  # check inputs
  if (is.null(variable)) {
    stop(
      "You need to specify a homophily variable. Either a categorical (=character, factor) or continous (=numeric, integer) variable."
    )
  } else{
    if (is.character(variable) | is.factor(variable)) {
      # check type
      if (type == 'absdiff') {
        stop("Please use type = categorical for chacacter/factor homophily variables.")
      }
      # check categories
      if (!is.null(these.categories.only)) {
        if(!all(these.categories.only %in% unique(variable))){
          stop(
            "The values you provided in these.categories.only are not subsets of the variable object."
          )
        }
      }
      # check length/number of categories
      if (length(unique(variable)) > 50) {
	      ##TODO LB: wait, I think I got this wrong here?
        stop("Please recode your homophily variable. It cannot have more than 20 categories.")
      }
    }
  }
  if (is.null(nodes)) {
    stop(
      "You need to specify a nodes object. This object lists the nodes in the graph (e.g., their ID)."
    )
  }
  if (type != 'categorical' & type != 'absdiff') {
    stop("Homophily type can either be categorical or absolute difference")
  }
  
  
  ### Categorical homophily
  if (type == 'categorical') {
    ## create block IDs
    labels <- levels(factor(variable))
    blockids <- as.numeric(plyr::mapvalues(
      variable,
      from = labels,
      to = c(1, numbers::Primes(length(variable) * 50))[1:(length(unique(variable)))]
    ))
    #TODO: package dependence here!!
    #verbose = TRUE:
    #print(unique(blockids))
    #print(unique(labels))
    ## create block matrix
    blocks <- blockids %*% t(blockids)
    
    ## now sort the blocks (according to which values are specified in these.categories.only)
    if (is.null(these.categories.only)) {
      # create homophily matrix
      homophily_mat <-
        ifelse(matrix(blocks %in% blockids ^ 2, nrow = length(nodes)), 10, 1) #hardcoded 10, 1
    } else{
      #these.categories.only = specified
      # which labels should be selected?
      selectedIDs <-
        numbers::Primes(length(variable))[labels %in% these.categories.only]
      homophily_mat <-
        ifelse(matrix(blocks %in% selectedIDs ^ 2, nrow = length(nodes)), 10, 1) #hardcoded 10, 1
    }
    ## label matrix
    rownames(homophily_mat) <- nodes
    colnames(homophily_mat) <- nodes
    
  } else{
    # type = absdiff
    ### Absolute difference
    homophily_mat <-
      matrix(0, nrow = length(nodes), ncol = length(nodes))
    rownames(homophily_mat) <- nodes
    colnames(homophily_mat) <- nodes
    ##TODO: this is super inefficient (2-for loops), move this to Rcpp
    for (i in 1:nrow(homophily_mat)) {
      for (j in 1:ncol(homophily_mat)) {
        homophily_mat[i, j] <-
          abs(variable[nodes %in% nodes[i]] - variable[nodes %in% nodes[j]])
      }
    }
  }
  
  ## treat zero-values
  if (is.null(zero_values)) {
    homophily_mat[homophily_mat == 0] <-
      min(homophily_mat[homophily_mat > 0]) / 10
  } else{
    homophily_mat[homophily_mat == 0] <- zero_values
  }
  
  ## return matrix with 0.1/1 or absolute difference
  return(homophily_mat)
}
