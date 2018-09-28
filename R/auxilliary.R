updateModel <- function(model, adj){

  callname <- as.character(model$call[1])
  fixXi <- length(grep('xi', deparse(model$call, width.cutoff = 500)))>0
  xi <- NULL
  if(fixXi)
    xi <- model$xi

  newcall <- NULL
  if(length(grep('ghype', callname))>0){
    newcall <- call(name = callname, object=adj, directed=model$directed, selfloops=model$selfloops, xi=xi, unbiased=all(model$omega==1))
  } else{
    if(length(grep('Block', callname))>0){
      newcall <- call(name = callname, adj=adj, labels=model$labels, directed=model$directed, selfloops=model$selfloops, xi=xi)
    }
  }
  return(newcall)
}
