get_node_state_probs<-function(x) {
  hisse.results <- x
  
  
  ConvertToBinaryState <- function(x) {
  x.trimmed <- x[,-1]
  state.0.indices <- which(grepl("0", colnames(x.trimmed)))
  state.1.indices <- which(grepl("1", colnames(x.trimmed)))
  x0.trimmed <- x.trimmed[, state.0.indices]
  if (is.null(dim(x0.trimmed)[1])) { #is a vector
    x0.trimmed <- matrix(data=x0.trimmed, nrow=length(x0.trimmed), ncol=1)
  }
  x1.trimmed <- x.trimmed[, state.1.indices]
  if (is.null(dim(x1.trimmed)[1])) { #is a vector
    x1.trimmed <- matrix(data=x1.trimmed, nrow=length(x1.trimmed), ncol=1)
  }
  result.vector.0 <- apply(x0.trimmed, 1, sum)
  result.vector.1 <- apply(x1.trimmed, 1, sum)
  result.vector <- result.vector.1 / (result.vector.0 + result.vector.1)
  names(result.vector) <- x[,1]
  return(result.vector)
}
  
  ConvertManyToBinaryState <- function(hisse.results, which.element) {
    AIC.weights <- GetAICWeights(hisse.results)
    storage.matrix <- matrix(nrow=dim(hisse.results[[1]][[which.element]])[1], ncol=0)
    for (i in sequence(length(hisse.results))) {
      storage.matrix <- cbind(storage.matrix, ConvertToBinaryState(x=hisse.results[[i]][[which.element]]))	
    }
    final.results <- apply(storage.matrix, 1, weighted.mean, w=AIC.weights)
    return(final.results)
  }
  
  GetAICWeights <- function(hisse.results) {
    AIC.vector <- sapply(hisse.results, "[[", "aic")
    delta.AIC.vector <- AIC.vector - min(AIC.vector)
    rel.likelihood <- exp(-0.5 * delta.AIC.vector)
    AIC.weight.vector <- rel.likelihood / sum(rel.likelihood)
    return(AIC.weight.vector)
  }
  
  
  states.internal <- ConvertManyToBinaryState(hisse.results, "node.mat")
  
  return(states.internal)
}

