
### uniform abundance
uniform <- function(vect){
  vect/sum(vect)  
}

### species counting
numberof <- function(vect){
  length(which(vect > 0))
}

### shannon statistical
shannon <- function(vect){
  vect1 <- vect/sum(vect)
  vect1 <- vect1[vect1 > 0]
  vect2 <- -log2(vect1)
  vect3 <- vect2 * vect1
  sum(vect3)
}

### merge low abundance 
merge_low_abundance <- function(x, vector_name){
  others_ind <- order(-x)[-(1:top)]
  others <- sum(x[others_ind])
  x[others_ind] <- 0
  x <- c(x, others = others)
  x
}

###remove value==0 and sort by rowSum
rm_sort <- function(data){
  data <- data[which(rowSums(data) > 0),]
  table <- as.matrix(data[order(-rowSums(data)), ])
  list(data,table)
}


###Statistical Frac dist
dist.Frac <- function(inMatrix, ...) {
  unifrac <- function(x, y){
    1 - length(which(x>0 & y>0))/length(which(x>0 | y>0))
  }
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=unifrac(as.vector(inMatrix[,i]),
                                 as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

###Statistical pam cluster
### x is a distance matrix and k the number of clusters
pam.clustering=function(x,k) { 
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

