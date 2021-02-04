


# Other functions ---------------------------------------------------------

# @description: compute the inices of the n lowest number in a vector
#               used to compute the indices of the n lowest p-values
#
# @parameters: 
#   - x: a numeric vector
#   - n: the number of indices with lowest values in x to return
#
# @return:
#   - indices of the n lowest number in the vector x
which.minn <- function(x,n=1){
  if (n==1)
    which.min(x)
  else
  {
    if (n>1){
      ii <- order(x,decreasing=FALSE)[1:min(n,length(x))]
      ii[!is.na(x[ii])]
    }
    else {
      stop("n must be >=1")
    }
  }
}

which.maxn <- function(x,n=1){
  if (n==1)
    which.max(x)
  else
  {
    if (n>1){
      ii <- order(x,decreasing=T)[1:min(n,length(x))]
      ii[!is.na(x[ii])]
    }
    else {
      stop("n must be >=1")
    }
  }
}

replace_neg_values <- function(x){
  x[x<0] <- 0 
  return(x)
} 
