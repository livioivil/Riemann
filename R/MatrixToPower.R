##' @description The functions compute the p-th power, the exponential and the natural logarithm of the input matrix.
##' @export MatrixToPower
##' @export MatrixExp
##' @export MatrixLog
##' @author Livio Finos, Marco Congedo

MatrixToPower <- function(M,p){
  ei=eigen(M)
  ei$vectors%*%diag(ei$values^p)%*%t(ei$vectors)
}


MatrixExp <- function(M){
  ei=eigen(M)
  ei$vectors%*%diag(exp(ei$values))%*%t(ei$vectors)
}

MatrixLog <- function(M){
  ei=eigen(M)
  ei$vectors%*%diag(log(ei$values))%*%t(ei$vectors)
}
