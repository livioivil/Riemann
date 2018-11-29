##' @description Validity check (and fix) for files path
##' @export
##' @author Livio Finos, Marco Congedo

MatrixToPower <- function(M,p){
  ei=eigen(M)
  ei$vectors%*%diag(ei$values^p)%*%t(ei$vectors)
}
