##' @description Computes the Riemannian distance using the Fisher (affine-invariant) metric between two Symmetric Positive Definite (SPD) matrices.
##' The distance is the length of the geodesic connecting \code{M1} to \code{M2}.
##' @name dist_Riemann
##' @title Riemannian distance 
##' @param M1 a SPD Matrix
##' @param M2 a SPD matrix (same dims of \code{M1}).
##' @param eps smallest eigenvalues. The count of eigenvalues smaller than \code{eps} is reported (just a quality control).
##' @value a distance (non negative scalar)
##' @references \href{Barachant A, Bonnet S, Congedo M, Jutten C (2012) Multi-Class Brain Computer Interface Classification by Riemannian Geometry
##' IEEE Transactions on Biomedical Engineering 59(4), 920-928.}{https://hal.archives-ouvertes.fr/hal-00681328}
##'  (see equation 4)
##' @export
##' @author Livio Finos, Marco Congedo


dist_Riemann <- function(M1,M2,eps=1E-9){
  ei=eigen(M1%*%solve(M2))
  D=sum(log(ei$values)^2)
  attr(D,"small_eigenvalues")=sum(ei$values<eps)
  sqrt(D)
}

#' @description computes the Riemannian norm of a SPD matrix; 
#' it is defined as the distance (\code{dist_Riemann}) between \code{M1} and the identity matrix.
##' @name norm_Riemann
##' @inheritParams dist_Riemann
##' @title Riemannian norm
##' @value a distance (non negative scalar)
##' @references \href{Barachant A, Bonnet S, Congedo M, Jutten C (2012) Multi-Class Brain Computer Interface Classification by Riemannian Geometry
##' IEEE Transactions on Biomedical Engineering 59(4), 920-928.}{https://hal.archives-ouvertes.fr/hal-00681328}
##'  (see equation 4)
##' @export
##' @author Livio Finos, Marco Congedo

norm_Riemann <- function(M1,eps=1E-9){
  ei=eigen(M1)
  D=sum(log(ei$values)^2)
  attr(D,"small_eigenvalues")=sum(ei$values<eps)
  sqrt(D)
}


##' @description Computes the Riemannian geodesic using the Fisher (affine-invariant) metric between two Symmetric Positive Definite (SPD) matrices.
##' \code{mean_geometric_Riemann} is the geodesic with \code{t=.5}.
##' @name geodesic_Riemann
##' @aliases mean_geometric_Riemann
##' @title Riemannian geodesic 
##' @inheritParams dist_Riemann
##' @param t (real number) arc-length of the geodesic from \code{M1} to \code{M2} 
##' (e.g., t=0 yields \code{M1}, t=1 yields \code{M2}, 
##' t=.5 yields the geometric mean of \code{M1} and \code{M2}).
##' @value a distance (non negative scalar)
##' @references \href{Barachant A, Bonnet S, Congedo M, Jutten C (2012) Multi-Class Brain Computer Interface Classification by Riemannian Geometry
##' IEEE Transactions on Biomedical Engineering 59(4), 920-928.}{https://hal.archives-ouvertes.fr/hal-00681328}
##'  (see equation 4)
##' @export
##' @author Livio Finos, Marco Congedo

geodesic_Riemann <- function(M1,M2,t){
  M1sqrt=MatrixToPower(M1,.5)
  M1sqrtinv=solve(M1sqrt)
  M1sqrt%*%MatrixToPower(M1sqrtinv%*%M2%*%M1sqrtinv,t)%*%M1sqrt
  
}

mean_geometric_Riemann <- function(M1,M2){
  geodesic_Riemann(M1,M2,.5)
}



########################
##' @description Computes the Riemannian distance using the Fisher (affine-invariant) metric between two Symmetric Positive Definite (SPD) matrices.
##' \code{norm_Riemann} computes the Riemannian norm of a SPD matrix; it is defined as the distance (\code{dist_Riemann}) between \code{M1} and the identity matrix.
##' The distance is the length of the geodesic connecting \code{M1} to \code{M2}.
##' @aliases exp_Map_Riemann log_Map_Riemann projection_on_tangent_space_Riemann
##' @param M1 a SPD Matrix
##' @param P a SPD matrix (same dims of \code{M1}). Point defining the tangent space.
##' @value \code{exp_Map_Riemann} returns a SPD matrix,  
##' \code{log_Map_Riemann} returns a symmetric matrix, 
##' \code{projection_on_tangent_space_Riemann} returns a vector.
##' @references Barachant A, Bonnet S, Congedo M, Jutten C (2012) Multi-Class Brain Computer Interface Classification by Riemannian Geometry
##' IEEE Transactions on Biomedical Engineering 59(4), 920-928. 
##' https://hal.archives-ouvertes.fr/hal-00681328
##' @export exp_Map_Riemann
##' @export log_Map_Riemann 
##' @export projection_on_tangent_space_Riemann
##' @author Livio Finos, Marco Congedo

exp_Map_Riemann <- function(M1,P){
  Psqrt=MatrixToPower(P,.5)
  Psqrtinv=solve(Psqrt)
  Psqrt%*%MatrixExp(Psqrtinv%*%M1%*%Psqrtinv)%*%Psqrt
  
}

log_Map_Riemann <- function(M1,P){
  Psqrt=MatrixToPower(P,.5)
  Psqrtinv=solve(Psqrt)
  Psqrt%*%MatrixLog(Psqrtinv%*%M1%*%Psqrtinv)%*%Psqrt
}

projection_on_tangent_space_Riemann <- function(M1,P){
  Psqrt=MatrixToPower(P,.5)
  Psqrtinv=solve(Psqrt)
  temp=MatrixLog(Psqrtinv%*%M1%*%Psqrtinv)
  temp=temp*(1-diag(nrow(P)))*sqrt(2)
  temp[upper.tri(temp,diag=TRUE)]
}
