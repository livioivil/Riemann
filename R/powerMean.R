##' @description Poer Mean algorithm, as published on IEEE Transactions on Signal Processing (2017), 65 (9) : 2211 - 2220.
##' \link{https://hal.archives-ouvertes.fr/hal-01500514/document}
##' @title powerMean_PDM
##' @export
##' @author Livio Finos, Marco Congedo

# Mat_list=list(matrix(sample(9),3,3),matrix(sample(9),3,3),matrix(sample(9),3,3))
# for(i in 1:3)
#   Mat_list[[i]]=Mat_list[[i]]%*%t(Mat_list[[i]])

powerMean_PDM <- function(Mat_list,p,w=NULL,zeta=1E-5,initialX=NULL){
  
  if(p==0){
    Pplus=powerMean_PDM(Mat_list,p=0.1,w=w,zeta=zeta,initialX=initialX)
    Pminus=powerMean_PDM(Mat_list,p=-0.01,w=w,zeta=zeta,initialX=initialX)
    #modificare: seconda call con initialX dalla prima call
    P=(Pplus+Pminus)/2
    attr(P,"errs")=cbind(attr(Pplus,"errs"),attr(Pminus,"errs"))
  } 
    
  N=ncol(Mat_list[[1]])
  
  phi=0.375/abs(p) 
  # p(-1, 1)\{0}, K positive weights w={w1,…,wK} such that
  # Σkwk=1 
  if(is.null(w)) w=rep(1,length(Mat_list))
  w=w/sum(w)
  
  # C*:
  if(p<0) Mat_list=lapply(Mat_list,MatrixToPower,-1)
  
  # initialize CHECK IT!!!
  # Initialize X as the principal square root inverse of (13) if p(0,1] or as
# its principal square root if p[-1,0)
  
  if(is.null(X)){
    X <- w[1]*MatrixToPower(Mat_list[[1]],p)
    invisible(lapply(2:length(Mat_list),
           function(i){ 
      X <<- X + w[i]*MatrixToPower(Mat_list[[i]],p)}))
    # X
    X=MatrixToPower(X,1/p)  
    
    X=MatrixToPower(X,-.5*sign(p))
  } 
  
##################
  check=zeta*2
  step=1
  while(check[step]>=zeta){  
    
    #compute H
    H=X%*%Mat_list[[1]]%*%t(X)
    invisible(
      lapply(2:length(Mat_list),function(i){ 
        XCXto_p=MatrixToPower(X%*%Mat_list[[i]]%*%t(X),abs(p))
        H <<- H + w[i]*XCXto_p
        }))
    # H
    
    # compute X
    X <- MatrixToPower(H,-phi)%*%X
    
    #check
    check <- c(check, 1/sqrt(N)*norm(H-diag(N),type="F"))
    step=step+1
  }
  
  if(p>0){
    P=solve(X)%*%t(solve(X))
  } else { #p<0
    P=t(X)%*%X
  }
  attr(P,"errs")=check
    P
}

H=(Mat_list[[1]])
invisible(
  lapply(2:length(Mat_list),function(i){ 
    H <<- H + w[i]*(Mat_list[[i]])
  }))
H-powerMean_PDM(Mat_list ,1,zeta=1E-7)


## mi pare non funzioni

H=solve(Mat_list[[1]])
invisible(
  lapply(2:length(Mat_list),function(i){ 
    H <<- H + w[i]*solve(Mat_list[[i]])
  }))
H-powerMean_PDM(Mat_list ,-1,zeta=1E-8)


powerMean_PDM(Mat_list ,0,zeta=1E-8)
