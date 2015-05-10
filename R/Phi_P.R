##************************************************************##
##***************** The creiteria function *******************##
##**********************Phi_p criterion **********************##
##***************** Type1: Phi_{p} Criteria ******************##

##********************* Distance Function ********************##
#'distance function
#'@include Phi_P.R
distance<-function(Design,t){
  n<-dim(Design)[1] #get the size of NLHDs
  dist_matrix<-matrix(rep(0,n*n),ncol=n)
  for(i in seq(from=1,to=n)){
    for(j in seq(from=1,to=n)){
      
      dist_matrix[i,j]<-(sum(abs(Design[i,]-Design[j,])^t))^(1/t)
      
    }
  }
  
  return(dist_matrix)
}
###***************** Phi_{p} value **********************###
#'Phi_p function
Phi_p<-function(Design,t,p){
  p<-p
  t<-t
  Design<-Design
  n<-dim(Design)[1] #get the size of GNLHDs
  dist<-distance(Design,t)
  distinverse<-rep(0,n*(n-1)/2)
  for(i in 2:n){
    for(j in 1:(i-1)){
      distinverse[(i-2)*(i-1)/2+j]<-1/dist[i,j]
    }
  }
  Phi_p<-(sum((distinverse)^(p)))^(1/p)
  return(Phi_p)
}
###******************************************************###
#' @export