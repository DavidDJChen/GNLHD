##****************************************************************************##
##********************Draw the density of GNLHD with structure ***************##
#'Phi_p_distribution function
#'@include Plot_Grid.R
Phi_p_distribution<-function(structure,dimension,number,p_value,t_value){
  structure<-structure
  dimension<-dimension
  p_value<-p_value
  t_value<-t_value
  u2<-GNLHD$new(s=structure,q=dimension)
  phi_p2<-rep(0,number)
  for(i in 1:number){
     a1<-u2$StandGNLHD()
     phi_p2[i]<-Phi_p(a1,t=t_value,p=p_value)
  }
  hist(phi_p2,freq=FALSE,main="s=structure,p=p_value and q=q_value",xlab="Phi_p value")   
}
##*************************************************************************************##
#'satter_grids funtion
##********************** Draw grids on scatter-plot ***********************************##
scatter_grids<-function(Design_Full,structure,lcm,k,wid){
  m<-length(structure)
  lcm<-lcm
  wid<-wid
  k<-k
  s<-structure
  Design_Full<-Design_Full
  X<-list()
  if(m==1){
    a1<-rep(0,s[m]*2)
    X<-matrix(a1,ncol=2)
    X[,1]<-(Design_Full[1:s[m],1]-0.5)/lcm
    X[,2]<-(Design_Full[1:s[m],2]-0.5)/lcm
  }else{
     X[[1]]<-matrix(rep(0,s[1]*2),ncol=2)
     X[[1]][1:s[1],1]<-(Design_Full[1:s[1],1]-0.5)/lcm
     X[[1]][1:s[1],2]<-(Design_Full[1:s[1],2]-0.5)/lcm
     for(i in 2:m){
       a<-rep(0,(s[i]-s[i-1])*2)
       X[[i]]<-matrix(a,ncol=2)
       X[[i]][1:(s[i]-s[i-1]),1]<-(Design_Full[(s[i-1]+1):s[i],1]-0.5)/lcm
       X[[i]][1:(s[i]-s[i-1]),2]<-(Design_Full[(s[i-1]+1):s[i],2]-0.5)/lcm
     }
  }
  for(i in 1:m){
    plot(X[[i]][,1],X[[i]][,2],xlim=c(0.04,0.962),ylim=c(0.04,0.962),xlab="x1",ylab="x2",pch=k[i])
    par(new=TRUE)
  }
  if(lcm==s[m]){
    for(i in 1:m){
       grid(s[i],s[i],col="gray",lwd=wid[i],lty=3)
    }
  }else{
    for(i in 1:m){
       grid(s[i],s[i],col="black",lwd=wid[i],lty=3)
    }
       grid(lcm,lcm,col="gray",lwd=wid[m+1],lty=3)
  }
}
##************************************************************##
#'@export