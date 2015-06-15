##*******************************************************************##
##****************** Project Name: Generating GNLH ******************##
##***************** Description: Generating the GNLH ****************##
##*********** Author: Daijun Chen ** Date: April-8th-2015 ***********##  
##*******************************************************************##
library(R6)
library(numbers)
##********************* Part 1: Basic Class---GNLHD *****************##
#'GGNLHD class
#'@include GNLHD_Construction.R
GGNLHD<-R6Class("GGNLHD",
                public=list(
                  s=NA,    #the structure of GNLHD
                  q=NA,    #the dimension of GNLHD
                  Lcm=NA,  #the least common multiplier of s
                  t=NA,    #the Lcm/n_i 
                  initialize=function(s,q){    #initialization
                    self$s<-s
                    self$q<-q
                    if(length(self$s)==1){ 
                      self$Lcm<-self$s
                    }else self$Lcm<-mLCM(self$s)
                    self$t<-self$Lcm/self$s
                  }    
                )
)

##******************************************************************##
#'LHD class
##******************** Part 2: Sub-Class 1---LHD *******************##
LHD<-R6Class("LHD",
             inherit=GGNLHD,
             public=list(
               LH=function(){#Latin Hypercube Sampling Method
                 num<-self$s #number of samples
                 dim<-self$q #dimension of samples
                 LH<-matrix(rep(0,num*dim),ncol=dim)
                 for(i in 1:dim){ #generating LH
                   LH[,i]<-sample(seq(1,num),num)
                 }
                 return(LH)
               },
               Swap=function(LH,column){#swap operation on LH
                 LH_old<-LH
                 LH_new<-LH_old
                 col<-column
                 n<-length(LH_old[,col])
                 if(n==1){ #special case: n==1
                   LH_new<-LH_old
                 }else{ 
                   swap_position<-sample(n,2)
                   i<-swap_position[1]
                   j<-swap_position[2]
                   LH_new[i,col]<-LH_old[j,col]
                   LH_new[j,col]<-LH_old[i,col]
                 }
                 return(LH_new)
               },
               StandLHD=function(){
                 n<-self$s[length(self$s)]#number of samples
                 StandLHD<-(self$LH()-runif(1,min=0,max=1))/n
                 return(StandLHD)
               }
             )
)
##*******************************************************************##
#'NLHD class
##******************* Part 3: Sub-Class 2---NLHD ********************##
NLHD<-R6Class("NLHD",
              inherit=GGNLHD,
              public=list(
                ##************************ NLH_permutation **************************##
                NLH_permutation=function(){#generating NLH 
                  n<-self$s[length(self$s)]#number of samples
                  t<-self$t
                  s<-self$s
                  t[length(s)+1]=1
                  for(i in 1:length(s)){                    
                    if(i==1){
                      permnp<-rep(0,len=n)
                      Z<-seq(from=1,to=s[i])
                      perm<-sample(Z)
                      for(j in seq(1,s[1])){
                        permnp[j]<-sample(seq(from=(perm[j]-1)*t[i]+1,
                                              to=(perm[j])*t[i]),1)    
                      }
                      C<-ceiling(permnp[1:s[i]]/t[2])
                    }
                    else{
                      if(i!=length(s)){
                        Z<-seq(from=1,to=s[i])
                        perm<-sample(setdiff(Z,C))    
                        for(j in seq(s[i-1]+1,s[i])){
                          permnp[j]<-sample(seq(from=(perm[j-s[i-1]]-1)*t[i]+1,
                                                to=(perm[j-s[i-1]])*t[i]),size=1)
                        }
                        C<-ceiling(permnp[1:s[i]]/t[i+1])
                      }
                      else{
                        Z<-seq(from=1,to=s[i])
                        perm<-sample(setdiff(Z,C))    
                        for(j in seq(s[i-1]+1,s[i])){
                          permnp[j]<-seq(from=(perm[j-s[i-1]]-1)*t[i]+1,
                                         to=(perm[j-s[i-1]])*t[i])
                        }
                        C<-ceiling(permnp[1:s[i]]/t[i+1]) 
                      }    
                    }
                    
                  }
                  return(permnp)
                },
                ##********************************* NLH *********************************##
                NLH=function(){#Nested Latin Hypercube Sampling
                  num<-self$s[length(self$s)]#number of samples
                  dim<-self$q #dimension of samples
                  NLH<-matrix(rep(0,num*dim),ncol=dim)
                  for(i in 1:dim){#generating NLH
                    NLH[,i]<-self$NLH_permutation()                     
                  }
                  return(NLH)
                },
                ##*************************** Swap operation ****************************##
                Swap=function(NLH,structure,column,Swap_type,Swap_layer){
                  NLH_old<-NLH
                  NLH_new<-NLH_old
                  col<-column
                  type<-Swap_type
                  layer<-Swap_layer
                  s<-structure
                  m<-length(s)
                  legal_set<-function(candidate,NLH_old_col,old_row,structure,layer){
                    NLH_old_col<-NLH_old_col
                    candidate<-candidate
                    len<-length(candidate)
                    legal_set<-rep(0,len)
                    s<-structure
                    m<-length(s)
                    t<-s[m]/s
                    layer<-layer
                    old_row<-old_row
                    for(i in 1:len){
                      candidatei<-candidate[i]
                      positioni<-which(NLH_old_col==candidatei)
                      for(j in 2:length(s)){
                        if((positioni>s[j-1])&&(positioni<=s[j])){
                          positionlayeri<-j
                        }  
                      }
                      judge<-rep(NA,positionlayeri-layer)
                      p=1
                      for(k in (layer):(positionlayeri-1)){
                        set_layer<-ceiling(setdiff(NLH_old_col[1:s[k]],old_row)/t[k])
                        if(length(setdiff(set_layer,ceiling(candidatei/t[k])))==length(set_layer)){
                          judge[p]<-TRUE
                        }else{
                          judge[p]<-FALSE
                        }
                        p=p+1
                      }
                      if(all(judge==TRUE)){
                        legal_set[i]=i
                      }else{
                        legal_set[i]=FALSE
                      }
                      
                    }
                    legal_set<-candidate[legal_set]
                    return(legal_set)
                  }
                  if(layer==1){
                    num<-s[layer]
                  }else{
                    num<-s[layer]-s[layer-1]
                  }
                  if(type=="in"){
                    do1.this<-"T1" # turn to layer_in swap
                  }else{
                    do1.this<-"T2" # turn to layer_between swap
                  }
                  if(layer!=m){
                    do2.this<-"T3"
                  }else{
                    do2.this<-"T4"
                  }
                  switch(do1.this,
                         T1={
                           if(num==1){
                             NLH_new<-NLH_old
                           }else{
                             exchange<-sample(num,2)
                             i<-exchange[1]
                             j<-exchange[2]
                             if(layer==1){
                               NLH_new[i,col]=NLH_old[j,col]
                               NLH_new[j,col]=NLH_old[i,col]
                             }else{
                               start<-s[layer-1]
                               starti<-start+i
                               startj<-start+j
                               NLH_new[starti,col]<-NLH_old[startj,col]
                               NLH_new[startj,col]<-NLH_old[starti,col]
                             }
                           }     
                         },
                         T2={
                           switch(do2.this,
                                  T3={
                                    if(layer==1){
                                      i<-sample(num,1)
                                      element<-NLH_old[i,col]
                                      t<-s[m]/s[layer]
                                      tag<-ceiling(element/t)
                                      domain<-seq(from=(tag-1)*t+1,to=tag*t)
                                      candidate<-setdiff(domain,element)
                                      legal_set<-legal_set(candidate,NLH_old[,col],element,s,layer)
                                      if(length(legal_set)==1){
                                        exchange<-legal_set
                                      }else{
                                        exchange<-sample(legal_set,1)
                                      }
                                      j<-which(NLH_old[,col]==exchange)
                                      NLH_new[i,col]=NLH_old[j,col]
                                      NLH_new[j,col]=NLH_old[i,col]
                                    }else{
                                      i<-sample(num,1)+s[layer-1]
                                      element<-NLH_old[i,col]
                                      t<-s[m]/s[layer]
                                      tag<-ceiling(element/t)
                                      domain<-seq(from=(tag-1)*t+1,to=tag*t)
                                      candidate<-setdiff(domain,element)
                                      legal_set<-legal_set(candidate,NLH_old[,col],element,s,layer)
                                      if(length(legal_set)==1){
                                        exchange<-legal_set
                                      }else{
                                        exchange<-sample(legal_set,1)
                                      }
                                      j<-which(NLH_old[,col]==exchange)
                                      NLH_new[i,col]=NLH_old[j,col]
                                      NLH_new[j,col]=NLH_old[i,col]                                               
                                    }
                                  },
                                  T4={
                                    if(num==1){
                                      NLH_new<-NLH_old
                                    }else{
                                      exchange<-sample(num,2)
                                      i<-exchange[1]
                                      j<-exchange[2]
                                      start<-s[layer-1]
                                      starti<-start+i
                                      startj<-start+j
                                      NLH_new[starti,col]<-NLH_old[startj,col]
                                      NLH_new[startj,col]<-NLH_old[starti,col]
                                    }
                                  }
                           )
                         }
                  )
                  return(NLH_new)
                },
                ##************************ StandNLHD ************************##
                StandNLHD=function(){
                  n<-self$s[length(self$s)]#number of samples
                  StandNLHD<-(self$NLH()-0.5)/n
                  return(StandNLHD)
                }
              )  
)

##************************************************************##
#'GNLHD class
##************** Part 4: Sub-Class 3---GNLHD ***************##
GNLHD<-R6Class("GNLHD",
               inherit=GGNLHD,
               public=list(
                 GNLH_illegal_set=function(){ #find out
                   s<-self$s  #the illegal sets in each 
                   lcm<-self$Lcm#layer
                   t<-self$t
                   n<-length(self$s)
                   illegal_set<-list()
                   for(i in seq(from=1,to=n-1)){
                     if(s[i+1]%%s[i]==0){
                       illegal_set[[i]]<-matrix(rep(0,t[i+1]),nrow=1)
                     }else{
                       Z<-seq(from=1,to=lcm)
                       X<-matrix(Z,ncol=t[i+1],byrow=TRUE)
                       is_bad<-ceiling(X/t[i])
                       bad<-rep(NA,n=s[i+1])
                       for(j in 1:s[i+1]){
                         if(length(unique(is_bad[j,]))!=1){
                           bad[j]<-1
                         }else bad[j]<-0
                       }
                       X_is_bad<-matrix(X[as.vector(which(bad==1)),],ncol=t[i+1])
                       illegal_set[[i]]<-X_is_bad
                     }
                   }
                   illegal_set[[n]]<-matrix(rep(0,t[n]),nrow=1)#do not judge the last layer
                   return(illegal_set) 
                 },
                 ##******************** Switch to exclude the illegal_set ************##
                 diagnose=function(i,illegal_set,k,permnp,C){##diagnose the ith layer
                   n<-self$s[length(self$s)]
                   s<-self$s #the structure of GNLH
                   lcm<-self$Lcm #least common multiplier
                   m<-length(s) # total m layers
                   t<-self$t #the multipliers of each layer
                   t[length(s)+1]=1
                   C<-C
                   if(i==1){
                     do1.this="T1"
                   }else{ 
                     do1.this="T2"
                   }
                   if(t[m]!=1){
                     do2.this="T3"
                   }else{
                     do2.this="T4"
                   }
                   if(i!=m){
                     do3.this="T5"
                   }else{
                     do3.this="T6"
                   }
                   switch(do1.this,
                          T1={
                            Z<-seq(from=1,to=s[i])
                            if(length(Z)!=1){
                              perm<-sample(Z,size=length(Z))
                            }else{
                              perm<-Z
                            }
                            for(j in seq(1,s[1])){ 
                              candidate<-seq(from=(perm[j]-1)*t[i]+1,
                                             to=(perm[j])*t[i])
                              bad_candidate<-rep(0,t[i])
                              v=1
                              repeat{
                                permnp[j]<-sample(setdiff(candidate,bad_candidate),1) 
                                judge<-rep(0,k)
                                p=1
                                for(w in 1:(length(s))){
                                  if(dim(illegal_set[[w]])[1]==0){
                                    
                                  }else{
                                    for(l in 1:dim(illegal_set[[w]])[1]){
                                      len<-length(setdiff(illegal_set[[w]][l,],permnp[1:j]))
                                      judge[p]<-(len>=(length(unique(illegal_set[[w]][l,]))-1))
                                      p=p+1
                                    }
                                  }
                                }
                                if(all(judge==TRUE)){
                                  break
                                }
                                bad_candidate[v]<-permnp[j]
                                v=v+1
                              }
                            }
                            C<-ceiling(permnp[1:s[1]]/t[2]) 
                          },
                          T2={
                            switch(do2.this,
                                   T3={
                                     Z<-seq(from=1,to=s[i])
                                     if(length(setdiff(Z,C))!=1){
                                       perm<-sample(setdiff(Z,C),size=length(setdiff(Z,C)))
                                     }else{
                                       perm<-setdiff(Z,C)
                                     }
                                     for(j in seq((s[i-1]+1),s[i])){
                                       candidate<-seq(from=(perm[j-s[i-1]]-1)*t[i]+1,
                                                      to=(perm[j-s[i-1]])*t[i])
                                       bad_candidate<-rep(0,t[i])
                                       v=1 
                                       repeat{
                                         permnp[j]<-sample(setdiff(candidate,bad_candidate),size=1)
                                         if(i!=m){
                                           judge<-rep(0,k)
                                           p=1
                                           for(w in 1:(length(s))){
                                             for(l in 1:(dim(illegal_set[[w]])[1])){
                                               if(dim(illegal_set[[w]])[1]==0){
                                                 
                                               }else{
                                                 len<-length(setdiff(illegal_set[[w]][l,],permnp[1:j]))
                                                 judge[p]<-(len>=(length(unique(illegal_set[[w]][l,]))-1))
                                                 p=p+1
                                               }
                                             }
                                           }
                                           if(all(judge==TRUE)){
                                             break
                                           }
                                           bad_candidate[v]<-permnp[j]
                                           v=v+1
                                         }else{
                                           
                                           break
                                         } 
                                       }
                                     }
                                     C<-ceiling(permnp[1:s[i]]/t[i+1])
                                   },
                                   T4={
                                     switch(do3.this,
                                            T5={
                                              Z<-seq(from=1,to=s[i])
                                              if(length(setdiff(Z,C))!=1){
                                                perm<-sample(setdiff(Z,C),size=length(setdiff(Z,C)))
                                              }else{
                                                perm<-setdiff(Z,C)
                                              }
                                              for(j in seq((s[i-1]+1),s[i])){
                                                candidate<-seq(from=(perm[j-s[i-1]]-1)*t[i]+1,
                                                               to=(perm[j-s[i-1]])*t[i])
                                                bad_candidate<-rep(0,t[i])
                                                v=1
                                                repeat{
                                                  permnp[j]<-sample(setdiff(candidate,bad_candidate),size=1)
                                                  judge<-rep(0,k)
                                                  p=1
                                                  for(w in 1:(length(s))){
                                                    for(l in 1:dim(illegal_set[[w]])[1]){
                                                      if(dim(illegal_set[[w]])[1]==0){
                                                        
                                                      }else{
                                                        len<-length(setdiff(illegal_set[[w]][l,],permnp[1:j]))
                                                        judge[p]<-(len>=(length(unique(illegal_set[[w]][l,]))-1))
                                                        p=p+1
                                                      }
                                                    }
                                                  }
                                                  if(all(judge==TRUE)){
                                                    break
                                                  }
                                                  bad_candidate[v]<-permnp[j]
                                                  v=v+1
                                                }
                                              }
                                              C<-ceiling(permnp[1:s[i]]/t[i+1])
                                              
                                            },
                                            T6={
                                              Z<-seq(from=1,to=s[i])
                                              if(length(setdiff(Z,C))!=1){
                                                perm<-sample(setdiff(Z,C),size=length(setdiff(Z,C)))
                                              }else{
                                                perm<-setdiff(Z,C)
                                              }
                                              for(j in seq((s[i-1]+1),s[i])){
                                                permnp[j]<-seq(from=(perm[j-s[i-1]]-1)*t[i]+1,
                                                               to=(perm[j-s[i-1]])*t[i])
                                              }
                                              C<-ceiling(permnp[1:s[i]]/t[i+1])
                                            }
                                     )
                                   }
                            )
                          }
                   )
                   return(list(permnp,C))
                 },
                 ##******************** GNLH_permutation() ***************************##
                 GNLH_permutation=function(){
                   n<-self$s[length(self$s)]
                   s<-self$s #the structure of GNLH
                   lcm<-self$Lcm #least common multiplier
                   m<-length(s) # total m layers
                   t<-self$t #the multipliers of each layer
                   t[length(s)+1]=1
                   illegal_set<-self$GNLH_illegal_set()
                   k=0
                   for(w in 1:(length(s))){
                     k<-k+dim(illegal_set[[w]])[1]##row numbers in illegal_set
                   }
                   permnp<-rep(0,len=lcm)
                   C<-list()
                   for(i in 1:m){
                     consequence<-(self$diagnose(i,illegal_set,k,permnp,C))
                     permnp<-consequence[[1]]
                     C<-consequence[[2]]
                   }
                   if((t[m]!=1)){
                     C<-permnp[1:s[m]]
                     Z<-seq(from=1,to=lcm)
                     perm<-sample(setdiff(Z,C),size=length(setdiff(Z,C)))  
                     permnp[(s[m]+1):lcm]<-perm[1:(lcm-s[m])]
                   }
                   return(permnp)
                 },
                 ##********************** GNLH_FULL *****************************##
                 GNLH_Full=function(){
                   num<-self$s[length(self$s)]
                   dim<-self$q
                   lcm<-self$Lcm
                   GNLH_Full<-matrix(rep(0,lcm*dim),ncol=dim)
                   GNLH<-matrix(rep(0,num*dim),ncol=dim)
                   for(i in 1:dim){
                     GNLH_Full[,i]<-self$GNLH_permutation()
                   }
                   return(GNLH_Full)
                 },
                 ##********************** GNLH **********************************##
                 GNLH=function(){
                   num<-self$s[length(self$s)]
                   GNLH_Full<-self$GNLH_Full()
                   GNLH<-GNLH_Full[1:num,]
                   return(GNLH)
                 }
                 ,
                 ##********************* Swap operation ************************##
                 Swap=function(GNLH_Full,structure,column,Swap_type,Swap_layer,lcm){
                   GNLH_Full_old<-GNLH_Full
                   GNLH_Full_new<-GNLH_Full_old
                   col<-column
                   type<-Swap_type
                   layer<-Swap_layer
                   s<-structure
                   m<-length(s)
                   lcm<-lcm
                   legal_set<-function(candidate,GNLH_Full_old_col,old_row,structure,layer,lcm){
                     GNLH_Full_old_col<-GNLH_Full_old_col
                     candidate<-candidate
                     len<-length(candidate)
                     legal_set<-rep(0,len)
                     s<-structure
                     m<-length(s)
                     lcm<-lcm
                     t<-lcm/s
                     layer<-layer
                     old_row<-old_row
                     for(i in 1:len){
                       candidatei<-candidate[i]
                       positioni<-which(GNLH_Full_old_col==candidatei)
                       if(positioni<=s[length(s)]){
                         for(j in 2:length(s)){
                           if((positioni>s[j-1])&&(positioni<=s[j])){
                             positionlayeri<-j
                           }  
                         }
                         judge<-rep(NA,positionlayeri-layer)
                         p=1
                         for(k in (layer):(positionlayeri-1)){
                           set_layer<-ceiling(setdiff(GNLH_Full_old_col[1:s[k]],old_row)/t[k])
                           if(length(setdiff(set_layer,ceiling(candidatei/t[k])))==length(set_layer)){
                             judge[p]<-TRUE
                           }else{
                             judge[p]<-FALSE
                           }
                           p=p+1
                         }
                         if(all(judge==TRUE)){
                           legal_set[i]=i
                         }else{
                           legal_set[i]=FALSE
                         }
                       }else{
                         judge<-rep(NA,length(s)-layer+1)
                         p=1 
                         for(k in (layer):(length(s))){
                           set_layer<-ceiling(setdiff(GNLH_Full_old_col[1:s[k]],old_row)/t[k])
                           if(length(setdiff(set_layer,ceiling(candidatei/t[k])))==length(set_layer)){
                             judge[p]<-TRUE
                           }else{
                             judge[p]<-FALSE
                           }
                           p=p+1
                         }
                         if(all(judge==TRUE)){
                           legal_set[i]=i
                         }else{
                           legal_set[i]=FALSE
                         }
                       }
                     }
                     legal_set<-candidate[legal_set]
                     return(legal_set)
                   }
                   if(layer==1){ #get the number of elements waiting for swapping
                     num<-s[layer]
                   }else{
                     num<-s[layer]-s[layer-1]
                   }
                   if(type=="in"){
                     do1.this<-"T1" # turn to layer_in swap
                   }else{
                     do1.this<-"T2" # turn to layer_between swap
                   }
                   if(layer!=m){
                     do2.this<-"T3"
                   }else{
                     do2.this<-"T4"
                   }
                   if(lcm==s[length(s)]){
                     do3.this<-"T5"
                   }else{
                     do3.this<-"T6"
                   }
                   switch(do1.this,
                          T1={
                            if(num==1){
                              GNLH_Full_new<-GNLH_Full_old
                            }else{
                              exchange<-sample(num,2)
                              i<-exchange[1]
                              j<-exchange[2]
                              if(layer==1){
                                GNLH_Full_new[i,col]=GNLH_Full_old[j,col]
                                GNLH_Full_new[j,col]=GNLH_Full_old[i,col]
                              }else{
                                start<-s[layer-1]
                                starti<-start+i
                                startj<-start+j
                                GNLH_Full_new[starti,col]<-GNLH_Full_old[startj,col]
                                GNLH_Full_new[startj,col]<-GNLH_Full_old[starti,col]
                              }
                            }     
                          },
                          T2={
                            switch(do2.this,
                                   T3={
                                     if(layer==1){
                                       i<-sample(num,1)
                                       element<-GNLH_Full_old[i,col]
                                       t<-lcm/s[layer]
                                       tag<-ceiling(element/t)
                                       domain<-seq(from=(tag-1)*t+1,to=tag*t)
                                       candidate<-setdiff(domain,element)
                                       legal_set<-legal_set(candidate,GNLH_Full_old[,col],element,s,layer,lcm)
                                       if(length(legal_set)==0){
                                         GNLH_Full_new<-GNLH_Full_old
                                       }else{
                                         if(length(legal_set)==1){
                                           exchange<-legal_set
                                         }else{
                                           exchange<-sample(legal_set,1)
                                         }
                                         j<-which(GNLH_Full_old[,col]==exchange)
                                         GNLH_Full_new[i,col]=GNLH_Full_old[j,col]
                                         GNLH_Full_new[j,col]=GNLH_Full_old[i,col]
                                       }
                                     }else{
                                       i<-sample(num,1)+s[layer-1]
                                       element<-GNLH_Full_old[i,col]
                                       t<-lcm/s[layer]
                                       tag<-ceiling(element/t)
                                       domain<-seq(from=(tag-1)*t+1,to=tag*t)
                                       candidate<-setdiff(domain,element)
                                       legal_set<-legal_set(candidate,GNLH_Full_old[,col],element,s,layer,lcm)
                                       if(length(legal_set)==0){
                                         GNLH_Full_new<-GNLH_Full_old
                                       }else{
                                         if(length(legal_set)==1){
                                           exchange<-legal_set
                                         }else{
                                           exchange<-sample(legal_set,1)
                                         }
                                         j<-which(GNLH_Full_old[,col]==exchange)
                                         GNLH_Full_new[i,col]=GNLH_Full_old[j,col]
                                         GNLH_Full_new[j,col]=GNLH_Full_old[i,col]
                                       }                                               
                                     }
                                   },
                                   T4={
                                     switch(do3.this,
                                            T5={
                                              if(num==1){
                                                GNLH_Full_new<-GNLH_Full_old
                                              }else{
                                                exchange<-sample(num,2)
                                                i<-exchange[1]
                                                j<-exchange[2]
                                                start<-s[layer-1]
                                                starti<-start+i
                                                startj<-start+j
                                                GNLH_Full_new[starti,col]<-GNLH_Full_old[startj,col]
                                                GNLH_Full_new[startj,col]<-GNLH_Full_old[starti,col]
                                              }
                                            },
                                            T6={
                                              i<-sample(num,1)+s[layer-1]
                                              element<-GNLH_Full_old[i,col]
                                              t<-lcm/s[layer]
                                              tag<-ceiling(element/t)
                                              domain<-seq(from=(tag-1)*t+1,to=tag*t)
                                              candidate<-setdiff(domain,element)
                                              legal_set<-legal_set(candidate,GNLH_Full_old[,col],element,s,layer,lcm)
                                              if(length(legal_set)==0){
                                                GNLH_Full_new<-GNLH_Full_old
                                              }else{
                                                if(length(legal_set)==1){
                                                  exchange<-legal_set
                                                }else{
                                                  exchange<-sample(legal_set,1)
                                                }
                                                j<-which(GNLH_Full_old[,col]==exchange)
                                                GNLH_Full_new[i,col]=GNLH_Full_old[j,col]
                                                GNLH_Full_new[j,col]=GNLH_Full_old[i,col]
                                              }
                                            }
                                     )
                                   }
                            )
                          }
                   )
                   return(GNLH_Full_new)
                 },
                 ##******************************* StandGNLHD ***********************************##
                 StandGNLHD=function(){
                   lcm<-self$Lcm#number of samples
                   StandGNLHD<-(self$GNLH()-0.5)/lcm
                   return(StandGNLHD)
                 }
                 ##******************************************************************************##
                 
               )
)
#'@export

