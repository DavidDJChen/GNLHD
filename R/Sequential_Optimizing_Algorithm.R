##***********************************************************************************##
##*************************** An efficient sequential algorithm *********************##
#'Is_repeat function
#'@include Sequential_Optimizing_Algorithm.R
##***********************************************************************************##
###*********** the same layer operation ****************###
samelayer_exchange <- function(column,J1){
  J1 <- J1
  column <- column
  l <- length(column)
  if(J1 <= (l*(l-1)/2)){
    J1 <- J1
  }else J1 <- l*(l-1)/2 #Judging whether J1 is larger than 
  #possible exchanging number
  exchange_output <- matrix(rep(column,J1),nrow=l)
  count <- 1
  if(count == 1){
    i <- sample(seq(from=1,to=l-1),1)
    if(i == l-1){
      j = l
    }else{
      j <- sample(seq(from = i+1,to = l), 1)
    }
    column_exchange <- column
    column_exchange[i] <- column[j]
    column_exchange[j] <- column[i]#exchanging the two elements 
    exchange_output[,count] <- column_exchange
    count <- count+1
  }
  while(count <= J1){
    i <- sample(seq(from = 1,to = l-1),1)
    if(i == l-1){
      j = l
    }else{
      j <- sample(seq(from = i+1,to = l),1)
    }
    column_exchange <- column
    column_exchange[i] <- column[j]
    column_exchange[j] <- column[i]#exchanging the two elements 
    discard <- rep(0,count-1)
    for(k in seq(1,count-1)){
      if(all(exchange_output[,k] == column_exchange)){
        discard[k] <- 1
      }else discard[k] <- 2
    }#Checking whether any exchanged column is equal to lastest column  
    
    if(all(discard == 2)){ 
      exchange_output[,count] <- column_exchange
      count <- count+1
    }
  }
  return(exchange_output) 
}
###*****************************************************###

###**************** Inner Loop of ESE for GNLHD******************###
innerloop <- function(GNLHD_best,GNLHD,s,M,J,t,p,T_h,lcm){#The inner loop function 
  n_acpt <- 0
  n_imp <- 0
  M <- M #The number of upmost column operation
  J <- J #The number of total new GNLHDs
  s <- s #The structure of GNLHD
  T_h <- T_h
  GNLHD <- GNLHD
  column_num <- dim(GNLHD)[2]
  number_in <- J #The number of layer-in swap operations 
  
  
  for(i in seq(1,M-1)){
    column_index <- (i%%column_num)+1 #The index of operating column
    column_operate <- GNLHD[,column_index] #The column operating
    column_exchange_in1 <- samelayer_exchange(column_operate[seq(1,s[1])],number_in) #The layer-in operation
    column_exchange_in <- matrix(rep(column_operate,number_in),
                                 nrow=length(column_operate))
    column_exchange_in[seq(1,s[1]),seq(1,number_in)] <- column_exchange_in1 #The new GNLHDs by layer-in operation
    
    
    column_exchange <- column_exchange_in #The new GNLHDs by two kinds of operations 
    t = t #variable parameter for Phi_p value
    p = p #variable parameter for Phi_p value
    Phi_p_value<-rep(0,dim(column_exchange)[2])
    GNLHD_exchange<-GNLHD
    for(j in seq(1,dim(column_exchange)[2])){ #get all design Phi_p_value
      GNLHD_exchange_prep <- GNLHD_exchange
      GNLHD_exchange_prep[,column_index] <- column_exchange[,j]
      Phi_p_value[j] <- Phi_p((GNLHD_exchange_prep[1:s[1],]-0.5)/lcm,t,p)
    } 
    try_value <- min(Phi_p_value)#It is possible to get more than one samllest try_index
    try_index <- min(which(Phi_p_value==try_value))#we get the smallest coordiate
    GNLHD_try <- GNLHD
    GNLHD_try[,column_index] <- column_exchange[,try_index]
    if((Phi_p_value[try_index]-Phi_p((GNLHD[1:s[1],]-0.5)/lcm,t,p)) <= T_h*runif(min=0,max=1,1)){
      GNLHD <- GNLHD_try
      n_acpt <- n_acpt+1
      if(Phi_p((GNLHD[1:s[1],]-0.5)/lcm,t,p) < Phi_p((GNLHD_best[1:s[1],]-0.5)/lcm,t,p)){
        GNLHD_best <- GNLHD
        n_imp <- n_imp+1
      }
    }  
  }
  return(GNLHD_best)
}
###********************************************************************###


###***************** outerloop of ESE *********************************###
outerloop<-function(GNLHD_initial,s,T_h_initial=0.1,M=100,J=6,t=2,p=50,
                    tolerance=0.1,alpha=c(0.8,0.9,0.7),lcm){
  M <- M
  J <- J 
  t <- t # The Phi_p variable
  p <- p # The Phi_p variable
  tolerance <- tolerance
  alpha <- alpha
  lcm <- lcm
  
  GNLHD <- GNLHD_initial
  GNLHD_best <- GNLHD
  T_h <- T_h_initial # Initial value is a small value
  for(l in 1:10){ # The convergence criteria
    GNLHD_oldbest <- GNLHD_best
    i <- 0
    n_acpt <- 0
    n_imp <- 0
    GNLHD_best <- innerloop(GNLHD_best,GNLHD,s,M,J,t,p,T_h,lcm)
    if((Phi_p((GNLHD_oldbest[1:s[1],]-0.5)/lcm,t,p)-Phi_p((GNLHD_best[1:s[1],]-0.5)/lcm,t,p)) > tolerance){
      flag_imp <- 1
    }else{
      flag_imp <- 0
    }
    if(flag_imp == 1){ # The improvement process
      if((n_acpt/M >= 0.1)&((n_imp/M) < (n_acpt/M))){
        T_h <- alpha[1]*T_h
      }else{
        if((n_acpt/M>=0.1)&((n_imp/M) == (n_acpt/M))){
          T_h <- T_h
        }else{
          T_h <- T_h/alpha[1]
        }
      }
    }else{
      if((n_acpt/M) < 0.1){
        T_h <- T_h/alpha[3]
      }else{
        T_h <- T_h*alpha[2]
      }
    }
  }
  return(GNLHD_best)
}

###****************************************************************###  


##***********************************************************************************##
##*************************** An efficient sequential algorithm *********************##

Is_repeat=function(Design){
  Is_repeat<-0
  row<-dim(Design)[1]
  col<-dim(Design)[2]
  dis_mat<-distance(Design,2)#L2 di#the distance matrix
  for(i in 2:row){
    for(j in 1:(i-1)){
      if(dis_mat[i,j]==0){
        Is_repeat<-Is_repeat+1
      }
    }
  }
  return(Is_repeat)
}
##**************************************************************##
#'Optimal_GNLHD_SequentialAlg function
##**************************************************************##
Optimal_GNLHD_SequentialAlg=function(GNLHD_in,GNLH_Full,iteration,w,T_h_initial=0.1,M=100,J=6,t=2,p=50,
                                     tolerance=0.1,alpha=c(0.8,0.9,0.7)){
  GNLH_Full<-GNLH_Full
  s<-GNLHD_in$s
  m<-length(s) #the number of layers
  t<-GNLHD_in$t
  lcm<-GNLHD_in$Lcm
  q<-GNLHD_in$q
  w<-w
  M<-M
  J<-J
  t<-t
  p<-p
  tolerance<-tolerance
  alpha<-alpha
  
  GNLH_Full<-outerloop(GNLHD_initial=GNLH_Full,s,T_h_initial=0.01,M=30,J=4,t=2,p=50,
                       tolerance=0.1,alpha=c(0.8,0.9,0.7),lcm=720) # Optimze the first layer by using ESE algorithm
  for(k in seq(from=2,to=m)){
    if(Is_repeat(ceiling(GNLH_Full[1:s[k],]/t[k-1]))==0){
      for(l in 1:iteration){
        col_swap<-sample(q,1)
        GNLH_try_Full<-GNLHD_in$Swap(GNLH_Full,s,col_swap,"between",k,lcm)
        f_value_try<-0
        f_value<-0
        for(v in 1:length(s)){
          f_value_try<-f_value_try+w[v]*Phi_p((GNLH_try_Full[1:s[v],]-0.5)/lcm,2,50)
        }
        for(v in 1:length(s)){
          f_value<-f_value+w[v]*Phi_p((GNLH_Full[1:s[v],]-0.5)/lcm,2,50)
        }
        if(f_value_try<f_value){
          GNLH_Full<-GNLH_try_Full
        }
      }
    }else{
      for(e in 1:iteration){
        dis_mat<-distance(ceiling(GNLH_Full[1:s[k],]/t[k-1]),2)
        dis_mat_1<-dis_mat+diag(s[k])
        k_row<-(which(dis_mat_1==0))%%s[k]
        k_row[which(k_row==0)]<-s[k]
        k_col<-ceiling(which(dis_mat_1==0)/s[k])
        repeat_row<-as.vector(k_row[which(k_row>s[k-1])]) # the repeated rows 
        col_swap<-sample(q,1) # the swapping column 
        if(length(repeat_row)==1){
          i_row<-repeat_row
        }else{
          i_row<-sample(repeat_row,1)
        }
        if(length(setdiff(seq(s[k-1]+1,s[k]),i_row))==1){
          j_row<-setdiff(seq(s[k-1]+1,s[k]),i_row)
        }else{
          j_row<-sample(setdiff(seq(s[k-1]+1,s[k]),i_row),1)
        }
        GNLH_try_Full<-GNLH_Full
        GNLH_try_Full[i_row,col_swap]<-GNLH_Full[j_row,col_swap]
        GNLH_try_Full[j_row,col_swap]<-GNLH_Full[i_row,col_swap]
        if(Is_repeat(ceiling(GNLH_try_Full[1:s[k],]/t[k-1]))<Is_repeat(ceiling(GNLH_Full[1:s[k],]/t[k-1]))){
          GNLH_Full<-GNLH_try_Full
        }
        if(Is_repeat(ceiling(GNLH_try_Full[1:s[k],]/t[k-1]))==0){
          GNLH_Full<-GNLH_try_Full
          break
        }
      }
    }
    for(l in 1:iteration){
      col_swap<-sample(q,1)
      GNLH_try_Full<-GNLHD_in$Swap(GNLH_Full,s,col_swap,"between",k,lcm)
      f_value_try<-0
      f_value<-0
      for(v in 1:length(s)){
        f_value_try<-f_value_try+w[v]*Phi_p((GNLH_try_Full[1:s[v],]-0.5)/lcm,2,50)
      }
      for(v in 1:length(s)){
        f_value<-f_value+w[v]*Phi_p((GNLH_Full[1:s[v],]-0.5)/lcm,2,50)
      }
      if(Phi_p((GNLH_try_Full[1:s[k],]-0.5)/lcm,2,50)<Phi_p((GNLH_Full[1:s[k],]-0.5)/lcm,2,50)){
        GNLH_Full<-GNLH_try_Full
      }
    }
  }
  return(GNLH_Full)
}
##***************************************************************************##
##***************************************************************************##
#'@export
