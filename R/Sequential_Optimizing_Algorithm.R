##***********************************************************************************##
##*************************** An efficient sequential algorithm *********************##
#'Is_repeat function
#'@include Sequential_Optimizing_Algorithm.R
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
Optimal_GNLHD_SequentialAlg=function(GNLHD,GNLH_Full,iteration){
  GNLH_Full<-GNLH_Full
  s<-GNLHD$s
  m<-length(s) #the number of layers
  t<-GNLHD$t
  lcm<-GNLHD$Lcm
  q<-GNLHD$q
  for(k in seq(from=2,to=m)){
    if(Is_repeat(ceiling(GNLH_Full[1:s[k],]/t[k-1]))==0){
      for(l in 1:iteration){
        col_swap<-sample(q,1)
        GNLH_try_Full<-GNLHD$Swap(GNLH_Full,s,col_swap,"between",k,lcm)
        if(Phi_p(GNLH_try_Full[1:s[m],],2,50)<Phi_p(GNLH_Full[1:s[m],],2,50)){
          GNLH_Full<-GNLH_try_Full}
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
        if(Is_repeat(ceiling(GNLH_try_Full[1:s[k],]/t[k-1]))==0){
          GNLH_Full<-GNLH_try_Full
          break
        }
      }
    }
    for(l in 1:iteration){
      col_swap<-sample(q,1)
      GNLH_try_Full<-GNLHD$Swap(GNLH_Full,s,col_swap,"between",k,lcm)
      if(Phi_p(GNLH_try_Full[1:s[m],],2,50)<Phi_p(GNLH_Full[1:s[m],],2,50)){
        GNLH_Full<-GNLH_try_Full}
    }
  }
  return(GNLH_Full)
}
##***************************************************************************##
#'@export