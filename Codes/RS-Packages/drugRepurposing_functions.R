
#### Function to determine the optimal rank using SVD
find_optimal_rank = function(A, threshold) {
  svd_result = svd(A)
  singular_values = svd_result$d
  
  # Calculate the cumulative explained variance
  cumulative_variance = cumsum(singular_values^2) / sum(singular_values^2)
  
  # Find the rank where the cumulative variance exceeds a certain threshold (e.g., 0.95 => 95%)
  optimal_rank = which(cumulative_variance >= threshold)[1]
  return(optimal_rank)
}


#### Non-negative Matrix Factorization (NMF) model using Lee-Seung update rules
NMF_ELee<-function(V, k, M=NULL, itnum = 100, eps=1e-9, epsDiff=1e-4, seed=myseed){
  # V: input matrix
  # k: selected rank
  # M: masked matrix (use in future)
  # itNum: number of iterations
  
  n=nrow(V)
  m=ncol(V)
  W<-matrix(runif(n*k), nrow = n)
  H<-matrix(runif(k*m), nrow = k)
  
  # W<-matrix(rnorm(n*k), nrow = n)
  # H<-matrix(rnorm(k*m), nrow = k)
  
  V1=V
  errDiff=err= NULL
  e1=1000
  for(i in 1:itnum){
    u_H<-t(W)%*%V1
    d_H<-t(W)%*%W%*%H + 1e-9
    H<-H*(u_H/d_H)
    
    u_W<-V1%*%t(H)
    d_W<-W%*%H%*%t(H) + 1e-9
    W<-W*(u_W/d_W)
    
    e= sum((V1-W%*%H)^2)
    ed= abs(e-e1)/e1
    
    err=c(err,e)
    errDiff=c(errDiff,ed)
    
    if (e < eps |  ed < epsDiff) break()
    
    e1=e
    
    if (!is.null(M)) V1= M*V + (1-M) * V1 # Trang's proposal
    #if (!is.null(M)) V1= (1-M)*V + M * V1
  }
  return(list(W = W, H = H, err = err, errDiff=errDiff))
}



#### Sort matrix
sort_matrix<-function(score_matrix,interact_matrix){
  #For each disease in columns, sort the score values (predicted) and change the values of the interaction matrix accordingly. In downstream step, go down to each row, compute performance metrics for top 1, 2, 3, .. to summarize the final results.
  
  sort_index=score_matrix
  for ( i in 1:ncol(sort_index)){
    x=-score_matrix[,i]
    #x=score_matrix[,i]
    o=order(x, decreasing = FALSE)
    sort_index[,i]=o
  }
  
  score_sorted = score_matrix
  y_sorted = interact_matrix
  for (i in 1:ncol(interact_matrix)){
    score_sorted[,i] = score_matrix[,i][sort_index[,i]]
    y_sorted[,i] = interact_matrix[,i][sort_index[,i]]
  }
  rownames(score_sorted)=NULL #row names is not meaningful, so remove them to avoid potential misunderstanding
  rownames(y_sorted)=NULL
  
  res=list(y_sorted=y_sorted,score_sorted=score_sorted,sort_index=sort_index)
  return(res)
  
}


##### Fold - cross validation

#n: number of element
#k:  number of folds
makeFolds<-function(k, n, seedVal){
  # to do LOO-CV, just set k equal to n
  set.seed(seedVal)
  s=trunc(n/k)
  foldId=rep(c(1:k),s+1)[1:n]
  sampleIDx=sample(seq(n))
  kFolds=split(sampleIDx,foldId)
  return(kFolds)
}


convert_element = function(nr, nc, k){
  # nr: number of row of original matrix
  # nc: number of column of original matrix
  # k: index of element
  # convert to the position on matrix for one element
  cl = trunc((k-1)/nr) + 1 ## column of element
  ro = k-(cl-1)*nr  ## row of element 
  return (c(ro,cl))
}
