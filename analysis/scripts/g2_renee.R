
#
# Function: get.g2stat
# Computes the G2 test statistic given two data matrices
# This is used internally by the G2 function
# Inputs: Z, W: both square, symmetric matrices with an equal number of rows
# Output: test statistic (single value)
# This is the most efficient way I found of computing this, making use of a trace property
get.g2stat.vec <- function(W,Z)
{
 g2tstat <- sum( crossprod( matrix(W,nrow=nrow(W)*ncol(W),ncol=1), matrix(Z,nrow=nrow(Z)*ncol(Z),ncol=1)  ))
 g2tstat
}


#
# Function: get.g2stat.multin
# Computes the G2 test statistic given two data matrices, under a multinomial distribution
# Here we write the test statistic as a function of W, with no particular structure given
# We also compute it in a more efficient way than we used to - using a trace property
# This is used internally by the G2 function
# Inputs: 
#  U = U1 U1^t, where U1 = (I-H)Y, a n*K matrix where n=number obs and K=number multinomial responses possible
#  tau.mat = X' W X, a n*n matrix 
# Output: test statistic (single value)
# 
get.g2stat.multin <- function(U, tau.mat)
{
  g2tstat <- sum( U * tau.mat )
  g2tstat
}



#
# Function: get.g2stat
# Computes the G2 test statistic given two data matrices
# This is used internally by the G2 function
# Inputs: Z, W: both square, symmetric matrices with an equal number of rows
# Output: test statistic (single value)
# This is the most efficient way I found of computing this, making use of a trace property
get.g2stat.had <- function(W,Z)
{
 g2tstat <- sum( W * Z )
 g2tstat
}

#
# Function: G2
### Input 
### dep data and indep data with samples on the rows and genes on the columns
### grouping: Either a logical value = F or A matrix with a single column and same number of rows as samples. 
###         Column name should be defined.
###         Contains clinical information of the samples. 
###         Should have two groups only. 
### stand : scaling the columns of the data, logical value
### nperm : number of permutations 

### Output
### A list containing G2 p.values and G2 test statistics

### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$$Sg


G2 <- function(dep.data,indep.data,stand = TRUE, nperm=100,grouping=F){
      
    nperm = nperm
      ## check for the number of samples in dep and indep data

      if (nrow(dep.data)!=nrow(indep.data)){
         cat("number of samples not same in dep and indep data","\n")
      }

     

      ### check for centering and standardizing the data
       
      dep.data = scale(dep.data,center=T,scale=stand)
      indep.data = scale(indep.data,center=T,scale=stand)

        

     
#################################################################################

   ### If grouping is T divide the dataset into two groups and get G2 pval for each group.
   
   if(length(grouping) > 1 ){

    ## check the number of groups
  
      if (length(table(grouping))>2){
       cat("number of groups cannot exceed 2, please check grouping variable","\n")
      }
   ## Group 1
   dep = dep.data[grouping == names(table(grouping))[1],]
   indep = indep.data[grouping == names(table(grouping))[1],]
   
   ## Calculate Z = XX' and W = YY'
   Z = tcrossprod(indep)
   W = tcrossprod(dep)
   ### G2 for Group1
   samp_names = rownames(indep)
   Sg1 = get.g2stat.had(W,Z)
  

   ### Permutations

   perm_samp = matrix(0,nrow=nrow(indep),ncol=nperm)   ## generate the permutation matrix
       for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
        }

   ## permutation starts
   for (perm in 1:nperm){
        permX = Z[perm_samp[,perm],]
        permX = permX[,perm_samp[,perm]]

        Sg1 = c( Sg1, get.g2stat.had(W,permX) )
      }
 
      perm_samp.g1 = perm_samp  ## permutation matrix for output
##############################################################################

   ## Group 2
   dep = dep.data[grouping == names(table(grouping))[2],]
   indep = indep.data[grouping == names(table(grouping))[2],]

   ## Calculate Z = XX' and W = YY'
   Z = tcrossprod(indep)
   W = tcrossprod(dep)
   
   ### G2 for Group2
   samp_names = rownames(indep)
   Sg2 = get.g2stat.had(W,Z)
  

   ### Permutations

   perm_samp = matrix(0,nrow(indep),nperm)   ## generate the permutation matrix
       for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
        }

   ## permutation starts
   for (perm in 1:nperm){
        permX = Z[perm_samp[,perm],]
        permX = permX[,perm_samp[,perm]]

        Sg2 = c( Sg2, get.g2stat.had(W,permX) )
      }

       perm_samp.g2 = perm_samp   #### permutation matrix for output

########################################################################

    #### G2 test statistic
    Sg1 = t(as.matrix(Sg1))
    Sg2 = t(as.matrix(Sg2))
    colnames(Sg1) = rep("g1",length(Sg1))
    colnames(Sg2) = rep("g2",length(Sg2))

    ### calculate the p-values for g1 and g2

    G2p1 = mean(Sg1[1]<= c(Inf , Sg1[2:(nperm+1)]))
    G2p2 = mean(Sg2[1]<= c(Inf , Sg2[2:(nperm+1)]))
    
    ### prepare for output
    G2p = c(G2p1,G2p2)
    Sg = cbind(Sg1,Sg2)
    perm_samp = list(perm_samp.g1, perm_samp.g2)
    names(perm_samp) = c("g1.perm.mat","g2.perm.mat")
    names(G2p) = c(paste(colnames(grouping)," status",":",names(table(grouping))[1],sep=""),paste(colnames(grouping)," status",":",names(table(grouping))[2],sep=""))
  
   } else {
         #### No  grouping of the samples.
         
         samp_names = 1:nrow(indep.data)
         
         ## Calculate Z = XX' and W = YY'
         Z = tcrossprod(indep.data)
         W = tcrossprod(dep.data)
   
         ### G2 
         Sg = get.g2stat.had(W,Z)
  

         ### Permutations

         perm_samp = matrix(0,nrow(indep.data),nperm)   ## generate the permutation matrix
         for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
         }

         ## permutation starts
         for (perm in 1:nperm){
            permX = Z[perm_samp[,perm],]
            permX = permX[,perm_samp[,perm]]

            Sg = c(Sg, get.g2stat.had(W,permX) )
         }


########################################################################

         #### G2 test statistic
         Sg = t(as.matrix(Sg))

         ### Calculte G2 pval
         G2p = mean(Sg[1]<= c(Inf , Sg[2:(nperm+1)]))
         names(G2p) = "G2 pval"
      }
     
return (list(perm = perm_samp,G2p = G2p,Sg = Sg))
}

#
# Function: G2.vec
### Input 
### dep data and indep data with samples on the rows and genes on the columns
### grouping: Either a logical value = F or A matrix with a single column and same number of rows as samples. 
###         Column name should be defined.
###         Contains clinical information of the samples. 
###         Should have two groups only. 
### stand : scaling the columns of the data, logical value
### nperm : number of permutations 

### Output
### A list containing G2 p.values and G2 test statistics

### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$$Sg

# This version of G2 uses the vec of matrices, then takes the sum of their inner product
# This is faster than using sum(diag( XXT %*% YYT )), 
# but it is on average slower than sum( xxt * yyt ), the sum of the Hadamard product of the two matrices
# The latter is now implemented in the G2 function, above


G2.vec <- function(dep.data,indep.data,stand = TRUE, nperm=100,grouping=F){
      
    nperm = nperm
      ## check for the number of samples in dep and indep data

      if (nrow(dep.data)!=nrow(indep.data)){
         cat("number of samples not same in dep and indep data","\n")
      }

     

      ### check for centering and standardizing the data
       
      dep.data = scale(dep.data,center=T,scale=stand)
      indep.data = scale(indep.data,center=T,scale=stand)

        

     
#################################################################################

   ### If grouping is T divide the dataset into two groups and get G2 pval for each group.
   
   if(length(grouping) > 1 ){

    ## check the number of groups
  
      if (length(table(grouping))>2){
       cat("number of groups cannot exceed 2, please check grouping variable","\n")
      }
   ## Group 1
   dep = dep.data[grouping == names(table(grouping))[1],]
   indep = indep.data[grouping == names(table(grouping))[1],]
   
   ## Calculate Z = XX' and W = YY'
   Z = tcrossprod(indep)
   W = tcrossprod(dep)
   ### G2 for Group1
   samp_names = rownames(indep)
   Sg1 = get.g2stat.vec(W,Z)
  

   ### Permutations

   perm_samp = matrix(0,nrow=nrow(indep),ncol=nperm)   ## generate the permutation matrix
       for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
        }

   ## permutation starts
   for (perm in 1:nperm){
        permX = Z[perm_samp[,perm],]
        permX = permX[,perm_samp[,perm]]

        Sg1 = c( Sg1, get.g2stat.vec(W,permX) )
      }
 
      perm_samp.g1 = perm_samp  ## permutation matrix for output
##############################################################################

   ## Group 2
   dep = dep.data[grouping == names(table(grouping))[2],]
   indep = indep.data[grouping == names(table(grouping))[2],]

   ## Calculate Z = XX' and W = YY'
   Z = tcrossprod(indep)
   W = tcrossprod(dep)
   
   ### G2 for Group2
   samp_names = rownames(indep)
   Sg2 = get.g2stat.vec(W,Z)
  

   ### Permutations

   perm_samp = matrix(0,nrow(indep),nperm)   ## generate the permutation matrix
       for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
        }

   ## permutation starts
   for (perm in 1:nperm){
        permX = Z[perm_samp[,perm],]
        permX = permX[,perm_samp[,perm]]

        Sg2 = c( Sg2, get.g2stat.vec(W,permX) )
      }

       perm_samp.g2 = perm_samp   #### permutation matrix for output

########################################################################

    #### G2 test statistic
    Sg1 = t(as.matrix(Sg1))
    Sg2 = t(as.matrix(Sg2))
    colnames(Sg1) = rep("g1",length(Sg1))
    colnames(Sg2) = rep("g2",length(Sg2))

    ### calculate the p-values for g1 and g2

    G2p1 = mean(Sg1[1]<= c(Inf , Sg1[2:(nperm+1)]))
    G2p2 = mean(Sg2[1]<= c(Inf , Sg2[2:(nperm+1)]))
    
    ### prepare for output
    G2p = c(G2p1,G2p2)
    Sg = cbind(Sg1,Sg2)
    perm_samp = list(perm_samp.g1, perm_samp.g2)
    names(perm_samp) = c("g1.perm.mat","g2.perm.mat")
    names(G2p) = c(paste(colnames(grouping)," status",":",names(table(grouping))[1],sep=""),paste(colnames(grouping)," status",":",names(table(grouping))[2],sep=""))
  
   } else {
         #### No  grouping of the samples.
         
         samp_names = rownames(indep.data)
         
         ## Calculate Z = XX' and W = YY'
         Z = tcrossprod(indep.data)
         W = tcrossprod(dep.data)
   
         ### G2 
         samp_names = rownames(indep.data)
         Sg = get.g2stat.vec(W,Z)
  

         ### Permutations

         perm_samp = matrix(0,nrow(indep.data),nperm)   ## generate the permutation matrix
         for(i in 1:ncol(perm_samp)){
           perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
         }

         ## permutation starts
         for (perm in 1:nperm){
            permX = Z[perm_samp[,perm],]
            permX = permX[,perm_samp[,perm]]

            Sg = c(Sg, get.g2stat.vec(W,permX) )
         }


########################################################################

         #### G2 test statistic
         Sg = t(as.matrix(Sg))

         ### Calculte G2 pval
         G2p = mean(Sg[1]<= c(Inf , Sg[2:(nperm+1)]))
         names(G2p) = "G2 pval"
      }
     
return (list(perm = perm_samp,G2p = G2p,Sg = Sg))
}
 
#
# Function: get.pval.percol
# This function takes a vector containing the observed test stat as the first entry, followed by values generated by permutation,
# and computed the estimated p-value
# Input
# x: a vector with length nperm+1
# Output
# the pvalue computed
get.pval.percol <- function(x){
  pval = mean(x[1]<= c(Inf , x[2:length(x)]))
  pval
}



#
# Function: G2.multin
# This is to compute the G2 test statistic under the assumption that the response follows a multinomial distribution
### Input 
### dep.data and indep.data with samples on the rows and genes on the columns
### nperm : number of permutations 
### W : a matrix that represents the (expected) correlation sttructure between effects of different SNPs on the (same or different) exons
###    If given, it must be a square, symmetric matrix with as many rows as the number of genes/covariates/columns in indep.data
###
### Output
### A list containing G2 p.values and G2 test statistics

### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$$Sg


G2.multin <- function(dep.data,indep.data,stand = TRUE, nperm=100,W=NULL){
  
  ## check for the number of samples in dep and indep data
    
  if (nrow(dep.data)!=nrow(indep.data)){
    stop("number of samples not the same in dep and indep data","\n")
  }
  
  ## check if W has compatible dimensions
  if( !is.null(W) ){ if(nrow(W) != ncol(W)) { stop("W must be a square, symmetric matrix","\n") }
                     if( nrow(W) != ncol(indep.data) ) { stop("W must have as many rows as the number of columns in indep.data","\n") }
  } else { W <- diag(rep(1,ncol(indep.data))) }

  nresponses <- ncol(dep.data)
  ncovariates <- ncol(indep.data)
  ### centering and standardizing the data are not done in this case
  
#  dep.data = scale(dep.data,center=T,scale=stand)
#  indep.data = scale(indep.data,center=T,scale=stand)
  
    #### No  grouping of the samples.
    
    ## Calculate U=(I-H)Y and UU', where Y has observations on rows;  
    ##  also tau.mat=X*W*X', where X has observations on rows and variables on columns
    ## NOTE: this formulation uses X with n obs on the rows and m covariates on the columns, so it is the transpose of the first calculations
    nsamples <- nrow(dep.data)
    n.persample <- rowSums(dep.data)
    n.all <- sum(dep.data)
    H <- (1/n.all)*matrix( rep(n.persample,each=nsamples),nrow=nsamples,byrow=T)
    U1 <- (diag(rep(1,nsamples)) - H) %*% dep.data
    U <- tcrossprod(U1)
    tau.mat <- indep.data %*% W %*% t(indep.data)
    samp_names = 1:nsamples ## this was rownames(indep.data), but I now do this so that rownames do not have to be added to the array tau.mat
    Sg = get.g2stat.multin(U,tau.mat)
    #
    ### G2 
    ### Permutations
    # When using permutations: only the rows of tau.mat are permuted
    # To check how the permutations can be efficiently applied, see tests_permutation_g2_multin.R


    perm_samp = matrix(0, nrow=nsamples, ncol=nperm)   ## generate the permutation matrix
    for(i in 1:ncol(perm_samp)){
      perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
    }
    
    ## permutation starts - recompute tau.mat  (or recompute U each time)
    for (perm in 1:nperm){
      tau.mat.perm = tau.mat[perm_samp[,perm],,drop=FALSE]          # permute rows
      tau.mat.perm = tau.mat.perm[,perm_samp[,perm],drop=FALSE]     # permute columns
      
      Sg = c(Sg, get.g2stat.multin(U,tau.mat.perm) )
    }
    
    
    ########################################################################
    
    #### G2 test statistic
    # *** recompute for a vector of values for each case - just reformat the result with as many rows as permutations + 1,
    # and as many columns as combinations of values of rho and mu
    Sg = t(as.matrix(Sg)) 

    ### Calculte G2 pval
    G2p = mean(Sg[1]<= c(Inf , Sg[2:(nperm+1)]))
    names(G2p) = "G2 pval"
      
  return (list(perm = perm_samp,G2p = G2p,Sg = Sg))
}


#
# Function: G2.multin.factor
# This is to compute the G2 test statistic under the assumption that the response follows a multinomial distribution
# Here we test each separate column in the dep.data as though it represents successes, against the sum of all other rows, thus representing failures
# So this can be used to apply the G2 test per exon for each gene, running automatically for all exons avaialable for that gene
### Input 
### dep.data and indep.data with samples on the rows and genes on the columns
### nperm : number of permutations 
### W : a matrix that represents the (expected) correlation sttructure between effects of different SNPs on the (same or different) exons
###    If given, it must be a square, symmetric matrix with as many rows as the number of genes/covariates/columns in indep.data
###
### Output
### A list containing G2 p.values and G2 test statistics

### Example : G2T = G2(dep.data = cgh, indep.data = expr, grouping=F, stand=TRUE, nperm=1000)
### G2 p.values : G2T$G2p
### G2 TS : G2T$$Sg


G2.multin.factor <- function(dep.data,indep.data,stand = TRUE, nperm=100,W=NULL){
  
  ## check for the number of samples in dep and indep data
  
  if (nrow(dep.data)!=nrow(indep.data)){
    stop("number of samples not the same in dep and indep data","\n")
  }
  
  ## check if W has compatible dimensions
  if( !is.null(W) ){ if(nrow(W) != ncol(W)) { stop("W must be a square, symmetric matrix","\n") }
                     if( nrow(W) != ncol(indep.data) ) { stop("W must have as many rows as the number of columns in indep.data","\n") }
  } else { W <- diag(rep(1,ncol(indep.data))) }
  
  nsamples <- nrow(dep.data)
  nresponses <- 2
  ncovariates <- ncol(indep.data)
  ### centering and standardizing the data are not done in this case
  
  ## Calculate U=(I-H)Y and UU', where Y has observations on rows;  
  ##  also tau.mat=X*W*X', where X has observations on rows and variables on columns
  ## NOTE: this formulation uses X with n obs on the rows and m covariates on the columns, so it is the transpose of the first calculations
  ## Here we compute one test stat per dependent variable, considering also the sum of all other dependent variables
  all.sg <- NULL
  n.persample <- rowSums(dep.data)
  n.all <- sum(dep.data)
  H <- (1/n.all)*matrix( rep(n.persample,each=nsamples),nrow=nsamples,byrow=T)
  tau.mat <- indep.data %*% W %*% t(indep.data)
  samp_names = 1:nsamples ## this was rownames(indep.data), but I now do this so that rownames do not have to be added to the array tau.mat
  
  for(xj in 1:ncol(dep.data))
  {   
    dep.bin <- cbind( dep.data[,xj],rowSums(dep.data[,-xj]) )
    U1 <- (diag(rep(1,nsamples)) - H) %*% dep.bin
    U <- tcrossprod(U1)
    Sg = get.g2stat.multin(U,tau.mat)
    all.sg <- c(all.sg,Sg)
  }
  #
  ### G2 
  ### Permutations
  # When using permutations: only the rows of tau.mat are permuted
  # To check how the permutations can be efficiently applied, see tests_permutation_g2_multin.R
  
  
  perm_samp = matrix(0, nrow=nsamples, ncol=nperm)   ## generate the permutation matrix
  for(i in 1:ncol(perm_samp)){
    perm_samp[,i] = samp_names[sample(1:length(samp_names),length(samp_names))]
  }
  
  ## permutation starts - recompute tau.mat  (or recompute U each time)
  sg.perm.mat <- matrix(0,nrow=ncol(dep.data),ncol=nperm)
  for (perm in 1:nperm){
    tau.mat.perm = tau.mat[perm_samp[,perm],,drop=FALSE]          # permute rows
    tau.mat.perm = tau.mat.perm[,perm_samp[,perm],drop=FALSE]     # permute columns
    for(xk in 1:length(all.sg))  sg.perm.mat[ xk , perm ]  <- get.g2stat.multin(U,tau.mat.perm)
  }
  
  
  ########################################################################
  
  #### G2 test statistic
  # *** recompute for a vector of values for each case - just reformat the result with as many rows as permutations + 1,
  # and as many columns as combinations of values of rho and mu
  
  ### Calculte G2 pval
  G2p = rowMeans( matrix(all.sg,nrow=length(all.sg),ncol=nperm+1) <= cbind(rep(Inf,length(all.sg)),sg.perm.mat ) )
  
  return (list(perm = perm_samp,G2p = G2p,Sg = all.sg))
}

