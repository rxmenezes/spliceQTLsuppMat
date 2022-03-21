# Function cor.utri
# Computes the correlation matrix between the rows of a given matrix, 
# then extracts the upper triangle of the resulting matrix, excluding the diagonal
# Result: a vector of correlations, with as many entries as nrow(dmat)*(nrow(dmat)-1)/2
# If nrow(dmat) = nsnps and ncol(dmat) = nsamples, then cor.utri(dmat) = vector of correlations between distinct snps
cor.utri <- function(dmat, method=method)
{
    cmat <- cor(t(dmat), method=method)
    cor.u <- cmat[ upper.tri(cmat) ]
    cor.u
}


# Function test.spQTL.perPart
# This function is to call test.spliceQTL using sapply, varying the index of the partition
# As it is to be used with sfSapply(), we define all foreign functions called within the function as well
# Inputs
# xk : the partition to be used
# partition: a vector with as many entries as samples in the data, containing
#   integers indicating the partitions
# snp.data: a matrix with the snp data (explanatory variables)
# yl: response variables, forming a matrix
# gexp.stand: logical, indicating if the gene expression is a matrix with SNPs on rows and samples on columns (TRUE)
#    or vice-versa; needed as input for the global test function
# nperm: the number of permutations to be used
# multin: logical, indicating whether a multinomial distribution for the random effects should
#    be used (TRUE) or a normal distribution (FALSE, the default)
# use.w: logical, indicating if the correlation matrix between the covariates in snp.data
#   should be taken into account in the test statistic (TRUE) or not (FALSE, the default)
# Output
#  The p-value computed for the given data partition
test.spQTL.perPart <- function(xk, partition, snp.data, yl,  
                               gexp.stand = T, nperm = 1000, multin = FALSE, use.w = FALSE)
{
    dp <- diag(as.numeric(partition == xk))
    dp <- dp[, colSums(dp) > 0]
    Xp <-  snp.data %*% dp
    if(use.w) {W <- tcrossprod(Xp)} else {W <- NULL}
    Yp.l <- yl %*% dp 
    
    # The following function is defined in file Projects/g2_renee.R
    get.g2stat.had <- function(W, Z)
    {
        g2tstat <- sum( W * Z )
        g2tstat
    }
    
    # The following function is defined in file Projects/g2_renee.R
    get.g2stat.multin <- function(U, tau.mat)
    {
        g2tstat <- sum( U * tau.mat )
        g2tstat
    }
    
    # The following function is defined in file Projects/f_spliceQTLfit.R
    test.spliceQTL <- function(data.exp, data.snp, gexp.stand=NULL, nperm=100, multin=FALSE, W = NULL)
    {
        if(gexp.stand) { total.gexp <- matrix(colSums(data.exp),nrow=nrow(data.exp),ncol=ncol(data.exp),byrow=T)
        total.gexp[,colSums(data.exp)==0] <- 10^(-3)
        data.exp <- data.exp/total.gexp
        } else { data.exp <- data.exp  }
        comp.snp.obs <- apply(is.na(data.snp),2,sum)==0
        data.snp <- data.snp[, comp.snp.obs  ]
        if(ncol(data.snp) > 0){
            data.snp <- data.snp[, comp.snp.obs ] 
            if(!multin) { pval <- G2(t(data.exp),t(data.snp),nperm=nperm,stand=FALSE)$G2p 
            } else { pval <- G2.multin(t(data.exp),t(data.snp),nperm=nperm,stand=FALSE, W=W)$G2p
            }
            myres <- c(nrow(data.exp),nrow(data.snp),pval)
        }
    }
    # The following function is defined in file Projects/g2_renee.R
    
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
    
    # The following function is defined in file Projects/g2_renee.R
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
    
    
    pval.res <- test.spliceQTL(  Yp.l,  Xp, gexp.stand = gexp.stand, nperm = nperm, multin = multin, W = W) 
    pval.res[3]
}