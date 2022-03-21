# f_spliceQTLfit.R

#####
# Function: splice.eqtl.fit
# This function computes the G2 test statistic for a set of responses, following a multinomial distribution
# It yields as result a p-value, which is computed using permutations
# Inputs
# xi: integer indicating which element to choose from a gene list 
#     - or, in this case, the row of the selection matrices to choose
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# exp.sel: matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.exp corresponding to the gene
# snp.sel:  matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.snp corresponding to the gene
# gexp.stand: logical indicating whether or not the gene expression data must be standardized 
#             by, per gene, dividing each column by its total 
#             This is especially important if we do not wish to find eQTL effects, as it corrects for the total gene expression
# nperm: the number of permutations to be used
#
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)
####

splice.eqtl.fit <- function(xi,data.exp,data.snp,exp.sel,snp.sel,gexp.stand=TRUE,nperm=100)
{
  sel.gexp <- data.exp[exp.sel[xi,1] : exp.sel[xi,2], ]
  # Standardize gene expression data by dividing exon expressions per sample by its sample gene expression
  if(gexp.stand) { total.gexp <- matrix(colSums(sel.gexp),nrow=nrow(sel.gexp),ncol=ncol(sel.gexp),byrow=T)
                   sel.gexp <- sel.gexp/total.gexp
                 }
  if(sum(is.na( snp.sel[xi,] )) == 0){
    sel.snp <- data.snp[snp.sel[xi,1] : snp.sel[xi,2], ]
    comp.snp.obs <- apply(is.na(sel.snp),2,sum)==0
    sel.snp <- sel.snp[, comp.snp.obs  ]
    if(ncol(sel.snp) > 0){
      sel.gexp <- sel.gexp[, comp.snp.obs ] 
      pval <-G2(t(sel.gexp),t(sel.snp),nperm=nperm,stand=FALSE)$G2p
      myres <- c(nrow(sel.gexp),nrow(sel.snp),pval)
    }
  } else { myres <- c(nrow(sel.gexp),NA,NA) }
myres
}

#####
# Function: splice.eqtl.fit2
# This function computes the G2 test statistic for a set of responses, following a multinomial distribution
# It yields as result a p-value, which is computed using permutations
# It differs from splice.eqtl.fit only in that it selects the rows of the data.exp and of the selection matrices
# using the gene ids, rather than assuming that the gene ids and the data matrices are in the same order
# This can be used if the list of genes used is smaller than the data, so containing a selection of genes to be tested
# Inputs
# xi: integer indicating which element to choose from a gene list 
#     - or, in this case, the row of the selection matrices to choose
# genes.used : the ids of genes to be used - these ids are contained in the rownames of the selection  matrices
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# exp.sel: matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.exp corresponding to the gene
#          The rownames of this selection matrix must be equal to the gene IDs
# snp.sel:  matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.snp corresponding to the gene
#          The rownames of this selection matrix must be equal to the gene IDs, and the same as those of exp.sel
# gexp.stand: logical indicating whether or not the gene expression data must be standardized 
#             by, per gene, dividing each column by its total 
#             This is especially important if we do not wish to find eQTL effects, as it corrects for the total gene expression
# nperm: the number of permutations to be used
# W: a square matrix with as many rows as there are covariates in the independent data set
#    It represents the correlation structure expected in the independent data
# w.type: if W=NULL and w.type is a given string, take the type: "cov" is the only one allowed. Then the 
# inner product of the indep data is taken
# multin: LOGICAL: indicate whether the response is assumed to follow a normal distribution & identity link (default, multin=F) 
#                  or a multinomial distribution and multinomial logistic link (multin=T)
#
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)
####

splice.eqtl.fit2 <- function(xi,genes.used,data.exp,data.snp,exp.sel,snp.sel,gexp.stand=TRUE,nperm=100,multin=FALSE, W=NULL ,w.type=NULL)
{
  exp.sel <- exp.sel[ rownames(exp.sel) == genes.used[xi] , ]
  snp.sel <- snp.sel[ rownames(snp.sel) == genes.used[xi] , ]
  sel.gexp <- data.exp[exp.sel[1] : exp.sel[2], ]
  # Standardize gene expression data by dividing exon expressions per sample by its sample gene expression
  if(gexp.stand) { total.gexp <- matrix(colSums(sel.gexp),nrow=nrow(sel.gexp),ncol=ncol(sel.gexp),byrow=T)
                   sel.gexp <- sel.gexp/total.gexp
  }
  if(sum(is.na( snp.sel )) == 0){
    sel.snp <- data.snp[snp.sel[1] : snp.sel[2], ]
    comp.snp.obs <- apply(is.na(sel.snp),2,sum)==0
    sel.snp <- sel.snp[, comp.snp.obs  ]
    if(ncol(sel.snp) > 0){
      sel.gexp <- sel.gexp[, comp.snp.obs ] 
      if(is.null(W)) { if(is.null(w.type)) { W <- diag(nrow(sel.snp)) } else {
        if(w.type=="cov") { W <- tcrossprod(sel.snp) }
      }  }
      if(!multin) {pval <- G2(t(sel.gexp),t(sel.snp),nperm=nperm,stand=FALSE)$G2p 
      } else {pval <- G2.multin(t(sel.gexp),t(sel.snp),nperm=nperm,stand=FALSE, W=W)$G2p
      }
      myres <- c(nrow(sel.gexp),nrow(sel.snp),pval)
    }
  } else { myres <- c(nrow(sel.gexp),NA,NA) }
  myres
}

#####
# Function: get.cor.test
# This function computes the correlation p-value between expression for a single probe, and a single covariate
# It has been written to be used with sapply
# It is called by splice.eqtl.lreg
# Inputs
# xi: integer indicating which element to choose from the covariates
# sel.gexp: the vector of responses, of length = sample size
# sel.snp: the matrix of covariates, with one covariate per row and as many rows as the sample size
# what: either "lr", in which case a linear regression is fitted, or "cor", in which case a correlation test is performed
# method: the correlation method to be used
# Output
# the p-value for the correlation test
####

get.cor.test <- function(xi,sel.gexp,sel.snp,what=c("lr","cor"),method=NULL)
{
  what <- match.arg(what)
  mydata.snp <- sel.snp[xi,]
  if(what=="cor") { method <- match.arg(method,c("pearson", "kendall", "spearman"))  
                    myp <- cor.test(sel.gexp,mydata.snp,method=method)$p.value
  } else {
    myp <- summary(lm(sel.gexp~mydata.snp))$coef[2,4]
  }
  myp
}


#####
# Function: splice.eqtl.lreg
# This function computes p-values for each pair of variables in the expression set and in the snp set
# It extracts the p-value by performing either a linear regression or a correlation test
# Inputs
# xi: integer indicating which element to choose from a gene list 
#     - or, in this case, the row of the selection matrices to choose
# genes.used : the ids of genes to be used - these ids are contained in the rownames of the selection  matrices
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# exp.sel: matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.exp corresponding to the gene
#          The rownames of this selection matrix must be equal to the gene IDs
# snp.sel:  matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.snp corresponding to the gene
#          The rownames of this selection matrix must be equal to the gene IDs, and the same as those of exp.sel
# gexp.stand: logical indicating whether or not the gene expression data must be standardized 
#             by, per gene, dividing each column by its total 
#             This is especially important if we do not wish to find eQTL effects, as it corrects for the total gene expression
#             I use here the same default as for the spliceQTL test, for consistency
# what: what is to be used to test for association. Options are: "lr", linear regression, or "cor", a correlation test
# method: if the correlation is chosen, the correlation method to be used 
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)
####

splice.eqtl.lreg <- function(xi,genes.used,data.exp,data.snp,exp.sel,snp.sel,gexp.stand=TRUE,what=c("lr","cor"),method=NULL,mydir)
{
  what=match.arg(what)
  exp.sel <- exp.sel[ rownames(exp.sel) == genes.used[xi] , ]
  snp.sel <- snp.sel[ rownames(snp.sel) == genes.used[xi] , ]
  sel.gexp <- data.exp[exp.sel[1] : exp.sel[2], ]
  if(gexp.stand) { total.gexp <- matrix(colSums(sel.gexp),nrow=nrow(sel.gexp),ncol=ncol(sel.gexp),byrow=T)
                   sel.gexp <- sel.gexp/total.gexp
  }
  if(sum(is.na( snp.sel )) == 0){
    sel.snp <- data.snp[snp.sel[1] : snp.sel[2], ]
    comp.snp.obs <- apply(is.na(sel.snp),2,sum)==0
    sel.snp <- sel.snp[, comp.snp.obs  ]
    if(ncol(sel.snp) > 0){
      sel.gexp <- sel.gexp[, comp.snp.obs ] 
      pvals.mat <- matrix(1,nrow=nrow(sel.gexp),ncol=nrow(sel.snp))
      for(xj in 1:nrow(sel.gexp))  
        {
        mygexp <- sel.gexp[xj,]
        pvals.mat[xj,] <- sapply(1:nrow(sel.snp),get.cor.test,mygexp,sel.snp,what=what,method=method)
      }
    }
  } else { pvals.mat <- NA }
  save(pvals.mat,file=paste(mydir,"/pvalues_mat_",genes.used[xi],"_",what,".RData",sep="")) # pvals.mat
}

#####
# Function: call.gt
# This function is used by splice.qtl.factor to call the gt function
# It assumes that the response follows a normal distribution
# It applies gt to the response chosen, using the gene set given
# Inputs:
# xi: the row of the response to be used
# dep.data: the dependent data matrix, with variables on rows and samples on columns
# indep.data: the independent data matrix containing the variables to compose the linear predictor under the alternative
#             it has variables on rows and samples on columns
#             the number of columns in indep.data must be the same as the one in dep.data
# what: what the function is supposed to return - either one of "z" (standardized test statistics for the responses), 
#                                                 "p" (p-values for the responses), or 
#                                                 "s" for the z-scores of the covariates involved -- note that the latter produces a vector per response
# makePlot: logical, indicating whether or not the covariates plot is to be produced
# Output: either a test statistic or a p-value, resulting from the global test

call.gt <- function(xi, dep.data, indep.data, what=c("z","p","s"),makePlot=FALSE)
{
  mystat <- what
  yvar <- dep.data[xi,]
  xmat <- t(as.matrix(indep.data))
  colnames(xmat) <- paste("c",1:ncol(xmat))
  mytest <- gt(yvar,xmat)
  if(is.null(what)) what == "z"
  if(what == "p") {myres <- p.value(mytest)
                  } else if(what == "z") {
                   myres <- z.score(mytest)
                  } else { myres <- result(covariates(mytest,what="z-score",plot=makePlot,sort=FALSE,cluster=FALSE)) }
  if(makePlot) covariates(mytest,what=what)
  myres
}
  
#####
# Function: call.gt.cov
# This function is used by splice.qtl.pairs to call the gt function
# It assumes that the response follows a normal distribution
# It applies gt to the response chosen, for the set of covariates given, and returns
# results for each covariate individually
# Inputs:
# xi: the row of the response to be used
# dep.data: the dependent data matrix, with variables on rows and samples on columns
# indep.data: the independent data matrix containing the variables to compose the linear predictor under the alternative
#             it has variables on rows and samples on columns
#             the number of columns in indep.data must be the same as the one in dep.data
# what: what the function is supposed to return - either one of "z" (standardized test statistics for the responses), or
#                                                 "p" (p-values for the responses), or
#                                                 "s" for the z-scores of the covariates involved -- note that the latter produces a vector per response
# makePlot: logical, indicating whether or not the covariates plot is to be produced
# Output: either a test statistic or a p-value, resulting from the global test

call.gt.cov <- function(xi, dep.data, indep.data, what=c("z","p","s"),makePlot=FALSE)
{
  yvar <- dep.data[xi,]
  xmat <- t(as.matrix(indep.data))
  colnames(xmat) <- paste("c",1:ncol(xmat))
  myres <- result( covariates(gt(yvar,xmat),what=what,plot=makePlot,sort=FALSE,cluster=FALSE) ) 
  if(what == "p") { myres <- myres[,1] } else {
    if(what == "s") { myres <- myres[,2] } else {
      myres <- myres[,2]/myres[,4]
    }
  }
  myres
}

  

#####
# Function: splice.eqtl.factor
# This function computes test statistics per response, so factorizing the G2 test statistic
# It involves computing the global test per response, yielding either test statistics or p-values
# Here it is assumed that the responses have a normal distribution
# Inputs
# xi: integer indicating which element to choose from a gene list 
#     - or, in this case, the row of the selection matrices to choose
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# exp.sel: matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.exp corresponding to the gene
# snp.sel:  matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.snp corresponding to the gene
# gexp.stand: logical indicating whether or not the gene expression data must be standardized 
#             by, per gene, dividing each column by its total 
#             This is especially important if we do not wish to find eQTL effects, as it corrects for the total gene expression
# what : what the function is supposed to return - either one of "t" (test statistics) or "p" (p-values)
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)
####

splice.eqtl.factor <- function(xi,data.exp,data.snp,exp.sel,snp.sel,gexp.stand=TRUE, what ="z",
                               makePlot=FALSE,file.name = paste("gt_covariates_Dep",xi,".pdf"))
{
  require(globaltest) || stop("this function requires the globaltest package but it is not available")
  if(missing(what)) { what <- "z"} else {
    match.arg(tolower(what), c("z", "p","s"))
  } # closes else of "what"

  sel.gexp <- data.exp[exp.sel[xi,1] : exp.sel[xi,2], ]
  # Standardize gene expression data by dividing exon expressions per sample by its sample gene expression
  if(gexp.stand) { total.gexp <- matrix(colSums(sel.gexp),nrow=nrow(sel.gexp),ncol=ncol(sel.gexp),byrow=T)
                   sel.gexp <- sel.gexp/total.gexp
                 }
  if(sum(is.na( snp.sel[xi,] )) == 0){  # this "if" closes with "else { myres <- rep(NA,nrows(sel.gexp)) }"
    sel.snp <- data.snp[snp.sel[xi,1] : snp.sel[xi,2], ]
    comp.snp.obs <- apply(is.na(sel.snp),2,sum)==0
    sel.snp <- sel.snp[, comp.snp.obs  ]
    
    if(ncol(sel.snp) > 0){
      sel.gexp <- sel.gexp[, comp.snp.obs ] 
      if(!makePlot){ 
        myres <- sapply( 1:nrow(sel.gexp) , call.gt , dep.data = sel.gexp , indep.data = sel.snp , what=what)
      } else {
        pdf(file.name,width=10,height=6)
        myres <- sapply( 1:nrow(sel.gexp) , call.gt , dep.data = sel.gexp , indep.data = sel.snp , what=what,
                         makePlot=TRUE,file.name = file.name)
        dev.off()
      }  #v closes "else" of "!makePlot"
     } # closes "if" of "ncol(sel.snp>0)"
   
  } else { myres <- rep(NA,nrow(sel.gexp)) }
myres
}

#####
# Function: splice.eqtl.factor.w
# This function computes test statistics per response, so factorizing the G2 test statistic
# It involves computing the multin G2 test per response, yielding either test statistics or p-values
# Here it is assumed that the responses have a multinomial distribution
# Inputs
# xi: integer indicating which element to choose from a gene list 
#     - or, in this case, the row of the selection matrices to choose
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# exp.sel: matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.exp corresponding to the gene
# snp.sel:  matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.snp corresponding to the gene
# gexp.stand: logical indicating whether or not the gene expression data must be standardized 
#             by, per gene, dividing each column by its total 
#             This is especially important if we do not wish to find eQTL effects, as it corrects for the total gene expression
# W : a covariance matrix between the effects (optional)
# w.type : either NULL (no covariance is taken into account if W is also NULL), or "cov", in which case
#          the covariance between effects is taken as the covariance between covariates
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)
####


splice.eqtl.factor.w <- function(xi,data.exp,data.snp,exp.sel,snp.sel,gexp.stand=TRUE, W=NULL ,w.type=NULL,nperm=100)
{
#  require(globaltest) || stop("this function requires the globaltest package but it is not available")
  
  sel.gexp <- data.exp[exp.sel[xi,1] : exp.sel[xi,2], ]
  # Standardize gene expression data by dividing exon expressions per sample by its sample gene expression
  if(gexp.stand) { total.gexp <- matrix(colSums(sel.gexp),nrow=nrow(sel.gexp),ncol=ncol(sel.gexp),byrow=T)
                   sel.gexp <- sel.gexp/total.gexp
  }
  if(sum(is.na( snp.sel[xi,] )) == 0){  # this "if" closes with "else { myres <- rep(NA,nrows(sel.gexp)) }"
    sel.snp <- data.snp[snp.sel[xi,1] : snp.sel[xi,2], ]
    comp.snp.obs <- apply(is.na(sel.snp),2,sum)==0
    sel.snp <- sel.snp[, comp.snp.obs  ]
    
    if(ncol(sel.snp) > 0){
      sel.gexp <- sel.gexp[, comp.snp.obs ] 
      if((w.type=="cov")&(is.null(W)) ) { W <- tcrossprod(sel.snp) }      
          myres <- G2.multin.factor(t(sel.gexp),t(sel.snp),nperm=nperm,stand=FALSE, W=W)$G2p
                   
    } # closes "if" of "ncol(sel.snp>0)"  
  } else { myres <- rep(NA,nrow(sel.gexp)) }
  myres
}




#####
# Function: splice.eqtl.factor2
# This function does the same as splice.eqtl.factor, but instead of selecting the SNPs that correspond to the gene,
# it loads a workspace with that SNP set
# For details about input and output,see splice.eqtl.factor
#####

splice.eqtl.factor2 <- function(xi,data.exp,exp.sel,gexp.stand=TRUE, what =c("z","p"),
                               makePlot=FALSE,file.name = paste("gt_covariates_Dep",xi,".pdf"))
{
  require(globaltest) || stop("this function requires the globaltest package but it is not available")
  if(missing(what)) { what <- "z"} else {
    match.arg(tolower(what), c("z", "p"))
  }
  sel.gexp <- data.exp[exp.sel[xi,1] : exp.sel[xi,2], ]
  # Standardize gene expression data by dividing exon expressions per sample by its sample gene expression
  if(gexp.stand) { total.gexp <- matrix(colSums(sel.gexp),nrow=nrow(sel.gexp),ncol=ncol(sel.gexp),byrow=T)
                   sel.gexp <- sel.gexp/total.gexp
  }
  load(paste(mydir.output,"/data_snp_geneN_",xi,".RData",sep="")) # sel.snp
  if(!makePlot){ 
        myres <- sapply( 1:nrow(sel.gexp) , call.gt , dep.data = sel.gexp , indep.data = sel.snp , what=what)
      } else {
        pdf(file.name,width=10,height=6)
        myres <- sapply( 1:nrow(sel.gexp) , call.gt , dep.data = sel.gexp , indep.data = sel.snp , what=what,
                         makePlot=TRUE,file.name = file.name)
        dev.off()
      }
    
  myres
}


#####
# Function: splice.eqtl.pairs
# This function computes test statistics per pair of response and covariate
# It is meant to be used after G2 or spliceQTL test, and after that has been applied per exon for all covariates
# It involves computing the global test per pair of response and covariate, yielding either test statistics or p-values
# Here it is assumed that the responses have a normal distribution
# Inputs
# xi: integer indicating which element to choose from a gene list 
#     - or, in this case, the row of the selection matrices to choose
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# exp.sel: matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.exp corresponding to the gene
# snp.sel:  matrix of integers, with as many rows as genes, and two columns
#          Per row, the indices for the first and last rows of data.snp corresponding to the gene
# results: list containing results from a previous run, possibly with another test, for each gene considering each exon separately
#         In general, this is the result of tests run per response set, and considering one response at a time
#         This indicates which of the responses drive the results for the response set - only those responses are going to be considered
#         We expect that each element of the list is a set of p-values, one p-value per response in each response set
# what : what the function is supposed to return - either one of "t" (test statistics) or "p" (p-values)
# cutoff: scalar - the cutoff to be used to declare the individual responses in each set significant
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)
####

splice.eqtl.pairs <- function(xi,data.exp,data.snp,exp.sel,snp.sel,snp.ann=NULL,results, what ="z",cutoff=0.01,
                               makePlot=FALSE,file.name = paste("gt_covariates_Dep",xi,".pdf"))
{
  require(globaltest) || stop("this function requires the globaltest package but it is not available")
  if(missing(what)) { what <- "z"} else {
    match.arg(tolower(what), c("z", "p","s"))
  } # closes else of "what"
  # Check if there are any exons that are individually significant; if none, finish
  myresults <- results[[xi]]
  sel.exon <- myresults <= cutoff.exon 
  if( sum(sel.exon) > 0 ) {  
    exp.sel.exon <- (exp.sel[xi,1] : exp.sel[xi,2])[sel.exon]
    sel.gexp <- data.exp[exp.sel.exon, ,drop=FALSE]
  if( sum(is.na( snp.sel[xi,] )) == 0 ){  # this "if" closes with "else { myres <- rep(NA,nrows(sel.gexp)) }"
    sel.snp <- data.snp[ snp.sel[xi,1] : snp.sel[xi,2] ,]
    comp.snp.obs <- apply(is.na(sel.snp),2,sum)==0
    sel.snp <- sel.snp[, comp.snp.obs  ]
    
    if(ncol(sel.snp) > 0){
      sel.gexp <- sel.gexp[, comp.snp.obs , drop=FALSE ] 
      if(!makePlot){ 
        myres <- sapply( 1:nrow(sel.gexp) , call.gt.cov , dep.data = sel.gexp , indep.data = sel.snp , what=what)
      } else {
        pdf(file.name,width=10,height=6)
        myres <- sapply( 1:nrow(sel.gexp) , call.gt.cov , dep.data = sel.gexp , indep.data = sel.snp , what=what,
                         makePlot=TRUE,file.name = file.name)
        dev.off()
      }  #v closes "else" of "!makePlot"
      
    } # closes "if" of "ncol(sel.snp>0)"
  }
  } else { myres <- NA }
  # Putting row and column names from the annotation
  if(is.null(dim(myres))) { myres <- matrix(myres,ncol=1) }
  colnames(myres) <- rownames(sel.gexp)                                  
#  if(!is.null(snp.ann)) { 
#    snp.ann.sel <- snp.ann[ snp.sel[xi,1] : snp.sel[xi,2] , ] 
#    rownames(myres) <- snp.ann.sel$ID
#  }  
  myres
}

#
# test.spliceQTL
# Function to test for spliceQTL
# Can take the normal or the multinomial
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# gexp.stand: logical indicating whether or not the gene expression data must be standardized 
#             by, per gene, dividing each column by its total 
#             This is especially important if we do not wish to find eQTL effects, as it corrects for the total gene expression
# nperm: the number of permutations to be used
#        observed for different covariates (SNPs); only used for the multinomial
# multin: LOGICAL: indicate whether the response is assumed to follow a normal distribution & identity link (default, multin=F) 
#                  or a multinomial distribution and multinomial logistic link (multin=T)
#
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)

test.spliceQTL <- function(data.exp,data.snp,gexp.stand=NULL,nperm=100,multin=FALSE, W=NULL)
  {
    if(gexp.stand) { total.gexp <- matrix(colSums(data.exp),nrow=nrow(data.exp),ncol=ncol(data.exp),byrow=T)
                     total.gexp[,colSums(data.exp)==0] <- 10^(-3)
                     data.exp <- data.exp/total.gexp
          } else { data.exp <- data.exp  }
  comp.snp.obs <- apply(is.na(data.snp),2,sum)==0
  data.snp <- data.snp[, comp.snp.obs  ]
  if(ncol(data.snp) > 0){
    data.snp <- data.snp[, comp.snp.obs ] 
    if(!multin) {pval <- G2(t(data.exp),t(data.snp),nperm=nperm,stand=FALSE)$G2p 
    } else {pval <- G2.multin(t(data.exp),t(data.snp),nperm=nperm,stand=FALSE, W=W)$G2p
    }
    myres <- c(nrow(data.exp),nrow(data.snp),pval)
  }
}

#
# test.spliceQTL.old
# Here we still use the old formulation with mu and rho. We now leave those out, and use a general W instead of W.rho
# Function to test for spliceQTL
# Can take the normal or the multinomial
# data.exp: gene expression data matrix, with one row per probe (in the RP3 case, exons) and one column per sample
#           No annotation columns
# data.snp: snp genotype data matrix (numeric), with one column per sample
#           Here it is assumed that the samples (columns) are in the same order as those of data.exp
# gexp.stand: logical indicating whether or not the gene expression data must be standardized 
#             by, per gene, dividing each column by its total 
#             This is especially important if we do not wish to find eQTL effects, as it corrects for the total gene expression
# nperm: the number of permutations to be used
# rho: constant (scalar or vector) representing the correlation between effects on different responses (exons); only used for the multinomial
# mu: constant (scalar or vector) representing the proportional change in the correlation between effects 
#        observed for different covariates (SNPs); only used for the multinomial
# multin: LOGICAL: indicate whether the response is assumed to follow a normal distribution & identity link (default, multin=F) 
#                  or a multinomial distribution and multinomial logistic link (multin=T)
#
# # Note: Only samples with complete observations for the SNPs are considered (no imputation at all used)

test.spliceQTL.old <- function(data.exp,data.snp,gexp.stand=NULL,nperm=100,multin=FALSE, rho = 0, mu = 0)
{
  if(gexp.stand) { total.gexp <- matrix(colSums(data.exp),nrow=nrow(data.exp),ncol=ncol(data.exp),byrow=T)
                   total.gexp[,colSums(data.exp)==0] <- 10^(-3)
                   data.exp <- data.exp/total.gexp
  } else { data.exp <- data.exp  }
  comp.snp.obs <- apply(is.na(data.snp),2,sum)==0
  data.snp <- data.snp[, comp.snp.obs  ]
  if(ncol(data.snp) > 0){
    data.snp <- data.snp[, comp.snp.obs ] 
    if(!multin) {pval <- G2(t(data.exp),t(data.snp),nperm=nperm,stand=FALSE)$G2p 
    } else {pval <- G2.multin(t(data.exp),t(data.snp),nperm=nperm,stand=FALSE, rho = rho, mu = mu)$G2p
    }
    myres <- c(nrow(data.exp),nrow(data.snp),pval)
  }
}

####
# Function my.rmultinom
# It is written to facilitate generating multinomial responses with the data in the current format
# We can thus call it with apply and let it be applied to the columns of the expected probs
# Input
# prob: vector of probabilities, as long as the number of levels in the multinomial
#       This is standardized internally by rmultinom
# size: total number of observations
# n: number of random vectors to draw
# Output
# A vector of the same length as prob
#####
my.rmultinom <-  function(prob,size,n) 
{
  resp <- rmultinom(n=n,size=size,prob=prob)
  resp 
}

#####
## function gen.response
# We use this to generate responses that are discrete, given a set of explanatory variables and effects
# We multiply effect by data, compute probability per response
# Note that only some samples display the effect
# Input
# x : matrix of covariates (nrow=ncovariates, ncol=nsamples)
# b : matrix of effects (nrow=ncovariates, ncol=nresponses)
# v.aff : vector of column indicators (integers) defining the subset of affected samples 
# total.gexp : scalar giving the sum (fixed) of the responses to be used
# Output
# matrix with multinomial responses
####
gen.response <- function(x, b, v.aff, total.gexp)
{
  nsamples <- ncol(x)
  nresponses <- ncol(b)
  v.unaff <- (1:nsamples)[!((1:nsamples) %in% v.aff)]
  exp.prob <- matrix(0,nrow=nresponses,ncol=nsamples)
  exp.prob[,v.aff] <- 1/( 1+exp( -1* t(b) %*% x[,v.aff] ) )
  exp.prob[,v.unaff] <- matrix(1,nrow=nresponses,ncol=length(v.unaff))
  resp <- apply(exp.prob,2,my.rmultinom,size=total.gexp,n=1)
}


  
