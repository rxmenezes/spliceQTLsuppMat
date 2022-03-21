#
# This function takes a dataset (genes on rows, cases/samples on columns) and
# applies the Wilcoxon test per gene to check if there are differentially expressed genes
# Inputs
# dataset: a matrix with ngenes rows and nsamples columns
# yvar: a variable with the same length as nsamples with two values, defining the groups to be compared
# To use with genetic subtypes

wpvals.gen <- function(dataset,yvar)
{
n.yvar <- yvar[!is.na(yvar)]
f.yvar <- factor(n.yvar)
n.drug <- length(f.yvar)
subset.drug <- dataset[,!is.na(yvar)]
drug.wpvals <- rep(0,nrow(subset.drug))
for(i in 1:nrow(subset.drug))
 {
  drug.wpvals[i] <- wilcox.test(as.numeric(subset.drug[i,f.yvar==levels(f.yvar)[1]]),as.numeric(subset.drug[i,f.yvar==levels(f.yvar)[2]]))$p.value
 }
drug.wpvals
}
#
# This function takes a dataset (genes on rows, cases/samples on columns) and
# applies the Wilcoxon test per gene to check if there are differentially expressed genes
# Inputs
# dataset: a matrix with ngenes rows and nsamples columns
# yvar: a variable with the same length as nsamples with three values, defining the groups to be compared
# The intermediate value is abandoned
# To use with drug resistance status
# 30
wpvals.drugs <- function(dataset,yvar)
{
n.yvar <- yvar[!is.na(yvar)]
f.yvar <- factor(n.yvar)
f.drug2 <- f.yvar[ f.yvar!=levels(f.yvar)[2] ]
n.drug <- length(f.drug2)
subset.drug1 <- dataset[,!is.na(yvar)]
subset.drug <- subset.drug1[,f.yvar!=levels(f.yvar)[2] ]
drug.wpvals <- rep(0,n.drug)
for(i in 1:n.drug)
 {
  drug.wpvals[i] <- wilcox.test(as.numeric(subset.drug[i,f.drug2==levels(f.drug2)[1]]),as.numeric(subset.drug[i,f.drug2==levels(f.drug2)[3]]))$p.value
 }
drug.wpvals
}

#
# This function takes a dataset (genes on rows, cases/samples on columns) and
# applies the Wilcoxon test per gene to check if there are differentially expressed genes
# Given a factor, it compares samples from one level to the remaining samples
# Inputs
# dataset: a matrix with ngenes rows and nsamples columns
# yvar: a variable with the same length as nsamples with at least two values, defining the groups to be compared
# val: the value of yvar to define the group, all samples not associated with this value are in the "other" group
# To use with genetic subtypes

wpvals.pairs <- function(dataset,yvar,val)
{
n.yvar <- yvar[!is.na(yvar)]
if(is.factor(n.yvar)==F)  f.yvar <- factor(n.yvar)
   else f.yvar <- n.yvar
n.drug <- length(f.yvar)
subset.drug <- dataset[,!is.na(yvar)]
drug.wpvals <- rep(0,nrow(subset.drug))
for(i in 1:nrow(subset.drug))
 {
  drug.wpvals[i] <- wilcox.test(as.numeric(subset.drug[i,f.yvar==val]),as.numeric(subset.drug[i,f.yvar!=val]))$p.value
 }
drug.wpvals
}

#
# This function takes a dataset (genes on rows, cases/samples on columns) and
# applies the Wilcoxon test per gene to check if there are differentially expressed genes
# Given a factor, it compares samples from one level to the samples from the second level
# Inputs
# dataset: a matrix with ngenes rows and nsamples columns
# yvar: a variable with the same length as nsamples with at least two values, defining the groups to be compared
# val1: the value of yvar to define the first group
# val2: the value of yvar to define the second group 
# To use with genetic subtypes

wpvals.2groups <- function(dataset,yvar,val1,val2)
{
n.yvar <- yvar[!is.na(yvar)]
if(is.factor(n.yvar)==F)  f.yvar <- factor(n.yvar)
   else f.yvar <- n.yvar
n.drug <- length(f.yvar)
subset.drug <- dataset[,!is.na(yvar)]
drug.wpvals <- rep(0,nrow(subset.drug))
for(i in 1:nrow(subset.drug))
 {
  drug.wpvals[i] <- wilcox.test(as.numeric(subset.drug[i,f.yvar==val1]),as.numeric(subset.drug[i,f.yvar==val2]))$p.value
 }
drug.wpvals
}


#
# This function takes a dataset (genes on rows, cases/samples on columns) and
# applies the t-test test per gene to check if there are differentially expressed genes
# Given a factor, it compares samples from one level to the remaining samples
# Inputs
# dataset: a matrix with ngenes rows and nsamples columns
# yvar: a variable with the same length as nsamples with at least two values, defining the groups to be compared
# val: the value of yvar to define the group, all samples not associated with this value are in the "other" group
# To use with genetic subtypes

tpvals.pairs <- function(dataset,yvar,val)
{
n.yvar <- yvar[!is.na(yvar)]
if(is.factor(n.yvar)==F)  f.yvar <- factor(n.yvar)
   else f.yvar <- n.yvar
n.drug <- length(f.yvar)
subset.drug <- dataset[,!is.na(yvar)]
drug.wpvals <- rep(0,nrow(subset.drug))
for(i in 1:nrow(subset.drug))
 {
  drug.wpvals[i] <- t.test(as.numeric(subset.drug[i,f.yvar==val]),as.numeric(subset.drug[i,f.yvar!=val]))$p.value
 }
drug.wpvals
}
#
# This function takes a dataset (genes on rows, cases/samples on columns) and
# applies the Student-t test per gene to check if there are differentially expressed genes
# Inputs
# dataset: a matrix with ngenes rows and nsamples columns
# yvar: a variable with the same length as nsamples with two values, defining the groups to be compared
# To use with genetic subtypes

tpvals.gen <- function(dataset,yvar)
{
n.yvar <- yvar[!is.na(yvar)]
f.yvar <- factor(n.yvar)
n.drug <- length(f.yvar)
subset.drug <- dataset[,!is.na(yvar)]
drug.wpvals <- rep(0,nrow(subset.drug))
for(i in 1:nrow(subset.drug))
 {
  drug.wpvals[i] <- t.test(as.numeric(subset.drug[i,f.yvar==levels(f.yvar)[1]]),as.numeric(subset.drug[i,f.yvar==levels(f.yvar)[2]]))$p.value
 }
drug.wpvals
}

#
# This function takes a dataset (genes on rows, cases/samples on columns) and
# applies the Student-t test per gene to check if there are differentially expressed genes
# Inputs
# dataset: a matrix with ngenes rows and nsamples columns
# yvar: a variable with the same length as nsamples with two values, defining the groups to be compared
# To use with genetic subtypes

logit.pvals.gen <- function(dataset,yvar)
{
n.yvar <- yvar[!is.na(yvar)]
f.yvar <- factor(n.yvar)
n.drug <- length(f.yvar)
subset.drug <- dataset[,!is.na(yvar)]
drug.wpvals <- rep(0,nrow(subset.drug))
#for(i in 1:nrow(subset.drug))
# {
#  drug.wpvals[i] <- anova( glm(as.numeric( subset.drug[i,f.yvar==levels(f.yvar)[1]] ),as.numeric( subset.drug[i,f.yvar==levels(f.yvar)[2]] ),test="Chisq" )$P[2]
# }
drug.wpvals
}

#
# This function does the BH-FDR multiple testing adjustment and produces graphs
# It requires the multtest library to be loaded
# It uses lines for the both graphs with just a selection of the first few points for the second, 
# so it makes for smaller pdfs than adj.bhfdr2 and adj.bhfdr3
# The points displayed are in another colour, and are only for features below the given threshold
# Also, the number of features under the FDR threshold is shown
# Input
# pvals: a vector containing the p-values of the features
# myvar: a string to be used in the graphs titles
# cutoff: the cutoff for the FDR to be used

adj.bhfdr <- function(pvals,myvar,cutoff=0.05)
{
  adj.pvals <- p.adjust(pvals,method="BH")
  ngenes <- sum(!is.na(pvals))
  n.sel <- sum(adj.pvals<=cutoff,na.rm=T)

# Making graphs to visualise effect

par(mfrow=c(1,3))
hist(pvals,main="Raw p-values",xlab="")
plot(1:ngenes,sort(pvals),main="Sorted raw p-values",xlab="Genes",ylab="",type="l",col="blue")
segments(0,0,ngenes,1,lty="dashed")
plot(1:ngenes,sort(adj.pvals),main=paste("FDR",myvar),xlab="Genes",ylab="sorted p-values",type="l",col="blue",ylim=c(0,1))
if(n.sel>0) points(1:n.sel,sort(adj.pvals)[1:n.sel],pch=20,col="red")
segments(0,cutoff,ngenes,cutoff,lty="dashed")
text(0.8*ngenes,0.5,labels=paste("N. selected",n.sel),cex=.7)
if(sum(is.na(pvals))>0) text(0.8*ngenes,0.6,labels=paste("N. NAs",sum(is.na(pvals))),cex=.7)

# Output

adj.pvals
}

#
# This function does the BH-FDR multiple testing adjustment and produces graphs
# It requires the multtest library to be loaded
# It uses lines for the first graph and points for the second, 
# so it makes for smaller pdfs than adj.bhfdr3, which uses points for both
# Input
# pvals: a vector containing the p-values of the features

adj.bhfdr2 <- function(pvals,myvar)
{
adj.pvals <- p.adjust(pvals,method="BH")
ngenes <- length(pvals)

# Making graphs to visualise effect

par(mfrow=c(1,3))
hist(pvals,main="Raw p-values",xlab="")
plot(1:ngenes,sort(pvals),main="Sorted raw p-values",xlab="Genes",ylab="",type="l",col="blue")
segments(0,0,ngenes,1,lty="dashed")
plot(1:ngenes,sort(adj.pvals),main=paste("FDR",myvar),xlab="Genes",ylab="sorted p-values",pch=20,col="blue")
segments(0,0.05,ngenes,0.05,lty="dashed")

# Output

adj.pvals
}

#
# This function does the BH-FDR multiple testing adjustment and produces graphs
# It requires the multtest library to be loaded
# It uses points for both graphs, so it makes heavier files
# Input
# pvals: a vector containing the p-values of the features

adj.bhfdr3 <- function(pvals,myvar)
{
adj.pvals <- p.adjust(pvals,method="BH")
ngenes <- length(pvals)

# Making graphs to visualise effect

par(mfrow=c(1,3))
hist(pvals,main="Raw p-values",xlab="")
plot(1:ngenes,sort(pvals),main="Sorted raw p-values",xlab="Genes",ylab="",pch=20,col="blue")
segments(0,0,ngenes,1,lty="dashed")
plot(1:ngenes,sort(adj.pvals),main=paste("FDR",myvar),xlab="Genes",ylab="sorted p-values",pch=20,col="blue")
segments(0,0.05,ngenes,0.05,lty="dashed")

# Output

adj.pvals
}


#
# This function does the BH-FDR multiple testing adjustment and produces graphs
# It requires the multtest library to be loaded
# Input
# pvals: a vector containing the p-values of the features
# same as adj.bhfdr, but with no graphs

adj.bhfdr.noG <- function(pvals)
{
adj.pvals <- p.adjust(pvals,method="BH")
# Output
adj.pvals
}

#
# Function to compute a table with numbers of genes selected for various levels of error
#
# Input
# pvals: a vector with (multiple-testing adjusted) p-values

table.npvals <- function(pvals)
{
error.levels <- c(0.01,0.025,0.05,0.075,0.10)
sel.id <- (pvals<=error.levels[1])
for(i in 2:length(error.levels))
 {
  sel.id <- cbind(sel.id, (pvals<=error.levels[i]) )
 }
colnames(sel.id) <- c(0.01,0.025,0.05,0.075,0.10)
n.selected <- apply(sel.id,2,sum)
}

#
# Function to make boxplots of data matrix with clones (rows) as factor and an
# extra factor from the data (group, say)
# It is mostly useful if gene expression levels are comparable.
# See a simpler way of doing this, by putting graphs per clone separately, in balike_getPS_ttest.R
#
# Inputs:
# mytopdata: data matrix with clones on rows and samples on columns
# --- we assume that the rownames of this matrix contain the clone ids,
#     unless explictly declared
# myfactor: a factor or vector with the same number of elements as the number of columns in
#     the data matrix
# mychosenlevel: factor level (or label, or value of the vector if not a factor)
#    I want to highlight in the graphs, so new factor
#    will have two levels, namely "mychosenlevel" and "the other level" 
# --- if not declared we assume that all values of myfactor are to be considered
#     separately  (in later version)
# mylabels: factor labels the final factor should have, to be included in the plot
#    

boxplot.factor <- function(mytopdata,myfactor,mychosenlevel,mylevels)
{
mytopdata <- data.ttop
probes.t <- factor(rep(rownames(mytopdata),ncol(mytopdata)))
data.top.vec <- array(mytopdata,dim=c(nrow(mytopdata)*ncol(mytopdata),1))[,1]
fclass.b.vec <- rep(0, nrow(mytopdata)*ncol(mytopdata) )
for(i in 1:ncol(mytopdata))
{
 fclass.b.vec[ ((i-1)*nrow(mytopdata)+1):(i*nrow(mytopdata)) ] <- rep(myfactor[i],nrow(mytopdata))
}
fclass.ba <- factor(fclass.b.vec==mychosenlevel,labels=mylevels)
boxplot( as.numeric(data.top.vec) ~ probes.t + fclass.ba )
}
