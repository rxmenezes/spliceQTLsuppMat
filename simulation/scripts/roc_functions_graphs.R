# roc_functions_graphs.R

###
# freq.fun
# Function to compute sensitivity values, written to be used with lapply
# It computes the frequency with which a vector of p-values is smaller than a scalar sign level alpha
# Input
# xi: the significance level to be chosen
# p : a vector of pvalues; 
#     this vector is already a selection of p-values for features known to be significant (to give the sensitivity)
#     or is a selection of p-values for features known to not be significant (to give the specificity)
# alpha: a vector with significance levels to be used to declare p-values significant (if p <= alpha ) or not (if p > alpha)
# Output
# a vector (same length as alpha) of values with the frequency at which the chosen p-value is smaller than or equal to each alpha value
###

freq.fun <- function(xi,p,alpha){
  sum( p <= alpha[xi] )
}

###
# roc.comp
# Function to compute values to produce a ROC curve
# Input
# p : a vector of pvalues
# effect: a binary vector with the same length as p, indicating which entries are under H0 (=0), and which are under Ha (=1)
# alpha: a vector with significance levels to be used to declare p-values significant (if p <= alpha ) or not (if p > alpha)
###

roc.comp <- function(p,effect,alpha=seq(0,1,by=0.001)){
  if (length(p) != length(effect)){
    stop("the number of features in p must be the same as in effect","\n")
  }
  total.sign <- sum(effect)
  total.nsign <- length(effect) - sum(effect)
  sens <- sapply(1:length(alpha),freq.fun, p=p[effect==1], alpha=alpha) / sum(effect==1)
  fpr <- sapply(1:length(alpha),freq.fun, p=p[effect==0], alpha=alpha) / sum(effect==0)
  mat.res <- cbind(sens,fpr)
  mat.res <- mat.res[order(mat.res[,2]),]
  mat.res
}


###
# roc.graph
# Function to produce ROC curves based on many different tests
# Input
# list.ss : a list of results from roc.comp, each element of which containing a matrix with 2 columns:  sens, spec
#           the names of the list elements are used to produce an informative legend
# mycols: colours to be used for the curves, must have the same length as list.ss
# pdf: logical indicating whether or not a pdf is to be produced
# file.name: if pdf is true, this is used as the file name
# Output
# The graph produced
###

roc.graph <- function(list.ss,mycols=NULL,my.xlim=c(0,1),my.ylim=c(0,1),main.title=NULL, pdf=TRUE,file.name="ROC curve.pdf"){
  ncurves <- length(list.ss)
  if((!is.null(mycols)) & (ncurves!=length(mycols))) {
    stop("The number of colours given must be the same as the number of curves")
  }
  if(!is.logical(pdf)){
    stop("The argument pdf must be logical")
  }
  if(is.null(mycols)) mycols <- rainbow(ncurves,start=0.1,end=0.9)
  mynames <- names(list.ss)
  if(is.null(main.title)) main.title="ROC curve"
  if(pdf) pdf(file.name,width=8,height=6)
    xi <- 1
   plot(list.ss[[xi]][,2],list.ss[[xi]][,1],type="l",col=mycols[xi],xlim=my.xlim,ylim=my.ylim,main=main.title,xlab="1-specificity",ylab="sensitivity")
   for(xi in 2:ncurves) lines( list.ss[[xi]][,2],list.ss[[xi]][,1],col=mycols[xi] )
   segments(0,0,1,1,lty="dashed",col="gray",lwd=2)
   legend("bottomright",legend=mynames,lty="solid",col=mycols)
  if(pdf) dev.off()
}
# 