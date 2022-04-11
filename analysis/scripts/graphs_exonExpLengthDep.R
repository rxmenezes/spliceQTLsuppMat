
# Define variables needed

exon.length <- exon.ann1$exon.length
gene <- factor(as.character(data.exp1$Gene_Symbol))
data.only <- as.matrix(data.exp1[, -c(1:3,ncol(data.exp1))])

# Make matrices into vectors

data.vec <- matrix(data.only, nrow = nrow(data.only)*ncol(data.only), ncol=1)[, 1]
elength.vec <- rep(exon.length, ncol(data.only))
gene.vec <- factor(rep(as.character(gene), ncol(data.only)))

# Fit linear mixed-effects model

res <- lmer(data.vec ~ elength.vec + (1|gene.vec))

# Make a graph of the dependence of exon expression on exon length for one gene
# This is supplementary figure 5

pdf(paste0(mydir.output, "/SFigure6.pdf"), width = 12, height = 5)

xi <-  "ENSG00000048707.8"
mydata <- data.only[gene == xi, ]
mylength <- rep( (exon.length)[gene == xi], ncol(mydata))
myexp <- matrix(as.matrix(mydata), nrow=nrow(mydata)*ncol(mydata), ncol=1)[, 1]
# compute intercept so as to make mean(yhat) equal to mean(y)
ahat.star <- mean(myexp)-summary(res)$coef[2,1]*mean(mylength)
denscols <- densCols(mylength, myexp)
plot(mylength, myexp, pch=20, col=denscols, main=paste(xi),
     xlab="exon length", ylab="exon expression (vst)")
lines(mylength[order(mylength)], 
      ahat.star + mylength[order(mylength)]*summary(res)$coef[2, 1],
      col="red")

dev.off()
