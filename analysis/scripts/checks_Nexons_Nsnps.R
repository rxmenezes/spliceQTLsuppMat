# checks_Nexons_Nsnps.R

max(mat.snp.sel[,2]-mat.snp.sel[,1]+1,na.rm=T) # 7303 using a big window, but after leaving out SNPs with MAF<0.05
 # 13563 -- 31381 on 2014-04-29, with a bigger window, and 26904 when NAs were left out  - on 20150428 I got 6111, possibly because of re-coding

# According to the selection matrices I computed, some genes have a v. large number of exons
mat.gexp.sel[mat.gexp.sel[,2]-mat.gexp.sel[,1]>100,]
#ENSG00000127481.9  5443 5547  # this leads to 105
#ENSG00000127603.17 5563 5665  # this leads to 103 
strange.genes <- rownames(mat.gexp.sel[mat.gexp.sel[,2]-mat.gexp.sel[,1]>100,])
sum(as.character(data.exp1$Gene_Symbol)==strange.genes[1])  # 105
sum(as.character(data.exp1$Gene_Symbol)==strange.genes[2])  # 103
# So these seem to be genuinely genes with many exons

results.altspl <- data.frame(results.altspl,n.exons.new = mat.gexp.sel[,2]-mat.gexp.sel[,1]+1, n.snps = mat.snp.sel[,2]-mat.snp.sel[,1]+1)

pdf(paste0(mydir.output, "/number_SNPs_per_gene_all.pdf"),width=10,height=6)
  plot(results.altspl$N.exons,results.altspl$n.exons.new,pch=20,col="blue",main="N.exons computed in two ways",xlab="number of rows with GeneID",
       ylab="number of cons rows")
  par(mfrow=c(1,2))
  hist(results.altspl$n.snps,col="blue",main="Number of SNPs per gene",breaks=50)
  hist(results.altspl$N.exons,col="red",main="Number of exons per gene",breaks=50)
  par(mfrow=c(1,1))
  mycols <- densCols(results.altspl$n.snps,results.altspl$n.exons.new)
  plot(results.altspl$n.snps,results.altspl$n.exons.new, pch=20,
   col=mycols,main="Number of exons and SNPs mapping to each gene (chr 1)", xlab="number of SNPs",ylab="number of exons")
  segments(0,0,4000,4000,lty="dashed",col="gray")
  plot(results.altspl$n.snps,results.altspl$n.exons.new, pch=20,log="xy",
   col="blue",main="Number of exons and SNPs mapping to each gene (chr 1)", xlab="number of SNPs",ylab="number of exons")
  segments(1,1,10000,10000,lty="dashed",col="gray")
dev.off()

save(results.altspl, file=paste0(mydir.output,"/Gene_nexons_nsnps_chr1_all.RData")) # results.altspl

