# graphs_allResults_PerGene_short_multinW_1graphPerPage.R

# This script is called by main_fit_eqtl_factors.R, itself a part of the analysis in main_analysis_spliceQTL_forVM.R
# The script with more graphs per gene can be found in graphs_allResults_PerGene.R

#load(paste(mydir.output,"/densities_correlations_betweenSNPs_matrix_perGene_x.RData",sep="")) # dens.x
#load(paste(mydir.output,"/densities_correlations_betweenSNPs_matrix_perGene_y.RData",sep="")) # dens.y
# Computing mean density too
#mean.dens <- apply(dens.y,1,mean)

for(xi in list.sel)
{
 ###
 # Selecting data corresponding to gene
 ###
  sel.gexp <- mat.gexp.sel[ xi ,]
  sel.snp <- mat.snp.sel[ xi ,]
  gexp.sel <- data.exp.only[ sel.gexp[1] : sel.gexp[2] , ]
  exon.ann.sel <- exon.ann1[ sel.gexp[1] : sel.gexp[2] , ]
  snp.sel <- mat.snp.only[ sel.snp[1] : sel.snp[2] , ]
  my.ann <- data.snp.ann[ sel.snp[1] : sel.snp[2] , ]
  g2.pval <- aspl.overl[xi,]
  pvals.exons <- results.aspl.factors[[xi]]
  cols.sign <- c("blue","gray")
  col.pvals.exons <- as.character( cut(pvals.exons,breaks=c(0,cutoff.exon,1),labels= cols.sign ) )
  sel.exons <- pvals.exons <= cutoff.exon
  ###
  pdf(paste0(mydir.output.perGene,"/graphs_gene_",xi,"_p1_exonExps.pdf"),width=10,height=6)
 ###
 # Exon expression data, per gene
 ###
  mylength <- factor( rep( rownames(gexp.sel) ,ncol(gexp.sel) ) )
  mypos <- rep( exon.ann.sel$exon.start, ncol(gexp.sel) )
  myexp <- matrix( gexp.sel,nrow=nrow(gexp.sel)*(ncol(gexp.sel)),ncol=1)[,1]
  denscols <- densCols(as.numeric(mylength),myexp)
 # equally spaced
  plot(as.numeric(mylength),myexp,pch=20,col=denscols,main=paste("Gene",xi,"G2 p-value=",round(g2.pval[3],5)),xlab="",ylab="exon expression (vst)") #,cex.axis=.4)
  points( (1:nrow(gexp.sel))[sel.exons] ,rep(max(myexp),sum(sel.exons)),pch='*',col="red")
 dev.off()
 pdf(paste0(mydir.output.perGene,"/graphs_gene_",xi,"_p2_exonCors.pdf"),width=10,height=6)
 ###
 # Correlations between exons, no dendrogram, so exons are displayed in genomic order
 ###
  heatmap.2(cor(t(gexp.sel)),symm=T,Rowv=F,Colv=F, col=bluered,density.info="density", RowSideColors=col.pvals.exons, dendrogram="none",
            trace="none",main=paste("Exons' Pearson correlations Gene",xi),cexRow=0.7,labRow=exon.ann.sel$exon.start,labCol=1:nrow(gexp.sel), # rep("",nrow(gexp.sel)), 
            symkey=T,symbreaks=T) #breaks=seq(from=-1,to=1,by=0.01) ) 
 legend("bottomleft",legend=c("sign","not sign"),fill=cols.sign,cex=.6)
 
 dev.off()
 #source("plots_maf_and_density_SNPcor.R")
 pdf(paste0(mydir.output.perGene,"/graphs_gene_",xi,"_p3_exonAndSNPCors.pdf"),width=10,height=6)
 
 ###
 # heatmap correlations between SNPs and exons
 ###
 # I want to make a variable as long as my.ann$start, indicating which start values are closest to the exon.ann values
 # alternatively, indicate which snps are within the gene

  gene.pos <- range(exon.ann.sel$exon.start)
  snp.in.gene <- (my.ann$POS >= gene.pos[1]) & (my.ann$POS <= gene.pos[2])
  # snp.in.gene should have only 0,1, 
  snp.in.gene[snp.in.gene] <- "black"
  snp.in.gene[snp.in.gene != "black"] <- "gray"

  snp.exon.cor <- cor(t(gexp.sel),t(snp.sel),method="s",use="complete.obs")
  #  colnames(snp.exon.cor) <- rep("",ncol(snp.exon.cor))
  #  rownames(snp.exon.cor) <- exon.ann.sel$exon.start
  heatmap.2(snp.exon.cor,Rowv=F,Colv=F, col=bluered,density.info="density",RowSideColors=col.pvals.exons,dendrogram="none",
            ColSideColors= snp.in.gene, trace="none",main=paste("Gene",xi, ", Spearman correlations"),cexRow=0.5,labCol=rep("",ncol(snp.exon.cor)),
            labRow=exon.ann.sel$exon.start, symkey=T,symbreaks=T) #breaks=seq(from=-1,to=1,by=0.01) )   
 legend("bottomleft",legend=c("sign","not sign"),fill=cols.sign,cex=.6)
 legend("topright",legend=c("within gene","outside"),fill=c("black","gray"),cex=.6)
 dev.off() 
}




