# graphs_spliceQTLfactor_results_W_forPresentation.R

# script called by main_analysis_spliceQTL.R

myalpha <- 0.05

my.cut <- myalpha/( nrow(res.all)-sum(is.na(res.all[,3])) )  # 3.633721e-05 for multin_W, very similar to the previous 3.628447e-05
res.nonas <- res.all[ !(is.na(res.all[,3])), ]
list.sel <- rownames( res.nonas[ res.nonas[,3] <= my.cut ,] )

# Extracting the table from the spliceQTL results that corresponds to these genes
aspl.overl <- res.all[list.sel,]
dim(aspl.overl) # 64   3 for normal, 34  3 for multin.W

# new number of tests: number already tested + the total number of exons corresponding to the selected genes
cutoff.exon <- myalpha/( sum(aspl.overl[,1]) + sum(!is.na(res.all[,3])) ) # 2.3 x 10-5 after shift, 2.788622e-05 for multin.W
# results.aspl.factors, chr 1 only, no 1-exon genes, my.nperm=10000 - 55 genes with no shift, 66 with shift

# We now make 1 heatmap per gene with ID in the following list
list.sel <- c("ENSG00000007341.12", "ENSG00000131236.11") 
# Other gene ids used previously: "ENSG00000008128.15","ENSG00000116171.11"
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
    g2.pval <- res.all[xi, ]
    pvals.exons <- results.aspl.factors[[xi]]
    cols.sign <- c("blue","gray")
    if(!is.null(pvals.exons)) {
        col.pvals.exons <- as.character( cut(pvals.exons,breaks=c(0,cutoff.exon,1),labels= cols.sign ) )
        sel.exons <- pvals.exons <= cutoff.exon
    } else {
        col.pvals.exons <- rep(cols.sign[2], g2.pval[1])
        sel.exons <- rep(FALSE, g2.pval[1])
    }
    
    pdf(paste(mydir.output.perGene,"/heatmap_exonAndSNPCors_gene_",xi,".pdf",sep=""),width=10,height=6)
    
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


