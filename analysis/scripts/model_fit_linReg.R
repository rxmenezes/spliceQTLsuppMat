# model_fit_linReg.R


data.exp1$Gene_Symbol <- factor(as.character(data.exp1$Gene_Symbol))  
genes.used <- levels(data.exp1$Gene_Symbol)
mydir.output.cov <- paste(mydir.output,"/covMatrices_perGene",sep="")

# tests to check function splice.eqtl.fit2 in script tests_g2multin_rhovec_muvec.R
results.aspl <- sapply(1:length(genes.used),splice.eqtl.lreg,genes.used = genes.used, data.exp = data.exp.only, data.snp = mat.snp.only, 
                         exp.sel = mat.gexp.sel, snp.sel =  mat.snp.sel,what="lr",mydir=mydir.output.cov)

# Note that splice.eqtl.lreg standardizes the expression of each exon per gene, whereby
# the model is robust to eQTL effects

# In case it is not possible to update v.results, use
#source("spliceQTL_update_v_results.R")

# For nperm=100, 1389; for nperm=1000, 511; for n=10^5, 349

# This object is generated with 3 rows (N.exons, N.SNPs and p.value) and as many columns as levels in f.gene
# So transpose to assign colnames to it

#save(v.results,file=paste(mydir.output,"/spliceQTL_fit_VariousNperm.RData",sep="")) # chr 1 only, no 1-exon genes
#save(v.results,file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_withShift.RData",sep="")) # v.results, chr 1 only, no 1-exon genes
save(v.results,file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_W_withShift.RData",sep="")) # v.results, chr 1 only, no 1-exon genes

#table(is.na(results.aspl[,2]),is.na(results.aspl[,3]))
#        FALSE TRUE
#  FALSE  1378    0
#  TRUE      0   11
# So for 11 genes we have NAs because these do not have any SNPs in their windows

# source("selection_signResults_spliceQTL_old.R") # this was used in model_fit_eqtl.R, but it no longer suits this modified script
# MTC and selection of significant results are now done in graphs_spliceQTL_results.R



