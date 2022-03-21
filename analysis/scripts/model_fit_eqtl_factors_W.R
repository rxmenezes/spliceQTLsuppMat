# model_fit_eqtl_factors_W.R

# Here we fit a global test per exon of each gene, thus factorizing the results
# This should be done only for a selection of the genes, such as the significant ones
# So, we run the function to factorize the test (run gt per response variable),
# but using only a selection of rows of the mat.gexp.sel and the corresponding selection from mat.snp.sel
# Automatically, only the genes corresponding to the rows of these matrices will be used

#sel.genes <- list.sel
mat.gexp.sel2 <- mat.gexp.sel[list.sel,]
mat.snp.sel2 <- mat.snp.sel[list.sel,]
# Since the 11 genes with empty windows, or NAs in the snp selection matrices, have been left out, we should have no NAs in the mat.snp.sel2

date()
results.aspl.factors <- sapply(1:length(list.sel),splice.eqtl.factor.w,data.exp = data.exp.only, data.snp = mat.snp.only, 
                          exp.sel = mat.gexp.sel2, snp.sel =  mat.snp.sel2,w.type="cov",nperm=10^5)
date()
names(results.aspl.factors) <- list.sel
# This is now a list, as for each gene we have a different number of exons/responses
save(results.aspl.factors, 
     file=paste0(mydir.output,"/spliceQTL_fit_factorized_multin_W.RData",sep="")) # results.aspl.factors, chr 1 only, no 1-exon genes, 1378 genes, my.nperm=10000 - 56 genes


# To extract the data for the selected genes, use for example
# data.exp1.sel <- data.exp1[as.character(data.exp1$Gene_Symbol) %in% sel.genes,]
# data.exp.only.sel <- data.exp.only[as.character(data.exp1$Gene_Symbol) %in% sel.genes,]

