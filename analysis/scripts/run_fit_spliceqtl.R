# run_fit_eqtl_fixedNperm.R

# This is based on the script in file "model_fit_eqtl.R", not available via the manuscript

# Here we compute the multinomial test per gene
# We leave out features with too large p-values after a number of permutations
# This is not very efficient yet, as the permutations are repeated from the beginning,
# The function runs with sapply(), so processing can be made faster by using
# sfSapply() from package snowfall, which makes use of parallel computing

# Here we test for eqtl per gene
# Objects: 
#     data.exp.only, mat.snp.only - data matrices
#     mat.gexp.sel, mat.snp.sel - to select rows from data.exp1, mat.snp.only
# Output: results.eqtl, with annotation columns for genes

# v.nperm contains the numbers of permutations each time
# v.cutoff contains the cut-offs for p-values each time - vector with the same length as v.nperm
# v.genes.used = list with the same length as v.nperm, containing the list of genes used each time - is updated
# v.results = list with the same length as v.nperm, with the results matrix computed each time
# The strategy is to let the function select the feature each time, using the geneID
# Once the row corresponding to the feature is determined, so is the row of the selection matrices too
# So, there is no need to make smaller versions of the data
# For this, we make a new splice.eqtl.fit2 function, where the list of genes used is also an input
# This is then used to select the rows of the selection matrices to be used (as many as there are genes, in the same order)
#

data.exp1$Gene_Symbol <- factor(as.character(data.exp1$Gene_Symbol))  
genes.used <- levels(data.exp1$Gene_Symbol)
v.results <- vector("list",length(v.nperm))
# Use w.type = "cov" as below to compute the multinomial test with W estimated from the data
# Use W as the identity matrix with as many rows as there are SNPs in the data (per gene)
# to compute the test with W=I
# Use multin = FALSE to use the test based on the normal distribution

for(xj in 1:length(v.nperm))
{
    my.nperm <- v.nperm[ xj ]
    results.aspl <- sapply(1:length(genes.used),splice.eqtl.fit2,genes.used = genes.used, 
                           data.exp = data.exp.only, data.snp = mat.snp.only, 
                           exp.sel = mat.gexp.sel, snp.sel =  mat.snp.sel,
                           nperm=my.nperm, multin = TRUE, w.type = "cov")
    results.aspl <- t(results.aspl)
    colnames(results.aspl) <- c("N.exons","N.SNPs","p.value") # add the extra columns, one per combination of alpha and mu values
    rownames(results.aspl) <- genes.used
    v.results[[xj]] <- results.aspl
    names(v.results)[xj] <- my.nperm
    save(results.aspl,file=paste0(mydir.output,"/spliceQTL_fit_Nperm",my.nperm,".RData"))
    # Results with NAs are not of interest for further calculations, so we make sure they are left out
    results.aspl <- results.aspl[ !is.na(results.aspl[,3]) , ]
    genes.used <- rownames(results.aspl)[ results.aspl[,3] <= v.cutoff[xj] ]
}

# This object is generated with 3 rows (N.exons, N.SNPs and p.value) and as many columns as levels in f.gene
# So transpose to assign colnames to it

save(v.results, file=paste0(mydir.output,"/spliceQTL_fit_VariousNperm_Multin_W.RData")) # v.results, chr 1 only, no 1-exon genes
