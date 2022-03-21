# model_fit_eqtl_factors_perExonAndSNP.R

# This script is called by main_analysis_Multin_W.R -- last modified in May 2015, but in 2016 I added this comment ONLY

# Now we perform the global test per exon, for the significant exons, extracting p-values for each SNP
# This yields a matrix per gene, which might have a single row of only one exon is significant
# Then results for all genes considered (the ones that were significant) are gathered into a list
# It is possible that none of the exons is significant, in which case we do not go further for that gene.

# Computing the cut-off to be used
aspl.overl <- res.all[list.sel,]
dim(aspl.overl) # 64   3 for normal, 34  3 for multin.W
# new number of tests: number already tested + the total number of exons corresponding to the selected genes
cutoff.exon <- myalpha/( sum(aspl.overl[,1]) + sum(!is.na(res.all[,3])) ) # 2.3 x 10-5 after shift, 2.788622e-05 for multin.W

# We first select rows of selection matrices for those genes that were significant in the first run (multivariate multinomial with W)
# We also use as input the list with results per exon calculated with the multinomial test with W=cov(X),
# and the cut-off to be used per exon

mat.gexp.sel2 <- mat.gexp.sel[list.sel,]
mat.snp.sel2 <- mat.snp.sel[list.sel,]

results.aspl.snps <- sapply(1:length(list.sel), splice.eqtl.pairs,
                            data.exp = data.exp.only, data.snp = mat.snp.only, 
                            exp.sel = mat.gexp.sel2, snp.sel =  mat.snp.sel2,
                            snp.ann=data.snp.ann,
                            results=results.aspl.factors, what="p", cutoff = cutoff.exon)
  
names(results.aspl.snps) <- list.sel

save(results.aspl.snps,file=paste0(mydir.output,"/spliceQTL_fit_pvalues_perPairExonSNP.RData")) # results.aspl.snps, chr 1 only, no 1-exon genes, 1378 genes, my.nperm=10000 - 56 genes
# Note that this is a list with one element per gene - for each gene this is a matrix with SNPS on rows and exons on columns
# Ids for the exons are colnames
# There are no ids for the SNPs - these we need to get from the annotation data snp.data.ann, for that gene - get mat.snp.sel2


