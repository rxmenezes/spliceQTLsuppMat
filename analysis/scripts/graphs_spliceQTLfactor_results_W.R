# graphs_spliceQTLfactor_results_W.R

# list.sel: a vector of gene IDs corresponding to selected genes, created in graphs_spliceQTL_results_W.R
# Extracting the table from the spliceQTL results that corresponds to these genes
aspl.overl <- res.all[list.sel,]
dim(aspl.overl) # 64   3 for normal, 34  3 for multin.W
# new number of tests: number already tested + the total number of exons corresponding to the selected genes
cutoff.exon <- myalpha/( sum(aspl.overl[,1]) + sum(!is.na(res.all[,3])) ) # 2.3 x 10-5 after shift, 2.788622e-05 for multin.W
# Results we obtained can be loaded below:
# load(paste0(mydir.output,"/spliceQTL_fit_factorized_multin_W.RData")) # results.aspl.factors
# results.aspl.factors, chr 1 only, no 1-exon genes, 1378 genes, my.nperm=10000 - 55 genes with no shift, 66 with shift
summary.pvals <- lapply(results.aspl.factors,summary) # for what summary.pvals looks like, see output/summary_pvals_results_withShift_multin_W.txt
mat.snp.sel.exons <- mat.snp.sel[list.sel,]

source(paste0(mydir.scripts, "/graphs_allResults_PerGene_short_multinW.R"))



