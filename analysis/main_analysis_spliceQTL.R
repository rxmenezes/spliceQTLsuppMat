# 
# This script preprocesses the data and applies the multinomial test, 
# with W estimated from the data  
# The GEUVADIS data for chromosome 1 is used

mydir <- "~/Documents/NKI/projects/spliceQTL/data and code/analysis"
mydir <- "~/Documents/NKI/projects/spliceQTL/man 2020/BiomJ/codeTest/spliceQTLsuppMat/analysis"

#
mydir.scripts <- file.path(mydir, "scripts")
mydir.data <- file.path(mydir,  "data")
mydir.output <-  file.path(mydir,  "output")
mydir.output.perGene <- mydir.output

myname <- "ObsData_Multin_W"

library(globaltest)
library(gplots)
library(lme4)

source(paste0(mydir.scripts, "/load_objects_spliceQTL.R")) 
# In the script we load objects, and leave out YRI samples
# We also define variables to be used, including the number of permutations, defined as a vector

source(paste0(mydir.scripts, "/make_selection_matrices_strand.R"))  # produce selection matrices  
# Note: here exon.ann and data.exp1 have their rows re-ordered by genomic position

source(paste0(mydir.scripts, "/checks_Nexons_Nsnps.R"))

### Run model, per gene, using all exons and all SNPs
# The following script can be skipped for a test run
#source(paste0(mydir.scripts, "/run_fit_spliceqtl.R"))

# Results presented in the article can be loaded directly:
load(paste0(mydir.output, "/spliceQTL_fit_VariousNperm_Multin_W.RData"))

source(paste0(mydir.scripts, "/merge_results_tables.R"))
# Merging multiple results into a single matrix 
# This script needs to be adapted depending on the results to be merged
# We suggest loading the results res.all as below:

source(paste0(mydir.scripts, "/graphs_spliceQTL_results_W.R")) # here establish threshold and select genes, leave graphs out
## Select some genes to work further

### Run tests per exon within each selected gene - this may take a while to run
source(paste0(mydir.scripts, "/model_fit_eqtl_factors_W.R")) 
# here we extract p-values per exon, for each selected gene 
# per run, exon chosen represents successes, and the sum of all other exons represent failures


# The script below may take a while to run, as it produces graphs per gene
# Each time, it computes the Spearman correlations between all pairs of SNPs and exons,
# which amounts to a large number of correlations

source(paste0(mydir.scripts, "/graphs_spliceQTLfactor_results_W.R")) # includes comparisons with the normal results via the ids of the genes 

source(paste0(mydir.scripts, "/model_fit_eqtl_factors_perExonAndSNP.R")) # only for significant exons, test one SNP at a time and get results
# this leads to a list of pairs of exons and SNPs

# selected with one but not with the other test
# res.all gives the number of exons and the number of SNPs per gene, so compute the total number of tests
# then get bonfcut = alpha/ntests = 0.05/ntests
# Then, per gene that is not selected (1376-34=1344), check how many p-values are below the threshold
# Add these all up
# For genes that are significant: repeat this for the exons that are NOT significant, adding them up

# Results we obtained can be loaded below:
# load(paste0(mydir.output,"/spliceQTL_fit_factorized_multin_W.RData")) # results.aspl.factors
# This object is a result of analysing 1378 genes,with my.nperm=10000 - 56 genes selected

# In the script below graphs per gene can be produced
# Change the gene ids as desired
###
# Below figures 2 and 3 are produced
###
source(paste0(mydir.scripts, "/graphs_spliceQTLfactor_results_W_forPresentation.R"))


