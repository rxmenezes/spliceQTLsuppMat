# load_objects_spliceQTL.R

load(paste0(mydir.data, "/workspace_loadedObjects_spliceQTL.RData"))

# Contained in the above workspace are functions which are stored in:

#source(paste0(mydir.scripts, "/g2_renee.R"))
#source(paste0(mydir.scripts, "/f_spliceQTLfit.R"))
#source(paste0(mydir.scripts, "/functions_wilcoxon_t_logistic_fdr.R"))

# Settings for numbers of permutations and significance cut-off are defined

#v.nperm <- c(10^2,10^3,10^5)
#v.cutoff <- c(0.2,0.1,0.05)

# Also included are the following objects:

# data.exp1: data.frame with exon level data, including annotation columns
# exon.ann1: exon data annotation, rows in the same order as data.exp1, 5 columns
# data.ann: exon data annotation, rows in the same order as data.exp1, 3 columns
# data.exp.only: data matrix containing the (exon-length normalized) exon-level data
# geneid.strand: annotation for chr 1, including exon strand
# mat.snp.only: data matrix with the SNP data (number of minor alleles per SNP),
#   where the SNPs correspond to rows and samples correspond to columns
#   Columns are assumed to be in the same order as those in data.exp.only
# data.snp.ann: data.frame with annotation for SNP data, with as many rows (and in the same
#   order) as in mat.snp.only

# Notes:

# exon-level data
# The exon-level data is normalized by exon length, per gene.
# This involved fitting a mixed-effects model to all exon observations per gene with 
# the exon having a random effect and the length having a fixed effect, after transforming
# the data using a hyperbolic arc-sine. 
# Residuals of this model fit are then used as the data. 
# As this may occasionally generates negative values, the minimum observed value plus epsilon
# was added to all observations to make them non-negative
# Genes with a single exon were left out of the data, as it makes no sense 
# to talk about alternative splicing in such cases

# Loading the exon expression data object data.exp1 
# load(paste0(mydir.data,"/data_exp1_no1exonGenes_chr1_vst_sorted_withAnn_NoYRIsamples_exonSorted_allPos.RData"))  # data.exp1

# The annotation for chr 1, including exon strand, can also be loaded using
# load(paste0(mydir.data,"/gencode_v12_chr1_GeneID_Strand.RData")) # geneid.strand

# data.ann is obtained using
data.ann <- data.exp1[, c(1:3)]

# data.exp.only is obtained using
data.exp.only <- data.exp1[, -c(1:3)]

# SNP data
# Here all observations are equal to either 0, 1 or 2, representing the number of minor alleles
# SNPs with MAF < .05, or with zero variance, were left out of the data

load(paste0(mydir.data,"/data_snp.RData")) # data.snp.ann, mat.snp.only

# The SNP data and annotation can also be loaded using:
# load(paste0(mydir.data,"/mat_snp_only_No0VarSNPs_MAFfrom5pc.RData")) # mat.snp.only
# load(paste0(mydir.data,"/mat_snp_ann_No0VarSNPs_MAFfrom5pc.RData")) # data.snp.ann

# Checking that columns are in the same order
# sum(colnames(mat.snp.only) == colnames(data.exp.only)) # 373, so OK!











