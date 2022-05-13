# make_suppTable4.R
# Internal reference: this comes from scripts server/main_table_results_2tests.R

### Supplementary table 4
# This table is made by cross-tabulating results of the spliceQTL test 
# and of the linear regression tests summarised per gene, using the GEUVADIS dataset.
###

load(paste(mydir.output,"/spliceQTL_fit_VariousNperm_asSingleMatrix_W_withShift.RData",sep="")) # res.all, a matrix with 3 columns and one row per gene
myalpha <- 0.05
my.cut <- myalpha/( nrow(res.all) - sum(is.na(res.all[, 3])) )  # 3.633721e-05 for multin_W, very similar to the previous 3.628447e-05
sum( res.all[, 3] <= my.cut ,na.rm=T) # 34 when SNPs with MAF<0.05 were left out, 
res.nonas <- res.all[ !(is.na(res.all[, 3])), ]
list.sel <- rownames( res.nonas[ res.nonas[, 3] <= my.cut , ] )
list.genes.notsel <- rownames( res.nonas[ res.nonas[, 3] > my.cut , ] )
# length(list.genes.notsel) # 1342

# Results G2_binomial_W, as a list, per exon, only for the significant genes
# results.aspl.factors, chr 1 only, no 1-exon genes, my.nperm=10000 
aspl.overl <- res.nonas[list.sel, ]
aspl.notsel <- res.nonas[list.genes.notsel, ]
# dim(aspl.overl) # 34  3 
# new number of tests: number already tested + the total number of exons corresponding to the selected genes
ntests.g2.exons <- ( sum(aspl.overl[, 1]) + nrow(res.nonas) ) # 1796
cutoff.exon <- myalpha/ntests.g2.exons # 2.783964e-05 new, 2.788622e-05 if ntests.g2.exons = 1793

# Results global test per pair of (significant) exon and SNP
# results.aspl.snps, chr 1 only, no 1-exon genes, my.nperm=10000 - 34 genes
names(results.aspl.snps) <- list.sel
ntests.g2.pairs <-  crossprod( aspl.overl[, 1], aspl.overl[, 2])[1, 1] # gives 1625539, before it was 1617287
cutoff.pairs <- myalpha/( ntests.g2.exons + ntests.g2.pairs   )  # gives 3.072508e-08, before it was 3.088174e-08

# Define Bonferroni significance cut-off to be used
# For this, we need the total number of tests (linreg)
# Note that the lists of genes, sel and not sel, do not include genes with NAS already

ntests.pairs <- crossprod(res.nonas[,1] , res.nonas[,2])[1, 1] # gives 53286710
cut.bonf <- myalpha/ntests.pairs   # gives 9.383203e-10

# Counting number of NS tests with spliceQTL test

table.linreg.butG2notSign <- matrix(0, nrow=length(list.genes.notsel), ncol=2)
rownames(table.linreg.butG2notSign ) <- list.genes.notsel
colnames(table.linreg.butG2notSign) <- c("notSign", "Sign")
what = "lr"
mydir.output.cov <- paste0(mydir.data, "/covMatrices_perGene")

# Load linear regression results for tests found to not be significant
# These results were produced by running the script 'model_fit_linreg.R', provided
# Here we provide the output directly for efficiency, and read these in below

for(xj in list.genes.notsel)
{  
    load(paste0(mydir.output.cov, "/pvalues_mat_", xj, "_", what, ".RData")) # pvals.mat
    table.linreg.butG2notSign[xj, ] <- c(sum( pvals.mat > cut.bonf ), sum( pvals.mat <= cut.bonf ))
}   

# Now count genes sign with linreg but not significant with GT at the exon level
# table.linreg.butG2BinNotSign contains sums of not.sign and sign testsfor linreg, 
# that were definitely notsign for G2Bin
table.linreg.butG2BinNotSign <- matrix(0, nrow=length(list.sel), ncol=2)
rownames(table.linreg.butG2BinNotSign ) <- list.sel
colnames(table.linreg.butG2BinNotSign) <- c("notSign", "Sign")
tables.sum <- matrix(0, nrow=2, ncol=2)
rownames(tables.sum) <- c(paste(what, "NS"), paste(what, "S"))
colnames(tables.sum) <- c("gt NS", "gt S")

for(xj in list.sel)
{
    res.exons <- results.aspl.factors[[ xj ]]
    load(paste0(mydir.output.cov, "/pvalues_mat_", xj, "_", what, ".RData")) # pvals.mat
    pvals.mat.exonNsel <- pvals.mat[res.exons > cutoff.exon, , drop=FALSE]
    table.linreg.butG2BinNotSign[xj, ] <- c(sum( pvals.mat.exonNsel > cut.bonf ), 
                                            sum( pvals.mat.exonNsel <= cut.bonf )) 
    # number of pairs found sign by linreg, for exons not sign
    if(sum( (res.exons <= cutoff.exon) )  > 0) {
        pvals.mat.exonSel <- pvals.mat[res.exons <= cutoff.exon, , drop=FALSE] # here nrow=nexons, ncol=nSNPs
        res.esnps <- results.aspl.snps[[ xj ]]  # here nrow=nSNPs, ncol=nexons
        for(xk in 1:sum( res.exons <= cutoff.exon ) ) 
        {  
            sign.linreg <- pvals.mat.exonSel[xk, ] <= cut.bonf
            sign.gt <- res.esnps[, xk] <= cutoff.pairs
            mytab <- matrix(c( sum((!sign.linreg) & (!sign.gt)), 
                               sum((sign.linreg) & (!sign.gt)), 
                               sum((!sign.linreg) & (sign.gt)), 
                               sum((sign.linreg) & (sign.gt)) ), nrow=2, ncol=2)
            tables.sum <- tables.sum + mytab
        }
    }
}

# For genes that are G2.NS - table.linreg.butG2notSign
# For genes that are G2.S but have no exons G2Bin sign - table.linreg.butG2BinNotSign
# Now join the two tables
table.genesNS.OrExonsNS <- rbind(table.linreg.butG2notSign, table.linreg.butG2BinNotSign)
# So to the table we have to add, to the first column (gt.NS), the values in this row
# colSums(table.genesNS.OrExonsNS)
#  notSign     Sign 
# 53021050    12790

# So we add this to tables.sum[,1]
tables.sum[,1] <- tables.sum[,1] + colSums(table.genesNS.OrExonsNS)
tables.sum
# gt NS gt S
# lr NS 53257624 3452
# lr S     15737 9897
