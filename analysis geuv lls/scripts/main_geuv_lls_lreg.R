###
# Script to analyse two datasets using tests: per exon-SNP pair (here a linear regression)
# and per gene, for all exons and SNPs at once (here the spliceQTL test)

source("lreg.R") # provided

# Per chr and dataset, the code yields a list of matrices, one per gene
for(chr in 10:22)
{ 
    for(data in c("Geuvadis", "LLS"))
    { 
        cat("Analysing", data, chr, ":", as.character(Sys.time()), "\n")
        #   rm(list = setdiff(ls(), c("data", "chr", "path"))); gc(); cat(".")
        load(file.path(path, paste0("temp.", data, ".chr", chr, ".RData"))); cat(".") # NOT provided
        # the object loaded above contains map (a data.frame), exons (a list) and snps (a list)
        myres <- lreg(Y = exons, X = snps, map = map, data = data)
        save(myres, file = paste( paste0("lregResults_", data, "_chr", chr, ".RData") ))
        rm(snps)
        rm(exons)
        rm(map)
        gc()
        pvalue <- test.multiple(Y = exons, X = snps, map = map, spe = 16, w.type = "cov"); cat(".")
        save(object=pvalue, file=file.path(opath, paste0("pvalNew.", data, ".chr", chr, ".RData"))); cat("\n")
    }
}

