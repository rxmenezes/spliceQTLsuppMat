# merge_results_tables.R

res.all <- v.results[[length(v.nperm)]]
for(xi in (length(v.nperm)-1):1)
{
 res.next <- v.results[[xi]] 
 res.next <- res.next[!(rownames(res.next) %in% rownames(res.all)),]
 res.all <- rbind(res.all,res.next)
}

mat.gexp.sel <- cbind(mat.gexp.sel,1:nrow(mat.gexp.sel))
mat.gexp.sel <- mat.gexp.sel[ order(rownames(mat.gexp.sel)) ,]
res.all <- res.all[ order(rownames(res.all)) ,]
res.all <- res.all[ order(mat.gexp.sel[,3]), ]
mat.gexp.sel <- mat.gexp.sel[order(mat.gexp.sel[,3]),]
# Checking that the two objects have rows in the same order
sum(rownames(mat.gexp.sel) == rownames(res.all)) # 1389, so OK!


# mat.gexp.sel was used to re-order rows in res.all, results from spliceQTL tests
# Checking that this order is the same as in mat.snp.sel, for consistency
sum(rownames(mat.gexp.sel) == rownames(mat.snp.sel)) # 1389, so OK!
# In case this script is followed by others, fix mat.gexp.sel
mat.gexp.sel <- mat.gexp.sel[,-3]
# There is no need to save mat.gexp.sel or mat.snp.sel since they are not altered after all

#save(res.all, file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_asSingleMatrix.RData",sep="")) # res.all, chr 1 only, no 1-exon genes
#write.table(res.all, file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_asSingleMatrix.txt",sep=""),sep="\t",col.names=NA) # chr 1 only, no 1-exon genes

# Results test stat with Normal dist
#save(res.all, file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_asSingleMatrix_withShift.RData",sep="")) # res.all, chr 1 only, no 1-exon genes
#write.table(res.all, file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_asSingleMatrix_withShift.txt",sep=""),sep="\t",col.names=NA) # chr 1 only, no 1-exon genes

# Results test stat with multinomial dist, W=cov(X)
save(res.all, file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_asSingleMatrix_W_withShift.RData",sep="")) # res.all, chr 1 only, no 1-exon genes
write.table(res.all, file=paste(mydir.output,"/spliceQTL_fit_VariousNperm_asSingleMatrix_W_withShift.txt",sep=""),sep="\t",col.names=NA) # chr 1 only, no 1-exon genes


