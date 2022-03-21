# make_partition_20_equalSetSizes.R

# Make a partition of the data creating 5 subsets of ~70 samples, and check which SNPs have to be left out because all their obs are zero.
# Use only the subset of SNPs that is present in all subsets. 
# Then select only that SNP subset from the large dataset

# Object expected: snp.data.l, a list with length equal to the number of genes
# One element in snp.data.l is a data.matrix containing the SNP data (number of minor alleles)
# corresponding to SNPs mapping to a gene, with SNPs on rows and samples on columns

set.seed <- 543678
npart <- 5      # the data set is divided into 5 approx equally-sized subsets
nshuffles <- 20 # the number of times the data set is subdivided - so the total number of subsets is npart*nshuffles
part.mat <- matrix(0,nrow=ncol(snp.data.l[[1]]),ncol=nshuffles)
nsam <- floor( ncol(snp.data.l[[1]])/npart ) # number of samples per subset - here 373/5 gives approximately 74
for(xj in 1:ncol(part.mat))
  {
   partition <- c(rep(1:npart, nsam), sample(1:npart, 3)) # This has 373 entries, all {1,2,3,4,5}
   partition <- sample(partition) # This creates a permutation of the vector above
   part.mat[, xj] <- partition
  }
# Now check if there are SNPs constant in either one of the partitions, per gene
# Eliminate those from the main data - there is no need to save the individual partitioned data matrices

rows.out <- vector("list", length(snp.data.l))
for(xi in 1:length(snp.data.l))
{
  all.data <- snp.data.l[[xi]]
  data.ann <- snp.ann.l[[xi]]
  for(xj in 1:ncol(part.mat))
  {  
    partition <- part.mat[, xj]
    for(xk in 1:npart)
   {
     dp <- diag(as.numeric(partition==xk))
     dp <- dp[,colSums(dp)>0]
     Xp <-  all.data %*% dp
     var.0 <- apply(Xp,1,var)==0
     if(sum(var.0)>0) { rows.out[[xi]] <- c(rows.out[[xi]],(1:nrow(Xp))[var.0] )  }
   }
   print(xj)
   print(xi)
   print(length(rows.out[[xi]]))
    }
  if(length(rows.out[[xi]])>0) { 
    all.data <- all.data[-rows.out[[xi]],]
    data.ann <- data.ann[-rows.out[[xi]],]
    snp.data.l[[xi]] <- all.data
    snp.ann.l[[xi]] <- data.ann
  }
}
# In this case, for gene 2 only 2 SNPs had to be left out.
#  So we shorten all gene sets to get them to have equal size
rows.out
lapply(snp.data.l, nrow) # 3286 for the first and third matrices, but 3284 for the second
new.nrow <- min(unlist(lapply(snp.data.l,nrow)))
leave.out <- (1:length(snp.data.l))[ unlist(lapply(snp.data.l,nrow)) == new.nrow]
for(xi in (1:length(snp.data.l))[-leave.out] ) snp.data.l[[xi]] <- snp.data.l[[xi]][-rows.out[[leave.out]],]
lapply(snp.data.l, nrow) # 3284 for all 3 matrices

# If something had changed, we would update the objects with the two rows below
save(snp.data.l,file="/snp_data_equalSetSizes.RData") # snp.data.l
save(snp.ann.l,file="/snp_ann_equalSetSizes.RData") # snp.ann.l


# Result: one matrix X per gene, with M SNPs/rows and n samples/columns
# I have made this into a  single list, with length = number of genes
# In addition, a vector "partition" defines the partitions: each entry  is one of {1,2,3,4,5}
# This vector is then used to partition using diag: 
#   first compute dp <- diag(as.numeric(partition==xi))
#   Then eliminate columns of dp that contain all zeroes: dp <- dp[,colSums(dp)>0]; check: ncol(dp) should be data size
#   Finally, compute the selection: Xp <-  X %*% dp, for each partition
# Note: the same partition is used in all cases

