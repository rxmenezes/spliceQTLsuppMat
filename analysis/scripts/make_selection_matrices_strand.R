# make_selection_matrices_strand.R

# This script then creates selection matrices, which are integer matrices:
# -- for the mRNA data, this is  matrix mat.gexp.sel 
# involving as many rows as there are genes in the data, and each row contains
# the first and last indices of rows in data.exp.only corresponding to exons
# mapping to that gene
# -- for the SNP data, this is the snp.rows.sel matrix, with as many rows as genes in the data,
# and each row contains the first and last indices of rows in mat.snp.only
# of SNPs mapping to that gene
# 

# Also save a matrix with gene ids, as well as start and end positions for genes

# Window: 1 Mb to each side of the start, where start = min(exon_start) if strand == "+"
#                                          and  start = max(exon_end) if strand == "-"

# Create a matrix with one row per gene, and one column for each of: n.exons, start, end, geneloc, gt-pvalue
#  where we here define n.exons as the number of exons found in the data for that gene,
#  start = min(start(all exons for that gene)), end = max(end(all exons for that gene)), geneloc = gene location as in "coord"column
genes.used <- levels(data.exp1$Gene_Symbol)
results.altspl <- matrix(0,nrow=length(genes.used),ncol=5)
rownames(results.altspl) <- genes.used
colnames(results.altspl) <- c("N.exons","StartEx","EndEx","Coord","StartGene")
# Creating matrices that will contain start and end rows corresponding to specific genes
mat.gexp.sel <- mat.snp.sel <- matrix(0,nrow=length(genes.used),ncol=2)
rownames(mat.gexp.sel) <- rownames(mat.snp.sel) <- genes.used
mystrand1 <- NULL
# Note that "N.exons" represents the number of exons in the file,
# "StartEx" is the start position of the first exon, "EndEx" is the end position of the last exon,
# and "Coord" is the gene coordinate given in the original gene annotation (downloaded with the data)
for(xi in 1:length(genes.used))
{
  sel.gexp <- (1:nrow(data.exp1))[ as.character(data.exp1$Gene_Symbol)==genes.used[xi] ]
  mystrand <- as.character(geneid.strand$strand)[geneid.strand$GeneID == genes.used[xi] ]
  mystrand1 <- c(mystrand1, mystrand)
  StartEx <- min(exon.ann1$exon.start[sel.gexp])
  EndEx <-  max(exon.ann1$exon.stop[sel.gexp])
  StartGene <- StartEx  # we keep this value if mystrand == "+" or if the strand is not known
  if(length(mystrand) > 0)  if( mystrand == "-"){ StartGene <- EndEx } # modified if negative strand
  N.exons <- length(sel.gexp)
  Coord <- data.exp1$Coord[sel.gexp][1]
  # Save data per gene as row in results.altspl
  results.altspl[xi,] <- c(N.exons,StartEx,EndEx,Coord,StartGene)
  # Also save rows in data matrix corresponding to exons for each gene - this is as before, in make_selection_matrices.R
  mat.gexp.sel[xi,] <- c(min(sel.gexp),max(sel.gexp))
  # Now for the SNP data, we find rows mapping to the window
  # See make_selection_matrices.R for selection of SNPs mapping to the gene only
  snp.rows.sel <- ( data.snp.ann$POS >= (StartGene-10^6) )& ( data.snp.ann$POS <= (StartGene+10^6) )
  if(any(snp.rows.sel)) { mat.snp.sel[xi,] <- range( (1:nrow(data.snp.ann))[snp.rows.sel] ) } else {
                          mat.snp.sel[xi,] <- c(NA,NA) }
}

results.altspl <- data.frame(results.altspl,Strand=mystrand1) # check what was done to get the strand info right, here I have only 
#   1389 genes that are in geneid.strand and in the data set
#load(paste(mydir.output,"/Gene_nexons_nsnps_chr1_all.RData",sep="")) # results.altspl, produced with script checks_Nexons_Nsnps.R


save(mat.gexp.sel, file=paste0(mydir.data,"/mat_gexp_sel_all.RData")) # chr 1 only, no 1-exon genes
save(mat.snp.sel, file=paste0(mydir.data,"/mat_snp_sel_all.RData"))   # chr 1 only, no 1-exon genes

# There were some NAs found

sum(is.na(mat.snp.sel)) # 26 in total after recoding
table(apply(is.na(mat.snp.sel),1,sum))
#   0    2 
#1376   13 
# Before re-coding we had only 11 here, now 13: the 2 extra must have been left out before - so they were probably recoded 
# AND they major allele must have a freq larger than 0.95
# So NAs are in both columns
nas.snp.sel <- rownames(mat.snp.sel)[is.na(mat.snp.sel[,2])]
nas.snp.sel
#nas.snp.sel
# [1] "ENSG00000135747.7"  "ENSG00000153207.10" "ENSG00000162711.12"
# [4] "ENSG00000162714.7"  "ENSG00000171161.7"  "ENSG00000171163.10"
# [7] "ENSG00000175137.9"  "ENSG00000185220.7"  "ENSG00000188295.9" 
#[10] "ENSG00000196418.6"  "ENSG00000197472.9"  "ENSG00000232274.1" 
#[13] "ENSG00000232336.1" # The last two are extra, compared to the previous round - they probably are located
# in a region where all existing SNPs have the newly calculated MAF <= 0.05

data.exp.out <- data.exp1[as.character(data.exp1$Gene_Symbol) %in% nas.snp.sel,]
range(data.exp.out$Coord)
# 143188828 249200395
range(data.exp.out$Coord)/10^6 # (in Mb)
# 143.1888 249.2004
names(table(data.exp.out$Coord/10^6))
# [1] "143.188828" "143.21107"  "247.09528"  "247.171395" "247.242113"
# [6] "247.267674" "247.33531"  "247.495148" "247.579458" "249.120832"
#[11] "249.132409" "249.153343" "249.200395"

range(data.snp.ann$POS)
#[1]     57952 245999953
# So most of these genes (11 of them) start downstream from the last SNP
# The two others have the two last IDs in the list

# So the 11 genes left out involve a window that starts downstream from the last SNP

# What about the other two?
data.exp.out[as.character(data.exp.out$Gene_Symbol) %in% nas.snp.sel[12:13],1:3]
#ENSG00000232274.1_143183876_143184895 ENSG00000232274.1   1 143211070
#ENSG00000232274.1_143188795_143189487 ENSG00000232274.1   1 143211070
#ENSG00000232274.1_143209430_143211070 ENSG00000232274.1   1 143211070
#ENSG00000232336.1_143188828_143189024 ENSG00000232336.1   1 143188828
#ENSG00000232336.1_143189423_143189695 ENSG00000232336.1   1 143188828
sum( (data.snp.ann$POS >= (143211070-10^6)) & (data.snp.ann$POS <= (143211070+10^6)) ) # 0
sum( (data.snp.ann$POS >= (143188828-10^6)) & (data.snp.ann$POS <= (143188828+10^6)) ) # 0

# So indeed, 11 genes are upstream from the last SNP, and two have no SNPs around in their window 
# It is possible these SNPs existed before but were left out in this new run





