# graphs_spliceQTL_results_W.R

# Selection: FWER <= 0.05 = alpha
# Total n.tests = nrow(res.all)
myalpha <- 0.05
# Computing the FWER threshold
my.cut <- myalpha/( nrow(res.all)-sum(is.na(res.all[,3])) )  # 3.633721e-05 for multin_W, very similar to the previous 3.628447e-05
sum( res.all[,3] <= my.cut ,na.rm=T) # 56 when all SNPs were used, 55 when SNPs with MAF<0.05 were left out, 
#                                      66 with the shift before recoding, 64 after recoding
#                                      34 with W in the multin test

res.nonas <- res.all[ !(is.na(res.all[,3])), ]
list.sel <- rownames( res.nonas[ res.nonas[,3] <= my.cut ,] )
length(list.sel) # 34, so OK!
nsel <- length(list.sel)




