# make_suppFigure3.R
# This plot was originally made in main_analysis.Rmd in folder "scripts vm new"

set.seed(54321)
nsnps.sel <- 1000
snpsel <- sample(1:nrow(mat.snp.only), nsnps.sel)
scors <- cor(t(mat.snp.only[snpsel, ]), method = "s")
#sum(abs(scors) > 0.5)

# Consecutive snps

snpsel <- sample(1:(nrow(mat.snp.only) - nsnps.sel), 1)
scors2 <- cor(t(mat.snp.only[snpsel:(snpsel+nsnps.sel), ]), method = "s")

pdf(paste0(mydir.output, "/suppFigure3.pdf"), width=8, height=6)
plot(density(scors2), type = "l", col = "blue", lwd = 2, xlab = "Spearman correlation between consecutive SNPs", ylab = "density", main = "")
dev.off()
