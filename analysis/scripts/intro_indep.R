# intro_indep.R


source(paste0(mydir.scripts, "/g2_renee.R"))
source(paste0(mydir.scripts, "/f_spliceQTLfit.R"))
source(paste0(mydir.scripts, "/functions_wilcoxon_t_logistic_fdr.R"))

v.nperm <- c(10^2,10^3,10^5)
v.cutoff <- c(0.2,0.1,0.05)




