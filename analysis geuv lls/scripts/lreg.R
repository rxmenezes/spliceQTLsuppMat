# lreg.R

lreg.single <- function(Y, X, table.ann, i, map = map, data = data){
    if(!is.null(table.ann)) {  gid <- rownames(table.ann)[i] 
                              } else { 
                               gid <- as.character(map$genes$gene_id)[i] }
    id.map <- which(as.character(map$genes$gene_id) == gid)
    # extract data
    ys <- map$exons[[id.map]]
    y <- Y[, ys, drop=FALSE]
    sdy <- apply(y, 2, sd)
    y <- y[, sdy > 0, drop = FALSE]
    xs <- seq(from=map$snps$from[id.map], to=map$snps$to[id.map], by=1)
    x <- X[, xs, drop=FALSE]
    sdx <- apply(x, 2, sd)
    x <- x[, sdx > 0, drop = FALSE]
    rm(Y,X); silent <- gc()
    myp <- NULL
    if( (ncol(y) > 0) & (ncol(x) > 0) ){
    for(xr in 1:ncol(y))
        for(xc in 1:ncol(x))
        {
            myp <- c(myp, summary(lm(y[, xr] ~ x[, xc]))$coef[2,4])
        }
    mat.p <- matrix(myp, nrow = ncol(y), ncol = ncol(x), byrow = TRUE)
    } else {
        mat.p <- NA
    }
    mat.p
}

# To compute the linear regression (and extract test result) per pair of response and covariate, 
# where table.ann contains a selection table with as rownames ids of genes to be used
# If table.ann is not supplied, this will apply the test to all genes in the data
lreg <- function(Y, X, map, table.ann = NULL, data = NULL){
    if(!is.null(table.ann)) { p <- nrow(table.ann)
    } else {
        p <- nrow(map$genes)
    }
    hm <- lapply(seq_len(p),
                 FUN=function(i) lreg.single(Y=Y, X=X, map=map, i=i, 
                                             table.ann = table.ann, data = data)
    )
    names(hm) <- as.character(map$genes$gene_id)
    hm
}


hmap.s1 <- function(Y, X, map, i, table.ann = table.ann, data = data){
    
    # extract data
    ys <- map$exons[[i]]
    y <- Y[, ys, drop=FALSE]
    xs <- seq(from=map$snps$from[i], to=map$snps$to[i], by=1)
    x <- X[, xs, drop=FALSE]
    rm(Y,X); silent <- gc()
    list(y= y, x = x)
}

hm1 <- function(Y, X, map, table.ann = NULL, data = NULL){
    
    p <- nrow(map$genes)
    hm <- lapply(X=seq_len(p),
                 FUN=function(i) hmap.s1(Y=Y, X=X, map=map, i=i, 
                                         table.ann = table.ann, data = data)
    )
    save(hm, file = "hm_test.RData")
}

hm2 <- function(Y, X, map, table.ann = NULL, data = NULL){
    
    g <- map$genes
    save(g, file = "ges_test.RData") # check if this has the ensembl IDs that I have in tb
}
