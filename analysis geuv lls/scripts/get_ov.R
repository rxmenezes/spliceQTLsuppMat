## Function getov
## On the basis of a gene id, selects items from each of 2 lists with that gene id as name,
## and computes the number of sub-items with value smaller than a threshold
##
## Inputs
## id: names of items within res1, all equal to names of items in res2
## res1, res2: two lists with the same object types, same length, same names - here p-values
## th: the p-value threshold cut-off, can also be a vector
## Output
## list with two elements per id, for genes with the same number of pairs for both datasets:
## cor: correlation between -log(p-values)
## ov.th: overlap between results, i.e. average number of pairs for which pair was selected in both cases 
## Genes without the same number of pairs get NA
getov.cor.old <- function(id, res1, res2)
{
    r1 <- res1[[ id ]]
    r2 <- res2[[ id ]]
    # Check if objects have same length/dimensions
    if( length(r1) == length(r2) ){ if((!is.null(dim(r1))) & (!is.null(dim(r2))) ){
        if( (nrow(r1) == nrow(r2)) & (ncol(r1) == ncol(r2))){ 
            ov.cor <- cor(-log(matrix(r1, nrow=nrow(r1)*ncol(r1), ncol = 1)[, 1]),
                          -log(matrix(r2, nrow=nrow(r2)*ncol(r2), ncol = 1)[, 1]))
        }
    }} else { ov.cor <- NA }
    ov.cor
}

getov.cor <- function(id, res1, res2)
{
    r1 <- res1[[ id ]]
    r2 <- res2[[ id ]]
    # Check if objects have same length/dimensions
    if( length(r1) == length(r2) ){ if((!is.null(dim(r1))) & (!is.null(dim(r2))) ){
        if( (nrow(r1) == nrow(r2)) & (ncol(r1) == ncol(r2))){ 
            ov.cor <- cor(-log(matrix(r1, nrow=nrow(r1)*ncol(r1), ncol = 1)[, 1]),
                          -log(matrix(r2, nrow=nrow(r2)*ncol(r2), ncol = 1)[, 1]))
        } else { ov.cor <- NA }
    }} else { ov.cor <- NA }
    ov.cor
}


getov.th <- function(id, res1, res2, th)
{
    r1 <- res1[[ id ]]
    r2 <- res2[[ id ]]
    if( length(r1) == length(r2) ){ if((!is.null(dim(r1))) & (!is.null(dim(r2))) ){
        if( (nrow(r1) == nrow(r2)) & (ncol(r1) == ncol(r2))){ 
            ov.v <- NULL
            for(myt in seq_along(th))
            {
                r1s <- r1 <= th[ myt ]
                r2s <- r2 <= th[ myt ]
                ov.v <- c(ov.v, sum(r1s == r2s)/length(r1s))
            }
        } else { 
            ov.v <- rep(NA, length(th))
        }
    }} else {
        ov.v <- rep(NA, length(th))
    }
    names(ov.v) <- th
    ov.v
}

# Function similar to the one above, but now we compute only the average number
# of selections per dataset and threshold, so this returns results for all genes
getov2 <- function(id, res1, res2, th)
{
    r1 <- res1[[ id ]]
    r2 <- res2[[ id ]]
    ov.v <- NULL
    for(myt in seq_along(th))
    {
        r1s <- r1 <= th[ myt ]
        r2s <- r2 <= th[ myt ]
        ov.v <- rbind(ov.v, c(sum(r1s)/length(r1s), sum(r2s)/length(r2s)))
    }
    rownames(ov.v) <- th
    colnames(ov.v) <- c("res1", "res2")
    ov.v
}
