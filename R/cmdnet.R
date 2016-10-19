#======================================================================================================
# cmdnet: an algorithm for identification of preserved networks in multiple datasets
# author: Nicholas Wisniewski
# date: October 17, 2016
#======================================================================================================



library(psych)
library(matrixcalc)


norm_vec <- function(x) sqrt(sum(x^2))


#' Correlation matrix distance
#' 
#' Computes the correlation matrix distance between two correlation matrices.
#' @param r1 The first correlation matrix
#' @param r2 The second correlation matrix
#' @param type whether to compute statistics at the network or node level
#' @param metric whether to use the angular distance or the cosine distance
#' @return A number or vector of distances, depending on the type of computation specified.
#' @export
cmdist <- function(r1, r2, type=c("network","node"), metric=c("angle", "cosine")){
    # sum(r1 * r2) is a faster way to compute tr(r1 %*% r2) for symmetric matrices
    if(type=="network" & metric=="cosine") return( 1 - sum(r1 * r2)/(norm(r1,"F")*norm(r2,"F")) ) 
    if(type=="network" & metric=="angle") return(2/pi * acos( round (sum(r1 * r2)/(norm(r1,"F")*norm(r2,"F")) ,13)   )  )
    if(type=="node" & metric=="cosine") return( 1 - diag(r1 %*% r2)/(apply(r1, 1, norm_vec)*apply(r2, 1, norm_vec))  )  
    if(type=="node" & metric=="angle") return( 1/pi * acos( round( diag(r1 %*% r2)/(apply(r1, 1, norm_vec)*apply(r2, 1, norm_vec)), 13) ) ) 
} 



# find genes common across all datasets, and align the rows
#' Clean the list of datasets
#' 
#' Takes a list of datasets, finds the nodes that are common across all, and matches the rows so they are in the same format.
#' @param datalist A list of datasets
#' @param geneset An optional list of row names to reduce the datasets by
#' @return A list of datasets that have matching rows
#' @export
cleanDatalist <- function(datalist, geneset=NULL){
    if(is.null(geneset)) geneset <- Reduce(intersect, mclapply(datalist, rownames))
    geneset_matches <- as.data.frame(lapply(datalist, FUN=function(x) match(geneset, rownames(x)) ))
    geneset_common <- na.omit(geneset_matches)  # remove genes that don't exist across all 
    datalist_common <- list()
    datalist_common <- vector("list",length(datalist))
    for (i in 1:length(datalist)) datalist_common[[i]] <- datalist[[i]][geneset_common[,i],]
    names(datalist_common) <- names(datalist)
    datalist_common
}


# compute full statistics
#' Correlation matrix distance across networks
#' 
#' Computes correlation matrix distance statistics for a list of datasets. 
#' @param datalist A list of datasets
#' @param genelist An optional list of row names to reduce the datasets by
#' @param mode Which of 3 computations to perform. Defaults to all 3.
#' @param corFUN A function specifying which correlation measure to use. Defaults to "cor"
#' @param xform Whether to transform corrlations using Fisher's z transform. Defaults to "none".
#' @param r0 An optional consensus matrix to compute deviations around. Defaults to NULL, and computes average matrix.
#' @param metric Whether to use angular distance or cosine distance. Defaults to "angle".
#' @return A list of correlation matrix distance statistics
#' @export
cmdnet <- function( datalist, genelist=NULL, mode=c("cmd","node","net"), corFUN=cor, xform=c("none", "fisher"), 
                    r0=NULL, metric=c("angle", "cosine")) {
    data.list <- cleanDatalist(datalist, genelist)
    groups <- names(data.list)
    ngroups <- length(data.list)
    nu <- unlist(lapply(data.list, FUN=function(x) ncol(x)))
    n <- sum(nu)  # total sample size
    
    r <- rw <- vector("list", ngroups)
    for (i in 1:ngroups) {  # loop over datasets and compute correlations
        r[[i]] <- corFUN(t(data.list[[i]]))
        if(is.null(r0) & xform[1]=="none") rw[[i]] <- nu[i] * r[[i]] # weighted by sample size
        if(is.null(r0) & xform[1]=="fisher") rw[[i]] <- nu[i] * fisherz(r[[i]]) # weighted by sample size and variance stabilized
    }
    if(is.null(r0)){
        r0 <- Reduce("+", rw) / n  # weighted mean correlation matrix
        if(xform[1]=="fisher") { r00 <- apply(r0, c(1:2), fisherz2r); diag(r00) <-1; r0 <- r00 }
    }
    
    # initialize some variables
    dev.matrix <- numeric()
    dev.vector <- matrix()
    dev.vector.cor1 <- matrix()
    dev.vector.cor2 <- matrix()
    meandev.matrix <- numeric()
    meandev.vector <- numeric()
    net <- matrix()
    
    # which parts of the anlysis are we gonna run?
    if("cmd" %in% mode){
        # Compute deviations about the mean correlation matrix
        dev.matrix <- sapply(r, function(ri) cmdist(r0, ri, type="network", metric=metric[1]))  # deviations
        names(dev.matrix) <- as.character(groups)
        meandev.matrix <- weighted.mean(dev.matrix, nu, na.rm=T)  # mean deviations
    }
    if("node" %in% mode){
        # Compute deviations for each row of the correlation matrix about its mean correlation vector
        dev.vector <- sapply(r, function(ri) cmdist(r0, ri, type="node", metric=metric[1])) # vector of deviations    
        colnames(dev.vector) <- as.character(groups)
        dev.vector.cor1 <- cor(dev.vector)
        dev.vector.cor2 <- cor(t(dev.vector))
        meandev.vector <- apply(dev.vector, 1, function(dev) weighted.mean(dev, nu, na.rm=T))
    }
    if("net" %in% mode){
        # Compute group network (similarity matrix)
        net <- matrix(NA, ngroups, ngroups)
        for (i in 1:ngroups) for(j in 1:ngroups) net[i, j] <- 1 - cmdist(r[[i]], r[[j]], type="network", metric=metric[1])
        rownames(net) <- colnames(net) <- groups
    }
    
    # Return
    list(cmd.dev = dev.matrix,
         cmd.meandev = meandev.matrix, 
         node.dev = dev.vector,  
         node.meandev = meandev.vector, 
         node.dev.cor1 = dev.vector.cor1,
         node.dev.cor2 = dev.vector.cor2,
         net = net,
         group.names = groups,
         group.r = r, 
         group.n = nu, 
         total.n = n,
         r0 = r0
    )
}


