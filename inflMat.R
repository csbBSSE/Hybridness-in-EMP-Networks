library(tidyverse)
options(stringsAsFactors = F)
library(compiler)

topo_to_int_mat <- function(topo_file) {
    # print(topo_file)
    df <- read.delim(topo_file, sep = " ", stringsAsFactors = F)
    if (ncol(df) != 3) {
        df <- read.delim(topo_file, stringsAsFactors = F)
    }
    # browser()
    colnames(df) <- c("Source", "Target", "Type")
    df <- df %>% 
        mutate(Source = str_remove_all(Source, "\\s")) %>%
        mutate(Target = str_remove_all(Target, "\\s")) %>%
        mutate(Type = ifelse(Type == 2, -1, 1))
    
    nodes <- unique(c(df$Source, df$Target)) %>% 
        sort(decreasing = T)
    n_nodes <- length(nodes)
    intmat <- rep(0, n_nodes * n_nodes) %>% 
        matrix(ncol = n_nodes)
    df1 <- df %>% 
        mutate(Source = sapply(Source, function(x) {which(nodes == x)})) %>% 
        mutate(Target = sapply(Target, function(x) {which(nodes == x)}))
    # browser()
    dummy <- apply(df1, 1, function(x) {
        # browser()
        i <- x[1]
        j <- x[2]
        k <- x[3]
        intmat[i,j] <<- k
    })
    return(list(intmat, nodes))
}

compute_power_matrix <- function(mat, power) {
    res <- mat
    if (power == 1)
    {
        return(res)
    }
    for (i in 2:power) {
        res <- res %*% mat
    }
    return(res)
}

influence_matrix <- function(topoFile, lmax = 10) {
    ls <- topo_to_int_mat(topoFile)
    intmat <- ls[[1]]
    intmax <- intmat
    intmax[which(intmax == -1)] <- 1
    res <- 0
    for (l in 1:lmax) {
        intM <- compute_power_matrix(intmat, l)
        maxM <- compute_power_matrix(intmax, l)
        r1 <- intM / maxM
        r1[is.nan(r1)] <- intM[is.nan(r1)]
        res <- res + r1
    }
    res <- res / lmax
    
    nodes <- ls[[2]] %>% str_replace_all(regex("\\W+"), "")
    
    influence_mat <- res
    colnames(influence_mat) <- rownames(influence_mat) <- nodes
    signal <- which(apply(intmat, 2, function(x){all(x==0)}))
    output <- which(apply(intmat, 1, function(x){all(x==0)}))
    secondary_signal <- which(apply(intmat, 2, function(x){
        if (length(signal) !=0)
            all(x[-signal] == 0)
        else
            F
    }))
    secondary_output <- which(apply(intmat, 1, function(x){
        if (length(output) != 0)
            all(x[-output] == 0)
        else
            F
    }))
    nonEssentials <- c(signal, output, secondary_signal, secondary_output)
    if(length(nonEssentials))
    {
        influence_reduced <- influence_mat[-nonEssentials, 
                                           -nonEssentials]
        nodes_reduced <- nodes[-nonEssentials] %>% c  %>% str_replace_all(regex("\\W+"), "")
    }
    else
    {
        influence_reduced <- influence_mat
        nodes_reduced <- nodes %>% str_replace_all(regex("\\W+"), "")
    }
    if (length(nodes_reduced) < 2)
        return()
    rownames(influence_reduced) <- colnames(influence_reduced) <- nodes_reduced
    
    
    influence_reduced <- influence_reduced[order(influence_reduced[,1]), order(influence_reduced[1,])]
    net <- str_remove(topoFile, ".topo") %>% paste0("Influence/",.)
    if (!dir.exists("Influence"))
        dir.create("Influence")
    write.csv(influence_reduced, paste0(net, ".csv"))
    write.csv(influence_mat, paste0(net, "_fullMat.csv"))
}
influence_matrix <- cmpfun(influence_matrix)

sigOut <- function(topoFile)
{
    ls <- topo_to_int_mat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]] %>% str_replace_all(regex("\\W+"), "")
    signal <- which(apply(intmat, 2, function(x){all(x==0)}))
    output <- which(apply(intmat, 1, function(x){all(x==0)}))
    secondary_signal <- which(apply(intmat, 2, function(x){
        if (length(signal) !=0)
            all(x[-signal] == 0)
        else
            F
    }))
    secondary_output <- which(apply(intmat, 1, function(x){
        if (length(output) != 0)
            all(x[-output] == 0)
        else
            F
    }))
    return(list(nodes[c(signal, secondary_signal)], nodes[c(output, secondary_output)]))
}
