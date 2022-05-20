library(tidyverse)
options(stringsAsFactors = F)
library(compiler)

topo_to_int_mat <- function(df, weights) {
    # print(topo_file)
    
    df <- df %>% 
        mutate(Type = weights)
    
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

influence_matrix <- function(df, weights, lmax = 10) {
    ls <- topo_to_int_mat(df, weights)
    intmat <- ls[[1]]
    intmax <- intmat
    intmax[which(intmax < 0)] <- -1*intmax[which(intmax < 0)]
    res <- 0
    for (l in 1:lmax) {
        intM <- compute_power_matrix(intmat, l)
        maxM <- compute_power_matrix(intmax, l)
        r1 <- intM / maxM
        r1[is.nan(r1)] <- intM[is.nan(r1)]
        res <- res + r1
    }
    res <- res / lmax
    
    nodes <- ls[[2]]
    
    influence_mat <- res
    colnames(influence_mat) <- rownames(influence_mat) <- nodes
    return(influence_mat)
}
influence_matrix <- cmpfun(influence_matrix)


InflLoopCalc <- function(df, weights, Enodes, Mnodes, lmax = 10)
{
    inflMat <- influence_matrix(df, weights, lmax)
    nodes <- c(Enodes, Mnodes)
    # browser()
    inflMat <- inflMat[nodes, nodes]
    df <- inflMat
    InflLoops <- lapply(1:length(nodes), function(i){
        sapply(i:length(nodes), function(j){
            n1 <- nodes[i]
            n2 <- nodes[j]
            mult <- df[n1, n2]*df[n2,n1]
            loop <- abs(df[n1,n2]) + abs(df[n2,n1])
            if(mult >0)
                if (df[n1,n2] > 0)
                    n <- 1
            else
                n <- 2
            else
                n <- 3
            c(ifelse(mult < 0, -1, 1), loop,n, i,j)
        }) %>% t %>% data.frame %>% set_names(c("Nature", "Strength", "Type","Node1", "Node2"))
    }) %>% reduce(rbind.data.frame) %>% mutate(Nature = ifelse(Nature == -1, "N", "P")) %>%
        mutate(Node1 = nodes[Node1], Node2 = nodes[Node2]) %>% unite("Node", Node1, Node2, sep = ",")
    typeCast <- c("DA", "TS", "NF")
    InflLoops$Type <- typeCast[InflLoops$Type]
    inflTS <- InflLoops %>% filter(Type == "TS") %>% select(Strength) %>% unlist %>% sum
    inflDA <- InflLoops %>% filter(Type == "DA") %>% select(Strength) %>% unlist %>% sum
    inflnegLoops <- InflLoops %>% filter(Type == "NF") %>% select(Strength) %>% unlist %>% sum
    inflposLoops <- inflTS + inflDA
    return(c(inflnegLoops, inflposLoops))
}
