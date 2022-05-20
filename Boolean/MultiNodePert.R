library(tidyverse)
library(stringr)
library(future)
library(future.apply)
library(data.table)
library(compiler)
source("D:\\Teams\\TeamsPaperFigures\\figCore.R")
source("D:\\Teams\\Final_Results\\codes\\inflMat.R")
source("D:\\Teams\\Final_Results\\codes\\stateLabeller.R")


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
hamming <- function(x)
{
    x <- str_split(x, "")
    sum(x[[1]]!=x[[2]])/(length(x[[1]])-2)
} 
hamming <- cmpfun(hamming)

dec2bin <- function(x, l=32) 
{
    paste(rev(rev(as.integer(rev(intToBits(x))))[1:l]), collapse = "")
} %>% cmpfun


state_gen <- function(n){
    s_list <- 1:(2^n)
    sapply(s_list, dec2bin, n)
} %>% cmpfun

inwards <- function(x,df){
    colnames(df) <- c("State", "Step")
    df %>% filter(Step == x) %>% select(State) %>% unlist
} %>% cmpfun

asyncInit <- function(sInit, update_matrix, maxT = 1000, init = NULL)
{#browser()
    f <- 0 #flag
    N_nodes <- length(sInit)
    up_ind_list <- sample(1:N_nodes, maxT, replace = T) # generate the update indices 
    s <- sInit
    for (k in 1:maxT){
        s_dummy <- s%*%update_matrix %>% sign # update the state
        if (all(s_dummy == s)) 
            f <- 1 # flag 1 means it is a proper steady state
        
        up_ind <- up_ind_list[k]
        s[up_ind] <- s_dummy[up_ind]
        
        if (f == 1) 
            break
    }
    if (is.null(init))
        init <- paste0("'", 
                       paste0(
                           replace(sInit, which(s==-1),0), 
                           collapse = ""),
                       "'")
    fin <- paste0("'", 
                  paste0(
                      replace(s, which(s==-1),0), 
                      collapse = ""),
                  "'")
    return(c(init, fin, f))
}

asyncInit <- cmpfun(asyncInit)

coherenceIter <- function(sInit, update_matrix, maxT = 1000, nNode = 1,
                          randChoice = F, nState = 100, nIter = 10)
{#browser()
    init <- sInit
    sInit <- str_split(sInit, "") %>% unlist
    sInit <- as.integer(sInit[2:(length(sInit)-1)])
    sInit <- ifelse(sInit == 0, -1, 1)
    if (nNode>1)
        randChoice <- T
    nNodes <- length(sInit)
    pertChoice <- sapply(1:nState, function(x){
        sample(1:nNodes, nNode)
    })
    if (nNode > 1)
        pertChoice <- t(pertChoice)
    else
        pertChoice <- matrix(pertChoice, ncol = 1)
    
    stateList <- apply(pertChoice, 1, function(x){
        s <- sInit
        s[x] <- -1*s[x]
        s
    })
    
    if(!randChoice)
        stateList <- sapply(1:nNodes, function(x){
            s <- sInit
            s[x] <- -1*s[x]
            s
        })
    # browser()
    df <- apply(stateList,2, function(x){
        sapply(1:nIter, function(i){
            asyncInit(x, update_matrix, maxT, init)
        }) %>% t %>% data.frame %>% set_names(c("init", "fin", "flag"))
    }, simplify = F) %>% reduce(rbind.data.frame) %>% 
        group_by(init, fin, flag) %>%
        summarise(Freq = n(), .groups = 'drop') %>% ungroup
    df$Freq <- round(df$Freq/sum(df$Freq),2)
    
    df
    
}

coherenceIter <- cmpfun(coherenceIter)

coherence <- function(net, maxT = 1000, nNode = 1,
                      randChoice = F, nState = 100, nIter = 10, sMax = 25, write=T)
{
    dataDir <- ""
    topoFile <- paste0(dataDir, net, ".topo")
    nodeOrder <- readLines(paste0(dataDir, net, "_nodes.txt")) %>% str_replace_all(regex("\\W+"), "")
    ls <- topo_to_int_mat(topoFile)
    update_matrix <- ls[[1]]
    nodes <- ls[[2]] %>% str_replace_all(regex("\\W+"), "")
    colnames(update_matrix) <- rownames(update_matrix) <- nodes
    update_matrix <- update_matrix[nodeOrder, nodeOrder]
    update_matrix <- 2*update_matrix + diag(length(nodes))
    
    states <- read_csv(paste0(dataDir, net, "_finFlagFreq.csv"), col_types = cols(), lazy = F) %>%
        filter(!is.na(Avg0), flag == 1)
    if (nrow(states) == 0)
        return(NA)
    if (nrow(states) > 2*sMax)
    {
        states <- states %>% arrange(frust0) %>% slice(c(1:sMax, (nrow(.)-(sMax-1)):nrow(.)))
    }
    states <- states$states
    # browser()
    df <- lapply(states, function(x){
        coherenceIter(x, update_matrix, maxT,nNode,
                      randChoice, nState, nIter)
    }) %>% reduce(rbind.data.frame) 
    df$flag <- as.integer(df$flag)
    df$nNode <- nNode
    df$Hamming <- df %>% select(init, fin) %>% apply(1, hamming)
    df$initPhen <- df$init %>% labeller_states(topoFile = paste0(net, ".topo"))
    df$finPhen <- df$fin %>% labeller_states(topoFile = paste0(net, ".topo"))
    if (write)
        write_csv(df, paste0(net, "_coherence.csv"), 
                  quote = "none")
    else
        return(df)
}

coherence <- cmpfun(coherence)


coherenceAllNode <- function(net, maxT = 1000, randChoice = T, 
                             nState = 100, nIter = 10, sMax = 25,
                             dataDir = NULL, reps = 1)
{
    dataDir <- paste0(randRaw, "/", net, "/")
    if (str_detect(net, "rand"))
        dataDir <- paste0(randRaw, "/", str_remove(net, "_rand.*"), "/")
    setwd(dataDir)
    l <- readLines(paste0(dataDir, net, "_nodes.txt")) %>% length
    # browser()
    dfList <- lapply(1:reps, function(rep){
        plan(multisession, workers = 6)
        df <- future_lapply(1:l, function(x){
            coherence(net, maxT, x, randChoice, nState, nIter, sMax, write = F)
        }) %>% reduce(rbind.data.frame) %>% mutate(Rep = rep)
        future:::ClusterRegistry("stop")
        return(df)
    })
    if (reps == 1)
        df <- dfList[[1]]
    else
        df <- dfList %>% reduce(rbind.data.frame)
    write_csv(df, paste0(WTcoherence, "/", net, "_allNodeCoherence_nPert", 
                         nState,"_nIter",nIter,"_reps", reps, ".csv"), quote = "none")
}
coherenceAllNode <- cmpfun(coherenceAllNode)

coherenceSingleNode <- function(nCores = 6)
{
    # setwd(paste0(randRaw, "/", net))
    topoFiles <- list.files(".", ".topo")
    nets <- str_remove(topoFiles, ".topo")
    coherences <- list.files(".", "coherence")
    netsDone <- coherences %>% str_remove("_coherence.csv")
    nets <- nets[!(nets %in% netsDone)]
    # dir.create("dummy")
    # setwd("dummy")
    
    plan(multisession, workers = nCores)
    dummy <- future_lapply(nets, function(x){
        # browser()
        print(x)
        df <- coherence(x, randChoice = F, nIter = 50, write = F)
        if(is.na(df))
            return()
        df <- df[df$init == df$fin,]%>% mutate(states = init, coherence = Freq) %>%
            # ungroup %>% 
            select(states, coherence)
        write_csv(df, paste0(x, "_coherence.csv"), quote = "none")
        
    })
}
coherenceSingleNode <- cmpfun(coherenceSingleNode)

folders <- c("grhl", "ovol", "oct", "nrf")
dummy <- sapply(folders, function(x){
    setwd(x)
    coherenceSingleNode()
    setwd("..")
})