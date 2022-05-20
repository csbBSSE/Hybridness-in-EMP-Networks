library(tidyverse)
library(philentropy)
library(Hmisc)
library(corrplot)
source("groupPlotter.R")
source("D:\\Github\\Projects\\Ongoing\\Mubasher\\Figures\\figCore.R")
source("D:\\Github\\Projects\\Ongoing\\Mubasher\\Figures\\utils.R")
options(stringsAsFactors = F, lazy = F)

eNodes <- list(c("GRHL2", "miR200"), c("GRHL2", "miR200", "OVOL2"), 
               c("GRHL2", "OVOL2", "Ecad", "miR34"), 
               c("GRHL2", "miR200", "OVOL2", "miR145"), 
               c("GRHL2", "miR200", "OVOL2", "miR34", "miR145", "Ecad"),
               c("Ecadherin", "miR34"), c("miR200", "miR145"), c("OVOL2", "miR200"))
mNodes <- list(c("ZEB"), c("ZEB"), c("ZEB", "KEAP1", "NRF2", "SNAIL"),
               c("ZEB", "OCT4"), c("ZEB", "KEAP1", "NRF2", "SNAIL", "OCT4"),
               c("ZEB", "KEAP1", "NRF2", "SNAIL"), c("OCT4", "ZEB"), c("ZEB"))

names(eNodes) <- c("grhl", "grhl-ovol", "grhl-ovol-nrf", "grhl-ovol-oct",
                   "grhl-ovol-oct-nrf", "nrf", "oct", "ovol")
names(mNodes) <- c("grhl", "grhl-ovol", "grhl-ovol-nrf", "grhl-ovol-oct",
                   "grhl-ovol-oct-nrf", "nrf", "oct", "ovol")

naZero <- function(df)
{
    df[is.na(df)] <- 0
    df
}

signalNode <- function()
{
    topoFile <- "wild.topo"
    ls <- topo_to_int_mat(topoFile)
    intmat <- ls[[1]]
    nodes <- ls[[2]]
    signal <- nodes[which(apply(intmat, 2, function(x){all(x==0)}))]
    signal
}

signalRemover <- function(state, nodes)
{
    signal <- signalNode()
    r <- which(nodes %in% signal)
    state <- str_split(state, "") %>% unlist
    state <- state[-r] %>% paste0(collapse = "")
    state
}

jsdDat <- function()
{
    jsdFile <- list.files(".", "JsdDat.csv")
    if (is_empty(jsdFile))
    {
        topoFiles <- list.files(".", ".topo")
        wtFile <-  "wild_finFlagFreq.csv"
        nodeFile <- "wild_nodes.txt"
        
        wtDfO <- read_csv(wtFile, col_types = cols()) %>%
            filter(!is.na(Avg0), flag == 1) %>% select(states, Avg0) %>%
            set_names(c("states", "WT"))
        nodes <- readLines(nodeFile) %>% 
            str_replace_all(regex("\\W+"), "")
        
        wtDf <- wtDfO %>% 
            separate(states, c("dummy", "dummy1",all_of(nodes), "dummy2"), sep = "") %>%
            select(-dummy)
        
        JSDs <- sapply(topoFiles, function(x){
            delDf <- read_csv(str_replace(x, ".topo", "_finFlagFreq.csv"), col_types = cols()) %>% 
                filter(!is.na(Avg0), flag == 1) %>% select(states, Avg0)
            nds <- readLines(str_replace(x, ".topo", "_nodes.txt")) %>% 
                str_replace_all(regex("\\W+"), "")
            wt <- wtDfO
            if (!all(nodes %in% nds))
                wt <- wtDf %>% select(dummy1, all_of(nds), dummy2, WT) %>% 
                unite("states", all_of(c("dummy1", nds, "dummy2")), sep = "") %>%
                group_by(states) %>% summarise(WT = sum(WT))
            df <- merge(wt, delDf, by = "states", all = T)
            df[is.na(df)] <- 0
            JSD(df %>% select(-states) %>% t)
        })
        df <- data.frame(Network = str_remove(topoFiles, ".topo"), 
                         JSD = JSDs)
        write.csv(df, "JsdDat.csv", row.names = F)
    }
    else
        df <- read.csv("JsdDat.csv")
    df
    
}

# cohDat <- function()
# {
#     coherencefiles <- list.files(".", "coherence")
#     coherenceDat <- sapply(coherencefiles, function(x){
#         d <- read_csv(x, col_types = cols()) %>% select(coherence) %>% unlist
#         c(min(d, na.rm = T), mean(d, na.rm = T), max(d, na.rm = T))
#     }) %>% t %>% data.frame %>% set_names(c("MinCoh", "MeanCoh", "MaxCoh")) %>%
#         mutate(Network = coherencefiles %>% str_remove("_coherence.csv"))
# }

frustDat <- function()
{
    freqFiles <- list.files(".", "finFlagFreq")
    nets <- freqFiles %>% str_remove("_finFlagFreq.csv")
    frustDat0 <- sapply(freqFiles, function(x){
        d <- read_csv(x, col_types = cols()) %>% filter(flag == 1, !is.na(Avg0)) %>%
            select(frust0) %>% unlist 
        c(min(d, na.rm = T), mean(d, na.rm = T), max(d, na.rm = T))
        
    }) %>% t %>% data.frame %>% set_names(c("MinFrust", "MeanFrust", "MaxFrust")) %>%
        mutate(Network = nets)
    frustDat0
}


JDat <- function()
{
    topoFiles <- list.files(".", ".topo")
    nets <- topoFiles %>% str_remove(".topo")
    topoFile <- "wild.topo"
    signal <- signalNode()
    groupMets <- topoFiles %>% sapply(function(x){
        corFile <- str_replace(x, ".topo", "_cor.csv")
        if(!file.exists(corFile))
            return(NA)
        corDat <- read.csv(corFile, row.names = 1) %>% as.matrix
        JMetric <- corDat[upper.tri(corDat)] %>% abs %>% sum
        JMetric <- JMetric - diag(corDat) %>% sum
        JMetric
    })
    data.frame(Network = nets, Jmetric = groupMets)
}
HybridDat <- function(eNodes, mNodes)
{
    topoFiles <- list.files(".", ".topo")
    nets <- topoFiles %>% str_remove(".topo")
    
    hybridness <- sapply(nets, function(net){
        print(net)
        df <- paste0(net, "_finFlagFreq.csv") %>% read_csv(col_types = cols())
        if ("EMTScore" %in% colnames(df))
        {
            df$phenotype <- ifelse(abs(df$EMTScore) <= 0.5, "Hybrid", "Terminal")
        }
        nodes <- readLines(paste0(net, "_nodes.txt")) %>% 
            str_replace_all(regex("\\W+"), "")
        dfStates <- df %>% select(states) %>% 
            separate(states, into = c("dummy1","dummy2",  nodes, "dummy3"), sep = "")
        df$E <- dfStates %>%
            select(all_of(eNodes)) %>% apply(1, function(x){
            x %>% as.integer %>% sum
        })
        df$M <- dfStates %>%
            select(all_of(mNodes)) %>% apply(1, function(x){
                x %>% as.integer %>% sum
            })
        df$Score <- df$E/length(eNodes) - df$M/length(mNodes)
        df$phenotype <- ifelse(abs(df$Score) <= 0.5, "Hybrid", "Terminal")
        df %>% filter(phenotype == "Hybrid") %>% 
            select(Avg0) %>%
            unlist %>% sum(na.rm = T)
    })
    data.frame(Network = nets, Hybridness = hybridness) %>% naZero
}

ssDatCompileBool <- function(folder)
{
    cwd <- getwd()
    setwd(boolRaw)
    setwd(folder)
    n <- str_remove(folder, "./")
    if(n %in% EMPNets)
    {
        groups <- groupCalc1("wild.topo")[[1]]
        eNode <- groups$E
        mNode <- groups$M
    }
    else
    {
        eNode <- eNodes[[n]]
        mNode <- mNodes[[n]]
    }
    
    dfJSD <- jsdDat()
    dfFrust <- frustDat()
    # dfCoh <- cohDat()
    # dfGS <- groupStrength()
    dfJ <- JDat()
    # dfloops <- loopsDat()
    dfHybrid <- HybridDat(eNode, mNode)
    # dfCons <- consistency(n)
    df <- list(dfJSD, dfFrust,dfJ,dfHybrid) %>% 
        reduce(merge, by = "Network", all = T) %>% mutate(Network = str_replace_all(Network, "-", "_"))
    
    
    setwd(bool)
    nam <- folder %>% str_remove("./") 
    if (!(folder %in% EMPNets))
        nam <- nam %>% str_remove_all("-")
    write.csv(df, paste0(nam, "_all.csv"), row.names = F)
    
    
}
folders <- list.dirs(".", recursive = F)
corDat <- lapply(folders, ssDatCompileBool)
# corMatPlot("ovol")