library(philentropy)

statesplit <- function(s){
    s <- s %>% str_split("") %>% unlist %>% as.integer
    s[-c(1, length(s))]
}
statesplit <- Vectorize(statesplit)

correlationMatBool <- function(freqFile, Net)
{#browser()
    net <- str_remove(freqFile, "_finFlagFreq.csv")
    print(net)
    corFile <- paste0("Correlations/", net, "_cor.csv")
    if (!file.exists(corFile)) {
        topoFile <- paste0(net, ".topo")
        # browser()
        nodes <- readLines(paste0(Net, "_nodes.txt"))
        corMat <- read_csv(freqFile, show_col_types = F) %>% 
            select(states, Avg0) %>% drop_na %>%
            mutate(Avg0 = Avg0*10000 %>% round) %>% apply(1, function(x){
                s <- x[1] %>% str_remove_all("'")
                n <- x[2] %>% as.integer
                rep(s, n)
            }) %>% unlist 
        if (length(corMat) < 3)
            return(NA)
        corMat <- corMat %>% lapply(function(x){
            str_split(x, "") %>% unlist %>% as.integer
        }) %>% reduce(rbind.data.frame) %>% 
            set_names(nodes) %>% cor
        colnames(corMat) <- colnames(corMat) %>% str_replace_all(regex("\\W+"), "")
        rownames(corMat) <- rownames(corMat) %>% str_replace_all(regex("\\W+"), "")
        if(!dir.exists("Correlations"))
            dir.create("Correlations")
        # browser()
        write.csv(data.frame(corMat), paste0("Correlations/", net, "_cor.csv"))
        # setwd("..")
    }
    # browser()
    df <- read.csv(corFile, row.names = 1) %>% as.matrix
    return(df)
}
correlationMatBool <- cmpfun(correlationMatBool)
corDat <- function(net) {
    freqFiles <- list.files(".", "_finFlagFreq")
    Jmetric <- sapply(freqFiles, function(x) {
        corDat <- correlationMatBool(x, net)
        JMe<- corDat[upper.tri(corDat)] %>% sum
        JMe <- JMe - diag(corDat) %>% sum
        JMe
    })
    corDf <- data.frame(Network = freqFiles %>% str_remove("_finFlagFreq.csv") %>%
                            str_replace(net, "rand"), 
                        Jmetric = Jmetric)
    corDf
}

JSDProper <- function(net, compute = F)
{
    jsdFile <- paste0(net, "_JSDdata.csv")
    if(!file.exists(jsdFile) || compute)
    {
        freqFiles <- list.files(".", "_finFlagFreq")
        wtFile <- paste0(boolRaw, "/", net, "/wild_finFlagFreq.csv")
        nodeFile <- paste0(boolRaw, "/", net, "/wild_nodes.txt")
        if (net == "silviera")
        {
            wtFile <- paste0(boolRaw, "/SIL/wild_finFlagFreq.csv")
            nodeFile <- paste0(boolRaw, "/SIL/wild_nodes.txt")
        }
        if (net == "silviera2")
        {
            wtFile <- paste0(boolRaw, "/SIL2/wild_finFlagFreq.csv")
            nodeFile <- paste0(boolRaw, "/SIL2/wild_nodes.txt")
        }
        if (net == "EMT_MET")
        {
            wtFile <- paste0(boolRaw, "/EMT_MET_reduced/wild_finFlagFreq.csv")
            nodeFile <- paste0(boolRaw, "/EMT_MET_reduced/wild_nodes.txt")
        }
        
        wtDfO <- read_csv(wtFile) %>%
            filter(!is.na(Avg0), flag == 1) %>% select(states, Avg0) %>%
            set_names(c("states", "WT"))
        
        
        JSDs <- sapply(freqFiles, function(x){
            delDf <- read_csv(x) %>% 
                filter(!is.na(Avg0), flag == 1) %>% select(states, Avg0)
            wt <- wtDfO
            df <- merge(wt, delDf, by = "states", all = T)
            df[is.na(df)] <- 0
            JSD(df %>% select(-states) %>% t)
        })
        df <- data.frame(Network = str_remove(freqFiles, "_finFlagFreq.csv") %>%
                             str_replace(net, "rand"), 
                         JSD = JSDs)
        write.csv(df, paste0(net, "_JSDdata.csv"), row.names = F)
    }
    df <- read.csv(jsdFile)
    return(df)
}
