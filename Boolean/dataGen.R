library(tidyverse)
library(compiler)

### Correlation Matrices ----

statesplit <- function(s){
    s <- s %>% str_split("") %>% unlist %>% as.integer
    s[-c(1, length(s))]
}
statesplit <- Vectorize(statesplit)

cormat <- function(topoFile)
{
    print(topoFile)
    net <- topoFile %>% str_remove(".topo")
    freqDf <- read.csv(paste0(net, "_finFlagFreq.csv")) %>% #filter(flag == 1) %>% 
        select(states, Avg0) %>% drop_na() %>% mutate(Avg0 = Avg0*10000 %>% round)
    nodes <- readLines(paste0(net, "_nodes.txt"))
    states <- rep(freqDf$states, freqDf$Avg0) %>% statesplit %>% t %>% data.frame %>%
        set_names(nodes)
    corr <- cor(states) %>% data.frame
    write.csv(corr, paste0(net, "_cor.csv"))
    return(corr)
}
cormat <- Vectorize(cmpfun(cormat))

correturn <- function(folder)
{
    cwd <- getwd()
    setwd(folder)
    topoFiles <- list.files(".", ".topo")
    corMats <- cormat(topoFiles)
    setwd(cwd)
    return(corMats)
}

folders <- list.dirs(".", recursive = F) %>% str_remove("./")
dummy <- sapply(folders, correturn)

### Loops ----

loopDatGen <- function(folder){
    setwd(folder)
    command <- paste0("python ../signConsistency.py undirectedLoops")
    system(command)
    command <- paste0("python ../Positive_loops.py undirectedLoops")
    system(command)
    setwd("..")
}
dummy <- sapply(EMPNets, loopDatGen)


### Influence Metrics ----

InflLoopCalc <- function(topoFile)
{
    inflFile <- paste0("Influence/", str_replace(topoFile, ".topo", ".csv"))
    inflMat <- read.csv(inflFile)
    
    nodes <- inflMat[[1]]
    rownames(inflMat) <- nodes
    df <- as.matrix(inflMat[,-1])
    row.names(df) <- nodes
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
    
    return(InflLoops)
}

InflLoopGen <- function(folder)
{
    setwd(boolRaw)
    setwd(folder)
    if(!dir.exists("InflLoops"))
        dir.create("InflLoops")
    topoFiles <- list.files(".", ".topo")
    sapply(topoFiles, function(x){
        InflLoopCalc(x) %>% write.csv(paste0("InflLoops/", x %>% str_replace(".topo", "_inflLoops.csv")), 
                                      row.names = F)
    })
    setwd("..")
}
folders <- list.dirs(".", recursive = F) %>% str_remove("./")
sapply(folders, InflLoopGen)
