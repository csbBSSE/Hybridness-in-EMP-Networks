### This file contains the code to generate raw data required to calculate network topology
## based metrics used in the manuscript. Please make sure to position the 'inflmat.R',
## 'signconsistency.py' and 'Positive_loops.py' scripts in the right folders before 
## running this script. 


library(tidyverse)
library(compiler)
source("inflMat.R") 

### Correlation Matrices ----

statesplit <- function(s){
    s <- s %>% str_split("") %>% unlist %>% as.integer
    s[-c(1, length(s))]
}
statesplit <- Vectorize(statesplit)


### Loops ----
# This function calculates directed an undirected loops for all topo files in the given folder
# It requires the python scripts to be in the parent directory of the current folder.
loopDatGen <- function(folder){
    setwd(folder)
    command <- paste0("python ../signConsistency.py undirectedLoops")
    system(command)
    command <- paste0("python ../Positive_loops.py directedLoops")
    system(command)
    setwd("..")
}
dummy <- sapply(folders, loopDatGen)


### Influence Metrics ----
# This function calculates the influence-based loops for a given topo file. 
# This requires the loading of inflMat.R script, available on the repo.

InflLoopCalc <- function(topoFile)
{
    inflFile <- paste0("Influence/", str_replace(topoFile, ".topo", ".csv"))
    if (!file.exists(inflFile)) {
        influence_matrix(topoFile)
    }
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

# This function calculates influence based loops for all topo files in the folder. 
InflLoopGen <- function(folder)
{
    # setwd(boolRaw)
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
