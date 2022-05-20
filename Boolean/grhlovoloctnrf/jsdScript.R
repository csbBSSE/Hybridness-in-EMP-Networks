library(tidyverse)
library(philentropy)

wtDat <- read.csv("wild_finFlagFreq.csv") %>% filter(flag == 1) %>%
    select(states, Avg0) %>% drop_na
wtNodes <- readLines("wild_nodes.txt")

stateReorder <- function(state, wtnodes, nodes)
{
    if (all(state == nodes))
        return(state)
    if (length(state) == length(nodes))
    {
        state <- str_split(state, "") %>% unlist
        names(state) <- nodes
        state <- paste0(state[wtnodes], collapes = "") %>% unname
        return(state)
    }
    else
    {
        state <- str_split(state, "") %>% unlist
        names(state) <- nodes
        state <- state[wtnodes]
        names(state) <- wtnodes
        state[is.na(state)] <- "0"
        state <- paste0(state, collapse = "") %>% unname
        return(state)
    }
}

jsdCalc <- function(topoFile, wtDat)
{
    nodes <- readLines(str_replace(topoFile, ".topo", "_nodes.txt"))
    df <- read.csv(str_replace(topoFile, ".topo","_finFlagFreq.csv")) %>% 
        filter(flag == 1) %>% select(states, Avg0) %>% drop_na %>% 
        set_names(c("states", "Net")) %>% 
        mutate(states = sapply(states, stateReorder, wtnodes = wtNodes, nodes = nodes)) %>%
        merge(wtDat, by = "states", all = T)
    df[is.na(df)] <- 0
    df %>% select(-states) %>% t %>% JSD
}
topoFiles <- list.files(".", ".topo")
jsds <- sapply(topoFiles, jsdCalc, wtDat = wtDat)
jsdDf <- data.frame(Net = topoFiles %>% str_remove(".topo"), JSD = jsds)
write.csv(jsdDf, "JSD.csv", row.names = F)
