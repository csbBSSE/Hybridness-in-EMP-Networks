library(tidyverse)
library(philentropy)

cleaner <- function(topoFile)
{
    topoDat <- readLines(topoFile) %>% str_squish()
    writeLines(topoDat, topoFile)
}

dirs <-list.dirs(".")[-1]

clean <- sapply(dirs, function(x){
    setwd(x)
    topoFiles <- list.files(".", ".topo$")
    d <- sapply(topoFiles, cleaner)
    setwd("..")
})

jsdCalc <- function(topoFile, WTfreqDf)
{
    net <- topoFile %>% str_remove(".topo")
    freqDf <- read.csv(paste0(net, "_finFlagFreq.csv"))
    freqDf <- freqDf %>% filter(flag == 1) %>%
        select(states, Avg0)
    if(nrow(freqDf) == 0)
        return(NA)
    freqDf <- freqDf[complete.cases(freqDf), ] %>% select(states, Avg0)%>%
        set_names(c("State", "randFreq"))
    df <- merge(WTfreqDf, freqDf, by = "State", all = T)
    # browser()
    df[is.na(df)] <- 0
    JSD(t(df[,-1]))
}

clean <- sapply(dirs, function(x){
    setwd(x)
    topoFiles <- list.files(".", ".topo$")
    WTfreqDf <- read.csv("wild_finFlagFreq.csv") %>%
        filter(flag == 1) %>%
        select(states, Avg0) %>% set_names(c("State", "WTFreq"))
    jsds <- sapply(topoFiles, jsdCalc, WTfreqDf = WTfreqDf)
    df <- data.frame(Net = topoFiles %>% str_remove(".topo"), 
                     JSD = jsds)
    net <- x %>% str_remove("./") %>% str_remove("topo.*")
    setwd("..")
    write.csv(df, paste0(net, "_JSDdat.csv"), row.names = F)
    
})
