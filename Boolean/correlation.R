library(tidyverse)
library(compiler)

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

