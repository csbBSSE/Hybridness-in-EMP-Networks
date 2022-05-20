library(tidyverse)

getToggles <- function(df)
{
    df <- df %>% filter(Type == 2) 
}

hiLoopData <- function(fol)
{
    setwd(fol)
    net <- fol
    topoDat <- 
    setwd("..")
}