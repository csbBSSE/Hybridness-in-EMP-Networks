library(tidyverse)

plasticityCalc <- function(folder)
{
    cwd <- getwd()
    setwd(folder)
    fils <- list.files(".", ".dat")
    nets <- str_remove(fils, "_solution.dat")
    plasticity <- sapply(fils, function(x){
        df <- read.delim(x, header = F)[[2]]
        sum(df > 1)/length(df)
    })
    setwd(cwd)
    df <- data.frame(Network = nets, Plasticity = plasticity)
    write.csv(df, paste0(folder, "_plasticity.csv"))
    return(df)
}