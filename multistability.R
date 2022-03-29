
#this code has network solution folder as input and output is plasticity score as excel file.

@author: Mubasher

library(tidyverse)
multistabilityCalc <- function(folder)
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

#type athe following command in the console and hit enter
multistabilityCalc("sol_files")   #folder_name is the name of the folder containing RACIPE solution files
