#This code has two files - excel file containing steady states (SS's) of the network and topology file of the network - as input and it calculates 
#frustration of each SS as output. Frustration of  a SS is calculated by taking the ratio of total number of frustrated edges in a SS to the total
#number of edges in that SS. Network frustration is then obtained by taking the maximum/minimum/mean of frustrations of all the SS's.

@author: Kishore


library(readxl)

frustCalc <- function(stateFile, topoFile)
{
    dfStates <- read_excel(stateFile, col_names = F)  # reading excel file containing network states (generated from MATLAB)
    dfTopo <- read.delim2(topoFile, sep = "")         # reading topo file of the network
    nodes <- unique(c(dfTopo$Source, dfTopo$Target))
    colnames(dfStates) <- nodes
    dfTopo$Type <- ifelse(dfTopo$Type == 2, -1, 1)
    frustration <- apply(dfStates, 1, function(s){
        s <- as.integer(s)
        names(s) <- nodes
        f <- apply(dfTopo,1,function(x){
        ifelse(s[x[1]]*s[x[2]]*as.integer(x[3]) == -1, 1, 0)
    })
        
        f <- sum(unlist(f))/nrow(dfTopo)
    })
    dfStates$frustration <- frustration
    write.csv(dfStates, paste0(topoFile, "_frustration.csv"), row.names = F)
}

#automatically reading excel and topo files and calculating .csv frustration files by calling frustCal function

exl_files <- list.files(pattern = "xlsx")
topo_files <- list.files(pattern = "topo")

for (i in 1:length(topo_files))
{
    topoFile <- topo_files[i]
    excelFil <- exl_files[i]
    frustCalc(excelFil, topoFile)
}


#read frust data files (.csv format) and calculate mean/min/max
df1 <- list.files(pattern = "csv") 
frust <- vector()
for (i in df1) {
    ldf <- read.csv(i)
    C <- as.data.frame(na.omit(ldf))
    frust[which(df1==i)] <- mean(C$frustration) #c(min(C$frustration), max(C$frustration),mean(C$frustration))
    print(frust[which(df1==i)])
}

names(frust) <- gsub(".topo_frustration.csv","",df1)
output <- as.data.frame(frust)
output
output <- cbind(" "=rownames(output), output)
write_xlsx(output, "frust_mean.xlsx")


