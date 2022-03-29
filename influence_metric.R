options(stringsAsFactors = F)
library(compiler)
library(stringr)


@author: Kishore

# This code has input as "topo_file_name.topo" & lmax =10. Its output is Influence Matrix in .csv format. This matrix is needed to calculate Influence values 

topo_to_int_mat <- function(topo_file) {
    # print(topo_file)
    df <- read.delim(topo_file, sep = "", stringsAsFactors = F)  # sep = ""/" " for 2/1 tab spaces between the columns in the topo file
    if (ncol(df) != 3) {
        df <- read.delim(topo_file, stringsAsFactors = F)
    }
    # browser()
    colnames(df) <- c("Source", "Target", "Type")
    df$Source = str_remove_all(df$Source, "\\s")
    df$Target = str_remove_all(df$Target, "\\s")
    df$Type = ifelse(df$Type == 2, -1, 1)
    
    nodes <- sort(unique(c(df$Source, df$Target)), decreasing = T) #decreasing order of first letter in names
    n_nodes <- length(nodes)
    intmat <- matrix(rep(0, n_nodes * n_nodes), ncol = n_nodes)
    df1 <- df #%>% 
    df1$Source = sapply(df1$Source, function(x) {which(nodes == x)})
    df1$Target = sapply(df1$Target, function(x) {which(nodes == x)})
    # browser()
    dummy <- apply(df1, 1, function(x) {
        # browser()
        i <- x[1]
        j <- x[2]
        k <- x[3]
        intmat[i,j] <<- k
    })
    return(list(intmat, nodes))
}

compute_power_matrix <- function(mat, power) {
    res <- mat
    if (power == 1)
    {
        return(res)
    }
    for (i in 2:power) {
        res <- res %*% mat
    }
    return(res)
}

influence_matrix <- function(topoFile, lmax = 10) {
    ls <- topo_to_int_mat(topoFile)     # Calling "topo_to_int_mat" function
    intmat <- ls[[1]]
    intmax <- intmat
    intmax[which(intmax == -1)] <- 1
    res <- 0
    for (l in 1:lmax) {
        intM <- compute_power_matrix(intmat, l)  # Calling "compute_power_matrix" function
        maxM <- compute_power_matrix(intmax, l)
        r1 <- intM / maxM
        r1[is.nan(r1)] <- intM[is.nan(r1)]
        res <- res + r1
    }
    res <- res / lmax
    
    nodes <- ls[[2]]
    
    influence_mat <- res
    colnames(influence_mat) <- rownames(influence_mat) <- nodes
    
    net <- str_remove(topoFile, ".topo")
    write.csv(influence_mat, paste0(net, "_fullMat.csv")) #writes INf. Matrix in CSV file
}

influence_matrix <- cmpfun(influence_matrix)


#reading all topo files at once and calculating influence matrices
listt <- list.files(pattern = "topo")
for (i in listt) {
    influence_matrix(i,lmax = 10)
}

#calculating influence values
#This code has input as Influence Matrix as CSV file and output as Influence VAlues

df1 <- list.files(pattern = "csv") #loading list of Infl Matrix files
Inf_Val <- vector()
for (i in df1) {
    ldf <- read.csv(i)
    C <- as.data.frame(ldf)
    C$X <- NULL
    C[C > 0] <- 0    # making positive values zeros
    CC <- abs(C) 
    Inf_Val[which(df1==i)] <- sum(CC)   # summing up absolute values of negative values only
    print(Inf_Val[which(df1==i)])
}
names(Inf_Val) <- gsub("_fullMat.csv","",df1)
output <- as.data.frame(Inf_Val)
output
output <- cbind(" "=rownames(output), output)
write_xlsx(output, "Infl_Mat2.xlsx")



