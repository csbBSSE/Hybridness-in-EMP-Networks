#This code takes RACIPE solution file (which is a matrix) as input and it gives J metric as output. It finds J metric in the following steps:
#Firstly, it finds the correlation (Pearson) matrix of the solution matrix.
#Secondly, it calculates the upper traingular matrix of the correlation matrix.
#Finally, it sums up the elements of the upper triangular matrix (without diagonal elements) by taking the absolute values of the elements.

library(tidyverse)
library(dplyr)
library(magrittr)
library(philentropy)
library(readxl)
library(xlsx)
library(data.table)
library(corrplot)

options(max.print=100000)

#reading solution files of networks in a folder
df1 <- list.files(pattern = "dat")
J_metric <- vector()  #defining an empty vector
for (i in df1) {
  ldf <- read.table(i)
  df2 <- as.data.frame(na.omit(ldf))
  #colnames(df2) <- NULL
  #rownames(df2) <- NULL
  data <- df2[, 4:(ncol(df2)-1)]    #4:6 for GRHL network  #drop snail column in every network
  
#calcualting J metric
  C <- cor(data, method = "pearson")  
  C[upper.tri(C)] <- NA
  C[is.na(C)] <- 0
  CC <- abs(C) 
  J_metric[which(df1==i)] <- sum(CC) - ncol(CC)
  print(J_metric[which(df1==i)])
}

names(J_metric) <- gsub("_solution.dat","",df1)
output <- as.data.frame(J_metric)
output 
output <- cbind(" "=rownames(output), output)
write_xlsx(output, "j_metric_run3.xlsx")
