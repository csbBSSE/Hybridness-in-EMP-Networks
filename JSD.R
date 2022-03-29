#This code uses JSD function to compare the steady state frequency distribution of any wild-type network with the steady state frequency 
#distribution of its perturbations.

@author: Mubasher

library(tidyverse)
library(dplyr)
library(magrittr)
library(philentropy)
library(readxl)
library(writexl)
library(xlsx)
#library(naturalsort)

#setwd("F:/NPDF-SERB-DST/RACIPE-ANALYSIS") ###this can also be done using session menu at the top###

options(max.print=100000)

#loading data files containing solution files of networks
df1 <- list.files(pattern = "xlsx")
jsdd <- vector()
for (i in df1) {
  ldf <- read_excel(i)
  df <- as.data.frame(na.omit(ldf))
  #colnames(df) <- NULL
  df2 <- df[, c((ncol(df)-1):ncol(df))] #Taking only last two columns viz. frequencies of wild-type and perturbed network 
  P <- t(df2[1])
  Q <- t(df2[2])
  x <- rbind(P, Q)
  jsdd[which(df1==i)] <- JSD(x, unit = "log2", est.prob = "empirical") #calculating JSD between the two frequency distributions
  print(jsdd[which(df1==i)])
}

names(jsdd) <- gsub(".xlsx","",df1)
output <- as.data.frame(jsdd)
output 
output <- cbind(" "=rownames(output), output)  #Binding rownames (which are network names) to the jsdd values. 
write_xlsx(output, "jsdd_run_3.xlsx")   
