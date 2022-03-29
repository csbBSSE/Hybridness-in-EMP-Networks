#This code implements higher order clustering method on the correlation matrix (calculated from the solution matrix of the wild-type network) 
#to identify epithelial (E) and mesenchymal (M) nodes in the wild-type network. Using this characteristic defintion of E and M nodes, it then 
#calculates Epithelial-Mesenchymal Transition (EMT) Score of the wild-type network and its perturbations by subtrating the average expressions 
#of all the epithelial nodes from the average expressions of all the mesenchymal nodes. Based on the EMT score, lying between -1 and +1, it 
#characterizes the SS as the hybrid state if its EMT score lies between -0.5 and +0.5.

@author: Kishore

library(readxl)
library(tidyverse)
library(magrittr)
library(dplyr) 

#reading the solution file of wild-type network file and doing correlations between the genes
ldf <- read.table("wild_solution.dat")  
df2 <- as.data.frame(na.omit(ldf))
#colnames(df2) <- NULL
#rownames(df2) <- NULL
data <- df2[, 4:ncol(df2)]              

C <- cor(data, method = "pearson")  

#colnames(C) <- c("miR200", "ZEB", "GRHL2", "OVOL2","miR145","OCT4", "Ecad", "NRF2", "SNAIL", "KEAP1", "miR34", "Iext")
#rownames(C) <- c("miR200", "ZEB", "GRHL2", "OVOL2","miR145","OCT4", "Ecad", "NRF2", "SNAIL", "KEAP1", "miR34", "Iext")

#categorizing genes as epithelial and mesenchmal markers
hc <- hclust(dist(C))
plot(hc)
groupOrder <- cutree(hc, 2)
group1 <- colnames(C)[groupOrder == 1]
group2 <- colnames(C)[groupOrder == 2]

gro1 <- shQuote(group1, type = "cmd")
gr1 <- paste(gro1, collapse = ", ")
cat(gr1)

gro2 <- shQuote(group2, type = "cmd")
gr2 <- paste(gro2, collapse = ", ")
cat(gr2)

#Function to calculate EMT score
EMscore <- function(state, EpithelialNodes, MesenchymalNodes, NodeOrder)
{
    epi <- state[which(NodeOrder %in% EpithelialNodes)] %>% as.integer %>%
        sum
    Mes <- state[which(NodeOrder %in% MesenchymalNodes)] %>% as.integer %>%
        sum
    score <- epi/length(EpithelialNodes) - Mes/length(MesenchymalNodes)
    score
}

hybridFreqScoreBased <- function(EpithelialNodes, MesenchymalNodes, NodeOrder)
{
    filz <- list.files(".", "xlsx")
    hybrid <- lapply(filz, function(x){
        net <- x %>% str_remove(".xlsx")
        print(net)
        df <- read_excel(x, col_names = F)kish
        if(ncol(df) == (length(NodeOrder)+1))
        {
            df[[ncol(df) + 1]] <- df[[ncol(df)]]
        }
            df <- set_names(df, c(NodeOrder, "WT", "Mut"))
        df1 <- df %>% select(all_of(NodeOrder))
        df1$Score <- df1 %>% apply(1, EMscore, EpithelialNodes = EpithelialNodes, 
                                   MesenchymalNodes = MesenchymalNodes, 
                                   NodeOrder = NodeOrder)
        df1 <- df1 %>% unite("State", NodeOrder, sep = "")
        
        df <- df %>% unite("State", NodeOrder, sep = "") %>%
            merge(df1, by = "State", all = T) %>%
            mutate(Hybrid = ifelse(abs(Score)<0.5, "Hybrid", "Terminal")) %>%
            group_by(Hybrid) %>%
            summarise(WT = sum(WT), Mut = sum(Mut)) %>%
            mutate(WT = WT/sum(WT), Mut = Mut/sum(Mut)) %>%
            mutate(Network = net)
        
        df
    }) %>% reduce(rbind.data.frame)
    write.csv(hybrid, "HybridFreq.csv", row.names = F)
    hybrid
}

#Calling the function to calculte EMT scores, network-wise

#Run for GRHL/OVOL/OCT/NRF/grhl-ovol et al.

hybridFreqScoreBased(c("miR200", "OVOL2"), c("ZEB"), c("miR200", "ZEB", "OVOL2")) 
                     

#hybridFreqScoreBased(c("miR200", "GRHL2", "OVOL2", "miR145", "Ecad", "miR34"), c("ZEB", "OCT4", "NRF2", "SNAIL", "KEAP1"),
                     #c("miR200", "ZEB", "GRHL2", "OVOL2","miR145","OCT4", "Ecad", "NRF2", "SNAIL", "KEAP1", "miR34", "Iext")) 
                     
#Run for SIL2 
hybridFreqScoreBased(c("MIR-34", "MIR-200", "CDH1", "MIR-1199", "GRHL2", "ESRP1", "MIR-190", "MIR-340", "CD24", "MIR-129"), 
                     c("EX-TGFB", "SNAIL1", "VIM", "ZEB1", "HAS2", "TGFB", "EMT", "HA", "CD44s", "TWIST1"), 
                     c("EX-TGFB", "SNAIL1", "MIR-34", "MIR-200", "CDH1", "VIM", "ZEB1", "MIR-1199", "GRHL2", "ESRP1",  "HAS2", 
                       "MIR-190", "MIR-340", "TGFB", "CD24", "EMT", "HA", "CD44s", "TWIST1", "MIR-129"))


#Run for EMT
hybridFreqScoreBased(c("CDH1", "miR-101", "miR-141", "miR-200a", "miR-200b", "miR-200c", "miR-34a"), 
                     c("FOXC2", "ZEB1", "ZEB2", "SNAI1", "TGF-beta", "SNAI2", "TWIST2", "VIM", "TWIST1"), 
                     c("FOXC2", "ZEB1", "KLF8", "CDH1", "miR-101", "ZEB2", "SNAI1", "miR-141", "TGF-beta", "miR-200a",  "miR-200b", 
                       "miR-200c", "miR-205", "miR-30c", "SNAI2", "miR-34a", "TWIST2", "miR-9", "VIM", "TWIST1", "TCF3", "GSC"))

#run for sil
hybridFreqScoreBased(c("MIR-34", "MIR-200", "CDH1", "MIR-1199", "GRHL2", "ESRP1", "MIR-190", "MIR-340", "CD24"), 
c("EX-TGFB", "SNAIL1", "VIM", "ZEB1", "HAS2", "TGFB", "EMT", "HA", "CD44s"), 
c("EX-TGFB", "SNAIL1", "MIR-34", "MIR-200", "CDH1", "VIM", "ZEB1", "MIR-1199", 
               "GRHL2", "ESRP1", "HAS2", "MIR-190", "MIR-340", "TGFB", "CD24", "EMT", "HA", "CD44s"))

#run for EMT_RACIPE2
hybridFreqScoreBased(c("CDH1", "miR-101", "miR-141", "miR-200a", "miR-200b", "miR-200c", "miR-205", "miR-30c", 
                       "miR-34a", "miR-9", "OVOL2", "GRHL2", "NP63", "CLDN7"), 
                     c("FOXC2", "ZEB1", "KLF8", "ZEB2", "SNAI1", "TGF-beta", "SNAI2", "TWIST2", "VIM", "TWIST1", "TCF3", "GSC"), 
                     c("FOXC2", "ZEB1", "KLF8", "CDH1", "miR-101", "ZEB2", "SNAI1",	"miR-141", "TGF-beta", "miR-200a", "miR-200b", "miR-200c", 
                       "miR-205", "miR-30c", "SNAI2",	"miR-34a", "TWIST2", "miR-9", "VIM", "TWIST1", "TCF3", "GSC", "OVOL2",	"GRHL2", "NP63", "CLDN7"))

#run for EMT_MET
hybridFreqScoreBased(c("GSK3", "KLF4", "Ecadherin", "cateninmemb", "SUFU", "RKIP", "Patched", "TrCP", "miR200"), 
                     c("ILK","AKT","PI3K","AXIN2","TCFLEF","CD44","TGFR","CDC42","NOTCHic","Csl","NFB","Csn","RAS", "DELTA",
                       "Frizzled","DSH","cateninnuc", "Destcompl","TWIST1","SNAI2","SNAI1","ZEB1","FOXC2","HEY1","cfos","EGR1", 
                       "MEK","ERK","SMO","FUS","Wnt","GLI","SMAD","IKK","Jagged","STAT","ZEB2","LIV1","RAF", "SOSGRB2","NOTCH",
                       "PAK1","SHH","SRC","LOXL23","cMet","JAK","TGF"), 
                     c("ILK","AKT","PI3K","AXIN2","TCFLEF","CD44","TGFR","CDC42","NOTCHic","Csl","NFB","Csn","RAS","DELTA","Frizzled",
                       "DSH","cateninnuc","Destcompl", "GSK3","KLF4","Ecadherin","cateninmemb","TWIST1","SNAI2","SNAI1","ZEB1","FOXC2",
                       "HEY1","cfos","EGR1","MEK","ERK","SMO","FUS","Wnt","GLI","SUFU","SMAD","IKK","Jagged","STAT","ZEB2","LIV1","RAF",
                       "SOSGRB2","RKIP","NOTCH","PAK1","SHH","Patched","SRC","LOXL23","TrCP","cMet","JAK","TGF","miR200"))


#Segregating the hybrid state frequencies from the terminal state (epithelial and mesenchymal) frequencies
csv_file <- read.csv("HybridFreq.csv")
colnames(csv_file) <- NULL
toDelete <- seq(0, nrow(csv_file), 2)     #index on even rows to be deleted
hybrid_filtered <-  csv_file[-toDelete, ] #deleting even rows
write.csv(hybrid_filtered, "HybridFreq_filtered_run3.csv", row.names = F)
hybrid_filtered
