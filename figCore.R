# This file initiates the analysis by setting the folder structure. Please create 
# corresponding folders and place the data in the folders according to the instructions
# given below.

library(tidyverse)
library(compiler)
library(ggrepel)
library(grid)
library(ggcorrplot)
options(stringsAsFactors = F)
mainFolder <- "D:/Github/Projects/Ongoing/Mubasher" 
cwd <- paste0(mainFolder, "/Figures")
dataFolder <- paste0(mainFolder, "/FinalResults") #This folder contains the data used to generate figures. Uploded on github
rac <- paste0(dataFolder, "/RACIPE") #This folder contains the summary of results for RACIPE simulations
bool <- paste0(dataFolder, "/Boolean") # Same as before but for boolean simulations
netMetrics <- paste0(dataFolder, "/NetMetrics") # For all edge perturbation networks, the topological metrics are store here. Can be calculated using provided codes
racipeRaw <- paste0(mainFolder, "/RACIPE/RawData") 
boolRaw <- paste0(mainFolder, "/Boolean")
topos <- paste0(dataFolder, "/TopoFiles") # Contains all topo Files. 
# hiLoops <- paste0(mainFolder, "/RACIPE/hiLoops") 
hiLoopCode <- paste0(mainFolder, "/HiLoop/countmotifs.py") # the code used to calculate hiloops. 

source(paste0(cwd, "/utils.R")) # some useful functions. The script is uploaded in the main folder of the github repo. 

nets <- c("grhl", "ovol", "oct", "nrf","grhlovol",
          "grhlovoloct", "grhlovolnrf", "grhlovoloctnrf","SIL", "SIL2",
          "EMT_RACIPE", "EMT_RACIPE2", 
          "EMT_MET_reduced")
netKey <- sapply(paste0(nets, ".topo"), NodeEdges, topoPath = topos)
names(netKey) <- nets
EMPNets <- nets[9:13]

netKey
source(paste0(cwd,"/plot_theme.R"))
source(paste0(cwd,"/utils.R"))


colN <- c("Network", "JSD", "Jmetric", "normJ", "InfVal", "InfVal2",
          "normI", "minFrust", "maxFrust", "meanFrust", "posLoops", "normPL", "negLoops",
          "normNL", "Plasticity", "PlasticityHybrid", "HOC","Hybridness", "PFL", "Type1", "Type2", 
          "MISSA", "Gs", "InflVal", "InflVal2", "MaxCoh", "MaxFrust", "MinCoh", "MinFrust",
          "posWeighted", "negWeighted", "posLoopsUnd", "negLoopsUnd", "posLoopsUndWeighted", 
          "negLoopsUndWeighted", "inconsistency", "inflts", "inflda", "inflposloops", 
          "inflnegloops", "normplweighted", "normplund", "normplinfl", "normplundweighted",
          "predFrust", "smallPFL")
colNKey <- c("Network", "JSD", "J Metric", "J Metric (normalized)",
             "Influence Metric", "Influence Metric (negative)", "Influence Metric (normalized)",
             "Minimum Frustration", "Maximum Frustration", "Mean Frustration", "Positive Feedback Loops (PFLs)", 
             "Fraction of PFLs", "Negative Feedback Loops (NFLs)", "Fraction NFLs", "Multistability", 
             "Hybrid Plasticity", 
             "HOC", "Hybridness", "PFL", "Type1", "Type2", "MISSA", "Groupstrength", 
             "Influence Metric", "Influence Metric (negative)",
             "Maximum Coherence", "Maximum Frustration", "Minimum Coherence",
             "Minimum Frustration", 
             "Weighted PFL", "Weighted NFL", "Undirected PFL", "Undirected NFL", "Undirected WPFL",
             "Undirected WNFL", "Inconsistency", "Toggle Switches (Influence)", 
             "Double Activation (Influence)","Positive Loops (Influence)", "Negative Loops (Influence)",
             "Fraction of Weighted PFLs", "Fraction of undirected PFLs", 
             "Fraction of Influence PFLs", "Fraction of undirected weighted PFLs", 
             "Predicted Frustration", "Small PFL")
names(colNKey) <- colN %>% tolower

