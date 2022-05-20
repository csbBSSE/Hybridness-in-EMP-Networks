library(tidyverse)
library(miscFuncs)

networks <- c("grhl", "ovol", "oct", "nrf", "grhlovol", "grhlovoloct", "grhlovolnrf",
              "grhlovoloctnrf", "silviera", "silviera2", "EMT_RACIPE", "EMT_RACIPE2", "EMT_MET")
EnodeList <- list(c("GRHL2", "miR200"), c("OVOL2", "miR200"), c("miR200", "miR145"), 
                  c("Ecadherin", "miR34"), c("GRHL2", "miR200", "OVOL2"),
                  c("GRHL2", "miR200", "OVOL2", "miR145"), c("GRHL2", "OVOL2", "Ecad", "miR34"),
                  c("GRHL2", "miR200", "OVOL2", "miR34", "miR145", "Ecad"),
                  c("MIR-34", "MIR-200", "CDH1", "MIR-1199", "GRHL2", "ESRP1", "MIR-190", 
                    "MIR-340", "CD24"), 
                  c("MIR-34", "MIR-200", "CDH1", "MIR-1199", "GRHL2", "ESRP1", "MIR-190", 
                    "MIR-340", "CD24", "MIR-129"),              
    c("CDH1", "miR-101", "miR-141", "miR-200a", "miR-200b", "miR-200c", "miR-34a"),
                  c("CDH1", "miR-101", "miR-141", "miR-200a", "miR-200b", "miR-200c", "miR-205", 
                    "miR-30c", "miR-34a", "miR-9", "OVOL2", "GRHL2", "NP63", "CLDN7"),
                  c("GSK3", "KLF4", "Ecadherin", "cateninmemb", "SUFU", "RKIP", "Patched", 
                    "TrCP", "miR200"))

MnodeList <- list(c("ZEB"), c("ZEB"), c("ZEB", "OCT4"),c("ZEB", "KEAP1", "NRF2", "SNAIL"),
                  c("ZEB"),c("OCT4", "ZEB"), c("ZEB", "KEAP1", "NRF2", "SNAIL"),
                   c("ZEB", "KEAP1", "NRF2", "SNAIL", "OCT4"),
    c("EX-TGFB", "SNAIL1", "VIM", "ZEB1", "HAS2", "TGFB", "EMT", "HA", "CD44s"),
                  c("EX-TGFB", "SNAIL1", "VIM", "ZEB1", "HAS2", "TGFB", "EMT", "HA", "CD44s", "TWIST1"),
                  c("FOXC2", "ZEB1", "ZEB2", "SNAI1", "TGF-beta", "SNAI2", "TWIST2", "VIM", "TWIST1"), 
                  c("FOXC2", "ZEB1", "KLF8", "ZEB2", "SNAI1", "TGF-beta", "SNAI2", "TWIST2", "VIM", 
                    "TWIST1", "TCF3", "GSC"), 
                  c("ILK","AKT","PI3K","AXIN2","TCFLEF","CD44","TGFR","CDC42","NOTCHic","Csl","NFB",
                    "Csn","RAS", "DELTA","Frizzled","DSH","cateninnuc", "Destcompl","TWIST1","SNAI2",
                    "SNAI1","ZEB1","FOXC2","HEY1","cfos","EGR1","MEK","ERK","SMO","FUS","Wnt",
                    "GLI","SMAD","IKK","Jagged","STAT","ZEB2","LIV1","RAF", "SOSGRB2","NOTCH",
                    "PAK1","SHH","SRC","LOXL23","cMet","JAK","TGF")
                  )

names(EnodeList) <- names(MnodeList) <- networks

netNames <- c("4N 7E", "4N 9E", "5N 10E", "8N 16E", "18N 33E", 
              "5N 11E", "7N 16E", "10N 22E", "12N 27E",
              "20N 40E", "22N 82E", "26N 100E", "57N 113E")
names(netNames) <- networks


colN <- c("Network", "JSD", "Jmetric", "normJ", "InfVal", "InfVal2",
          "normI", "minFrust", "maxFrust", "meanFrust", "posLoops", "normPL", "negLoops",
          "normNL", "Plasticity", "PlasticityHybrid", "HOC","Hybridness", "PFL", "Type1", "Type2", 
          "MISSA", "Gs", "InflVal", "InflVal2", "MaxCoh", "MaxFrust", "MinCoh", "MinFrust",
          "posWeighted", "negWeighted", "posLoopsUnd", "negLoopsUnd", "posLoopsUndWeighted", 
          "negLoopsUndWeighted", "inconsistency", "inflts", "inflda", "inflposloops", 
          "inflnegloops", "normplweighted", "normplund", "normplinfl", "normplundweighted",
          "predFrust")
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
             "Predicted Frustration")
names(colNKey) <- colN %>% tolower

correlGrob <- function(df, x, y, xPos = NULL, yPos = NULL, method = "pearson")
{
    corr <- cor.test(df[[x]] %>% as.numeric, as.numeric(df[[y]]), method = method)
    pVal <- ifelse(corr$p.value < 0.05, "*", "")
    xPos <- ifelse(!is.null(xPos), xPos, 0.5)
    yPos <- ifelse(!is.null(yPos), yPos, 0.9)
    grob <- grobTree(textGrob(paste0("\u03c1 : ", round(corr$estimate, 2), pVal), 
                              x=xPos,  y=yPos, hjust=0,
                              gp=gpar(col="black", fontsize=18, fontface="bold")))
    return(grob)
}

stringBreak <- function(s, maxPieces = 3, maxLength = 15) {
    sPieces <- s %>% str_split(" ") %>% unlist
    if (length(sPieces) == 2)
        return(str_replace(s, " ", "\n"))
    i <- 1
    pieceList <- sapply(1:(maxPieces-1), function(x){
        # browser()
        thisPiece <- ""
        while(nchar(paste(thisPiece, sPieces[i]))<maxLength && i <= length(sPieces)) {
            if (thisPiece == "")
                thisPiece <- sPieces[i]
            else
                thisPiece <- paste(thisPiece, sPieces[i])
            i <<- i + 1
        }
        thisPiece
    })
    j <- i-1
    if (j < length(sPieces)) {
        pieceList[maxPieces] <- paste(sPieces[i:(length(sPieces))], collapse = " ")
    }
    pieceList %>% paste(collapse = "\n")
}
stringBreak <- Vectorize(stringBreak, vectorize.args = "s")
