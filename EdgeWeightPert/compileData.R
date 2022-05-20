wd <- "D:/Github/Projects/Ongoing/Mubasher/EdgeWeightPert"

setwd(wd)
source("D:/Github/Projects/Ongoing/Mubasher/Figures/figCore.R")
source("setup.R")
source("influence.R")
source("loopDat.R")
source("corJSD.R")


phenotype <- function(stateList, Enodes, Mnodes, nodeOrder) {
    stateList <- stateList %>% str_remove_all("'") %>% 
        str_split("") %>% sapply(as.integer) %>% t %>% data.frame %>%
        set_names(nodeOrder)
    Escore <- stateList %>% select(all_of(Enodes)) %>% rowSums
    Mscore <- stateList %>% select(all_of(Mnodes)) %>% rowSums
    Score <- Escore/length(Enodes) - (Mscore/length(Mnodes))
    phenotype <- rep("H", length(Score))
    phenotype[(abs(Score) > 0.5 & Score < 0)] <- "M"
    phenotype[(abs(Score) > 0.5 & Score < 0)] <- "E"
    return(phenotype)
}

infMat <- function(Enodes, Mnodes) {
    fil <- list.files(".", "_new.csv")
    df <- read_csv(fil, col_types = cols())
    dfTopo <- df %>% select(Source, Target)
    df <- df %>% select(-Source, -Target)
    dfInfl <- sapply(df, function(x) {
        InflLoopCalc(dfTopo, x, Enodes, Mnodes)
    }) %>% t %>% data.frame %>% set_names(c("inflnegloops", "inflposloops"))
    dfInfl$Network <- paste0("rand_", colnames(df))
    dfInfl
}




compileDat <- function(net, Enodes, Mnodes) {
    setwd(paste0(net, "_rand"))
    nodes <- list.files(".", "nodes.txt") %>% readLines
    filz <- list.files(".", "finFlagFreq")
    hybridFrust <- sapply(filz, function(x) {
        print(x)
        d <- read_csv(x, col_types = cols()) %>% filter(flag == 1) %>% drop_na
        if (nrow(d) ==0)
            return(c(NA, NA, NA, NA))
        d <- d %>% 
            mutate(Phenotype = phenotype(states, Enodes, Mnodes, nodes))
        hybridness <- d %>% filter(Phenotype == "H") %>% select(Avg0) %>%
            unlist %>% sum
        frustration <- c(min(d$frust0), max(d$frust0), sum(d$frust0*d$Avg0))
        c(hybridness, frustration)
    }) %>% t %>% data.frame %>% set_names(c("Hybridness", "minFrust", "maxFrust", "meanFrust")) %>% 
        mutate(Network = filz %>% str_remove("_finFlagFreq.csv") %>% str_replace(net, "rand"))
    inflDat <- infMat(Enodes, Mnodes)
    loops <- loopsDat(net)
    jsd <- JSDProper(net)
    jmetric <- corDat(net)
    # browser()
    df <- merge(hybridFrust, inflDat, by = "Network", all = T) %>%
        merge(loops, by = "Network", all = T) %>%
        merge(jsd, by = "Network", all = T) %>%
        merge(jmetric, by = "Network", all = T)
    setwd("..")
    write_csv(df, paste0(net, "_all.csv"), quote = "none")
}
d <- sapply(networks[-(1:8)], function(net) {
    compileDat(net, EnodeList[[net]], MnodeList[[net]])
})
