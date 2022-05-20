source("D:/Teams/Figures/figCore.R")
library(Hmisc)
library(philentropy)
library(future.apply)
library(future)

JSDProper <- function(net, compute = F)
{
    setwd(edgeDel)
    setwd(net)
    jsdFile <- paste0(net, "_JSDdata.csv")
    if(!file.exists(jsdFile) || compute)
    {
        topoFiles <- c(list.files(".", "1_\\d+?.topo"), paste0(net, ".topo"))
        wtFile <- paste0(randRaw, "/", net, "/", net, "_finFlagFreq.csv")
        nodeFile <- paste0(randRaw, "/", net, "/", net, "_nodes.txt")
        if (net == "EMT_MET")
        {
            wtFile <- paste0(randRaw, "/EMT_MET_reduced/EMT_MET_reduced_finFlagFreq.csv")
            nodeFile <- paste0(randRaw, "/EMT_MET_reduced/EMT_MET_reduced_nodes.txt")
        }
            
        wtDfO <- read_csv(wtFile) %>%
            filter(!is.na(Avg0), flag == 1) %>% select(states, Avg0) %>%
            set_names(c("states", "WT"))
        nodes <- readLines(nodeFile) %>% 
            str_replace_all(regex("\\W+"), "")
        
        wtDf <- wtDfO %>% 
            separate(states, c("dummy", "dummy1",all_of(nodes), "dummy2"), sep = "") %>%
            select(-dummy)
        
        JSDs <- sapply(topoFiles, function(x){
            delDf <- read_csv(str_replace(x, ".topo", "_finFlagFreq.csv")) %>% 
                filter(!is.na(Avg0), flag == 1) %>% select(states, Avg0)
            nds <- readLines(str_replace(x, ".topo", "_nodes.txt")) %>% 
                str_replace_all(regex("\\W+"), "")
            wt <- wtDfO
            if (!all(nodes %in% nds))
                wt <- wtDf %>% select(dummy1, all_of(nds), dummy2, WT) %>% 
                unite("states", all_of(c("dummy1", nds, "dummy2")), sep = "") %>%
                group_by(states) %>% summarise(WT = sum(WT))
            df <- merge(wt, delDf, by = "states", all = T)
            df[is.na(df)] <- 0
            JSD(df %>% select(-states) %>% t)
        })
        df <- data.frame(Net = str_remove(topoFiles, ".topo"), JSD = JSDs)
        write.csv(df, paste0(net, "_JSDdata.csv"), row.names = F)
    }
    
}

loopCalc <- function(net, compute = F)
{
    setwd(edgeDel)
    setwd(net)
    loopFile <- paste0(net, "_loopData.csv")
    if (!file.exists(loopFile) || compute)
    {
        command <- paste0("python ../Positive_loops.py ", net, ".topo")
        system(command)
        topoFiles <- c(list.files(".", "1_\\d+?.topo"), paste0(net, ".topo"))
        nets <- topoFiles %>% str_remove(".topo")
        loopData <- sapply(nets, function(x){
            df <- read.csv(paste0(x, "_loops.csv"))
            posLoops <- c(sum(df$Nature == "P"),
                          df %>% filter(Nature == "P") %>% mutate(Edge_count = 1/Edge_count) %>% 
                              select(Edge_count) %>% unlist %>% sum)
            negLoops <- c(sum(df$Nature == "N"),
                          df %>% filter(Nature == "N") %>% mutate(Edge_count = 1/Edge_count) %>% 
                              select(Edge_count) %>% unlist %>% sum)
            posDf <- df %>% filter(Nature == "P") %>% mutate(ID = 1:nrow(.))
            cycles <- posDf$Cycles %>% str_split(",")
            nodes <- unlist(cycles) %>% unique
            hiLoops1 <- lapply(nodes, function(n){
                l <- posDf %>% filter(str_detect(Cycles, n)) %>% select(ID) %>%
                    unlist %>% sort
                data.frame(loopSet = paste0(l, collapse = ","), length = length(l)) 
            }) %>% reduce(rbind.data.frame) %>%
                filter(length > 2) %>% select(loopSet) %>% unlist 
            hiLoops1 <- sapply(1:length(hiLoops1), function(i){
                h <- hiLoops1[i]
                hl <- hiLoops1[-i]
                s <- sum(str_detect(hl, h))
                if (s == 0)
                    T
                else
                    F
            }) %>% sum
            pos2 <- posDf %>% filter(Edge_count == 2)
            
            hiLoops2 <- apply(pos2, 1, function(r){
                nodes <- r[1] %>% str_split(",") %>% unlist
                loops1 <- posDf %>% filter(str_detect(Cycles, nodes[1])) %>% 
                    select(ID) %>% unlist
                loops2 <- posDf %>% filter(str_detect(Cycles, nodes[2])) %>% 
                    select(ID) %>% unlist
                loops1 <- loops1[!(loops1 %in% loops2)]
                loops2 <- loops2[!(loops2 %in% loops1)]
                if (length(loops1) && length(loops2))
                    T
                else
                    F
            }) %>% sum
            c(posLoops, negLoops, hiLoops1, hiLoops2)
        }) %>% t %>% data.frame %>% 
            set_names(c("posLoops", "posWeighted","negLoops", "negWighted","hiLoops1", "hiLoops2")) %>%
            mutate(Net = nets)
        write.csv(loopData, paste0(net, "_loopData.csv"), row.names = F)
    }
    
}

consistencyFormatting <- function(net, edges, edgeNames)
{
    cycles <- read_csv(paste0(net, "_undirectedLoops.csv"), lazy = F) %>% select(Cycles) %>% unlist
    cycleDat <- lapply(cycles, function(cyc){
        nodes <- str_split(cyc, ",")[[1]]
        if(length(nodes) == 1)
        {
            ed <- data.frame(Edge1 = paste0(nodes, ",", nodes))
        }
        else
        {
            nodes <- c(nodes, nodes[1])
            ed <- lapply(1:(length(nodes)-1), function(i){
                e <- c(nodes[i], nodes[i+1])
                e <- list(e, rev(e)) %>% sapply(paste0, collapse = ",")
                e <- e[which(e %in% edgeNames)]
            }) %>% expand.grid %>% set_names(paste0("Edge", 1:ncol(.)))
        }
        
        ed %>% mutate(Nature = ed %>% apply(1, function(e){edges[e] %>% prod})) %>%
            mutate(Nature = ifelse(Nature == 1, "P", "N")) %>% 
            unite("EdgeList",contains("Edge"), sep = ";") %>%
            mutate(Cycles = cyc, Edge_count = str_count(EdgeList, ";") + 1)
    }) %>% reduce(rbind.data.frame)
    cycleDat
}

consistency <- function(net, compute = F)
{
    setwd(edgeDel)
    setwd(net)
    loopFile <- paste0(net, "_Consistency.csv")
    if (is_empty(loopFile) || compute)
    {
        # browser()
        command <- paste0("python ../signConsistency.py undirectedLoops ", net, ".topo")
        system(command)
        topoFiles <- c(list.files(".", "1_\\d+?.topo"), paste0(net, ".topo"))
        nets <- topoFiles %>% str_remove(".topo")
        loopData <- sapply(nets, function(x){
            topoDf <- read.delim(paste0(x, ".topo"), sep = "") %>% 
                unite(col = "Edge",sep = ",", Source, Target) %>% 
                mutate(Type = ifelse(Type == 2, -1, 1))
            edges <- topoDf$Type
            names(edges) <- topoDf$Edge
            edgeNames <- names(edges)
            setwd("undirectedLoops")
            df <- consistencyFormatting(x, edges, edgeNames)
            setwd("..")
            # browser()
            # d <- read.csv(paste0("undirectedLoops/",x, "_undirectedLoops.csv"))
            posLoops <- c(sum(df$Nature == "P"), 
                          df %>% filter(Nature == "P") %>% mutate(Edge_count = 1/Edge_count) %>% 
                              select(Edge_count) %>% unlist %>% sum)
            negLoops <- c(sum(df$Nature == "N"), 
                          df %>% filter(Nature == "N") %>% mutate(Edge_count = 1/Edge_count) %>% 
                              select(Edge_count) %>% unlist %>% sum)
            fracNeg <- negLoops[1]/(negLoops[1] + posLoops[1])
            print(x)
            write.csv(df, paste0("undirectedLoops/",x, "_undirectedLoopsF.csv"), row.names = F)
            df$ID <- 1:nrow(df)
            edgeDat <- lapply(edgeNames, function(e){
                df %>% filter(str_detect(EdgeList, e)) %>% group_by(Nature) %>%
                    summarise(Count = n(), loops = paste0(ID, collapse = ",")) %>% 
                    mutate(Edge = e)
            }) %>% reduce(rbind.data.frame) %>% 
                mutate(Nature = factor(Nature, levels = c("P", "N"))) %>% 
                spread(key = Nature, value = Count)
            if (!("P" %in% colnames(edgeDat)))
                edgeDat$P <- 0
            if (!("N" %in% colnames(edgeDat)))
                edgeDat$N <- 0
            edgeDat[is.na(edgeDat)] <- 0
            negDat <- edgeDat %>% filter(P == 0) %>% arrange(-N)
            posDat <- edgeDat %>% filter(N == 0)
            mixDat <- edgeDat %>% filter(P != 0, N != 0)
            nLoops <- df$ID[df$Nature == "N"]
            coveredLoops <- c()
            swapEdges <- 0
            if(nrow(negDat) != 0)
            {
                dummy <- sapply(negDat$loops, function(l){
                    l <- l %>% str_split(",") %>% unlist %>% as.integer
                    i <- intersect(coveredLoops, l)
                    if (is.null(i) || length(i) < length(l))
                    {
                        coveredLoops <<- c(coveredLoops, l) %>% unique
                        swapEdges <<- swapEdges + 1
                    }
                })
                if (length(coveredLoops) < length(nLoops))
                {
                    remaining <- nLoops[!(nLoops %in% coveredLoops)]
                    dat <- edgeDat %>% 
                        filter(sapply(loops, function(l){
                            any(str_detect(l, as.character(remaining)))
                        }))
                    loops <- paste0(dat$loops, collapse = ",") %>% str_split(",") %>%
                        unlist %>% as.integer
                    pLoops <- loops[!(loops %in% remaining)]
                    dummy <- sapply(dat$loops, function(l){
                        l <- l %>% str_split(",") %>% unlist %>% as.integer
                        l <- l[!(l %in% pLoops)]
                        i <- intersect(coveredLoops, l)
                        if (is.null(i) || length(i) < length(l))
                        {
                            coveredLoops <<- c(coveredLoops, l) %>% unique
                            swapEdges <<- swapEdges + 1
                            return(T)
                        }
                        return(F)
                    })
                    ed <- dat$Edge[dummy]
                    dat <- edgeDat %>% 
                        filter(sapply(loops, function(l){
                            any(str_detect(l, as.character(pLoops)))
                        })) %>% filter(!(Edge %in% ed))
                    coveredLoops <- c()
                    dummy <- sapply(dat$loops, function(l){
                        l <- l %>% str_split(",") %>% unlist %>% as.integer
                        i <- intersect(coveredLoops, l)
                        if (is.null(i) || length(i) < length(l))
                        {
                            coveredLoops <<- c(coveredLoops, l) %>% unique
                            swapEdges <<- swapEdges + 1
                            return(T)
                        }
                        return(F)
                    })
                    
                }
            }
            c(posLoops, negLoops, swapEdges)
        }) %>% t %>% data.frame %>% 
            set_names(c("posLoopsUnd", "posLoopsUndWeighted","negLoopsUnd", "negLoopsUndWeighted",
                        "inconsistency")) %>%
            mutate(Net = nets)
        write.csv(loopData, paste0(net, "_Consistency.csv"), row.names = F)
    }
}

correlationMatBool <- function(topoFile)
{#browser()
    print(topoFile)
    net <- topoFile %>% str_remove(".topo")
    nodes <- readLines(topoFile %>% str_replace(".topo", "_nodes.txt"))
    corMat <- read_csv(topoFile %>% str_replace(".topo", "_finFlagFreq.csv"), 
                       show_col_types = F) %>% 
        filter(flag == 1) %>% select(states, Avg0) %>% drop_na %>%
        mutate(Avg0 = Avg0*10000 %>% round) %>% apply(1, function(x){
            s <- x[1] %>% str_remove_all("'")
            n <- x[2] %>% as.integer
            rep(s, n)
        }) %>% unlist 
    if (length(corMat) < 3)
        return(NA)
    corMat <- corMat %>% lapply(function(x){
        str_split(x, "") %>% unlist %>% as.integer
    }) %>% reduce(rbind.data.frame) %>% 
        set_names(nodes) %>% cor
    colnames(corMat) <- colnames(corMat) %>% str_replace_all(regex("\\W+"), "")
    rownames(corMat) <- rownames(corMat) %>% str_replace_all(regex("\\W+"), "")
    if(!dir.exists("Correlations"))
        dir.create("Correlations")
    # browser()
    write.csv(data.frame(corMat), paste0("Correlations/", net, "_cor.csv"))
}
correlationMatBool <- cmpfun(correlationMatBool)

InflCorrelHybrid <- function(net)
{
    setwd(edgeDel)
    setwd(net)
    topoFiles <- c(list.files(".", "1_\\d+?.topo"), paste0(net, ".topo"))
    df <- sapply(topoFiles, function(x){
        n <- x %>% str_remove(".topo")
        topoFile <- paste0(n, ".topo")
        InflFile <- paste0("Influence/", n, "_fullMat.csv")
        if(!file.exists(InflFile))
        {
            influence_matrix(topoFile)
        }
        InflMat <- read.csv(InflFile) %>% select(-1) %>% as.matrix %>% as.vector
        InflMetric <- InflMat[which(InflMat<0)] %>% sum
        CorFile <- paste0("Correlations/", n, "_cor.csv")
        if(!file.exists(CorFile))
            correlationMatBool(topoFile)
        corDat <- read.csv(CorFile, row.names = 1) %>% as.matrix
        JMetric <- corDat[upper.tri(corDat)] %>% sum
        JMetric <- JMetric - diag(corDat) %>% sum
        Hybrid <- read.csv(paste0(n, "_finFlagFreq.csv")) %>% 
            filter(!is.na(Avg0), flag == 1, abs(EMTScore) <0.5) %>% 
            select(Avg0) %>% unlist %>% sum
        c(InflMetric, JMetric, Hybrid)
    }) %>% t %>% 
        data.frame %>% set_names(c("InflMetric", "JMetric", "hybridFreq")) %>% 
        mutate(Net = topoFiles %>% str_remove(".topo"))
    write.csv(df, paste0(net, "_InflCorHybrid.csv"), row.names = F)
}

metricData <- function(net, method = "pearson", compute = F)
{
    if (compute)
    {
        JSDProper(net)
        loopCalc(net)
        consistency(net)
        InflCorrelHybrid(net)
    }
    
    setwd(edgeDel)
    setwd(net)
    frustCoh <- read_csv(paste0("../", net, "_ALL.csv")) %>% 
        select(minFrust, maxFrust, minFrust, minCoh, maxCoh, meanCoh, Net)
    inflCor <- read.csv(paste0(net, "_InflCorHybrid.csv"))
    loopDat <- read.csv(paste0(net, "_loopData.csv"))
    cons <- read.csv(paste0(net, "_Consistency.csv"))
    jsdDat <- read.csv(paste0(net, "_JSDdata.csv"))
    nets <- jsdDat$Net
    gs <- read.csv(paste0("../", ifelse(net == "EMT_MET", "EMT_MET_reduced", net), 
                          "_GroupMetrics2.csv")) %>%
        mutate(Gs = (abs(G11) + abs(G12) + abs(G22) + abs(G21))/4) %>% 
        select(Net, Gs)
    df <- list(frustCoh, inflCor, loopDat, jsdDat, cons, gs) %>% 
        reduce(merge, by = "Net", all = T)
    write.csv(df, paste0(net, "_mubasherAll.csv"), row.names = F)
    net1 <- ifelse(net == "EMT_MET", "EMT_MET_reduced", net)
    wd <- "D:\\Github\\Projects\\Ongoing\\Mubasher\\FinalResults\\Boolean"
    write.csv(df, paste0(wd, "/", net1, "_all.csv"), row.names = F)
    cordat <- df %>% select(-Net) %>% as.matrix %>% rcorr(type = method)
    corDf <- cordat$r %>% data.frame %>% 
        mutate(M1 = rownames(.)) %>% 
        gather(key = "Metric", value = "Correlation", -M1)
    corP <- cordat$P %>% data.frame %>%
        mutate(M1 = rownames(.)) %>%
        gather(key = "Metric", value = "pValue", -M1) %>%
        mutate(Significance = ifelse(pValue < 0.05, "", "X")) %>% 
        select(Metric, M1, Significance)
    df <- merge(corDf, corP, by = c("Metric", "M1"), all = T)
    ggplot(df, aes(x = Metric, y = M1, fill = Correlation)) +
        geom_tile() + geom_text(aes(label = Significance)) +
        theme_Publication() + 
        theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
              axis.text = element_text(size = rel(1)),
              legend.position = "right",
              legend.direction = "vertical", legend.key.height = unit(0.8, "cm")) +
        scale_fill_gradient2(low = "red", high = "blue", limits = c(-1,1)) + 
        labs(x ="", y = "")
    ggsave(paste0(net, "_metricCorrel.png"), width = 6, height = 5)
    print(net)
}
sapply(c(EMPNets[1:4],"EMT_MET"), metricData, compute = T)
