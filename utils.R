source("inflMat.R")

stringBreak <- function(s, maxPieces = 3, maxLength = 15) {
    sPieces <- s %>% str_split(" ") %>% unlist
    if (length(sPieces) <= 2)
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

naZero <- function(df)
{
    df[is.na(df)] <- 0
    df
}

racReader <- function(f)
{
    net <- str_remove(f, "_all.csv")
    df <- read_csv(f, lazy = F, col_types = cols()) %>% set_names(tolower(colnames(.))) %>%
        select(network, jsd, jmetric, hybridness, plasticity, type1, type2, missa)
    rem <- ifelse(net == "EMT_RACIPE", "EMT", net)
    df2 <- read_csv(paste0(rac, "/Frustration/", net, "_PlasticityFrustHybrid.csv"), lazy = F, 
                    col_types = cols()) %>% set_names(tolower(colnames(.))) %>%
        select(network, minfrust, maxfrust, meanfrust) %>% 
        mutate(network = str_replace_all(network, "-", "_") %>% str_remove(paste0(rem, "_")))
    df <- merge(df, df2, by = "network", all = T)
    dfMetrics <- read_csv(paste0(netMetrics, "/", 
                                 str_replace(f, "_all", "_netMetrics")), lazy = F, col_types = cols())
    colnames(dfMetrics) <- dfMetrics %>% colnames %>% tolower
    dfMetrics <- dfMetrics %>%
        mutate(normpl = posloops/(posloops + negloops), 
               normplweighted = posweighted/(posweighted+negweighted),
               normplund = posloopsund/(posloopsund + negloopsund),
               normplundweighted = posloopsundweighted/(posloopsundweighted + negloopsundweighted),
               normplinfl = inflposloops/(inflposloops + inflnegloops))
    merge(df, dfMetrics, by = "network", all = T)
}
boolReader <- function(f)
{
    df <- read_csv(f, lazy = F, col_types = cols()) %>% set_names(tolower(colnames(.)))
    df2 <- read_csv(paste0(rac,"/", f), col_types = cols()) %>%
        set_names(tolower(colnames(.))) %>%
        select(network,type1, type2, missa)
    df <- merge(df, df2, by = "network", all = T)
    dfMetrics <- read_csv(paste0(netMetrics, "/", 
                                 str_replace(f, "_all", "_netMetrics")), lazy = F, col_types = cols())
    colnames(dfMetrics) <- dfMetrics %>% colnames %>% tolower
    dfMetrics <- dfMetrics %>%
        mutate(normpl = posloops/(posloops + negloops), 
               normplweighted = posweighted/(posweighted+negweighted),
               normplund = posloopsund/(posloopsund + negloopsund),
               normplundweighted = posloopsundweighted/(posloopsundweighted + negloopsundweighted),
               normplinfl = inflposloops/(inflposloops + inflnegloops))
    merge(df, dfMetrics, by = "network", all = T)
}

groupCalc1 <- function(topoFile, remNodes = NULL)
{#browser()
    # print(topoFile)
    inflFile <- paste0("Influence/", str_replace(topoFile, ".topo", ".csv"))
    if (!file.exists(inflFile))
        influence_matrix(topoFile)
    if (!file.exists(inflFile))
        return(list(NA, NA))
    inflMat <- read.csv(inflFile, stringsAsFactors = F)
    
    # nodes <- readLines(str_replace(topoFile, ".topo", "_nodes.txt"))
    nodes <- rownames(inflMat) <- colnames(inflMat)[-1]
    df <- as.matrix(inflMat[,-1])
    if (!is.null(remNodes))
    {
        ids <- which(nodes %in% remNodes)
        nodes <- nodes[-ids]
        df <- df[-ids, -ids]
    }
    df1 <- apply(df, 2, function(x){
        ifelse(x > 0, 1, -1)
    })
    df1 <- cbind(df1, t(df1))
    hc <- hclust(dist(df1))
    clust <- cutree(hc, 2)
    g1 <- nodes[clust == 1] %>% sort
    g2 <- nodes[clust == 2] %>% sort
    if(g1[length(g1)]> g2[1])
    {
        g0 <- g1
        g1 <- g2
        g2 <- g0
    }
    l <- list(g1, g2)
    mirdetect <- c(sum(str_detect(g1 %>% str_to_upper, "MIR")), 
                   sum(str_detect(g2 %>% str_to_upper, "MIR")))
    egroup <- which.max(mirdetect)
    mgroup <- ifelse(egroup == 1, 2, 1)
    names(l)[c(egroup, mgroup)] <- c("E", "M")
    return(list(l, list(nodes)))
}
groupCalc1 <- cmpfun(groupCalc1)


AvgStdCalc <- function(df, key)
{
    df$Avg <- df %>% select(contains(key)) %>% apply(1, mean)
    df$Std <- df %>% select(contains(key)) %>% apply(1, sd)
    df %>% select(-contains(key))
}

xlsxTocsv <- function(fl)
{
    df <- readxl::read_xlsx(fl)
    fl <- str_replace(fl, "xslx", "csv")
    write_csv(df, fl, quote = "none")
}

NodeEdges <- function(topoFile, topoPath)
{
    df <- read.delim(paste0(topoPath, "/",topoFile), sep = "")
    nodes <- c(df$Source, df$Target) %>% unique %>% length
    edges <- nrow(df)
    return(paste0(nodes, "N ", edges, "E"))
}

bmCoeff <- function(vec)
{
    vec <- na.omit(vec)
    s <- sd(vec)
    m <- mean(vec)
    n <- length(vec)
    sk <- sum((vec-mean(vec))^3)/((length(vec)-1)*sd(vec)^3)
    ku <- sum((vec-mean(vec))^4)/((length(vec)-1)*sd(vec)^4) -3
    (sk^2 + 1)/(ku + 3*((n-1)^2)/((n-2)*(n-3)))
}


cfgAnalysis <- function(net)
{
    cfgFile <- paste0(net, ".cfg")
    cfgDat <- readLines(cfgFile)
    numEdges <- cfgDat[30] %>% str_extract("\\d+") %>% as.integer
    numNodes <- cfgDat[31] %>% str_extract("\\d+") %>% as.integer
    nodes <- cfgDat[32:(31+numNodes)] %>% str_remove("^\\d+?\t") %>%
        str_trim %>% str_replace_all(regex("\\W+"), "")
    topoDat <- cfgDat[(32+numNodes):length(cfgDat)] %>% str_split("\t") %>%
        reduce(rbind.data.frame) %>% set_names(c("ID", "S", "Tar", "Type")) %>%
        mutate(Type = ifelse(Type == "2", -1, 1), S = as.integer(S), 
               Tar = as.integer(Tar))
    if (!file.exists(paste0(net, ".topo")))
    {
        topoDf <- topoDat %>% select(-ID) %>% mutate(S = nodes[S], Tar = nodes[Tar]) %>%
            set_names(c("Source", "Target", "Type")) %>%
            mutate(Type = ifelse(Type == -1, 2, 1))
        write_delim(topoDf, paste0(net, ".topo"), delim = " ", quote = "none")
    }
    
    return(list(nodes, topoDat))
}

frustCalcRAC <- function(state, topoDat)
{
    sum(state[topoDat$S]*state[topoDat$Tar] == topoDat$Type)/nrow(topoDat)
}
frustCalcRAC <- frustCalcRAC %>% cmpfun


discretize <- function(net)
{
    setwd(paste0(RACIPE_WT, "/", net))
    ls <- cfgAnalysis(net)
    nodes <- ls[[1]]
    topoDat <- ls[[2]]
    solutionDf <- paste0(net, "_solution.dat") %>% read_delim(delim = "\t") %>%
        set_names(c("ParIndex", "nStates", "Count", nodes)) %>% 
         mutate(Count = Count/max(Count)) %>%
        select(all_of(nodes), Count)
    dots <- lapply(nodes, as.symbol)
    states <- solutionDf %>% select(-Count) %>% sapply(function(x){
        y <- (x-mean(x))/sd(x)
        ifelse(y>0, 1, -1)
    }) %>% data.frame %>% set_names(nodes) %>% mutate(Count = solutionDf$Count) %>% 
        group_by(across(nodes)) %>% summarise(Frequency = sum(Count)) 
    states$Frequency <- states$Frequency/sum(states$Frequency)
    states$Frustration <- states %>% select(all_of(nodes)) %>% 
                   apply(1, function(x){frustCalcRAC(x, topoDat)})
    statesDf <- states %>%
        unite(State, all_of(nodes), sep = "") %>% 
        mutate(State = paste0("'",State %>% str_replace_all("-1", "0"),"'"))
    write_csv(statesDf, paste0(net, "_discreteStates.csv"), quote = "none")
    print(net)
}

getEMSONodes <- function(net)
{
    wd <- getwd()
    topoFile <- paste0(net, ".topo")
    if (!file.exists(topoFile))
    {
        net1 <- str_remove(net, "_rand.*")
        setwd(paste0(randRaw, "/", net1))
    }
    
    
    ls <- topo_to_int_mat(topoFile)
    intMat <- ls[[1]]
    nodes <- ls[[2]] %>% str_replace_all(regex("\\W+"), "")
    colnames(intMat) <- rownames(intMat) <- nodes
    signal <- which(apply(intMat, 2, function(x){all(x==0)}))
    output <- which(apply(intMat, 1, function(x){all(x==0)}))
    secondary_signal <- which(apply(intMat, 2, function(x){
        if (length(signal) !=0)
            all(x[-signal] == 0)
        else
            F
    }))
    secondary_output <- which(apply(intMat, 1, function(x){
        if (length(output) != 0)
            all(x[-output] == 0)
        else
            F
    }))
    ls <- groupCalc1(topoFile)
    groups <- ls[[1]] %>% unlist
    sigs <- unique(nodes[c(signal, secondary_signal)])
    names(sigs) <- paste0("S", 1:length(sigs))
    outs <- unique(nodes[c(output, secondary_output)])
    if (length(outs) > 0)
        names(outs) <- paste0("O", 1:length(outs))
    nodes <- unique(c(groups, sigs, outs))
    nodesInGroups <- which(nodes %in% groups)
    dummy <- sapply(1:length(nodes), function(x){
        n <- nodes[x]
        if(n %in% groups)
            names(nodes)[x] <<- names(groups[groups == nodes[x]])
        if(n %in% sigs)
            names(nodes)[x] <<- names(sigs[sigs == nodes[x]])
        if(n %in% outs)
            names(nodes)[x] <<- names(outs[outs == nodes[x]])
    })
    
    setwd(wd)
    return(nodes)
}

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



multiFactorCorrelation <- function(df,x, y,z, label = T, method = "pearson")
{
    facts <- unique(df[[x]])
    s <- sapply(facts, function(f){
        d <- df[df[[x]] == f, ] %>% select(z,y)
        d <- cor.test(d[[y]], d[[z]], method = method)
        if(label)
            paste0("\u03c1: ",round(d$estimate,2), ifelse(d$p.value < 0.05, "*", ""))
        else
            c(d$estimate, d$p.value)
    })
    if (label)
    {
        names(s) <- facts
    }
    else
    {
        s <- s %>% t %>% data.frame %>% set_names(c("Correlation", "pValue")) %>%
            mutate(Factors = facts)
    }
    s
    
}

multiFactorCorrelationAnova <- function(df,x, y,z, label = T, method = "pearson")
{
    facts <- unique(df[[x]])
    s <- sapply(facts, function(f){
        d <- df[df[[x]] == f, ] %>% select(z,y)
        d <- cor(d[[y]], d[[z]], method = method)
        d1 <- df %>% select(all_of(c(x,y,z))) %>% set_names(c("X", "Y", "Z"))
        p <- 
        if(label)
            paste0("\u03c1: ",round(d$estimate,2), ifelse(d$p.value < 0.05, "*", ""))
        else
            c(d$estimate, d$p.value)
    })
    if (label)
    {
        names(s) <- facts
    }
    else
    {
        s <- s %>% t %>% data.frame %>% set_names(c("Correlation", "pValue")) %>%
            mutate(Factors = facts)
    }
    s
    
}


