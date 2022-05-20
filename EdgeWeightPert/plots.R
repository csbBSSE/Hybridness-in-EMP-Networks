plotFolder <- "D:/Github/Projects/Ongoing/Mubasher/Figures/edgeWeight"
wd <- "D:/Github/Projects/Ongoing/Mubasher/EdgeWeightPert"

setwd(wd)

source("setup.R")
library(grid)
library(gridExtra)

scatter <- function(net)
{
    setwd(wd)
    df <- read_csv(paste0(net, "_all.csv"), col_types = cols())
    colnames(df) <- tolower(colnames(df))
    # dfWT <- df %>% filter(network == "wild") %>% mutate(network = "WT")
    setwd(plotFolder)
    if(!dir.exists(("ScatterPlots")))
        dir.create("ScatterPlots")
    setwd("ScatterPlots")
    if(!dir.exists(net))
        dir.create(net)
    setwd(net)
    # colnames(df) <- tolower(colnames(df))
    metrics <- colnames(df)[-1]
    plots <- sapply(metrics, function(metric){
        # browser()
        # print(metric)
        df[[metric]] <- as.numeric(df[[metric]])
        # dfWT[[metric]] <- as.numeric(dfWT[[metric]])
        gr <- correlGrob(df, metric, "hybridness", method = "spearman")
        
        ggplot(df, aes_string(x = metric, y = "hybridness")) +
            geom_point(size = 2) + 
            annotation_custom(gr) + theme_Publication() + 
            theme(axis.text = element_text(size = rel(1.5))) +
            labs(x = colNKey[metric], y = "Hybridness") 
        ggsave(paste0(net, "_", metric, ".png"), width = 5.5, height = 5)
    })
}

sapply(networks, scatter)


scatterCustom <- function(net, m1, m2, RACIPE = T) {
    setwd(wd)
    print(net)
    m1 <- tolower(m1)
    m2 <- tolower(m2)
    
    df <- read_csv(paste0(net, "_all.csv"))
    colnames(df) <- tolower(colnames(df))
    setwd(plotFolder)
    if(!dir.exists(("ScatterPlots")))
        dir.create("ScatterPlots")
    setwd("ScatterPlots")
    if(!dir.exists(net))
        dir.create(net)
    setwd(net)
    df[[m1]] <- as.numeric(df[[m1]])
    df[[m2]] <- as.numeric(df[[m2]])
    gr <- correlGrob(df, m1, m2, method = "spearman")
    
    ggplot(df, aes_string(x = m1, y = m2)) +
        geom_point(size = 2) + 
        annotation_custom(gr) + theme_Publication() + 
        theme(axis.text = element_text(size = rel(1.5))) +
        labs(x = colNKey[m1], y = colNKey[m2]) 
    ggsave(paste0(net, "_", m1, "_", m2, ".png"), width = 5.5, height = 5)
}

frustMets <- c("minFrust", "maxFrust", "meanFrust")
sapply(frustMets, function(x){
    sapply(networks, scatterCustom, m1 = "predFrust", m2 = x)
})

heatMapMultDat <- function(net, metrics, racipe = T, m2 = "hybridness")
{
    print(net)
    df <- paste0(net, "_all.csv") %>% read_csv
    df <- df %>% 
        set_names(colnames(.) %>% tolower) %>% 
        select(all_of(metrics %>% tolower), all_of(m2 %>% tolower)) %>% drop_na %>% 
        mutate_all(as.numeric)
    corDf <- sapply(metrics, function(metric){
        print(metric)
        corRAC <- cor.test(df[[tolower(m2)]], df[[tolower(metric)]], method= "spearman")
        c(corRAC$estimate, corRAC$p.value)
    }) %>% t %>% data.frame %>% set_names(c("Correlation", "pValue")) %>% 
        mutate(Metrics = colNKey[tolower(metrics)], Network = netNames[net])
    corDf
}

heatMapMult <- function(nets, metrics, nam, racipe = T, m2 = "hybridness")
{
    setwd(wd)
    df <- lapply(nets, heatMapMultDat, metrics = metrics, racipe = racipe, m2 = m2) %>% 
        reduce(rbind.data.frame) %>%
        mutate(Significance = ifelse(pValue < 0.05, "", "X")) %>%
        mutate(Metrics = stringBreak(Metrics)) %>%
        mutate(Network = factor(Network, levels = netNames[nets] %>% unname),
               Metrics = factor(Metrics, levels = rev(colNKey[tolower(metrics)] %>% stringBreak)))
    ggplot(df, aes(x = Network, y = Metrics, fill = Correlation)) +
        geom_tile() + 
        geom_text(aes(label = Significance)) + theme_Publication() + labs(y = "") + 
        scale_fill_gradient2(low = "blue", high = "red", limits = c(-1,1)) +
        theme(legend.position = "top", legend.key.width = unit(0.8, "cm"),
              axis.text.x = element_text(angle= 60, hjust = 1, vjust = 1))
    setwd(plotFolder)
    if(!dir.exists("Heatmaps"))
        dir.create("Heatmaps")
    setwd("Heatmaps")
    ggsave(paste0(nam, ".png"), width = 10, height = 2.5 + length(metrics)*0.8)
}

heatMapMult(networks, c("jsd", "jmetric","meanfrust", "predfrust", "posloops", "negloops", 
                        "inflposloops", "inflnegloops"), "all")
