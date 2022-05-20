wd <- "D:/Github/Projects/Ongoing/Mubasher/RandomNets"
plotFolder <- "D:/Github/Projects/Ongoing/Mubasher/Figures/Rand"
source("D:/Github/Projects/Ongoing/Mubasher/Figures/figCore.R")
source("D:/Github/Projects/Ongoing/Mubasher/Figures/utils.R")

setwd(netMetrics)
compiledData <- lapply(nets, function(net) {
    read.csv(paste0(net, "_randNetMetrics.csv")) %>% 
        mutate(Net = net)
}) %>% reduce(rbind.data.frame) %>%
    set_names(tolower(colnames(.))) %>%
    mutate(normpl = posloops/(posloops + negloops),
           normplweighted = posweighted/(posweighted + negweighted))
WTData <- compiledData %>% filter(network == "wild")

metrics <- c("predfrust", "normpl", "normplweighted", "inflposloops", "inconsistency", "type1", "missa")

percentiles <- sapply(metrics, function(x){
    df <- compiledData %>% select(all_of(x), net) %>% 
        set_names(c("Metric", "Network")) %>% 
        mutate(Network = factor(netKey[Network], levels = netKey[nets]))
    dfWT <- WTData %>% select(all_of(x), net) %>% 
        set_names(c("Metric", "Network")) %>% 
        mutate(Network = factor(netKey[Network], levels = netKey[nets]))
    ggplot(df, aes(x = Network, y = Metric)) + geom_boxplot() + 
        geom_point(data = dfWT, color = "red", size = 1, shape = 3, stroke = 3) +
        labs(y = colNKey[tolower(x)]) + theme_Publication() +
        theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
    ggsave(paste0(plotFolder, "/", x, "_violins.png"), width = 7, height = 4)
    df %>% split(df$Network) %>% sapply(function(d) {
        sum(d$Metric < dfWT$Metric[dfWT$Network == as.character(unique(d$Network))])/nrow(d)
    })
}) %>% data.frame %>% set_names(tolower(metrics)) %>%
    mutate(Network = factor(netKey, levels = netKey)) %>%
    gather(key = "Metric", value = "Value", -Network) %>% 
    mutate(Metric = stringBreak(colNKey[Metric])) %>% 
    mutate(Metric = factor(Metric, levels = colNKey[metrics] %>% stringBreak %>% rev))

ggplot(percentiles, aes(x = Network, y = Metric, fill = Value)) +
    geom_tile() + scale_fill_viridis_c() + theme_Publication() + 
    theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
          legend.position = "top", legend.key.width = unit(0.8, "cm")) +
    labs(x = "", y = "", fill = "Percentile")
ggsave(paste0(plotFolder, "/percentiles.png"), width = 13, height = 2.5 + length(metrics)*0.6)

    
    