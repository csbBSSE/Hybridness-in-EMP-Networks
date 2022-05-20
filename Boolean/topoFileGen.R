edelFuncs <- function(net)
{
    topoFile <- paste0(net, ".topo")
    if(!dir.exists(net))
        dir.create(net)
    setwd(net)
    file.copy(paste0("../", net, ".topo"), "./wild.topo")
    topoDf <- read.delim("wild.topo",sep = "")
    edges <- topoDf %>% mutate(edgeNature = ifelse(Type == 1, 2, 1)) %>% 
        unite("Edge", Source, Target, Type, sep ="_") %>%
        mutate(ID = 1:nrow(.))
    delTopos <- apply(edges, 1, function(x){
        id <- x[3] %>% as.integer
        eN <- x[2]
        write_delim(topoDf[-id, ], paste0(x[1], "_0.topo"), delim = " ")
        t <- topoDf
        t[id, 3] <- eN %>% as.integer
        write_delim(t, paste0(x[1], "_", eN, ".topo"), delim = " ")
    })
    setwd("..")
}
sapply(c("grhl", "ovol", "oct", "nrf"), edelFuncs)
sapply(c("SIL", "SIL2", "EMT_RACIPE", "EMT_RACIPE2", "EMT_MET_reduced"), edelFuncs)
