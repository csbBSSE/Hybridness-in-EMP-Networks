loopsDat <- function(net)
{
    # browser()
    filz <- paste0("Loops_", net, ".csv")
    df <- read_csv(filz, col_types = cols()) 
    colnames(df) <- str_replace(colnames(df), "Net Weights-", "rand_")
    df <- df %>% 
        mutate(Cycles = str_remove(Cycles, "\\[") %>% str_remove("\\]") %>% str_remove_all("'")) %>%
        mutate(Edge_count = str_count(Cycles, ",") + 1) %>% 
        mutate(Nature = ifelse(rand_1 > 0, "P", "N"))
    colz <- colnames(df)[str_detect(colnames(df), "rand")]
    loops <- sapply(colz, function(x){
        d <- df %>% select(Edge_count, Nature, all_of(x)) %>%
            set_names("Edge_count", "Nature", "Weight") %>%
            group_by(Nature) %>%
            summarise(Weighted = sum(Weight)) %>%
            arrange(Nature)
        posLoops <- filter(d, Nature == "P") %>% select(-Nature) %>% unlist
        if(is_empty(posLoops))
            posLoops <- c(0)
        negLoops <- filter(d, Nature == "N") %>% select(-Nature) %>% unlist
        if(is_empty(negLoops))
            negLoops <- c(0)
        predFrust <- df %>% select(Edge_count, Nature, all_of(x)) %>%
            set_names("Edge_count", "Nature", "Weight") %>%
            filter(Edge_count <= 6, Edge_count > 1, Nature == "N") %>%
            select(Weight) %>% unlist %>% sum
        smallPFL <- df %>% select(Edge_count, Nature, all_of(x)) %>%
            set_names("Edge_count", "Nature", "Weight") %>%
            filter(Edge_count <= 6, Edge_count > 1, Nature == "P") %>%
            select(Weight) %>% unlist %>% sum
        c(posLoops, negLoops, predFrust, smallPFL)
    }) %>% t %>% data.frame %>% set_names(c( "posLoops", "negLoops", "predFrust",
                                             "smallPFL")) %>%
        mutate(Network = colz)
    loops
}

hiLoopData <- function(net) {
    topoDf <- read.delim(paste0(net, ".topo"), sep = "") %>% 
        mutate(Edge = paste0(Source, "_", Target))
    edges <- topoDf %>% select(Edge) %>% unlist
    nodes <- c(topoDf$Source, topoDf$Target) %>% unique
    filz <- paste0("Loops_", net, ".csv")
    df <- read_csv(filz, col_types = cols()) 
    colnames(df) <- str_replace(colnames(df), "Net Weights-", "rand_")
    
    ## Type1
    loops <- df %>% select(Cycles) %>% unlist %>%
        str_remove("\\[") %>% str_remove("\\]") %>% str_remove_all("'")
    loopEdges <- loops %>% str_split(",") %>% lapply(function(x) {
        x <- c(x, x[1])
        if (length(x) == 2)
        {
            edges <- paste0(x, collapse = "_")
        }
        sapply(2:length(x), function(i) {
            paste0(x[i-1], "_", x[i])
        })
    }) %>% sapply(paste0, collapse = ",")
    hiLoops <- matrix(nrow = 1, ncol = 3)
    hiLoopCount <- 0
    edgeLoops <- lapply(edges, function(x) {
        s <- which(str_detect(loopEdges, x))
        if (length(s)<3) {
            return(s)
        }
        else {
            hiLoops <<- rbind(hiLoops, combn(s, 3) %>% t)
        }
        return(s)
    })
    nodeLoops <- lapply(nodes, function(x) {
        s <- which(str_detect(loops, x))
        if (length(s)<3) {
            return(s)
        }
        else {
            hiLoopCount <<- hiLoopCount + choose(length(s), 3)
        }
        return(s)
    })
}
