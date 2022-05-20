library(miscFuncs)
library(tidyverse)

topoFiles <- list.files(".", ".topo")

randNetGen(topoFiles)


source("dataGen.R")

source("networkMetric.R")
