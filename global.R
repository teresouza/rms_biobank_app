reqPackages <- c("shiny", "shinyjs", "ggplot2", "ggthemes", "ggsci","tidyr", "dplyr", 
                 "data.table", "shinyWidgets", "plotly", "shinycssloaders", "tibble", "DT")
package.check <- lapply(
  reqPackages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
allData <- readRDS("./data/allData_March2022.RDS")
rna <- allData$rna
menu.data <- allData$menu
fusion <- allData$fusion
snvs <- allData$snvs
cna <- allData$cnas
metadata <- allData$metadata
drugs.cor <- allData$drugs
abvs <- read.delim("./data/RMS_biobank_abvs_March2021.txt", header = T, check.names = F)
