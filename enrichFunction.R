library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

f_T_CD_CE <- read.csv("D:/R program/proteomics/20210129 Tunicamycin MG132/proteinGroups/fullList/f_T_CD_CE.csv",
                      header = TRUE,
                      stringsAsFactors = FALSE)

setwd("D:/R program/proteomics/20210129 Tunicamycin MG132/proteinGroups")

mouse <- read.csv("mouseUniprot.csv",
                  header = TRUE,
                  stringsAsFactors = FALSE)

enrich1 <- function(proteinID){
  e1 <- enrichGO(gene = proteinID,
               universe = mouse$Entry,
               OrgDb = org.Mm.eg.db,
               keyType = "UNIPROT",
               ont = "CC",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
  
  e2 <- e1@result %>%
    separate(GeneRatio, into = c("xIn", "xQ"),
           sep = "/", convert = TRUE) %>%
    separate(BgRatio, into = c("bIn", "bQ"),
           sep = "/", convert = TRUE ) %>%
    mutate(foldEnrichment = xIn/xQ/(bIn/bQ),
         mlogP = -log10(pvalue)) %>%
    filter(p.adjust < 0.05 & xIn >= 5) %>%
    filter(foldEnrichment > 1.5) %>%
    dplyr::arrange(dplyr::desc(foldEnrichment)) %>%
    dplyr::select(ID, Description, foldEnrichment, pvalue, Count) %>%
    head(30)
  
  e2
}


f <- enrich1(T_ctrBinder$proteinID)


    