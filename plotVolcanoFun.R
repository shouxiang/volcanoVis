library(tidyverse)
f_T_CD_CE <- read.csv("D:/R program/proteomics/20210129 Tunicamycin MG132/proteinGroups/fullList/f_T_CD_CE.csv",
                      header = TRUE,
                      stringsAsFactors = FALSE)

plotVolcano <- function(df, fcTH = 1.5, pvalueTH = 0.04666594) {

  xMin <- min(df$log2fc) - .2
  xMax <- max(df$log2fc) + .2
  yMax <- max(-log10(df$pvalue)) + .2

  hit <- df %>% filter(log2fc >= log2(fcTH) & pvalue <= pvalueTH)

  ggplot(df, aes(log2fc, log1P)) +
    geom_point(colour = "grey", alpha = .4, size = 3.5) +
    geom_point(data = hit,
               colour = "red", size = 3.5) +
    scale_x_continuous(limits = c(xMin, xMax)) +
    scale_y_continuous(limits = c(NA, yMax)) +
    geom_segment(aes(x = log2(fcTH), y = -log10(pvalueTH),
                     xend = xMax, yend = -log10(pvalueTH)),
               linetype = "dashed") +
    geom_segment(aes(x = log2(fcTH), y = -log10(pvalueTH),
                     xend = log2(fcTH), yend = yMax),
               linetype = "dashed") +
    theme_classic() +
    theme(
      strip.background = element_rect(colour = "white"),
      axis.text.x = element_text(margin = margin(t = 2.5, b = 2.5),
                                 size = 15,
                                 face = "bold",
                                 color = "black"),
      axis.text.y = element_text(margin = margin(r = 2.5, l = 2.5),
                                 size = 15,
                                 face = "bold",
                                 color = "black"),
      axis.title.x = element_text(size = 20,
                                  face = "bold",
                                 color = "black"),
      axis.title.y = element_text(size = 20,
                                  face = "bold",
                                 color = "black"),
      axis.ticks.length = unit( 1.5, "mm"),
      plot.title = element_text(face = 'bold',
                                colour = 'red',
                                size = 26,
                                hjust = 0.5)) +
    labs(x = "log2(FC)", y = "-logP") +
    ggtitle(hit %>% nrow())
}



OmicsXZ::plotVolcano(f_T_CD_CE)
