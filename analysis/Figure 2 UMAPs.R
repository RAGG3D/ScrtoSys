library(WGCNA)
library(reshape2)
library(stringr)
library(tidyverse)
library(tidybulk)
library(org.Hs.eg.db)
library(Hmisc)
library(msigdbr)
library(survival)
library(ComplexHeatmap)
require(vegan)
cancer_list = c("BLCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "LIHC", "LUAD",
                "PRAD", "READ", "STAD", "THCA")
scrto <- read_csv("C:/Rstudio/Chapter 3/data/comprehensive secreted protein gene list.csv")

TCGAun_mat <- map_dfr(list("BLCA", "COAD", "ESCA", "HNSC", "KICH", "KIRC", "LIHC", "LUAD",
                           "PRAD", "READ", "STAD", "THCA"),
                      function(cancer0){
                        x <- read_csv(paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/", cancer0,".csv")) %>%
                          dplyr::filter(symbol %in% unique(scrto$symbol))})

TCGAumap <- TCGAun_mat %>%
  dplyr::select(sample_ct, symbol, raw_count_scaled_adjusted, cancer) %>%
  reduce_dimensions(.element = sample_ct, .feature = symbol, .abundance = raw_count_scaled_adjusted, method="UMAP", top = 2329)

Figure2 <- TCGAumap %>%
  ggplot(aes(x = `UMAP1`, y = `UMAP2`, color = cancer)) + 
  geom_point(size = 0.5)+
  stat_ellipse() +
  theme_bw() + 
  scale_color_ucscgb() + 
  theme(text=element_text(size=20, family="sans"),
        legend.position = "bottom") +
  labs(title = "") + 
  guides(colour = guide_legend(override.aes = list(size=10)),
         fill=guide_legend(nrow=2,byrow=TRUE))

ggsave(plot = Figure2, "C:/Rstudio/Chapter 3/Figures/Figure 2 (secretome T+N).png", 
       height = 15, width = 15)
