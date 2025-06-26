library(tidyverse)
library(tidybulk)
library(org.Hs.eg.db)
library(Hmisc)
library(msigdbr)
library(survival)
library(survminer)
library(ComplexHeatmap)

CIBERSORT <- function(x, c1, c2, c3, ref, meth, act) {
  as_tibble(setDT(x)[, lapply(.SD, sum), by = c(c1, c2), .SDcols = c3]) %>%
    deconvolve_cellularity(
      .sample = !!sym(c1),
      .transcript = !!sym(c2),
      .abundance = !!sym(c3),
      reference = ref,
      method = meth,
      action = act
    )
}

ref = readRDS("C:/Rstudio/Chapter 3/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("C:/Rstudio/Chapter 3/data/cellname.csv")))$`Cell Type`
scrto <- read_csv("C:/Rstudio/Chapter 3/data/comprehensive secreted protein gene list.csv")
h_df = msigdbr(species = "Homo sapiens", category = "H")
cancerlist = list("KICH", "READ", 
                  "LIHC", "BLCA", "STAD", "COAD", 
                  "PRAD", "THCA", "HNSC", "KIRC", "LUAD", 
                  "ESCA")
cancerlist = list("LUAD", 
                  "ESCA")
####*functions*####


ggsurvp_cat_item <- function(x, p, nrow, palette, xlab, ylab, title, size0){
  ggsurvplot(
    survfit(
      Surv(total_living_days, vital_status) ~ item,
      data = x
    ),
    data = x,
    facet.by = c("cat"),
    scales = "free_x",
    conf.int = T,
    risk.table = F,
    conf.int.alpha = 0.15,
    pval = p,
    nrow = nrow,
    legend.title = " ",
    short.panel.labs = T,
    palette = palette
  ) + theme_classic() +
    theme(text=element_text(size=size0, family="sans"),
          legend.position = "bottom") +
    labs(x = xlab, y = ylab, title = title) +
    guides(linetype = FALSE)
}

three_heatmaps_of_pth_scrto_TN <- function(cancerlist){
  map(cancerlist, function(cancer0){
    dir.create(paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/"))
    #### normal hallmark score and scrtos 
    x2 <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene normal pathway+scrto corr.csv")) %>%
      filter(!grepl("HALLMARK_", X2)) %>%
      filter(p.value <= 0.05) %>%
      dplyr::select(X1, X2, value) %>% 
      pivot_wider(names_from = "X2", values_from = value) %>% 
      replace(is.na(.), 0) %>%
      pivot_longer(-X1, names_to = "Secretome_Normal", values_to = "value") 
    
    
    p2 <- x2 %>%
      tidybulk::rename("Pathway" = "X1") %>%
      tidyHeatmap::heatmap(
        .column = Secretome_Normal,
        .row = Pathway,
        .value = value,   
        palette_value = circlize::colorRamp2(
          seq(min(x2$value), max(x2$value), length.out = 5), 
          rev(RColorBrewer::brewer.pal(5, "RdBu"))
        ),
        # Arguments passed to ComplexHeatmap 
        clustering_distance_columns = "manhattan",
        clustering_method_columns = "ward.D",
        clustering_distance_rows = "manhattan",
        clustering_method_rows = "ward.D",
        column_title = paste0(cancer0, " Normal")
      ) 
    
    #### tumor hallmark score and scrtos
    x1 <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene pathway+scrto corr.csv")) %>%
      filter(!grepl("HALLMARK_", X2)) %>%
      filter(p.value <= 0.05) %>%
      dplyr::select(X1, X2, value) %>% 
      pivot_wider(names_from = "X2", values_from = value) %>% 
      replace(is.na(.), 0) %>%
      pivot_longer(-X1, names_to = "Secretome_Tumor", values_to = "value") 
    
    
    p1 <- x1 %>%
      tidybulk::rename("Pathway" = "X1") %>%
      tidyHeatmap::heatmap(
        .column = Secretome_Tumor,
        .row = Pathway,
        .value = value,   
        palette_value = circlize::colorRamp2(
          seq(min(x1$value), max(x1$value), length.out = 5), 
          rev(RColorBrewer::brewer.pal(5, "RdBu"))
        ),
        
        clustering_distance_columns = "manhattan",
        clustering_method_columns = "ward.D",
        clustering_distance_rows = "manhattan",
        clustering_method_rows = "ward.D",
        column_title = paste0(cancer0, " Tumour")
      ) 
    
    #### tumor-normal hallmark score and scrtos
    x3 <- x1 %>% mutate(tumor = value, Secretome = Secretome_Tumor) %>%
      dplyr::select(X1, Secretome, tumor) %>%
      left_join(x2 %>% mutate(normal = value, Secretome = Secretome_Normal) %>%
                  dplyr::select(X1, Secretome, normal)) %>%
      mutate(value = tumor - normal, Secretome_Complement = Secretome)%>%
      dplyr::select(X1, Secretome_Complement, value) %>% 
      replace(is.na(.), 0)
    
    p3 <- x3 %>%
      tidybulk::rename("Pathway" = "X1") %>%
      tidyHeatmap::heatmap(
        .column = Secretome_Complement,
        .row = Pathway,
        .value = value,   
        palette_value = circlize::colorRamp2(
          seq(min(x3$value), max(x3$value),  length.out = 5), 
          rev(RColorBrewer::brewer.pal(5, "RdBu"))
        ),
        
        clustering_distance_columns = "manhattan",
        clustering_method_columns = "ward.D",
        clustering_distance_rows = "manhattan",
        clustering_method_rows = "ward.D",
        column_title = paste0(cancer0, " T-N")
      ) 
    
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/heatmap Tumor Normal.pdf"), height = 8, width = 10)
    draw(p1+p2)
    dev.off()
    
    
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/heatmap Tumor (Tumor-Normal).pdf"), height = 8, width = 10)
    draw(p1+p3)
    dev.off()
    
  })
}

clinical_combine <- function(x, cancer) {
  x %>% 
    inner_join(read.csv(paste0("C:/Rstudio/Chapter 3/data/TCGA clinical/clinical_", tolower(cancer), ".csv")), by = c("sample" = "bcr_patient_barcode"))%>%
    mutate(total_living_days = as.numeric(as.character(days_to_death)), age = -as.numeric(as.character(days_to_birth))/365) %>% 
    mutate(na = is.na(total_living_days)) %>%
    mutate(total_living_days = ifelse(na == "TRUE", as.numeric(as.character(last_contact_days_to)), total_living_days)) %>%
    mutate(vital_status = ifelse(vital_status == "Dead", 1, 0))}



####*Running*####
three_heatmaps_of_pth_scrto_TN(cancerlist)


map(cancerlist, function(cancer0){
  pthscore <- map_dfr(list.files(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/", cancer0, "/")), function(name0){
    read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/", cancer0, "/", name0))%>%
      mutate(sample_ct = gsub(".csv", "", name0))
  }) %>%
    dplyr::select(sample_ct, ID,enrichmentScore, pvalue) 
  
  x <- pthscore %>%
    mutate(sample = gsub("-Tumor", "", sample_ct)) %>%
    clinical_combine(cancer0) 
  
  coxpth <- map_dfr(as.list(unique(x$ID)), function(pth0){
    a <- summary(coxph(Surv(total_living_days, vital_status) ~
                         enrichmentScore,
                       data=x %>% filter(ID == pth0) ))
    
    data.frame(ID = pth0,
               HR = a$coefficients[2],
               P = a$coefficients[5]) 
    
  })

    
    IDs = (coxpth %>% filter(P <= 0.05))$ID 
    
    ####***Tumor (Heatmap + KM)***#####
    pth_surv <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene pathway+scrto corr.csv")) %>%
      filter(p.value <= 0.05) %>%
      replace(is.na(.), 0) %>%
      filter(X1 %in% IDs) %>%
      tidybulk::rename("ID" = "X1", "Scrto_pth Tumour corr" = "value", "corr p" = p.value) %>%
        inner_join(coxpth) %>%
      tidybulk::rename("Hallmark HR" = "HR") %>%
      write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/Pthscore/", cancer0, " pthscore_surv.csv"), row.names = F)
    
    
  })

####*Pathways and secretomes in normal tissues*####
map(cancerlist, function(cancer0){
  pthscore <- map_dfr(list.files(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/", cancer0, "/")), function(name0){
    read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/", cancer0, "/", name0))%>%
      mutate(sample_ct = gsub(".csv", "", name0))
  }) %>%
    dplyr::select(sample_ct, ID,enrichmentScore, pvalue) 
  
  x <- pthscore %>%
    mutate(sample = gsub("-Normal", "", sample_ct)) %>%
    clinical_combine(cancer0) 
  
  coxpth <- map_dfr(as.list(unique(x$ID)), function(pth0){
    a <- summary(coxph(Surv(total_living_days, vital_status) ~
                         enrichmentScore,
                       data=x %>% filter(ID == pth0) ))
    
    data.frame(ID = pth0,
               HR = a$coefficients[2],
               P = a$coefficients[5]) 
    
  })
  
  
  IDs = (coxpth %>% filter(P <= 0.05))$ID 
  
  ####***Tumor (Heatmap + KM)***#####
  pth_surv <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene pathway+scrto corr.csv")) %>%
    filter(p.value <= 0.05) %>%
    replace(is.na(.), 0) %>%
    filter(X1 %in% IDs) %>%
    tidybulk::rename("ID" = "X1", "Scrto_pth Normal corr" = "value", "corr p" = p.value) %>%
    inner_join(coxpth) %>%
    tidybulk::rename("Hallmark HR" = "HR") %>%
    mutate(Tissue = "Normal") %>%
    write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/Pthscore/", cancer0, " pthscore_surv.csv"), row.names = F)
  
  
})

######**Heatmap + KM**#####
map(cancerlist, function(cancer0){
  pthgenes <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/Pthscore/", cancer0, " pthscore_surv.csv")) %>%
  mutate(value = `Scrto_pth Tumour corr`, HR = `Hallmark HR`) %>%
    mutate(ID = gsub("HALLMARK_", "", ID)) %>%
    filter(X2 %in% (as.data.frame(hubs) %>% filter(hubs %in% scrto$symbol))$hubs)
  
  mods <- as.data.frame(hubs) %>% rownames_to_column("Color") %>%
    mutate(yii = "ME") %>%
    unite("Var1", c("yii", "Color"), sep = "")
    
    
    
 x1 <- pthgenes %>% mutate(hubs = X2) %>% dplyr::select(ID, hubs, value, HR) %>%
   inner_join(mods) %>%
   unite("hubs", c("Var1", "hubs"), sep = "_") 
 
 x2 <- x1 %>%
   rbind(mod_trait %>%  
           tidybulk::rename("ID" = "Var2" ) %>%
           inner_join(mods) %>%
           filter(hubs %in% pthgenes$X2) %>%
           unite("hubs", c("Var1", "hubs"), sep = "_") %>%
           tidybulk::rename("value" = "Cor") %>%
           mutate(HR = 1) %>%
           dplyr::select(hubs, ID, value, HR)) %>%
   filter(!ID %in% x1$ID)

 
 p1 <- x1 %>% dplyr::select(ID, hubs, value) %>%
   pivot_wider(names_from = ID, values_from = value) %>%
   mutate_all(~replace(., is.na(.), 0)) %>%
   pivot_longer(-hubs, names_to = "ID", values_to = "value") %>%
   inner_join(unique(x1 %>% dplyr::select(ID, HR))) %>%
   tidybulk::rename("Pathway" = "ID", "Secretome" = "hubs") %>%
  tidyHeatmap::heatmap(Pathway, Secretome, value, scale = "none",
                       palette_value = circlize::colorRamp2(
                         seq(min(x1$value), max(x1$value), length.out = 11), 
                         RColorBrewer::brewer.pal(11, "BrBG")
                       )) %>%
  tidyHeatmap::add_tile(
    HR, 
palette = circlize::colorRamp2(
      seq(0, 2, length.out = 3), 
      c("Green", "White", "Red")
    )
  )%>%
    tidyHeatmap::as_ComplexHeatmap() %>%
    ComplexHeatmap::draw(heatmap_legend_side = "left"   ) 
 pdf(file=paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/", cancer0, " Heatmap HR tumor pthscore_surv.pdf"), height = 5, width = 10)
 draw(p1)
 dev.off()
 
 p2 <- x2 %>% dplyr::select(ID, hubs, value) %>%
   pivot_wider(names_from = ID, values_from = value) %>%
   mutate_all(~replace(., is.na(.), 0)) %>%
   pivot_longer(-hubs, names_to = "ID", values_to = "value") %>%
   inner_join(unique(x2 %>% dplyr::select(ID, HR))) %>%
   tidybulk::rename("Pathway" = "ID", "Secretome" = "hubs") %>%
   tidyHeatmap::heatmap(Pathway, Secretome, value, scale = "none",
                        palette_value = circlize::colorRamp2(
                          seq(min(x1$value), max(x1$value), length.out = 11), 
                          RColorBrewer::brewer.pal(11, "BrBG")
                        )) %>%
   tidyHeatmap::as_ComplexHeatmap() %>%
   ComplexHeatmap::draw(heatmap_legend_side = "left"   ) 
 pdf(file=paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/", cancer0, " Heatmap HR tumor clinical.pdf"), height = 5, width = 10)
 draw(p2)
 dev.off()
})

map(cancerlist, function(cancer0){
  gene <- read_csv(paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/", cancer0, ".csv")) 
  deg <- gene %>%
    test_differential_abundance(~tissue, 
                                .sample = sample_ct,
                                .abundance = raw_count_scaled_adjusted,
                                .transcript = symbol,
                                action = "get",
                                test_above_log2_fold_change = T) %>%
    write.csv(paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/", cancer0, " DEG.csv"), row.names = F)
  
  x <- read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/Pthscore/", cancer0, " pthscore_surv.csv")) %>%
    mutate(value = `Scrto_pth Tumour corr`, HR = `Hallmark HR`) %>%
    mutate(ID = gsub("HALLMARK_", "", ID)) %>%
    filter(X2 %in% uniqeu(deg$symbol))
  p <- x %>% dplyr::select(ID, X2, value) %>%
    pivot_wider(names_from = ID, values_from = value) %>%
    mutate_all(~replace(., is.na(.), 0)) %>%
    pivot_longer(-X2, names_to = "ID", values_to = "value") %>%
    inner_join(unique(x %>% dplyr::select(ID, HR))) %>%
    tidybulk::rename("Pathway" = "ID", "Secretome" = "X2") %>%
    tidyHeatmap::heatmap(Pathway, Secretome, value,    scale = "row",
                         palette_value = circlize::colorRamp2(
                           seq(min(x$value), max(x$value), length.out = 11), 
                           RColorBrewer::brewer.pal(11, "RdBu")
                         )) %>%
    tidyHeatmap::add_tile(
      HR, 
      palette = circlize::colorRamp2(
        seq(0, 2, length.out = 3), 
        c("Green", "White", "Red")
      )
    )%>%
    tidyHeatmap::as_ComplexHeatmap() %>%
    ComplexHeatmap::draw(heatmap_legend_side = "left"   ) 
  pdf(file=paste0("C:/Rstudio/Chapter 3/output/Figures/", cancer0, "/", cancer0, " DEGs Heatmap HR tumor pthscore_surv.pdf"), height = 5, width = 10)
  draw(p)
  dev.off()
})

cancer0 = "BLCA"
gene <- read_csv(paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/", cancer0, ".csv")) %>%
  filter(tissue == "Tumor") 

cell <- read_csv(paste0("C:/Rstudio/Chapter 3/data/", cancer0, ".csv")) %>%
  unite("sample_ct", c("sample", "tissue"), sep = "-") 



ggplot(cell, aes(x=celltype, y=log2(abundance+1))) + 
  geom_boxplot(color="red", fill="orange", alpha=0.2) +
  theme_classic()
trait <- paste0("C:/Rstudio/Chapter 3/data/TCGA clinical/clinical_", tolower(cancer0)," WGCNA.csv")
MEs_col <- read_csv(paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " Module scores.csv"))
x <- cell %>%
  inner_join(MEs_col %>%
  rownames_to_column("sample") %>%
  dplyr::select(sample, MEgreen) %>%
  mutate(sample_ct = read_csv(trait)$sample) %>%
    dplyr::select(-sample))

a <- reshape::melt(cor(x %>% dplyr::select(-sample_ct))[-c(1:24), c(1:24)])
a.p <- rcorr(as.matrix(x %>% dplyr::select(-sample_ct)))
a <- a %>%
  rownames_to_column("X1") %>%
  inner_join(as.data.frame(a.p$P[-c(1:24), c(1:24)]) %>% 
                   mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:24), c(1:24)]))) %>% 
                   gather(X2, p.value, -X1)) %>%
  dplyr::select(-X2)

yellow_cells <- (a %>% filter(p.value <= 0.05))$X1
split_median <- function(x, cat, item){
  x %>%
    group_by(!!as.name(cat)) %>%
    mutate(item = factor(Hmisc::cut2(!!as.name(item), g = 2), labels = c(1:nlevels(Hmisc::cut2(!!as.name(item), g = 2))))) %>%
    ungroup() %>%
    mutate(item = gsub("1", "L", item)) %>%
    mutate(item = gsub("2", "H", item))
}

z <- cell %>% pivot_longer(-sample_ct, values_to = "portion", names_to = "celltype") %>%
  mutate(sample = gsub("-Tumor", "", sample_ct)) %>%
  split_median("celltype", "portion") %>%
  clinical_combine("BLCA") %>%
  mutate(cat = celltype) %>%
  mutate(item = fct_relevel( item, "L", "H")) %>%
  as.data.frame()
ggsurvp_cat_item(z%>% filter(celltype %in% c("IL2NK","mDC")), T, 3, "jco", "", "", cancer0, 16)

x <- gene %>% filter(symbol == "COL6A3") %>%
  mutate(sample = gsub("-Tumor", "", sample_ct)) %>%
  split_median("symbol", "raw_count_scaled_adjusted") %>%
  clinical_combine("BLCA") %>%
  mutate(cat = symbol) %>%
  mutate(item = fct_relevel( item, "L", "H")) %>%
  as.data.frame()
coxph(Surv(total_living_days, vital_status) ~ raw_count_scaled_adjusted, x)  

p0 <- ggsurvp_cat_item(x, T, 1, "jco", "", "", cancer0, 16)

y <- x %>% inner_join(cell %>% pivot_longer(-sample_ct, names_to = "celltype", values_to = "fraction") %>%
                        filter(celltype %in% yellow_cells) %>%
                   split_median("celltype", "fraction") %>%
                   tidybulk::rename("cell" = "item")) %>%
  unite("item", c("item", "cell"), sep = "/") %>%
  unite("cat", c("cat", "celltype"), sep = "/")  %>%
  mutate(item = fct_relevel(item, "L/L", "L/H", "H/L", "H/H")) %>%
  as.data.frame()
ggsurvp_cat_item(y%>%
                   filter(item %in% c("L/L", "L/H"))%>%
                   filter(cat %in% c("COL6A3/IL2NK", "COL6A3/mDC")), T, 3, c("#0073C2FF","#EFC000FF"), "", "", cancer0, 16)
p1 <- ggsurvp_cat_item(y%>%
                         filter(item %in% c("H/L", "H/H")), T, 3, c("#868686FF","#A73030FF"), "", "", cancer0, 16)
p1a <- ggsurvp_cat_item(y %>%
                          filter(item %in% c("H/L", "H/H")) %>%
                          filter(cat %in% c("COL6A3/IL2NK", "COL6A3/mDC")), T, 1, c("#868686FF","#A73030FF"), "", "", cancer0, 16)

pthscore <- map_dfr(list.files(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/", cancer0, "/")), function(name0){
  read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/", cancer0, "/", name0))%>%
    mutate(sample_ct = gsub(".csv", "", name0))
}) %>%
  dplyr::select(sample_ct, ID,enrichmentScore, pvalue) 

z <- cell %>% dplyr::select(sample_ct, Epi, Fibro, IL2NK, mDC, ReNK, SPANK, 
                            contains("GD T"), contains("Naive CD8 T")) %>%
  inner_join(pthscore %>%
               filter(ID %in% c("HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE")) %>%
               dplyr::select(-pvalue) %>%
               pivot_wider(names_from = "ID", values_from = "enrichmentScore")) 
a <- reshape::melt(cor(z %>% dplyr::select(-sample_ct))[-c(1:8), c(1:8)])
p2 <- as_tibble(a ) %>%
  ggplot(aes(x = X2, y = X1, fill = value)) +
  geom_tile(stat = "identity") +
  theme_classic() +
  theme(text=element_text(size=16, family="sans"),
        panel.spacing.x = unit(0,"line"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "", title = "") + 
  scale_fill_gradient2(low = "#006699", high = "#b30000", mid = "white",
                       midpoint = 0, 
                       name = "Correlation",
                       limits=range(-0.4, 0.5),
                       breaks = seq(-0.4, 0.4, 0.5))

pthcore_HR <- read_csv("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/Pthscore/BLCA pthscore_surv.csv") %>%
  inner_join(read_csv("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/BLCA whole gene pathway+scrto corr.csv") %>%
               tidybulk::rename("ID" = "X1"))

x <- gene %>%
  filter(symbol %in% c("COL6A3", "CD226", "CD244", "CRTAM", "KIR2DL4", "KLRC1", 
                       "KLRC2", "KLRC3", "KLRC4", "KLRK1", "NCR1", "NCR2", "NCR3",
                       "CD80", "CD83","ITGAX")) %>%
  dplyr::select(sample_ct, symbol, raw_count_scaled_adjusted) %>%
  pivot_wider(names_from = "symbol", values_from = "raw_count_scaled_adjusted") %>%
  dplyr::select(sample_ct, COL6A3, CD226, CD244, CRTAM, KIR2DL4, KLRC1, 
                KLRC2, KLRC3, KLRC4, KLRK1, NCR1, NCR2, NCR3,
                CD80, CD83,ITGAX)

a <- reshape::melt(cor(x %>% dplyr::select(-sample_ct))[-c(1:1), c(1:1)]) %>%
  rownames_to_column("symbol")
p3 <- as_tibble(a )  %>%
  mutate(COL6A3 = "COL6A3") %>%
  mutate(symbol = fct_relevel(symbol, "CD80", "CD83","ITGAX",
                              "CD226", "CD244", "CRTAM", "KIR2DL4", "KLRC1", 
                              "KLRC2", "KLRC3", "KLRC4", "KLRK1", "NCR1", "NCR2", "NCR3")) %>%
  ggplot(aes(y = COL6A3, x = symbol, fill = value)) +
  geom_tile(stat = "identity") +
  theme_classic() +
  theme(text=element_text(size=16, family="sans"),
        panel.spacing.x = unit(0,"line"),
        axis.text.x = element_text(angle = -90, hjust = 0.5, vjust = 0.5),
        legend.position = "bottom") +
  labs(x = "", y = "", title = "") + 
  scale_fill_gradient2(low = "#006699", high = "#b30000", mid = "white",
                       midpoint = 0, 
                       name = "Correlation",
                       limits=range(-0.2, 0.2),
                       breaks = seq(-0.2, 0.1, 0.2))

ggsave(plot=p0, "C:/Rstudio/Chapter 3/Figures/Figure 4 COL6A3 BLCA.pdf", device = "pdf", height = 4, width = 4)
ggsave(plot=p1a, "C:/Rstudio/Chapter 3/Figures/Figure 4 COL6A3+cell BLCA.pdf", device = "pdf", height = 4, width = 8)
ggsave(plot=p2, "C:/Rstudio/Chapter 3/Figures/Figure 4 COL6A3+pth BLCA.pdf", device = "pdf", height = 4, width = 10)

x <- read_csv("C:/Rstudio/Chapter 3/output/WGCNA/BLCA/BLCA.module_eipgengene.xls.csv")%>%
  pivot_longer(-MEs, names_to = "sample_ct", values_to = "item") %>%
  split_median("MEs", "item") %>%
  filter(grepl("-Tumor", sample_ct)) %>%
  mutate(sample = gsub("-Tumor", "", sample_ct)) %>%
  clinical_combine("BLCA") %>%
  mutate(cat = MEs) %>%
  as.data.frame()
p4 <- ggsurvp_cat_item(x, T, 5, "jco", "","","", 10)

x <- cell %>% dplyr::select(sample_ct, mDC) %>%
  mutate(cat = "mDC", item = mDC) %>%
  split_median("cat", "item") %>%
  filter(grepl("-Tumor", sample_ct)) %>%
  mutate(sample = gsub("-Tumor", "", sample_ct)) %>%
  clinical_combine("BLCA")%>%
  mutate(item = fct_relevel(item,"L", "H")) %>%
  as.data.frame()
ggsurvp_cat_item(x, T, 1, "jco", "","","", 20)

x <- x %>% dplyr::select(sample_ct, pathologic_stage, total_living_days, vital_status, lymphovascular_invasion_present)
p <- ggsurvp_cat_item(x %>%
                   mutate(cat = "Stage", item = pathologic_stage) %>%
                   filter(item!=""), T, 1, "jco", "","","", 20)


