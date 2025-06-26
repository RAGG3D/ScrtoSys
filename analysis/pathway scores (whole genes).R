if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(tidyverse)
library(tidybulk)
library(org.Hs.eg.db)
library(Hmisc)
ref = readRDS("C:/Rstudio/Chapter 3/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("C:/Rstudio/Chapter 3/data/cellname.csv")))$`Cell Type`
scrto <- read_csv("C:/Rstudio/Chapter 3/data/comprehensive secreted protein gene list.csv")

install.packages("msigdbr")
library(msigdbr)

h_df = msigdbr(species = "Homo sapiens", category = "H")
###**###
#####*C2*#####
map(list("KICH", "LIHC"),
    function(cancer0){
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      sample0 = unique(a$sample_ct)
      
      #samplelist = sample0[(which(sample0 == "TCGA-67-6217-Tumor", arr.ind=TRUE)+1):length(sample0)]
      samplelist = sample0
      map(as.list(samplelist), function(sample1){
        x <- (a) %>%
          filter(sample_ct == sample1) %>%
          filter(!symbol %in% scrto$symbol) %>%
          inner_join(toTable(org.Hs.egSYMBOL)) 
        
        
        entrez_rank_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){
          my_gene_collection = msigdbr::msigdbr(species = species) %>%  filter( tolower(gs_cat) %in% tolower(gene_collections) )
          my_gene_collection |>
            nest(data = -gs_cat) |>
            mutate(fit =
                     map(
                       data,
                       ~ 	clusterProfiler::GSEA(
                         my_entrez_rank,
                         TERM2GENE=.x %>% dplyr::select(gs_name, entrez_gene),
                         pvalueCutoff = 1
                       )
                       
                     )) |>
            mutate(test =
                     map(
                       fit,
                       ~ .x |>
                         # ggplot2::fortify(showCategory=Inf) %>%
                         as_tibble() |>
                         rowid_to_column(var = "idx_for_plotting")
                       #%>%
                       #	mutate(plot = future_imap(ID, ~ enrichplot::gseaplot2(fit, geneSetID = .y, title = .x)))
                       
                     )) |>
            dplyr::select(-data)
          
        }
        
        
        .data = x
        .entrez = x$gene_id
        .arrange_desc = x$raw_count_scaled_adjusted
        species = "Homo sapiens"
        gene_sets = c("c2")
        .sample = x$sample_ct
        
        a <- x %>%
          pivot_transcript(.transcript = gene_id) %>%
          arrange(desc(raw_count_scaled_adjusted)) %>%
          dplyr::select(gene_id, raw_count_scaled_adjusted) %>%
          deframe() %>%
          entrez_rank_to_gsea(species, gene_collections = gene_sets) 
        
        as.data.frame(a$test)%>%
          write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/",cancer0,"/", sample1, ".csv"), row.names = F)
      }) 
    })

cancer0 = "ESCA"
map(list("KICH", "LIHC"),
    function(cancer0){
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      sample0 = unique(a$sample_ct)
      
      samplelist = sample0
      #samplelist = samplelist[1:10]
      b <- map_dfr(as.list(samplelist), function(sample1){
        read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/",cancer0,"/", sample1, ".csv")) %>% 
          mutate(sample_ct = sample1) %>%
          dplyr::select(sample_ct, ID, enrichmentScore) %>%
          pivot_wider(names_from = ID, values_from = enrichmentScore)
      })
      x <- a %>% filter(symbol%in% unique(scrto$symbol)) %>%
        dplyr::select(sample_ct, symbol, raw_count_scaled_adjusted) %>%
        pivot_wider(names_from = symbol, values_from = raw_count_scaled_adjusted) %>%
        inner_join(b)
      
      
      a <- reshape::melt(cor(x %>% dplyr::select(-sample_ct))[-c(1:4919), c(1:4919)]) 
      a.p <- rcorr(as.matrix((x %>% dplyr::select(-sample_ct)))) 
      a %>% inner_join(as.data.frame(a.p$P[-c(1:4919), c(1:4919)]) %>% 
                         mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:4919), c(1:4919)]))) %>% 
                         gather(X2, p.value, -X1)) %>%
        write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/", cancer0, " pathway+scrto corr.csv"), row.names = F)
      
    })




#####*H whole genes*#####"KICH", "READ", "LIHC", "BLCA", "STAD", 
map(list("PRAD", "THCA", "HNSC", "KIRC", "LUAD", 
         "ESCA", "BRCA"),
    function(cancer0){
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      sample0 = unique(a$sample_ct)
      
      #samplelist = sample0[(which(sample0 == "TCGA-67-6217-Tumor", arr.ind=TRUE)+1):length(sample0)]
      samplelist = sample0
      map(as.list(samplelist), function(sample1){
        x <- (a) %>%
          filter(sample_ct == sample1) %>%
          inner_join(toTable(org.Hs.egSYMBOL)) 
        
        
        entrez_rank_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){
          my_gene_collection = msigdbr::msigdbr(species = species) %>%  filter( tolower(gs_cat) %in% tolower(gene_collections) )
          my_gene_collection |>
            nest(data = -gs_cat) |>
            mutate(fit =
                     map(
                       data,
                       ~ 	clusterProfiler::GSEA(
                         my_entrez_rank,
                         TERM2GENE=.x %>% dplyr::select(gs_name, entrez_gene),
                         pvalueCutoff = 1
                       )
                       
                     )) |>
            mutate(test =
                     map(
                       fit,
                       ~ .x |>
                         # ggplot2::fortify(showCategory=Inf) %>%
                         as_tibble() |>
                         rowid_to_column(var = "idx_for_plotting")
                       #%>%
                       #	mutate(plot = future_imap(ID, ~ enrichplot::gseaplot2(fit, geneSetID = .y, title = .x)))
                       
                     )) |>
            dplyr::select(-data)
          
        }
        
        
        .data = x
        .entrez = x$gene_id
        .arrange_desc = x$raw_count_scaled_adjusted
        species = "Homo sapiens"
        gene_sets = c("h")
        .sample = x$sample_ct
        
        a <- x %>%
          pivot_transcript(.transcript = gene_id) %>%
          arrange(desc(raw_count_scaled_adjusted)) %>%
          dplyr::select(gene_id, raw_count_scaled_adjusted) %>%
          deframe() %>%
          entrez_rank_to_gsea(species, gene_collections = gene_sets) 
        
        as.data.frame(a$test)%>%
          write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/",cancer0,"/", sample1, ".csv"), row.names = F)
      }) 
    })


##### Secretome related pathway in Tumor
map(list("KICH", "READ", 
         "LIHC", "BLCA", "STAD", "COAD", 
         "PRAD", "THCA", "HNSC", "KIRC", "LUAD", 
         "ESCA"),
    function(cancer0){
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      
      sample0 = unique((a %>% filter(grepl("Tumor", sample_ct)))$sample_ct)
      samplelist = sample0
      
      b <- map_dfr(as.list(samplelist), function(sample1){
        read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/",cancer0,"/", sample1, ".csv")) %>% 
          mutate(sample_ct = sample1) %>%
          dplyr::select(sample_ct, ID, enrichmentScore) %>%
          pivot_wider(names_from = ID, values_from = enrichmentScore)
      })
      x <- a %>% filter(symbol%in% unique(scrto$symbol)) %>%
        dplyr::select(sample_ct, symbol, raw_count_scaled_adjusted) %>%
        pivot_wider(names_from = symbol, values_from = raw_count_scaled_adjusted) %>%
        inner_join(b)
      
      ###*Models package in tidyR*###
      a <- reshape::melt(cor(x %>% dplyr::select(-sample_ct))[-c(1:2329), c(1:2329)]) 
      a.p <- rcorr(as.matrix((x %>% dplyr::select(-sample_ct)))) 
      a %>% inner_join(as.data.frame(a.p$P[-c(1:2329), c(1:2329)]) %>% 
                         mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:2329), c(1:2329)]))) %>% 
                         gather(X2, p.value, -X1)) %>%
        write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene pathway+scrto corr.csv"), row.names = F)
      
    })

##### Secretome related pathway in Normal
map(list("KICH", "READ", 
         "LIHC", "BLCA", "STAD", "COAD", 
         "PRAD", "THCA", "HNSC", "KIRC", "LUAD", 
         "ESCA"),
    function(cancer0){
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      
      sample0 = unique((a %>% filter(grepl("Normal", sample_ct)))$sample_ct)
      samplelist = sample0
      
      b <- map_dfr(as.list(samplelist), function(sample1){
        read_csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/",cancer0,"/", sample1, ".csv")) %>% 
          mutate(sample_ct = sample1) %>%
          dplyr::select(sample_ct, ID, enrichmentScore) %>%
          pivot_wider(names_from = ID, values_from = enrichmentScore)
      })
      x <- a %>% filter(symbol%in% unique(scrto$symbol)) %>%
        dplyr::select(sample_ct, symbol, raw_count_scaled_adjusted) %>%
        pivot_wider(names_from = symbol, values_from = raw_count_scaled_adjusted) %>%
        inner_join(b)
      
      
      a <- reshape::melt(cor(x %>% dplyr::select(-sample_ct))[-c(1:2329), c(1:2329)]) 
      a.p <- rcorr(as.matrix((x %>% dplyr::select(-sample_ct)))) 
      a %>% inner_join(as.data.frame(a.p$P[-c(1:2329), c(1:2329)]) %>% 
                         mutate(X1 = rownames(as.data.frame(a.p$P[-c(1:2329), c(1:2329)]))) %>% 
                         gather(X2, p.value, -X1)) %>%
        write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/", cancer0, " whole gene normal pathway+scrto corr.csv"), row.names = F)
      
    })

#####*H whole gene set*#####
map(list("KICH", "READ", 
         "LIHC", "BLCA", "STAD", "COAD", 
         "PRAD", "THCA", "HNSC", "KIRC", "LUAD", 
         "ESCA", "BRCA"),
    function(cancer0){
      dir.create(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/",cancer0,"/"))
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      sample0 = unique(a$sample_ct)
      
      #samplelist = sample0[(which(sample0 == "TCGA-67-6217-Tumor", arr.ind=TRUE)+1):length(sample0)]
      samplelist = sample0
      map(as.list(samplelist), function(sample1){
        x <- (a) %>%
          filter(sample_ct == sample1) %>%
          inner_join(toTable(org.Hs.egSYMBOL)) 
        
        
        entrez_rank_to_gsea = function(my_entrez_rank, species, gene_collections  = NULL){
          my_gene_collection = msigdbr::msigdbr(species = species) %>%  filter( tolower(gs_cat) %in% tolower(gene_collections) )
          my_gene_collection |>
            nest(data = -gs_cat) |>
            mutate(fit =
                     map(
                       data,
                       ~ 	clusterProfiler::GSEA(
                         my_entrez_rank,
                         TERM2GENE=.x %>% dplyr::select(gs_name, entrez_gene),
                         pvalueCutoff = 1
                       )
                       
                     )) |>
            mutate(test =
                     map(
                       fit,
                       ~ .x |>
                         # ggplot2::fortify(showCategory=Inf) %>%
                         as_tibble() |>
                         rowid_to_column(var = "idx_for_plotting")
                       #%>%
                       #	mutate(plot = future_imap(ID, ~ enrichplot::gseaplot2(fit, geneSetID = .y, title = .x)))
                       
                     )) |>
            dplyr::select(-data)
          
        }
        
        
        .data = x
        .entrez = x$gene_id
        .arrange_desc = x$raw_count_scaled_adjusted
        species = "Homo sapiens"
        gene_sets = c("h")
        .sample = x$sample_ct
        
        a <- x %>%
          pivot_transcript(.transcript = gene_id) %>%
          arrange(desc(raw_count_scaled_adjusted)) %>%
          dplyr::select(gene_id, raw_count_scaled_adjusted) %>%
          deframe() %>%
          entrez_rank_to_gsea(species, gene_collections = gene_sets) 
        
        as.data.frame(a$test)%>%
          write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Whole gene sets/",cancer0,"/", sample1, ".csv"), row.names = F)
      }) 
    })

