library(tidyverse)
library(tidybulk)
library(org.Hs.eg.db)
library(Hmisc)
library(msigdbr)
ref = readRDS("C:/Rstudio/Chapter 3/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("C:/Rstudio/Chapter 3/data/cellname.csv")))$`Cell Type`
scrto <- read_csv("C:/Rstudio/Chapter 3/data/comprehensive secreted protein gene list.csv")



h_df = msigdbr(species = "Homo sapiens", category = "H")
####*Functions*####
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

H_scrto <- h_df %>% filter(gene_symbol %in% scrto$symbol)
#####*H whole genes*#####"KICH", "READ", "LIHC", "BLCA", "STAD", 
map(list("KICH", "BLCA", "READ", "LIHC", "STAD"),
    function(cancer0){
      a <- read_csv(
        paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",cancer0,".csv")
      )
      sample0 = unique(a$sample_ct)
      
      #samplelist = sample0[(which(sample0 == "TCGA-67-6217-Tumor", arr.ind=TRUE)+1):length(sample0)]
      samplelist = sample0
      map(as.list(unique(scrto$symbol)), function(scrto1){
        dir.create(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Single scrto removal/",cancer0,"/"))
        dir.create(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Single scrto removal/",cancer0,"/", scrto1, "/"))
        map(as.list(samplelist), function(sample1){
          if(scrto1 %in% H_scrto$gene_symbol) {
            x <- (a) %>%
          filter(sample_ct == sample1) %>%
          filter(symbol != scrto1) %>%
          inner_join(toTable(org.Hs.egSYMBOL)) 
          } else {
            x <- (a) %>%
              filter(sample_ct == sample1) %>%
              inner_join(toTable(org.Hs.egSYMBOL)) 
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
          write.csv(paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/Hallmarks/Single scrto removal/",cancer0,"/", scrto1, "/", sample1, ".csv"), row.names = F)
      })
      })
       
    })

