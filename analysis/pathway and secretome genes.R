source("E:/PhD Projects/Rprojects/SProtosurv/src/functions.R")
ref = readRDS("C:/Rstudio/Chapter 3/data/my_ref.rds")[,-c(1,4,11,17,21,24,28,31,32,33)]
colnames(ref) <- (data.frame(Var2 = colnames(ref)) %>%
                    inner_join(read_csv("C:/Rstudio/Chapter 3/data/cellname.csv")))$`Cell Type`
scrto <- read_csv("C:/Rstudio/Chapter 3/data/comprehensive secreted protein gene list.csv")



cancer0 = "KICH"

map(list( "PRAD",
          "ESCA"),
    function(cancer0){
      paste0("C:/Rstudio/Chapter 3/output/TCGA pathway scores/",cancer0,"/", sample1, ".csv")
      
    })
