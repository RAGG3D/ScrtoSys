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
cancerlist = list("KICH", "READ", 
                  "LIHC", "STAD", "COAD", 
                  "PRAD", "THCA", "HNSC", "KIRC", "LUAD", 
                  "ESCA")

cancer0 = "BLCA"
WGCNA_everything <- function(cancer0){
  dir.create(paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0,"/"))
  type = "unsigned"
  exprMat <- paste0("C:/Rstudio/Chapter 3/output/WGCNA/",cancer0, "/", cancer0)
  # 
  # biweight mid-correlation & bicor
  # corType: pearson or bicor
  corType = "bicor"
  
  corFnc = ifelse(corType=="pearson", cor, bicor)

  maxPOutliers = ifelse(corType=="pearson",1,0.05)
  
  robustY = ifelse(corType=="pearson",T,F)
  WGCNA::allowWGCNAThreads()
  gene <- read_csv(paste0("C:/Rstudio/Chapter 3/data/Adjusted counts all gene/",  cancer0, ".csv")) %>%
    mutate(TMM_adjusted_log2 = log2(raw_count_scaled_adjusted+1))
  
  dataExpr <- gene %>%
    dplyr::select(sample_ct, symbol, TMM_adjusted_log2) %>%
    filter(grepl("Tumor", sample_ct)) %>%
    pivot_wider(names_from = "sample_ct", values_from = "TMM_adjusted_log2") %>%
    column_to_rownames("symbol") 
  
  
  m.mad <- apply(dataExpr,1,mad)
  dataExprVar <- dataExpr[which(m.mad > 0),]
  
  dataExpr <- as.data.frame(t(dataExprVar))
  

  gsg = goodSamplesGenes(dataExpr, verbose = 3)
  
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
      printFlush(paste("Removing genes:", 
                       paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", 
                       paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  sampleTree = hclust(dist(dataExpr), method = "average")
  pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0,"/",cancer0,  ".hcluster.tree.pdf"), onefile=F, paper="special", 
      width=10, height=7, bg="white", pointsize=6)
  plot(sampleTree, main = paste0(cancer0, " Sample clustering to detect outliers"), sub="", xlab="")
  dev.off()
  
  # sample network based on squared Euclidean distance note that we
  # transpose the data
  A = adjacency(t(dataExpr), type = "distance")
  # this calculates the whole network connectivity
  k = as.numeric(apply(A, 2, sum)) - 1
  # standardized connectivity
  Z.k = scale(k)
  # Designate samples as outlying if their Z.k value is below the threshold
  thresholdZ.k = -5  # often -2.5
  
  # the color vector indicates outlyingness (red)
  outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")
  
  # calculate the cluster tree using flahsClust or hclust
  sampleTree = hclust(as.dist(1 - A), method = "average")
  # Convert traits to a color representation: where red indicates high
  # values
  # traitColors = data.frame(numbers2colors(datTraits, signed = FALSE))
  # dimnames(traitColors)[[2]] = paste(names(datTraits), "C", sep = "")
  # datColors = data.frame(outlierC = outlierColor, traitColors)
  # Plot the sample dendrogram and the colors underneath.
  plotDendroAndColors(sampleTree, groupLabels = names(outlierColor), 
                      colors = outlierColor, 
                      main = "Sample dendrogram and trait heatmap")
  pdf(file = paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, " plotDendroAndColorsAndTree.pdf"))
  plotDendroAndColors(sampleTree, groupLabels = names(outlierColor), 
                      colors = outlierColor, 
                      main = "Sample dendrogram and trait heatmap")
  dev.off()
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  
  sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                          networkType=type, verbose=5)
  power = sft$powerEstimate
  power = 4
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # R-square=0.85
  abline(h=0.85,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
       cex=cex1, col="red")
  
  pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " .softpower.pdf"), onefile=F, paper="special", 
      bg="white", pointsize=6)
  par(mfrow = c(1,2))
  cex1 = 0.9
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # R-square=0.85
  abline(h=0.85,col="red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
       cex=cex1, col="red")
  dev.off()
  if (is.na(power)){
    print("Using experience power since no suitable power found.")
    power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                   ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                 ifelse(type == "unsigned", 6, 12))       
                   )
    )
  }
  
  print(paste("Finally chooosed power is :", power))
  net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                         TOMType = type, minModuleSize = 25,
                         networkType = type,
                         reassignThreshold = 0, mergeCutHeight = 0.2,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=F, corType = corType, 
                         maxPOutliers=maxPOutliers, loadTOM=F,
                         TOMDenom = "min",  deepSplit = 1,
                         stabilityCriterion = "Individual fraction", 
                         saveTOMFileBase = paste0(exprMat, ".tom"),
                         verbose = 3, randomSeed=1117)
  # Convert labels to colors for plotting
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  # Plot the dendrogram and the module colors underneath
  pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " plotDendroAndColors.pdf"), onefile=F, paper="special", 
      bg="white", pointsize=6)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors,
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.5,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  dynamicColors <- labels2colors(net$unmergedColors)
  plotDendroAndColors(net$dendrograms[[1]], cbind(dynamicColors,moduleColors),
                      c("Dynamic Tree Cut", "Module colors"),
                      dendroLabels = FALSE, hang = 0.5,
                      addGuide = TRUE, guideHang = 0.05)

  gene_module <- data.frame(ID=colnames(dataExpr), module=moduleColors)
  gene_module = gene_module[order(gene_module$module),]
  write.table(gene_module,file=paste0(exprMat,".gene_module.xls"),
              sep="\t",quote=F,row.names=F)
  
  MEs = net$MEs

  MEs_col = MEs
  colnames(MEs_col) = paste0("ME", labels2colors(
    as.numeric(str_replace_all(colnames(MEs),"ME",""))))
  MEs_col = orderMEs(MEs_col)
  
  MEs_colt = as.data.frame(t(MEs_col))
  colnames(MEs_colt) = rownames(dataExpr)
  write.table(MEs_colt,file=paste0(exprMat,".module_eipgengene.xls"),
              sep="\t",quote=F)

  pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " plotEigengeneNetworks.pdf"), onefile=F, paper="special", 
      bg="white", pointsize=6)
  plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                        marDendro = c(3,3,2,4),
                        marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                        xLabelsAngle = 90)
  dev.off()
  
  hubs = chooseTopHubInEachModule(dataExpr, colorh=moduleColors, power=power, type=type)
  
########################################################################################################################
  
  
  trait <- paste0("C:/Rstudio/Chapter 3/data/TCGA clinical/clinical_", tolower(cancer0)," WGCNA.csv")
  if(trait != "") {
    sampleName = as.data.frame(rownames(dataExpr))%>%
      tidybulk::rename("sample" = "rownames(dataExpr)")%>%
      left_join(read_csv(trait)) 
      df <- sampleName
      df[] <- lapply(df, function(x) {x[grepl("[", x, fixed = TRUE)] <- NA; x})
      df <- as.data.frame(df)
      df[, c("bcr_patient_uuid" )] <- sapply(df[, c("bcr_patient_uuid" )], unclass)
      colname = colnames(df)
      x <- map_dfc(as.list(1:ncol(df)), function(i){
        y <- (df %>% dplyr::select(contains(colname[i])) %>%
                group_by(!!as.name(colname[i])) %>%
                summarise() %>%
                mutate(l = c(1:nrow(.))) %>%
                left_join(df %>% dplyr::select(contains(colname[i]))))[,2]
        colnames(y) <- colname[i]
      })

    
    traitData = read_csv(paste0("C:/Rstudio/Chapter 3/data/TCGA clinical/clinical_", tolower(cancer0)," WGCNA.csv")) %>%
      column_to_rownames("sample")
    #sampleTree2 = hclust(dist(dataExpr), method = "average")
    # Convert traits to a color representation: white means low, red means high, grey means missing entry
    traitColors = numbers2colors(traitData, signed = FALSE);
    # Plot the sample dendrogram and the colors underneath.
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " Sample dendrogram and trait heatmap.pdf"), onefile=F, paper="special", 
        bg="white", pointsize=6)
    plotDendroAndColors(sampleTree, traitColors,
                        groupLabels = names(traitData), 
                        main = "Sample dendrogram and trait heatmap")
    dev.off()

    MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " Eigengene adjacency heatmap.pdf"), onefile=F, paper="special", 
        bg="white", pointsize=6)
    plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                          marDendro = c(3,3,2,4), marHeatmap = c(3,4,2,2),
                          plotDendrograms = T, xLabelsAngle = 90)
    
    
    dev.off()

    if (corType=="pearsoon") {
      modTraitCor = cor(MEs_col, traitData, use = "p")
      modTraitP = corPvalueStudent(modTraitCor, nSamples)
    } else {
      modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
      modTraitCor = modTraitCorP$bicor
      modTraitP   = modTraitCorP$p
    }

    textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
    dim(textMatrix) = dim(modTraitCor)
    
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " labeledHeatmap.pdf"), onefile=F, paper="special", 
        pointsize=6)
    labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
                   yLabels = colnames(MEs_col), 
                   cex.lab = 0.5, 
                   ySymbols = colnames(MEs_col), colorLabels = FALSE, 
                   colors = blueWhiteRed(50), 
                   textMatrix = textMatrix, setStdMargins = FALSE, 
                   cex.text = 0.5, zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    dev.off()
    
    mod_trait <- as.data.frame(melt(modTraitCor)) %>%
      tidybulk::rename("Cor" = "value") %>%
      inner_join(as.data.frame(melt(modTraitP))%>%
                   tidybulk::rename("P" = "value")) %>% 
      filter(P <= 0.05) 
    
    mod_trait %>%
      write.csv(paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " Module_traits relationships.csv"), row.names = F)
    
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
    
    
    pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " Module-trait relationships.pdf"), onefile=F, paper="special", 
        bg="white", pointsize=6)
    labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
                   yLabels = colnames(MEs_col), 
                   cex.lab = 0.5, 
                   ySymbols = colnames(MEs_col), colorLabels = FALSE, 
                   colors = blueWhiteRed(50), 
                   textMatrix = textMatrix, setStdMargins = FALSE, 
                   cex.text = 0.5, zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    
    dev.off()
    modTraitCorMelt = as.data.frame(modTraitCor)
    write.table(modTraitCorMelt,file=paste0(exprMat,".module_trait_correlation.xls"),
                sep="\t",quote=F)
    modTraitCorMelt$ID = rownames(modTraitCor)
    modTraitCorMelt = melt(modTraitCorMelt)
    colnames(modTraitCorMelt) <- c("Module","Trait","PersonCorrelationValue")
    modTraitPMelt = as.data.frame(modTraitP)
    write.table(modTraitPMelt,file=paste0(exprMat,".module_trait_correlationPvalue.xls"),
                sep="\t",quote=F)
    modTraitPMelt$ID = rownames(modTraitP)
    modTraitPMelt = melt(modTraitPMelt)
    colnames(modTraitPMelt) <- c("Module","Trait","Pvalue")
    #modTraitCorP = cbind(modTraitCorMelt, Pvalue=modTraitPMelt$Pvalue)
    modTraitCorP = merge(modTraitCorMelt, modTraitPMelt, by=c("Module","Trait"))
    write.table(modTraitCorP,file=paste0(exprMat,".module_trait_correlationPvalueMelt.xls"),
                sep="\t",quote=F,row.names=F)
  }
  
  


if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


geneTraitCorMelt = as.data.frame(geneTraitCor)
write.table(geneTraitCorMelt,file=paste0(exprMat,".gene_trait_correlation.xls"),
            sep="\t",quote=F)
geneTraitCorMelt$ID = rownames(geneTraitCor)
geneTraitCorMelt = melt(geneTraitCorMelt)
colnames(geneTraitCorMelt) <- c("Gene","Trait","PersonCorrelationValue")
geneTraitPMelt = as.data.frame(geneTraitP)
write.table(geneTraitPMelt,file=paste0(exprMat,".gene_trait_correlationPvalue.xls"),
            sep="\t",quote=F)
geneTraitPMelt$ID = rownames(geneTraitP)
geneTraitPMelt = melt(geneTraitPMelt)
colnames(geneTraitPMelt) <- c("Gene","Trait","Pvalue")
#geneTraitCorP = cbind(geneTraitCorMelt, Pvalue=geneTraitPMelt$Pvalue)
geneTraitCorP = merge(geneTraitCorMelt, geneTraitPMelt, by=c("Gene","Trait"))
write.table(geneTraitCorP,
            file=paste0(exprMat,".gene_trait_correlationPvalueMelt.xls"),
            sep="\t",quote=F,row.names=F)

#plot_me_trat <- cbind(dynamicColors,moduleColors,geneTraitCor)
geneTraitCorColor <- numbers2colors(geneTraitCor)

plotDendroAndColors(net$dendrograms[[1]],
                    cbind(dynamicColors,moduleColors,geneTraitCorColor),
                    c("Dynamic Tree Cut", "Module colors", colnames(geneTraitCor)),
                    dendroLabels = FALSE, hang = 0.5,
                    addGuide = TRUE, guideHang = 0.05)

pdf(file=paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " Dynamic Tree Cut.pdf"), onefile=F, paper="special", 
    bg="white", pointsize=6)
plotDendroAndColors(net$dendrograms[[1]],
                    cbind(dynamicColors,moduleColors,geneTraitCorColor),
                    c("Dynamic Tree Cut", "Module colors", colnames(geneTraitCor)),
                    dendroLabels = FALSE, hang = 0.5,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

module = "green"
pheno = "pathologic_stage"
modNames = substring(colnames(MEs_col), 3)

module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))

moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
pdf(file = "verboseScatterplot.pdf")
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   abline = T,
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

module = "green"
if(T){
  gene <- colnames(dataExpr) 
  inModule <- moduleColors==module
  modgene <- gene[inModule]
  
  TOM <- TOMsimilarityFromExpr(dataExpr,power=power)
  modTOM <- TOM[inModule,inModule]
  dimnames(modTOM) <- list(modgene,modgene)
  
  nTop = 100
  IMConn = softConnectivity(dataExpr[, modgene]) 
  top = (rank(-IMConn) <= nTop) 
  filter_modTOM <- modTOM[top, top]
  
  # for visANT
  vis <- exportNetworkToVisANT(filter_modTOM,
                               file = paste("step8_visANTinput-",module,".txt",sep = ""),
                               weighted = T,threshold = 0)
  # for cytoscape
  cyt <- exportNetworkToCytoscape(filter_modTOM,
                                  edgeFile = paste(paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " step8_CytoscapeInput-edges-"), paste(module, collapse="-"), ".txt", sep=""),
                                  nodeFile = paste(paste0("C:/Rstudio/Chapter 3/output/WGCNA/", cancer0, "/",cancer0, " step8_CytoscapeInput-nodes-"), paste(module, collapse="-"), ".txt", sep=""),
                                  weighted = TRUE,
                                  threshold = 0.15,
                                  nodeNames = modgene[top], 
                                  nodeAttr = moduleColors[inModule][top])
}

}

cancer0 = "BLCA"
map(list("BLCA"), function(cancer0){
  WGCNA_everything(cancer0)
})




