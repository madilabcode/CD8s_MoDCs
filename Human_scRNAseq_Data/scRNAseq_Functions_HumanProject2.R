setwd("C:\\Kfir_Thesis_Asaf_Madi\\Spatial_Transcriptomics\\Cell_Cell_Communication")
library(findPC)
library(readr)
library(dplyr)
library(dbplyr)
library(tibble)
library(anndata)
library(Seurat)
#library(SeuratData)
#library(SeuratDisk)
library(rhdf5)
library(RColorBrewer)
library(future)
library(ggplot2)
library(reticulate)
library(tidyverse)
library(plyr)
library(FNN)
library(stringr)
library(RANN)
library(ggsignif)
library(patchwork)
library(sleepwalk)
library(gridExtra)
library(janitor)
library(purrr)
library(ggpubr)
library(grid)



#GENERIC FUNCTION FOR GENERATING BOX PLOT OF SIGNATURE SCORES

Box_Plot_Sig_Scores <- function(SlideRDS, up_list, down_list, Slide_Name, Cytokine_Name, saveRDS_flag = F) {
  Slide_Seurat <- readRDS(SlideRDS)
  
  Slide_Seurat@assays[["RNA"]] <- Slide_Seurat@assays[["SCT"]]
  
  Signature_Scores <- make_signature_old(Slide_Seurat, up = up_list, idents = c(0,1),down = down_list)
  
  if (saveRDS_flag == T) {
    saveRDS(Signature_Scores,file=paste0("C:\\Kfir_Thesis_Asaf_Madi\\Spatial_Transcriptomics\\Cell_Cell_Communication\\Signature_Scores\\Merging_All_Plots\\",Cytokine_Name,"\\",Slide_Name,"_Signature_Scores_List_",Cytokine_Name,".rds"))
  } 
  
  #Signature_Scores <- readRDS("C:\\Kfir_Thesis_Asaf_Madi\\Spatial_Transcriptomics\\Cell_Cell_Communication\\Signature_Scores\\IL10\\MPM01_Signature_Scores_List_IL10.rds")
  
  Slide_meta_data_sig_scores <- Signature_Scores[[2]]$data
  
  Slide_Intra_meta_data_sig_scores <- subset(Slide_meta_data_sig_scores, seurat_clusters == 1,select=c('seurat_clusters','SigUint'))
  
  Slide_Extra_meta_data_sig_scores <- subset(Slide_meta_data_sig_scores, seurat_clusters == 0,select=c('seurat_clusters','SigUint'))
  
  Slide_Intra_Extra_meta_data_sig_scores <- rbind(Slide_Intra_meta_data_sig_scores, Slide_Extra_meta_data_sig_scores)
  Slide_Intra_Extra_meta_data_sig_scores$seurat_clusters <- as.character(Slide_Intra_Extra_meta_data_sig_scores$seurat_clusters)
  
  
  Signature_Scores_t_test <- t.test(subset(Slide_Intra_Extra_meta_data_sig_scores, seurat_clusters==1,select='SigUint'), subset(Slide_Intra_Extra_meta_data_sig_scores, seurat_clusters==0, select='SigUint'), alternative = "two.sided", var.equal = FALSE)
  
  p_value <- Signature_Scores_t_test[3][[1]]
  
  if(p_value < 0.05){
    #BOXPLOT
    p <- ggplot(Slide_Intra_Extra_meta_data_sig_scores, aes(x=seurat_clusters, y=SigUint, fill = seurat_clusters)) + ggtitle(paste0("Box Plot of ", Cytokine_Name," Gene Set Signature Scores in Extra-tumorality vs Intra-tumorality CD8 Cells")) + theme(plot.title = element_text(size=10)) + scale_x_discrete(labels = c("Extra","Intra"))+
      geom_boxplot(width=0.3)+ stat_summary(fun.y=mean,  geom="point", color="red")+ xlab(paste0(Slide_Name, " Slide")) +ylab("SigUint Expression") + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())  +geom_signif(comparisons = list(c("0","1")),map_signif_level = TRUE, test = "t.test", tip_length = 0.02, vjust = 0) + scale_fill_manual(values=c("#9b62cb", "#dddae0"),name = "Tumorality Position", labels = c("Extra Mean","Intra Mean"))
  }
  else{
    #BOXPLOT
    p <- ggplot(Slide_Intra_Extra_meta_data_sig_scores, aes(x=seurat_clusters, y=SigUint, fill = seurat_clusters,alpha=1)) + ggtitle(paste0("Box Plot of ", Cytokine_Name," Gene Set Signature Scores in Extra-tumorality vs Intra-tumorality CD8 Cells")) + theme(plot.title = element_text(size=10)) + scale_x_discrete(labels = c("Extra","Intra"))+
      geom_boxplot(width=0.3,aes(alpha=1))+ stat_summary(fun.y=mean,  geom="point", color="red")+ xlab(paste0(Slide_Name, " Slide")) +ylab("SigUint Expression") + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())  +geom_signif(comparisons = list(c("0","1")),map_signif_level = TRUE, test = "t.test", tip_length = 0.02, vjust = 0) + scale_fill_manual(values=c("#9b62cb", "#dddae0"),name = "Tumorality Position", labels = c("Extra Mean","Intra Mean"), aes(alpha = 1))
  }
  
  #BOXPLOT
  #p <- ggplot(Slide_Intra_Extra_meta_data_sig_scores, aes(x=seurat_clusters, y=SigUint, fill = seurat_clusters)) + ggtitle(paste0("Box Plot of ", Cytokine_Name," Gene Set Signature Scores in Extra-tumorality vs Intra-tumorality CD8 Cells")) + theme(plot.title = element_text(size=13)) + scale_x_discrete(labels = c("Extra","Intra"))+
  #geom_boxplot(width=0.3)+ stat_summary(fun.y=mean,  geom="point", color="red")+ xlab(paste0(Slide_Name, " Slide")) +ylab("SigUint Expression") + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())  +geom_signif(comparisons = list(c("0","1")),map_signif_level = TRUE, test = "t.test", tip_length = 0.02, vjust = 0) + scale_fill_manual(values=c("#9b62cb", "#dddae0"),name = "Tumorality Position", labels = c("Extra Mean","Intra Mean"))
  
  #ggsave(p ,file = paste0(".\\Signature_Scores\\",Cytokine_Name,"\\",Slide_Name,"_Intra_vs_Extra_Signature_BoxPlot_",Cytokine_Name,".png"), dpi=300, width=10, height=7)
  
  return(p)
}


#Merging all plots

returnBoxPlot <- function(signatureScoresRDS, Cytokine_Name, basedir) {
  Slide_Name <- sub("_Signature_Score_*.*", "", signatureScoresRDS) 
  Signature_Scores <- readRDS(paste0(basedir,"/",signatureScoresRDS))
  
  Slide_meta_data_sig_scores <- Signature_Scores[[2]]$data
  
  Slide_Intra_meta_data_sig_scores <- subset(Slide_meta_data_sig_scores, seurat_clusters == 1,select=c('seurat_clusters','SigUint'))
  
  Slide_Extra_meta_data_sig_scores <- subset(Slide_meta_data_sig_scores, seurat_clusters == 0,select=c('seurat_clusters','SigUint'))
  
  Slide_Intra_Extra_meta_data_sig_scores <- rbind(Slide_Intra_meta_data_sig_scores, Slide_Extra_meta_data_sig_scores)
  Slide_Intra_Extra_meta_data_sig_scores$seurat_clusters <- as.character(Slide_Intra_Extra_meta_data_sig_scores$seurat_clusters)
  
  
  Signature_Scores_t_test <- t.test(subset(Slide_Intra_Extra_meta_data_sig_scores, seurat_clusters==1,select='SigUint'), subset(Slide_Intra_Extra_meta_data_sig_scores, seurat_clusters==0, select='SigUint'), alternative = "two.sided", var.equal = FALSE)
  
  p_value <- Signature_Scores_t_test[3][[1]]
  
  if(p_value < 0.05){
    #BOXPLOT
    p <- ggplot(Slide_Intra_Extra_meta_data_sig_scores, aes(x=seurat_clusters, y=SigUint, fill = seurat_clusters)) + ggtitle(paste0("Box Plot of ", Cytokine_Name," Gene Set Signature Scores in Extra-tumorality vs Intra-tumorality CD8 Cells")) + theme(plot.title = element_text(size=10)) + scale_x_discrete(labels = c("Extra","Intra"))+
      geom_boxplot(width=0.3)+ stat_summary(fun.y=mean,  geom="point", color="red")+ xlab(paste0(Slide_Name, " Slide")) +ylab("SigUint Expression") + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())  +geom_signif(comparisons = list(c("0","1")),map_signif_level = TRUE, test = "t.test", tip_length = 0.02, vjust = 0) + scale_fill_manual(values=c("#9b62cb", "#dddae0"),name = "Tumorality Position", labels = c("Extra Mean","Intra Mean"))
  }
  else{
    #BOXPLOT
    p <- ggplot(Slide_Intra_Extra_meta_data_sig_scores, aes(x=seurat_clusters, y=SigUint, fill = seurat_clusters,alpha=1)) + ggtitle(paste0("Box Plot of ", Cytokine_Name," Gene Set Signature Scores in Extra-tumorality vs Intra-tumorality CD8 Cells")) + theme(plot.title = element_text(size=10)) + scale_x_discrete(labels = c("Extra","Intra"))+
      geom_boxplot(width=0.3,aes(alpha=1))+ stat_summary(fun.y=mean,  geom="point", color="red")+ xlab(paste0(Slide_Name, " Slide")) +ylab("SigUint Expression") + theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())  +geom_signif(comparisons = list(c("0","1")),map_signif_level = TRUE, test = "t.test", tip_length = 0.02, vjust = 0) + scale_fill_manual(values=c("#9b62cb", "#dddae0"),name = "Tumorality Position", labels = c("Extra Mean","Intra Mean"), aes(alpha = 1))
  }
  return(p)
}



#GENERIC FUNCTION FOR GENERATING META DATA and PVAL OF SIGNATURE SCORES

metadata_pval_make_signature <- function(make_sig_obj, idents, ttest_flag = FALSE) {
  
  metadata <- make_sig_obj[[2]]$data
  
  sig_scores_binded <- data.frame()
  
  for (i in idents) {
    
    sig_scores_ <- subset(metadata, seurat_clusters == i,select=c('seurat_clusters','SigUint'))
    #assign(paste0("sig_scores_", i), sig_scores_)
    sig_scores_binded <- rbind(sig_scores_binded, sig_scores_)
  }
  
  sig_scores_binded$seurat_clusters <- as.character(sig_scores_binded$seurat_clusters)
  
  p_value <- "Undefined"
  
  if ( ttest_flag == TRUE) { #If TRUE then 'idents' must be of length 2
    
    Signature_Scores_t_test <- t.test(subset(sig_scores_binded, seurat_clusters==idents[1],select='SigUint'), subset(sig_scores_binded, seurat_clusters==idents[2], select='SigUint'), alternative = "two.sided", var.equal = FALSE)
    
    p_value <- Signature_Scores_t_test[3][[1]]
    
  }
  #Else ANOVA?
  
  return( list(sig_scores_binded, p_value) )
}



#GENERIC FUNCTION FOR GENERATING BOX PLOT OF SIGNATURE SCORES

Box_Plot_By_Make_Sig <- function(Seurat_Obj ,metadata_pval_obj, title = "", xlabs = NULL, xlabtitle = "", ttest_flag = FALSE, bicolors = c("#9b62cb", "#dddae0"), legend_name = "") {
  idents = c(unique(metadata_pval_obj[[1]]$seurat_clusters))
  if (is.null(xlabs)) { #Obtaining the assigned names of the numerical idents
    #xlabs <-  unique(subset(Seurat_Obj@meta.data,seurat_clusters %in% idents)$orig.ident)
    xlabs <- sort(as.vector(unique(as.data.frame(subset(Seurat_Obj@meta.data,seurat_clusters %in% idents)$orig.ident))[,1]))
  } else {
    xlabs <- xlabs[as.character(idents)]
  }
  if( metadata_pval_obj[[2]] < 0.05 & ttest_flag == TRUE ) {
    
    boxPlot <- 
      ggplot(metadata_pval_obj[[1]], aes(x=seurat_clusters, y=SigUint, fill = seurat_clusters)) +
      ggtitle(title) + 
      theme(plot.title = element_text(size=10)) + 
      scale_x_discrete(labels = xlabs) +
      geom_boxplot(width=0.3) + 
      stat_summary(fun.y=mean,  geom="point", color="red") + 
      xlab(xlabtitle) + 
      ylab("SigUint Expression") + 
      theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) + 
      geom_signif(comparisons = list(idents),map_signif_level = TRUE, test = "t.test", tip_length = 0.02, vjust = 0) + 
      scale_fill_manual(values=bicolors,name = legend_name, labels = c(paste0(xlabs[1]," Mean"), paste0(xlabs[2]," Mean")))
    
  } else if ( metadata_pval_obj[[2]] >= 0.05 & ttest_flag == TRUE ) {
    
    boxPlot <- 
      ggplot(metadata_pval_obj[[1]], aes(x=seurat_clusters, y=SigUint, fill = seurat_clusters, alpha = 1)) +
      ggtitle(title) +
      theme(plot.title = element_text(size=10)) +
      scale_x_discrete(labels = xlabs) +
      geom_boxplot(width=0.3, aes(alpha = 1)) +
      stat_summary(fun.y=mean,  geom="point", color="red") +
      xlab(xlabtitle) + ylab("SigUint Expression") +
      theme(axis.line = element_line(colour = "black"), panel.grid.major =   element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank()) + 
      geom_signif(comparisons = list(idents),map_signif_level = TRUE, test = "t.test", tip_length = 0.02, vjust = 0) + 
      scale_fill_manual(values=bicolors,name = legend_name, labels = c(paste0(xlabs[1]," Mean"), paste0(xlabs[2]," Mean")), aes(alpha = 1))
    
  }# else if ( metadata_pval_obj[[2]] >= 0.05 & ttest_flag == TRUE ) {
  #geom_signif with multiple comparisons??
  #}
  return(boxPlot)
}






Seurat_Single_Cell_Analysis <- function(Seurat_Obj, cluster_resolution = 0.8, prplxty = 2,npcs_value = 50) {
  
  #Applying normalization on the data - stored in @data
  Seurat_Obj = NormalizeData(Seurat_Obj)
  #Find differentially expressed genes between cells
  Seurat_Obj = FindVariableFeatures(Seurat_Obj)
  #Storing gene names
  all.genes = rownames(Seurat_Obj)
  #Scaling the raw data - stored in @scale
  Seurat_Obj = ScaleData(Seurat_Obj, vars.to.regress = "percent.mt")
  #Perform linear dimensional reduction - PCA - stored under @reductions$pca
  Seurat_Obj = RunPCA(Seurat_Obj, features = VariableFeatures(object = Seurat_Obj), npcs = npcs_value)
  #Determine the ‘dimensionality’ of the dataset - which PCs to take
  #ElbowPlot(Seurat_CD8A, ndims = 40) + theme_classic()
  
  #TAKING THE MAX number of PCs from the list of all methods to find the elbow point
  n_dims_PCA <- max(findPC(sort(Seurat_Obj@reductions$pca@stdev,decreasing = TRUE), method='all', figure = T, number=40))
  #Finding neighbors, we don't need FindClustering, because we have our own set of Idents
  Seurat_Obj = FindNeighbors(Seurat_Obj, dims = 1:n_dims_PCA)
  #Accordingly clustering the cells
  Seurat_Obj = FindClusters(Seurat_Obj, resolution = cluster_resolution)
  #Run tsne
  Seurat_Obj = RunTSNE(Seurat_Obj, dims = 1:n_dims_PCA, verbose = F)
  #TSNE plot
  TSNE_plot <- DimPlot(Seurat_Obj, reduction = "tsne", label= TRUE, pt.size=0.5)
  
  ggsave(plot = TSNE_plot,file = ".\\TSNE_SCRNA_Func.png", dpi=300, width=10, height=4.5)
  
  idents_flag = FALSE
  #Loading new idents names from user by letting them running FeaturePlot
  while (!idents_flag) {
    Genes_To_FeaturePlot = readline(prompt = "According to the tsne plot, which genes would you like to FeaturePlot (comma separated) ? :")
    Genes_To_FeaturePlot <- gsub(" ", "", Genes_To_FeaturePlot) #Removing whitespace from string
    Feature_plot <- FeaturePlot(Seurat_Obj, features = unlist(strsplit(Genes_To_FeaturePlot,",")), order=TRUE,pt.size=0.5, reduction="tsne")
    ggsave(plot = Feature_plot,file = ".\\Features_SCRNA_Func.png", dpi=300, width=10, height=4.5)
    
    
    decision = readline(prompt = "Have you decided what names to give to the tsne clusters according to the FeaturePlot ? (Yes/No) :")
    
    if ( decision == "Yes") {
      idents_flag = TRUE
    }
    
  }
  
  Idents_Names = readline(prompt = "Enter the new names of the clusters in order (from 0):")
  Idents_Names = unlist(strsplit(Idents_Names,","))
  Seurat_Obj@active.ident = Seurat_Obj$seurat_clusters
  
  for (i in 0:(length(Idents_Names) - 1)) {
    newIdent <- Idents_Names[i+1]
    names(newIdent) <- as.character(i)
    print(newIdent)
    Seurat_Obj = RenameIdents(Seurat_Obj, newIdent)
    
  }
  
  #Updating orig.ident in metadata - merging Idents(Seurat_Obj) with metadata by rownames(cell ids)
  Seurat_Obj@meta.data <- transform(merge(Seurat_Obj@meta.data,Idents(Seurat_Obj),by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  #Removing orig.ident column (it is old)
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data[,-which(colnames(Seurat_Obj@meta.data) == "orig.ident")]
  #Rename the colname y to orig.ident
  colnames(Seurat_Obj@meta.data)[colnames(Seurat_Obj@meta.data) == "y"] ="orig.ident"
  #Relocating the new orig.ident column to be the first column
  Seurat_Obj@meta.data <- Seurat_Obj@meta.data %>% relocate(orig.ident, .before = nCount_RNA)
  
  TSNE_Plot_Final <- DimPlot(Seurat_Obj, reduction = "tsne", label=FALSE, pt.size=0.5)
  ggsave(plot = TSNE_Plot_Final, file = ".\\TSNE_Named_SCRNA_Func.png", dpi=300, width=10, height=4.5)
  
  return (Seurat_Obj)
}


UpdateIdents <- function() {
  
  #Setting identity numbers - Intra = 1, Extra = 0
  
  #Intra df - cell_ID | Identity
  Intra_idents_df <- data.frame(Cell_IDs = c(names(Intra_tumoral_Intrst_neighbors_list)), Ident = 1)
  
  Intra_idents_df <- data.frame(Intra_idents_df[,-1], row.names=Intra_idents_df[,1])
  
  colnames(Intra_idents_df) <- "seurat_clusters"
  
  #Extra df - cell_ID | Identity
  Extra_idents_df <- data.frame(Cell_IDs = c(names(Extra_tumoral_Intrst_neighbors_list)), Ident = 0)
  
  Extra_idents_df <- data.frame(Extra_idents_df[,-1], row.names=Extra_idents_df[,1])
  
  colnames(Extra_idents_df) <- "seurat_clusters"
  
  #Binding both dfs
  Idents_df <- rbind(Intra_idents_df,Extra_idents_df)
  
  
  Subset_Cells@meta.data <- Subset_Cells@meta.data[,-which(colnames(Subset_Cells@meta.data) == "seurat_clusters")]
  Subset_Cells@meta.data <- transform(merge(Subset_Cells@meta.data,Idents_df,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  
  #Ordering of the cells is important because Idents(Seurat) is a factor object of (factor,integer) and it's not flexible
  order_of_cells <- as.data.frame(Idents(Subset_Cells))
  #Filtering for only cell_ids dataframe to get the order of the cell_ids in Idents
  order_of_cells <- tibble::rownames_to_column(order_of_cells, "cell_ids")
  order_of_cells <- as.data.frame(order_of_cells[,-2])
  colnames(order_of_cells) <- "cell_ids"
  
  #Reordering Idents_df according to Idents(Subset_Cells) order
  Idents_df <- Idents_df[order(match(rownames(Idents_df), order_of_cells$cell_ids)), , drop = FALSE]
  
  #SETTING IDENTS AS WE WANTED! :)
  Idents(Subset_Cells) <- Idents_df
  
  
  
  
}
