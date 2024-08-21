# Kfir Inbal Generic Function R Code
#
# DESeq2 Analysis Mouse Genome
#
# This function converts EMSMBEL ids into gene names;
# then applies DESeq2 function of each sample vs control using WALD test;
# also, it applies DESeq2 function of all samples using LRT test;
# Then normalizes the counts for PCA plot;
# Plotting dispersions plot;
# Then does pairwise comparison between all samples and for each pair plots MA plot and Volcano Plot 
# and generates DGE table;
# Finally, prepares two dataframes ready for heatmap:
#   1) Dataframe of the normalized counts table with gene names of LRT DEG table.
#   2) Dataframe of the normalized counts table with gene names that are common to all pairwise comparisons.
#
# @params:
#   [I] Counts - assumes counts table with ENSEMBL gene ids as rows and samples as columns.
#   [II] Samples_Groups_By_Order - group classification for each sample according to their order in the counts table. e.g. c("Treatment","Treatment","Control","Control")
#   [III] Samples_Colors_By_Order - color classification for each sample according to their order in the coutns table. e.g. c("blue","blue","green","green")
#   [IV] Normalization Method - possible options: "rlog" or "vst". Default "rlog".
#   [V] plots_path - a path to where to save the plots and csv tables.
#   [VI] Experiment_Name - a name for the experiment, used for saving plot and table names. Default the current time and date.
#   [VII] PCA_width - PCA plot image width. Default 5.
#   [VIII] PCA_height - PCA plot image height. Default 4.
#
# @return: a list of all data generated in this function for reuse outside the function.



Mouse_DESeq2_Analysis <- function(Counts, Normalization_Method = "rlog",Samples_Groups_By_Order, Samples_Colors_By_Order, plots_path = ".\\", Experiment_Name = Sys.time(), PCA_width = 10, PCA_height = 7, Comparisons = TRUE, ref = TRUE) {
  #############Preparing Metadata ##################
  
  Counts_metaData <- data.frame(row.names = colnames(Counts), group = factor(Samples_Groups_By_Order), color = Samples_Colors_By_Order) #Assumptions - control group is called 'Control' ; the control group is located after all treatment groups
  
  ############# Converting ensembl gene IDs to Gene Name ##########################
  Counts <- tibble::rownames_to_column(Counts, "Gene")
  
  #Annotate the Ensembl gene IDs to gene symbols:
  if (!exists("ensembl")) {   
    ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
  }
  
  annot <- getBM(
    attributes = c(
      'mgi_symbol',
      'ensembl_gene_id'),
    filters = 'ensembl_gene_id',
    values = Counts$Gene,
    mart = ensembl)
  
  annot <- annot[!(is.na(annot$mgi_symbol) | annot$mgi_symbol == ""), ]
  
  Counts <- merge(
    x = Counts,
    y =  annot,
    by.y = 'ensembl_gene_id',
    by.x = "Gene")
  
  
  
  #Removing column "Gene" which holds ensembl ids
  Counts <- Counts[,-which(colnames(Counts) %in% c("Gene"))]
  #Moving the column "mgi_symbol" to be the first column
  Counts <- Counts %>%
    relocate(mgi_symbol, .before = 1)
  #Renaming the column "mgi_symbol" to "Gene"
  names(Counts)[names(Counts) == 'mgi_symbol'] <- 'Gene'
  
  #Creating column "sumExpression" of all columns
  Counts$sumExpression <- rowSums(Counts[,-1])
  
  #Select duplicate Gene names by max of "sumExpression"
  Counts <- Counts %>% 
    group_by(Gene) %>% 
    slice_max(sumExpression, n = 1,with_ties = FALSE) %>%
    ungroup()
  
  #Removing column "sumExpression"
  Counts <- subset(Counts, select = -sumExpression)
  
  
  #First column to row.names
  Counts <- Counts %>% remove_rownames %>% column_to_rownames(var="Gene")
  
  #Make sure there are no duplicates (Should return 0 if there are no duplicates, >=1 if there are)
  sum(duplicated(row.names(Counts)))
  ############################################################################
  
  ########## Calculating gene length for each gene ##########################
  
  #Calculating gene length
  gene_coords = getBM(attributes=c("mgi_symbol","ensembl_gene_id", "start_position","end_position"), filters="mgi_symbol", values=row.names(Counts), mart=ensembl)
  
  
  sum(duplicated((gene_coords$ensembl_gene_id)))
  
  gene_coords <- gene_coords[!duplicated(gene_coords$mgi_symbol),]
  all(gene_coords$mgi_symbol %in% row.names(Counts))
  all(row.names(Counts)%in% gene_coords$mgi_symbol)
  
  #First column to row.names
  gene_coords <- gene_coords %>% remove_rownames %>% column_to_rownames(var="mgi_symbol")
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  gene_coords
  
  sum(duplicated(row.names(gene_coords)))
  
  gene_coords$effLength <- gene_coords$size - 203.7 + 1
  
  ############################################################################
  ###############################DESEQ 2#######################################
  
  #Make sure the row namees in metaData matches to col names in raw_data
  all(colnames(Counts) %in% rownames(Counts_metaData))
  
  #Make sure they are in the same order
  all(colnames(Counts) == rownames(Counts_metaData))
  
  
  Counts_dds <- DESeqDataSetFromMatrix(countData=Counts, 
                                       colData=Counts_metaData, 
                                       design=~group)
  
  #pre-filtering: removing rows with low gene counts
  #keeping rows that have at least 10 reads total
  #keep <- rowSums(counts(Counts_dds)) >= 10
  Counts_dds <- Counts_dds[which(rowSums(counts(Counts_dds)) >= 10), ]
  
  #set the factor level
  if (ref == TRUE) {
    Counts_dds$group <- relevel(Counts_dds$group, ref = "Control") #Without ref it will just order it alphabetically. the untreated group(control) should be first, that's the reference it compares it to
  }
  
  if (Normalization_Method == "rlog") {
    Counts_nmdata <- rlog(Counts_dds, blind = FALSE)
    #PCA Plot
    PCA_p <- plotPCA(Counts_nmdata, intgroup = "group") + geom_text_repel(aes(label=name)) + scale_color_manual(values = unique(Counts_metaData$color)) + theme_classic() + coord_fixed(xlim(-20,40)) #hjust="bottom"
    
    ggsave(plot = PCA_p,file = paste0(plots_path,"\\PCA_Plot_",Experiment_Name,"_Deseq_rlog",".png"), dpi=300, width=PCA_width, height=PCA_height)
    
  } else if (Normalization_Method == "vst") {
    Counts_nmdata <- vst(Counts_dds, blind = FALSE)
    #PCA Plot
    PCA_p <- plotPCA(Counts_nmdata, intgroup = "group") + geom_text_repel(aes(label=name)) + scale_color_manual(values = unique(Counts_metaData$color)) + theme_classic() + coord_fixed(xlim(-20,40)) #hjust="bottom"
    
    ggsave(plot = PCA_p,file = paste0(plots_path,"\\PCA_Plot_",Experiment_Name,"_Deseq_vst",".png"), dpi=300, width=PCA_width, height=PCA_height)
  }
  
  
  #Run DESeq
  Counts_dds_Wald <- DESeq(Counts_dds)
  
  #Run DESeq LRT
  Counts_dds_LRT <- DESeq(Counts_dds, test = "LRT", reduced = ~ 1)
  deg_res_LRT = results(Counts_dds_LRT)
  deg_res_LRT <- na.omit(deg_res_LRT)
  deg_res_sigs_LRT <- deg_res_LRT[deg_res_LRT$padj < 0.05, ]
  
  #"#672A70","#138808"
  
  #Dispersions
  jpeg(paste0(plots_path,"\\Dispersion_Plot_",Experiment_Name,"_Deseq",".jpg"))
  Disp_p <- plotDispEsts(Counts_dds_Wald)
  dev.off()
  
  if(Comparisons == FALSE) {
    ################################# HEATMAP########################################
    #Counts_VST_data@assays@data@listData[[1]] ????? CHECK 
    #LRT Heatmap
    df_for_Heatmap <- deg_res_sigs_LRT[order(deg_res_sigs_LRT$log2FoldChange, decreasing = TRUE),]
    #Obtaining normalized counts
    nm_counts_deseq_LRT <- rlog(Counts_dds_LRT, blind=FALSE)
    #Filtering to only genes from the DEG sigs
    mat_LRT <- assay(nm_counts_deseq_LRT)[rownames(df_for_Heatmap), ]
    
    
    list_to_return <- list(Metadata = Counts_metaData,Raw_Counts = Counts,Counts_Normalized_data = list(data = Counts_nmdata,Normalization_Method = Normalization_Method, DESeq_Object_Type = "On Counts_dds"), Gene_Lengths = gene_coords, dds = Counts_dds ,dds_LRT = Counts_dds_LRT ,dds_Wald = Counts_dds_Wald, PCA_plot = PCA_p ,Disp_plot = Disp_p, LRT = deg_res_sigs_LRT, Heatmap_df_LRT = mat_LRT)
    
    saveRDS(object = list_to_return,file = paste0(plots_path,"\\",Experiment_Name,".rds"))
    
    return(list_to_return)
    
  }
  
  for_merge_list <- list()
  Counts_res_list <- list()
  grps_done <- list()
  for(grp in unique(Counts_metaData[,1])) {
    for(grp2 in unique(Counts_metaData[,1])) {
      if (!(grp2 %in% grps_done) & grp != grp2) {
        #Counts_res <- paste("Counts_res", g, sep = "_")
        #assign(Counts_res, results(Counts_dds, contrast = c(colnames(Counts_metaData)[1], g, "Control")))
        Counts_res <- results(Counts_dds_Wald, contrast = c(colnames(Counts_metaData)[1], grp, grp2))
        #Summary
        summary(Counts_res)
        #let's look at the results table
        head(Counts_res) 
        #baseMean = average of the normalized counts taken over all the samples. 
        #The lower the baseMean, the noisier it's gonna be.
        #log2FoldChange = positive - upregulated in treated vs control ; 
        #negative - downregulated in treated vs control.
        #lfcSE = standard error estimates for the log2FC.
        #stat = Wald test value for each gene
        Counts_res <- na.omit(Counts_res)
        summary(Counts_res)
        
        
        #MA PLOT
        jpeg(paste0(plots_path,"\\MA_Plot_",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_Deseq",".jpg"))
        MA_p <- plotMA(Counts_res)
        dev.off()
        
        #Plot the most basic volcano plot
        #For the most basic volcano plot, only a single data-frame, data-matrix, or tibble of test results is required, containing point labels, log2FC, and adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
        
        EVolcano_Plot <- EnhancedVolcano(Counts_res,
                                         lab = rownames(Counts_res),
                                         labSize = 0.5,
                                         x = 'log2FoldChange',
                                         y = 'pvalue',
                                         xlim = c(-8,8)) + labs(title = paste0(grp," Versus ", grp2, " Volcano Plot"),
                                                                subtitle = "Default cutoffs: pvalue-cutoff 10e-6; log2FC-cutoff > |2|")+
          theme(plot.title = element_text(size = 10),
                plot.subtitle = element_text(size = 5)) + theme_classic() + geom_segment(aes(x = 2, y = -log10(min(Counts_res$pvalue)), xend = 5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15, length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp), x=(0 + 6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA) + geom_segment(aes(x = -2, y = -log10(min(Counts_res$pvalue)), xend = -5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15,length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp2), x=(0 + -6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA)
        
        ggsave(plot = EVolcano_Plot,file = paste0(plots_path,"\\Volcano_General_",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_Deseq",".png"), dpi=300, width=7, height=5)
        
        
        Counts_res_logFC_gt2 <- Counts_res[Counts_res$log2FoldChange > 2, ]
        Counts_res_logFC_lt_neg2 <- Counts_res[Counts_res$log2FoldChange < -2, ]
        
        top_10_res_genes_up <- c(row.names(head(Counts_res_logFC_gt2[order(Counts_res_logFC_gt2$pvalue),], 10)))
        top_10_res_genes_down <- c(row.names(head(Counts_res_logFC_lt_neg2[order(Counts_res_logFC_lt_neg2$pvalue),], 10)))
        
        
        print(paste("UP:",top_10_res_genes_up))
        print(paste("DOWN:",top_10_res_genes_down))
        
        EVolcano_Plot_Specific_Labels <- EnhancedVolcano(Counts_res,
                                                         lab = rownames(Counts_res),
                                                         selectLab = c("Il27","Il1b","Ifi205","Ifit3","Cd40","Tnfsf15","Ccl5","Ifit3b","Cxcl10","Cxcl2","Ifit2","Tnf","Il12b","Il12a","Isg20","Socs3","Ccl4","Cxcl1", top_10_res_genes_up, top_10_res_genes_down),
                                                         labSize = 2,
                                                         x = 'log2FoldChange',
                                                         y = 'pvalue',
                                                         xlim = c(-8,8),
                                                         drawConnectors = TRUE) +
          labs(title = paste0(grp," Versus ",grp2, " Volcano Plot"),
               subtitle = "Default cutoffs: pvalue-cutoff 10e-6; log2FC-cutoff > |2|") +
          theme(plot.title = element_text(size = 10),
                plot.subtitle = element_text(size = 5)) + theme_classic() + geom_segment(aes(x = 2, y = -log10(min(Counts_res$pvalue)), xend = 5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15, length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp), x=(0 + 6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA) + geom_segment(aes(x = -2, y = -log10(min(Counts_res$pvalue)), xend = -5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15, length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp2), x=(0 + -6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA) 
        
        
        
        ggsave(plot = EVolcano_Plot_Specific_Labels,file = paste0(plots_path,"\\Volcano_Spcfc_Labels",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_Deseq",".png"), dpi=300, width=7, height=5)
        
        colnames(Counts_res)[colnames(Counts_res) == "log2FoldChange"] = paste0("lgFC_", grp, "_vs_",grp2)
        #Filtering out insignificant DEG
        deg_res_sigs <- Counts_res[Counts_res$padj < 0.05, ]
        summary(deg_res_sigs)
        
        #Sort summary list by p-value
        deg_res_sigs <- deg_res_sigs[order(deg_res_sigs$padj),]
        write.csv(deg_res_sigs, file = paste0(plots_path,"\\DEG_",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_padj_lt_0.05.csv"))
        
        for_merge_list <- c(for_merge_list, list(c(deg_res_sigs)))
        
        data_list <- list(c(deg_res = Counts_res, deg_res_sig = deg_res_sigs, MA_plot = list(MA_p), Disp_plot = list(Disp_p), Volcano_Plot_General = list(EVolcano_Plot), Volcano_P_Specific_Labels =  list(EVolcano_Plot_Specific_Labels)))
        
        Counts_res_list <- c(Counts_res_list, `names<-`(data_list, paste(grp,"_vs_",grp2)))
        grps_done <- append(grps_done,grp)
        
        
      }
    }
  }
  
  
  
  if (length(unique(Counts_metaData[,1])) > 2 ) {
    
    deg_res_sigs_Merged <- Reduce(function(x, y) merge(x, y, by='row.names'), for_merge_list, accumulate=F)
    deg_res_sigs_Merged <- as.data.frame(deg_res_sigs_Merged)
    #First column to row.names
    deg_res_sigs_Merged <- deg_res_sigs_Merged %>% remove_rownames %>% column_to_rownames(var="Row.names")
    write.csv(deg_res_sigs_Merged, file = paste0(plots_path,"\\LogFC_DEG_",Experiment_Name,"_--Merged","--_padj_lt_0.05.csv"))
    
    #Merged Heatmap
    df_for_Heatmap <- deg_res_sigs_Merged[order(deg_res_sigs_Merged[, 2], decreasing = TRUE),]
    #Obtaining normalized counts
    nm_counts_deseq <- rlog(Counts_dds_Wald, blind=FALSE)
    #Filtering to only genes from the DEG sigs
    mat_merged <- assay(nm_counts_deseq)[rownames(df_for_Heatmap), ]
    
    
  } else if (length(unique(Counts_metaData[,1])) == 2 ) {
    deg_res_sigs_Merged <- data.frame()
    mat_merged <- data.frame()
  }
  
  ################################# HEATMAP########################################
  #Counts_VST_data@assays@data@listData[[1]] ????? CHECK 
  #LRT Heatmap
  df_for_Heatmap <- deg_res_sigs_LRT[order(deg_res_sigs_LRT$log2FoldChange, decreasing = TRUE),]
  #Obtaining normalized counts
  nm_counts_deseq_LRT <- rlog(Counts_dds_LRT, blind=FALSE)
  #Filtering to only genes from the DEG sigs
  mat_LRT <- assay(nm_counts_deseq_LRT)[rownames(df_for_Heatmap), ]
  
  
  list_to_return <- list(Metadata = Counts_metaData,Raw_Counts = Counts,Counts_Normalized_data = list(data = Counts_nmdata,Normalization_Method = Normalization_Method, DESeq_Object_Type = "On Counts_dds"), Gene_Lengths = gene_coords,dds = Counts_dds, dds_LRT = Counts_dds_LRT ,dds_Wald = Counts_dds_Wald, PCA_plot = PCA_p ,Disp_plot = Disp_p, Counts_res_list = Counts_res_list, forMerge = for_merge_list, merged = deg_res_sigs_Merged, LRT = deg_res_sigs_LRT, Heatmap_df_LRT = mat_LRT, Heatmap_df_merged = mat_merged)
  
  saveRDS(object = list_to_return,file = paste0(plots_path,"\\",Experiment_Name,".rds"))
  
  return(list_to_return)
  
}























Human_DESeq2_Analysis <- function(Counts, Normalization_Method = "rlog",Samples_Groups_By_Order, Samples_Colors_By_Order, plots_path = ".\\", Experiment_Name = Sys.time(), PCA_width = 5, PCA_height = 4, Comparisons = TRUE, ref = TRUE) {
  #############Preparing Metadata ##################
  
  Counts_metaData <- data.frame(row.names = colnames(Counts), group = factor(Samples_Groups_By_Order), color = Samples_Colors_By_Order) #Assumptions - control group is called 'Control' ; the control group is located after all treatment groups
  
  ############# Converting ensembl gene IDs to Gene Name ##########################
  Counts <- tibble::rownames_to_column(Counts, "Gene")
  
  #Annotate the Ensembl gene IDs to gene symbols:
  if (!exists("ensembl")) {   
    ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
  }
  
  annot <- getBM(
    attributes = c(
      'hgnc_symbol',
      'ensembl_gene_id'),
    filters = 'ensembl_gene_id',
    values = Counts$Gene,
    mart = ensembl)
  
  annot <- annot[!(is.na(annot$hgnc_symbol) | annot$hgnc_symbol == ""), ]
  
  Counts <- merge(
    x = Counts,
    y =  annot,
    by.y = 'ensembl_gene_id',
    by.x = "Gene")
  
  
  
  #Removing column "Gene" which holds ensembl ids
  Counts <- Counts[,-which(colnames(Counts) %in% c("Gene"))]
  #Moving the column "hgnc_symbol" to be the first column
  Counts <- Counts %>%
    relocate(hgnc_symbol, .before = 1)
  #Renaming the column "hgnc_symbol" to "Gene"
  names(Counts)[names(Counts) == 'hgnc_symbol'] <- 'Gene'
  
  #Creating column "sumExpression" of all columns
  Counts$sumExpression <- rowSums(Counts[,-1])
  
  #Select duplicate Gene names by max of "sumExpression"
  Counts <- Counts %>% 
    group_by(Gene) %>% 
    slice_max(sumExpression, n = 1,with_ties = FALSE) %>%
    ungroup()
  
  #Removing column "sumExpression"
  Counts <- subset(Counts, select = -sumExpression)
  
  
  #First column to row.names
  Counts <- Counts %>% remove_rownames %>% column_to_rownames(var="Gene")
  
  #Make sure there are no duplicates (Should return 0 if there are no duplicates, >=1 if there are)
  sum(duplicated(row.names(Counts)))
  ############################################################################
  
  ########## Calculating gene length for each gene ##########################
  
  #Calculating gene length
  gene_coords = getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "start_position","end_position"), filters="hgnc_symbol", values=row.names(Counts), mart=ensembl)
  
  
  sum(duplicated((gene_coords$ensembl_gene_id)))
  
  gene_coords <- gene_coords[!duplicated(gene_coords$hgnc_symbol),]
  all(gene_coords$hgnc_symbol %in% row.names(Counts))
  all(row.names(Counts)%in% gene_coords$hgnc_symbol)
  
  #First column to row.names
  gene_coords <- gene_coords %>% remove_rownames %>% column_to_rownames(var="hgnc_symbol")
  gene_coords$size=gene_coords$end_position - gene_coords$start_position
  gene_coords
  
  sum(duplicated(row.names(gene_coords)))
  
  gene_coords$effLength <- gene_coords$size - 203.7 + 1
  
  ############################################################################
  ###############################DESEQ 2#######################################
  
  #Make sure the row namees in metaData matches to col names in raw_data
  all(colnames(Counts) %in% rownames(Counts_metaData))
  
  #Make sure they are in the same order
  all(colnames(Counts) == rownames(Counts_metaData))
  
  
  Counts_dds <- DESeqDataSetFromMatrix(countData=Counts, 
                                       colData=Counts_metaData, 
                                       design=~group)
  
  #pre-filtering: removing rows with low gene counts
  #keeping rows that have at least 10 reads total
  #keep <- rowSums(counts(Counts_dds)) >= 10
  Counts_dds <- Counts_dds[which(rowSums(counts(Counts_dds)) >= 10), ]
  
  #set the factor level
  if (ref == TRUE) {
    Counts_dds$group <- relevel(Counts_dds$group, ref = "Control") #Without ref it will just order it alphabetically. the untreated group(control) should be first, that's the reference it compares it to
  }
  
  if (Normalization_Method == "rlog") {
    Counts_nmdata <- rlog(Counts_dds, blind = FALSE)
    #PCA Plot
    PCA_p <- plotPCA(Counts_nmdata, intgroup = "group") + geom_text_repel(aes(label=name)) + scale_color_manual(values = unique(Counts_metaData$color)) + theme_classic() + coord_fixed(xlim(-20,40)) #hjust="bottom"
    
    ggsave(plot = PCA_p,file = paste0(plots_path,"\\PCA_Plot_",Experiment_Name,"_Deseq_rlog",".png"), dpi=300, width=PCA_width, height=PCA_height)
    
  } else if (Normalization_Method == "vst") {
    Counts_nmdata <- vst(Counts_dds, blind = FALSE)
    #PCA Plot
    PCA_p <- plotPCA(Counts_nmdata, intgroup = "group") + geom_text_repel(aes(label=name)) + scale_color_manual(values = unique(Counts_metaData$color)) + theme_classic() + coord_fixed(xlim(-20,40)) #hjust="bottom"
    
    ggsave(plot = PCA_p,file = paste0(plots_path,"\\PCA_Plot_",Experiment_Name,"_Deseq_vst",".png"), dpi=300, width=PCA_width, height=PCA_height)
  }
  
  
  #Run DESeq
  Counts_dds_Wald <- DESeq(Counts_dds)
  
  #Run DESeq LRT
  Counts_dds_LRT <- DESeq(Counts_dds, test = "LRT", reduced = ~ 1)
  deg_res_LRT = results(Counts_dds_LRT)
  deg_res_LRT <- na.omit(deg_res_LRT)
  deg_res_sigs_LRT <- deg_res_LRT[deg_res_LRT$padj < 0.05, ]
  
  #"#672A70","#138808"
  
  #Dispersions
  jpeg(paste0(plots_path,"\\Dispersion_Plot_",Experiment_Name,"_Deseq",".jpg"))
  Disp_p <- plotDispEsts(Counts_dds_Wald)
  dev.off()
  
  if(Comparisons == FALSE) {
    ################################# HEATMAP########################################
    #Counts_VST_data@assays@data@listData[[1]] ????? CHECK 
    #LRT Heatmap
    df_for_Heatmap <- deg_res_sigs_LRT[order(deg_res_sigs_LRT$log2FoldChange, decreasing = TRUE),]
    #Obtaining normalized counts
    nm_counts_deseq_LRT <- rlog(Counts_dds_LRT, blind=FALSE)
    #Filtering to only genes from the DEG sigs
    mat_LRT <- assay(nm_counts_deseq_LRT)[rownames(df_for_Heatmap), ]
    
    
    list_to_return <- list(Metadata = Counts_metaData,Raw_Counts = Counts,Counts_Normalized_data = list(data = Counts_nmdata,Normalization_Method = Normalization_Method, DESeq_Object_Type = "On Counts_dds"), Gene_Lengths = gene_coords, dds = Counts_dds ,dds_LRT = Counts_dds_LRT ,dds_Wald = Counts_dds_Wald, PCA_plot = PCA_p ,Disp_plot = Disp_p, LRT = deg_res_sigs_LRT, Heatmap_df_LRT = mat_LRT)
    
    saveRDS(object = list_to_return,file = paste0(plots_path,"\\",Experiment_Name,".rds"))
    
    return(list_to_return)
    
  }
  
  for_merge_list <- list()
  Counts_res_list <- list()
  grps_done <- list()
  for(grp in unique(Counts_metaData[,1])) {
    for(grp2 in unique(Counts_metaData[,1])) {
      if (!(grp2 %in% grps_done) & grp != grp2) {
        #Counts_res <- paste("Counts_res", g, sep = "_")
        #assign(Counts_res, results(Counts_dds, contrast = c(colnames(Counts_metaData)[1], g, "Control")))
        Counts_res <- results(Counts_dds_Wald, contrast = c(colnames(Counts_metaData)[1], grp, grp2))
        #Summary
        summary(Counts_res)
        #let's look at the results table
        head(Counts_res) 
        #baseMean = average of the normalized counts taken over all the samples. 
        #The lower the baseMean, the noisier it's gonna be.
        #log2FoldChange = positive - upregulated in treated vs control ; 
        #negative - downregulated in treated vs control.
        #lfcSE = standard error estimates for the log2FC.
        #stat = Wald test value for each gene
        Counts_res <- na.omit(Counts_res)
        summary(Counts_res)
        
        
        #MA PLOT
        jpeg(paste0(plots_path,"\\MA_Plot_",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_Deseq",".jpg"))
        MA_p <- plotMA(Counts_res)
        dev.off()
        
        #Plot the most basic volcano plot
        #For the most basic volcano plot, only a single data-frame, data-matrix, or tibble of test results is required, containing point labels, log2FC, and adjusted or unadjusted P values. The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.
        
        EVolcano_Plot <- EnhancedVolcano(Counts_res,
                                         lab = rownames(Counts_res),
                                         labSize = 0.5,
                                         x = 'log2FoldChange',
                                         y = 'pvalue',
                                         xlim = c(-8,8)) + labs(title = paste0(grp," Versus ", grp2, " Volcano Plot"),
                                                                subtitle = "Default cutoffs: pvalue-cutoff 10e-6; log2FC-cutoff > |2|")+
          theme(plot.title = element_text(size = 10),
                plot.subtitle = element_text(size = 5)) + theme_classic() + geom_segment(aes(x = 2, y = -log10(min(Counts_res$pvalue)), xend = 5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15, length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp), x=(0 + 6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA) + geom_segment(aes(x = -2, y = -log10(min(Counts_res$pvalue)), xend = -5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15,length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp2), x=(0 + -6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA)
        
        ggsave(plot = EVolcano_Plot,file = paste0(plots_path,"\\Volcano_General_",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_Deseq",".png"), dpi=300, width=7, height=5)
        
        
        Counts_res_logFC_gt2 <- Counts_res[Counts_res$log2FoldChange > 2, ]
        Counts_res_logFC_lt_neg2 <- Counts_res[Counts_res$log2FoldChange < -2, ]
        
        top_10_res_genes_up <- c(row.names(head(Counts_res_logFC_gt2[order(Counts_res_logFC_gt2$pvalue),], 10)))
        top_10_res_genes_down <- c(row.names(head(Counts_res_logFC_lt_neg2[order(Counts_res_logFC_lt_neg2$pvalue),], 10)))
        
        
        print(paste("UP:",top_10_res_genes_up))
        print(paste("DOWN:",top_10_res_genes_down))
        
        EVolcano_Plot_Specific_Labels <- EnhancedVolcano(Counts_res,
                                                         lab = rownames(Counts_res),
                                                         selectLab = c("Il27","Il1b","Ifi205","Ifit3","Cd40","Tnfsf15","Ccl5","Ifit3b","Cxcl10","Cxcl2","Ifit2","Tnf","Il12b","Il12a","Isg20","Socs3","Ccl4","Cxcl1", top_10_res_genes_up, top_10_res_genes_down),
                                                         labSize = 2,
                                                         x = 'log2FoldChange',
                                                         y = 'pvalue',
                                                         xlim = c(-8,8),
                                                         drawConnectors = TRUE) +
          labs(title = paste0(grp," Versus ",grp2, " Volcano Plot"),
               subtitle = "Default cutoffs: pvalue-cutoff 10e-6; log2FC-cutoff > |2|") +
          theme(plot.title = element_text(size = 10),
                plot.subtitle = element_text(size = 5)) + theme_classic() + geom_segment(aes(x = 2, y = -log10(min(Counts_res$pvalue)), xend = 5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15, length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp), x=(0 + 6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA) + geom_segment(aes(x = -2, y = -log10(min(Counts_res$pvalue)), xend = -5, yend = -log10(min(Counts_res$pvalue))), arrow = arrow(angle = 15, length = unit(0.5, "cm"))) + geom_label(aes(label=paste("Up In ",grp2), x=(0 + -6)/2), y = -log10(min(Counts_res$pvalue)) + 5, size=2, label.size = NA, fill = NA) 
        
        
        
        ggsave(plot = EVolcano_Plot_Specific_Labels,file = paste0(plots_path,"\\Volcano_Spcfc_Labels",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_Deseq",".png"), dpi=300, width=7, height=5)
        
        colnames(Counts_res)[colnames(Counts_res) == "log2FoldChange"] = paste0("lgFC_", grp, "_vs_",grp2)
        #Filtering out insignificant DEG
        deg_res_sigs <- Counts_res[Counts_res$padj < 0.05, ]
        summary(deg_res_sigs)
        
        #Sort summary list by p-value
        deg_res_sigs <- deg_res_sigs[order(deg_res_sigs$padj),]
        write.csv(deg_res_sigs, file = paste0(plots_path,"\\DEG_",Experiment_Name,"_--",grp,"_Versus_",grp2,"--_padj_lt_0.05.csv"))
        
        for_merge_list <- c(for_merge_list, list(c(deg_res_sigs)))
        
        data_list <- list(c(deg_res = Counts_res, deg_res_sig = deg_res_sigs, MA_plot = list(MA_p), Disp_plot = list(Disp_p), Volcano_Plot_General = list(EVolcano_Plot), Volcano_P_Specific_Labels =  list(EVolcano_Plot_Specific_Labels)))
        
        Counts_res_list <- c(Counts_res_list, `names<-`(data_list, paste(grp,"_vs_",grp2)))
        grps_done <- append(grps_done,grp)
        
        
      }
    }
  }
  
  
  
  if (length(unique(Counts_metaData[,1])) > 2 ) {
    
    deg_res_sigs_Merged <- Reduce(function(x, y) merge(x, y, by='row.names'), for_merge_list, accumulate=F)
    deg_res_sigs_Merged <- as.data.frame(deg_res_sigs_Merged)
    #First column to row.names
    deg_res_sigs_Merged <- deg_res_sigs_Merged %>% remove_rownames %>% column_to_rownames(var="Row.names")
    write.csv(deg_res_sigs_Merged, file = paste0(plots_path,"\\LogFC_DEG_",Experiment_Name,"_--Merged","--_padj_lt_0.05.csv"))
    
    #Merged Heatmap
    df_for_Heatmap <- deg_res_sigs_Merged[order(deg_res_sigs_Merged[, 2], decreasing = TRUE),]
    #Obtaining normalized counts
    nm_counts_deseq <- rlog(Counts_dds_Wald, blind=FALSE)
    #Filtering to only genes from the DEG sigs
    mat_merged <- assay(nm_counts_deseq)[rownames(df_for_Heatmap), ]
    
    
  } else if (length(unique(Counts_metaData[,1])) == 2 ) {
    deg_res_sigs_Merged <- data.frame()
    mat_merged <- data.frame()
  }
  
  ################################# HEATMAP########################################
  #Counts_VST_data@assays@data@listData[[1]] ????? CHECK 
  #LRT Heatmap
  df_for_Heatmap <- deg_res_sigs_LRT[order(deg_res_sigs_LRT$log2FoldChange, decreasing = TRUE),]
  #Obtaining normalized counts
  nm_counts_deseq_LRT <- rlog(Counts_dds_LRT, blind=FALSE)
  #Filtering to only genes from the DEG sigs
  mat_LRT <- assay(nm_counts_deseq_LRT)[rownames(df_for_Heatmap), ]
  
  
  list_to_return <- list(Metadata = Counts_metaData,Raw_Counts = Counts,Counts_Normalized_data = list(data = Counts_nmdata,Normalization_Method = Normalization_Method, DESeq_Object_Type = "On Counts_dds"), Gene_Lengths = gene_coords,dds = Counts_dds, dds_LRT = Counts_dds_LRT ,dds_Wald = Counts_dds_Wald, PCA_plot = PCA_p ,Disp_plot = Disp_p, Counts_res_list = Counts_res_list, forMerge = for_merge_list, merged = deg_res_sigs_Merged, LRT = deg_res_sigs_LRT, Heatmap_df_LRT = mat_LRT, Heatmap_df_merged = mat_merged)
  
  saveRDS(object = list_to_return,file = paste0(plots_path,"\\",Experiment_Name,".rds"))
  
  return(list_to_return)
  
}

