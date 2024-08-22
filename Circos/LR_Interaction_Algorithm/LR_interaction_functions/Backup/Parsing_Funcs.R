LigandParser <- function(MS_data, cell_line, read_thresh){
  MS_dataSub <- MS_data[, c('Gene Symbol', grep(cell_line, MS_data %>% colnames, value = TRUE))] %>% na.omit()
  colnames(MS_dataSub) <- c('Gene Symbol', 'MS_Read')
  MS_dataSubOrd <- MS_dataSub %>% arrange(desc(MS_Read))
  MS_dataSubOrdUni <- MS_dataSubOrd[!(MS_dataSubOrd$`Gene Symbol` %>% duplicated()),]
  MS_dataSubOrdUni <- MS_dataSubOrdUni %>% filter(MS_Read > read_thresh) %>%
                        remove_rownames() %>% column_to_rownames('Gene Symbol') 
  return(MS_dataSubOrdUni)
}

ReceptorParser <- function(ReceptorDE_Res, cell_line, FC_thresh){
  #ReceptorDE_ResSub <- ReceptorDE_Res[,-grep("Significantly", ReceptorDE_Res %>% colnames())]
  #ReceptorDE_ResSub <- ReceptorDE_ResSub[, c('gene_name' , grep(cell_line, ReceptorDE_ResSub %>% colnames, value = TRUE))]
  #colnames(ReceptorDE_ResSub) <- c("gene_name", "log2FoldChange", "padj")
  ReceptorDE_ResSub <- ReceptorDE_Res
  ReceptorDE_ResSubOrd <- ReceptorDE_ResSub %>% arrange(desc(log2FoldChange))
  ReceptorDE_ResSubOrdUni <- ReceptorDE_ResSubOrd[!(ReceptorDE_ResSubOrd$gene_name %>% duplicated()),] %>% na.omit()
  DE_Receptors <- ReceptorDE_ResSubOrdUni %>% filter(log2FoldChange > FC_thresh & padj < 0.05) %>%
                        remove_rownames() %>% column_to_rownames('gene_name') 
  DE_Receptors <- DE_Receptors %>% select(log2FoldChange) 
  return(DE_Receptors)
}

RNAseqEXPParser <- function(RNAseqEXPMean, cell_line){
  RNAseqEXPMeanSub <- RNAseqEXPMean %>% select(c('gene_name',cell_line)) 
  RNAseqEXPMeanSubOrd <- RNAseqEXPMeanSub %>% arrange(desc(get(cell_line)))
  RNAseqEXPMeanSubOrd <- RNAseqEXPMeanSubOrd[!(RNAseqEXPMeanSubOrd$gene_name %>% duplicated()),] %>%
    remove_rownames() %>% column_to_rownames('gene_name') 
  return(RNAseqEXPMeanSubOrd)
}