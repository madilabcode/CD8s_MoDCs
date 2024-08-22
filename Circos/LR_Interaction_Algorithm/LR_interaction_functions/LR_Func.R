
#*** create the Ligand_Receptor table ***#

createLRtable <- function(fromCol, toCol, num2compare="topQ",
                          LigandReceptorTable = LigandReceptorTable,
                          EXP.Ligands, EXP.Receptors, enable_overlap = F){
  require(tidyverse)  
  
  #few cutoff options for comparing ligands and receptors
    
  if(num2compare == "Median"){ #above median
    fromMedian <- fromCol[,1] %>% median()
    fromtop    <- fromCol %>% subset(. > fromMedian)
    toMedian   <- toCol[,1]  %>% median()
    totop      <- toCol %>% subset(. > toMedian)
    
  }else if(num2compare == "Mean"){ #above mean
    fromMean <- fromCol[,1] %>% mean()
    fromtop  <- fromCol %>% subset(. > fromMean)
    toMean   <- toCol[,1]   %>% mean()
    totop    <- toCol %>% subset(. > toMean)
    
  }else if(num2compare == "topQ"){ #top Quarter
    fromTQ  <- fromCol[,1] %>% quantile() %>% .["75%"]
    fromtop <- fromCol %>% subset(. > fromTQ)
    toTQ    <- toCol[,1]   %>% quantile() %>% .["75%"]
    totop   <- toCol %>% subset(. > toTQ)
    
  }else if(num2compare == "topDecile"){ #top Decile
    fromTD  <- fromCol[,1] %>% quantile(prob = seq(0, 1, length = 11)) %>% .["90%"]
    fromtop <- fromCol %>% subset(. > fromTD)
    toTD    <- toCol[,1]   %>% quantile(prob = seq(0, 1, length = 11)) %>% .["90%"]
    totop   <- toCol %>% subset(. > toTD)
    
  }else if(num2compare == "All"){ #all
    fromtop <- fromCol
    totop   <- toCol
    
  }else if(num2compare %>% as.numeric() %>% is.numeric()){
    fromtop <- as.data.frame(fromCol %>% head(num2compare %>% as.numeric()))
    totop   <- as.data.frame(toCol   %>% head(num2compare %>% as.numeric()))
  }
  
  #subset LR table to candidates from fromtop and totop
  LRSub1 <- LigandReceptorTable[which(LigandReceptorTable$from %in% rownames(fromtop)),]
  LRSub2 <- LRSub1[which(LRSub1$to %in% rownames(totop)),]
  
  # extract the expression of the ligands and insert to LR table
  LRSub3 <- inner_join(LRSub2, rownames_to_column(EXP.Ligands),   by = c('from' = "rowname"))
  
  # extract the expression of the receptors and insert to LR table
  LR <- inner_join(LRSub3, rownames_to_column(EXP.Receptors), by = c('to' = "rowname"))
  
  colnames(LR) <- c("Ligand","Receptor","LigandAbundance","ReceptorExp")
  LR <- LR %>% dplyr::filter(ReceptorExp > 0)
  
  #remove overlaps between ligand and receptor
  if(enable_overlap == F){
    LR <- LR[!(LR$Ligand   %in% base::intersect(LR$Ligand, LR$Receptor)),]
    LR <- LR[!(LR$Receptor %in% base::intersect(LR$Ligand, LR$Receptor)),]  #reciprocally
  }else{
    overlap_genes <- base::intersect(LR$Ligand, LR$Receptor)
    
    overlap_genes_i_L <- which(LR$Ligand %in% overlap_genes)
    LR[overlap_genes_i_L,'Ligand'] <- paste0(LR[overlap_genes_i_L,'Ligand'], "_L")
    
    overlap_genes_i_R <- which(LR$Receptor %in% overlap_genes)
    LR[overlap_genes_i_R,'Receptor'] <- paste0(LR[overlap_genes_i_R,'Receptor'], "_R")
  }
  
  return(LR %>% unique())
  print("Done LR")
}