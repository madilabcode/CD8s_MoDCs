FClegenedNcolor <- function(LR, toCol){
  #create the color scheme and legend DE pval
  DE_Data_R <- toCol %>% dplyr::filter(toCol %>% row.names() %in% LR$Receptor)
  Max <- DE_Data_R$log2FoldChange %>% max(na.rm = T) %>% round(2)
  DEcol_fun = colorRamp2(c(0, Max), c("white", "orange"))
  lgd_DE = Legend(at =   c(0, Max),
                  col_fun = DEcol_fun,
                  title_position = "leftcenter-rot",
                  title = paste0("log FC Rec."),
                  legend_height = unit(2, "cm"))
  print("Done FC")
  return(list(DEcol_fun,lgd_DE,DE_Data_R))
}

FClegenedNcolor_Ligand <- function(LR, fromCol){
  #create the color scheme and legend DE pval
  DE_Data_R_L <- fromCol %>% dplyr::filter(fromCol %>% row.names() %in% LR$Ligand)
  Max <- DE_Data_R_L$log2FoldChange %>% max(na.rm = T) %>% round(2)
  DEcol_fun_L = colorRamp2(c(0, Max), c("white", "#28a10a"))
  lgd_DE_L = Legend(at =   c(0, Max),
                  col_fun = DEcol_fun_L,
                  title_position = "leftcenter-rot",
                  title = paste0("log FC Lig."),
                  legend_height = unit(2, "cm"))
  print("Done FC")
  return(list(DEcol_fun_L,lgd_DE_L,DE_Data_R_L))
}

DSAcalculator <- function(Gene ,RNAseqEXPCol, threshold = 100, directed_mouse_ppi){
  ReceptorTB <- directed_mouse_ppi[which(directed_mouse_ppi$`Input-node Gene Symbol` == Gene),]
  row.names(ReceptorTB) <- ReceptorTB$`Output-node Gene Symbol`
  ensName <- ProtNamesInfo[which(ProtNamesInfo$preferred_name == Gene),"protein_external_id"]  
  ReceptorTBString <- ProtActions[which(ProtActions$item_id_b == ensName),]
  ReceptorTBActive <- ReceptorTBString[which(ReceptorTBString$mode == "activation"),] %>% subset(.,score>threshold)
  if(ReceptorTBActive$mode %>% length() > 0){
    ReceptorTBActive$item_id_b <- Gene
    ReceptorTBActive$item_id_a <- ReceptorTBActive$item_id_a %>% as.character()
    for(ensName in ReceptorTBActive$item_id_a){
      ReceptorTBActive[which(ReceptorTBActive$item_id_a == ensName),"item_id_a"] <- ENSNamesInfo[ensName,"preferred_name"]}
    inter <- base::intersect(ReceptorTB$`Output-node Gene Symbol`, ReceptorTBActive$item_id_a) %>% as.vector()
    if(inter %>% length() > 0){
      DSA <- RNAseqEXPCol[base::intersect(RNAseqEXPCol %>% rownames, inter),] %>% mean %>% round(1)
      n <- length(inter)
      return(list('DSA_Score' = DSA, 'Number' = n, 'Gene_Names' = inter))
    }else{
      return(NA)}
  }else{
    return(NA)}
}
"
DSAcalculator <- function(Gene, RNAseqEXPCol, threshold = 100, LigandReceptorTable){
  # can the receptor act as a ligand?
  if(intersect(Gene, LigandReceptorTable$from) %>% length() > 0){
    ReceptorTB <- LigandReceptorTable[which(LigandReceptorTable$from == Gene),] %>% unique() 
    ensName <- ProtNamesInfo[which(ProtNamesInfo$preferred_name == Gene),'protein_external_id']  
    ReceptorTBActive <- ProtActions %>% subset(item_id_a == ensName & mode == 'activation' & score > threshold)
    # are there receptors of mode activation
    if(ReceptorTBActive %>% nrow() > 0){
      ReceptorTBActive$item_id_a <- Gene
      # translate ligand names from ENS to symbol
      ReceptorTBActive$item_id_b <- ReceptorTBActive$item_id_b %>% as.character()
      mrg <- merge(x = ReceptorTBActive, y = ENSNamesInfo, by.x = 'item_id_b', by.y = 0)
      # which genes are receptors to Gene from activation subset DF ReceptorTBActive
      inter <- intersect(ReceptorTB$to, mrg$preferred_name) %>% as.vector()
      #the they 
      if(inter %>% length() > 0){
          DSA <- RNAseqEXPCol[intersect(RNAseqEXPCol %>% rownames, inter),] %>% mean %>% round(3)
          n <- length(inter)
          return(list('DSA_Score' = DSA, 'Number' = n, 'Gene_Names' = inter))
      }else{
        return(NA)
      }
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}
"
allReceptorDSA_Creator <- function(LR ,threshold = 100, RNAseqEXPCol, directed_mouse_ppi){
  allReceptorDSA <- data.frame(matrix(vector(), 0, 5, dimnames=list(c(), c("Receptor","DSA","n","nLVL1","nLVL2"))),
                               stringsAsFactors=F) %>% 
                        dplyr::mutate(Receptor = as.character(Receptor), 
                               DSA = as.numeric(DSA),
                               n = as.numeric(n), 
                               nLVL1 = as.numeric(nLVL1), 
                               nLVL2 = as.numeric(nLVL2))  
  
  #*#*# LVL 1 #*#*#
  for(Receptor in LR$Receptor %>% unique){
    LVL1_DSA     <- NA           
    LVL1_n       <- NA
    LVL2_Genes   <-NA
    GeneDSAstats <- DSAcalculator(Receptor, RNAseqEXPCol, threshold, directed_mouse_ppi)
    if(GeneDSAstats %>% is.list){
      LVL1_DSA <- GeneDSAstats$DSA_Score        
      LVL1_n   <- GeneDSAstats$Number
      LVL2_Genes <- GeneDSAstats$Gene_Names
      #*#*# LVL 2 #*#*#
      LVL2DSAvec <- NA
      LVL2nVec   <- NA
      for(Gene2 in LVL2_Genes){
        Gene2DSAstats <- DSAcalculator(Gene2, RNAseqEXPCol, threshold, directed_mouse_ppi)
        if(Gene2DSAstats %>% is.list){
          LVL2DSAvec <- append(LVL2DSAvec, Gene2DSAstats$DSA_Score)
          LVL2nVec   <- append(LVL2nVec, Gene2DSAstats$Number)
          LVL3_Genes <- Gene2DSAstats$Gene_Names
        }
      }
      #RECORD DATA if LVL3 exists
      finalDSA <- sum(LVL1_DSA*2 , mean(LVL2DSAvec, na.rm = T), na.rm = T) %>% round(2)
      allReceptorDSA <- add_row(allReceptorDSA, Receptor = Receptor, DSA = finalDSA,
                                n = sum(LVL1_n, LVL2nVec, na.rm = T), nLVL1 = LVL1_n, nLVL2 = sum(LVL2nVec, na.rm = T))
    }
  }
  row.names(allReceptorDSA) <- allReceptorDSA$Receptor
  return(allReceptorDSA)
}

DSAlegenedNcolor <- function(LR, threshold = 100, RNAseqEXPCol, directed_mouse_ppi){
  #allReceptorDSA <- allReceptorDSA_Creator(LR, 100, RNAseqEXPCol, directed_mouse_ppi)
  #allReceptorDSA[sapply(allReceptorDSA$DSA, is.infinite),"DSA"] <- 0
  #allReceptorDSA$DSA <- round(log2(allReceptorDSA$DSA),2)
  #View(allReceptorDSA)
  allReceptorDSA <- read.csv("C:\\Kfir_Thesis_Asaf_Madi\\Project_2 - BulkRNAseq\\Circos_Code\\DSA_Flow_Calculation_DF.csv",row.names = 1)
  if(dim(allReceptorDSA)[1] != 0){
    MAX <- max(allReceptorDSA$DSA, na.rm = T) 
    DSAcol_fun = colorRamp2( c(0, MAX), c("white","#800080"))
    lgd_DSA    = Legend(at = c(0, MAX),
                        col_fun = DSAcol_fun,
                        title_position = "leftcenter-rot",
                        title = "DSA",
                        legend_height = unit(2, "cm")
                        )
    print(paste0("Done DSA Calculation"))
    return(list(DSAcol_fun, lgd_DSA, allReceptorDSA))
  }else{return(NULL)}
} 

EXPlegenedNcolorFrom <- function(LR, z_scores = F){
  GeneExpressionFrom <- LR %>% dplyr::select(Ligand, LigandAbundance) %>% unique()
  GeneExpressionFrom$LigandAbundance <- log2(GeneExpressionFrom$LigandAbundance) %>% round(2)
  rownames(GeneExpressionFrom) <- GeneExpressionFrom$Ligand 
  GeneExpressionFrom$Z <- scale(GeneExpressionFrom$LigandAbundance) %>% round(2)
  
  if(z_scores == T){
      max_from <- max(GeneExpressionFrom$Z, na.rm = T)
      min_from <- min(GeneExpressionFrom$Z, na.rm = T)
      GeneExpressionFrom <- GeneExpressionFrom %>% dplyr::select(-c(LigandAbundance))
  }else{
      max_from <- max(GeneExpressionFrom$LigandAbundance, na.rm = T)
      min_from <- 0
  }
  
  EXPRcol_funFrom = colorRamp2( c(min_from, max_from), c("white","blue"))
  lgd_EXPRFrom    = Legend(at = c(min_from, max_from),
                           col_fun = EXPRcol_funFrom, title_position = "leftcenter-rot", title = "Ligands Z Val.",
                           legend_height = unit(2, "cm"))
  
  print("Done Ligand Expression")
  return(list(EXPRcol_funFrom, lgd_EXPRFrom, GeneExpressionFrom))
}

EXPlegenedNcolorTo <- function(LR, z_scores = F){
  GeneExpressionTo <- LR %>% dplyr::select(Receptor, ReceptorExp) %>% unique()
  GeneExpressionTo$ReceptorExp <- log2(GeneExpressionTo$ReceptorExp) %>% round(2)
  rownames(GeneExpressionTo) <- GeneExpressionTo$Receptor 
  GeneExpressionTo$Z <- scale(GeneExpressionTo$ReceptorExp) %>% round(2)
  
  if(z_scores == T){
    max_to <- max(GeneExpressionTo$Z, na.rm = T)
    min_to <- min(GeneExpressionTo$Z, na.rm = T)
    GeneExpressionTo <- GeneExpressionTo %>% dplyr::select(-c(ReceptorExp))
  }else{
    max_to <- max(GeneExpressionTo$ReceptorExp, na.rm = T)
    min_to <- 0
  }
  
  EXPRcol_funTo = colorRamp2( c(min_to, max_to), c("white","red"))
  lgd_EXPRTo    = Legend(at = c(min_to, max_to),
                         col_fun = EXPRcol_funTo, title_position = "leftcenter-rot", title = "Receptors Z Val.",
                         legend_height = unit(2, "cm"))
  
  print("Done Receptor Expression")
  return(list(EXPRcol_funTo, lgd_EXPRTo, GeneExpressionTo))
}
