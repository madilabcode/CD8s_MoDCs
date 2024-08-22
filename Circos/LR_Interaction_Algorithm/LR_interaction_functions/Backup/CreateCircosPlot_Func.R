createCircosPlot <- function(LR = LR,
                             LigandReceptorTable,
                             toCol, 
                             threshold = 100,
                             RNAseqEXPCol = NA,
                             Figure_name){
  library(gridBase)
  library(circlize)
  library(ComplexHeatmap)
  
  LR_EXP.from <- EXPlegenedNcolorFrom(LR, z_scores = T)
  LR_EXP_fun.from <- LR_EXP.from[[1]]
  lgd_LR_EXP.from <- LR_EXP.from[[2]]
  EXP_Data_Sub.from <- LR_EXP.from[[3]]
  
  LR_EXP.to <- EXPlegenedNcolorTo(LR, z_scores = T)
  LR_EXP_fun.to <- LR_EXP.to[[1]]
  lgd_LR_EXP.to <- LR_EXP.to[[2]]
  EXP_Data_Sub.to <- LR_EXP.to[[3]]
  
  Activation_DSA <- DSAlegenedNcolor(LR, threshold, RNAseqEXPCol, directed_mouse_ppi)
  if(!is.null(Activation_DSA)){
    Activation_DSA_fun <- Activation_DSA[[1]]
    lgd_Activation_DSA <- Activation_DSA[[2]]
    Activation_DSA_Data_Sub <- Activation_DSA[[3]]  
  }
  
  FC <- FClegenedNcolor(LR, toCol)
  FC_fun <- FC[[1]]
  lgd_FC <- FC[[2]]
  FC_Data_Sub <- FC[[3]]
  
  # combine legends
  if(!is.null(Activation_DSA)){
    lgd_list_vertical = packLegend(lgd_LR_EXP.from, lgd_LR_EXP.to, lgd_FC, lgd_Activation_DSA, gap = unit(0.7, "cm")) 
  }else{
    lgd_list_vertical = packLegend(lgd_LR_EXP.from, lgd_LR_EXP.to, lgd_FC, gap = unit(0.7, "cm")) 
  }
  
  # receptors rectangles are grey
  grid.colPRE = cbind(as.character(LR$Receptor), c("grey")) %>% as.data.frame()
  # naming the sectors (slices) of the circle
  grid.col <- setNames(as.character(grid.colPRE$V2), grid.colPRE$V1)
  # vector of features that are in the plot
  factors <- unique(c(as.character(LR$Ligand),as.character(LR$Receptor)))
  # subset LR to just names, otherwise it affects the width of the arrows
  LR <- LR[,c("Ligand","Receptor")]
  
  #the plot itself
  pdf(Figure_name, width=10, height=7)
  circos.clear()
  #making the middle of the circle vertical as opposed to horizontal
  circosPlot <- circos.par(start.degree = 90, clock.wise = F, points.overflow.warning=FALSE)
  
  chordDiagram(LR,
               scale = F,
               big.gap = 20,
               grid.col = grid.col,
               directional = 1,
               annotationTrack = "grid",
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               link.sort = TRUE,
               preAllocateTracks = list(track.height = LR %>% dimnames() %>% unlist() %>% strwidth() %>% max())) 
  
                # Ligand Expression loop  
                for(Gene in EXP_Data_Sub.from %>% row.names()){circos.rect(
                  xleft = 0,
                  xright = c(grep(paste0("^",Gene,"$"), LR$Ligand)) %>% length(),
                  ybottom = 0.7, ytop = 0.9, 
                  col    = LR_EXP_fun.from(EXP_Data_Sub.from[Gene,2]), 
                  border = LR_EXP_fun.from(EXP_Data_Sub.from[Gene,2]),
                  sector.index = Gene, track.index = 1)}
                
                # Receptor Expression loop  
                for(Gene in EXP_Data_Sub.to %>% row.names()){circos.rect(
                  xleft = 0,
                  xright = c(grep(paste0("^",Gene,"$"), LR$Receptor)) %>% length(),
                  ybottom = 0.7, ytop = 0.9, 
                  col    = LR_EXP_fun.to(EXP_Data_Sub.to[Gene,2]), 
                  border = LR_EXP_fun.to(EXP_Data_Sub.to[Gene,2]),
                  sector.index = Gene, track.index = 1)}  
                
                # FC loop  
                for(receptor in intersect(FC_Data_Sub %>% rownames(), LR$Receptor)){circos.rect(
                  xleft = 0,
                  xright = c(grep(paste0("^",receptor,"$"),LR$Receptor)) %>% length()-0.1,
                  ybottom = 0.95, ytop = 1.15, 
                  col    = FC_fun(FC_Data_Sub[receptor,"log2FoldChange"]), 
                  border = FC_fun(FC_Data_Sub[receptor,"log2FoldChange"]),
                  sector.index = receptor, track.index = 1)}  
                
                if(!is.null(Activation_DSA)){
                #Activation DSA loop
                  for(receptor in intersect(Activation_DSA_Data_Sub$Receptor, LR$Receptor)){circos.rect(
                    xleft = 0,
                    xright = grep(paste0("^",receptor,"$"),LR$Receptor) %>% length-0.1,
                    ybottom = 1.2, ytop = 1.4+(Activation_DSA_Data_Sub[receptor,"n"]/50), 
                    col    = Activation_DSA_fun(Activation_DSA_Data_Sub[receptor,"DSA"]), 
                    border = Activation_DSA_fun(Activation_DSA_Data_Sub[receptor,"DSA"]),
                    sector.index = receptor, track.index = 1)}
                }
  
  circos.track(track.index = 1, panel.fun = function(x, y){
    circos.text(cex = 0.7, CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))}, bg.border = NA)
  
  circle_size = unit(1, "snpc")
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, height = circle_size, just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  upViewport()
  draw(lgd_list_vertical, just = c("left", "bottom"), x = unit(1, "cm"), y = unit(1, "cm"))
  while(!is.null(dev.list()))  dev.off()
  
  print(paste0("Done making Circos figure"))
}