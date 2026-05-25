
#================================================================================#
#### Análisis de RNA-seq de experimento de Olivo (expresion ~ luz + variedad) ####
#================================================================================#


#### FUNCIONES ####

# Definir función para generar diagramas de Venn (2 grupos)

VenniVidiVinci2 <- function(deseqOut1, deseqOut2) {
  vennList <- list(degs1 = rownames(deseqOut1[which(deseqOut1$padj < 0.05), ]), 
                   degs2 = rownames(deseqOut2[which(deseqOut2$padj < 0.05), ]))
  ggvenn(vennList, 
         fill_color = c("darkgrey", "lightblue"), 
         stroke_size = 0.5, 
         set_name_size = 4, 
         auto_scale = T, 
         show_percentage = F)
}

# Definir función para generar diagramas de Venn (3 grupos)

VenniVidiVinci3 <- function(deseqOut1, deseqOut2, deseqOut3) {
  vennList <- list(degs1 = rownames(deseqOut1[which(deseqOut1$padj < 0.05), ]), 
                   degs2 = rownames(deseqOut2[which(deseqOut2$padj < 0.05), ]),
                   degs3 = rownames(deseqOut3[which(deseqOut3$padj < 0.05), ]))
  ggvenn(vennList, 
         fill_color = c("lightgreen", "lightblue", "tomato"), 
         stroke_size = 0.5, 
         set_name_size = 4, 
         show_percentage = F)
}

# Definir función de calcular test hipergeométrico con dos salidas de DESeq2

autoPhyper <- function(deseqOut1, deseqOut2) {
  deglist1 <- rownames(deseqOut1[which(deseqOut1$padj <= 0.05),])
  group1 <- length(deglist1) # number of DEGs in group 1
  deglist2 <- rownames(deseqOut2[which(deseqOut2$padj <= 0.05),])
  group2 <- length(deglist2) # number of DEGs in group 2
  overlap <- length(intersect(deglist1, deglist2)) # number of overlapping DEGs between groups
  totalgenes <- length(unique(c(rownames(deseqOut1), rownames(deseqOut2)))) # all analyzed genes
  HGtest <- phyper(overlap, group2, totalgenes - group2, group1, lower.tail = FALSE, log.p = FALSE) # Hyper geometric test
  print(paste("Hypergeometric test p-value: ", HGtest)) # print the p-value
}

# Definir función del Volcano Plot

customVolcano <- function(deseqResults, dataTitle) {
  deseqResults$diffexpr <- "NO"
  deseqResults$diffexpr[deseqResults$log2FoldChange > 0 & deseqResults$padj < 0.05] <- "UP"
  deseqResults$diffexpr[deseqResults$log2FoldChange < 0 & deseqResults$padj < 0.05] <- "DOWN"
  ggplot(data = deseqResults, aes(x = log2FoldChange, y = -log10(padj), colour = diffexpr)) +	
    geom_point(alpha = 0.4, size = 1.75) + 
    theme_minimal() + labs(title = dataTitle) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "black", size = 1, linetype = "solid")) + 
    xlab("log2 fold change") + ylab("-log10 q-value") +
    scale_colour_manual(values = c("steelblue", "grey", "red")) + 
    scale_x_continuous(limits = c(min(deseqResults$log2FoldChange), max(deseqResults$log2FoldChange))) + 
    scale_y_continuous(limits = c(0, -log10(min(na.omit(deseqResults$padj)))))
}

# Función auxiliar para mapear columnas en topGO

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# Función para correr topGO en bucle (ACTUALIZAR NÚMERO DE NODOS POR EJECUCIÓN!!!)

topGOloop <- function(degList) {
  for (i in 1:length(degList)) {
    GOlist <- as.vector(degList[[i]]$padj)
    names(GOlist) <- rownames(degList[[i]])
    # Generate factor TRUE/FALSE for the genes in gene list present in whole gene set
    geneNames <- names(geneID2GO)
    myInterestingGenes <- names(GOlist)
    GOlist.list <- factor(as.integer(geneNames %in% myInterestingGenes))
    names(GOlist.list) <- geneNames
    ## BIOLOGICAL PROCESS ##
    # Create topGOdata object
    GOlist_BP <- new("topGOdata", 
                     description = "", 
                     ontology = "BP", # which GO term family to analyze
                     allGenes = GOlist.list, # input gene list
                     nodeSize = 5, # filter nodes (GO terms) with less than 5 seq counts
                     annot = annFUN.gene2GO, # function for extracting reference annotation from GO mapping file
                     gene2GO = geneID2GO) # mapping file
    # Classic Fisher exact test
    GOlist_BP.fet <- runTest(GOlist_BP, algorithm = "classic", statistic = "fisher")
    # Weighted Fisher exact test
    GOlist_BP.weight <- runTest(GOlist_BP, algorithm = "weight", statistic = "fisher")
    # Compile and save the results
    GOlist_BP.allRes <- GenTable(GOlist_BP, 
                                 classicFisher = GOlist_BP.fet, 
                                 weightFisher = GOlist_BP.weight,
                                 orderBy = "weightFisher", 
                                 ranksOf = "weightFisher", 
                                 numChar = 1000,
                                 topNodes = 2839) # BP: UniProt: 2839 nodes
    GOlist_BP.allRes$classicFisher[GOlist_BP.allRes$classicFisher == "< 1e-30"] <- "1e-30"
    GOlist_BP.classicFDR <- p.adjust(GOlist_BP.allRes$classicFisher, method = "fdr")
    GOlist_BP.allRes$classicFDR <- GOlist_BP.classicFDR
    GOlist_BP.allRes$weightFisher[GOlist_BP.allRes$weightFisher == "< 1e-30"] <- "1e-30"
    GOlist_BP.weightFDR <- p.adjust(GOlist_BP.allRes$weightFisher, method = "fdr")
    GOlist_BP.allRes$weightFDR <- GOlist_BP.weightFDR
    write.table(GOlist_BP.allRes, 
                file = paste0("topGO/", names(degList)[i], "_topGO_results_BP.txt"),
                sep = '\t', 
                quote = F)
    ## MOLECULAR FUNCTION ##
    # Create topGOdata object
    GOlist_MF <- new("topGOdata", 
                     description = "", 
                     ontology = "MF", # which GO term family to analyze
                     allGenes = GOlist.list, # input gene list
                     nodeSize = 5, # filter nodes (GO terms) with less than 5 seq counts
                     annot = annFUN.gene2GO, # function for extracting reference annotation from GO mapping file
                     gene2GO = geneID2GO) # mapping file
    # Classic Fisher exact test
    GOlist_MF.fet <- runTest(GOlist_MF, algorithm = "classic", statistic = "fisher")
    # Weighted Fisher exact test
    GOlist_MF.weight <- runTest(GOlist_MF, algorithm = "weight", statistic = "fisher")
    # Compile and save the results
    GOlist_MF.allRes <- GenTable(GOlist_MF, 
                                 classicFisher = GOlist_MF.fet, 
                                 weightFisher = GOlist_MF.weight,
                                 orderBy = "weightFisher", 
                                 ranksOf = "weightFisher", 
                                 numChar = 1000,
                                 topNodes = 895) # MF: UniProt: 895 nodes
    GOlist_MF.allRes$classicFisher[GOlist_MF.allRes$classicFisher == "< 1e-30"] <- "1e-30"
    GOlist_MF.classicFDR <- p.adjust(GOlist_MF.allRes$classicFisher, method = "fdr")
    GOlist_MF.allRes$classicFDR <- GOlist_MF.classicFDR
    GOlist_MF.allRes$weightFisher[GOlist_MF.allRes$weightFisher == "< 1e-30"] <- "1e-30"
    GOlist_MF.weightFDR <- p.adjust(GOlist_MF.allRes$weightFisher, method = "fdr")
    GOlist_MF.allRes$weightFDR <- GOlist_MF.weightFDR
    write.table(GOlist_MF.allRes, 
                file = paste0("topGO/", names(degList)[i], "_topGO_results_MF.txt"), 
                sep = '\t', 
                quote = F)
    ## CELLULAR COMPONENT ##
    # Create topGOdata object
    GOlist_CC <- new("topGOdata", 
                     description = "", 
                     ontology = "CC", # which GO term family to analyze
                     allGenes = GOlist.list, # input gene list
                     nodeSize = 5, # filter nodes (GO terms) with less than 5 seq counts
                     annot = annFUN.gene2GO, # function for extracting reference annotation from GO mapping file
                     gene2GO = geneID2GO) # mapping file
    # Classic Fisher exact test
    GOlist_CC.fet <- runTest(GOlist_CC, algorithm = "classic", statistic = "fisher")
    # Weighted Fisher exact test
    GOlist_CC.weight <- runTest(GOlist_CC, algorithm = "weight", statistic = "fisher")
    # Compile and save the results
    GOlist_CC.allRes <- GenTable(GOlist_CC, 
                                 classicFisher = GOlist_CC.fet, 
                                 weightFisher = GOlist_CC.weight,
                                 orderBy = "weightFisher", 
                                 ranksOf = "weightFisher", 
                                 numChar = 1000,
                                 topNodes = 465) # CC: UniProt: 465 nodes
    GOlist_CC.allRes$classicFisher[GOlist_CC.allRes$classicFisher == "< 1e-30"] <- "1e-30"
    GOlist_CC.classicFDR <- p.adjust(GOlist_CC.allRes$classicFisher, method = "fdr")
    GOlist_CC.allRes$classicFDR <- GOlist_CC.classicFDR
    GOlist_CC.allRes$weightFisher[GOlist_CC.allRes$weightFisher == "< 1e-30"] <- "1e-30"
    GOlist_CC.weightFDR <- p.adjust(GOlist_CC.allRes$weightFisher, method = "fdr")
    GOlist_CC.allRes$weightFDR <- GOlist_CC.weightFDR
    write.table(GOlist_CC.allRes, 
                file = paste0("topGO/", names(degList)[i], "_topGO_results_CC.txt"), 
                sep = '\t', 
                quote = F)
  } 
}

# Función para representar gráficamente subconjuntos de GO terms enriquecidos

goTermPlot <- function(goTermSubset) {
  mainPlot <- 
    ggplot(goTermSubset, 
           aes(x = fct_reorder(Term, weightFDR, .desc = TRUE), 
               y = -log10(weightFDR), 
               size = log2(Significant/Expected), 
               fill = -log10(weightFDR))) +
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5)) +
    scale_fill_continuous(low = 'steelblue', high = 'tomato') +
    xlab('') + ylab('Enrichment score') +
    labs(title = paste('Enriched GO terms')) +
    # Theme and general format
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      # Axis format
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      #Legend format
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    coord_flip()
  return(mainPlot)
}

