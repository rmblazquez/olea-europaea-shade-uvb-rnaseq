
#================================================================================#
#### Análisis de RNA-seq de experimento de Olivo (expresion ~ luz + variedad) ####
#================================================================================#

# Establecemos entorno de trabajo

setwd("/home/rmblazquez/Documentos/Resultados/RNAseq_acebuche/")

library(ggplot2)
library(gridExtra)
library("RColorBrewer")
library(tidyverse)
library(dplyr)

library(dunn.test)

library(tximport)
library(DESeq2)

library(edgeR)
library(ggvenn)
library(EnhancedVolcano)
library(pheatmap)
library(reshape2)

source('/home/rmblazquez/Documentos/Scripts/RNAseq_olivo_functions.R')


#### GENERAR MATRIZ DE DISEÑO EXPERIMENTAL ####

# Importamos la matriz de diseño experimental, y elegimos solo los factores que nos interesan 
# En este caso nos quedamos con Treatment y con Variant
# ¡Ojo! He incluído el % GC y el % de duplicados promedios de las lecturas PE por muestra, desde MultiQC

design.raw <- read.csv("RNAseq_olivo_design.csv", sep = ';') # MOVER EL FICHERO AL WORKING DIRECTORY
design <- data.frame(sample = design.raw$sample, 
                     treatment = design.raw$treatment,
                     variant = design.raw$subspp,
                     time = design.raw$time)
design$treatment <- as.factor(design$treatment)
levels(design$treatment) <- c("control", "shade", "uvb")
design$variant <- as.factor(design$variant)
design$time <- as.factor(design$time)
levels(design$time) <- c("t0", "t1")
rownames(design) <- design$sample
design$sample <- NULL


#### GENERAR MATRIZ DE EXPRESIÓN ####

# Generar una variable "file" con referencia de muestra y path al "quant.sf"
# Para esto hay que reconocer los archivos quant en carpetas nombradas como 
# la muestra, y generar un vector con los paths de cada "quant.sf"

samples <- rownames(design)
dir <- paste0(getwd())
files <- c()
for (i in 1:length(samples)) {
  files[i] <- paste0(dir, "/Salmon/", samples[i] , ".quant/quant.sf")
}
names(files) <- samples

# Generar objeto tximport con los conteos

data.tx <- tximport(files, type = "salmon", txOut = TRUE) # Sumariza cuantificaciones de transcritos


#### ANÁLISIS DE EXPRESIÓN DIFERENCIAL ####

# datos solo para el tiempo 1

data.t1 <- data.tx$counts[, 25:49]

## Comparar efecto de tratamientos dentro de cada variante:

# Para comparación por grupos individuales, combinamos tratamiento y variante 
# en la variable group, y usamos la variable expDesign

group <- as.factor(paste0(design$treatment, design$variant))

expDesign <- model.matrix(~ 0 + group)

table <- data.frame(sample = rownames(expDesign), condition = group)

# colnames(expDesign) nos da los grupos del vector de contrastes:
# Para usar un grupo de referencia/control, ponemos -1 en su posición
# Para usar un grupo como tratamiento a comparar, ponemos 1 en su posición
# Para ignorar grupos, ponemos 0 en el vector
# Por ejemplo, para comparar los grupos controleuropaea y shadeeuropaea,
# usaremos el vector de contrastes c(-1,0,0,1,0,0,0,0,0)

# Analizar con DESeq2 comparaciones por pares de grupos

dds.int <- DESeqDataSetFromMatrix(countData = round(data.t1, digits = 0), 
                                  colData = table[25:49,], 
                                  design = expDesign[25:49,])

dds.int <- dds.int[rowSums(counts(dds.int)) > 10, ]

dds.int <- dds.int[apply(counts(dds.int), 1, sd) != 0, ]

dds.int <- DESeq(dds.int)

# Comparar control vs. shade en europaea: c(-1,0,0,1,0,0,0,0,0)

res.int.shade.eur <- results(dds.int, 
                             contrast = c(-1,0,0,1,0,0,0,0,0), 
                             alpha = 0.05) 

# Comparar control vs. uvb en europaea: c(-1,0,0,0,0,0,1,0,0) 

res.int.uvb.eur <- results(dds.int, 
                           contrast = c(-1,0,0,0,0,0,1,0,0), 
                           alpha = 0.05) 

# Comparar control vs. shade en guanchica: c(0,-1,0,0,1,0,0,0,0) 

res.int.shade.gua <- results(dds.int, 
                             contrast = c(0,-1,0,0,1,0,0,0,0), 
                             alpha = 0.05) 

# Comparar control vs. uvb en guanchica: c(0,-1,0,0,0,0,0,1,0)

res.int.uvb.gua <- results(dds.int, 
                           contrast = c(0,-1,0,0,0,0,0,1,0), 
                           alpha = 0.05) 

# Comparar control vs. shade en maroccana: c(0,0,-1,0,0,1,0,0,0) 

res.int.shade.mar <- results(dds.int, 
                             contrast = c(0,0,-1,0,0,1,0,0,0), 
                             alpha = 0.05) 

# Comparar control vs. uvb en maroccana: c(0,0,-1,0,0,0,0,0,1)

res.int.uvb.mar <- results(dds.int, 
                           contrast = c(0,0,-1,0,0,0,0,0,1), 
                           alpha = 0.05) 

## Comparación entre control y tratamientos a nivel general

# Comparar control vs. sombra: c(-1,-1,-1,1,1,1,0,0,0)

res.int.shade <- results(dds.int, 
                         contrast = c(-1,-1,-1,1,1,1,0,0,0), 
                         alpha = 0.05) 

# Comparar control vs. UVB: c(-1,-1,-1,0,0,0,1,1,1)

res.int.uvb <- results(dds.int, 
                       contrast = c(-1,-1,-1,0,0,0,1,1,1), 
                       alpha = 0.05) 

## Comparación entre variantes a nivel general

# Comparar europaea vs. guanchica: c(-1,1,0,-1,1,0-1,1,0)

res.int.eur_gua <- results(dds.int,
                           contrast = c(-1,1,0,-1,1,0,-1,1,0),
                           alpha = 0.05)

# Comparar europaea vs. maroccana: c(-1,0,1,-1,0,1,-1,0,1)

res.int.eur_mar <- results(dds.int,
                           contrast = c(-1,0,1,-1,0,1,-1,0,1),
                           alpha = 0.05)


# Comparar guanchica vs. maroccana: c(0,-1,1,0,-1,1,0,-1,1)

res.int.gua_mar <- results(dds.int,
                           contrast = c(0,-1,1,0,-1,1,0,-1,1),
                           alpha = 0.05)

## Comparación de trasfondo de variante con los datos del t0

# datos solo para el tiempo 0

data.t0 <- data.tx$counts[, 1:24]

dds.bkg <- DESeqDataSetFromMatrix(countData = round(data.t0, digits = 0),
                                  colData = table[1:24,],
                                  design = expDesign[1:24,])

dds.bkg <- dds.bkg[rowSums(counts(dds.bkg)) > 10, ]

dds.bkg <- dds.bkg[apply(counts(dds.bkg), 1, sd) != 0, ]

dds.bkg <- DESeq(dds.bkg)

# Comparar europaea vs. guanchica: c(-1,1,0,-1,1,0-1,1,0)

res.bkg.eur_gua <- results(dds.bkg,
                           contrast = c(-1,1,0,-1,1,0,-1,1,0),
                           alpha = 0.05)

# Comparar europaea vs. maroccana: c(-1,0,1,-1,0,1,-1,0,1)

res.bkg.eur_mar <- results(dds.bkg,
                           contrast = c(-1,0,1,-1,0,1,-1,0,1),
                           alpha = 0.05)

# Comparar guanchica vs. maroccana: c(0,-1,1,0,-1,1,0,-1,1)

res.bkg.gua_mar <- results(dds.bkg,
                           contrast = c(0,-1,1,0,-1,1,0,-1,1),
                           alpha = 0.05)

## Comparación con todos los tiempos (para WGCNA)

dds.all <- DESeqDataSetFromMatrix(countData = round(data.tx$counts, digits = 0), 
                                  colData = table, 
                                  design = expDesign)

dds.all <- dds.all[rowSums(counts(dds.all)) > 10, ]

dds.all <- dds.all[apply(counts(dds.all), 1, sd) != 0, ]

dds.all <- DESeq(dds.all)


#### EXPLORACIÓN DE RESULTADOS ####

## Comprobar el sentido de la sobre-expresión

# Verificamos que el FC toma valores positivos para sobreexpresión en tratamiento 
# y negativos para sobreexpresión en control:

# res.int.shade.eur[which(res.int.shade.eur$padj < 0.05), ][1,]

# Oleur061Scf0009g06001.1: logFC -12, más expresión en control (denominador) que en tratamiento (numerador)

# data.t1[which(rownames(data.t1) == "Oleur061Scf0009g06001.1"), ]

# Oleur061Scf0009g06001.1: expresión muestras controleuropaea (43, 45, 48)= 623.276, 2419.757, 1.645   
# Oleur061Scf0009g06001.1: expresión muestras shadeeuropaea (46) = 0

## Imprimir número de DEGs por comparación

# tratamiento dentro de variante
nrow(res.int.shade.eur[which(res.int.shade.eur$padj < 0.05), ])
nrow(res.int.uvb.eur[which(res.int.uvb.eur$padj < 0.05), ])
nrow(res.int.shade.gua[which(res.int.shade.gua$padj < 0.05), ])
nrow(res.int.uvb.gua[which(res.int.uvb.gua$padj < 0.05), ])
nrow(res.int.shade.mar[which(res.int.shade.mar$padj < 0.05), ])
nrow(res.int.uvb.mar[which(res.int.uvb.mar$padj < 0.05), ])

# solo tratamiento
nrow(res.int.shade[which(res.int.shade$padj < 0.05), ])
nrow(res.int.uvb[which(res.int.uvb$padj < 0.05), ])

# solo variante (t1)
nrow(res.int.eur_gua[which(res.int.eur_gua$padj < 0.05), ])
nrow(res.int.eur_mar[which(res.int.eur_mar$padj < 0.05), ])
nrow(res.int.gua_mar[which(res.int.gua_mar$padj < 0.05), ])

# solo variante (t0)
nrow(res.bkg.eur_gua[which(res.bkg.eur_gua$padj < 0.05), ])
nrow(res.bkg.eur_mar[which(res.bkg.eur_mar$padj < 0.05), ])
nrow(res.bkg.gua_mar[which(res.bkg.gua_mar$padj < 0.05), ])

## MDS

# datos t1

colorTreatment <- design$treatment[25:49]
levels(colorTreatment) <- c('green', 'blue', 'red')

shapeVariant <- design$variant[25:49]
levels(shapeVariant) <- c(23, 24, 25)

plotMDS(dds.int, gene.selection = "common", 
        top = dim(counts(dds.int))[2], 
        dim.plot = c(1,2), 
        # labels = colnames(data.t1), 
        labels = NULL,
        pch = as.integer(shapeVariant),
        col = as.vector(colorTreatment))

# datos t0

colorTreatment2 <- design$treatment[1:24]
levels(colorTreatment2) <- c('seagreen', 'purple', 'orange')

shapeVariant2 <- design$variant[1:24]
levels(shapeVariant2) <- c(23, 24, 25)

plotMDS(dds.bkg, gene.selection = "common", 
        top = dim(counts(dds.bkg))[2], 
        dim.plot = c(1,2), 
        # labels = colnames(data.t0), 
        labels = NULL,
        pch = as.integer(shapeVariant2),
        col = as.vector(colorTreatment2))

# datos t0 y t1

colorTreatmentTime <- as.factor(paste0(design$treatment, design$time))
levels(colorTreatmentTime) <- c('seagreen', 'green', 'purple','blue', 'orange', 'red')

shapeVariantAll <- design$variant
levels(shapeVariantAll) <- c(23, 24, 25)

plotMDS(dds.all, gene.selection = "common", 
        top = dim(counts(dds.all))[2], 
        dim.plot = c(1,2), 
        # labels = design.raw$genotype,
        # labels = colnames(data.tx$counts),
        labels = NULL,
        pch = as.integer(shapeVariantAll),
        col = as.vector(colorTreatmentTime))


## Boxplot para normalización

par(mfcol = c(1,2))
boxplot(log(counts(dds.int),2), col = colorTreatment)
boxplot(log(counts(dds.int, normalize = TRUE), 2), col = colorTreatment)

## Venn Diagram

# bivariante

VenniVidiVinci2(res.int.shade, res.int.uvb)

VenniVidiVinci3(res.int.eur_gua, res.int.eur_mar, res.int.gua_mar)

VenniVidiVinci3(res.bkg.eur_gua, res.bkg.eur_mar, res.bkg.gua_mar)

# dentro de grupo

VenniVidiVinci3(res.int.shade.eur, res.int.shade.gua, res.int.shade.mar)

VenniVidiVinci3(res.int.uvb.eur, res.int.uvb.gua, res.int.uvb.mar)

# t0 vs t1 por variantes

VenniVidiVinci2(res.int.eur_gua, res.bkg.eur_gua)

VenniVidiVinci2(res.int.eur_mar, res.bkg.eur_mar)

VenniVidiVinci2(res.int.gua_mar, res.bkg.gua_mar)

# DEGs por variante t0

eur.t0 <-
unique(
  c(
    rownames(res.bkg.eur_gua[which(res.bkg.eur_gua$padj < 0.05 & res.bkg.eur_gua$log2FoldChange > 0), ]),
    rownames(res.bkg.eur_mar[which(res.bkg.eur_mar$padj < 0.05 & res.bkg.eur_gua$log2FoldChange > 0), ])
  )
)

gua.t0 <-
unique(
  c(
    rownames(res.bkg.gua_mar[which(res.bkg.gua_mar$padj < 0.05 & res.bkg.gua_mar$log2FoldChange > 0), ]),
    rownames(res.bkg.eur_gua[which(res.bkg.eur_gua$padj < 0.05 & res.bkg.eur_gua$log2FoldChange < 0), ])
  )
)

mar.t0 <-
unique(
  c(
    rownames(res.bkg.gua_mar[which(res.bkg.gua_mar$padj < 0.05 & res.bkg.gua_mar$log2FoldChange < 0), ]),
    rownames(res.bkg.eur_mar[which(res.bkg.eur_mar$padj < 0.05 & res.bkg.eur_mar$log2FoldChange < 0), ])
  )
)

ggvenn(list(eur.t0 = eur.t0, gua.t0 = gua.t0, mar.t0 = mar.t0), 
       fill_color = c('seagreen', 'purple', 'orange'), 
       stroke_size = 0.5, 
       set_name_size = 4, 
       show_percentage = F)

# DEGs por variante t1

eur.t1 <-
  unique(
    c(
      rownames(res.int.eur_gua[which(res.int.eur_gua$padj < 0.05 & res.int.eur_gua$log2FoldChange > 0), ]),
      rownames(res.int.eur_mar[which(res.int.eur_mar$padj < 0.05 & res.int.eur_gua$log2FoldChange > 0), ])
    )
  )

gua.t1 <-
  unique(
    c(
      rownames(res.int.gua_mar[which(res.int.gua_mar$padj < 0.05 & res.int.gua_mar$log2FoldChange > 0), ]),
      rownames(res.int.eur_gua[which(res.int.eur_gua$padj < 0.05 & res.int.eur_gua$log2FoldChange < 0), ])
    )
  )

mar.t1 <-
  unique(
    c(
      rownames(res.int.gua_mar[which(res.int.gua_mar$padj < 0.05 & res.int.gua_mar$log2FoldChange < 0), ]),
      rownames(res.int.eur_mar[which(res.int.eur_mar$padj < 0.05 & res.int.eur_mar$log2FoldChange < 0), ])
    )
  )

ggvenn(list(eur.t1 = eur.t1, gua.t1 = gua.t1, mar.t1 = mar.t1), 
       fill_color = c("lightgreen", "lightblue", "tomato"), 
       stroke_size = 0.5, 
       set_name_size = 4, 
       show_percentage = F)

# DEGs t0 vs DEGs t1 por variante

ggvenn(list(eur.t0 = eur.t0, eur.t1 = eur.t1), 
       fill_color = c("lightblue", "grey"), 
       stroke_size = 0.5, 
       set_name_size = 4, 
       show_percentage = F)

ggvenn(list(gua.t0 = gua.t0, gua.t1 = gua.t1), 
       fill_color = c("lightblue", "grey"), 
       stroke_size = 0.5, 
       set_name_size = 4, 
       show_percentage = F)

ggvenn(list(mar.t0 = mar.t0, mar.t1 = mar.t1), 
       fill_color = c("lightblue", "grey"), 
       stroke_size = 0.5, 
       set_name_size = 4, 
       show_percentage = F)


## Hypergeometric test

# bivariante

autoPhyper(res.int.shade, res.int.uvb)

autoPhyper(res.int.eur_gua, res.int.eur_mar)
autoPhyper(res.int.eur_gua, res.int.gua_mar)
autoPhyper(res.int.eur_mar, res.int.gua_mar)

# dentro de grupo

autoPhyper(res.int.shade.eur, res.int.shade.gua)
autoPhyper(res.int.shade.eur, res.int.shade.mar)
autoPhyper(res.int.shade.gua, res.int.shade.mar)

autoPhyper(res.int.uvb.eur, res.int.uvb.gua)
autoPhyper(res.int.uvb.eur, res.int.uvb.mar)
autoPhyper(res.int.uvb.gua, res.int.uvb.mar)

## Volcano plots

# bivariantes

grid.arrange(customVolcano(res.int.shade, "Control vs. Shade"),
             customVolcano(res.int.uvb, "Control vs. UVB"),
             ncol = 2)

grid.arrange(customVolcano(res.int.eur_gua, "europaea vs. guanchica"),
             customVolcano(res.int.eur_mar, "europaea vs. maroccana"),
             customVolcano(res.int.gua_mar, "guanchica vs. maroccana"),
             ncol = 3)

grid.arrange(customVolcano(res.bkg.eur_gua, "europaea vs. guanchica"),
             customVolcano(res.bkg.eur_mar, "europaea vs. maroccana"),
             customVolcano(res.bkg.gua_mar, "guanchica vs. maroccana"),
             ncol = 3)

# dentro de grupo

grid.arrange(customVolcano(res.int.shade.eur, "control vs. shade (europaea)"),
             customVolcano(res.int.shade.gua, "control vs. shade (guanchica)"),
             customVolcano(res.int.shade.mar, "control vs. shade (maroccana)"),
             customVolcano(res.int.uvb.eur, "control vs. UVB (europaea)"),
             customVolcano(res.int.uvb.gua, "control vs. UVB (guanchica)"),
             customVolcano(res.int.uvb.mar, "control vs. UVB (maroccana)"),
             ncol = 3, nrow = 2)


## Heatmap para comparación de sombra y ultravioleta

# Seleccionar DEGs de treatment

genes_sig.int.shade <- res.int.shade[which(res.int.shade$padj <= 0.05),]

genes_sig.int.uvb <- res.int.uvb[which(res.int.uvb$padj <= 0.05),]

# Extraer conteos normalizados del objeto DESeq2

mat_norm <- as.data.frame(counts(dds.int, normalized = TRUE))

# Subconjunto con genes significativos

genes_sig <- rbind(genes_sig.int.shade, genes_sig.int.uvb)

mat_heatmap <- log2(mat_norm[rownames(mat_norm) %in% unique(rownames(genes_sig)), ] + 1)

# Generar vector de grupos para tratamiento

treatmentGroup <- data.frame(Treatment = as.factor(design$treatment[25:49]),
                             Variant = as.factor(design$variant[25:49])
                             )

rownames(treatmentGroup) <- colnames(mat_heatmap)

treatmentSort <- c(
  25,28,31,34,37,40,43,45,48,# Control
  26,29,32,35,38,41,46,      # Shade
  27,30,33,36,39,42,44,47,49 # UVB
  ) # comentar si no queremos ordenar por tratamiento las muestras

treatmentSort <- treatmentSort - 24

treatmentGroup <- treatmentGroup[treatmentSort, ]

mat_heatmap <- mat_heatmap[, treatmentSort]

annotation_colors <- list(
  Treatment = c(
    control = "lightblue",
    shade = "black",
    uvb = "red"),
  Variant = c(
    europaea = "lightgreen",
    guanchica = "orange",
    maroccana = "blue")
)

# Heatmap

pheatmap(
  mat_heatmap,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  annotation_col = treatmentGroup,
  annotation_colors = annotation_colors,
  fontsize_row = 8,
  color = colorRampPalette(c("white", "blue", "darkblue"))(100),
  border_color = NA
)

## Expresión de genes individuales

# Listas de DEGs por sombra en común a todas las variedades (tratamiento por variedad):

DEGs.int.shade <-
  intersect(
    intersect(
      rownames(res.int.shade.eur[which(res.int.shade.eur$padj < 0.05), ]), 
      rownames(res.int.shade.gua[which(res.int.shade.gua$padj < 0.05), ])
    ), 
    rownames(res.int.shade.mar[which(res.int.shade.mar$padj < 0.05), ])
  )

# Listas de DEGs por UVB en común a todas las variedades (tratamiento por variedad):

DEGs.int.uvb <- 
  intersect(
    intersect(
      rownames(res.int.uvb.eur[which(res.int.uvb.eur$padj < 0.05), ]), 
      rownames(res.int.uvb.gua[which(res.int.uvb.gua$padj < 0.05), ])
    ), 
    rownames(res.int.uvb.mar[which(res.int.uvb.mar$padj < 0.05), ])
  )

# Barplot para comprobar niveles de expresión por DEG (convertir en función?)

DEGs.int.shade.counts <- 
  data.t1[which(rownames(data.t1) %in% DEGs.int.shade),]
colnames(DEGs.int.shade.counts) <- paste0(design$treatment[25:49], design$variant[25:49])

plotList.shade <- list()
for (i in 1:length(DEGs.int.shade)) {
  geneExprPlot.shade <- melt(t(DEGs.int.shade.counts[i,]))
  geneExprPlot.shade$Var1 <- NULL
  colnames(geneExprPlot.shade) <- c("treatment", "counts")
  geneExprPlot.shade$treatment <- factor(geneExprPlot.shade$treatment, levels = 
                                   c("controleuropaea", "controlguanchica", "controlmaroccana",
                                     "shadeeuropaea", "shadeguanchica", "shademaroccana", 
                                     "uvbeuropaea", "uvbguanchica", "uvbmaroccana"))
  plot <- 
    ggplot(data = geneExprPlot.shade, aes(x = treatment, y = counts, fill = treatment)) + 
      geom_bar(stat = "identity") +
      ggtitle(rownames(DEGs.int.shade.counts)[i])
  plotList.shade[[i]] <- plot
}

DEGs.int.uvb.counts <- 
  data.t1[which(rownames(data.t1) %in% DEGs.int.uvb),]
colnames(DEGs.int.uvb.counts) <- paste0(design$treatment[25:49], design$variant[25:49])

plotList.uvb <- list()
for (i in 1:length(DEGs.int.uvb)) {
  geneExprPlot.uvb <- melt(t(DEGs.int.uvb.counts[i,]))
  geneExprPlot.uvb$Var1 <- NULL
  colnames(geneExprPlot.uvb) <- c("treatment", "counts")
  geneExprPlot.uvb$treatment <- factor(geneExprPlot.uvb$treatment, levels = 
                                         c("controleuropaea", "controlguanchica", "controlmaroccana",
                                           "shadeeuropaea", "shadeguanchica", "shademaroccana", 
                                           "uvbeuropaea", "uvbguanchica", "uvbmaroccana"))
  plot <- 
    ggplot(data = geneExprPlot.uvb, aes(x = treatment, y = counts, fill = treatment)) + 
      geom_bar(stat = "identity") +
      ggtitle(rownames(DEGs.int.uvb.counts)[i])
  plotList.uvb[[i]] <- plot
}

