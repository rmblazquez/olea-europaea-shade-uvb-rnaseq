
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

library(reshape2)

library(topGO)
library(scales)
library(forcats)
source('/home/rmblazquez/Documentos/Scripts/RNAseq_olivo_functions.R')


#### ANÁLISIS GO TERMS #### 

## Código para generar fichero de GO terms a partir de UniProt.
# UNIPROT MUESTRA IDS QUE NO SON DE PLANTAS Y METEN RUIDO, DEPURAR!!!!

# Del excel Olea_europaea_v061_annotations.xlsx, exportar como CSV delimitado por comas.
# En linux, extraer UniProt IDs de con `cut -f4 -d$',' Olea_europaea_v061_annotations.csv | cut -d$'.' -f1 | sort | uniq | tr '\n' ' ' > UniProt_IDs_list.txt`
# Importamos los códigos UniProt del fichero de anotación del genoma del olivo, exportade en CSV del excel

oleurIDs <- read.csv("topGO/Annotations/Olea_europaea_v061_annotations.csv", sep = ",", header = TRUE)

# Los IDs de UniProt se han analizado con UniProtKB ID mapper, exportando la tabla con los GO term IDs

uniprotIDs <- read.csv("topGO/Annotations/Oleur_v061_UniProt_idmapping_2026_01_27.tsv", sep = "\t", header = TRUE)

# Generamos un data frame con los gene IDs de olivo y los uniprot IDs asociados, los ordenamos por uniprotID

oleurIDs.uniprot <- data.frame(uniprotID = oleurIDs$swp_term, oleurID = oleurIDs$Gene.ID)

oleurIDs.uniprot.t <- oleurIDs.uniprot[which(oleurIDs.uniprot$uniprotID != ""), ]

oleurIDs.uniprot.sort <- oleurIDs.uniprot.t[order(oleurIDs.uniprot.t$uniprotID), ]

# Generamos otro data frame con el ID uniprot y el ID GOterm, y se ordena alfabéticamente

uniprotIDs.GO <- data.frame(uniprotID = uniprotIDs$From, GO.terms = uniprotIDs$Gene.Ontology.IDs)

uniprotIDs.GO.sort <- uniprotIDs.GO[order(uniprotIDs.GO$uniprotID), ]

# Generamos una variable con el número de veces que se repiten los códigos UniProt en el fichero de anotación

uniprot.times <- table(oleurIDs.uniprot$uniprotID[which(oleurIDs.uniprot$uniprotID != "")]) # Filtrar nombres en blanco

# Se importan a R y se generan las repeticiones de GO terms con este código:

names(uniprot.times) <- NULL
uniprot.times <- as.vector(uniprot.times)

id.times.vec <- c()
up.times.vec <- c()
go.times.vec <- c()
for (i in 1:length(uniprot.times)) {
  id.times <- rep(oleurIDs.uniprot.sort[i,2], each = uniprot.times[i])
  up.times <- rep(uniprotIDs.GO.sort[i,1], each = uniprot.times[i])
  go.times <- rep(uniprotIDs.GO.sort[i,2], each = uniprot.times[i])
  id.times.vec <- c(id.times.vec, id.times)
  up.times.vec <- c(up.times.vec, up.times)
  go.times.vec <- c(go.times.vec, go.times)
}

write.table(unique(data.frame(GeneID = oleurIDs.uniprot.sort$oleurID,
			      UnirpotID.oleur = oleurIDs.uniprot.sort$uniprotID,
			      UniprotID.uniprot = up.times.vec,
			      GOterm = gsub(";", ",", go.times.vec))), 
            "topGO/Annotations/uniprot.GOreps.txt", 
            sep = '\t', quote = FALSE, row.names = FALSE)  




# PROBLEMAS A LA HORA DE GENERAR EL GO TERM UNIVERSE DE TAIR!!!! NO USAR HASTA ARREGLARLO!!!
# Se puede hacer lo mismo con los IDs de TAIR: 
# En linux, extraer TAIR IDs de con `cut -f2 -d$',' Olea_europaea_v061_annotations.csv | cut -d$'.' -f1 | sort | uniq | tr '\n' ' ' > TAIR_IDs_list.txt`
# Generamos un data frame con los gene IDs de olivo y los TAIR IDs asociados, los ordenamos por TAIR ID
#oleurIDs.araport <- data.frame(araportID = substr(oleurIDs$Atha_gene, 1, 9), oleurID = oleurIDs$Gene.ID)
#oleurIDs.araport.t <- oleurIDs.araport[which(oleurIDs.araport$araportID != ""), ]
#oleurIDs.araport.sort <- oleurIDs.araport.t[order(oleurIDs.araport.t$araportID), ]
# Importar a ID mapper de UniProt, seleccionar Araport -> Uniprot, añadir columna de GO term IDs, y exportar GO terms.
#araportIDs <- read.csv("topGO/Annotations/Oleur_v061_TAIR_idmapping_2026_01_27.tsv", sep = "\t", header = TRUE)
# Generamos otro data frame con el TAIR ID y el ID GOterm, y se ordena alfabéticamente
#araportIDs.GO <- data.frame(araportID = araportIDs$From, GO.terms = araportIDs$Gene.Ontology.IDs)
#araportIDs.GO.sort <- araportIDs.GO[order(araportIDs.GO$araportID), ]
# Generamos una variable con el número de veces que se repiten los códigos TAIR en el fichero de anotación
#araport.times <- table(oleurIDs.araport$araportID[which(oleurIDs.araport$araportID != "")]) # Filtrar nombres en blanco
# Se importan a R y se generan las repeticiones de GO terms con este código:
#names(araport.times) <- NULL
#araport.times <- as.vector(araport.times)
#id.times.vec <- c()
#at.times.vec <- c()
#go.times.vec <- c()
#for (i in 1:length(araport.times)) {
#  id.times <- rep(oleurIDs.araport.sort[i,2], each = araport.times[i])
#  at.times <- rep(araportIDs.GO.sort[i,1], each = araport.times[i])
#  go.times <- rep(araportIDs.GO.sort[i,2], each = araport.times[i])
#  id.times.vec <- c(id.times.vec, id.times)
#  at.times.vec <- c(at.times.vec, at.times)
#  go.times.vec <- c(go.times.vec, go.times)
#}
#write.table(unique(data.frame(GeneID = oleurIDs.araport.sort$oleurID, 
#                              araportID.oleur = oleurIDs.araport.sort$araportID,
#                              araportID.uniprot = at.times.vec,
#			      GOterm = gsub(";", ",", go.times.vec))),
#            "topGO/Annotations/TAIR.GOreps.txt", 
#            sep = '\t', quote = FALSE, row.names = FALSE) 


## Importar GO terms del genoma del olivo (opciones: UniProt o TAIR)

GOterms <- read.csv("topGO/Annotations/uniprot.GOreps.txt", sep = '\t', header = TRUE) # Comprobar que las dos columnas de UniProtIDs coinciden

GOterms.topGO <- data.frame(GeneID = GOterms$GeneID, GOterm = GOterms$GOterm) # Seleccionar columnas con los IDs de olivo y los GO terms

# Seleccionar solo los GO terms de genes expresados en nuestro transcriptoma

GOterms.int <- GOterms.topGO[which(GOterms.topGO$GeneID %in% rownames(counts(dds.int))), ]

# Exportar universo GO como tabla.

# write.table(GOterms.int, "topGO/Annotations/GOterm_universe_UniProt.csv", sep = '\t', quote = F, row.names = F)

# Importar con topGO::readMappings()

geneID2GO <- readMappings(file = "/home/rmblazquez/Documentos/Resultados/RNAseq_acebuche/topGO/Annotations/GOterm_universe_UniProt.csv")

# Seleccionar DEGs de treatment y variant (GENERADOS EN SCRIPT RNAseq_olivo_DGEA.R)

genes_sig.int.shade <- res.int.shade[which(res.int.shade$padj <= 0.05),]

genes_sig.int.uvb <- res.int.uvb[which(res.int.uvb$padj <= 0.05),]

genes_sig.int.eur_gua <- res.int.eur_gua[which(res.int.eur_gua$padj <= 0.05),]

genes_sig.int.eur_mar <- res.int.eur_mar[which(res.int.eur_mar$padj <= 0.05),]

genes_sig.int.gua_mar <- res.int.gua_mar[which(res.int.gua_mar$padj <= 0.05),]

genes_sig.bkg.eur_gua <- res.bkg.eur_gua[which(res.bkg.eur_gua$padj <= 0.05),]

genes_sig.bkg.eur_mar <- res.bkg.eur_mar[which(res.bkg.eur_mar$padj <= 0.05),]

genes_sig.bkg.gua_mar <- res.bkg.gua_mar[which(res.bkg.gua_mar$padj <= 0.05),]

# Seleccionar DEGs de treatment.variant

genes_sig.int.shade.eur <- res.int.shade.eur[which(res.int.shade.eur$padj <= 0.05),]

genes_sig.int.shade.gua <- res.int.shade.gua[which(res.int.shade.gua$padj <= 0.05),]

genes_sig.int.shade.mar <- res.int.shade.mar[which(res.int.shade.mar$padj <= 0.05),]

genes_sig.int.uvb.eur <- res.int.uvb.eur[which(res.int.uvb.eur$padj <= 0.05),]

genes_sig.int.uvb.gua <- res.int.uvb.gua[which(res.int.uvb.gua$padj <= 0.05),]

genes_sig.int.uvb.mar <- res.int.uvb.mar[which(res.int.uvb.mar$padj <= 0.05),]

# Generar lista de sets para analizar (cada lista es un objeto DESeq2 con los DEGs)

degList.treatment <- list(shade.up = genes_sig.int.shade[which(genes_sig.int.shade$log2FoldChange > 0),], 
                          shade.down = genes_sig.int.shade[which(genes_sig.int.shade$log2FoldChange < 0),], 
                          uvb.up = genes_sig.int.uvb[which(genes_sig.int.uvb$log2FoldChange > 0),],
                          uvb.down = genes_sig.int.uvb[which(genes_sig.int.uvb$log2FoldChange < 0),])

degList.variants <- list(eur_gua.eur = genes_sig.int.eur_gua[which(genes_sig.int.eur_gua$log2FoldChange > 0),], 
                         eur_gua.gua = genes_sig.int.eur_gua[which(genes_sig.int.eur_gua$log2FoldChange < 0),], 
                         eur_mar.eur = genes_sig.int.eur_mar[which(genes_sig.int.eur_mar$log2FoldChange > 0),],
                         eur_mar.mar = genes_sig.int.eur_mar[which(genes_sig.int.eur_mar$log2FoldChange < 0),],
                         gua_mar.gua = genes_sig.int.gua_mar[which(genes_sig.int.gua_mar$log2FoldChange > 0),],
                         gua_mar.mar = genes_sig.int.gua_mar[which(genes_sig.int.gua_mar$log2FoldChange < 0),])

degList.bkgrVar <- list(bkg.eur_gua.eur = genes_sig.bkg.eur_gua[which(genes_sig.bkg.eur_gua$log2FoldChange > 0),], 
                        bkg.eur_gua.gua = genes_sig.bkg.eur_gua[which(genes_sig.bkg.eur_gua$log2FoldChange < 0),], 
                        bkg.eur_mar.eur = genes_sig.bkg.eur_mar[which(genes_sig.bkg.eur_mar$log2FoldChange > 0),],
                        bkg.eur_mar.mar = genes_sig.bkg.eur_mar[which(genes_sig.bkg.eur_mar$log2FoldChange < 0),],
                        bkg.gua_mar.gua = genes_sig.bkg.gua_mar[which(genes_sig.bkg.gua_mar$log2FoldChange > 0),],
                        bkg.gua_mar.mar = genes_sig.bkg.gua_mar[which(genes_sig.bkg.gua_mar$log2FoldChange < 0),])

degList.withinVar <- list(shade.eur.up = genes_sig.int.shade.eur[which(genes_sig.int.shade.eur$log2FoldChange > 0),],
                          shade.eur.down = genes_sig.int.shade.eur[which(genes_sig.int.shade.eur$log2FoldChange < 0),],
                          shade.gua.up = genes_sig.int.shade.gua[which(genes_sig.int.shade.gua$log2FoldChange > 0),],
                          shade.gua.down = genes_sig.int.shade.gua[which(genes_sig.int.shade.gua$log2FoldChange < 0),],
                          shade.mar.up = genes_sig.int.shade.mar[which(genes_sig.int.shade.mar$log2FoldChange > 0),],
                          shade.mar.down = genes_sig.int.shade.mar[which(genes_sig.int.shade.mar$log2FoldChange < 0),],
                          uvb.eur.up = genes_sig.int.uvb.eur[which(genes_sig.int.uvb.eur$log2FoldChange > 0),],
                          uvb.eur.down = genes_sig.int.uvb.eur[which(genes_sig.int.uvb.eur$log2FoldChange < 0),],
                          uvb.gua.up = genes_sig.int.uvb.gua[which(genes_sig.int.uvb.gua$log2FoldChange > 0),],
                          uvb.gua.down = genes_sig.int.uvb.gua[which(genes_sig.int.uvb.gua$log2FoldChange < 0),],
                          uvb.mar.up = genes_sig.int.uvb.mar[which(genes_sig.int.uvb.mar$log2FoldChange > 0),],
                          uvb.mar.down = genes_sig.int.uvb.mar[which(genes_sig.int.uvb.mar$log2FoldChange < 0),])

topGOloop(degList.treatment)

topGOloop(degList.variants)

topGOloop(degList.bkgrVar)

topGOloop(degList.withinVar)

# Set de DEGs solapados entre sombra y UVB

res.int.overlap.names <-
  intersect(
    rownames(res.int.shade[which(res.int.shade$padj < 0.05), ]), 
    rownames(res.int.uvb[which(res.int.uvb$padj < 0.05), ]))

res.int.overlap.shade <- res.int.shade[rownames(res.int.shade) %in% res.int.overlap.names,]

res.int.overlap.uvb <- res.int.uvb[rownames(res.int.uvb) %in% res.int.overlap.names,]

res.int.overlap.shade2 <- res.int.overlap.shade$padj
names(res.int.overlap.shade2) <- rownames(res.int.overlap.shade)

res.int.overlap.uvb2 <- res.int.overlap.uvb$padj
names(res.int.overlap.uvb2) <- rownames(res.int.overlap.uvb)

res.int.overlap.df <- data.frame(shade.padj = res.int.overlap.shade2, uvb.padj = res.int.overlap.uvb2)
res.int.overlap.df$padj <- apply(res.int.overlap.df, 1, mean)

degList.int.overlap <- list(treatment.overlap = res.int.overlap.df)

topGOloop(degList.int.overlap)


# Hay que crear subcarpetas y mover manualmente los ficheros de topGO para seguir el análisis

# Representar gráficamente los 5 términos GO (BP) enriquecidos de cada comparación:
# Nos centramos en la ontología nde Biological Process, que es la más relevante

# GO terms (BP) de efecto general de tratamiento

int.shade.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/treatment/shade.up_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

int.shade.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/treatment/shade.down_topGO_results_BP.txt",
                                         sep = '\t',
                                         header = TRUE)

int.uvb.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/treatment/uvb.up_topGO_results_BP.txt",
                                     sep = '\t',
                                     header = TRUE)

int.uvb.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/treatment/uvb.down_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

goTermPlot(int.shade.up.enrichedGOs[1:5,])
goTermPlot(int.shade.down.enrichedGOs[1:5,])
goTermPlot(int.uvb.up.enrichedGOs[1:5,])
goTermPlot(int.uvb.down.enrichedGOs[1:5,])

int.overlap.enrichedGOs <- read.delim("topGO/UniProt_GOs/treatment/treatment.overlap_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

goTermPlot(two.overlap.enrichedGOs[1:5,])


# GO terms (BP) de efecto general de variante

int.eur_gua.eur.enrichedGOs <- read.delim("topGO/UniProt_GOs/variant/eur_gua.eur_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

int.eur_gua.gua.enrichedGOs <- read.delim("topGO/UniProt_GOs/variant/eur_gua.gua_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

int.eur_mar.eur.enrichedGOs <- read.delim("topGO/UniProt_GOs/variant/eur_mar.eur_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

int.eur_mar.mar.enrichedGOs <- read.delim("topGO/UniProt_GOs/variant/eur_mar.mar_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

int.gua_mar.gua.enrichedGOs <- read.delim("topGO/UniProt_GOs/variant/gua_mar.gua_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

int.gua_mar.mar.enrichedGOs <- read.delim("topGO/UniProt_GOs/variant/gua_mar.mar_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

goTermPlot(two.eur_gua.eur.enrichedGOs[1:5,])
goTermPlot(two.eur_gua.gua.enrichedGOs[1:5,])
goTermPlot(two.eur_mar.eur.enrichedGOs[1:5,])
goTermPlot(two.eur_mar.mar.enrichedGOs[1:5,])
goTermPlot(two.gua_mar.gua.enrichedGOs[1:5,])
goTermPlot(two.gua_mar.mar.enrichedGOs[1:5,])

# GO terms (BP) de efecto de fondo de variante (t0)

bkg.eur_gua.eur.enrichedGOs <- read.delim("topGO/UniProt_GOs/bkgr/bkg.eur_gua.eur_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

bkg.eur_gua.gua.enrichedGOs <- read.delim("topGO/UniProt_GOs/bkgr/bkg.eur_gua.gua_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

bkg.eur_mar.eur.enrichedGOs <- read.delim("topGO/UniProt_GOs/bkgr/bkg.eur_mar.eur_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

bkg.eur_mar.mar.enrichedGOs <- read.delim("topGO/UniProt_GOs/bkgr/bkg.eur_mar.mar_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

bkg.gua_mar.gua.enrichedGOs <- read.delim("topGO/UniProt_GOs/bkgr/bkg.gua_mar.gua_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

bkg.gua_mar.mar.enrichedGOs <- read.delim("topGO/UniProt_GOs/bkgr/bkg.gua_mar.mar_topGO_results_BP.txt",
                                          sep = '\t',
                                          header = TRUE)

goTermPlot(bkg.eur_gua.eur.enrichedGOs[1:5,])
goTermPlot(bkg.eur_gua.gua.enrichedGOs[1:5,])
goTermPlot(bkg.eur_mar.eur.enrichedGOs[1:5,])
goTermPlot(bkg.eur_mar.mar.enrichedGOs[1:5,])
goTermPlot(bkg.gua_mar.gua.enrichedGOs[1:5,])
goTermPlot(bkg.gua_mar.mar.enrichedGOs[1:5,])

# GO terms (BP) de efecto de tratamiento dentro de cada variante

shade.eur.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/shade.eur.up_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

shade.eur.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/shade.eur.down_topGO_results_BP.txt",
                                         sep = '\t',
                                         header = TRUE)

uvb.eur.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/uvb.eur.up_topGO_results_BP.txt",
                                     sep = '\t',
                                     header = TRUE)

uvb.eur.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/uvb.eur.down_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

shade.gua.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/shade.gua.up_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

shade.gua.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/shade.gua.down_topGO_results_BP.txt",
                                         sep = '\t',
                                         header = TRUE)

uvb.gua.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/uvb.gua.up_topGO_results_BP.txt",
                                     sep = '\t',
                                     header = TRUE)

uvb.gua.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/uvb.gua.down_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

shade.mar.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/shade.mar.up_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

shade.mar.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/shade.mar.down_topGO_results_BP.txt",
                                         sep = '\t',
                                         header = TRUE)

uvb.mar.up.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/uvb.mar.up_topGO_results_BP.txt",
                                     sep = '\t',
                                     header = TRUE)

uvb.mar.down.enrichedGOs <- read.delim("topGO/UniProt_GOs/withinVar/uvb.mar.down_topGO_results_BP.txt",
                                       sep = '\t',
                                       header = TRUE)

goTermPlot(shade.eur.up.enrichedGOs[1:5,])
goTermPlot(shade.gua.up.enrichedGOs[1:5,])
goTermPlot(shade.mar.up.enrichedGOs[1:5,])
goTermPlot(uvb.eur.up.enrichedGOs[1:5,])
goTermPlot(uvb.gua.up.enrichedGOs[1:5,])
goTermPlot(uvb.mar.up.enrichedGOs[1:5,])


#### GO TERM ANALYSIS PARA MÓDULOS ####

# Comprobar la significancia de cada módulo para cada comparación. Usar lista de geneIDs y FDRs con la mayor significancia.
# Cargar librerías y funciones de topGO.

# Importar universod de GO terms con topGO::readMappings()

geneID2GO <- readMappings(file = "/home/rmblazquez/Documentos/Resultados/RNAseq_acebuche/topGO/Annotations/GOterm_universe_UniProt.csv")

# Importar lista de genes por módulo de WGCNA

genesByModule <- read.table("WGCNA_genesByModule.txt", sep = "\t", header = TRUE)

# seleccionamos genes de módulos 3, 10 y 4 en los resultados de comparación control vs. sombra (3 y 10) y uvb (4)
FDRs.shade.M3 <- res.int.shade[rownames(res.int.shade) %in% rownames(genesByModule[which(genesByModule$X2 == "brown"),]),]
FDRs.shade.M10 <- res.int.shade[rownames(res.int.shade) %in% rownames(genesByModule[which(genesByModule$X2 == "purple"),]),]
FDRs.uvb.M4 <- res.int.uvb[rownames(res.int.uvb) %in% rownames(genesByModule[which(genesByModule$X2 == "yellow"),]),]

# Generamos la lista de comparaciones a analizar
goterm.list.m3 <- list(shade.M3 = FDRs.shade.M3)
goterm.list.m10 <- list(shade.M10 = FDRs.shade.M10)
goterm.list.m4 <- list(uvb.M4 = FDRs.uvb.M4)

# Realizamos los análisis con topGOloop
setwd("/home/rmblazquez/Documentos/Resultados/RNAseq_acebuche/")
topGOloop(goterm.list.m3)
topGOloop(goterm.list.m10)
topGOloop(goterm.list.m4)

# GO terms (BP) de efecto general de tratamiento
shade.M3.enrichedGOs <- read.delim("topGO/UniProt_GOs/modules/shade.M3_topGO_results_BP.txt",
                                   sep = '\t',
                                   header = TRUE)
shade.M10.enrichedGOs <- read.delim("topGO/UniProt_GOs/modules/shade.M10_topGO_results_BP.txt",
                                    sep = '\t',
                                    header = TRUE)
uvb.M4.enrichedGOs <- read.delim("topGO/UniProt_GOs/modules/uvb.M4_topGO_results_BP.txt",
                                 sep = '\t',
                                 header = TRUE)

goTermPlot(shade.M3.enrichedGOs[1:20,])
goTermPlot(shade.M10.enrichedGOs[1:20,])
goTermPlot(uvb.M4.enrichedGOs[1:20,])

