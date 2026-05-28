
#===========================================================#
#### Análisis de redes de co-expresión de genes en Olivo ####
#===========================================================#

library(WGCNA)
source('/home/rmblazquez/Documentos/Scripts/wgcna2igraph.R')
library(igraph)
library(dynamicTreeCut)
library(stats)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
library(pheatmap)

enableWGCNAThreads(nThreads = 12)

setwd("/home/rmblazquez/Documentos/Resultados/RNAseq_acebuche/WGCNA")

#### FUNCIONES ####

# Función para generar informes de media + SE

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm = TRUE),
      se = sd(x[[col]], na.rm = TRUE)/sqrt(length(x[[col]])))
  }
  data_sum <- ddply(data, groupnames, .fun = summary_func,
                    varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# Función para hacer plots de correlación GS vs GMM. Inputs:
#  gsdata (data frame con nombres de genes, asignación a módulo, y gene signification para una comparación), 
#  moduleNum (formato MX, X == número), 
#  moduleCol (formato = color)
corGS_GMM <- function(gsdata, moduleNum, moduleCol) {
  gsdata.m <- gsdata[which(gsdata$module == as.character(moduleNum)), ]
  exprGenes.m <- exprGenes[which(exprGenes$module == as.character(moduleCol)), ]
  i <- 1
  whichModule <- moduleEigen[, which(colnames(moduleEigen) == paste0("ME", as.character(moduleCol)))]
  gModMem <- c()
  while (i <= length(rownames(exprGenes.m))) {
    gModMem[i] <- cor(as.numeric(exprGenes.m[i, c(1:49)]), whichModule, method = "spearman") # PROBLEMA: moduleEigen$MEbrown!!!
    i <- i + 1
  }
  plot(as.numeric(gsdata.m[,2]), as.numeric(abs(gModMem)),
       xlab = paste0(as.character(moduleNum), " gene significance"), ylab = paste0(as.character(moduleNum), " |gene module membership|"), col = as.character(moduleCol))
}

# Función para calcular la gene significance de cada módulo para cada variable
# geneSig.list debe ser una variable de significación tipo "exprGenes$geneSig.xxx"
geneSig_calc <- function(geneSig.list) {
  gsdata <- as.data.frame(cbind(exprGenes$module, as.numeric(geneSig.list)))
  rownames(gsdata) <- rownames(exprGenes)
  colnames(gsdata) <- c("module", "significance")
  gsdata$module <- factor(gsdata$module, 
                          levels = blockColors)
  levels(gsdata$module) <- c("singlets", "M1", "M2", "M3", "M4", "M5", "M6", 
                             "M7", "M8", "M9", "M10", "M11", "M12",
                             "M13", "M14", "M15", "M16") # Ajustar al número de módulos
  gsdata$significance <- as.numeric(gsdata$significance)
  gsdata.dsum <- data_summary(gsdata, varname = "significance", groupnames = "module")[2:17,]
  colnames(gsdata.dsum)[2] <- deparse(substitute(geneSig.list))
  gsdata.dsum
}


#### IMPORTAR DATOS DE DESeq2 ####

# Importar conteos normalizados del script "RNAseq_olivo_expressionAnalysis.R"
# Usar los datos normalizados y filtrados del objeto `dds.all` (línea 452):

normCounts <- counts(dds.all, normalize = TRUE) # dds.all no está definido en este script, calcular con DESeq2

# Nos quedamos con los genes que tengan más de 2/3 de valores distintos a cero (o sea, los que tienen 16 ceros por gen)

normFilteredCounts <- normCounts[rowSums(normCounts == 0) <= 16, ]

# Eliminamos genes con varianza próxima a cero (0.05)

normFilteredCounts <- normFilteredCounts[apply(normFilteredCounts, 1, var) >= 0.05, ]

# Transformamos los conteos normalizados en log2 + 1 para atenuar las diferencias de magnitud

data.Expr <- t(log2(normFilteredCounts + 1))
data.Gene <- rownames(normFilteredCounts)
data.Sample <- colnames(normFilteredCounts) 
colnames(data.Expr) <- data.Gene
rownames(data.Expr) <- data.Sample

#### CALCULAR VALOR DE UMBRAL SUAVE ####

# Generamos un vector de valores de soft threshold (beta) a comprobar 

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(data.Expr, 
                         powerVector = powers, 
                         verbose = 5)
gc()

# Graficamos los valores 

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

# El valor de sft en la que la conectividad de la red y la independencia de la topología son óptimas está en torno a 8

beta1 <- 8

# Comprobamos cómo se correlaciona la expresión con la escala para ese valor de beta

k.dat.Expr <- softConnectivity(data.Expr, power = beta1) - 1	
scaleFreePlot(k.dat.Expr, main = paste("data Expression, power = ", beta1), truncated = F)
gc()


#### DETECCIÓN DE MÓDULOS EN LA RED: MÉTODO ESTÁTICO ####

# Filtramos genes con varianza de expresión de cero

kCut <- (dim(data.Expr)[[2]])
kRank <- rank(-k.dat.Expr)
vardat.Expr <- apply(data.Expr, 2, var) 
restk <- kRank <= kCut & vardat.Expr > 0 

# Calculamos la matriz de adyacencia
# Usar método de Spearman porque es menos sensible a los valores atípicos

ADJdat.Expr <- adjacency(datExpr = data.Expr[,restk], 
                         power = beta1, 
                         corOptions = "use = 'p', method = 'spearman'")
gc()

# Convertimos la matriz de adyacencia en matriz de disimilitud

dissTOMdat.Expr <- TOMdist(ADJdat.Expr)

# La representamos como un dendrograma con cluster jerárquico

hierTOMdat.Expr <- hclust(as.dist(dissTOMdat.Expr), method = "average")

#### DETECCIÓN DE MÓDULOS EN LA RED: MÉTODO HÍBRIDO ####

# Extraemos los módulos de la red
# Comprobamos cómo distintos valores de detectCutHeight, minModuleSize, 
# y deepSplit afectan a la generación de módulos
# maxBlockSize: idealmente debería ser el número total de genes, pero WGCNA limita el número máximo en función de la RAM
# Para 600 GB, el máximo es de 46360, por lo que si queremos usar 100 GB, podemos poner el límite en unos 7000 (https://support.bioconductor.org/p/84938/)

cor <- WGCNA::cor # Para evitar error provocado por conflicto de nombres entre stats::cor y WGCNA::cor, revertir tras análisis con blockwiseModules()

# Generar varias redes jugando con los parámetros "detectCutHeight", "deepSplit", y "minModuleSize"
# Matriz de combinaciones de parámetros

detectCutHeight.v <- rep(seq(0.8, 0.95, by = 0.05), each = 12)
deepSplit.v <- rep(seq(0, 4, by = 2), times = 4, each = 4)
minModuleSize.v <- rep(seq(50, 200, by = 50), times = 12)
bwm.params <- data.frame(dch = detectCutHeight.v, ds = deepSplit.v, mms = minModuleSize.v)

# Bucle que genera una red para cada combinación de parámetros especificada en la matriz de parámetros

data.ntwk.list <- list()
for (i in 1:nrow(bwm.params)){
  data.ntwk <- blockwiseModules(datExpr = data.Expr[,restk], 
                                maxBlockSize = 7000,          # Limitar número a 7000 por estima de uso de RAM
                                networkType = "unsigned",     # Unsigned para que la dirección de la sobre-expresión no se tenga en cuenta
                                power = beta1, 
                                detectCutHeight = bwm.params[i, 1], 
                                deepSplit = bwm.params[i, 2], 
                                minModuleSize = bwm.params[i, 3], 
                                saveTOMs = FALSE, 
                                verbose = F)
  data.ntwk.list[[i]] <- data.ntwk
}
gc()

# Incluir en el nombre de cada red (elemento de la lista) los parámetros utilizados)

list.names <- c()
for (i in 1:nrow(bwm.params)){
 list.names[i] <- paste0("detectCutHeight = ", 
			 bwm.params[i, 1], 
			 ", deepSplit = ", 
			 bwm.params[i, 2], 
			 ", minModuleSize = ", 
			 bwm.params[i,3])
}

names(data.ntwk.list) <- list.names

# Representación por método híbrido

module.numbers <- c()
for (i in 1:nrow(bwm.params)) {
  module.numbers <- c(module.numbers, length(table(data.ntwk.list[[i]]$colors)))
}

names(module.numbers) <- list.names

# Crear vector de colores del método estático con parámetros comparables al híbrido

colorhdat.Expr <- cutreeStaticColor(hierTOMdat.Expr, cutHeight = 0.95, minSize = 200)

# Graficar los módulos de las redes para compararlos en un dendrograma

pdf("WGCNA_modules_hybrid.pdf")
par(mfrow = c(3, 1), mar = c(2, 4, 1, 1))
plot(hierTOMdat.Expr, main = "Network Dendrogram", labels = F, xlab = "", sub = "")
plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = colorhdat.Expr)) #static
plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = data.ntwk.list[[48]]$colors)) #hybrid
# Usar este fragmento de código para plotear todos los módulos generados a partir de la matriz de parámetros
#for (i in 1:nrow(bwm.params)) {
#plotColorUnderTree(hierTOMdat.Expr, colors = data.frame(module = data.ntwk.list[[i]]$colors)) #hybrid
#}
dev.off() # los resultados con dch = 0.95 incluyen más genes por módulo, con deepSplit = 4 y mms = 200 saca 16 módulos de tamaños decentes [[48]]

# Relación de genes en cada módulo de la red seleccionada (modificar el índice i)

genesByModule <- data.frame(cbind(data.Gene[restk], data.ntwk.list[[48]]$colors))

write.table(genesByModule, file = "WGCNA_genesByModule.txt", sep = "\t") # en caso de querer exportar los datos

# Cargar desde tabla si no está cargado en el entorno:
# genesByModule <- read.table("WGCNA_genesByModule.txt", sep = "\t", header = TRUE)

# Representación de la red 
networkGraph <- 
wgcna2igraph(data.ntwk.list[[48]],      # red
             data.Expr[, restk],        # matriz de expresión
	     modules2plot = names(sort(table(genesByModule$X2), decreasing = TRUE))[-1], # Módulos a representar (eliminamos "grey" con el [-1])
             colors2plot = names(sort(table(genesByModule$X2), decreasing = TRUE))[-1],  # Colores de los nodos por pertenencia a módulo
	     #kME.threshold = 0.5,       # representar nodos com más de 0.5 de valor kME
	     #adjacency.threshold = 0.1, # representar ejes con adyacencia mayor a 0.1
             adj.power = 8,             # exponente de la transformación a adjyacencia
	     edge.alpha = 0.1,
	     node.size = 1
)

# Eliminar loops (edges que vuelven a su nodo de origen) con simplify(), eliminar nodos no conectados con delete_vertices(degree()==0)

networkGraph.clean <- delete_vertices(simplify(networkGraph), degree(networkGraph)==0) 

plot(networkGraph.clean)
# GENERA UN GRÁFICO GIGANTE Y NO INFORMATIVO!!! REPASAR COLORES O PROBAR CON https://ramellose.github.io/networktutorials/wgcna.html

#### EIGENGENE STATISTICS #### 

# Añadir a la tabla de eigengenes (MEs) las variables

data.ntwk <- data.ntwk.list[[48]] # Cambiar el índice según los módulos que se hayan seleccionado
# 48: "detectCutHeight = 0.95, deepSplit = 4, minModuleSize = 200"

data.ntwk$MEs$treatment <- as.factor(design$treatment)
data.ntwk$MEs$variant <- as.factor(design$variant)
data.ntwk$MEs$time <- as.factor(design$time)
write.table(data.ntwk$MEs, file = "MEs_values_for_stats.txt", sep = "\t", quote = F)

# Ejecutar esta orden cuando los eigengenes no están cargados en R
# data.ntwk <- list()
# data.ntwk$MEs <- read.table("MEs_values_for_stats.txt", sep = '\t', header = T, row.names = 1)

# Ejecutamos el PCA
rownames(data.ntwk$MEs) <- rownames(data.Expr)
MEpca <- prcomp(data.ntwk$MEs[,1:16])
s <- summary(MEpca) # Varianza explicada: PC1 = 27.0%, PC2 = 18.7%

# Generamos un vector de colores

color.group <- as.factor(paste0(data.ntwk$MEs$treatment, data.ntwk$MEs$variant, data.ntwk$MEs$time))
levels(color.group) <- c('seagreen', 'green', 'purple', 'blue', 'orange', 'red',
                         'seagreen', 'green', 'purple', 'blue', 'orange', 'red',
                         'seagreen', 'green', 'purple', 'blue', 'orange', 'red')

# Biplot de los eigengenes

plot(MEpca$x[,1], MEpca$x[,2], 
     xlab = paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), 
     ylab = paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), 
     pch = 21, col = "black", bg = color.group, cex = 1, las = 1, asp = 1)
abline(v = 0, lty = 2, col = "grey50")
abline(h = 0, lty = 2, col = "grey50")
l.x <- MEpca$rotation[,1] * 0.9
l.y <- MEpca$rotation[,2] * 0.9
arrows(x0 = 0, x1 = l.x, y0 = 0, y1 = l.y, col = "red", length = 0.15, lwd = 1.5)
l.pos <- l.y
lo <- which(l.y < 0)
hi <- which(l.y > 0)
l.pos <- replace(l.pos, lo, "1")
l.pos <- replace(l.pos, hi, "3")
text(l.x, l.y, labels = row.names(MEpca$rotation), col = "red", pos = l.pos)

# Cargamos los nombres ("MEcolores") de los módulos

MEcolors <- colnames(data.ntwk$MEs[1:16]) # range from 1 to total number of modules

# Hacemos pruebas de Kruskal-Wallis para cada módulo

MEcounter <- 1
kw.results <- c()
while (MEcounter <= length(MEcolors)){
  Treatment.kw <- kruskal.test(as.formula(paste0("data.ntwk$MEs$", MEcolors[MEcounter], " ~ data.ntwk$MEs$treatment")))
  Variant.kw <- kruskal.test(as.formula(paste0("data.ntwk$MEs$", MEcolors[MEcounter], " ~ data.ntwk$MEs$variant")))
  Time.kw <- kruskal.test(as.formula(paste0("data.ntwk$MEs$", MEcolors[MEcounter], " ~ data.ntwk$MEs$time")))
  MEvalues <- as.data.frame(data.ntwk$MEs[MEcounter])
  my_title <- MEcolors[MEcounter]
  my_vector <- c(as.character(my_title), 
		 as.character(Treatment.kw$data.name), 
		 as.numeric(Treatment.kw$statistic), 
		 as.numeric(Treatment.kw$p.value), 
		 as.numeric(Variant.kw$statistic), as.numeric(Variant.kw$p.value), 
		 as.numeric(Time.kw$statistic), as.numeric(Variant.kw$p.value))
  kw.results <- rbind(kw.results, my_vector)
  MEcounter <- MEcounter + 1
}

rownames(kw.results) <- kw.results[,1]
colnames(kw.results) <- c("module", "contrast", "treatment.stat", "treatment.pvalue", "variant.stat", "variant.pvalue", "time.stat", "time.pvalue")
write.table(kw.results, "ME_KWtest_results.txt", sep = "\t", quote = FALSE, rownames = FALSE)

# Generamos boxplots para todos los módulos

data.ntwk$MEs$group <- as.factor(paste0(data.ntwk$MEs$treatment, data.ntwk$MEs$variant, data.ntwk$MEs$time))

MEcounter <- 1
pdf(file = "MEs_boxplots.pdf")
while (MEcounter <= length(data.ntwk$MEs[1:16])){
  bxplt <- ggboxplot(data.ntwk$MEs, x = "treatment", y = MEcolors[MEcounter], color = "group", 
                     repel = T, font.label = list(size = 14, face = "plain"), 
                     add = "jitter", shape = "variant",
                     ggtheme = theme_gray())
  print(bxplt)
  MEcounter <- MEcounter + 1
}
dev.off()

# Correlación entre módulos

cor <- stats::cor # Volver a señalar cor como en el paquete de stats para no usar cor de WGCNA

MEcor <- cor(as.data.frame(data.ntwk$MEs[1:16])) # hay 16 módulos

pheatmap(MEcor, color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100))


#### GENE SIGNIFICANCE (-log(FDR)) & GENE MODULE MEMBERSHIP (cor(gene, ME)) ####

# Relacionamos conteos normalizados, identificadores de genes, y módulos
# Primero, cargamos los conteos normalizados de todos los genes analizados en WGCNA

# normFilteredCounts <- counts(dds.all, normalize = TRUE)

exprGenes <- as.data.frame(normFilteredCounts) # los conteos normalizados de dds.all deben estar aquí cargados

# Incluír columna de módulos en los conteos

exprGenes$module <- genesByModule$X2

# Ejecutar script RNAseq_olivo_expressionAnalysis.R en el entorno, 
# para obtener las variables de resultados de DESeq2 de los análisis dds.int, dds.bkg
# y seleccionar solo los genes incluidos en el WGCNA

int.shade <- res.int.shade[rownames(res.int.shade) %in% rownames(exprGenes), ]
int.uvb <- res.int.uvb[rownames(res.int.uvb) %in% rownames(exprGenes), ]

bkg.eur_gua <- res.bkg.eur_gua[rownames(res.bkg.eur_gua) %in% rownames(exprGenes), ]
bkg.eur_mar <- res.bkg.eur_mar[rownames(res.bkg.eur_mar) %in% rownames(exprGenes), ]
bkg.gua_mar <- res.bkg.gua_mar[rownames(res.bkg.gua_mar) %in% rownames(exprGenes), ]

int.shade.eur <- res.int.shade.eur[rownames(res.int.shade.eur) %in% rownames(exprGenes), ]
int.shade.gua <- res.int.shade.gua[rownames(res.int.shade.gua) %in% rownames(exprGenes), ]
int.shade.mar <- res.int.shade.mar[rownames(res.int.shade.mar) %in% rownames(exprGenes), ]

int.uvb.eur <- res.int.uvb.eur[rownames(res.int.uvb.eur) %in% rownames(exprGenes), ]
int.uvb.gua <- res.int.uvb.gua[rownames(res.int.uvb.gua) %in% rownames(exprGenes), ]
int.uvb.mar <- res.int.uvb.mar[rownames(res.int.uvb.mar) %in% rownames(exprGenes), ]


## Calcular la gene significance por módulo (la media del valor -log2(FDR) por módulo)

# Generamos listas de FDRs por contraste de hipótesis

exprGenes$geneSig.int.shade <- -log(as.numeric(int.shade$padj) + 0.001, 2) 
exprGenes$geneSig.int.uvb <- -log(as.numeric(int.uvb$padj) + 0.001, 2) 

exprGenes$geneSig.int.shade.eur <- -log(as.numeric(int.shade.eur$padj) + 0.001, 2) 
exprGenes$geneSig.int.shade.gua <- -log(as.numeric(int.shade.gua$padj) + 0.001, 2) 
exprGenes$geneSig.int.shade.mar <- -log(as.numeric(int.shade.mar$padj) + 0.001, 2) 
exprGenes$geneSig.int.uvb.eur <- -log(as.numeric(int.uvb.eur$padj) + 0.001, 2) 
exprGenes$geneSig.int.uvb.gua <- -log(as.numeric(int.uvb.gua$padj) + 0.001, 2) 
exprGenes$geneSig.int.uvb.mar <- -log(as.numeric(int.uvb.mar$padj) + 0.001, 2) 

# Generamos un vector de nombres de módulos ordenados por el número de genes que incluyen

blockColors <- names(sort(table(exprGenes$module), decreasing = T))
# blockColors[1:2] <- c("grey", "turquoise") # Ejecutar cuando M1 tenga más secuencias que los singlets (grey)

# Calculamos la significación de cada módulo para cada contraste:

# para sombra 
gsdata.int.shade.dsum <- geneSig_calc(exprGenes$geneSig.int.shade)

# para UVB
gsdata.int.uvb.dsum <- geneSig_calc(exprGenes$geneSig.int.uvb)


# Graficar los valores de Gene Significance por módulo

# para sombra
gsdata.int.shade.plot <- ggplot(gsdata.int.shade.dsum, aes(x = module, y = geneSig.int.shade, fill = module)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = geneSig.int.shade - se, ymax = geneSig.int.shade + se), width=.2,
                position = position_dodge(.9)) +
  scale_fill_manual("Module", values = blockColors[2:17])

# para UVB
gsdata.int.uvb.plot <- ggplot(gsdata.int.uvb.dsum, aes(x = module, y = geneSig.int.uvb, fill = module)) + 
  geom_bar(stat = "identity", color = "black", 
           position = position_dodge()) +
  geom_errorbar(aes(ymin = geneSig.int.uvb - se, ymax = geneSig.int.uvb + se), width=.2,
                position = position_dodge(.9)) +
  scale_fill_manual("Module", values = blockColors[2:17])

grid.arrange(gsdata.int.shade.plot, gsdata.int.uvb.plot,  ncol = 1, nrow = 2)


## Calcular "gene module membership" (correlación entre expresión de genes y valor del eigengen del módulo)

moduleEigen <- read.table("MEs_values_for_stats.txt", sep = '\t', header = T) # exportado de data.net$MEs, cargar de la tabla si no está ya cargado en el entorno

# Módulo 3 (brown): tiene alta gene significance en sombra

corGS_GMM(gsdata.int.shade, "M3", "brown")

# Módulo 10 (purple): alta gene significance en UVB

corGS_GMM(gsdata.int.uvb, "M10", "purple")

