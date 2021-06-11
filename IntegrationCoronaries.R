#BiocManager::install(c("ReactomePA","biomaRt","clusterProfiler"))

library(umap)
library(Seurat)
library(gplots)
library(limma)
#library(ReactomePA)
library(dplyr)
#library(biomaRt)
#library(clusterProfiler)
#library(org.Hs.eg.db)
#library(reactome.db)
#library(SingleR)
#library(ggplot2)
#library(ggfortify)
#library(survival)
#library(tidyverse)
#library(tableone)
#library(beeswarm)
#library(pROC)

#reticulate::py_install(packages ='umap-learn')
##############################################

copycat = 1 # copycat object on and off - it creates seurat objects with same content but different gene names (ensembl, HGCN), as it is difficult to change it in existing seurat object

#stamp for filenames from this analysis
a_name = "16_with_clint_rib_min_scaled_withGN_unscaled_PCAs_QN"

setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/analysis_draft")

raw.genecounts = read.table (file = "raw_counts.txt.minRib.txt.PC.txtt",header = T,sep = "\t",row.names = 1)

setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/Clint")

raw.genecounts2 = read.table (file = "genecountsraw_no_aorta.txt.PC.txt.minRib2.txt",header = T,sep = "\t",row.names = 2)
raw.genecounts2 = raw.genecounts2[,-1]
setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/analysis_draft")

#correct for UMI saturation the log = ln!!!! no log10 no log2 its natural log
raw.genecounts=round(-4096*(log(1-(raw.genecounts/4096))))

#exclude spike-ins
raw.genecounts=raw.genecounts[grep("ERCC",row.names(raw.genecounts),invert=TRUE),]
#exclude unvanted genes
raw.genecounts=raw.genecounts[grep("UGDH.AS1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM2P2",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("LOC100131257",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("KCNQ1OT1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("MALAT1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("MAB21L3",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("EEF1A1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]
#
#exclude spike-ins
raw.genecounts2=raw.genecounts2[grep("ERCC",row.names(raw.genecounts2),invert=TRUE),]
#exclude unvanted genes
raw.genecounts2=raw.genecounts2[grep("UGDH.AS1",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("PGM2P2",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("LOC100131257",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("KCNQ1OT1",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("MALAT1",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("PGM5P2",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("MAB21L3",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("EEF1A1",row.names(raw.genecounts2),invert=TRUE),]
raw.genecounts2=raw.genecounts2[grep("PGM5P2",row.names(raw.genecounts2),invert=TRUE),]
#


#Scale before quantile normalization
raw.genecounts=t(t(raw.genecounts)/colSums(raw.genecounts))*100000
mn = mean(colSums(raw.genecounts2))
raw.genecounts2=t(t(raw.genecounts2)/colSums(raw.genecounts2))*mn

#quantile normalize
raw.genecounts=round(limma::normalizeQuantiles(raw.genecounts))
raw.genecounts2=round(limma::normalizeQuantiles(raw.genecounts2))

#create seurat object
colon <- CreateSeuratObject(raw.genecounts, min.cells = 1,min.features = 1,
                            project = "AE")

raw.genecountsHG=raw.genecounts
namesHG=sapply( strsplit( row.names(raw.genecountsHG), "_ENSG"), "[", 1)
 row.names(raw.genecountsHG)=namesHG
  
colonHG <- CreateSeuratObject(raw.genecountsHG, min.cells = 1,min.features = 1,
                                 project = "AE")
  

colon2 <- CreateSeuratObject(raw.genecounts2, min.cells = 1,min.features = 1,
                            project = "Clint")


colon.combined <- merge(colonHG, y = colon2, add.cell.ids = c("AE", "Clint"), project = "AE_Clint")


#normalize the counts

colon.list <- SplitObject(colon.combined, split.by = "orig.ident")
colon.list <- colon.list[c("AE", "Clint")]


for (i in 1:length(colon.list)) {
  colon.list[[i]] <- NormalizeData(colon.list[[i]], verbose = FALSE)
  colon.list[[i]] <- FindVariableFeatures(colon.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

reference.list <- colon.list[c("AE", "Clint")]
colon.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30,k.filter=150)
colon.integrated <- IntegrateData(anchorset = colon.anchors, dims = 1:30)





#select fixed number of variable genes

colon.integrated<- FindVariableFeatures(colon.integrated, selection.method = "vst", nfeatures = 5000)

top10 <- head(VariableFeatures(colon.integrated), 10)



#scale data based on all genes
all.genes <- rownames(colon.integrated)
colon.integrated <- ScaleData(colon.integrated, features = c(all.genes))


#runPCAS
colon.integrated <- RunPCA(colon.integrated, features = VariableFeatures(object = colon.integrated))
VizDimLoadings(colon.integrated, dims = 1:2, reduction = "pca")
DimPlot(colon.integrated, reduction = "pca")
DimHeatmap(colon.integrated, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(colon.integrated, dims = 1:9, cells = 500, balanced = TRUE)
colon.integrated <- JackStraw(colon.integrated, num.replicate = 100)
colon.integrated <- ScoreJackStraw(colon.integrated, dims = 1:20)
JackStrawPlot(colon.integrated, dims = 1:20)
ElbowPlot(colon.integrated)



colon.integrated <- NormalizeData(object = colon.integrated, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
    
  



#find clusters and make tSNE
colon.integrated <- FindNeighbors(colon.integrated, dims = 1:12)
colon.integrated <- FindClusters(colon.integrated, resolution = 0.5)
head(Idents(colon.integrated), 10)
colon.integrated <- RunTSNE(object = colon.integrated, dims.use = 1:12, do.fast = TRUE)
colon.integrated <- RunUMAP(object = colon.integrated, dims = 1:12)

DimPlot(colon.integrated, reduction = "umap")
DimPlot(colon.integrated, reduction = "umap",split.by = "orig.ident" )
DimPlot(colon.integrated, reduction = "pca")
DimPlot(colon.integrated, reduction = "tsne")
TSNEPlot(object = colon.integrated)



#find marker genes
colon.markers <- FindAllMarkers(object = colon.integrated, only.pos = TRUE, min.pct = 0.2, 
                                thresh.use = 0.2)



cor.matrix = cor(as.matrix(colon.integrated@assays$integrated@data))
colon.integrated[["mean.cor"]] <- colMeans(cor.matrix)
aa <- cor.matrix[order(row(cor.matrix), -cor.matrix)]
second.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 2])
tenth.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 10])
names(second.best)=row.names(cor.matrix)
colon.integrated[["second.best.cor"]] <- second.best
names(tenth.best)=row.names(cor.matrix)
colon.integrated[["tenth.best.cor"]] <- tenth.best

write.table(cor.matrix, file = "cor.matrix.txt", sep = "\t", quote = F)

