
library(umap)
library(Seurat)
library(limma) 



##############################################
setwd("C:/Users/micha/OneDrive/Documents/R/atheroexpress/analysis_draft")

#read in read count matrix files
raw.genecounts = read.table (file = "raw_counts.txt.minRib.txt.PC.txtt",header = T,sep = "\t",row.names = 1)

#correct for UMI saturation the log = ln!!!! no log10 no log2 its natural log; there are 4096 possible UMIs
raw.genecounts=round(-4096*(log(1-(raw.genecounts/4096))))

#check
raw.genecounts[1:10,1:10]

#exclude spike-ins
raw.genecounts=raw.genecounts[grep("ERCC",row.names(raw.genecounts),invert=TRUE),]
#exclude unvanted genes - historically those genes have mapping issues
raw.genecounts=raw.genecounts[grep("UGDH.AS1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM2P2",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("LOC100131257",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("KCNQ1OT1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("MALAT1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("MAB21L3",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("EEF1A1",row.names(raw.genecounts),invert=TRUE),]
raw.genecounts=raw.genecounts[grep("PGM5P2",row.names(raw.genecounts),invert=TRUE),]


#Scale before quantile normalization
raw.genecounts=t(t(raw.genecounts)/colSums(raw.genecounts))*100000

#quantile normalize
raw.genecounts=round(limma::normalizeQuantiles(raw.genecounts))

#create seurat object
plaque <- CreateSeuratObject(raw.genecounts, min.cells = 1,min.features = 1,
                            project = "AE")

#add correlation coefficients into the object 
cor.matrix = cor(log2(raw.genecounts+10))
plaque[["mean.cor"]] <- colMeans(cor.matrix)
aa <- cor.matrix[order(row(cor.matrix), -cor.matrix)]
second.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 2])
tenth.best= (matrix(aa, nrow(cor.matrix), byrow = TRUE)[, 10])
names(second.best)=row.names(cor.matrix)
plaque[["second.best.cor"]] <- second.best
names(tenth.best)=row.names(cor.matrix)
plaque[["tenth.best.cor"]] <- tenth.best



#get percentatge of mito genes
plaque[["percent.mt"]] <- PercentageFeatureSet(plaque, pattern = "^MT-")
percent.mt=PercentageFeatureSet(plaque, pattern = "^MT-")
VlnPlot(plaque, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(plaque, features = c("second.best.cor", "tenth.best.cor", "mean.cor"), ncol = 3)

plaque<- FindVariableFeatures(plaque, selection.method = "vst", nfeatures = 5000)
top <- head(VariableFeatures(plaque), 35)

#plot variable genes
plot1 <- VariableFeaturePlot(plaque)
plot1
plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
plot2


#scale data based on all genes
all.genes <- rownames(plaque)
plaque <- ScaleData(plaque, features = c(all.genes))


#runPCAS

plaque <- RunPCA(plaque, features = VariableFeatures(object = plaque))
VizDimLoadings(plaque, dims = 1:2, reduction = "pca")
DimPlot(plaque, reduction = "pca")
DimHeatmap(plaque, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(plaque, dims = 1:9, cells = 500, balanced = TRUE)
plaque <- JackStraw(plaque, num.replicate = 100)
plaque <- ScoreJackStraw(plaque, dims = 1:20)
JackStrawPlot(plaque, dims = 1:20)
ElbowPlot(plaque)

plaque <- NormalizeData(object = plaque, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
  


#find clusters and make tSNE
plaque <- FindNeighbors(plaque, dims = 1:12)
plaque <- FindClusters(plaque, resolution = 0.5)
head(Idents(plaque), 10)
plaque <- RunTSNE(object = plaque, dims.use = 1:12, do.fast = TRUE)
plaque <- RunUMAP(object = plaque, dims = 1:12)

#DimPlot(plaque)
DimPlot(plaque, reduction = "umap")
DimPlot(plaque, reduction = "pca")
DimPlot(plaque, reduction = "tsne")
TSNEPlot(object = plaque)


plaque.markers <- FindAllMarkers(object = plaque, only.pos = TRUE, min.pct = 0.2, 
                                thresh.use = 0.2)


