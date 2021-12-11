#win PC
install.packages('Seurat')
library(Seurat)

#linux 


#usage
library(stats)
library(base)
library(dplyr)
#create pbmc project
pbmc.data <- Read10X(data.dir = "D:/undergrad_research/pbmc3k/pbmc3k/pbmc3k/filtered_matrices_mex/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
#check the expression of particular genes in 1:60 cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:60]

#calculate MT QC and add it to metadata
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#check the fisrt 5 metadata
head(pbmc@meta.data, 5)

#violin plot of the three metadata objects
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#scatter plot of the 3 metadata objects
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#â€œLogNormalizeâ€? that normalizes the feature expression measurements for each cell by the total expression, 
#multipies this by a scale factor (10,000 by default), and log-transforms the result.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

#Find variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 15 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 15)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaling the data to :
#Shifts the expression of each gene, so that the mean expression across cells is 0
#Scales the expression of each gene, so that the variance across cells is 1
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:4, reduction = "pca")




#the structure of pbmc object
View(pbmc@commands$FindVariableFeatures)
View(pbmc@assays$RNA@var.features)
View(pbmc@meta.data$nFeature_RNA)
