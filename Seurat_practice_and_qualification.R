library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggraph)
library(clustree)


{
setwd("D:/QSB-academic-courses/ESAM472-RNA-sequencing/Assignment5/opt5_data")
dir_filtered<-'D:/QSB-academic-courses/ESAM472-RNA-sequencing/Assignment5/opt5_data';
#gene.column=2 selects gene names from the file features.tsv.gz, instead of gene-ids
#Note: all 3 files in "filter" directory should be gzipped
sc.counts <-Read10X(data.dir=dir_filtered,gene.column=2 , unique.features = TRUE)
#the input of sc.counts should be dgc matrix, some formats of data attached in papers should be transformed to this 
sc<-CreateSeuratObject(counts = sc.counts , project = "marrow",min.cells = 3 , min.features = 100)
#
sc[["percent.mt"]]<-PercentageFeatureSet(sc,pattern ="^MT-")
sc[["percent.ribo"]]<-PercentageFeatureSet(sc,pattern ="^RP[SL][[:digit:]]")
VlnPlot(sc,features= c("nFeature_RNA" ,"nCount_RNA" ,"percent.mt" ,"percent.ribo"),ncol=4,log = TRUE)
#copy picture to .pdf
dev.copy2pdf(file ="vioplot_of_ribo_and_mt_rna.pdf" , width =10 , height =8)

#thresholds for filter out outiers
#Adjust the thresholding parameters until you feel like you have good values
min_genes <- 100
max_genes <- 10000
min_reads <- 100
max_reads <- 50000
#scatter plot of counts and genes in log scale, the for lines are thresholds, each dot presents a cell
FeatureScatter(sc, feature1="nCount_RNA" , feature2 = "nFeature_RNA") +
  scale_y_log10() + scale_x_log10() +
  geom_hline(yintercept = min_genes) +
  geom_hline(yintercept = max_genes) +
  geom_vline(xintercept = min_reads) +
  geom_vline(xintercept = max_reads)
#save picture as eps
dev.copy2eps(file ="scatter_counts_features_each_cell.eps" , width =10 , height =8)
#the mt threshold selected from the results above
mt.thresh=30
sc<-subset(x = sc, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes
           & nCount_RNA > min_reads & nCount_RNA < max_reads & percent.mt < mt.thresh)

}

#use sctransform to visualize data (sctransform use negativebinomial residual regression) with respect to 
#library size and mitochondrial RNA counts
library(sctransform)
sc <-SCTransform ( sc, vars.to.regress = c("percent.mt") , verbose = T )

#shuffle the PCs across cells (to distory the inner correlation) to select good PC numbers
{
  library(Matrix)
  nshuffle<-25
  n_evs=100
  evnull <- rep(0 , n_evs)
  library(irlba)
  shuffled_expression<-as.matrix(sc@assays$SCT@scale.data);
  L <-irlba(shuffled_expression ,n_evs)
  evs<-L$d
  for (each in 1:nshuffle) {
    cat(paste0(each ,"...\ n ") )
    sct_tmp<-t(apply(shuffled_expression, 1,sample,replace=FALSE));
    L_tmp<-irlba(sct_tmp,n_evs)
    evnull<-evnull+L_tmp$d
    
  }
  evnull=evnull/nshuffle
  
}
plot(1:n_evs,evs,xlab="Top PC index", ylab="Singular Values", 
     main="Singular Values of PCs in shuffled and unshuffled data", col="blue")
lines(1:n_evs, evnull, type="o", col="green")
dev.copy2pdf(file ="selectPCs1-100.pdf" , width =10 , height =8)


#PCA after data normalization
sc <- RunPCA(sc, verbose = T,npcs =100)
#Find nearest neibourhood,k.param is the number of nearest neighbors
resolution_<-c(0.1 ,0.15 ,0.2 ,0.25 ,0.3 ,0.4 ,0.5)
library(ggraph)
library(clustree)
sc<-FindNeighbors(sc,graph.name="mgraph20",dims=1:100,verbose=TRUE,k.param=20)
sc <- FindClusters(object=sc , graph.name ="mgraph20" , verbose = T ,
                   resolution = resolution_ , algorithm =2)
#plot the clustertree under different resolution
clustree(sc,prefix="mgraph20_res.")
dev.copy2pdf(file ="clustertree_kparam_20.pdf" , width =10 , height =8)

sc<-FindNeighbors(sc,graph.name="mgraph40",dims=1:100,verbose=TRUE,k.param=40)
sc <- FindClusters(object=sc , graph.name ="mgraph40" , verbose = T ,
                   resolution = resolution_ , algorithm =2)
clustree(sc,prefix="mgraph40_res.")
dev.copy2pdf(file ="clustertree_kparam_40.pdf" , width =10 , height =8)

sc<-FindNeighbors(sc,graph.name="mgraph80",dims=1:100,verbose=TRUE,k.param=80)
sc <- FindClusters(object=sc , graph.name ="mgraph80" , verbose = T ,
                   resolution = resolution_ , algorithm =2)
clustree(sc,prefix="mgraph80_res.")
dev.copy2pdf(file ="clustertree_kparam_80.pdf" , width =10 , height =8)

sc<-FindNeighbors(sc,graph.name="mgraph120",dims=1:100,verbose=TRUE,k.param=120)
sc <- FindClusters(object=sc , graph.name ="mgraph120" , verbose = T ,
                   resolution = resolution_ , algorithm =2)
clustree(sc,prefix="mgraph120_res.")
dev.copy2pdf(file ="clustertree_kparam_120.pdf" , width =10 , height =8)

sc<-FindNeighbors(sc,graph.name="mgraph200",dims=1:100,verbose=TRUE,k.param=200)
sc <- FindClusters(object=sc , graph.name ="mgraph200" , verbose = T ,
                   resolution = resolution_ , algorithm =2)
clustree(sc,prefix="mgraph200_res.")
dev.copy2pdf(file ="clustertree_kparam_200.pdf" , width =10 , height =8)

#Fix the resolution and k.param
Idents(object = sc) <- "mgraph20_res.0.5"



sc <- FindClusters(object=sc , graph.name ="mgraph" , verbose = T ,
                     resolution = resolution_ , algorithm =2)



#
DimPlot(sc , reduction ="pca" , label=TRUE , pt.size =0.75)
dev.copy2pdf(file ="pc1and2.pdf" , width =10 , height =8)
DimPlot(sc,reduction ="pca" , dims = c (2,3) , label = TRUE , pt.size =0.75)
dev.copy2pdf(file ="pc2and3.pdf" , width =10 , height =8)
DimPlot(sc,reduction ="pca" , dims = c (3,4) , label = TRUE , pt.size =0.75)
dev.copy2pdf(file ="pc3and4.pdf" , width =10 , height =8)






#Selection of perplexity by visualization
{
sc <- RunTSNE (sc , reduction ="pca" , dims =1:100 , perplexity=30, method ="FIt-SNE" ,
                 max_iter =1000 , verbose = TRUE , reduction.name="tsne30" , reduction.key ="tSNE30_")
sc <- RunTSNE (sc , reduction ="pca" , dims =1:100 , perplexity=60, method ="FIt-SNE" ,
               max_iter =1000 , verbose = TRUE , reduction.name="tsne60" , reduction.key ="tSNE60_")
sc <- RunTSNE (sc , reduction ="pca" , dims =1:100 , perplexity=125, method ="FIt-SNE" ,
               max_iter =1000 , verbose = TRUE , reduction.name="tsne125" , reduction.key ="tSNE125_")
sc <- RunTSNE (sc , reduction ="pca" , dims =1:100 , perplexity=250, method ="FIt-SNE" ,
               max_iter =1000 , verbose = TRUE , reduction.name="tsne250" , reduction.key ="tSNE250_")
}
plot1 <- DimPlot (sc , reduction ="tsne30" , label = TRUE , pt.size =0.75)
plot2 <- DimPlot (sc , reduction ="tsne60" , label = TRUE , pt.size =0.75)
plot3 <- DimPlot (sc , reduction ="tsne125" , label = TRUE , pt.size =0.75)
plot4 <- DimPlot (sc , reduction ="tsne250" , label = TRUE , pt.size =0.75)
plot_grid ( plot1 , plot2 , plot3 , plot4 , labels = c ( 'perp 30'  , 'perp 60'  , 'perp 125' , 'perp 250') ,
            label_size = 12 , ncol =2)
dev.copy2pdf(file ="tsne_perplexity.pdf" , width =10 , height =8)


#identify marker genes
sc.markers.roc <- FindAllMarkers (sc,only.pos=TRUE , min.pct = 0.25 ,
                                       logfc.threshold = 0.25 , test.use ="roc")
sc.markers.roc.thresholded <- filter(sc.markers.roc , myAUC >=0.7)
sc.markers.roc.thresholded.grouped<- sc.markers.roc.thresholded %>% group_by (cluster)
sc.markers.roc.thresholded.grouped.sorted<-sc.markers.roc.thresholded.grouped %>%arrange(desc(myAUC),.by_group=TRUE)
sc.markers.roc.df <-as.data.frame (sc.markers.roc.thresholded.grouped.sorted)
write.csv (sc.markers.roc.thresholded.grouped.sorted , file ="marrow_cluster_markers.csv")

#test the performance of marker gene selection by testing 2 marker genes
FeaturePlot(sc, reduction ="tsne60" , features = c ("S100A9" ,"MALAT1") ,
              pt.size = 0.75 , ncol = 2 , order = TRUE,label = TRUE)
dev.copy2pdf(file ="tsne60_cluster0_2genes.pdf" , width =10 , height =8)
VlnPlot(sc, features=c("S100A9" ,"MALAT1") , ncol =2 , log=TRUE )
dev.copy2pdf(file ="vioplot_tsne60_cluster0_2genes.pdf" , width =10 , height =8)

FeaturePlot(sc, reduction ="tsne60" , features = c ("PPT1" ,"CST3") ,
            pt.size = 0.75 , ncol = 2 , order = TRUE,label = TRUE)
dev.copy2pdf(file ="tsne60_cluster11_2genes.pdf" , width =10 , height =8)
VlnPlot(sc, features=c("PPT1" ,"CST3") , ncol =2 , log=TRUE )
dev.copy2pdf(file ="vioplot_tsne60_cluster11_2genes.pdf" , width =10 , height =8)


FeaturePlot(sc, reduction ="tsne60" , features = c ("MZB1" ,"IGLC3") ,
            pt.size = 0.75 , ncol = 2 , order = TRUE,label = TRUE)
dev.copy2pdf(file ="tsne60_cluster14_2genes.pdf" , width =10 , height =8)
VlnPlot(sc, features=c("MZB1" ,"IGLC3") , ncol =2 , log=TRUE )
dev.copy2pdf(file ="vioplot_tsne60_cluster14_2genes.pdf" , width =10 , height =8)
