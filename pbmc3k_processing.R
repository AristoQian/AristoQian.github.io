{
library(Matrix)
matrix_dir = "/Users/cqq7269/Downloads/filtered_gene_bc_matrices/GRCh38/pbmc3k/filtered_gene_bc_matrices/hg19/"   
barcode.path <- paste0(matrix_dir, "barcodes.tsv")
features.path <- paste0(matrix_dir, "genes.tsv")
matrix.path <- paste0(matrix_dir, "matrix.mtx")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
}

pca_pbmc3k<-prcomp(as.matrix(t(pbmc3k[which(rowSums(pbmc3k)>0),])),center =TRUE, scale.=TRUE)
#calculating PC variance contribution
std<-pca_pbmc3k$sdev
val<-std^2
proportion_val<-val/sum(val)
culmulative_proportion<-cumsum(proportion_val)
{
  par(mar=c(6,6,2,2))
  plot(val[1:200],type="b",
       cex=1,
       cex.lab=1,
       cex.axis=1,
       lty=1,
       lwd=1,
       xlab="PC index",
       ylab="Variance Contribution")
}

#seleting 100 PCs
data<-pca_pbmc3k$x[,1:100]
wss<-sapply(1:50, function(x){
  return(kmeans(data,centers = x)$tot.withinss)
})

plot(1:50,wss)
#set 20 as clustering number
km_pbmc3k<-kmeans(as.matrix(pca_pbmc3k$x[,1:100]),20,iter.max = 20,nstart = 20)
library(Rtsne)
pbmc3k_tsne<-Rtsne(data,dims = 2,pca = FALSE,max_iter = 1200,theta = 0.15,perplexity = 30,verbose = T)

#visualization
pbmc3k_kmdata<-cbind(pbmc3k_tsne$Y,km_pbmc3k$cluster)
vis_pbmc3k_km<-data.frame(cbind(pbmc3k_tsne$Y,str_c("cluster",as.character(km_pbmc3k$cluster))))
colnames(vis_pbmc3k_km)<-c("tSNE1","tSNE2","Cluster")
library(tidyverse)
write_csv(vis_pbmc3k_km,"vis_pbmc3k_km.csv");vis_pbmc3k_km<-read_csv("vis_pbmc3k_km.csv")
plot_km_tsne_pbmc3k<-ggplot(vis_pbmc3k_km,aes(x=tSNE1, y=tSNE2,color=Cluster))+geom_point()
plot_km_tsne_pbmc3k

#---------------------------------------------------------------------------------------
#ratio analysis
kmeans_original_data<-data.frame(t(pbmc3k));kmeans_original_data$cluster<-as.factor(km_pbmc3k$cluster)
total.umi<-apply(pbmc3k,2,sum)
kmeans_original_data$total.umi<-total.umi
std.cell.pbmc3k<-apply(pbmc3k[which(rowSums(pbmc3k)>0),], 2, function(x){return(sd(x/sum(x)))})
kmeans_original_data$std.cell<-std.cell.pbmc3k
library(LICORS)
library(stringr)
normalized_pbmc3k<-data.frame(normalize(t(pbmc3k),byrow = TRUE))
colnames(normalized_pbmc3k)<-str_c(colnames(normalized_pbmc3k),"_probability")
pbmc3k_merge<-cbind(kmeans_original_data,normalized_pbmc3k)
row.names(pbmc3k_merge)<-row.names(kmeans_original_data)
library(ggplot2)
ggplot(pbmc3k_merge%>%select(c(ENSG00000181163,ENSG00000181163_probability,total.umi,cluster))%>%filter(ENSG00000181163_probability<0.005),
       aes(x=total.umi,y=ENSG00000181163))+geom_point(aes(color=ENSG00000181163_probability))

ggplot(pbmc3k_merge%>%select(c(ENSG00000120129,ENSG00000120129_probability,total.umi,cluster))%>%filter(ENSG00000120129_probability<0.01),
       aes(x=total.umi,y=ENSG00000120129))+geom_point(aes(color=ENSG00000120129_probability))

ggplot(pbmc3k_merge%>%select(c(ENSG00000113732,ENSG00000113732_probability,total.umi,cluster))%>%filter(ENSG00000113732_probability<0.005),
       aes(x=total.umi,y=ENSG00000113732))+geom_point(aes(color=ENSG00000113732_probability))
