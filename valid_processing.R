total.umi<-apply(valid,2,sum);plot(total.umi,as.matrix(valid["GAPDH",]))
kmeans_original_data$total.umi<-total.umi
std.cell.valid<-apply(valid[which(rowSums(valid)>0),], 2, function(x){return(sd(x/sum(x)))})
kmeans_original_data$std.cell<-std.cell.valid
library(LICORS)
library(stringr)
normalized_valid<-data.frame(normalize(t(valid),byrow = TRUE))
colnames(normalized_valid)<-str_c(colnames(normalized_valid),"_probability")
{
  kmeans_original_cluster1<-kmeans_original_data%>%filter(cluster==1)
  kmeans_original_cluster2<-kmeans_original_data%>%filter(cluster==2)
  kmeans_original_cluster3<-kmeans_original_data%>%filter(cluster==3)
  kmeans_original_cluster4<-kmeans_original_data%>%filter(cluster==4)
}

ggplot(kmeans_original_data,aes(x=total.umi,y=GAPDH))+geom_point(aes(color=cluster))
ggplot(kmeans_original_data,aes(x=total.umi,y=LMNA))+geom_point(aes(color=cluster))
ggplot(kmeans_original_data%>%filter(cluster==4),aes(x=total.umi,y=CD40))+geom_point()+geom_smooth(method="loess")
ggplot(kmeans_original_data,aes(x=total.umi,y=CD40))+geom_point(aes(color=cluster))
ggplot(kmeans_original_data,aes(x=total.umi,y=FHL3))+geom_point(aes(color=cluster))
ggplot(kmeans_original_data,aes(x=total.umi,y=FOS))+geom_point(aes(color=cluster))

union(valid_merge$CD40_probability,valid_merge$CD40_probability)
0.0044406851
#prob std vs total cell UMI

ggplot(kmeans_original_data,aes(x=total.umi,y=std.cell.valid))+geom_point()

valid_merge<-cbind(kmeans_original_data,normalized_valid)
row.names(valid_merge)<-row.names(kmeans_original_data)
#prob of CD40 vs cell total UMI
plot(as.matrix(kmeans_original_data$total.umi),as.matrix(kmeans_original_data[,"CD40"]/sum(kmeans_original_data[,"CD40"])))

#
ggplot(valid_merge%>%select(c(CD40,CD40_probability,total.umi,cluster))%>%filter(CD40_probability>quantile(valid_merge$CD4_probability,0.98)),
       aes(x=total.umi,y=CD40))+geom_point(aes(color=cluster))
#plot the total.umi vs GAPDH with color indicating the estimated prob of GAPDH in each cell
ggplot(valid_merge%>%select(c(GAPDH,GAPDH_probability,total.umi,cluster)),
       aes(x=total.umi,y=GAPDH))+geom_point(aes(color=GAPDH_probability))

ggplot(valid_merge%>%select(c(GAPDH,GAPDH_probability,total.umi,cluster))%>%filter(GAPDH_probability<0.015),
       aes(x=total.umi,y=GAPDH))+geom_point(aes(color=GAPDH_probability))




ggplot(valid_merge%>%select(c(LMNA,LMNA_probability,total.umi,cluster)),
       aes(x=total.umi,y=LMNA))+geom_point(aes(color=LMNA_probability))

ggplot(valid_merge%>%select(c(GAPDH,GAPDH_probability,total.umi,cluster))%>%filter(GAPDH_probability<0.001),
       aes(x=total.umi,y=GAPDH))+geom_point(aes(color=GAPDH_probability))


ggplot(valid_merge%>%select(c(CD8B,CD8B_probability,total.umi,cluster))%>%filter(CD8B_probability<0.001),
       aes(x=total.umi,y=CD8B))+geom_point(aes(color=CD8B_probability))




ggplot(valid_merge%>%select(c(FTH1,FTH1_probability,total.umi,cluster))%>%filter(FTH1_probability<0.01),
       aes(x=total.umi,y=FTH1))+geom_point(aes(color=FTH1_probability))




ggplot(valid_merge%>%select(c(BAD,BAD_probability,total.umi,cluster)),
       aes(x=total.umi,y=BAD))+geom_point(aes(color=BAD_probability))



ggplot(valid_merge%>%select(c(CYS1,CYS1_probability,total.umi,cluster)),
       aes(x=total.umi,y=CYS1))+geom_point(aes(color=CYS1_probability))
