
# title: "Proximal proteomics reveals a landscape of human nuclear condensates"
# author: "Ruofei Li (李若菲)"
# date: "05/13/2025"

  
# Load Required libraries

library(pheatmap)  
library(corrplot)
library(psych)

#### Fig 1d GO signature analysis


#The file "GO_AllLists.csv" contains GO analysis results of a bait or nuclear condensate (NC) MS proteome list, which was downloaded from Metascape (https://metascape.org/)

setwd("/data")
data<-read.csv("GO_AllLists.csv", sep = ",",  header = TRUE) 
data<-data[,c(3,4,1,17)] 
data<-data[which(data$Log.q.value.<= -2),] ## Screening for GO terms with q value <= 1×10⁻² 

#These two files manually curated GO terms into biological processes, for subsequent organization of bait/NC lists
GO_list<-read.csv("GO term.txt", sep = "\t",  header = TRUE) 
GO_category_order<-read.csv("GO category.txt", sep = "\t",  header = TRUE)


GO_category<-merge(GO_list,data,by=c(1:3),all.x = T)
GO_category$Log.q.value.<-abs(GO_category$Log.q.value.)

GO_category[is.na(GO_category[]==T)]<-0
GO_category$CategoryID<-factor(GO_category$CategoryID,levels =GO_category_order$CategoryID )

a_map<-GO_category[order(GO_category$CategoryID),c(1,2,4,5)]
a_map[4][a_map[4]>=20 & a_map[4]<30]<-24
a_map[4][a_map[4]>=30]<-30
rownames(a_map)<-a_map$GO

## Set the name of the bait/NC lists
names(a_map)[4]<-"my protein"

bk<-c(seq(0,19.9,by=0.1),seq(20,30,by=0.1))

p_col<-c(colorRampPalette(c("#0b68ff","#a0abfd"))(3),
         colorRampPalette(c("#6e50ec","#cf9fff"))(4),
         colorRampPalette(c("#c53207","#ffc16d"))(4),
         colorRampPalette(c("#f9f0a7","#9dd085","#42b4af"))(4)
)
names(p_col)<-GO_category_order$CategoryID
ann_col<-list(CategoryID=p_col)

##Using lg(q value) to generate a "barcode" of the bait/NC
p=pheatmap(a_map[4],cluster_rows = F,cluster_cols = F,
         color = c(colorRampPalette(c("white","skyblue1","slateblue2","purple4"))(length(bk)*2/3),colorRampPalette(c("purple4","firebrick3"))(length(bk)*1/3)),
         show_rownames = F,annotation_names_row =F,
         annotation_row = a_map[3],breaks = bk,
         annotation_colors = ann_col,width=7, height=11,
         filename = "/results/GO_barcode_heatmap.pdf"
)
write.table(a_map,file="/results/GO_barcode.txt",sep="\t",quote=FALSE)

#### Fig 2a,c and Extended Fig 3c, RA matrix
### Marker proteins of 18 NCs were performed PhastID, the enriched interactors with normalized abundance were documented in "NC_marker_matrix.txt"
NC_matrix<-read.csv("NC_marker_matrix.txt", sep = "\t",  header = TRUE)

#Calculate pairwise RA scores between baits to generate a matrix 
discoPlot1<-function(x){
  x<-x[,-1]
  x<-x^0.5   ##Abundance transformation to compress high-abundance variance and enhance sensitivity for low-abundance signals
  co<-cor(x)  #Pearson correlation calculation
  co[co[]<0]<-0  # Set negative correlations to zero
  co<-co
  # Association Strength Normalization
  b<-colSums(co)-1 #Calculate adjusted column sums (excluding self-correlations)
  co[co[]==1]<-NA
  for (i in 1:ncol(co)) { #Column-wise normalization，scale to [0,1] range maintaining relative proportions
    # Mathematical transformation:
    # RA_ij = (r_ij/Σr_ik) * min(1/RA_ij) 
    co[,i]<-co[,i]/b[i]   
    co[,i]<-co[,i]*min(1/co[,i],na.rm = T)
  }
  
  co<-as.matrix(co)^2 # Non-linear Enhancement
 # Visualization
  col=colorRampPalette(c("white","orange","firebrick3"))
  a<-corrplot(co,tl.cex = 0.5, col = col(200),na.label = "square",
              na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F,
              )           
  return(a)
}

RA_matrix<-discoPlot1(NC_matrix)

col=colorRampPalette(c("white","orange","firebrick3"))
pdf("/results/RA_bait_matrix.pdf", width = 8, height = 8)
corrplot(RA_matrix$corr,tl.cex = 0.8, col = col(200),na.label = "square",
          na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F,
        ) 
dev.off() 

#RA matrix of 18 NCs, calculated from above matrix
bait_matrix<-as.data.frame(RA_matrix$corr)
write.table(bait_matrix,file="/results/RA_bait_matrix.txt",sep="\t",quote=FALSE)
bait_matrix[is.na(bait_matrix)==T]<-1
NC_list<-list(TC=bait_matrix[1],
              NC_FC=bait_matrix[2],
              NC_DFC=bait_matrix[6:7],
              NC_GC=bait_matrix[3:5],
              CB=bait_matrix[8:11],
              GEM=bait_matrix[12:16],
              HLB=bait_matrix[17],
              CEN=bait_matrix[18],
              PNC=bait_matrix[19:20],
              PS=bait_matrix[21:24],
              CLB=bait_matrix[25],
              NS=bait_matrix[26:28],
              IK=bait_matrix[29],
              SNB=bait_matrix[30],
              PML_NB=bait_matrix[31:34],
              HET=bait_matrix[35:36],
              TEC=bait_matrix[37],
              NPC=bait_matrix[38],
              NL=bait_matrix[39:40],
              S68B=bait_matrix[41])

for (i in 1:length(NC_list)) {
  NC_list[[i]]$total<-rowSums(NC_list[[i]])/ncol(NC_list[[i]])
  names(NC_list[[i]])[ncol(NC_list[[i]])]<-names(NC_list)[i]
  NC_list[[i]]<-NC_list[[i]][ncol(NC_list[[i]])]
}

NC_crosstalk<-data.frame(list=names(NC_list))
rownames(NC_crosstalk)<-names(NC_list)

for (i in 1:length(NC_list)) {
  a<-c(NC_list[[i]][1,1],
       NC_list[[i]][2,1],
       sum(NC_list[[i]][6:7,1])/2,
       sum(NC_list[[i]][3:5,1])/3,
       sum(NC_list[[i]][8:11,1])/4,
       sum(NC_list[[i]][12:16,1])/5,
       NC_list[[i]][17,1],
       NC_list[[i]][18,1],
       sum(NC_list[[i]][19:20,1])/2,
       sum(NC_list[[i]][21:24,1])/4,
       NC_list[[i]][25,1],
       sum(NC_list[[i]][26:28,1])/3,
       NC_list[[i]][29,1],
       NC_list[[i]][30,1],
       sum(NC_list[[i]][31:34,1])/4,
       sum(NC_list[[i]][35:36,1])/2,
       NC_list[[i]][37,1], 
       NC_list[[i]][38,1],
       sum(NC_list[[i]][39:40,1])/2,
       NC_list[[i]][41,1])
  NC_crosstalk<-data.frame(NC_crosstalk,a)
  names(NC_crosstalk)[i+1]<-names(NC_list)[i]
}

normalize<-function(x){
  for (i in 1:ncol(x)) {
    x[,i]<-x[,i]*min(1/x[,i])
  }
  test<-x
  a<-1
  while(a<= ncol(x)){
    for (b in 1:nrow(x)) {
      test[b,a]<-mean(x[b,a]+x[a,b])
    }
    a=a+1
  }
  corrplot(test,tl.cex = 1.2, col = COL1('Purples'),na.label = "square",method = "square",
           na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F)
  
}  

n<-normalize(as.matrix(NC_crosstalk[-1]))
n<-as.matrix(n$corr)
write.table(n,file="/results/RA_NC_matrix.txt",sep="\t",quote=FALSE)


col=colorRampPalette(c("white","#7da9f0","#4e26b4"))
pdf("/results/RA_NC_matrix.pdf", width = 8, height = 8)
corrplot(n,tl.cex = 1, col = col(200),na.label = "square",method = "square",
         na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F, type = 'upper',diag = FALSE,
         col.lim = c(0,1))
dev.off()

#### Fig 3a,c prediction of protein localization and potential novel NC

UK_protein1<-read.csv("ESRP1.txt", sep = "\t",  header = TRUE)
UK_protein2<-read.csv("BUD13.txt", sep = "\t",  header = TRUE)

NC_matrix<-read.csv("NC_marker_matrix.txt", sep = "\t",  header = TRUE)

prediction_matrix<-merge(NC_matrix,UK_protein1,by=1,all = T)
prediction_matrix<-merge(prediction_matrix,UK_protein2,by=1,all = T)

names(prediction_matrix)[1]<-"Majority.protein.IDs"
prediction_matrix[is.na(prediction_matrix)]<-0


discoPlot2<-function(x){  ##x as a matrix
  x<-x[,-1]
  x<-x^0.5
  co<-cor(x)
  co[co[]<0]<-0
  co<-co[1:41,]
  b<-c(colSums(co[,1:41])-1,colSums(co[,42:ncol(co)]))
  co[co[]==1]<-NA
  for (i in 1:ncol(co)) {
    co[,i]<-co[,i]/b[i]
    co[,i]<-co[,i]*min(1/co[,i],na.rm = T)
  }
  
  co<-as.matrix(co)^2
  col=colorRampPalette(c("white","orange","firebrick3"))
  a<-corrplot(co,tl.cex = 0.6, col = col(200),na.label = "square",
              na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F)
  return(a)
}

UK_RA<-discoPlot2(prediction_matrix)
UK_RA_matrix<-UK_RA[["corr"]]

col=colorRampPalette(c("white","orange","firebrick3"))
pdf("/results/RA_unknown_protein_matrix.pdf", width = 5, height = 8)
corrplot(UK_RA_matrix[,42:43],tl.cex = 0.8, col = col(200),na.label = "square",
              na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F)
dev.off()

write.table(UK_RA_matrix[,42:43],file="/results/RA_unknown_protein_matrix.txt",sep="\t",quote=FALSE)

#calcualte the pearson's r
Prsn_r<-UK_RA_matrix
Prsn_r[is.na(Prsn_r)==T]<-1
Prsn_r<-cor(Prsn_r,method = "pearson")
Prsn_r<-Prsn_r[,42:ncol(Prsn_r)]
write.table(Prsn_r,file="/results/unknown_protein_pearson_R.txt",sep="\t",quote=FALSE)
#calcualte the p-value

p<-UK_RA_matrix
p[is.na(p)==T]<-1
p_matrix<-corr.test(p)
p_matrix<-p_matrix$p
p_matrix<-p_matrix[,42:ncol(p_matrix)]
write.table(p_matrix,file="/results/unknown_protein_p_value.txt",sep="\t",quote=FALSE)

#######NOVA map for Fig 4a, 4b, 5d
library(dplyr)
library(patchwork)
library(Hmisc)
library(ggplot2)
library(openxlsx)
library(harmony)
library(readxl)
library(corrplot)

discoPlotALL_t<-function(x){ ### Using same method to calculate a serises of RA scores of each bait
  x<-x[,-1]
  x<-x^0.5
  co<-cor(x)
  co[co[]<0]<-0
  co<-co[1:41,]
  b<-c(colSums(co[,1:41])-1,colSums(co[,42:ncol(co)]))
  co[co[]==1]<-NA
  for (i in 1:ncol(co)) {
    co[,i]<-co[,i]/b[i]
    co[,i]<-co[,i]*min(1/co[,i],na.rm = T)
  }
  
  co<-as.matrix(co)^2
  col=colorRampPalette(c("white","orange","firebrick3"))
  a<-corrplot(co,tl.cex = 0.4, col = col(200),na.label = "square",
              na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F)
  return(a)
}

#Loading datasets "362baits.txt", which includes PhastID datasets from this study and other published BioID datasets, to prepare the RA matrix for NOVA map clustering

NOVA_baits<-read.csv("362baits.txt", sep = "\t",  header = TRUE)
NOVA_RA_matrix<-discoPlotALL_t(NOVA_baits)
test<-NOVA_RA_matrix$corr
write.table(test,file="/results/NOVA_RA_matrix.txt",sep="\t",quote=FALSE,col.names = T,row.names =T)

####using UMAP to profile the NOVA map
library(Seurat)
test[is.na(test)==T]<-1.1 
bait_data<-as.matrix(test)  
bait <- CreateSeuratObject(counts = bait_data)
bait <- NormalizeData(bait, normalization.method = "LogNormalize", scale.factor = 10)
bait <- FindVariableFeatures(bait, selection.method = "vst", nfeatures = 41)
all.pros <- colnames(bait)
NC_vn <- ScaleData(bait, features = all.pros)
NC_vn <- RunPCA(NC_vn, features =  VariableFeatures(object = NC_vn))
print(NC_vn[["pca"]], dims = 1:5, nfeatures = 5)

##Clustering
NC_vn <- FindNeighbors(NC_vn, dims = 1:13)
NC_vn <- FindClusters(NC_vn, resolution = 1.5) 
NC_vn <- RunUMAP(NC_vn, dims = 1:7, label = T,seed.use =1)

# Extract UMAP coordinates
head(NC_vn@reductions$umap@cell.embeddings) 
p1 <- DimPlot(NC_vn, reduction = "umap",pt.size=1)
load("NOVA.rData")
mydata=NOVA@reductions$umap@cell.embeddings

## Visualizing NOVA map
##1. NOVA map
mydata=read.csv("UMAP-xy.txt", sep = "\t", header = TRUE)
color= read_excel("color label.xlsx",sheet = 1)
plot<-merge(mydata,data.frame(color[-1],row.names = color$node),by='row.names',all=F)

pdf("/results/NOVA.pdf", width = 8, height = 8)
NOVA<-ggplot(data = plot,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(size=plot$size,colour=as.factor(plot$color),alpha = 0.5)+
  theme_bw()+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())
print(NOVA)
dev.off()
# 2. NOVA of DAXX mutants
color= read_excel("color label.xlsx",sheet = 2)
plot<-merge(mydata,data.frame(color[-1],row.names = color$node),by='row.names',all=F)

highlight_points <- subset(plot, Row.names %in% c("DAXX","DAXX-A297P","DAXX-E465X","DAXX-L130R","DAXX-R306X"))

pdf("/results/NOVA-DAXX.pdf", width = 8, height = 8)
DAXX <- ggplot() +
  geom_point(data = plot,aes(x = UMAP_1, y = UMAP_2),
             size = plot$size,colour = as.factor(plot$color),alpha = 0.5) +
  geom_point(data = highlight_points,aes(x = UMAP_1, y = UMAP_2), color = "#CC0A0A",
             size = highlight_points$size * 1.2,alpha = 0.8) +
  ggrepel::geom_text_repel(data = highlight_points,
                           aes(x = UMAP_1, y = UMAP_2, label = Row.names),
                           color = "black", box.padding = 0.5, min.segment.length = 0, segment.color = "grey50",size = 4) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank())
print(DAXX)
dev.off()

###3. NOVA of HSF1 before and after heat shock
color= read_excel("color label.xlsx",sheet = 3)
plot<-merge(mydata,data.frame(color[-1],row.names = color$node), by='row.names',all=F)
highlight_HSF1 <- subset(plot, Row.names %in% c("HSF1-HS","HSF1+HS"))

pdf("/results/NOVA-HSF1.pdf", width = 8, height = 8)
HSF1 <- ggplot() +
  geom_point(data = plot,aes(x = UMAP_1, y = UMAP_2),
             size = plot$size,colour = as.factor(plot$color),alpha = 0.5) +
  geom_point(data = highlight_HSF1,aes(x = UMAP_1, y = UMAP_2), color = "#CC0A0A",
             size = highlight_HSF1$size * 1.2,alpha = 0.8) +
  ggrepel::geom_text_repel(data = highlight_HSF1,
                           aes(x = UMAP_1, y = UMAP_2, label = c("HSF1","HSF1+HS")),
                           color = "black", box.padding = 0.5, min.segment.length = 0, segment.color = "grey50",size = 4) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank())
print(HSF1)
dev.off()
