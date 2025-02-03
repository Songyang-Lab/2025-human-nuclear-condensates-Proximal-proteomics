---
# title: "Proximal proteomics reveals a landscape of human nuclear condensates"
# author: "Ruofei Li"
# date: "03/02/2025"
---
  
# Load Required libraries

library(pheatmap)  
library(corrplot)
library(psych)

#### Fig 1d GO signature analysis


#The file "GO_AllLists.csv" contains GO analysis results downloaded from Metascape (https://metascape.org/)

setwd("...")

data<-read.csv("GO_AllLists.csv", sep = ",",  header = TRUE) 
data<-data[,c(3,4,1,17)] 
data<-data[which(data$Log.q.value.<= -2),]

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

## Set the name of the protein
names(a_map)[4]<-"my protein"

bk<-c(seq(0,19.9,by=0.1),seq(20,30,by=0.1))

p_col<-c(colorRampPalette(c("#0b68ff","#a0abfd"))(3),
         colorRampPalette(c("#6e50ec","#cf9fff"))(4),
         colorRampPalette(c("#c53207","#ffc16d"))(4),
         colorRampPalette(c("#f9f0a7","#9dd085","#42b4af"))(4)
)
names(p_col)<-GO_category_order$CategoryID
ann_col<-list(CategoryID=p_col)


pheatmap(a_map[4],cluster_rows = F,cluster_cols = F,
         color = c(colorRampPalette(c("white","skyblue1","slateblue2","purple4"))(length(bk)*2/3),colorRampPalette(c("purple4","firebrick3"))(length(bk)*1/3)),
         show_rownames = F,annotation_names_row =F,
         annotation_row = a_map[3],breaks = bk
         ,annotation_colors = ann_col
)


#### Fig 2a,c and Extended Fig 3c, RA matrix

NC_matrix<-read.csv("NC_marker_matrix.txt", sep = "\t",  header = TRUE)

#RA matrix of 41 NC marker proteins
discoPlot1<-function(x){
  x<-x[,-1]
  x<-x^0.5
  co<-cor(x)
  co[co[]<0]<-0
  co<-co
  b<-colSums(co)-1
  co[co[]==1]<-NA
  for (i in 1:ncol(co)) {
    co[,i]<-co[,i]/b[i]
    co[,i]<-co[,i]*min(1/co[,i],na.rm = T)
  }
  
  co<-as.matrix(co)^2
  col=colorRampPalette(c("white","orange","firebrick3"))
  a<-corrplot(co,tl.cex = 0.5, col = col(200),na.label = "square",
              na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F)
  return(a)
}

RA_matrix<-discoPlot1(NC_matrix)


#RA matrix of 18 NCs

bait_matrix<-as.data.frame(RA_matrix$corr)

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
       NC_matrix[[i]][37,1], 
       NC_matrix[[i]][38,1],
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

col=colorRampPalette(c("white","#7da9f0","#4e26b4"))
corrplot(n,tl.cex = 0.6, col = col(200),na.label = "square",method = "square",
         na.label.col="grey",tl.col = "black",outline =F,addrect  = 0.5,is.corr = F, type = 'upper',diag = FALSE,
         col.lim = c(0,1))


#### Fig 3a,c prediction of protein localization and potential novel NC

UK_protein1<-read.csv("ESRP1.txt", sep = "\t",  header = TRUE)
UK_protein2<-read.csv("BUD13.txt", sep = "\t",  header = TRUE)

NC_matrix<-read.csv("NC_marker_matrix.txt", sep = "\t",  header = TRUE)

prediction_matrix<-merge(NC_matrix,UK_protein1,by=1,all = T)
prediction_matrix<-merge(prediction_matrix,UK_protein2,by=1,all = T)

names(prediction_matrix)[1]<-"Majority.protein.IDs"
prediction_matrix[is.na(prediction_matrix)]<-0


discoPlot2<-function(x){  ##x为矩阵
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

#calcualte the pearson's r
Prsn_r<-UK_RA_matrix
Prsn_r[is.na(Prsn_r)==T]<-1
Prsn_r<-cor(Prsn_r,method = "pearson")
Prsn_r<-Prsn_r[,42:ncol(Prsn_r)]

#calcualte the p-value

p<-UK_RA_matrix
p[is.na(p)==T]<-1
p_matrix<-corr.test(p)
p_matrix<-p_matrix$p
p_matrix<-p_matrix[,42:ncol(p_matrix)]



