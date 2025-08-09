rm(list=ls())
options(stringsAsFactors = F)
.libPaths(c("~/SeuratV4", .libPaths()))
library(dior)
library(Seurat)
library(ggplot2)
library(copykat)
library(infercnv)
library(AnnoProbe)
setwd("/Volumes/WD/ccRCC_Epi/Epithelial/Epi_copykat/copykat_result/")

load("/Volumes/WD/ccRCC_Epi/Epithelial/EpithelialRaw.RData")
table(experiment.aggregate@meta.data[["orig.ident"]],experiment.aggregate@meta.data[["dataset"]])

experiment.aggregate<-subset(experiment.aggregate,dataset=="GSE156632")
experiment.aggregate@meta.data[["dataset"]]<-as.factor(as.character(experiment.aggregate@meta.data[["dataset"]]))
experiment.aggregate@meta.data[["group"]]<-as.factor(as.character(experiment.aggregate@meta.data[["group"]]))
experiment.aggregate@meta.data[["orig.ident"]]<-as.factor(as.character(experiment.aggregate@meta.data[["orig.ident"]]))
for (pi in c("1","2","3","4","5")) {
  scT<-subset(experiment.aggregate,orig.ident==paste("GSE156632T",pi,sep = ""))
  scN<-subset(experiment.aggregate,orig.ident==paste("GSE156632N",pi,sep = ""))
  sc2<-merge(scT,scN)
  normal_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("GSE156632N",pi,sep = ""))]
  normalMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=normal_cells),assay = "RNA",slot = "counts"))
  tumor_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("GSE156632T",pi,sep = ""))]
  tumorMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=tumor_cells),assay = "RNA",slot = "counts"))
  dat=cbind(normalMat,tumorMat)
  results <- copykat(rawmat=dat, id.type="S",cell.line = "no",ngene.chr = 1,
                     win.size=25, KS.cut=0.1, sam.name=paste(unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),sep = ""),
                     distance="euclidean", norm.cell.names=normal_cells,plot.genes = "TRUE",genome = "hg20",
                     n.cores=10,output.seg="FLASE")
  save(results,file = paste("CopyKAT_",unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),".RData",sep = ""))
}

load("/Volumes/WD/ccRCC_Epi/Epithelial/EpithelialRaw.RData")
experiment.aggregate<-subset(experiment.aggregate,dataset=="GSE159115")
experiment.aggregate@meta.data[["dataset"]]<-as.factor(as.character(experiment.aggregate@meta.data[["dataset"]]))
experiment.aggregate@meta.data[["group"]]<-as.factor(as.character(experiment.aggregate@meta.data[["group"]]))
experiment.aggregate@meta.data[["orig.ident"]]<-as.factor(as.character(experiment.aggregate@meta.data[["orig.ident"]]))
for (pi in c("2","3","5","6")) {
  scT<-subset(experiment.aggregate,orig.ident==paste("GSE159115P",pi,"T",sep = ""))
  scN<-subset(experiment.aggregate,orig.ident==paste("GSE159115P",pi,"N",sep = ""))
  sc2<-merge(scT,scN)
  normal_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("GSE159115P",pi,"N",sep = ""))]
  normalMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=normal_cells),assay = "RNA",slot = "counts"))
  tumor_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("GSE159115P",pi,"T",sep = ""))]
  tumorMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=tumor_cells),assay = "RNA",slot = "counts"))
  dat=cbind(normalMat,tumorMat)
  results <- copykat(rawmat=dat, id.type="S",cell.line = "no",ngene.chr = 1,
                     win.size=25, KS.cut=0.1, sam.name=paste(unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),sep = ""),
                     distance="euclidean", norm.cell.names=normal_cells,plot.genes = "TRUE",genome = "hg20",
                     n.cores=10,output.seg="FLASE")
  save(results,file = paste("CopyKAT_",unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),".RData",sep = ""))
}

load("/Volumes/WD/ccRCC_Epi/Epithelial/EpithelialRaw.RData")
experiment.aggregate<-subset(experiment.aggregate,dataset=="GSE178481")
experiment.aggregate@meta.data[["dataset"]]<-as.factor(as.character(experiment.aggregate@meta.data[["dataset"]]))
experiment.aggregate@meta.data[["group"]]<-as.factor(as.character(experiment.aggregate@meta.data[["group"]]))
experiment.aggregate@meta.data[["orig.ident"]]<-as.factor(as.character(experiment.aggregate@meta.data[["orig.ident"]]))
for (pi in c("BM1","PR2","PR4","PR6")) {
  scT<-subset(experiment.aggregate,orig.ident==paste("RCC.",pi,".PTumor",sep = ""))
  scN<-subset(experiment.aggregate,orig.ident==paste("RCC.",pi,".Normal",sep = ""))
  sc2<-merge(scT,scN)
  normal_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("RCC.",pi,".Normal",sep = ""))]
  normalMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=normal_cells),assay = "RNA",slot = "counts"))
  tumor_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("RCC.",pi,".PTumor",sep = ""))]
  tumorMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=tumor_cells),assay = "RNA",slot = "counts"))
  dat=cbind(normalMat,tumorMat)
  results <- copykat(rawmat=dat, id.type="S",cell.line = "no",ngene.chr = 1,
                     win.size=25, KS.cut=0.1, sam.name=paste(unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),sep = ""),
                     distance="euclidean", norm.cell.names=normal_cells,plot.genes = "TRUE",genome = "hg20",
                     n.cores=10,output.seg="FLASE")
  save(results,file = paste("CopyKAT_",unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),".RData",sep = ""))
}

load("/Volumes/WD/ccRCC_Epi/Epithelial/EpithelialRaw.RData")
experiment.aggregate<-subset(experiment.aggregate,dataset=="GSE178481")
experiment.aggregate@meta.data[["dataset"]]<-as.factor(as.character(experiment.aggregate@meta.data[["dataset"]]))
experiment.aggregate@meta.data[["group"]]<-as.factor(as.character(experiment.aggregate@meta.data[["group"]]))
experiment.aggregate@meta.data[["orig.ident"]]<-as.factor(as.character(experiment.aggregate@meta.data[["orig.ident"]]))
for (pi in c("PR3.PTumor1","PR3.PTumor2","PR3.PTumor3")) {
  scT<-subset(experiment.aggregate,orig.ident==paste("RCC.",pi,sep = ""))
  scN<-subset(experiment.aggregate,orig.ident==paste("RCC.","PR3.Normal",sep = ""))
  sc2<-merge(scT,scN)
  normal_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("RCC.","PR3.Normal",sep = ""))]
  normalMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=normal_cells),assay = "RNA",slot = "counts"))
  tumor_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("RCC.",pi,sep = ""))]
  tumorMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=tumor_cells),assay = "RNA",slot = "counts"))
  dat=cbind(normalMat,tumorMat)
  results <- copykat(rawmat=dat, id.type="S",cell.line = "no",ngene.chr = 1,
                     win.size=25, KS.cut=0.1, sam.name=paste(unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),sep = ""),
                     distance="euclidean", norm.cell.names=normal_cells,plot.genes = "TRUE",genome = "hg20",
                     n.cores=10,output.seg="FLASE")
  save(results,file = paste("CopyKAT_",unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),".RData",sep = ""))
}

load("/Volumes/WD/ccRCC_Epi/Epithelial/EpithelialRaw.RData")
experiment.aggregate<-subset(experiment.aggregate,dataset=="GSE178481")
experiment.aggregate@meta.data[["dataset"]]<-as.factor(as.character(experiment.aggregate@meta.data[["dataset"]]))
experiment.aggregate@meta.data[["group"]]<-as.factor(as.character(experiment.aggregate@meta.data[["group"]]))
experiment.aggregate@meta.data[["orig.ident"]]<-as.factor(as.character(experiment.aggregate@meta.data[["orig.ident"]]))
for (pi in c("PR5.PTumor1","PR5.PTumor2","PR5.PTumor3")) {
  scT<-subset(experiment.aggregate,orig.ident==paste("RCC.",pi,sep = ""))
  scN<-subset(experiment.aggregate,orig.ident==paste("RCC.","PR5.Normal",sep = ""))
  sc2<-merge(scT,scN)
  normal_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("RCC.","PR5.Normal",sep = ""))]
  normalMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=normal_cells),assay = "RNA",slot = "counts"))
  tumor_cells  <- row.names(sc2@meta.data)[which(sc2@meta.data$orig.ident==paste("RCC.",pi,sep = ""))]
  tumorMat=as.data.frame(GetAssayData(subset(experiment.aggregate, cells=tumor_cells),assay = "RNA",slot = "counts"))
  dat=cbind(normalMat,tumorMat)
  results <- copykat(rawmat=dat, id.type="S",cell.line = "no",ngene.chr = 1,
                     win.size=25, KS.cut=0.1, sam.name=paste(unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),sep = ""),
                     distance="euclidean", norm.cell.names=normal_cells,plot.genes = "TRUE",genome = "hg20",
                     n.cores=10,output.seg="FLASE")
  save(results,file = paste("CopyKAT_",unique(scT@meta.data[["dataset"]]),unique(scT@meta.data[["orig.ident"]]),".RData",sep = ""))
}