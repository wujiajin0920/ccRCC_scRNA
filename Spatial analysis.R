library(SpaCET)
library(Seurat)
library(RColorBrewer)
library(ggplot2)
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
# filter low-quality spots and calculate the QC metrics
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=10)
# plot the QC metrics
SpaCET.visualize.spatialFeature(
  SpaCET_obj,colors = brewer.pal(5, 'Oranges'),
  spatialType = "QualityControl", 
  spatialFeatures=c("UMI"),
  imageBg = TRUE
)
ggsave("QC_UMI.pdf",width = 5.2,height = 5)
SpaCET.visualize.spatialFeature(
  SpaCET_obj, colors = brewer.pal(5, 'Reds'),
  spatialType = "QualityControl", 
  spatialFeatures=c("Gene"),
  imageBg = TRUE
)
ggsave("QC_Gene.pdf",width = 5.2,height = 5)
# deconvolve ST data
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType="KIRC", coreNo=8)
# show the ST deconvolution results
SpaCET_obj@results$deconvolution$propMat[,1:6]
# show the spatial distribution of malignant cells and macrophages.
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant")
)
ggsave("Malignant.pdf",width = 5.2,height = 5)

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Macrophage M2")
)
ggsave("Macrophage_M2.pdf",width = 5.2,height = 5)

SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="Hallmark")
rownames(SpaCET_obj@results$GeneSetScore)
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
)
ggsave("Hallmark_EMT.pdf",width = 5.2,height = 5)

gmt1 <- list(
  Plasitic = c("CA9","NPTX2","CD70","CCDC146","AXL","CAV1","IFI44L","MT3","PHLDA1","RNASET2","BIRC3","VCAN","RGS5"))
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = gmt1)
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures=c("Plasitic")
)
ggsave("Plasitic_GeneSetScore.pdf",width = 5.2,height = 5)

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneExpression", 
  spatialFeatures=c("AXL")
)
ggsave("AXL_GeneExpression.pdf",width = 5.2,height = 5)

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant","Macrophage","CAF","Neutrophil")
)
# show the spatial distribution of all cell types.
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures="All", 
  sameScaleForFraction = TRUE,
  pointSize = 0.1, 
  nrow=5
)
# calculate the cell-cell colocalization.
SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)

# visualize the cell-cell colocalization.
SpaCET.visualize.colocalization(SpaCET_obj)
ggsave("GSM7974884_colocalization.pdf",width = 16,height = 7)
# calculate the L-R network score across ST spots.
SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj,coreNo=8)

# visualize the L-R network score.
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "LRNetworkScore", 
  spatialFeatures=c("Network_Score","Network_Score_pv")
)
# Ligand-Receptor analysis for a co-localized cell-type pair
SpaCET_obj <- SpaCET.CCI.cellTypePair(SpaCET_obj, cellTypePair=c("CAF","Macrophage"))
## [1] "CAF and Macrophage M2 have potential intercellular interaction in the current tissue."

# Visualize the interaction analysis of a co-localized cell-type pair.
SpaCET.visualize.cellTypePair(SpaCET_obj, cellTypePair=c("CAF","Macrophage"))
# Identify the Tumor-Stroma Interface
SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)

# Visualize the Interface
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface")
# Combine the interface and interaction spots
SpaCET_obj <- SpaCET.combine.interface(SpaCET_obj, cellTypePair=c("CAF","Macrophage M2"))

# Visualize the Interface. The spatialFeature should be formated as "Interface&celltype1_celltype2". Celltype 1 and 2 is in the order of alphabet.
SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", spatialFeature = "Interface&CAF_Macrophage M2")
# Compute the distance of CAF-M2 to tumor border
SpaCET.distance.to.interface(SpaCET_obj, cellTypePair=c("CAF", "Macrophage M2"))
# further deconvolve malignant cell states
SpaCET_obj <- SpaCET.deconvolution.malignant(SpaCET_obj, coreNo = 8)

# show cancer cell state fraction of the first five spots
SpaCET_obj@results$deconvolution$propMat[c("Malignant cell state A","Malignant cell state B"),1:6]

SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant","Malignant cell state A","Malignant cell state B"), 
  nrow=1
)

# run gene set calculation
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets="TLS")

# visualize TLS
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "GeneSetScore", 
  spatialFeatures = c("TLS")
)
