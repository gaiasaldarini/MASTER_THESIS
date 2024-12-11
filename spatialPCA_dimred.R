library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggsignif)
library(readr)
library(spacexr)
library(scatterpie)
library(pheatmap)
library(SpatialPCA)
library(glmGamPoi)
library(harmony)
#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')
library(presto)
library(igraph)
library(plotly)
library(ggpubr)

path = "/beegfs/scratch/ric.hsr/Gaia_saldarini/FINAL/"

samples_infected = c("C1", "D1", "C2", "D2", "B3", "C3", "D3")
col_samples = brewer.pal(7, "Reds")
names(col_samples) = samples_infected

tot_coor = readRDS(paste0(path, "data/tot_coor.rds"))
slices_merged = readRDS(paste0(path, "data/slices_merged.rds"))
host_genes = read.csv(paste0(path, "/data/host_genes.csv"))[[1]]
infected_aligned = slices_merged[host_genes, slices_merged$orig.ident %in% samples_infected]
beads = read.csv(paste0(path, "data/beads.csv"))[[1]]
filtered_spatialPCA = infected_aligned[, ! (colnames(infected_aligned) %in% c(beads))]
dim(filtered_spatialPCA)


# FATTO IN BACKGROUND CON 128G 8 CPU
z_values = c(0, 5, 35, 45, 65, 70, 80)

coor_2d = tot_coor[colnames(filtered_spatialPCA), c("tissue_x", "tissue_y")]
for ( i in 1:7){
  coor_2d$z[endsWith(rownames(coor_2d), samples_infected[i])] = z_values[i]
}

combined_coords = coor_2d
combined_coords$z = combined_coords$z / 10 # per avere x y z nella stessa unità di misura
colnames(combined_coords) = c("x", "y", "z")
dim(combined_coords)

count_sub = filtered_spatialPCA@assays$Spatial$counts
xy_coords = as.matrix(combined_coords[, 1:3])
dim(count_sub)
dim(xy_coords)

LIBD = CreateSpatialPCAObject(counts=count_sub, location=xy_coords, project = "SpatialPCA",gene.type="spatial",sparkversion="sparkx",numCores_spark=8,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20)
LIBD = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
LIBD = SpatialPCA_EstimateLoading(LIBD,fast=FALSE,SpatialPCnum=20) 
LIBD = SpatialPCA_SpatialPCs(LIBD, fast=FALSE)

saveRDS(LIBD, paste0(path, "data/LIBD.rds"))

dim(LIBD@counts)
dim(LIBD@normalized_expr)
spatially_variable_genes = rownames(LIBD@normalized_expr)
saveRDS(spatially_variable_genes, paste0(path, "data/spatially_variable_genes.rds"))
spatialPCs = LIBD@SpatialPCs 
saveRDS(spatialPCs, paste0(path, "data/spatialPCs.rds"))