
library(Seurat)
library(spacexr)

path = "/beegfs/scratch/ric.hsr/Gaia_saldarini/FINAL/"

## METHOD DECONVOLUTION - RCTD

reference = readRDS("/beegfs/scratch/ric.cirillo.transcriptomics/ric.cirillo.transcriptomics/Progetti_Lore/Xu_ref.rds")
cell_types = as.factor(reference$cell_type)
names(cell_types) = colnames(reference)
reference = Reference(counts = reference@assays$RNA@counts, cell_types = cell_types)

slices_list = readRDS(paste0(path, "data/slices_list_unfiltered.rds"))

n_slices = 10
for (i in 1:n_slices){
  
  print(i)
  slices_list[[i]] = SCTransform(slices_list[[i]],assay='Spatial',vst.flavor='v2')
  puck<-spacexr::SpatialRNA(slices_list[[i]]@images[[1]]@coordinates[,c('imagerow','imagecol')], slices_list[[i]]@assays$SCT@counts)
  
  
  spacexr::create.RCTD(puck,reference,max_cores = 2,test_mode = F, CELL_MIN_INSTANCE = 6)->myRCTD
  
  
  myRCTD <- spacexr::run.RCTD(myRCTD, doublet_mode = 'full')
  
  
  
  results <- myRCTD@results
  
  norm_weights = as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))
  norm_weights <- t(norm_weights)
  ncol <- length(colnames(slices_list[[i]])[!(colnames(slices_list[[i]]) %in% colnames(norm_weights))])
  if (ncol > 0) {
    tmp <- matrix(0, nrow = nrow(norm_weights), ncol = ncol)
    colnames(tmp) <- colnames(slices_list[[i]])[!(colnames(slices_list[[i]]) %in% colnames(norm_weights))]
    norm_weights <- cbind(norm_weights, tmp)[,colnames(slices_list[[i]])]
  }
  
  slices_list[[i]][['RCTD']]<-CreateAssayObject(data = norm_weights)
  
  
}

saveRDS(slices_list, paste0(path, "data/deconvolution_list.rds"))