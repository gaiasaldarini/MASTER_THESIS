# MASTER_THESIS
## Master degree in Bioinformatics for computational genomics, Università degli studi di Milano, Politecnico di Milano
## Internship at COSR - IRCCS San Raffaele

Rmd scripts for analyses of spatial transcriptomic data (Visium):
1. Quality.Rmd: data import, quality filters (spots/genes/mt), beads identification, MABS info as covariates in Seurat object, RCTD deconvolution pipeline
2. Ctrl_Infected.Rmd: SC-like pipeline on 10 slices, comparison between 2 biological conditions, focus on shared clusters
3. Deconvolution.Rmd: deconvolution results on 12 cell types (3 macro-groups), correlation for co-localization, BALT structures (T and B cells interaction, Ccr6), macrophages signatures (M1, M2)
4. 3D_workflow.Rmd: 3D pipeline on Infected slices (slices alignment, SCT, SVG, spatialPCA, Harmony), focus on reproducibility (cfr cluster 6, 9), merge to define macro-domains, characterization of macro-domains, 3D visualization, MABS distribution and INOS pathway along macro-domains
5. SC_vs_3D_workflow.Rmd: SC_like pipeline on 7 Infected, comparison with 3D workflow results at the level of variable features, clusters markers, clusters spatial distribution (Moran index)
6. Granuloma.Rmd: granuloma identification, SC-like pipeline, biological characterization of sub-clusters, DESeq model on granuloma size on z-axis, trajectory analysis-pseudotime, 3D visualization, MABS distribution
7. Bulk_MABS.Rmd: public bulk transcriptomic data, analysis with DESeq, WGCNA, focus on 2 bacterial transcripts (rpoB, lsr2)
8. Pathogen.Rmd: spatial distribution, DEG MABS+ MABS-, 3D graph on MABS+ spots, Leiden clustering of the graph, DeSeq analysis on rpoB/lsr2 signal on the graph's clusters, study of cell types distribution along diffusion pattern from the beads to MABS+ spots 
  
