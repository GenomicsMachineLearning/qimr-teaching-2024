library(dplyr)
library(Voyager)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(ggplot2)
library(scater)
library(rlang)
library(scran)
library(scuttle)
library(terra)
library(sf)
library(rmapshaper)
library(scran)
library(stringr)
library(EBImage)
library(patchwork)
library(bluster)
library(rjson)
theme_set(theme_bw())

# Layout
custom_theme <- function() {
  theme_bw() +
    theme(
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10, face = "bold"),
      legend.position = "right",
      legend.box.just = "right"
    )
}

## Override functions for notebook uste - ignore for command line
# custom_theme <- function() {
#   theme_bw() +
#     theme(
#       legend.text = element_text(size = 14),
#       legend.title = element_text(size = 16, face = "bold"),
#       axis.text = element_text(size = 12),
#       axis.title = element_text(size = 14, face = "bold"),
#       legend.position = "right",
#       legend.box.just = "right"
#     )
# }

# Display plot over X-Windows
# show_plot <- function(plot, width = 800, height = 400, res = 100) {
#   temp_filename <- tempfile(fileext = ".png")
#   while (dev.cur() > 1) {
#     dev.off()
#   }
#   grDevices::png(filename = temp_filename, width = width, height = height, res = res)
#   print(plot)
#   dev.off()
#   img <- png::readPNG(temp_filename)
#   grid::grid.raster(img)
# }
# Display plot normally
show_plot <- function(plot, width = 800, height = 400, res = 100) {
    return(plot)
}

# Setup plotting resolution.
options(repr.plot.width = 20, repr.plot.height = 16)

#Original way to load Visium Data
#raw_sfe <- SpatialFeatureExperiment::read10xVisiumSFE(dirs = skin_data_dir, samples = ".", type = "sparse", data = "raw")
#Voyager::plotImage(raw_sfe)
#transposed_raw_sfe <- SpatialFeatureExperiment::transpose(raw_sfe)

# Setup Data Directory
data_dir <- R.utils::getAbsolutePath('../data')
mouse_dir <- glue::glue("{data_dir}/Visium_Mouse_Olfactory_Bulb/outs")
skin_A2_dir <- glue::glue("{data_dir}/Visium_Skin_A2/outs")
raw_sfe <- readRDS(glue::glue("{data_dir}/Visium_Mouse_Olfactory_Bulb.rds"))
transposed_raw_sfe <- raw_sfe

# Set Random Number Seed.
set.seed(29)

# Get number of genes detected and counts
colData(transposed_raw_sfe)$nGenes <- colSums(counts(transposed_raw_sfe) > 0)
colData(transposed_raw_sfe)$nCounts <- colSums(counts(transposed_raw_sfe))

# Find mitochdonrial genes
is_mt <- str_detect(rowData(transposed_raw_sfe)$symbol, "^MT-")
qc_sfe <- scuttle::addPerCellQCMetrics(transposed_raw_sfe, subsets = list(mito = is_mt))

# Normally MITO % is set to 20.
processed_sfe <- transposed_raw_sfe[, qc_sfe$subsets_mito_percent < 20]

# Remove genes with zero counts
processed_sfe <- processed_sfe[rowSums(counts(processed_sfe)) > 0,]
colData(processed_sfe)$nCounts <- colSums(counts(processed_sfe))
colData(processed_sfe)$nGenes <- colSums(counts(processed_sfe) > 0)
processed_sfe <- scuttle::logNormCounts(processed_sfe)

# Create an in tissue sample
sfe_tissue <- processed_sfe[, colData(processed_sfe)$in_tissue]
sfe_tissue <- sfe_tissue[rowSums(counts(sfe_tissue)) > 0,]
rowData(sfe_tissue)$means <- rowMeans(counts(sfe_tissue))
rowData(sfe_tissue)$vars <- rowVars(counts(sfe_tissue))

# Coefficient of variance
rowData(sfe_tissue)$cv2 <- rowData(sfe_tissue)$vars/rowData(sfe_tissue)$means^2
sfe_tissue <- scuttle::logNormCounts(sfe_tissue)

# Run PCAs
processed_dec <- scran::modelGeneVar(processed_sfe, lowess = FALSE)
processed_hvgs <- scran::getTopHVGs(processed_dec, n = 2000)
processed_sfe <- BiocSingular::runPCA(processed_sfe, ncomponents = 30, subset_row = processed_hvgs, scale = TRUE)
tissue_dec <- scran::modelGeneVar(sfe_tissue, lowess = FALSE)
tissue_hvgs <- scran::getTopHVGs(tissue_dec, n = 2000)
sfe_tissue <- BiocSingular::runPCA(sfe_tissue, ncomponents = 30, subset_row = tissue_hvgs, scale = TRUE)

# Cluster
colData(processed_sfe)$cluster <- bluster::clusterRows(
    reducedDim(processed_sfe, "PCA")[,1:3],
    BLUSPARAM = SNNGraphParam(cluster.fun = "leiden",
    cluster.args = list(resolution_parameter = 0.5, objective_function = "modularity"))
)
colData(sfe_tissue)$cluster <- bluster::clusterRows(
    reducedDim(sfe_tissue, "PCA")[,1:3],
    BLUSPARAM = SNNGraphParam(cluster.fun = "leiden",
    cluster.args = list(resolution_parameter = 0.5, objective_function = "modularity"))
)

# Run UMAP
processed_sfe <- scater::runUMAP(processed_sfe, dimred = "PCA", n_dimred = 3)
sfe_tissue <- scater::runUMAP(sfe_tissue, dimred = "PCA", n_dimred = 3)

# Find Markers for in Tissue
markers <- scran::findMarkers(sfe_tissue, groups = colData(sfe_tissue)$cluster,
                       test.type = "wilcox", pval.type = "all", direction = "up")
genes_use <- vapply(markers, function(x) rownames(x)[1], FUN.VALUE = character(1))

# Annotate Tissue with Polygon
sp <- SpatialFeatureExperiment::spotPoly(sfe_tissue)
SpatialFeatureExperiment::dimGeometry(sfe_tissue, "spotPoly", MARGIN = 2) <- sp
(tb <- SpatialFeatureExperiment::annotGeometry(sfe_tissue, "tissueBoundary"))

# Run all for spatial graphics.
(g_all <- SpatialFeatureExperiment::findSpatialNeighbors(processed_sfe, MARGIN = 2, method = "tri2nb"))
(g_specific <- SpatialFeatureExperiment::findSpatialNeighbors(sfe_tissue, MARGIN = 2, method = "tri2nb"))
spatialGraph(processed_sfe, "graph1", MARGIN = 2) <- g_all
spatialGraph(sfe_tissue, "graph1", MARGIN = 2) <- g_specific
colGraph(sfe_tissue, "visium") <- SpatialFeatureExperiment::findVisiumGraph(sfe_tissue, zero.policy = TRUE)
colGraph(sfe_tissue, "visium_B") <- SpatialFeatureExperiment::findVisiumGraph(sfe_tissue, style = "B", zero.policy = TRUE)
sfe_tissue <- Voyager::runUnivariate(sfe_tissue, type = "localG_perm", features = "Ptgds", colGraphName = "visium_B",
                                     swap_rownames = "symbol")
dec <- scran::modelGeneVar(sfe_tissue)
hvgs <- scran::getTopHVGs(dec, n = 2000)
sfe_tissue <- Voyager::colDataUnivariate(sfe_tissue, features = c("nCounts", "nGenes"), colGraphName = "graph1", type = "moran")
sfe_tissue <- Voyager::runMoransI(sfe_tissue, features = hvgs)

neg_moran <- rownames(sfe_tissue)[order(rowData(sfe_tissue)$moran_sample01, decreasing = FALSE)[1:6]]
sfe_tissue <- Voyager::runUnivariate(sfe_tissue, "moran.mc", neg_moran, colGraphName = "graph1", nsim = 200, alternative = "less")
colGraph(sfe_tissue, "knn5") <- SpatialFeatureExperiment::findSpatialNeighbors(sfe_tissue, method = "knearneigh", dist_type = "idw", k = 5, style = "W")
sfe_tissue <- Voyager::runMultivariate(sfe_tissue, "multispati", colGraphName = "knn5", nfposi = 20, nfnega = 20)
sfe_tissue$clusts_nonspatial <- scran::clusterCells(sfe_tissue, use.dimred = "PCA",
                                             BLUSPARAM = NNGraphParam(
                                              cluster.fun = "leiden",
                                              cluster.args = list(
                                                  objective_function = "modularity",
                                                  resolution_parameter = 1
                                              )
                                          ))
sfe_tissue$clusts_multispati <- bluster::clusterRows(SingleCellExperiment::reducedDim(sfe_tissue, "multispati")[,1:20],
                                            BLUSPARAM = NNGraphParam(
                                                cluster.fun = "leiden",
                                                cluster.args = list(
                                                    objective_function = "modularity",
                                                    resolution_parameter = 1
                                                )
                                            )
                                           )
