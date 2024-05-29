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
library(png)
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

show_plot <- function(plot, filename, width = 800, height = 400, res = 100) {
  png(filename = filename, width = width, height = height, res = res)
  print(plot)
  dev.off()
  img <- readPNG(filename)
  grid::grid.raster(img)
}

options(repr.plot.width = 20, repr.plot.height = 16)

data_dir <- R.utils::getAbsolutePath('../../data')
mouse_dir <- glue::glue("{data_dir}/Visium_Mouse_Olfactory_Bulb/outs")

raw_sfe <- readRDS(glue::glue("{data_dir}/Visium_Mouse_Olfactory_Bulb.rds"))
transposed_raw_sfe <- raw_sfe

is_mt <- str_detect(rowData(transposed_raw_sfe)$symbol, "^mt-")
sum(is_mt)
colData(transposed_raw_sfe)$nGenes <- colSums(counts(transposed_raw_sfe) > 0)
colData(transposed_raw_sfe)$nCounts <- colSums(counts(transposed_raw_sfe))
qc_sfe <- scuttle::addPerCellQCMetrics(transposed_raw_sfe, subsets = list(mito = is_mt))

processed_sfe <- transposed_raw_sfe[, qc_sfe$subsets_mito_percent < 20]
processed_sfe <- processed_sfe[rowSums(counts(processed_sfe)) > 0,]

colData(processed_sfe)$nCounts <- colSums(counts(processed_sfe))
