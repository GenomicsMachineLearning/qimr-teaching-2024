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
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.box.just = "right"
    )
}
options(repr.plot.width = 20, repr.plot.height = 16)

data_dir <- R.utils::getAbsolutePath('../../data')
mouse_dir <- glue::glue("{data_dir}/Visium_Mouse_Olfactory_Bulb/outs")

raw_sfe <- readRDS(glue::glue("{data_dir}/Visium_Mouse_Olfactory_Bulb.rds"))
transposed_raw_sfe <- raw_sfe

