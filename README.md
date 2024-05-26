# Materials

* 001 - Clustering/Cell Typing/Label Transfer/Deconvolution (Andrew Causer)
  * Voyager - Cell Segmentation, Visualising QC, Clustering, Moran-I, etc.
* 002 - Spatial Statistics (Andrew Newman)
  * Voyager - Cell Segmentation, Visualising QC, Moran-I, etc.
* 003 - Single cell visualization (Levi & Xiao)
  * Naspari for CosMX (Levi)
  * stLearn/Xenium Explorer for Xenium/CODEX (Xiao)
* 004 - Multimodalities 
  * Inferred Copy Number Variation (Prakrithi)
  * Spatial metabolomics + transcriptomics (Andrew)
  * CODEX + Visium (Xiao)
* 005 - Community Analysis [Feng]
  *
* 006 - Spatial Proteomics [Xiao & Quan]
  * 
* 007 - cell-cell interaction analysis [Onkar - Levi]
  * 
* 008 - Deep Learning [Xiao & Quan]
  * 

# Data

* One Visium, one Xenium, one CosmiX
  * Melanoma - single cell, Visium, Xenium, and CosmiX.
* Instructions on how to download

# Installation

Using Conda/Mamba/Micromamba:
* MacOS
* Windows

How many environments:
* One for Voyager/Clustering/etc in R.
* 

## Setup Micromamba
```
$ /working/joint_projects/P3903/teaching2024-winter-qimr/micromamba/micromamba shell init
```

## Mamba Configuration

Create ~/.condarc

```
channels:
  - conda-forge
  - bioconda
  - defaults
channel_priority: flexible
```

## Install Dependencies
This installs:
 * R and dependencies for Seurat v5, Voyager and CARD, and
 * Python and dependencies for Jupyter.

MacOS:

```
$ conda create --name qimr-teaching --subdir osx-64 python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching
```

Others:
```
$ conda create --name qimr-teaching python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching
```

HPC:
```
$ micromamba create -p /working/joint_projects/P3903/teaching2024-winter-qimr/conda-envs/base-python-R python=3.10 r-base=4.3 r-devtools -y
$ micromamba activate /mnt/lustre/working/joint_projects/P3903/teaching2024-winter-qimr/conda-envs/base-python-R
```

For OSes other than HPC replace all calls to micromamba to conda.

Dependencies:
```
$ micromamba install -c bioconda bioconductor-scater bioconductor-splatter bioconductor-edgeR bioconductor-bluster bioconductor-BiocFileCache bioconductor-glmGamPoi bioconductor-SingleCellExperiment bioconductor-SummarizedExperiment bioconductor-ScaledMatrix bioconductor-BiocParallel bioconductor-ebimage r-anndata -y
$ micromamba install -c conda-forge r-scico r-ggnewscale r-magick r-rjson r-ragg r-units r-stringi r-sf r-s2 r-reticulate r-stringi r-tidyverse r-r.utils r-Seurat r-SeuratObject r-sctransform r-proj r-rcpptoml r-spdep r-lme4 r-ggrastr r-dbscan r-hdf5r r-optparse r-memuse r-sfheaders r-zeallot r-rmapshaper -y
$ micromamba install -c bioconda presto r-presto bioconductor-dropletutils r-MuSiC -y 
$ micromamba install -c conda-forge r-mcmcpack r-fields r-concaveman r-scatterpie r-ggcorrplot r-nnls r-pbmcapply r-NMF -y
```

R Dependencies to Install Directly:
```
$ R
> install.packages("remotes", dependencies = FALSE)
> remotes::install_github("drighelli/SpatialExperiment", dependencies=FALSE)
> remotes::install_github("pachterlab/SpatialFeatureExperiment", ref="devel", dependencies=FALSE)
> remotes::install_github("pachterlab/Voyager", ref="devel", dependencies=FALSE)
> remotes::install_version('wrMisc', dependencies = FALSE)
> remotes::install_github('YingMa0107/CARD', dependencies = FALSE)
```

Python Dependencies:
```
$ micromamba install -c conda-forge jupyter pandas fontconfig freetype libtiff r-irkernel scanpy -y
$ python -m pip install ucdeconvolve igraph leidenalg
```

### Export Environment File

For each operating system save the dependencies as:

```
$ conda env export > environments-[linux|windows|macos].yml
```
