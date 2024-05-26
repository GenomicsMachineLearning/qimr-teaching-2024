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

# Package Managers

We recommend Conda for most users and operating systems. Micromamba is included here for use on the HPC and people
who are used to using package managers. 

## Setup Conda

** Add instructions here for MacOS and Windows**

## Setup Micromamba
```
$ /working/joint_projects/P3903/teaching2024-winter-qimr/micromamba/micromamba shell init
```

### Configuration

Create ~/.condarc

```
channels:
  - conda-forge
  - bioconda
  - defaults
channel_priority: flexible
```

## Installing Managed Dependencies (Required)
This installs:
 * R and its dependencies for Seurat v5,
 * Python and dependencies for Jupyter.

### Using a package manager with and an environment.yaml (Recommended)

Windows:
```
C:\> conda env create --name qimr-teaching-2024 --file=environment-windows.yml -y
```

MacOS:
```
$ conda env create --subdir osx-64 --name qimr-teaching-2024 --file=environment-macos.yml -y
```

Linux:
```
$ conda env create --name qimr-teaching-2024 --file=environment-linux.yml -y
```

### Using a package manager (manual installation)

#### Create and Activate a New Environment
MacOS:
```
$ conda create --name qimr-teaching-2024 --subdir osx-64 python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching-2024
```

Others:
```
$ conda create --name qimr-teaching-2024 python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching-2024
```

HPC:
```
$ micromamba create -p /working/joint_projects/P3903/teaching2024-winter-qimr/conda-envs/base-python-R python=3.10 r-base=4.3 r-devtools -y
$ micromamba activate /mnt/lustre/working/joint_projects/P3903/teaching2024-winter-qimr/conda-envs/base-python-R
```

For OSes other than HPC replace all calls to micromamba to conda.

#### Install Dependencies

R Dependencies:
```
$ conda install -c bioconda bioconductor-scater bioconductor-splatter bioconductor-edgeR bioconductor-bluster bioconductor-BiocFileCache bioconductor-glmGamPoi bioconductor-SingleCellExperiment bioconductor-SummarizedExperiment bioconductor-ScaledMatrix bioconductor-BiocParallel bioconductor-ebimage r-anndata -y
$ conda install -c conda-forge r-scico r-ggnewscale r-magick r-rjson r-ragg r-units r-stringi r-sf r-s2 r-reticulate r-stringi r-tidyverse r-r.utils r-Seurat r-SeuratObject r-sctransform r-proj r-rcpptoml r-spdep r-lme4 r-ggrastr r-dbscan r-hdf5r r-optparse r-memuse r-sfheaders r-zeallot r-rmapshaper -y
$ conda install -c bioconda presto r-presto bioconductor-dropletutils r-MuSiC -y 
$ conda install -c conda-forge r-mcmcpack r-fields r-concaveman r-scatterpie r-ggcorrplot r-nnls r-pbmcapply r-NMF -y
```

Python Dependencies:
```
$ conda install -c conda-forge jupyter pandas fontconfig freetype libtiff r-irkernel scanpy -y
```

## Install Environment Specific Dependencies (Required)
This install dependencies that aren't managed by packages and need to be installed directly from source.

R Dependencies to install directly:
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
$ python -m pip install ucdeconvolve igraph leidenalg
```

### Export Environment File (not required)

In order to use the environments again, save the conda environment file and remove the fingerprint information:
```
$ conda env export > environments-[linux|windows|macos].yml
```
