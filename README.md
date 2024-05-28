# QIMR Spatial and Machine Learning Teaching Material 2024

An overview of the materials to be covered in this course:
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

# Acccessing the QIMR HPC

## Log into the HPC with Windows

* Download and Install MobaXterm
  * https://mobaxterm.mobatek.net/
* Use MobaXterm to connect to hpcpbs01

## Log into the HPC with MacOS

* Download and install XQuartz https://www.xquartz.org/

In a terminal, use SSH locally:
```
$ ssh -Y [username]@hpcpbs01
```

SSH into from VPC:
```
$ ssh -Y [username]@hpcpbs01.adqimr.ad.lan
```

Go to the teaching materials directory:
```
git clone https://github.com/GenomicsMachineLearning/qimr-teaching.git
cd qimr-teaching
```

## Start a new runner node

```
$ qsub -IX -l ncpus=2,mem=32GB,walltime=8:00:00
```

## Run IPython in a Conda environment

```
$ micromamba /path/conda-env/XXX
$ ipython
```

# Data

* One Visium, one Xenium, one CosmiX
  * Melanoma - single cell, Visium, Xenium, and CosmiX.
* Instructions on how to download

# Installation

To use the notebooks we require using a Conda package manager (or compatible tool like Mamba or Micromamba),
installing the packages with Conda and then installing additional packages directly into R and Python.

For each environment, the steps are:
* Download and install a package manager (Conda),
* Install packages with conda (using an environment YAML file), and
* Installing R and Python packages using R and Python (pip).

There are XXX environments that need to be setup (as some versions of some software are not compatible with others):
* qimr-teaching-2024 - for use with materials 001 and 002 (Voyager/Clustering/etc in R).
* 

# Package Managers

We recommend Conda for most users and operating systems. Micromamba is included here for use on the HPC and people
who are used to using package managers. 

## Setup Conda

Follow the instructions to download and install Conda: https://www.anaconda.com/download/success

## Setup Micromamba

Follow the instructions to install Micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

Once downloaded, ensure Micromamba is initialised:
```
$ micromamba shell init
```

Create ~/.condarc

```
channels:
  - conda-forge
  - bioconda
  - defaults
  - r
channel_priority: flexible
```

## Installing Managed Dependencies
This installs:
 * R and its dependencies for Seurat v5,
 * Python and dependencies for Jupyter.

### Using a package manager with and an environment YAML file (Recommended)

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

## Install Environment Specific Dependencies
This install dependencies that aren't managed by packages and need to be installed directly from source.

R Dependencies:
```
$ R
> install.packages("remotes", dependencies = FALSE)
> remotes::install_github("drighelli/SpatialExperiment", dependencies=FALSE)
> remotes::install_github("pachterlab/SpatialFeatureExperiment", ref="devel", dependencies=FALSE)
> remotes::install_github("pachterlab/Voyager", ref="devel", dependencies=FALSE)
> remotes::install_version('wrMisc', dependencies = FALSE)
> remotes::install_github('YingMa0107/CARD', dependencies = FALSE)
> remotes::install_github('jinworks/CellChat', dependencies = FALSE)
```

Python Dependencies:
```
$ python -m pip install poetry
$ poetry install
```

# Additional Information (Not Required for Course Participant)

## Manual Installation/Recreating Environment YAML

### Create and Activate a New Environment
MacOS:
```
$ conda create --name qimr-teaching-2024 --subdir osx-64 python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching-2024
```

Windows/MacOS:
```
$ conda create --name qimr-teaching-2024 python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching-2024
```

HPC:
```
$ micromamba create -p [some-directory]/conda-envs/qimr-teaching-2024 python=3.10 r-base=4.3 r-devtools -y
$ micromamba activate [some-directory]/conda-envs/qimr-teaching-2024
```

Replace the calls below with "micromamba" instaed of "conda".

### Install Dependencies

R Dependencies:
```
$ conda install -c bioconda bioconductor-scater bioconductor-scran bioconductor-splatter bioconductor-edgeR bioconductor-bluster bioconductor-BiocFileCache bioconductor-glmGamPoi bioconductor-SingleCellExperiment bioconductor-SummarizedExperiment bioconductor-ScaledMatrix bioconductor-BiocParallel bioconductor-ebimage r-anndata -y
$ conda install -c conda-forge r-ggalluvial r-svglite r-sna r-ggpubr r-ggnetwork r-matrix=1.6-3 r-scico r-ggnewscale r-magick r-rjson r-ragg r-units r-stringi r-sf r-s2 r-reticulate r-stringi r-tidyverse r-r.utils r-Seurat r-SeuratObject r-sctransform r-proj r-rcpptoml r-spdep r-lme4 r-ggrastr r-dbscan r-hdf5r r-optparse r-memuse r-sfheaders r-zeallot r-rmapshaper -y
$ conda install -c bioconda presto r-presto bioconductor-dropletutils r-MuSiC bioconductor-hoodscanr -y 
$ conda install -c conda-forger-mcmcpack r-fields r-concaveman r-scatterpie r-ggcorrplot r-nnls r-pbmcapply r-NMF -y
```

Python Dependencies:
```
$ conda install -c conda-forge jupyter pandas fontconfig freetype libtiff r-irkernel scanpy -y
```

**Follow the "Install Environment Specific Dependencies" section.**

### Export Environment File

In order to use the environments again, save the conda environment file and remove the fingerprint information:
```
$ conda env export > environments-[linux|windows|macos].yml
```
