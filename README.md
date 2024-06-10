![Genomics and Machine Learning Lab, QIMR and UQ Logo](images/logo.png)

# QIMR Spatial and Machine Learning Teaching Material 2024

An overview of the materials to be covered in this course:
* [Module 001 - Introducing Spatial Analysis](#Module-001---Clustering-and-Cell-Types)
  * Preprocessing - Andrew C
  * Clustering and Cell Typing - Andrew C
  * Deconvolution and Label Transfer - Andrew C
* [Module 002 - Spatial Statistics](#Module-002---Spatial-Statistics)
  * Tissue Segmentation - Andrew N
  * Spatial Statistics with Voyager - Andrew N
* [Module 003 - Downstream Analysis](#Module-003---Downstream-Analysis)
  * Inferring Malignant Cells using CNV Profiles - Prakrithi
  * Community Analysis - Feng
  * Cell-Cell Interaction Analysis - Onkar & Levi
* [Module 004 - Spatial Single Cell Visualisation](#Module-004---Spatial-Single-Cell)
  * Single Cell Spatial Data Visualisation - Levi
  * stLearn/Xenium Explorer - Xiao
* [Module 005 - Spatial Proteomics](#Module-005---Spatial-Proteomics)
  * Xiao & Quan
* [Module 006 - Deep Learning](#Module-006---Deep-Learning)
  * Xiao & Quan

# Data

* Module 1
  * 1.1 - Cell Typing Tutorial
    * [visium.RDS](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/visium.RDS)
    * [xenium.RDS](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/xenium.RDS)
  * 1.2 - Cell Typing Example
    * [scRNA_processed.RDS](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/scRNA_processed.RDS)
    * [visium_processed.RDS](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/visium_processed.RDS)
    * [xenium_processed.RDS](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/xenium_processed.RDS) 
* Module 2
  * 2.1 - Tissue Segmentation
    * [Visium_Mouse_Olfactory_Bulb.tar.gz](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/Visium_Mouse_Olfactory_Bulb.tar.gz)
    * [Visium_Skin_A2.tar.gz](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/Visium_Skin_A2.tar.gz) (Soon) 
  * 2.2 - Spatial Statistics
    * [Visium_Mouse_Olfactory_Bulb.rds](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/Visium_Mouse_Olfactory_Bulb.rds)
* Module 3
  * 3.1 - CNV Profiling
  * 3.2 - Community Analysis
    * [CosMx_Skin_Melanoma.RData](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/CosMx_Skin_Melanoma.RData)
  * 3.3 - Neighborhood Coordination
    * [spatial.csv](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/spatial.csv)
  * 3.4 - Cell-Cell Interaction
    * [Visium_Skin_A2_cellchat.rds](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/Visium_Skin_A2_cellchat.rds)
    * [visium_decon.csv](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/visium_decon.csv)
    * [scalefactors_json.json](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/scalefactors_json.json)
* Module 4
  * [Xenium_with_labels.zarr.tar.gz](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/Xenium_with_labels.zarr.tar.gz)
  * [cosmx.h5ad](https://workshop-2024.s3.ap-southeast-2.amazonaws.com/cosmx.h5ad)
* Module 5
* Module 6

# Accessing the Training Materials on the QIMR HPC

The steps you will need to perform:
* Log into the HPC,
* Use a runtime environment, and
* Run iPython or R.

## Log into the HPC with Windows

* Download and Install MobaXterm
  * https://mobaxterm.mobatek.net/
* Use MobaXterm to connect to hpcpbs01

## Log into the HPC with MacOS

* Download and install XQuartz https://www.xquartz.org/

In a terminal, use SSH locally:
```
$ ssh -CY [username]@hpcpbs01
```

SSH into from VPC:
```
$ ssh -CY [username]@hpcpbs01.adqimr.ad.lan
```

## Start a new runner node

```
$  qsub -IX -q spatial -l ncpus=2,mem=24GB,walltime=8:00:00
```

## Using a Conda Environment

There is one Conda environments that will be used for all training materials:
* python-r

```
$ /scratch/qimr-spatial-teaching-2024/micromamba/micromamba shell init
$ source ~/.bashrc
$ micromamba activate /scratch/qimr-spatial-teaching-2024/conda-envs/python-r
```


# Module 001 - Clustering and Cell Types

From the workshop directory start an R prompt:
```bash
$ R
```

View the [notebook here](001-clustering-cell-typing/1.2_Workshop_Cell_Typing_Example.ipynb).

# Module 002 - Spatial Statistics

From the workshop directory go to the notebooks directory and being an R prompt:
```bash
$ cd /scratch/qimr-spatial-teaching-2024/002-spatial-statistics
$ R
```

You can view the notebooks here:
* [Tissue Segmentation](./002-spatial-statistics/2.1-tissue_segmentation.ipynb).
* [Spatial Statistics](./002-spatial-statistics/2.2-spatial-statistics.ipynb).

# Module 003 - Downstream Analysis 

You can view the notebooks here:
* [CNV Profiling](./003-downstream-analysis/3.1_QIMR_CNV_profiling.ipynb).
* [Cell community identification](./003-downstream-analysis/3.2_hoodscanR.ipynb).
* [Neighborhood Coordination and Cell Community Identification](./003-downstream-analysis/3.3_neighborhood.ipynb).
* [Cell-Cell Interactions](./003-downstream-analysis/3.4_Cell_Cell_Interaction.html).

# Module 004 - Spatial Single Cell 

```bash
cd /scratch/qimr-spatial-teaching-2024/004-spatial-single-cell
ipython
```
You can view the notebook here:
* [Spatial Data Visualisation](./004-spatial-single-cell/single_cell_visualisations.ipynb).

# Module 005 - Spatial Proteomics 

```bash
cd /scratch/qimr-spatial-teaching-2024/005-spatial-proteomics
ipython
```
You can view the notebooks here:
* [Mapping CODEX in Visium Data](./005-spatial-proteomics/mapping_CODEX_Visium.ipynb).
* [Spatial Proteomics Analysis for CODEX Data](./005-spatial-proteomics/CODEX_analysis.ipynb).

# Module 006 - Deep Learning 

```bash
cd /scratch/qimr-spatial-teaching-2024/006-deep-learning
ipython
```

You can view the notebook here:
* [Deep Learning](./006-deep-learning/Deep_learning_tutorial.ipynb).

# Additional Information (Not Required for Course Participants)

The following is not required for participation in the course but is described here to document how the
requirements were set up and to help you install your own  environment in the future. 

## Package Managers

We recommend Conda package management for most users and operating systems. Micromamba is included here for use on the 
HPC and people who are used to using Conda or other package managers. 

Download and install either:
* Conda: https://www.anaconda.com/download/success
* Micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html

### Using a package manager with and an environment YAML file (Recommended)

The steps are:
* Create a new environment with the OS specific YAML file,
* Install any extra software directly from source (see "Installing Dependencies from Source").

MacOS:
```
$ conda env create --subdir osx-64 --name qimr-teaching-2024 --file=environment-macos.yml -y
```

In a custom directory:
```
$ conda create --subdir osx-64 --prefix [some-directory]/conda-envs/qimr-teaching-2024 --file=environment-macos.yml -y
```

Currently, all dependencies are only available for x86 (Intel). If you are running an M-series CPU 
(M1, M2, M3, M4, etc.) it will run these under emulation.

Linux:

In your default conda directory:
```
$ conda env create --name qimr-teaching-2024 --file=environment-linux.yml -y
```

In a custom directory:
```
$ conda create --prefix [some-directory]/conda-envs/qimr-teaching-2024 --file=environment-linux.yml
```

### Recreating the Conda Environment YAML

If the above YAML fails to install, or you wish to upgrade a dependency, or you wish to recreate it. Run the following
to recreate the Conda environment. 

The steps are:
* Configure the package manager,
* Create a base environment and install packages with conda/micromamba, and
* Installing any further source based dependencies in R and Python (using pip or poetry) directly.

#### Configuring the Package Manager

Firstly, we need to configure conda/micromamba in how it resolves dependencies. We assume you've already installed the
package manager previously.

Create ~/.condarc

```
channels:
  - conda-forge
  - bioconda
  - defaults
  - r
channel_priority: flexible
```

### Create and Activate a New Environment

Next, we need to create the base environment. We're going to assume using "conda" in the following steps, but you can
replace the calls to "conda" with "micromamba" (if that's what you're using).

MacOS:
```
$ conda create --name qimr-teaching-2024 --subdir osx-64 python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching-2024
```

Linux:
```
$ conda create --name qimr-teaching-2024 python=3.10 r-base=4.3 r-devtools -y
$ conda activate qimr-teaching-2024
```

HPC:
```
$ micromamba create -p [some-directory]/conda-envs/qimr-teaching-2024 python=3.10 r-base=4.3 r-devtools -y
$ micromamba activate [some-directory]/conda-envs/qimr-teaching-2024
```

### Installing Managed Dependencies

R Dependencies:
```
$ conda install -c bioconda bioconductor-ebimage bioconductor-biocneighbors bioconductor-scater bioconductor-scran bioconductor-splatter bioconductor-edgeR bioconductor-bluster bioconductor-BiocFileCache bioconductor-glmGamPoi bioconductor-SingleCellExperiment bioconductor-SummarizedExperiment bioconductor-ScaledMatrix bioconductor-BiocParallel bioconductor-ebimage r-anndata -y
$ conda install -c conda-forge r-dendextend r-parallelDist r-dlm r-mixtools r-ggalluvial r-svglite r-sna r-ggpubr r-ggnetwork r-matrix=1.6-3 r-scico r-ggnewscale r-magick r-rjson r-ragg r-units r-stringi r-sf r-s2 r-reticulate r-stringi r-tidyverse r-r.utils r-Seurat r-SeuratObject r-sctransform r-proj r-rcpptoml r-spdep r-lme4 r-ggrastr r-dbscan r-hdf5r r-optparse r-memuse r-sfheaders r-zeallot r-rmapshaper -y
$ conda install -c bioconda presto r-presto bioconductor-dropletutils r-MuSiC bioconductor-hoodscanr bioconductor-infercnv -y 
$ conda install -c conda-forge r-mcmcpack r-fields r-concaveman r-scatterpie r-ggcorrplot r-nnls r-pbmcapply r-NMF -y
```

Python Dependencies:
```
$ conda install -c conda-forge jupyter pandas fontconfig freetype libtiff r-irkernel scanpy -y
```

### Installing Unmanaged Dependencies from Source
This install dependencies that aren't managed by packages and need to be installed directly from source.

To install R dependencies run:
```
$ Rscript dependencies.R
```

To install Python dependencies run:
```
$ python -m pip install -r requirements.txt
```

### Cloning an Environment
An example of copying it from group directory (P3903) to a temporary directory on a local scratch:

```
$ conda create --prefix [some-directory]/qimr-spatial-teaching-2024/conda-envs/qimr-teaching-2024 --clone [some-directory]qimr-spatial-teaching-2024/conda-envs/source-conda-dir
```

### Export Environment File

In order to use the environments again, save the conda environment file and remove the fingerprint information:
```
$ conda env export > environment-[linux|windows|macos].yml
```

