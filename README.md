![Genomics and Machine Learning Lab, QIMR and UQ Logo](images/logo.png)

# QIMR Spatial and Machine Learning Teaching Material 2024

An overview of the materials to be covered in this course:
* Module 001 - Introducing Spatial Analysis
  * Preprocessing - Andrew C
  * Clustering and Cell Typing - Andrew C
  * Deconvolution and Label Transfer - Andrew C
* Module 002 - Spatial Statistics
  * Tissue Segmentation - Andrew N
  * Spatial Statistics with Voyager - Andrew N
* Module 003 - Downstream Analysis
  * Inferring Malignant Cells using CNV Profiles - Prakrithi
  * Community Analysis - Feng
  * Cell-Cell Interaction Analysis - Onkar & Levi
* Module 004 - Spatial Single Cell Visualisation
  * Xenium/CosMx Visualisation - Levi
  * stLearn/Xenium Explorer - Xiao
* Module 005 - Spatial Proteomics
  * Xiao & Quan
* Module 006 - Deep Learning
  * Xiao & Quan

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
$ ssh -YC [username]@hpcpbs01
```

SSH into from VPC:
```
$ ssh -YC [username]@hpcpbs01.adqimr.ad.lan
```

## Start a new runner node

```
$  qsub -IX -q spatial -l ncpus=2,mem=24GB,walltime=8:00:00
```

## Using a Conda Environment

There is one Conda environments that will be used for all training materials:
* qimr-teaching-2024

```
$ /scratch/qimr-teaching-2024/micromamba/micromamba shell init
$ source ~/.bashrc
$ micromamba activate /scratch/qimr-teaching-2024/conda-envs/qimr-teaching-2024
```

# Data

* On scratch on HPC075 (all data is here on the training day)
* Also available via link to public data on AWS - tar.gz for Visium and Xenium (for you to work at another time).

# Module 001 - Clustering and Cell Types

View the [notebook here](./001-clustering-cell-typing/notebooks/1.2_Workshop_Cell_Typing_Example.ipynb).

From the workshop directory go to the notebooks directory and being an R prompt:
```
$ R
```

# Module 002 - Spatial Statistics

View the [notebook here](./002-spatial-statistics/notebooks/voyager.ipynb).

From the workshop directory go to the notebooks directory and being an R prompt:
```
$ cd 002-spatial-statistics/notebooks
$ R
```

You should be able to open the "voyager.ipyn" file in your IDE or web browser and follow the instructions.

# Module 003 - Spatial Single Cell Visualization 

```
$ cd 003-spatial-single-cell
$ ipython
```

# Module 004 - Multimodalities 

# Module 005 - Community Analysis 

# Module 006 - Spatial Proteomics

# Module 007 - Cell-Cell Interaction
   
# Module 008 - Deep Learning

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
$ conda create --subdir osx-64 --prefix [some-directory]/conda-envs/qimr-teaching-2024 --file=environment-macos.yml
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

R Dependencies:
```
$ R
> install.packages("remotes", dependencies = FALSE)
> remotes::install_github("drighelli/SpatialExperiment", dependencies=FALSE)
> remotes::install_github("pachterlab/SpatialFeatureExperiment", ref="devel", dependencies=FALSE)
> remotes::install_github("pachterlab/Voyager", ref="devel", dependencies=FALSE)
> remotes::install_version('wrMisc', dependencies = FALSE)
> remotes::install_github('YingMa0107/CARD', dependencies = FALSE)
> remotes::install_github('renozao/NMF', ref="0.30.4.900", dependencies = FALSE)
> remotes::install_github('jinworks/CellChat', dependencies = FALSE)
> remotes::install_github('navinlabcode/copykat', dependencies = FALSE)
```

Python Dependencies:
```
$ python -m pip install poetry
$ poetry install
```

OR:
```
pip install "spatialdata[extra]>=0.1.2" "torch>=2.2.2" "torchvision>=0.17.2" "lightning>=2.2.2" "pyro-ppl>=1.9.0" "squidpy>=1.4.1" "monai>=1.3.0" "plotly>=5.22.0" "sopa[baysor,cellpose,snakemake,tangram]>=1.0.14" "jupyterlab>=4.2.1" "ipywidgets>=8.1.2" "monkeybread>=1.0.3" "harmonypy>=0.0.9" "matplotlib<3.9.0" "netgraph>=4.13.2" "python-louvain>=0.16"
pip install multimodal_cci
pip install "git+https://github.com/xiao233333/stLearn.git@update-dependency" 
```

### Cloning an Environment
An example of copying it from group directory (P3903) to a temporary directory on a local scratch:

```
$ conda create --prefix /scratch/qimr-teaching-2024/conda-envs/qimr-teaching-2024 --clone [some-directory]/teaching2024-winter-qimr/conda-envs/xiao_pyr2
```

### Export Environment File

In order to use the environments again, save the conda environment file and remove the fingerprint information:
```
$ conda env export > environment-[linux|windows|macos].yml
```

