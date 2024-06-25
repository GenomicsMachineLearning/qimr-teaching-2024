![Genomics and Machine Learning Lab, QIMR and UQ Logo](images/logo.png)

# GML Spatial and Machine Learning Teaching Material 2024

An overview of the materials to be covered in this course:
* [Module 000 - Single Cell RNAseq Data Analysis](#Module-000---Single-Cell-RNAseq-Data-Analysis)
  * Single Cell Data Analysis - Quan
* [Module 001 - Spatial Single Cell Visualisation](#Module-001---Spatial-Single-Cell)
  * Single Cell Spatial Data Visualisation - Levi
  * stLearn/Xenium Explorer - Xiao
* [Module 002 - Introducing Spatial Analysis](#Module-002---Clustering-and-Cell-Types)
  * Preprocessing - Andrew C
  * Clustering and Cell Typing - Andrew C
  * Deconvolution and Label Transfer - Andrew C
* [Module 003 - Downstream Analysis](#Module-003---Downstream-Analysis)
  * Inferring Malignant Cells using CNV Profiles - Prakrithi
  * Community Analysis - Feng
  * Cell-Cell Interaction Analysis - Onkar & Levi
* [Module 004 - Spatial Statistics](#Module-004---Spatial-Statistics)
  * Tissue Segmentation - Andrew N
  * Spatial Statistics with Voyager - Andrew N
* [Module 005 - Spatial Proteomics](#Module-005---Spatial-Proteomics)
  * Xiao & Quan
* [Module 006 - Deep Learning](#Module-006---Deep-Learning)
  * Xiao & Quan

# Data
* Module 0 scRNAseq (not about spatial omics)
  * [PBMC 3K scRNAseq data pbmc3k_filtered_gene_bc_matrices.tar.gz](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) (7.6 MB)
* Module 1
  * [Xenium_with_labels.zarr.tar.gz](https://downloads.gmllab.com/qimr-teaching-2024/Xenium_with_labels.zarr.tar.gz) (1.0 GB)
  * [cosmx.h5ad](https://downloads.gmllab.com/qimr-teaching-2024/cosmx.h5ad) (1.3 GB)
* Module 2
  * 2.1 - Cell Typing Tutorial
    * [visium.RDS](https://downloads.gmllab.com/qimr-teaching-2024/visium.RDS) (7.4 MB)
    * [xenium.RDS](https://downloads.gmllab.com/qimr-teaching-2024/xenium.RDS) (33 MB)
  * 2.2 - Cell Typing Example
    * [scRNA_processed.RDS](https://downloads.gmllab.com/qimr-teaching-2024/scRNA_processed.RDS) (364 MB)
    * [visium_processed.RDS](https://downloads.gmllab.com/qimr-teaching-2024/visium_processed.RDS) (36 MB)
    * [xenium_processed.RDS](https://downloads.gmllab.com/qimr-teaching-2024/xenium_processed.RDS) (98 MB)
* Module 3
  * 3.1 - CNV Profiling
    * [infercnv.tar.gz](https://downloads.gmllab.com/qimr-teaching-2024/infercnv.tar.gz) (191 MB)
  * 3.2 - Community Analysis
    * [CosMx_Skin_Melanoma.RData](https://downloads.gmllab.com/qimr-teaching-2024/CosMx_Skin_Melanoma.RData) (7.2 MB)
  * 3.3 - Neighborhood Coordination
    * [spatial.csv](https://downloads.gmllab.com/qimr-teaching-2024/spatial.csv) (2.6 MB)
  * 3.4 - Cell-Cell Interaction
    * [Visium_Skin_A2_cellchat.rds](https://downloads.gmllab.com/qimr-teaching-2024/Visium_Skin_A2_cellchat.rds) (14 MB)
    * [visium_decon.csv](https://downloads.gmllab.com/qimr-teaching-2024/visium_decon.csv) (375 KB)
    * [scalefactors_json.json](https://downloads.gmllab.com/qimr-teaching-2024/scalefactors_json.json) (204 Bytes)
* Module 4
  * 4.1 - Tissue Segmentation
    * [Visium_Mouse_Olfactory_Bulb.tar.gz](https://downloads.gmllab.com/qimr-teaching-2024/Visium_Mouse_Olfactory_Bulb.tar.gz) (31 MB)
    * [Visium_Skin_A2.tar.gz](https://downloads.gmllab.com/qimr-teaching-2024/Visium_Skin_A2.tar.gz) (Coming Soon) 
  * 4.2 - Spatial Statistics
    * [Visium_Mouse_Olfactory_Bulb.rds](https://downloads.gmllab.com/qimr-teaching-2024/Visium_Mouse_Olfactory_Bulb.rds) (28 MB)
* Module 5
  * [CODEX.tar.gz](https://downloads.gmllab.com/qimr-teaching-2024/CODEX.tar.gz) (447 MB)
* Module 6

# Running the Learning Materials

Copy and paste each of the following lines into your terminal once you have logged into the workshop server:
* ```/software/bin/micromamba shell init```
* ```source ~/.bashrc```
* ```micromamba activate /software/conda-envs/winter_school_2024```
* ```git clone https://github.com/GenomicsMachineLearning/qimr-teaching-2024```
* ```~/qimr-teaching-2024/runme.sh```

The output will look something like:
```bash
Port 3502 is available

Command to create ssh tunnel:
ssh -N -L 3502:10.10.10.10:3502 foo@10.10.10.10
Use a Browser on your local machine to go to:
localhost:3502  (prefix w/ https:// if using password)

[I 2024-06-20 05:57:41.633 ServerApp] Extension package jupyter_lsp took 0.1372s to import
[I 2024-06-20 05:57:44.647 ServerApp]     http://127.0.0.1:3502/tree?token=abc123
```

* Copy the line beginning with "ssh" into a new terminal, on your local computer, and hit [Enter].
* Copy the text beginning with "http://127.0.0.1" into a new tab in your browser, and hit [Enter].

# Module 000 - Single Cell RNAseq Data Analysis
* [scRNAseq Data Analysis](./000-single-cell-RNAseq/000_Single_Cell_RNAseq_Analysis_2024.ipynb)

# Module 001 - Spatial Single Cell 

* [Spatial Data Visualisation](./001-spatial-single-cell/single_cell_visualisations.ipynb).

# Module 002 - Clustering and Cell Types

* [Clustering and Cell Typing](002-clustering-cell-typing/2.2_Workshop_Cell_Typing_Example.ipynb).

# Module 003 - Downstream Analysis 

* [CNV Profiling](./003-downstream-analysis/3.1_QIMR_CNV_profiling.ipynb).
* [Cell community identification](./003-downstream-analysis/3.2_hoodscanR.ipynb).
* [Neighborhood Coordination and Cell Community Identification](./003-downstream-analysis/3.3_neighborhood.ipynb).
* [Cell-Cell Interactions CellChat](./003-downstream-analysis/3.4_CCI_CellChat.ipynb).
* [Cell-Cell Interactions stLearn](./003-downstream-analysis/3.5_CCI_stLearn_MMCCI.ipynb).

# Module 004 - Spatial Statistics

* [Tissue Segmentation](004-spatial-statistics/4.1-tissue_segmentation.ipynb).
* [Spatial Statistics](004-spatial-statistics/4.2-spatial-statistics.ipynb).

# Module 005 - Spatial Proteomics 

* [Mapping CODEX in Visium Data](./005-spatial-proteomics/mapping_CODEX_Visium.ipynb).
* [Spatial Proteomics Analysis for CODEX Data](./005-spatial-proteomics/CODEX_analysis.ipynb).

# Module 006 - Deep Learning 

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
$ conda env create --subdir osx-64 --prefix [some-directory]/conda-envs/qimr-teaching-2024 --file=environment-macos.yml -y
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
$ conda install -c conda-forge r-dendextend r-parallelDist r-mixtools r-ggalluvial r-svglite r-sna r-ggpubr r-ggnetwork r-matrix=1.6-3 r-scico r-ggnewscale r-magick r-rjson r-ragg r-units r-stringi r-sf r-s2 r-reticulate r-stringi r-tidyverse r-r.utils r-Seurat r-SeuratObject r-sctransform r-proj r-rcpptoml r-spdep r-lme4 r-ggrastr r-dbscan r-hdf5r r-optparse r-memuse r-sfheaders r-zeallot r-rmapshaper -y
$ conda install -c bioconda presto r-presto bioconductor-dropletutils r-MuSiC bioconductor-hoodscanr bioconductor-infercnv -y 
$ conda install -c conda-forge r-mcmcpack r-fields r-concaveman r-scatterpie r-ggcorrplot r-nnls r-pbmcapply r-NMF r-terra -y
```

Python Dependencies:
```
$ conda install -c conda-forge jupyter pandas fontconfig freetype libtiff r-irkernel numpy==1.26.4 scanpy -y
```

### Installing Unmanaged Dependencies from Source
This install dependencies that aren't managed by packages and need to be installed directly from source.

To install R dependencies run:
```
$ Rscript dependencies.R
```

To install Python dependencies run:
```
$ python -m pip install --use-pep517 -r requirements.txt
```

For MacOS, ensure you have Rust and rust-up installed in order to install Python dependencies 
(run ```rustup target add x86_64-apple-darwin``` before installing). This is for
gseapy and fastremap.

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

