#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://mirror.aarnet.edu.au/pub/CRAN/"))

install.packages("remotes", dependencies = FALSE)
remotes::install_version('dlm', dependencies = FALSE)
# Install version 1.15.2
remotes::install_github("drighelli/SpatialExperiment", dependencies=FALSE)
# Install version 1.7.0
remotes::install_github("pachterlab/SpatialFeatureExperiment", ref="devel", dependencies=FALSE)
# Install version 1.7.0
remotes::install_github("pachterlab/Voyager", ref="devel", dependencies=FALSE)
remotes::install_version('wrMisc', dependencies = FALSE)
# Install version 1.1
remotes::install_github('YingMa0107/CARD', dependencies = FALSE)
remotes::install_github('renozao/NMF', ref="0.30.4.900", dependencies = FALSE)
# Install version 2.1.2
remotes::install_github('jinworks/CellChat', dependencies = FALSE)
# Install version 1.1.0
remotes::install_github('navinlabcode/copykat', dependencies = FALSE)