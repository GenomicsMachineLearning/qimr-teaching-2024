#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://mirror.aarnet.edu.au/pub/CRAN/"))

install.packages("remotes", dependencies = FALSE)
remotes::install_github("drighelli/SpatialExperiment", dependencies=FALSE)
remotes::install_github("pachterlab/SpatialFeatureExperiment", ref="devel", dependencies=FALSE)
remotes::install_github("pachterlab/Voyager", ref="devel", dependencies=FALSE)
remotes::install_version('wrMisc', dependencies = FALSE)
remotes::install_github('YingMa0107/CARD', dependencies = FALSE)
remotes::install_github('renozao/NMF', ref="0.30.4.900", dependencies = FALSE)
remotes::install_github('jinworks/CellChat', dependencies = FALSE)
remotes::install_github('navinlabcode/copykat', dependencies = FALSE)