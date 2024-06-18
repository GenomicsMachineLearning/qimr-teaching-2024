set.seed(29)
(g_all <- SpatialFeatureExperiment::findSpatialNeighbors(processed_sfe, MARGIN = 2, method = "tri2nb"))
(g_specific <- SpatialFeatureExperiment::findSpatialNeighbors(sfe_tissue, MARGIN = 2, method = "tri2nb"))
spatialGraph(processed_sfe, "graph1", MARGIN = 2) <- g_all
spatialGraph(sfe_tissue, "graph1", MARGIN = 2) <- g_specific
colGraph(sfe_tissue, "visium") <- SpatialFeatureExperiment::findVisiumGraph(sfe_tissue, zero.policy = TRUE)
colGraph(sfe_tissue, "visium_B") <- SpatialFeatureExperiment::findVisiumGraph(sfe_tissue, style = "B", zero.policy = TRUE)
sfe_tissue <- Voyager::runUnivariate(sfe_tissue, type = "localG_perm", features = "Pcp4", colGraphName = "visium_B",
                                     swap_rownames = "symbol")
dec <- scran::modelGeneVar(sfe_tissue)
hvgs <- scran::getTopHVGs(dec, n = 2000)
sfe_tissue <- Voyager::colDataUnivariate(sfe_tissue, features = c("nCounts", "nGenes"), colGraphName = "graph1", type = "moran")
sfe_tissue <- Voyager::runMoransI(sfe_tissue, features = hvgs)
# colFeatureData(sfe_tissue)[c("nCounts", "nGenes"),]
pos_moran <- rownames(sfe_tissue)[order(rowData(sfe_tissue)$moran_sample01, decreasing = TRUE)[1:6]]
neg_moran <- rownames(sfe_tissue)[order(rowData(sfe_tissue)$moran_sample01, decreasing = FALSE)[1:6]]
sfe_tissue <- Voyager::runUnivariate(sfe_tissue, "moran.mc", neg_moran, colGraphName = "graph1", nsim = 200, alternative = "less")
colGraph(sfe_tissue, "knn5") <- SpatialFeatureExperiment::findSpatialNeighbors(sfe_tissue, method = "knearneigh", dist_type = "idw", k = 5, style = "W")
sfe_tissue <- Voyager::runMultivariate(sfe_tissue, "multispati", colGraphName = "knn5", nfposi = 20, nfnega = 20)
sfe_tissue$clusts_nonspatial <- scran::clusterCells(sfe_tissue, use.dimred = "PCA",
                                             BLUSPARAM = NNGraphParam(
                                              cluster.fun = "leiden",
                                              cluster.args = list(
                                                  objective_function = "modularity",
                                                  resolution_parameter = 1
                                              )
                                          ))
sfe_tissue$clusts_multispati <- bluster::clusterRows(SingleCellExperiment::reducedDim(sfe_tissue, "multispati")[,1:20],
                                            BLUSPARAM = NNGraphParam(
                                                cluster.fun = "leiden",
                                                cluster.args = list(
                                                    objective_function = "modularity",
                                                    resolution_parameter = 1
                                                )
                                            )
                                           )