#!/usr/bin/env Rscript

library(SPOTlight)
library(SingleCellExperiment)
library(Seurat)

args <- commandArgs(trailingOnly=TRUE)
dwn_smplng <- as.numeric(args[1])
i <- as.numeric(args[2])

tech <- "Chromiumv3"; tissue <- "pbmc"
org <- "hs"
source("misc/paths_vrs.R")
print(id_comp)

clust_vr <- "seurat_clusters"
cl_n <- 20
hvg <- 0
ntop <- 100
transf <- "uv"
method <- "nsNMF"

id_nmf <- sprintf("cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s", 
                  cl_n, hvg, ntop, transf, method)

#### Load data ####
se_obj <- readRDS(file = sprintf("%s/%s/downsampling/seurat_umap_%s.RDS", 
                                 an_00, robj_dir, id_comp))

# Remove clusters 11-13 to be able to compare the different datasets
keep_id_downsampling <- rownames(se_obj@meta.data[! se_obj$seurat_clusters %in% c("11", "12", "13"), ])
se_obj@meta.data[, clust_vr] <- paste0("clust_", se_obj@meta.data[, clust_vr])

se_obj <- se_obj[, rownames(se_obj@meta.data) %in% keep_id_downsampling]

# Benchmarking time
start_time <- Sys.time()
out_ls <- spatial_decon_syn_assessment_nmf_fun(se_obj = se_obj,
                                               clust_vr = clust_vr,
                                               verbose = TRUE,
                                               cl_n = cl_n,
                                               hvg = hvg,
                                               ntop = ntop,
                                               transf = transf,
                                               method = method)

total_time <- difftime(Sys.time(), start_time, units = "mins")
out_ls[[2]]["time_mins"] <- total_time

out_dir <- sprintf("%s/benchmark", an_nmf)
dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(object = out_ls,
        file = sprintf("%s/dwnsmpl_benchmark_ls_iter-%s_%s_%s.RDS", 
                       out_dir, i, id_nmf, id_comp))

# for (dwn_smplng in c("5000", "10000", "15000", "20000", "50000")) {
#   tissue <- 'pbmc'
#   org <- "hs"
#   tech <- "Chromiumv3"
#   source('misc/paths_vrs.R')
#   print(id_comp)
#   clust_vr <- "seurat_clusters"
#   
#   id_nmf <- sprintf("cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s", 
#                     cl_n, hvg, ntop, transf, method)
#   se_obj <- readRDS(file = sprintf("%s/%s/downsampling/seurat_umap_%s.RDS", 
#                                    an_00, robj_dir, id_comp))
#   # Remove clusters 11-13 to be able to compare the different datasets
#   keep_id_downsampling <- rownames(se_obj@meta.data[! se_obj$seurat_clusters %in% c("11", "12", "13"), ])
#   
#   #### Discard Megakaryocytes ####
#   # Discard them bc they are not present in all the cells
#   se_obj <- se_obj[, rownames(se_obj@meta.data) %in% keep_id_downsampling]
#   
#   se_obj <- SCTransform(se_obj)
# 
#   Seurat::Idents(object = se_obj) <- se_obj@meta.data[, clust_vr]
#   cluster_markers_all <- Seurat::FindAllMarkers(object = se_obj,
#                                                 verbose = TRUE,
#                                                 only.pos = TRUE,
#                                                 assay = "SCT",
#                                                 slot = "data",
#                                                 min.pct = 0.9,
#                                                 max.cells.per.ident = 100)
# 
#   saveRDS(object = cluster_markers_all,
#           file = sprintf("%s/%s/cluster_markers_%s.RDS",
#                          an_mouse, robj_dir, id_comp))
# 
# }
