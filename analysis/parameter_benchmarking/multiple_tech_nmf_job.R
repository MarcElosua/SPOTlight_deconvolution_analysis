#!/usr/bin/env Rscript

library(SPOTlight)
library(dplyr)
library(SingleCellExperiment)
library(Seurat)
library(purrr)

args <- commandArgs(trailingOnly=TRUE)
tech <- args[1]
i <- args[2]

dwn_smplng <- 'both'; 
tissue <- 'pbmc'
org <- "hs"
# tech <- "ICELL8"
source('misc/paths_vrs.R')
print(id_comp)

clust_vr <- "nnet2"
cl_n <- 20
hvg <- 3000
ntop <- 100
transf <- "uv"
method <- "nsNMF"
logFC <- 0
pct1 <- 0.5

id_nmf <- sprintf("cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s", 
                  cl_n, hvg, ntop, transf, method)

#### Load data ####
se_obj <- readRDS(file = sprintf("%s/%s/techs/seurat_umap_%s.RDS", 
                                 an_00, robj_dir, id_comp))

#### Discard Megakaryocytes ####
# Discard them bc they are not present in all the cells
se_obj <- se_obj[, se_obj$nnet2 != "Megakaryocytes"]

# Benchmarking time
start_time <- Sys.time()
out_ls <- spatial_decon_syn_assessment_nmf_fun(se_obj = se_obj,
                                               clust_vr = clust_vr,
                                               verbose = TRUE,
                                               cl_n = cl_n,
                                               hvg = hvg,
                                               ntop = ntop,
                                               transf = transf,
                                               method = method,
                                               logFC = logFC,
                                               pct1 = pct1)

total_time <- difftime(Sys.time(), start_time, units = "mins")
out_ls[[2]]["time_mins"] <- total_time

out_dir <- sprintf("%s/benchmark", an_nmf)
dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(object = out_ls,
        file = sprintf("%s/tech_benchmark_ls_iter-%s_%s_%s.RDS", 
                       out_dir, i, id_nmf, id_comp))

# for (tech in c("Chromium", "inDrop", "C1HT-medium", "C1HT-small", "CEL-Seq2", "ddSEQ", "Drop-Seq", "ICELL8", "MARS-Seq", "Chromium (sn)", "Quartz-Seq2", "SCRB-Seq", "Smart-Seq2")) {
#   dwn_smplng <- 'both';
#   tissue <- 'pbmc'
#   org <- "hs"
#   # tech <- "ICELL8"
#   source('misc/paths_vrs.R')
#   print(id_comp)
#   clust_vr <- "nnet2"
# 
#   se_obj <- readRDS(file = sprintf("%s/%s/techs/seurat_umap_%s.RDS",
#                                    an_00, robj_dir, id_comp))
# 
#   #### Discard Megakaryocytes ####
#   # Discard them bc they are not present in all the cells
#   se_obj <- se_obj[, se_obj$nnet2 != "Megakaryocytes"]
# 
#   se_obj <- SCTransform(se_obj)
# 
#   Seurat::Idents(object = se_obj) <- se_obj@meta.data[, clust_vr]
#   cluster_markers_all <- Seurat::FindAllMarkers(object = se_obj,
#                                                 verbose = TRUE,
#                                                 only.pos = TRUE,
#                                                 assay = "SCT",
#                                                 slot = "data")
# 
#   saveRDS(object = cluster_markers_all,
#           file = sprintf("%s/%s/cluster_markers_%s.RDS",
#                          an_mouse, robj_dir, id_comp))
# 
# }


# for (i in seq_len(10)) {
#   tech <- "Chromium (sn)"
# 
#   dwn_smplng <- 'both'; 
#   tissue <- 'pbmc'
#   org <- "hs"
#   source('misc/paths_vrs.R')
#   print(id_comp)
#   
#   clust_vr <- "nnet2"
#   cl_n <- 20
#   hvg <- 3000
#   ntop <- 100
#   transf <- "uv"
#   method <- "nsNMF"
#   logFC <- 0
#   pct1 <- 0.5
#   
#   id_nmf <- sprintf("cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s", 
#                     cl_n, hvg, ntop, transf, method)
#   
#   #### Load data ####
#   se_obj <- readRDS(file = sprintf("%s/%s/techs/seurat_umap_%s.RDS", 
#                                    an_00, robj_dir, id_comp))
#   
#   #### Discard Megakaryocytes ####
#   # Discard them bc they are not present in all the cells
#   se_obj <- se_obj[, se_obj$nnet2 != "Megakaryocytes"]
#   
#   # Benchmarking time
#   start_time <- Sys.time()
#   out_ls <- spatial_decon_syn_assessment_nmf_fun(se_obj = se_obj,
#                                                  clust_vr = clust_vr,
#                                                  verbose = TRUE,
#                                                  cl_n = cl_n,
#                                                  hvg = hvg,
#                                                  ntop = ntop,
#                                                  transf = transf,
#                                                  method = method,
#                                                  logFC = logFC,
#                                                  pct1 = pct1)
#   
#   total_time <- difftime(Sys.time(), start_time, units = "mins")
#   out_ls[[2]]["time_mins"] <- total_time
#   
#   out_dir <- sprintf("%s/benchmark", an_nmf)
#   dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
#   saveRDS(object = out_ls,
#           file = sprintf("%s/tech_benchmark_ls_iter-%s_%s_%s.RDS", 
#                          out_dir, i, id_nmf, id_comp))
# }
