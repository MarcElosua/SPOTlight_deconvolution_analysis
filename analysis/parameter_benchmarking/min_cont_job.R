#!/usr/bin/env Rscript

library(NMF)
library(nnls)
library(topicmodels)
library(dplyr)
library(ggplot2)
library(purrr)
library(Seurat)
library(SPOTlight)
library(edgeR)

#### Setting vrs ####
tech <- "sc"
tissue <- "allen_ref_14k"
dwn_smplng <- "both"
org <- "mm"
source("misc/paths_vrs.R")

# clust_vr <- "subclass_label"
clust_vr <- "subclass"
cl_n <- 100
method <- "nsNMF"
transf <- "uv"
hvg <- 3000
logFC <- 1
pct1 <- 0.9

#### Load downsampled data ####
# se_obj <- readRDS(file = "data/MusMusculus/allen_reference/allen_ref_70k_processed.RDS")
se_obj <- readRDS(file = "analysis/mouse_brain/R_obj/allen_cortex_processed.RDS")

args <- commandArgs(trailingOnly=TRUE)
min_cont <- as.numeric(args[1])
i <- args[2]

# set.seed(234)
# Iterate n times over n_clust
id_nmf <- sprintf("cln-%s_hvg-%s_transf-%s_method-%s-%s_min_cont", 
                  cl_n, hvg, transf, method, min_cont)

se_obj@meta.data[, clust_vr] <- gsub(pattern = "[[:punct:]]|[[:blank:]]",
                                     ".",
                                     x = se_obj@meta.data[, clust_vr],
                                     perl = TRUE)

##########################
### Generate test data ###
##########################
test_spot_ls <- test_spot_fun(se_obj = se_obj,
                              clust_vr = clust_vr,
                              n = 1000)

test_spot_counts <- test_spot_ls[[1]]
colnames(test_spot_counts) <- paste("mixt", 1:ncol(test_spot_counts), sep = "_")
test_spot_metadata <- test_spot_ls[[2]]

############################
#### Preprocessing data ####
############################
### Marker genes
#### Extract the top marker genes from each cluster ####
# Seurat::Idents(object = se_obj) <- se_obj@meta.data[, clust_vr]
# cluster_markers_all <- Seurat::FindAllMarkers(object = se_obj,
#                                               verbose = TRUE,
#                                               only.pos = TRUE,
#                                               assay = "SCT",
#                                               slot = "data",
#                                               min.pct = 0.9
#                                               # max.cells.per.ident = 100
# )

print("Loading markers")
# cluster_markers_all <- readRDS(file = "data/MusMusculus/allen_reference/markers_allen_ref_70k.RDS")
# saveRDS(cluster_markers_all, file = "data/MusMusculus/allen_reference/markers_allen_ref_14k.RDS")
cluster_markers_all <- readRDS(file = "data/MusMusculus/allen_reference/markers_allen_ref_14k.RDS")

cluster_markers_filt <- cluster_markers_all %>%
  filter(avg_logFC > logFC & pct.1 > pct1)

####################
### Downsampling ###
####################
# Downsample number of genes and number of samples
print("Downsampling genes and cells")
se_obj_down <- downsample_se_obj(se_obj = se_obj,
                                 clust_vr = clust_vr,
                                 cluster_markers = cluster_markers_filt,
                                 cl_n = cl_n,
                                 hvg = hvg)
###################
#### Train NMF ####
###################
print("Train NMF")
nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers_filt,
                        se_sc = se_obj_down,
                        mtrx_spatial = test_spot_counts,
                        transf = transf,
                        clust_vr = clust_vr,
                        method = method)
# dir.create(sprintf("%s/%s", an_nmf, robj_dir), showWarnings = FALSE, recursive = TRUE)
saveRDS(object = nmf_mod_ls,
        file = sprintf("%s/%s/nmf_mod_ls_%s_%s", an_nmf, robj_dir, id_nmf, i))
nmf_mod <- nmf_mod_ls[[1]]

# get matrix H
h <- coef(nmf_mod)

#################################
#### Get mixture composition ####
#################################
# Run test spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit.
ct_topic_profiles <- topic_profile_per_cluster_nmf(h = h,
                                                   train_cell_clust = nmf_mod_ls[[2]],
                                                   clust_vr = clust_vr)
print("Deconvolute synthetic spots")
pred_comp <- mixture_deconvolution_nmf(nmf_mod = nmf_mod,
                                       mixture_transcriptome = test_spot_counts,
                                       transf = transf,
                                       reference_profiles = ct_topic_profiles, 
                                       min_cont = min_cont)

################################
#### Performance statistics ####
################################
ct_cols <- colnames(pred_comp)[which(colnames(pred_comp) != "res_ss")]
raw_statistics_ls <- test_synthetic_performance(test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
                                                spot_composition_mtrx = pred_comp[, ct_cols])

out_dir <- sprintf("%s/min_cont", an_nmf)
dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(object = raw_statistics_ls,
        file = sprintf("%s/min_cont_benchmark_ls_iter-%s_%s.RDS", 
                       out_dir, i, id_nmf))
