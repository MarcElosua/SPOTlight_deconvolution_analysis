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

clust_vr <- "subclass"
cl_n <- 20
# hvg <- 0
ntop <- 100
transf <- "uv"
method <- "nsNMF"

#### Load data ####
allen_reference <- readRDS(file = "analysis/mouse_brain/R_obj/allen_cortex_processed.RDS")

args <- commandArgs(trailingOnly=TRUE)
hvg_str <- args[1]
i <- args[2]

if (hvg_str == "uns") {
    hvg <- hvg_str
} else {
    hvg <- as.numeric(hvg_str)
}


# set.seed(234)
# Iterate n times over n_clust
id_nmf <- sprintf("cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s", 
                  cl_n, hvg_str, ntop, transf, method)

start_t <- Sys.time()
out_ls <- spatial_decon_syn_assessment_nmf_fun(se_obj = allen_reference,
                                               clust_vr = clust_vr,
                                               verbose = TRUE,
                                               cl_n = cl_n,
                                               hvg = hvg,
                                               ntop = ntop,
                                               transf = transf,
                                               method = method)

total_time <- difftime(Sys.time(), start_t, units = "mins")
out_ls[[2]]["time_mins"] <- total_time

out_dir <- sprintf("%s/benchmark", an_nmf)
dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(object = out_ls,
        file = sprintf("%s/hvg_benchmark_ls_iter-%s_%s.RDS", 
                       out_dir, i, id_nmf))
