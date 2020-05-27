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
# cl_n <- 20
hvg <- 0
ntop <- 100
transf <- "uv"
method <- "nsNMF"
logFC <- 1
pct1 <- 0.9

args <- commandArgs(trailingOnly=TRUE)
cl_n <- as.numeric(args[1])
i <- as.numeric(args[2])

id_nmf <- sprintf("cln-%s_hvg-%s_ntop-%s_transf-%s_method-%s", 
                  cl_n, hvg, ntop, transf, method)

#### Load data ####
allen_reference <- readRDS(file = "analysis/mouse_brain/R_obj/allen_cortex_processed.RDS")

start_t <- Sys.time()
out_ls <- spatial_decon_syn_assessment_nmf_fun(se_obj = allen_reference,
                                               clust_vr = clust_vr,
                                               verbose = TRUE,
                                               cl_n = cl_n,
                                               hvg = hvg,
                                               ntop = ntop,
                                               transf = transf,
                                               method = method, 
                                               logFC = logFC,
                                               pct1 = pct1)

total_time <- difftime(Sys.time(), start_t, units = "mins")
out_ls[[2]]["time"] <- total_time

out_dir <- sprintf("%s/benchmark", an_nmf)
dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(object = out_ls,
        file = sprintf("%s/cln_benchmark_ls_iter-%s_%s.RDS", 
                       out_dir, i, id_nmf))
