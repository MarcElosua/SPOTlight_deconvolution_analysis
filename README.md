# SPOTlight: Seeded NMF regression to Deconvolute Spatial Transcriptomics Spots with Single-Cell Transcriptomes

This repository contains all the scripts notebooks and reports to reproduce the analysis of our paper SPOTlight: Seeded NMF regression to Deconvolute Spatial Transcriptomics Spots with Single-Cell Transcriptomes. Here we describe which data we used, how to access it, the most important packages and versions used, and how to navigate the directories and files in this repository.

## Data

All the data used for this study was publicly available at the Gene Expression Omnibus or through the pertinent web locations. The data in this project can be divided by tissue as well as by modality.

### Mus Musculus Brain
* [Cortical + Hippocampal scRNAseq](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq): Downloaded from the Allen Brain Map, includes 76,533 cells. Specifications on how the data was generated can be found [here](https://portal.brain-map.org/atlases-and-data/rnaseq/protocols-mouse-cortex-and-hippocampus).

* [Anterior + Posterior ST](https://support.10xgenomics.com/spatial-gene-expression/datasets/): Downloaded from the 10X website, includes 2 Anterior and 2 Posterior slices of adult mice.

### PDAC
* [PDAC scRNAseq/ST](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111672): PDAC scRNAseq and ST data was obtained from the paper [Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas](https://pubmed.ncbi.nlm.nih.gov/31932730/). Data for PDAC-A can be downloaded through GSM3036909, GSM3036910, GSM3405527, GSM3405528, and PDAC-B GSM3405531, GSM3405532, GSM3405533. ST data used in this paper can be accessed through GSM3036911 and GSM4100723.
* [PDAC immune scRNAseq](https://bigd.big.ac.cn/bioproject/browse/PRJCA001063): The immune reference dataset used in this study was obtained from the paper [ Single-cell RNA-seq Highlights Intra-Tumoral Heterogeneity and Malignant Progression in Pancreatic Ductal Adenocarcinoma ](https://pubmed.ncbi.nlm.nih.gov/31273297/). It contains 41,986 cells from which 10,623 are immune cells and can be accessed through the Genome Sequence Archive under project [PRJCA001063](https://bigd.big.ac.cn/bioproject/browse/PRJCA001063).

### PBMC
* [Benchmarking scRNAseq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133549): This data was obtained from the paper [Benchmarking single-cell RNA-sequencing protocols for cell atlas projects](https://www.nature.com/articles/s41587-020-0469-4?draft=marketing). It contains scRNAseq data collected with different protocols and sequencing depths.

## Package versions
These are the versions of the most important packages used throughout all the analysis:
CRAN:
* [tidyverse 1.3.0](https://cran.r-project.org/web/packages/tidyverse/vignettes/paper.html)
* [NMF 0.22.0](https://cran.r-project.org/web/packages/NMF/index.html)
* [Seurat 3.1.4.9903](https://satijalab.org/seurat/v3.1/spatial_vignette.html)

## Docker environments
Docker environments with all the pertinent dependencies are available at dockerhub.

To download the R environment image

    # R environment
    docker pull marcelosua/spotlight_env_r:latest

To download the R studio image

    # Rstudio environment
    docker pull marcelosua/spotlight_env_rstudio:latest

## File system
This repository contains 3 main analysis directories, each corresponding to the 3 main analysis of the paper. Scripts within each folder are numbered in the order they are to be run.
>>>>>>> 671522b044bb568d9c1a44fd72583574ad8e5d33
