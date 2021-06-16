# Variables
ver <- "2019-11-26"

# Set ID
id_comp <- sprintf("%s_%s_%s_%s_%s",
                   ver, tech, org, tissue, dwn_smplng)

# Hyperparameters
al <- 0.01

## PATHS and common variables
version_dir <- sprintf("%s",ver)
dir.create(path = version_dir, 
           showWarnings = F, 
           recursive = T)

an <- "analysis"
an_00 <- sprintf("%s/00_QC", an)
an_01 <- sprintf("%s/01_train_LDA", an)
an_02 <- sprintf("%s/02_topic_profiles_cluster", an)
an_03 <- sprintf("%s/03_lda_prediction", an)
an_04 <- sprintf("%s/04_prediction_assessment", an)
an_mouse <- sprintf("%s/mouse_brain", an)
an_nmf <- sprintf("%s/NMF_NNLS_approach", an)
an_pdac <- sprintf("%s/pancreas_PDAC", an)
an_tools <- sprintf("%s/tool_benchmarking", an)
plt_dir <- sprintf("%s/plots_%s", ver, ver) 
dir.create(path = plt_dir, 
           showWarnings = FALSE, 
           recursive = TRUE)

robj_dir <- sprintf("%s/R_objects_%s", ver, ver)
dir.create(path = robj_dir, 
           showWarnings = FALSE, 
           recursive = TRUE)

dir.create(path = sprintf("%s/%s",robj_dir,tech), 
           showWarnings = F, 
           recursive = T)

# Load fonts
extrafont::loadfonts(quiet = TRUE)

# Set arial as default font
theme_set <- function() {
  ggplot2::update_geom_defaults(
    "text",
    list(
      family = "Arial",
      # size = 5,
      fontface = "plain",
      color = "#000000", # black
      hjust = 0.5,
      vjust = -0.5)
  )
  ggplot2::theme(text = ggplot2::element_text(family = "Arial"))
                # add your other theme things here)
        }

ggplot <- function(...) {
  ggplot2::ggplot(...) + theme_set()
}
