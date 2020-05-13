# COVID-19
The R/Python scripts for the analysis of single-cell RNA-seq data from COVID-19 patients

## 1. Requirement
We analyzed the scRNA-seq data in a Linux system with R (version 3.6.1) and Python (version 3.6.8) enviroment. The following software and packages are also required:

software|version|enviroment
-|-|-
Cellranger|3.1.0|Linux
Seurat|3.1.4|R
dplyr|0.8.4|R
patchwork|1.0.0|R
scrublet|0.2.1|Python
scipy|1.0.0|Python
pandas|0.24.2|Python
matplotlib|3.0.3|Python
seaborn|0.9.0|Python

## 2. Installation
Users need to copy the scripts to the same path as the raw data folders (i.e. the path containing the "P1-1r1/", "P1-1r2/", "P1-2r1/", "P1-2r2/", "P2-1/", "P2-2/", and "P2-3/" folders), and run these scripts in R studio or Jupyter notebook.

## 3. Step by step analysis

### 3.1 Map the raw sequencing data to the genome reference

    bash step0_map_raw_data_using_cellranger.sh

### 3.2 Detect doublet by using Scrublet
Run **step1_detect_doublets_using_scrublet.ipynb** in Jupyter notebook. This script decects doublets from the cells of patients and healthy controls respectively.

### 3.3 Integrate patients with healthy controls
Run **step2_integrate_patients_and_healthy_controls.R** in R studio. This script also cluster cells and perform UMAP analysis on the scRNA-seq data.

### 3.4 Plot UMAP diagram and marker gene expression violinplot
Run **step3_plot_umap_and_marker_gene_expression.ipynb** in Jupyter notebook. This script illustrates UMAP diagram for the cells and violinplot for the expression of marker genes. All diagrams will be automatically presented in the jupyter interface.

### 3.5 Generate DEGs between different disease stages 
Run **step4_find_DEGs_between_disease_stages.ipynb** in Jupyter notebook. This script generates DEGs between different disease stages, for CD14 monocytes and effector CD8 T cells respectively. It also present the expression heatmap of DEGs in PNG format file.
