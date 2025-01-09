# MOSAIC Pilot Study
<img src="mosaic.png" width="150" align="right"/>

This repository hosts the code to reproduce the analysis steps for all technologies across all indications, including the visualization code for creating Figure and Supplementary Figures in our manuscript, titled 

<b> Transcriptome Analysis of Archived Tumor Tissues by Visium, GeoMx DSP, and Chromium Methods Reveals Inter- and Intra-Patient Heterogeneity.</b> 

# Installation
To start with replicating our result, clone the repository to your working directory.

``` bash
git clone https://github.com/bdsc-tds/mosaic_pilot_study.git
```

# Singularity container
The computational environment can be recreated by installing the singularity container from the .def file. The R version we use for this workflow is 4.3.2, and renv is used to track specific versions of packages. Please find files related to renv in reproducibility/r/metadata, and use r.def in reproducibility/r to build the corresponding container:

``` bash
# Note: the current working directory is the root of this repo
cd reproducibility/r
singularity build --fakeroot --force /path/to/the/built/container r.def
```

Launch the container from your terminal.

``` bash
singularity shell --bind /scratch,/users,/work /path/to/container/environment.sif
```

# Analysis
The folder structure of this repository is detailed as follow
```
    CHUV
        ├── Chromium
        ├── Visium
        ├── GeoMx
        └── Manuscript_Figure
            ├── Figure 1 - 6
            └── SuppFig
```

The analysis for each technology should be run in sequential order based on the naming of files. For bash files, only scripts titled main_.sh should be submitted as jobs to high computing clusters. Note that path to container should be modified in _.sh files. 

For visualization, the following code links are mapped to the creation of each figure/supplement figure. 

# Visualization
## Figures
| Script                         | Figures         |
|--------------------------------|-----------------|
|  Data Characteristics ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig1*)) | Fig. 1 |
|  Visium GeoMx Deconvolution Specificity ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig2*)) | Fig. 2 |
|  Visium GeoMx Registration ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig3_12_Vis_Geo_Mapped)) | Fig. 3 |
|  Visium Annotation Gallery ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig4_Vis_patho_decon_gallery)) | Fig. 4 |
|  GeoMx Annotation Gallery ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig4_Vis_patho_decon_gallery)) | Fig. 4 |
|  Intra-patient heterogeneity in Visium ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig5_Visium_clustering_biology)) | Fig. 5 |
|  Intra-patient heterogeneity in GeoMx ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig5_Visium_clustering_biology)) | Fig. 5 |
|  Inter-patient heterogeneity across all technologies ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig6_Three_Tech_Dotplot)) | Fig. 6|

## Supplementary Figures
| Script                         | Supplementary Figures         |
|--------------------------------|-------------------------------|
|  Data Characteristics ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Descriptive)) | Fig. S1 |
|  Chromium UMAP ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Chrom_pt_spec_tumor_marker_dotplot)) | Fig. S2 |
|  Chromium Characteristics ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Descriptive)) | Fig. S3 |
|  Chromium Characteristics ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Chrom_pt_spec_tumor_marker_dotplot)) | Fig. S3 dotplot |
|  GeoMx Descriptive ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/GeoMx_marker_exp_heatmap))) | Fig. S4 |
|  Visium Sample Gallery ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Visium_Sample_Gallery)) | Fig. S5 |
|  Visium Pathology Deconvolution Agreement ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/SuppFig/Descriptive_Heatmap/Vis_Heatmap_per_Patho_Decon_avgfraction.R)) | Fig. S6|
|  Visium GeoMx Immune Abundance ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/GeoMx_Visium_Immune_RedDim)) | Fig. S7 a,b|
|  Visium Chromium Intra-patient heterogeneity ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/B3_Chrom_DE)) | Figs. S7 c-g |
|  Visium Integrated Annotation ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Visium_Integration)) | Fig. S8 |
|  Visium Ridge Plot ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/Fig6_Three_Tech_Dotplot/visium_prep_level1_5_level4_pt_spec.R)) | Fig. S9 a,b|
|  GeoMx Ridge Plot ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/Fig6_Three_Tech_Dotplot/Geo_dotplot.R)) | Fig. S9 c,d|
|  Decon Assisted Inter-patient heterogeneity across all technologies ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/Fig6_Three_Tech_Dotplot/e_Three_Tech_Dotplot_with_decon_chrom_level2_vis_level1_5_level4_pt_spec.R)) | Fig. S9 e|