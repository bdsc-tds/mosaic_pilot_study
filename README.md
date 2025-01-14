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
The computational environment can be recreated by installing the singularity container from the .def file. The R version we use for this workflow is 4.3.2, Bioconductor version 3.18, and renv is used to track specific versions of packages. Please find files related to renv in `/reproducibility/r/metadata`, and use `r.def` in `/reproducibility/r` to build the corresponding container:

``` bash
# Note: the current working directory is the root of this repo
cd reproducibility/r
singularity build --fakeroot --force /path/to/the/built/container r.def
```

The install time for this container on a regular computer or a computing cluster is approximately 2 hours with 8GB memory size. Launch the container from your terminal.

``` bash
singularity shell --bind /scratch,/users,/work /path/to/container/environment.sif
```

For Visium deconvolution with cell2location, a separate conda environment was used and can be installed with: 
```bash
conda env create -f /Owkin/C2L_code/env.yml
```

Activate the deconvolution environment:
```bash
conda activate cell2loc
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

The analysis for each technology should be run in sequential order based on the naming of files. For bash files, only scripts titled main_.sh should be submitted as jobs to high computing clusters. Note that path to container should be modified in _.sh files. The computing time of each analysis step follows its method. For instance, as detailed in cell2location deconvolution, "a Visium sample with 4039 locations and 10241 genes should take about 17-40 minutes, depending on the computer's GPU hardware".

For visualization, the following code links are mapped to the creation of each figure/supplement figure. 

# Visualization
## Figures
| Script                         | Figures         |
|--------------------------------|-----------------|
|  Data Characteristics ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig1*)) | Fig. 1 |
|  Visium GeoMx Deconvolution Specificity ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig2*)) | Fig. 2 |
|  Visium GeoMx Registration ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig3_12_Vis_Geo_Mapped)) | Fig. 3 |
|  Visium Annotation Gallery ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig4_Vis_patho_decon_gallery); [Revision](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/revision/Fig4_Vis)) | Fig. 4 a,b, g-i |
|  GeoMx Annotation Gallery ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/revision/Fig4_Geo)) | Fig. 4 c,d, g-i |
|  Visium GeoMx Pathology Deconvolution Heatmap ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/revision/Fig4_patho_decon_heatmap)) | Fig. 4 e,f |
|  TLS in L1 Visium ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig5_Visium_clustering_biology/L1)) | Fig. 5 a |
|  TLS in L1 GeoMx ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig5_Geo/Geo_L1_Overlay.R)) | Fig. 5 a |
|  Intra-patient Heterogeneity in B3 Visium ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig5_Visium_clustering_biology/B3); [Revision](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig5_B3_DE_Pathways/volcano_final/vis_volcano_B3.R)) | Fig. 5 b,d |
|  Intra-patient Heterogeneity in B3 GeoMx ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig5_Geo/Fig5_B3_Decon_Pie.R); [Revision](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig5_B3_DE_Pathways/volcano_final/geo_volcano_B3.R)) | Fig. 5 c,e |
|  UpSet Plot ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig5_B3_DE_Pathways/UpSet_plot.R)) | Fig. 5 f |
|  Inter-patient heterogeneity across all technologies ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/Fig6_Three_Tech_Dotplot)) | Fig. 6|

## Supplementary Figures
| Script                         | Supplementary Figures         |
|--------------------------------|-------------------------------|
|  Data Characteristics ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Descriptive)) | Fig. S1 |
|  DV200 Blockage ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/DV200_nExpressedGene.R)) | Fig. S1 b-e |
|  Disease Biology Frequency ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Descriptive)) | Fig. S2 |
|  Chromium UMAP ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Chrom_pt_spec_tumor_marker_dotplot)) | Fig. S3 |
|  Chromium Cell Type Frequency ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Descriptive)) | Fig. S4 |
|  Chromium Dotplot ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Chrom_pt_spec_tumor_marker_dotplot)) | Fig. S4 dotplot |
|  GeoMx Descriptive ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/GeoMx_marker_exp_heatmap))) | Fig. S5 |
|  Visium Sample Gallery ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Visium_Sample_Gallery)) | Fig. S6 b,c; Fig. S7 a,b |
|  Visium GeoMx Immune Abundance ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/GeoMx_Visium_Immune_RedDim)) | Fig. S7 c,d|
|  Visium Pathology Deconvolution Agreement ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig4_patho_decon_heatmap_supp/Vis_Heatmap_per_Patho_Decon_avgfraction_final_level4.R)) | Fig. S8 a|
|  GeoMx Pathology Deconvolution Agreement ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig4_patho_decon_heatmap_supp/Geo_Heatmap_per_AOI_Decon_avgfraction_level4.R)) | Fig. S8 b|
|  L1 GeoMx ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig5_Geo/Geo_L1_Overlay.R)) | Figs. S9 a |
|  L4 GeoMx ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/Fig5_Geo_supp/Geo_L4_Overlay_Final.R)) | Figs. S9 b-c |
|  Visium Chromium Intra-patient heterogeneity ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/B3_Chrom_DE)) | Figs. S9 d-f |
|  Visium UMAP by Sample ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/SuppFig/Visium_Integration/Visium_Integration_UMAPs_sample.R)) | Figs. S10 |
|  Visium UMAP by Biology ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/SuppFig/Visium_Integration/Visium_Integration_UMAPs_patho_decon.R)) | Figs. S11 |
|  GeoMx TSNE Batch Effect ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/10_manuscript_figure_helper.R)) | Figs. S12 |
|  Visium Integrated Annotation ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/tree/main/CHUV/Manuscript_Figure/SuppFig/Visium_Integration)) | Fig. S13 |
|  Ranked Gene Selection ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/FC_healthyvstumor_3techs.R)) | Fig. S14 |
|  Visium Ridge Plot ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/Fig6_Three_Tech_Dotplot/visium_prep_level1_5_level4_pt_spec.R)) | Fig. S15 a,b|
|  GeoMx Ridge Plot ([Code](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/CHUV/Manuscript_Figure/Fig6_Three_Tech_Dotplot/Geo_dotplot.R)) | Fig. S15 c,d|
|  Decon Assisted Inter-patient heterogeneity across all technologies ([Revision](https://github.com/bdsc-tds/mosaic_pilot_study/blob/main/revision/new_Supp_Fig6.R)) | Fig. S15 e|