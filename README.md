# NPC_scRNA_ST_code
Data Analysis Process of "Single-cell and spatial transcriptomics reveal the molecular and histological mechanisms of radioresistance and immune escape in recurrent nasopharyngeal carcinoma".

For data access, please visit: https://ngdc.cncb.ac.cn/, the BioProject ID is PRJCA015229.

Below is the specific content introduction of the code file：

RNA-seq:

01: Deconv for bulk RNA-seq

02: Survival analysis

scRNA-seq:

01: Load scRNA-seq data data and preprocess

02: Integation by harmony and clustering/annotation for scRNA-seq data

03: NMF anaylsis for single malignant cells

04: DEG analysis

05: Crosstalk (ligand and receptor analysis) using liana

06.1: Trajectory analysis using Monocle3 for CAF (Cancer-associated fibroblast)

06.2: Trajectory analysis using Monocle3 for CD8T

07: CytoTRACE (Cellular (Cyto) Trajectory Reconstruction Analysis using gene Counts and Expression) analysis for CAF

Spatial transcriptomics:

01: Progeny (Pathway RespOnsive GENes for activity inference) for spatial data

02: Deconvolution for ST data by using SpaCET (Spatial Cellular Estimator for Tumors)

03: Colocalization score

04: Stemness scores for Malignant near/not near CAFs

05: Spatial interactions using COMMOT

06: Neighbor enrichment using Suqidpy

07: Distance from myeloid to malignant on ST data
