# Single-cell sequencing of the human midbrain reveals glial activation and a neuronal state specific to Parkinson's disease

Data analysis workflow to reproduce the findings of the manuscript:  "Single-cell sequencing of the human midbrain reveals glial activation and a neuronal state-specific to Parkinson's disease". Raw-data is publically available in the Gene Expression Omnibus (GEO) with the accession number [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783). 

## General data analysis workflow

### Single-cell RNA sequencing

We investigated 11 human midbrain sections using the 10X scRNAseq solution: 5 IPD patients and 6 sex- and age-matched controls. The CellRanger 3.0 pipeline was used to count the transcripts. Scrublet was used to identify and filter out dublets. Seurat3 and monocle3 R (version > 4.0) packages were used for data normalization, sample integration, clustering, and trajectory inference. Individual functions for each of these steps are called from bash workflow scripts. 

#### Reads mapping and UMI count

UMI count matrix and cell metadata are available as supplementary files in the GEO accession number [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783). Reads were aligned on the 10X pre-indexed reference human genome (hg19, GRCh38). 

#### Preprocessing, duplet-scoring, cell-filtering, sample normalization & integration, clustering & cell-type annotation

`./scrnaseq/run_analysis_midbrain.sh`

### Genetic enrichment



### Image analysis
