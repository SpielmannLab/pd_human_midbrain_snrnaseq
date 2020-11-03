# PD variants enrichment on cell-type marker genes using MAGMA

Scripts to perform MAGMA analysis for various genetic and transcriptomics projects.

MAGMA version 1.07 downloaded from https://ctg.cncr.nl/software/magma

MAGMA is a gene set enrichment analysis method, which tests the joint association of all SNPs in a gene with the phenotype, while accounting for LD structure between SNPs. Competitive gene set analysis was performed on SNP p-values from the GWAS summary statistics of Nalls et al. (excluding 23andMe) and the publicly available European subset of 1000 Genomes Phase 3 was used as a reference panel to estimate LD between SNPs.
We used MAGMA here to assess the association of cell type-specific expressed genes with genetic risk of Parkinson's disease (PD), in order to identify disease-relevant cell types in the midbrain.

Preparation of MAGMA input files:

MAGMA input files are all in text formats. It requires:
- PD GWAS summary statistics file: ‘nallsEtAl2019_excluding23andMe_allVariants.tab’ downloaded from https://drive.google.com/file/d/1FZ9UL99LAqyWnyNBxxlx6qOUlfAnublN/view?usp=sharing and including an extra column 'N_all' for sample size per SNP.
- LOC_SNPs file: contains SNP information extracted from GWAS summary statistics. It should contain three columns: SNP ID, chromosome and base pair position.
- gene.loc file: the gene location file (NCBI GRCh37 build) download from https://ctg.cncr.nl/software/magma
- g1000_eur: 1000 Genomes data in plink format used as reference for LD
- gene_set file: with each row corresponding to a gene set: name of the gene set followed by the gene IDs, separated by whitespace. (see example genes_set_input_file)

MAGMA analysis:

A basic analysis in MAGMA consists of three steps: first, an annotation step to map SNPs onto genes; second, a gene analysis step to compute gene p-values; and three, a gene-level analysis step: either a generalized gene-set analysis, a gene property analysis, or both. (Please see MAGMA manual_v1.07b.pdf available in https://ctg.cncr.nl/software/magma for more details)

1. The annotation step is a pre-processing step prior to the actual analysis, in which SNPs are mapped to genes. The mapping is based on genomic location, assigning a SNP to a gene if the SNP’s location falls inside the region provided for each gene. Gene boundaries were defined as the transcribed region of each gene. An extended window of 35 kb upstream and 10 kb downstream of each gene was added to the gene boundaries. The annotation step will produce an output file with the .genes.annot suffix.
2. In the gene analysis step the gene p-values are computed. The gene analysis results plus gene correlations are output into a formatted output file with .genes.raw suffix. Here, we are using the default gene analysis model (snp-wise=mean) which uses the sum of -log(SNP p-value) as test statistic. The gene analysis step will produce an output file with .genes.raw suffix.
3. The gene-level (Gene-set) analysis: The competitive gene-set analysis is implemented as a linear regression model on a gene-level data matrix.
    - Primary output for the gene-level analysis is written to a .gsa.out file (see example gene_sets_analaysis_output.gsa.out)
    - Per-gene output for every gene-set is written to a .gsa.sets.genes.out file (see example gene_sets_analaysis_output.gsa.sets.genes.out)
    - Per-gene output for all genes in the analysis is written to a .gsa.genes.out file (see example gene_sets_analaysis_output.gsa.genes.out)

Example: 
Evaluation of enrichment of PD-associated variants with cell-type-specific differentially expressed genes (DEG).
- Input files: (genes_set_input_file) 11 genes-sets corresponding to the 11 cell-type.
- Output files:
    - gene_sets_analaysis_output.gsa.out : statistical results for every gene-sets.
    - gene_sets_analaysis_output.gsa.sets.genes.out : statistical results for every gene in gene-sets.
    - gene_sets_analaysis_output.gsa.genes.out : statistical results for all genes used in the analysis.

