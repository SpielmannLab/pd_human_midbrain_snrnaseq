#
# Estimate expression difference given a colData variable and a clustering-resolution or a cell-ontology (Seurat3 objects only)
#

" Estimate expression difference given a colData variable and a clustering-resolution or a cell-ontology (Seurat3 objects only)

Usage: diff_exp_seurat2monocle.R --jobname=<value> [--mcafile=<file>] --groupvar=<value> --level2test=<string> [--samplevar=<value>] --resolution=<value> [--distribution=<string>] [--det_thr=<value>] [--ncores=<value>] --specie=<value> --infolder=<folder> --outfolder=<folder>

Options:
  -h --help              Show this screen.
  --version              00.99.01
  --jobname=<file>       Descriptive name for your experiment.
  --mcafile=<file>       rds file with mca (Seurat3) object.
  --groupvar=<value>     Variable to do comparison.
  --level2test=<string>  groupvar level to calculate the differential gene expression effect-size.
  --samplevar=<value>    Sample variable. Ref for pseudo-bulk.
  --resolution=<value>   Clustering resolution.
  --distribution=<value> Distribution to use for diff. exp. e.g. quasipoisson, negbinomial, poisson, binomial, gaussian, zipoisson, or zinegbinomial.
  --det_thr=<value>      Minimum percentage of expressing cells of a given gene to be considered.
  --ncores=<value>       Number of cores to be used.  
  --infolder=<file>      Path to the single_cell_data .rds files.  
  --outfolder=<file>     Path to results folder.

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)
	
# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
	  'Rtsne', 'MASS', 'htmlwidgets', 'monocle3')

if (!require("BiocManager", character.only = TRUE)) {
	install.packages("BiocManager")
	BiocManager::install()
} else {
	ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
	if (any(!ipkgs)) {
		BiocManager::install(pkgs[!ipkkgs])
	} else {
		message("\n\nCool! your machine has everything is needed.\n\n")
	}
}

# Functions
seurat2monocleDE <- function(mca, jobname, groupvar, level2test, det_thr, distribution, ncores) {
	
	## Monocle diff. exp. framework
	cds <- new_cell_data_set(expression_data = mca@assays$RNA@counts,
			 cell_metadata = mca@meta.data,
			 gene_metadata = row_df[rownames(mca@assays$RNA@counts), ])

	cds <- preprocess_cds(cds, method = "PCA", 
		      num_dim = 25,
		      norm_method = "log", 
		      use_genes = NULL,
		      scaling = TRUE,  
		      verbose = TRUE, 
		      cores = 10)

	print(dim(cds))

	# Sub-setting genes that are detected in at least det_thr % of the cells of at least one var level
	print(table(mca@meta.data[[groupvar]]))
	
	genes <- lapply(unique(mca@meta.data[[groupvar]]), function(condition) {
	
		       scds <- cds[, which(mca@meta.data[[groupvar]] == condition)]
		       scds <- scds[which(rowSums(exprs(scds) > 0) > det_thr*ncol(exprs(scds))), ]
		       rownames(scds)
	
		  })

	igenes <- unique(unlist(genes))
	print(head(igenes))
	print(length(igenes))

	sub_cds <- cds[igenes, ]

	cgene_fits <- fit_models(sub_cds,
				 model_formula_str = paste0("~", groupvar),
				 expression_family = distribution,
				 cores = ncores)

	fit_coefs <- coefficient_table(cgene_fits)

	colnames(fit_coefs)
	table(fit_coefs$term)
	unique(fit_coefs$term)

	coef_res <- fit_coefs %>% 
	       filter(term == paste0(groupvar, level2test)) %>% 
	       filter (q_value < 0.05) %>%
	       select(gene, gene_short_name, term, q_value, estimate) %>% 
	       arrange(q_value) %>%
	       data.frame

	print(head(coef_res, 10))

	write.table(coef_res, paste0(outputfolder, "/monocle3_diff_exp_", jobname, "_qp.tsv"),
		   sep = "\t", quote = FALSE, row.names = FALSE)

	pdf(paste0(outputfolder, "/monocle3_diff_exp_", jobname, "_qp.pdf"))

	plot(

		plot_genes_violin(sub_cds[head(arrange(coef_res, q_value),10)[["gene"]], ], 
				 group_cells_by="condition", ncol=5) 

		)

	dev.off()

}

# --- Input parameters

jobname <- arguments$jobname
groupvar <- arguments$groupvar
samplevar <- arguments$samplevar
res <- arguments$resolution
distribution <- arguments$distribution
specie <- arguments$specie
inputfolder <- arguments$infolder
outputfolder <- arguments$outfolder
det_thr <- arguments$det_thr
ncores <- arguments$ncores
level2test <- arguments$level2test

# --- Read cell_data_set object with merged samples

if(!is.null(arguments$mcafile)) {

	mca_file <- paste0(inputfolder, '/', arguments$mcafile)
	jobname <- gsub("\\.rds", "", arguments$mcafile)
	
} else {

	mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')
	
}

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# --- Setting default cell identity (cell-type annotation)
# --- Check if numeric (resolution) or character (cell-ontology-variable)

if(is.na(as.numeric(`resolution`))) {

	cluster_var <- `resolution`
	print(cluster_var)
	print(table(mca@meta.data[[cluster_var]]))
	Idents(mca) <- cluster_var
	
	} else {
	
	cluster_var <- colnames(mca@meta.data)[grep(paste0("snn_res.", `resolution`), colnames(mca@meta.data))][1]
	print(cluster_var)
	print(table(mca@meta.data[[cluster_var]]))
	Idents(mca) <- cluster_var

}

# Checking cell metadata variable
print(table(mca@meta.data[[groupvar]]))

# Checking sample metadata
print(table(mca@meta.data[[samplevar]]))

# Gene annotation
if(specie == "hg") {

	refbiomart <- data.table::fread("./data/ensg_biomart_20191206.tsv")

} else if(specie == "mm"){

	refbiomart <- data.table::fread("./data/biomart_mm_20191206.tsv")

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}


# Counts
raw_eset <- mca@assays$RNA@counts 
print(raw_eset[1:5, 1:10])
print(dim(raw_eset))

row_df <- mca@assays$RNA@meta.features
row_df[["gene"]] <- rownames(row_df)
print(dim(row_df))
row_df <- merge(row_df, refbiomart, by.x ="gene", by.y = "Gene stable ID", all.x = TRUE) 
row_df <- filter(row_df, !duplicated(`gene`))
print(dim(row_df))
colnames(row_df)[7] <- "gene_short_name"
rownames(row_df) <- row_df[["gene"]]
print(head(row_df, 2))


message("Single-cell gene differential expression modelling.")
lapply(unique(mca@meta.data[[ident_col]]), function(cell) {

	       print(cell)

	       cells <- rownames(mca@meta.data)[which(mca@meta.data[[ident_col]] == cell)]
	       
	       sub_mca <- mca[, cells]
	       
	       seurat2monocleDE(mca=sub_mca, 
	       			jobaname=paste0(gsub("\\+| ", "_", cell), '_', "DEG"),
	       			groupvar=groupvar,
	       			leveltest=level2test,
	       			det_thr=det_thr,
	       			distribution=distribution,
	       			ncores=ncores)
})

message("Single-cell gene differential expression modelling is done *.*")
