#
# Estimate cell trajectory for every single cluster. Input seurat3 objects only. Monocle3 for trajectory inference 
#

" Estimate cell trajectory for every single cluster. Input seurat3 objects only. Monocle3 for trajectory inference 

Usage: subtrajectory_seurat2monocle.R --jobname=<value> [--mcafile=<file>] --groupvar=<value> [--samplevar=<value>] [--alignvar=<string>] --resolution=<value> [--nhvg=<value>] --npcs=<value> --specie=<value> --ncores=<value> --infolder=<folder> --outfolder=<folder>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     rds file with mca (Seurat3) object.
  --groupvar=<value>   Variable to do comparison. eg. condition: IPD, CO.
  --samplevar=<value>  Sample variable. Ref for pseudo-bulk.
  --alignvar=<string>  Variable o regress out before trajectory inference.
  --resolution=<value> Clustering resolution.
  --nhvg=<value>       Number of highly variable genes to be considered for trajectory reconstruction.
  --npcs=<value>       Number of PCs to consider for trajectory reconstruction.
  --ncores=<value>     Number of cores to be used.
  --infolder=<file>    Path to the single_cell_data .rds files.
  --outfolder=<file>   Path to results folder.

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- 
	
# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
	  'Rtsne', 'MASS', 'htmlwidgets', 'monocle3', '')

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

jobname <- arguments$jobname
groupvar <- arguments$groupvar
samplevar <- arguments$samplevar
resolution <- arguments$resolution
alignvar <- arguments$alignvar
nhvg <- arguments$nhvg
npcs <- arguments$npcs
ncores <- arguments$ncores
specie <- arguments$specie
inputfolder <- arguments$infolder
outputfolder <- arguments$outfolder

# --- Read cell_data_set 
if(!is.null(arguments$mcafile)) {

	mca_file <- paste0(inputfolder, '/', arguments$mcafile)
	jobname <- gsub("\\.rds", "", arguments$mcafile)
	
} else {

	mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')
	
}

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# Setting default cell identity based on clustering resolution

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

# gene metadata
row_df <- mca@assays$RNA@meta.features
row_df[["gene"]] <- rownames(row_df)
print(dim(row_df))
row_df <- merge(row_df, refbiomart, by.x ="gene", by.y = "Gene stable ID", all.x = TRUE) 
row_df <- filter(row_df, !duplicated(`gene`))
print(dim(row_df))
colnames(row_df)[7] <- "gene_short_name"
rownames(row_df) <- row_df[["gene"]]
print(head(row_df))

print(unique(mca@meta.data[[cluster_var]]))

lapply(unique(mca@meta.data[[cluster_var]]), function(cell) {

	       # Sub-setting dataset
	       cells <- rownames(mca@meta.data)[which(mca@meta.data[[cluster_var]] == cell)]
	       sub_mca <- mca[, cells]
	       DefaultAssay(sub_mca) <- "RNA"
	       
	       # Highly variable genes (hvg)
	       sub_mca <- FindVariableFeatures(sub_mca,
				    selection.method = "vst",
				    nfeatures = nhvg)
	       hvg <- sub_mca@assays$RNA@var.features
	       
	       # Trajectory inference
	       seurat2monocleTR(mca=sub_mca, 
	       			jobname=paste0(gsub("\\+| ", "_", cell), '_', "trajectory"), 
	       			ncores=ncores,
	       			npcs=npcs,
	       			hvg=hvg, 
	       			alignvar=alignvar,
	       			det_thr=0.10)

})

seurat2monocleTR <- function(mca, jobname, ncores=1, npcs=25, hvg=NULL, alignvar=NULL, det_thr=0.10) {
	
	## Monocle diff. exp. framework
	cds <- new_cell_data_set(expression_data = mca@assays$RNA@counts,
			 cell_metadata = mca@meta.data,
			 gene_metadata = row_df[rownames(mca@assays$RNA@counts), ])

	cds <- preprocess_cds(cds, method = "PCA", 
		      num_dim = npcs,
		      norm_method = "log", 
		      use_genes = hvg,
		      scaling = TRUE,  
		      verbose = TRUE, 
		      cores = ncores)

	if(!is.null(alignvar)) {

		cds <- align_cds(cds, 
				 alignment_group = alignvar)

		cds <- reduce_dimension(cds,
					reduction_method = "UMAP")

	} else {

		cds <- reduce_dimension(cds,
			       reduction_method = "UMAP", 
			       preprocess_method = "PCA")
	
	}

       cds <- cluster_cells(cds, 
			    preprocess_method = "PCA")

       print(table(cds@clusters$UMAP$partitions))

       cds <- learn_graph(cds,
			  use_partition = TRUE)

       if(grepl("Microglia", jobname)) {

	       cds <- order_cells(cds,
				  reduction_method = "UMAP",
				  "Y_25"		  
	       )

       } else {
	       
	       cds <- order_cells(cds,
				  reduction_method = "UMAP",
				  get_earliest_principal_node(cds)		  
		)

       }

       pdf(paste0(outputfolder, '/', jobname, "_", "_monocle3_trajectory.pdf"))

       plot(
	       plot_cells(cds,
			  color_cells_by = samplevar,
			  label_groups_by_cluster=TRUE,
			  label_leaves=TRUE,
			  label_branch_points=TRUE)
	       )
       plot(
	       plot_cells(cds,
			  color_cells_by = groupvar,
			  label_groups_by_cluster=FALSE,
			  label_leaves=FALSE,
			  label_branch_points=FALSE)
	       )
       plot(
	       plot_cells(cds,
			  color_cells_by = "pseudotime",
			  label_groups_by_cluster=FALSE,
			  label_leaves=FALSE,
			  label_branch_points=FALSE)
	       )
      dev.off()

	# Differential expression accross the trajectory
	genes <- lapply(unique(pData(cds)[[groupvar]]), function(condition) {

			      scds <- cds[, which(pData(cds)[[groupvar]] == condition)]
			      scds <- scds[which(rowSums(exprs(scds) > 0) > det_thr*ncol(exprs(scds))), ]
			      rownames(scds)
	       })

        igenes <- unique(unlist(genes))

        igenes <- igenes[which(igenes %in% hvg)]

	pseudotime_genes <- graph_test(cds[igenes, ], 
				       neighbor_graph="principal_graph", 
				       cores=ncores)

	pseudotime_genes <- pseudotime_genes %>%
		filter(q_value < 0.05) %>%
		arrange(q_value)

	write.table(pseudotime_genes, paste0(outputfolder, "/monocle3_diff_exp_", jobname, "_pseudotime.tsv"),
		   sep = "\t", quote = FALSE, row.names = FALSE)

       pdf(paste0(outputfolder, '/', jobname, "_", "_trajectory_genes_pseudotime.pdf"), height = 50)
	       plot(
	       
	       	plot_genes_in_pseudotime(cds[rowData(cds)$gene %in% head(arrange(pseudotime_genes, q_value), 10)[["gene"]], ],
                        color_cells_by="condition",
			min_expr=0.5)
	       )
       dev.off()

       saveRDS(cds, paste0(outputfolder, '/', jobname, "_", "monocle3_trajectory_pseudotime.rds"))

}

