#
# Cluster a seurat3 object using several resolution points
#

"  Cluster a seurat3 object using several resolution points

Usage: cluster_seurat3.R --jobname=<value> [--mcafile=<file>] --npcs=<value> [--assay=<value>] --testpcs=<logical> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     Optional. If !is.null, do not use jobname to guess mcafilename.
  --npcs=<value>       Number of PCs to consider for clustering.
  --assay=<value>      Assay to use for similarity network.
  --testpcs=<logical>  Wheather the mca input object contains JackStraw-test results.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.
  --ncores=<value>     Number of processors to use

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
	  'Rtsne', "monocle3")

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

# ---------- Functions ---------

getcluster <- function(mca, res) {

	mca <- FindClusters(mca, 
			    graph.name = NULL,
			    modularity.fxn = 1, 
			    initial.membership = NULL, 
			    weights = NULL,
			    node.sizes = NULL, 
			    resolution = res, 
			    algorithm = 1, 
			    n.start = 10,
			    n.iter = 10, 
			    random.seed = 0, 
			    group.singletons = TRUE,
			    temp.file.location = NULL, 
			    edge.file.name = NULL, 
			    verbose = TRUE)

	res_cl <- colnames(mca@meta.data)[grep("_snn_res.", colnames(mca@meta.data))]
	print(res_cl)

	g <- DimPlot(mca, 
		     reduction = "umap", 
		     pt.size = 0.01,
		     group.by = res_cl) + 
	  	ggplot2::theme(legend.position = "none") +
		ggtitle(res)
	
	print(tail(colnames(mca@meta.data)))

	clstr <- select(mca@meta.data, matches(paste0("_snn_res.", res)))
	print(head(clstr, 2))
	rownames(clstr) <- rownames(mca@meta.data)

	return(list('g' = g, 
		    'cluster' = clstr))
}


# ------------------------------

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
npcs <- as.numeric(arguments$npcs)
testpcs <- as.logical(arguments$testpcs)

if(!is.null(arguments$assay)) {
	`assay` <- arguments$assay
} else {
	`assay` <- "SCT"
} 


# Memory limits
options(future.globals.maxSize=91943040000000)
# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

# --- Read cell_data_set object with merged samples
if(!is.null(arguments$mcafile)) {
	mca_file <- paste0(inputfolder, '/', arguments$mcafile)
	jobname <- gsub("\\.rds", "", arguments$mcafile)
} else {
	mca_file <- paste0(inputfolder, '/', jobname, '_seurat3_merged_clustered.rds')
}

mca <- readRDS(mca_file)

print(mca)

# Remove legacy clustering metadata
new_meta <- mca@meta.data[-grep("_snn_res", colnames(mca@meta.data))]
rownames(new_meta) <- rownames(mca@meta.data)
mca@meta.data <- new_meta[colnames(mca), ]
#

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# Vis PCs contribution
if(testpcs) {

	# check informative PCs
	pdf(paste0(outputfolder, '/', jobname, "_jackStraw.pdf"))
		plot(
			JackStrawPlot(mca, dims=seq(npcs))
		)
	dev.off()

	message("########################################################################################")
	message(paste("Now check the", paste0(outputfolder, '/', jobname, "_jackStraw.pdf"), "file and", 
		      "decide an informative number of PCs to consider for clustering. Then, re-run this",
		      "function using --testpcs=FALSE."))
	message("########################################################################################")

} else {

	
	# --- Rounds of clustering

	DefaultAssay(mca) <- `assay`

	mca <- RunPCA(mca,
			  assay = `assay`,
			  npcs = npcs,
			  verbose = TRUE)

	mca <- RunUMAP(mca,
		       assay = `assay`,
		       dims = seq(npcs),
		       reduction = "pca",
		       n.neighbors = 30,
		       metric = "cosine",
		       min.dist = 0.3, spread = 1, 
		       set.op.mix.ratio = 1,
		       local.connectivity = 1L, 
		       repulsion.strength = 1)

	mca <- FindNeighbors(mca,
			     assay = `assay`,
			     dims = seq(npcs), 
			     verbose = TRUE)


	res <- c(0.01, 0.02, 0.04, 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.7, 1, 1.5, 2, 3)
	
	cls <- mclapply(res, function(r) getcluster(mca, r), mc.cores = ncores)
 
	cluster_jobname <- paste0(jobname, '_', ncol(mca), '_cells_and_', npcs, 'PCs_from_', nrow(mca), "_genes")  

	# --- Write QC plots
	title <- ggdraw() + 
	  draw_label(
	    cluster_jobname,
	    fontface = 'bold',
	    hjust = 0.5
	  ) +
	  theme(
	    plot.margin = margin(0, 0, 0, 0)
	  )

	rfile <- paste0(outputfolder, '/', cluster_jobname, '_', 'clustering.pdf')

	print(rfile)

	cl_gplots <- lapply(cls, function(cl) cl[['g']])

	pdf(rfile, height = 11.69, width = 8.27) 
	    plot(
		 plot_grid(title, plot_grid(plotlist = cl_gplots, ncol = 2), 
		      ncol = 1, rel_heights = c(0.02, 1))
		 )
	dev.off()

	jpeg(paste0(outputfolder, '/', cluster_jobname, "_", "clustering.jpeg"), 
			   height = 1980, width = 2980, pointsize = 74, quality = 100)
	    plot(
	    	plot_grid(title, plot_grid(plotlist = cl_gplots, ncol = 4), 
		      ncol = 1, rel_heights = c(0.02, 1))
	    )
	dev.off()

	meta_cls <- Reduce(cbind, lapply(cls, function(cl) cl[['cluster']]))

	new_meta <- cbind(mca@meta.data, meta_cls[rownames(mca@meta.data), ])
	rownames(new_meta) <- rownames(mca@meta.data)
	
	mca@meta.data <- new_meta[colnames(mca), ]
	print(head(mca@meta.data,2))

	mca_file <- paste0(outputfolder, '/', jobname, '_', `assay`, '_seurat3_clustered_resrange.rds')
	print(mca)
	saveRDS(mca, file = mca_file)
	message(paste("The cell_data_set object (Seurat3)", 
		      mca_file, "has been created."))
}
