#
# Read an integrated mca object (Seurat3 - active assay: integrated) and return the same object with normalized RNA and SCT assays
#
"  Read an integrated mca object (Seurat3 - active assay: integrated) and return the same object with normalized RNA and SCT assays

Usage: merge_sc_seurat3_step2.R --jobname=<value> [--mcafile=<file>] --npcs=<value> [--nhvg=<value>] [--nneigh=<value>] --colgroup=<value> --ncores=<value> --testpcs=<logical> --infolder=<folder> --outfolder=<folder> [--TMP=<logical>]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     Seurat3 object containing several samples integrated with merge_sc_seurat3_step1_5.R.
  --npcs=<value>       Number of principal components to use.
  --nhvg=<value>       Number of highly variable genes to consider.
  --nneigh=<value>     Number of nearest neighbors for umap estimation.
  --colgroup=<value>   Cell metadata variable to color cells.
  --ncores=<value>     Number of threads to use.
  --testpcs=<logical>  Wheather to test pc explanation.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.

"-> doc


suppressMessages(library(docopt))
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
	  'Rtsne', "phateR", "Rmagic", "monocle3", "")

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

# --------------------------------
# Functions
# --------------------------------

mca2tsne <- function(mca, npcs, colgroup, ncores) {

	pca <- Embeddings(object = mca@reductions$pca)

	tsne_res <- Rtsne::Rtsne(pca[, seq(npcs)], 
			 dims = 2, pca = FALSE,
			 check_duplicates = FALSE, 
			 num_threads = ncores)

	tsne_data <- data.frame(tsne_res$Y)
	row.names(tsne_data) <- rownames(pca)
	colnames(tsne_data) <- c("tSNE1", "tSNE2")
	tsne_data$id <- rownames(tsne_data)
	pdata <- mca@meta.data
	pdata$id <- rownames(pdata)
	tsne_data <- merge(tsne_data, pdata, by = "id")

	g_tsne <- ggplot(tsne_data, aes(x = tSNE1, 
					y = tSNE2, 
					color = !!rlang::sym(colgroup))) +
		    geom_point(size = 0.01) +
		    theme_classic() +
		    theme(legend.position="bottom")

	return(g_tsne)
}


mca2phate <- function(mca, `assay`, npcs, colgroup) {

	FM <- mca[[`assay`]]@scale.data # ngenes x cells

	phate_obj <- phate(t(FM),
			  ndim = 2, 
			  knn = 5, 
			  decay = 40, 
			  n.landmark = 500,
			  gamma = 1, 
			  t = "auto", 
			  knn.dist.method = "euclidean",
			  t.max = 100, 
			  npca = npcs,
			  n.jobs = 1 # no parallel here (...=1)
			  )

	phate_df <- data.frame(phate_obj$embedding)
	phate_df[['barcode']] <- rownames(phate_df) 

	cluster_df <- data.frame(mca@meta.data)
	cluster_df[['barcode']] <- rownames(cluster_df)

	phate_df <- merge(phate_df, cluster_df, by = "barcode")

	g_phate <- ggplot(data.frame(phate_df)) +
		geom_point(aes(PHATE1, PHATE2, 
			       color=!!rlang::sym(colgroup)), shape = ".") +
		theme_classic() +
		theme(legend.position="bottom")

	return(g_phate)
}

test_dim <- function(mca, `assay`, npcs) {

		g_elbow <- ElbowPlot(mca, 
				     ndims = npcs,
				     reduction = "pca")

		mca <- JackStraw(mca, 
				 reduction = "pca", 
				 assay = `assay`, 
				 dims = npcs,
				 num.replicate = 100, 
				 prop.freq = 0.01, 
				 verbose = TRUE,
				 maxit = 1000)

		mca <- ScoreJackStraw(mca, 
				      dims = seq(npcs), 
				      reduction = "pca",
				      score.thresh = 1e-05)
		
		g_js <- JackStrawPlot(mca, 
				      dims = seq(npcs),
)

		g_pcafeat <- DimHeatmap(mca, 
					dims = ceiling(seq(npcs)/4), 
					cells = 500, 
					balanced = TRUE,
					assays = `assay`,
					fast = FALSE)

		return(list('elbow'=g_elbow,
			    'js'=g_js,
			    'pcafeat'=g_pcafeat))
	}


plot_pca_umap <- function(mca, colgroup) {

	g_pca <- DimPlot(mca, 
			 reduction = "pca", 
			 pt.size = 0.01, 
			 group.by=colgroup) +
		ggplot2::theme(legend.position = "bottom")


	g_umap <- DimPlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  group.by=colgroup) +
		ggplot2::theme(legend.position = "bottom")

	return(list('pca'=g_pca,
		    'umap'=g_umap))
}


vis_meta <- function(mca) {
	metavar <- colnames(mca@meta.data)
	
	mt <- lapply(metavar, function(colgroup) {

		if(class(mca@meta.data[[colgroup]]) %in% c("integer", "numeric")) {

			mca@meta.data[[colgroup]] <- as.numeric(mca@meta.data[[colgroup]])
		
			print("Num:")
			print(colgroup)
			print(class(mca@meta.data[[colgroup]]))

			FeaturePlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  features=colgroup) +
		 		#ggplot2::theme(legend.position = "none") +
				ggplot2::theme(legend.position = "none") + 
				ggtitle(colgroup) +
				scale_colour_gradient(low = "grey", high = "blue")
	
		} else if(class(mca@meta.data[[colgroup]]) %in% c("character", "factor")) {
			
			mca@meta.data[[colgroup]] <- as.factor(mca@meta.data[[colgroup]])
			
			print("Fac:")
			print(colgroup)
			print(class(mca@meta.data[[colgroup]]))
	
			DimPlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  group.by=colgroup) +
		 		#ggplot2::theme(legend.position = "none") +
				ggplot2::theme(legend.position = "none") + 
				ggtitle(colgroup)
		}
	})

	return(mt)
}

norm_fs_sct <- function(mca, method, nhvg, npcs, nneigh, testpca, colgroup) {

	if(method == "effectsize") {

		DefaultAssay(mca) <- "RNA"

		mca <- NormalizeData(mca,
			     assay = "RNA",
			     normalization.method = "LogNormalize", 
			     scale.factor = 10000,
			     margin = 1, 
			     verbose = TRUE)

		mca <- FindVariableFeatures(mca,
				    assay = "RNA",
				    selection.method = "vst", 
				    nfeatures = nhvg, 
				    verbose = FALSE)

		mca <- ScaleData(mca,
			 features = mca@assays$RNA@var.features, 
			 assay = "RNA",
			 vars.to.regress = NULL, 
			 split.by = NULL, 
			 model.use = "linear",
			 use.umi = FALSE, 
			 do.scale = TRUE, 
			 do.center = TRUE,
			 scale.max = 10, 
			 block.size = 1000, 
			 min.cells.to.block = 3000,
			 verbose = TRUE)

		tmp_mca <- RunPCA(mca,
		      features = mca@assays$RNA@var.features,
		      assay = "RNA",
		      npcs = npcs, 
		      verbose = TRUE,
		      seed.use = 444)

		tmp_mca <- RunUMAP(tmp_mca,
		       assay = "RNA",
		       dims = seq(npcs),
		       reduction = "pca",
		       n.neighbors = nneigh,
		       metric = "cosine",
		       min.dist = 0.3, spread = 1, 
		       set.op.mix.ratio = 1,
		       local.connectivity = 1L, 
		       repulsion.strength = 1)

		pca_umap <- plot_pca_umap(tmp_mca, colgroup)

		tsne <- mca2tsne(tmp_mca, npcs, colgroup, ncores) 

		`phate` <- mca2phate(tmp_mca, "RNA", npcs, colgroup)

		mt <- vis_meta(tmp_mca)
		nc <- ceiling(sqrt(length(mt)))
	
		png(paste0(outputfolder, "/", jobname, "_integrated_npcs_", npcs, "_nvgs_", nhvg,
			   "_Seurat3_factor_size_non_corrected.png"),
		    width = 500*nc, height = 500*nc)
			plot(
			     plot_grid(plotlist = mt, ncol = nc)
			)
		dev.off()

		if(testpca) {

			dim_ggs <- test_dim(tmp_mca, "RNA", npcs)

			return(list("mca" = mca,
				    "gg" = list("pca_umap" = pca_umap,
						"tsne" = tsne,
						"phate" = `phate`),
				    "testdim" = dim_ggs))

		} else {

			return(list("mca" = mca,
				    "gg" = list("pca_umap" = pca_umap,
						"tsne" = tsne,
						"phate" = `phate`)))

		}

	} else if(method == "sct") {
	
		mca <- SCTransform(mca,
			   assay = "RNA",
			   vars.to.regress = NULL, 
			   variable.features.n = nhvg,
			   verbose = TRUE,
			   return.only.var.genes = TRUE,
			   seed.use = TRUE)

		tmp_mca <- RunPCA(mca,
		      features = mca@assays$SCT@var.features,
		      assay = "SCT",
		      npcs = npcs, 
		      verbose = TRUE,
		      seed.use = 444)

		tmp_mca <- RunUMAP(tmp_mca,
		       assay = "SCT",
		       dims = seq(npcs),
		       reduction = "pca",
		       n.neighbors = nneigh,
		       metric = "cosine",
		       min.dist = 0.3, 
		       spread = 1, 
		       set.op.mix.ratio = 1,
		       local.connectivity = 1L, 
		       repulsion.strength = 1)

		pca_umap <- plot_pca_umap(tmp_mca, colgroup)

		tsne <- mca2tsne(tmp_mca, npcs, colgroup, ncores) 

		`phate` <- mca2phate(tmp_mca, "RNA", npcs, colgroup)

		mt <- vis_meta(tmp_mca)

		nc <- ceiling(sqrt(length(mt)))

		png(paste0(outputfolder, "/", jobname, "_integrated_npcs_", npcs, "_nvgs_", nhvg,
			   "_Seurat3_SCT_non_corrected.png"),
		    width = 500*nc, height = 500*nc)
			plot(
			     plot_grid(plotlist = mt, ncol = nc)
			)
		dev.off()

		if(testpca) {

			dim_ggs <- test_dim(tmp_mca, "SCT", npcs)

			return(list("mca" = mca,
				    "gg" = list("pca_umap" = pca_umap,
						"tsne" = tsne,
						"phate" = `phate`),
				    "testdim" = dim_ggs))

		} else {

			return(list("mca" = mca,
				    "gg" = list("pca_umap" = pca_umap,
						"tsne" = tsne,
						"phate" = `phate`)))

		}
	}
}

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
colgroup <- arguments$colgroup
jobname <- arguments$jobname
npcs <- as.numeric(arguments$npcs)
testpcs <- as.logical(arguments$testpcs)

if(!is.null(arguments$nhvg)) {

	nhvg <- as.numeric(arguments$nhvg)

} else {

	nhvg <- 500

}

if(!is.null(arguments$mcafile)) {

	mca_file <- paste0(outputfolder, '/', arguments$mcafile)
	jobname <- gsub("\\.rds", "", arguments$mcafile)

} else {

	mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_merged.rds')

}

if(!is.null(arguments$nneigh)) {

	nneigh <- arguments$nneigh

} else {

	nneigh <- 30

}

# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

# ------------
# Reading merge mca file  
# ------------

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# Maximum vector size whithin the session env
options(future.globals.maxSize=41943040000000)

# assays
print(mca)
print(sapply(names(mca), function(a) dim(mca[[a]])))

# cell metadata
print(head(mca@meta.data,2))

# Run

# Effect size
method <- "effectsize"
rna_mca <- norm_fs_sct(mca, method, nhvg, npcs, nneigh, testpcs, colgroup)

# SCT
method <- "sct"
sct_mca <- norm_fs_sct(rna_mca[["mca"]], method, nhvg, npcs, nneigh, testpcs, colgroup)

# integrated
pca_umap <- plot_pca_umap(mca, colgroup)
tsne <- mca2tsne(mca, npcs, colgroup, ncores) 
`phate` <- mca2phate(mca, "integrated", npcs, colgroup)
int_mca <- list("gg" = list("pca_umap" = pca_umap,
		 "tsne" = tsne,
		 "phate" = `phate`))

rna_title <- ggdraw() + 
		draw_label("Factor size", fontface = 'bold', hjust = 0.5) +
		theme(plot.margin = margin(0, 0, 0, 0))

sct_title <- ggdraw() + 
		draw_label("SCT", fontface = 'bold', hjust = 0.5) +
		theme(plot.margin = margin(0, 0, 0, 0))

int_title <- ggdraw() + 
		draw_label("Integrated", fontface = 'bold', hjust = 0.5) +
		theme(plot.margin = margin(0, 0, 0, 0))


png(paste0(outputfolder, "/", jobname,
	   "_Seurat3_integration_sample_nhvg_", nhvg, "_npcs_", npcs, "_allassays.png"),
    width = 2500, height = 1500)

	plot(plot_grid(# rna
		       plot_grid(rna_title,
		      		plot_grid(rna_mca[["gg"]][["pca_umap"]][['pca']], rna_mca[["gg"]][['tsne']],
					rna_mca[["gg"]][["pca_umap"]][['umap']], rna_mca[["gg"]][['phate']], ncol = 4),
		      		ncol = 1, rel_heights = c(0.1, 1)),
		       # sct
		       plot_grid(sct_title,
		      		plot_grid(sct_mca[["gg"]][["pca_umap"]][['pca']], sct_mca[["gg"]][['tsne']],
					sct_mca[["gg"]][["pca_umap"]][['umap']], sct_mca[["gg"]][['phate']], ncol = 4),
		      		ncol = 1, rel_heights = c(0.1, 1)),
		       # integrated
		       plot_grid(int_title,
		      		plot_grid(int_mca[["gg"]][["pca_umap"]][['pca']], int_mca[["gg"]][['tsne']],
					int_mca[["gg"]][["pca_umap"]][['umap']], int_mca[["gg"]][['phate']], ncol = 4),
		      		ncol = 1, rel_heights = c(0.1, 1)),

		       ncol = 1)
	)

dev.off()

if(testpcs) {

	# Integrated
	int_ggs <- test_dim(mca, "integrated", npcs)

	png(paste0(outputfolder, "/", jobname,
		   "_Seurat3_integration_sample_nhvg_", nhvg, "_npcs_", npcs, "_allassays_dimensions.png"),
	    width = 2500, height = 1500)

		plot(plot_grid(# rna
		       plot_grid(rna_title,
		      		plot_grid(rna_mca[["testdim"]][["elbow"]], rna_mca[["testdim"]][['js']], ncol = 2),
		      		ncol = 1, rel_heights = c(0.1, 1)),
		       # sct
		       plot_grid(sct_title,
		      		plot_grid(sct_mca[["testdim"]][["elbow"]], sct_mca[["testdim"]][['js']], ncol = 2),
		      		ncol = 1, rel_heights = c(0.1, 1)),
		       # integrated
		       plot_grid(int_title,
		      		plot_grid(int_ggs[["elbow"]], int_ggs[["testdim"]][['js']], ncol = 2),
		      		ncol = 1, rel_heights = c(0.1, 1)),

		       ncol = 1)
		)
	dev.off()

}

# --- Write cell_data_set object with merged samples

mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_integrated_norm.rds')
out_mca <- sct_mca[["mca"]]
DefaultAssay(out_mca) <- "integrated"

print(out_mca)

saveRDS(out_mca, file = mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been created."))
