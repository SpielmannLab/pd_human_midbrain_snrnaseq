#
# Read and merge several single-cell (Seurat3) objects following the Seurat approach
#

" Read and merge several single-cell (Seurat3) objects following the Seurat approach

Usage: merge_sc_seurat3_step1.R --jobname=<value> --npcs=<value> [--nhvg=<value>] [--nneigh=<value>] [--method=<value>] --ncores=<value> --testpcs=<logical> [--specie=<value>] [--write=<logical>] [--covars=<value>] --infolder=<folder> --outfolder=<folder> [--explore=<logical>] [<file> <file>...]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --npcs=<value>       Number of principal components to use.
  --nhvg=<value>       Number of highly variable genes to consider.
  --nneigh=<value>     Number nearest neighbors for umap estimation.
  --method=<value>     Norm. methods. (SCT, or factor.size = RNA)
  --ncores=<value>     Number of threads to use.
  --testpcs=<logical>  Wheather PCs contribution to the global gene expression variance should be tested.
  --specie=<value>     Either mm or hg. Cell cycle gene reference.
  --write=<logical>    Wheather or not to write the mca preprocessed object.
  --covars=<value>     String separated by /// with cell metavariables to regress out.
  --infolder=<file>    Path to the single_cell_data .rds files.  
  --outfolder=<file>   Path to results folder.

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
	  'Rtsne', "phateR", "Rmagic", "monocle3")

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

pymagic <- reticulate::import("magic")

# ------- Functions

cds2mca <- function(cds) {
	mca <- CreateSeuratObject(counts = exprs(cds),
			  meta.data = data.frame(pData(cds)),
			  project = jobname)
	mca@assays$RNA@meta.features <- data.frame(fData(cds))
	mca
}


mca2tsne <- function(mca, npcs, colgroup, ncores) {

	pca <- Embeddings(object = mca@reductions$pca)

	print(head(pca))

	tsne_res <- Rtsne::Rtsne(pca[, seq(npcs)], 
			 dims = 2, pca = FALSE,
			 check_duplicates = FALSE, 
			 num_threads = ncores)

	print("tsne run!")

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
	
		g_pcafeat <- DimHeatmap(mca, 
					dims = ceiling(seq(npcs)/4), 
					cells = 500, 
					balanced = TRUE,
					assays = `assay`,
					fast = FALSE)

		if(`assay` != "SCT") {

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
				      dims = seq(npcs))

	
			return(list('elbow'=g_elbow,
				    'js'=g_js,
				    'pcafeat'=g_pcafeat))

		} else {

			return(list('elbow'=g_elbow,
				    'js'=NULL,
				    'pcafeat'=g_pcafeat))

		}
}

plot_pca_umap <- function(mca, colgroup) {

	print(head(mca@meta.data,2))

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


			if(length(which(is.na(mca@meta.data[[colgroup]]))) == nrow(mca@meta.data)) {

				mca@meta.data[[colgroup]][which(is.na(mca@meta.data[[colgroup]]))] <- 0
			
			}

			mca@meta.data[[colgroup]] <- as.numeric(mca@meta.data[[colgroup]])
		
			print("Num:")
			print(colgroup)
			print(class(mca@meta.data[[colgroup]]))

			FeaturePlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  features=colgroup) +
				ggplot2::theme(legend.position = "none") + 
				ggtitle(colgroup) +
				scale_colour_gradient(low = "grey", high = "blue")
	
		} else if(class(mca@meta.data[[colgroup]]) %in% c("character", "factor")) {

			if(length(which(is.na(mca@meta.data[[colgroup]]))) == nrow(mca@meta.data)) {

				mca@meta.data[[colgroup]][which(is.na(mca@meta.data[[colgroup]]))] <- "missing"
			
			}

			mca@meta.data[[colgroup]] <- as.factor(mca@meta.data[[colgroup]])
			
			print("Fac:")
			print(colgroup)
			print(class(mca@meta.data[[colgroup]]))
	
			DimPlot(mca, 
			  reduction = "umap", 
			  pt.size = 0.01, 
			  group.by=colgroup) +
				ggplot2::theme(legend.position = "none") + 
				ggtitle(colgroup)
		}
	})

	return(mt)
}

norm_fs_sct <- function(mca, name, method, nhvg, npcs, nneigh, testpca, colgroup, cc.markers, ncores, covars) {

	if(grepl("ENSG|ENSMUSG", rownames(mca)[1])) {

		   genename <- "Gene stable ID"
	} else {
		   genename <- "gene"
	}
		

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

		print("Cell cycle...")

		g2test <- rownames(mca)[1]
		print(g2test)

		if(grepl("\\.\\d+$", g2test)) {

			g_key <- rownames(mca)
			names(g_key) <- gsub("\\.\\d+$", "", g_key)
			print(head(g_key))

			s.feat <- g_key[names(g_key) %in% filter(cc.markers, phase == "s")[["Gene stable ID"]]]
			names(s.feat) <- NULL
			print(mca[s.feat, ])
			g2m.feat <- g_key[names(g_key) %in% filter(cc.markers, phase == "g2m")[["Gene stable ID"]]]
			names(g2m.feat) <- NULL
			print(mca[g2m.feat, ])

			mca <- CellCycleScoring(mca, 
					s.features = s.feat,
					g2m.features = g2m.feat,
					set.ident = FALSE)

			print("... done!")

		} else {

			print(filter(cc.markers, phase == "s")[[genename]])

			mca <- CellCycleScoring(mca, 
					s.features = filter(cc.markers, phase == "s")[[genename]],
					g2m.features = filter(cc.markers, phase == "g2m")[[genename]],
					set.ident = FALSE)

		}

		sex <- data.frame(t(mca[["RNA"]][grep("Xist|XIST|ENSMUSG00000086503|ENSG00000229807", rownames(mca[["RNA"]])), rownames(mca@meta.data)]))

		if(ncol(sex) == 0) {
			print("No Xist detected:!?")
			mca@meta.data[["sex"]] <- 0
		} else {
			mca@meta.data[["sex"]] <- sex[[1]]
		}

		hk <- data.frame(t(mca[["RNA"]][grep("Actb|ACTB|ENSMUSG00000029580|ENSG00000075624", rownames(mca[["RNA"]])), rownames(mca@meta.data)]))
		mca@meta.data[["housekeeping"]] <- hk[[1]]
		print(head(mca@meta.data,2))

		mca <- ScaleData(mca,
			 features = mca@assays$RNA@var.features, 
			 assay = "RNA",
			 vars.to.regress = covars, 
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
		       
		tmp_mca <- FindNeighbors(tmp_mca,
					dims = seq(npcs), 
					reduction = "pca",
					verbose = TRUE)

		tmp_mca <- FindClusters(tmp_mca,
					verbose = TRUE)
	
		print("Ploting...")
		pca_umap <- plot_pca_umap(tmp_mca, colgroup)

		print("Getting tSNE embedings...")
		tsne <- mca2tsne(tmp_mca, npcs, colgroup, ncores) 

		print("Getting phate embedings...")
		`phate` <- mca2phate(tmp_mca, "RNA", npcs, colgroup)

		mt <- vis_meta(tmp_mca)
		nc <- ceiling(sqrt(length(mt)))
	
		if(is.null(covars)) {

			gg_file <- paste0(outputfolder, "/", name, "_integrated_npcs_", npcs, "_nvgs_", nhvg, 
					  "_Seurat3_factor_size_non_corrected.png")

		} else {

			gg_file <- paste0(outputfolder, "/", name, "_integrated_npcs_", npcs, "_nvgs_", nhvg, 
					  "_Seurat3_factor_size_corrected_", 
					  paste(covars, collapse = "_"), ".png")

		}

		png(gg_file,
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

	} else if(method == "SCT") {

		print("Cell cycle...")

		g2test <- rownames(mca)[1]
		print(g2test)
		if(grepl("\\.\\d+$", g2test)) {

			g_key <- rownames(mca)
			names(g_key) <- gsub("\\.\\d+$", "", g_key)
			print(head(g_key))

			print(head(g_key[names(g_key) %in% filter(cc.markers, phase == "s")[["Gene stable ID"]]]))
			print(head(g_key[names(g_key) %in% filter(cc.markers, phase == "g2m")[["Gene stable ID"]]]))
			print(head(filter(cc.markers, phase == "g2m")[["Gene stable ID"]]))

			mca <- CellCycleScoring(mca, 
					s.features = g_key[names(g_key) %in% filter(cc.markers, phase == "s")[["Gene stable ID"]]],
					g2m.features = g_key[names(g_key) %in% filter(cc.markers, phase == "g2m")[["Gene stable ID"]]],
					set.ident = FALSE)

		} else {

			print(filter(cc.markers, phase == "s")[[genename]])

			mca <- CellCycleScoring(mca,
					s.features = filter(cc.markers, phase == "s")[[genename]],
					g2m.features = filter(cc.markers, phase == "g2m")[[genename]],
					set.ident = FALSE)

		}

		mca@meta.data[["cc.difference"]] <- mca@meta.data[["S.Score"]] - mca@meta.data[["G2M.Score"]]

		sex <- data.frame(t(mca[["RNA"]][grep("Xist|XIST|ENSMUSG00000086503|ENSG00000229807", rownames(mca[["RNA"]])), rownames(mca@meta.data)]))
		if(ncol(sex) == 0) {
			print("No Xist detected:!?")
			mca@meta.data[["sex"]] <- 0
		} else {
			mca@meta.data[["sex"]] <- sex[[1]]
		}

		hk <- data.frame(t(mca[["RNA"]][grep("Actb|ACTB|ENSMUSG00000029580|ENSG00000075624", rownames(mca[["RNA"]])), rownames(mca@meta.data)]))
		mca@meta.data[["housekeeping"]] <- hk[[1]]
		print(head(mca@meta.data,2))

		mca <- SCTransform(mca,
			   assay = "RNA",
			   vars.to.regress = covars, 
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
		       
		tmp_mca <- FindNeighbors(tmp_mca,
					dims = seq(npcs), 
					reduction = "pca",
					verbose = TRUE)

		tmp_mca <- FindClusters(tmp_mca,
					verbose = TRUE)

		pca_umap <- plot_pca_umap(tmp_mca, colgroup)

		tsne <- mca2tsne(tmp_mca, npcs, colgroup, ncores) 

		`phate` <- mca2phate(tmp_mca, "SCT", npcs, colgroup)

		mt <- vis_meta(tmp_mca)

		nc <- ceiling(sqrt(length(mt)))

		if(is.null(covars)) {

			gg_file <- paste0(outputfolder, "/", name, "_integrated_npcs_", npcs, "_nvgs_", nhvg, 
					  "_Seurat3_SCT_non_corrected.png")

		} else {

			gg_file <- paste0(outputfolder, "/", name, "_integrated_npcs_", npcs, "_nvgs_", nhvg, 
					  "_Seurat3_SCT_corrected_", 
					  paste(covars, collapse = "_"), ".png")

		}

		png(gg_file,
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

get_norm_single_mca <- function(mca, mca_name, nhvg, npcs, nneigh, testpcs, colgroup, cc.markers, ncores, w, covars, explore) {

	# assays
	print(mca_name)
	print(mca)

	# cell metadata
	print(head(mca@meta.data,2))
	# cell-cycle markers
	print(head(cc.markers,2))
	# Covariaties to regress
	print(covars)

	# Run
	# SCT
	method <- "SCT"
	sct_mca <- norm_fs_sct(mca, mca_name, method, nhvg, npcs, nneigh, testpcs, colgroup, cc.markers, ncores, covars)

	rna_title <- ggdraw() + 
		draw_label("Factor size", fontface = 'bold', hjust = 0.5) +
		theme(plot.margin = margin(0, 0, 0, 0))

	if(explore) {

		# Effect size
		method <- "effectsize"
		rna_mca <- norm_fs_sct(mca, mca_name, method, nhvg, npcs, nneigh, testpcs, colgroup, cc.markers, ncores, covars)

		sct_title <- ggdraw() + 
			draw_label("SCT", fontface = 'bold', hjust = 0.5) +
			theme(plot.margin = margin(0, 0, 0, 0))

		png(paste0(outputfolder, "/", mca_name,
			   "_Seurat3_nhvg_", nhvg, "_npcs_", npcs, "_allassays.png"),
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
				       ncol = 1)
			)
		dev.off()

		print("Test PCs:")
		print(testpcs)
	
		if(testpcs) {

			png(paste0(outputfolder, "/", mca_name,
				   "_Seurat3_nhvg_", nhvg, "_npcs_", npcs, "_allassays_dimensions.png"),
			    width = 2500, height = 1500)
	
			plot(
				plot_grid(# rna
				       plot_grid(rna_title,
			      		plot_grid(rna_mca[["testdim"]][["elbow"]], rna_mca[["testdim"]][['js']], ncol = 2),
			      		ncol = 1, rel_heights = c(0.1, 1)),
			       # sct
			       plot_grid(sct_title,
			      		plot_grid(sct_mca[["testdim"]][["elbow"]], sct_mca[["testdim"]][['js']], ncol = 2),
			      		ncol = 1, rel_heights = c(0.1, 1)),
			       ncol = 1)
			)
			dev.off()
			
		}

	}

	if(w) {
		# --- Write cell_data_set object with merged samples
		mca_file <- paste0(outputfolder, '/', mca_name, '_', npcs, 'PCs_', nhvg, 'g_seurat3.rds')
		out_mca <- sct_mca[["mca"]]
		DefaultAssay(out_mca) <- "SCT"
		
		out_mca <- RunPCA(out_mca,
			      npcs = npcs, 
			      verbose = TRUE,
			      seed.use = 123)

		out_mca <- RunUMAP(out_mca,
			       dims = seq(npcs),
			       reduction = "pca",
			       n.neighbors = nneigh,
			       metric = "cosine",
			       min.dist = 0.3)
		       
		out_mca <- FindNeighbors(out_mca,
					dims = seq(npcs), 
					reduction = "pca",
					verbose = TRUE)

		out_mca <- FindClusters(out_mca,
					verbose = TRUE)
		       
		print(out_mca)
	
		saveRDS(out_mca, file = mca_file)
	
		message(paste("The cell_data_set object (Seurat3)", 
			      mca_file, "has been created."))
	}
		      
}

# -----------------

if(is.null(arguments$nhvg)) {

	nhvg <- 500

} else {

	nhvg <- as.numeric(arguments$nhvg)

}

if(!is.null(arguments$nneigh)) {

	nneigh <- arguments$nneigh

} else {

	nneigh <- 30

}

if(is.null(arguments$method)) {

	method <- "SCT"

} else {

	method <- arguments$method

}

testpcs <- as.logical(arguments$testpcs)
if(is.null(arguments$specie)) {

	specie <- "mm"

} else {

	specie <- arguments$specie

}

if(is.null(arguments$write)) {

	w <- FALSE

} else {

	w <- as.logical(arguments$write)

}

if(is.null(arguments$explore)) {

	explore <- TRUE

} else {

	explore <- as.logical(arguments$explore)

}
if(is.null(arguments$covars)) {

	covars <- NULL

} else {

	covars <- unlist(strsplit(arguments$covars, "///"))

}

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
cdsfiles <- arguments$file
npcs <- as.numeric(arguments$npcs)

# Maximum vector size whithin the session env
# Mb*1024^2
# 100000 * 1024^2 = 104857600000
options(future.globals.maxSize=91943040000000)

# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

# ------------
# Reading merge mca file  
# ------------
cdss <- lapply(cdsfiles, function(...) readRDS(paste0(inputfolder, '/', ...)))

cdsnames <- gsub("\\.rds", "", cdsfiles)
names(cdss) <- cdsnames
print("This files will be SCT and Factor size normalized (fill up RNA and SCT assays):")
print(names(cdss))

if(any(!sapply(cdss, class)=="cell_data_set")) {
	cl <- sapply(cdss, class)
	print(cl)
	#stop("Please make sure that all .rds files are cell_data_set objects.")
	mcas <- cdss
} else {
	# From monocle3 to Seurat3 object
	mcas <- lapply(cdss, cds2mca)
}

# refbiomart
# ------------------------------------------------------------------------------
# Marker expression visualization (umap)
# ------------------------------------------------------------------------------

if(specie == "hg") {

	refbiomart <- data.table::fread("./data/ensg_biomart_20191206.tsv")

	sc.markers <- merge(data.frame("gene" = cc.genes.updated.2019$s.genes,
				       "phase" = "s"),
			    refbiomart, by.x = "gene", 
			    by.y = "Gene name", all.x = TRUE) 

	g2mc.markers <- merge(data.frame("gene" = cc.genes.updated.2019$g2m.genes,
				       "phase" = "g2m"),
			    refbiomart, by.x = "gene", 
			    by.y = "Gene name", all.x = TRUE) 

	cc.markers <- rbind(sc.markers, g2mc.markers)

} else if(specie == "mm"){

	refbiomart <- data.table::fread("./data/human_to_mouse_biomart_20200519.tsv")

	sc.markers <- merge(data.frame("gene" = cc.genes.updated.2019$s.genes,
				       "phase" = "s"),
			    refbiomart, by.x = "gene", 
			    by.y = "Gene name", all.x = TRUE) 

	g2mc.markers <- merge(data.frame("gene" = cc.genes.updated.2019$g2m.genes,
				       "phase" = "g2m"),
			    refbiomart, by.x = "gene", 
			    by.y = "Gene name", all.x = TRUE) 

	cc.markers <- rbind(sc.markers, g2mc.markers)
	cc.markers <- cc.markers %>%
		filter(`Mouse gene name` != '') %>%
		select(-c(`Gene stable ID`, `gene`))
	colnames(cc.markers)[c(2:3)] <- c('Gene stable ID', 'gene')
	print(head(cc.markers,2))

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}

# Visualizing blood
if(!is.null(arguments$blood)) {

	lapply(names(mcas), function(mca_name) {

		mca <- mcas[[mca_name]]

		mca <- NormalizeData(mca,
			     assay = "RNA",
			     normalization.method = "LogNormalize", 
			     scale.factor = 10000,
			     margin = 1, 
			     verbose = TRUE)

		meta <- mca@meta.data
		meta[['barcode']] <- rownames(meta)

		m <- data.matrix(mca@assays$RNA@data)
	
		g <- c("ENSMUSG00000052187", "ENSMUSG00000052217", "ENSMUSG00000029580")
		g <- g[which(g %in% rownames(m))]
		print(g)
		sub_mat <- data.frame(t(m[unique(g), ]))
		sub_mat[['barcode']] <- rownames(sub_mat)
		sub_eset <- merge(sub_mat, meta, by = "barcode")

		if(any(!grepl("_RAW", colnames(sub_eset)))) {
			   colnames(sub_eset)[grep("detected", colnames(sub_eset))] <- paste(colnames(sub_eset)[grep("detected", colnames(sub_eset))], "RAW", sep = "_") 
		}
		print(head(sub_eset))

		gg <- ggplot(sub_eset, aes(x = cell_depth_counts, y = cell_depth_genes_detected_RAW)) +
			geom_point(aes(colour = ENSMUSG00000052187), alpha = .6, size = 0.6) +
			theme_classic() +
			scale_colour_gradient2(low = "grey40", mid = "white", high = "red",
					       midpoint = quantile(sub_eset$ENSMUSG00000052187, 0.5)) +
			ggtitle("Hbb-y blood")

		gg_actb <- ggplot(sub_eset, aes(x = cell_depth_counts, y = cell_depth_genes_detected_RAW)) +
			geom_point(aes(colour = ENSMUSG00000029580), alpha = .6, size = 0.6) +
			scale_colour_gradient2(low = "grey40", mid = "white", high = "red",
					       midpoint = quantile(sub_eset$ENSMUSG00000029580, 0.5)) +
			theme_classic() +
			ggtitle("Actb housekeeping")

		gg_log <- ggplot(sub_eset, aes(x = log10(cell_depth_counts), y = log10(cell_depth_genes_detected_RAW))) +
			geom_point(aes(colour = ENSMUSG00000052187), alpha = .6, size = 0.6) +
			theme_classic() +
			scale_colour_gradient2(low = "grey40", mid = "white", high = "red",
					       midpoint = quantile(sub_eset$ENSMUSG00000052187, 0.5)) +
			ggtitle("Hbb-y blood")

		gg_actb_log <- ggplot(sub_eset, aes(x = log10(cell_depth_counts), y = log10(cell_depth_genes_detected_RAW))) +
			geom_point(aes(colour = ENSMUSG00000029580), alpha = .6, size = 0.6) +
			scale_colour_gradient2(low = "grey40", mid = "white", high = "red",
					       midpoint = quantile(sub_eset$ENSMUSG00000029580, 0.5)) +
			theme_classic() +
			ggtitle("Actb housekeeping")

		png(paste0(outputfolder, "/", mca_name, ".png"),
		    width = 900, height = 400)

		plot(
	
		       plot_grid(gg, gg_actb, gg_log, gg_actb_log, ncol = 2)

		       )

		dev.off()

		pdf(paste0(outputfolder, "/", mca_name, ".pdf"),
		    width = 15, height = 8)

		plot(
	
		       plot_grid(gg, gg_actb, gg_log, gg_actb_log, ncol = 2)

		       )

		dev.off()


	})
}


# Normalizing

for(sample_name in names(mcas)) {

	print(sample_name)
	get_norm_single_mca(mcas[[sample_name]], sample_name, nhvg, npcs, nneigh, testpcs, colgroup, cc.markers, ncores, w, covars, explore)
	
}
