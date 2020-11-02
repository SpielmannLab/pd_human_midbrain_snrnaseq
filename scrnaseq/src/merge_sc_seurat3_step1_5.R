#
# Integration protocols for several already pre-processed (Seurat3) objects 
#

" Integration protocols for several already pre-processed (Seurat3) objects

Usage: merge_sc_seurat3_step1_5.R --jobname=<value> --npcs=<value> [--nhvg=<value>] --metafile=<file> --samplevar=<value> --method=<value> --normethod=<value> [--nneigh=<value>] [--refsamples=<value>] --ncores=<value> --testpcs=<logical> --infolder=<folder> --outfolder=<folder> [<file> <file>...]

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --npcs=<value>       Number of principal components to use.
  --nhvg=<value> Optional. Number of genes to use during integration.
  --metafile=<file>    Metadata samples x features.
  --samplevar=<value>  Var. name with sample ids.
  --method=<value>     Integration strategy: full_cca, ref_cca, ref_rpc
  --normethod=<value>  Norm. method. either fs (Factor.size) or sct (SCT)
  --nneigh=<value> Optional. Number of neirest neighbors fro umap embeding.
  --refsamples=<value> Reference samples for ref_cca and ref_rpc methods. Save comp. power.
  --ncores=<value>     Number of threads to use.
  --testpcs=<logical>  Wheather to test pc explanation.
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

# -------------------------------
pymagic <- reticulate::import("magic")
# ------------------------------

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


integrate_samples <- function(mcas, method, nm, nhvg, npcs, nneigh, sample_meta, refsample) {

	if(nm == "sct") {

		mcas <- lapply(mcas, function(x) {
			       DefaultAssay(x) <- "SCT"
			       x
		})


		if(method == "full_cca") {
			

			mca.features <- SelectIntegrationFeatures(object.list = mcas,
								  nfeatures = nhvg)

			print(head(mca.features))

			mcas <- PrepSCTIntegration(object.list = mcas, 
				   anchor.features = mca.features,
				   verbose = TRUE)

			print(names(mcas))

			mca.anchors <- FindIntegrationAnchors(object.list = mcas, 
				      normalization.method = "SCT", 
				      anchor.features = mca.features)

			mca <- IntegrateData(anchorset = mca.anchors,
					     normalization.method = "SCT",
					     verbose = TRUE)

			mca <- RunPCA(mca,
			      features = mca@assays$integrated@var.features,
			      assay = "integrated",
			      npcs = npcs, 
			      verbose = TRUE,
			      seed.use = 444)

			mca <- RunUMAP(mca,
			       assay = "integrated",
			       dims = seq(npcs),
			       reduction = "pca",
			       n.neighbors = nneigh,
			       metric = "cosine",
			       min.dist = 0.3, 
			       spread = 1, 
			       set.op.mix.ratio = 1,
			       local.connectivity = 1L, 
			       repulsion.strength = 1)
	
			# Sample metadata
			meta <- mca@meta.data
			meta[["barcode_idx"]] <- rownames(meta)
		        
			meta %>%
			        mutate('sample' = gsub(paste0("_filtered_FILTERED_qc_", 
							      npcs, "PCs_", nhvg, 
							      "g_seurat3|", 
							      nhvg, "g_seurat3.rds"), "", `sample`)) -> meta
                                                                                                                                                                                                                                                            
			print(head(meta,2))
			print(table(meta$sample))

			mt <- vis_meta(mca)
			nc <- ceiling(sqrt(length(mt)))

			png(paste0(outputfolder, "/", jobname, "_npcs_", npcs, "_nvgs_", nhvg,
				   "_Seurat3_SCT_integrated_", method, '_', nm, ".png"),
			    width = 600*nc, height = 600*nc)
	
				plot(
				     plot_grid(plotlist = mt, ncol = nc)
				)
		
			dev.off()

			return(mca)

		} else if(method == "ref_cca") {

			mca.features <- SelectIntegrationFeatures(object.list = mcas, 
						  nfeatures = nhvg)

			mcas <- PrepSCTIntegration(object.list = mcas, 
					   anchor.features = mca.features,
					   verbose = TRUE)


			refs <- grep(refsamples, names(mcas))
			print(names(mcas)[refs])

			mca.anchors <- FindIntegrationAnchors(object.list = mcas, 
					      normalization.method = "SCT", 
					      anchor.features = mca.features, 
					      reference = refs)

			mca <- IntegrateData(anchorset = mca.anchors, 
			     normalization.method = "SCT",
			     verbose = TRUE)
	
			mca <- RunPCA(mca,
			      features = mca@assays$integrated@var.features,
			      assay = "integrated",
			      npcs = npcs, 
			      verbose = TRUE,
			      seed.use = 444)

			mca <- RunUMAP(mca,
			       assay = "integrated",
			       dims = seq(npcs),
			       reduction = "pca",
			       n.neighbors = nneigh,
			       metric = "cosine",
			       min.dist = 0.3, 
			       spread = 1, 
			       set.op.mix.ratio = 1,
			       local.connectivity = 1L, 
			       repulsion.strength = 1)
	
			# Sample metadata
			meta <- mca@meta.data
			meta[["barcode_idx"]] <- rownames(meta)                                    
	
			meta %>%
			        mutate('sample' = gsub(paste0("_filtered_FILTERED_qc_", npcs, "PCs_", nhvg, "g_seurat3"), "", `sample`)) -> meta
		                                                                                                                                                                                                                                                                    
		
			print(head(meta,2))
			print(table(meta$sample))
                                                                                                                                                                                                                                    
			print(head(sample_meta,2))
			print(table(sample_meta$sample))
		
			new_meta <- merge(meta, sample_meta, by = "sample")
			rownames(new_meta) <-   new_meta[["barcode_idx"]]
			table(new_meta$sample) 
			print(head(new_meta,2))
	                                                                                                                                                                                                                                                                                          
			mca@meta.data <- new_meta[colnames(mca), ]
	
			mt <- vis_meta(mca)
			nc <- ceiling(sqrt(length(mt)))
		
			png(paste0(outputfolder, "/", jobname, "_npcs_", npcs, "_nvgs_", nhvg,
				   "_Seurat3_SCT_integrated_", method, '_', nm, ".png"),
			    width = 600*nc, height = 600*nc)
																					
				plot(
				     plot_grid(plotlist = mt, ncol = nc)
				)
		
			dev.off()
	
			return(mca)

		}
	
	} else if(nm == "fs") {

		mcass <- lapply(seq(mcas), function(mca_n) {
			  meta <- mcas[[mca_n]]@meta.data
			  rownames(meta) <- paste0(rownames(meta), '_', mca_n)
			  mca <- mcas[[mca_n]][["RNA"]]
			  print(mca)
			  mca <- RenameCells(object = mca,
					     new.names = paste0(colnames(mca), '_', mca_n))
			  CreateSeuratObject(counts = mca@counts,
					     assay = "RNA",
					     meta.data = meta[colnames(mca), ])
		})

		names(mcass) <- names(mcas)
		
		mcas <- lapply(mcass, function(x) {

				       DefaultAssay(x) <- "RNA"

				       x <- NormalizeData(x,
							  assay = "RNA",
							  normalization.method = "LogNormalize",
							  scale.factor = 10000,
							  margin = 1,
							  verbose = TRUE)

				       x <- FindVariableFeatures(x,
								 assay = "RNA",
								 selection.method = "vst",
								 nfeatures = nhvg,
								 verbose = FALSE)

				       x <- ScaleData(x,
						      features = x@assays$RNA@var.features,
						      assay = "RNA",
						      verbose = TRUE)

				       x
			})


		if(method == "full_cca") {

			mca.features <- SelectIntegrationFeatures(object.list = mcas, 
								  nfeatures = nhvg)


			mca.anchors <- FindIntegrationAnchors(object.list = mcas,
							      anchor.features = mca.features,							      
							      reduction = "cca",
							      dims = seq(npcs))

			mca <- IntegrateData(anchorset = mca.anchors,
					     normalization.method = "LogNormalize",
					     features = mca.features,
					     verbose = TRUE)

			mca <- ScaleData(mca,
					 features = mca@assays$integrated@var.features, 
					 assay = "integrated",
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

			mca <- RunPCA(mca,
				      features = mca@assays$integrated@var.features,
				      assay = "integrated",
				      npcs = npcs, 
				      verbose = TRUE,
				      seed.use = 444)

			mca <- RunUMAP(mca,
			       assay = "integrated",
			       dims = seq(npcs),
			       reduction = "pca",
			       n.neighbors = nneigh,
			       metric = "cosine",
			       min.dist = 0.3, spread = 1, 
			       set.op.mix.ratio = 1,
			       local.connectivity = 1L, 
			       repulsion.strength = 1)
	
			# Sample metadata
			meta <- mca@meta.data
			meta[["barcode_idx"]] <- rownames(meta)
		                                                                                                                                                                                                                                                                    
			meta %>%
			        mutate('sample' = gsub(paste0("_filtered_FILTERED_qc_", npcs, "PCs_", nhvg, "g_seurat3"), "", `sample`)) -> meta
	
			print(head(meta,2))
			print(table(meta$sample))

			mt <- vis_meta(mca)
			nc <- ceiling(sqrt(length(mt)))

			png(paste0(outputfolder, "/", jobname, "_npcs_", npcs, "_nvgs_", nhvg,
				   "_Seurat3_factor_size_integrated_", method, '_', nm, ".png"),
			    width = 600*nc, height = 600*nc)
	
				plot(
				     plot_grid(plotlist = mt, ncol = nc)
				)
		
			dev.off()

			return(mca)
	
		} else if(method == "full_rpca") {

			mca.features <- SelectIntegrationFeatures(object.list = mcas,
								  nfeatures = nhvg,
								  assay = rep("RNA", length(mcas)))

			scale.mcas <- lapply(mcas, function(x) {

			       x <- ScaleData(x,
					      assay = "RNA",
					      features = mca.features, 
					      verbose = TRUE)
			       
			       x <- RunPCA(x,
					   assay = "RNA",
					   features = mca.features, 
					   npcs = npcs,
					   verbose = FALSE)

			       return(x)
			})

			print("These samples will be full rpca integrated:")
			print(names(scale.mcas))
	
			mca.anchors <- FindIntegrationAnchors(object.list = scale.mcas,
							      assay = rep("RNA", length(scale.mcas)),
							      anchor.features = mca.features,
							      reduction = "rpca",
							      dims = seq(npcs))

			mca <- IntegrateData(anchorset = mca.anchors,
					     normalization.method = "LogNormalize",
					     dims = seq(npcs), 
					     verbose = TRUE)

			mca <- ScaleData(mca,
					 assay = "integrated",
					 verbose = TRUE)

			mca <- RunPCA(mca,
			      features = mca@assays$integrated@var.features,
			      assay = "integrated",
			      npcs = npcs, 
			      verbose = TRUE,
			      seed.use = 444)

			mca <- RunUMAP(mca,
			       assay = "integrated",
			       dims = seq(npcs),
			       reduction = "pca",
			       n.neighbors = nneigh,
			       metric = "cosine",
			       min.dist = 0.3, 
			       spread = 1, 
			       set.op.mix.ratio = 1,
			       local.connectivity = 1L, 
			       repulsion.strength = 1)

			# Sample metadata >> cell metadata
			meta <- mca@meta.data 
			meta[["barcode_idx"]] <- rownames(meta)
	
			meta %>%
			        mutate('sample' = gsub(paste0("_filtered_FILTERED_qc_", npcs, "PCs_", nhvg, "g_seurat3"), "", `sample`)) -> meta

			mt <- vis_meta(mca)
			nc <- ceiling(sqrt(length(mt)))
		
			png(paste0(outputfolder, "/", jobname, "_npcs_", npcs, "_nvgs_", nhvg,
				   "_Seurat3_factor_size_integrated_", method, '_', nm, ".png"),
			    width = 600*nc, height = 600*nc)
																					
				plot(
				     plot_grid(plotlist = mt, ncol = nc)
				)
		
			dev.off()
		
			return(mca)

		} else if(method == "ref_rpca") {

		
			mca.features <- SelectIntegrationFeatures(object.list = mcas,
								  nfeatures = nhvg,
								  assay = rep("RNA", length(mcas)))

			scale.mcas <- lapply(mcas, function(x) {

			       x <- ScaleData(x,
					      assay = "RNA",
					      features = mca.features, 
					      verbose = TRUE)
			       
			       x <- RunPCA(x,
					   assay = "RNA",
					   features = mca.features, 
					   npcs = npcs,
					   verbose = FALSE)

			       return(x)
			})

			refs <- grep(refsamples, names(scale.mcas))
			print(names(scale.mcas)[refs])
	
			mca.anchors <- FindIntegrationAnchors(object.list = scale.mcas,
							      assay = rep("RNA", length(scale.mcas)),
							      anchor.features = mca.features,
							      reduction = "rpca",
							      reference = refs,
							      dims = seq(npcs))

			mca <- IntegrateData(anchorset = mca.anchors,
					     normalization.method = "LogNormalize",
					     dims = seq(npcs), 
					     verbose = TRUE)

			mca <- ScaleData(mca,
					 assay = "integrated",
					 verbose = TRUE)

			mca <- RunPCA(mca,
			      features = mca@assays$integrated@var.features,
			      assay = "integrated",
			      npcs = npcs, 
			      verbose = TRUE,
			      seed.use = 444)

			mca <- RunUMAP(mca,
			       assay = "integrated",
			       dims = seq(npcs),
			       reduction = "pca",
			       n.neighbors = nneigh,
			       metric = "cosine",
			       min.dist = 0.3, 
			       spread = 1, 
			       set.op.mix.ratio = 1,
			       local.connectivity = 1L, 
			       repulsion.strength = 1)

			# Sample metadata >> cell metadata
			meta <- mca@meta.data
			meta[["barcode_idx"]] <- rownames(meta)
		                                                                                                                                                                                                                                                                    
			meta %>% 
			        mutate('sample' = gsub(paste0("_filtered_FILTERED_qc_", npcs, "PCs_", nhvg, "g_seurat3"), "", `sample`)) -> meta
	
			print(head(meta, 2))
			print(head(sample_meta, 2))

			new_meta <- merge(meta, sample_meta, by = "sample")
			rownames(new_meta) <- new_meta[["barcode_idx"]]
			print(head(new_meta, 2))
		                                                                                                                                                                                                                                                                                          
			mca@meta.data <- new_meta[colnames(mca), ]     

			mt <- vis_meta(mca)
			nc <- ceiling(sqrt(length(mt)))
		
			png(paste0(outputfolder, "/", jobname, "_npcs_", npcs, "_nvgs_", nhvg,
				   "_Seurat3_factor_size_integrated_", method, '_', nm, ".png"),
			    width = 600*nc, height = 600*nc)
																					
				plot(
				     plot_grid(plotlist = mt, ncol = nc)
				)
		
			dev.off()
		
			return(mca)
		}
	}
}
# ----

method <- arguments$method
nm <- arguments$normethod
nneigh <- arguments$nneigh
refsamples <- arguments$refsamples

if(!is.null(arguments$nhvg)) {

	nhvg <- as.numeric(arguments$nhvg)

} else {

	nhvg <- 500
}

if(!is.null(arguments$nneigh)) {

	nneigh <- as.numeric(arguments$nneigh)

} else {

	nneigh <- 30
}

if(!is.null(arguments$refsamples)) {

	refsamples <- arguments$refsamples

} else {

	refsamples <- NA
}


nhvg <- as.numeric(arguments$nhvg)
metafile <- arguments$metafile
samplevar <- arguments$samplevar
cdsfiles <- arguments$file
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
npcs <- as.numeric(arguments$npcs)

testpcs <- as.logical(arguments$testpcs)

# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

options(future.globals.maxSize=61943040000000)

# ------------
# Reading cds 
# ------------

mcas <- lapply(cdsfiles, function(x) readRDS(paste0(inputfolder, '/', x)))
names(mcas) <- gsub("\\.rds", "", cdsfiles)

message(paste("The cell_data_set objects (Seurat3)", 
	      paste(cdsfiles, collapse = " "), "has been read and will be integrated."))

print(names(mcas))

# Adding sample id column ~ rds file name. 
mcas <- lapply(names(mcas), function(s) {
		       mca <- mcas[[s]]
		       mca@meta.data[["sample"]] <- s
		       mca
	})

names(mcas) <- gsub("\\.rds", "", cdsfiles)
print(names(mcas))
print(lapply(mcas, dim))
print(lapply(mcas, function(x) head(rownames(x))))
print(lapply(mcas, function(x) head(which(is.na(rownames(x))))))
print(lapply(mcas, function(x) head(x@meta.data,2)))


# -----

# Read sample metadata
print(metafile)
sample_meta <- data.table::fread(metafile)

print("Samples to be integrated:")
print(mcas)
print("Using the strategy:")
print(method)
print("And the normalization:")
print(nm)
print("Using this # of variable genes:")
print(nhvg)
print("Considering these # of PCs")
print(npcs)
print("# neighbors for UMAP embedings:")
print(nneigh)
print("Using these samples as reference:")
print(refsamples)
print("Sample metadata:")
print(head(sample_meta,2))

i.mca <- integrate_samples(mcas, method, nm, nhvg, npcs, nneigh, sample_meta, refsamples)

# integrate data and keep full geneset
mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_integrated_', "npcs_", npcs, "_nvgs_", nhvg, '_', method, '_', nm, '.rds')
print(i.mca)
print(mca_file)
saveRDS(i.mca, file = mca_file)
message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been created."))
