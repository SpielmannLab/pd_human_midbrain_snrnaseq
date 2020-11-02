#
# Identify marker genes in a seurat3 object for given cell-phenodata variable
#

" Identify marker genes in a seurat3 object for given cell-phenodata variable

Usage: get_mkr_genes_seurat3.R --jobname=<value> [--mcafile=<file>] [--assay=<value>] [--mincells=<value>] [--min_cell_prop=<value>] [--min_logfc=<value>] --resolution=<value> --specie=<value> --infolder=<folder> --outfolder=<folder> --ncores=<value> 

Options:
  -h --help               Show this screen.
  --version               00.99.01
  --jobname=<file>        Descriptive name for your experiment.
  --mcafile=<file>        Seurat3 object after clustering.
  --assay=<value>         Assay whose data slot will be used.
  --mincells=<value>      Min. number of cell in clusters.
  --min_cell_prop=<value> Min. proprtion of cells of a given cluster expressing a given gene to be considered.
  --min_logfc=<value>     Min. logFC among cell clusters of a given gene to be considered.
  --resolution=<value>    Clustering resolution. 	
  --specie=<value>        Either hg or mm.
  --infolder=<file>       Path to the single_cell_data .rds files.  
  --outfolder=<file>      Path to results folder.
  --ncores=<value>        Number of processors to use

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "Seurat", "future", "sctransform",
	  'Rtsne', 'data.table', "monocle3")

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

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
res <- as.numeric(arguments$resolution)
specie <- arguments$specie

if(!is.null(arguments$mincells)) {
	mincells <- as.numeric(arguments$mincells)
} else {
	mincells <- 50
}

if(!is.null(arguments$assay)) {
	`assay` <- arguments$assay
} else {
	`assay` <- "RNA"
}

if(!is.null(arguments$min_cell_prop)) {
	min_cell_prop <- as.numeric(arguments$min_cell_prop)
} else {
	min_cell_prop <- 0.25
}

if(!is.null(arguments$min_logfc)) {
	min_logfc <- as.numeric(arguments$min_logfc)
} else {
	min_logfc <- 0.25
}

# ---------------------------------------------------
# Setting multicore env for some Seurat applications
# ---------------------------------------------------
ncores <- as.numeric(arguments$ncores)
plan("multiprocess", workers = ncores)
plan()

# ----------------------------------------------------
#  Read cell_data_set object with cluster metadata
# ----------------------------------------------------

# --- Read cell_data_set object with merged samples
if(!is.null(arguments$mcafile)) {

	mca_file <- paste0(inputfolder, '/', arguments$mcafile)
	jobname <- gsub("\\.rds", "", arguments$mcafile)
	
} else {

	mca_file <- paste0(inputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')
	
}


print(mca_file)

mca <- readRDS(mca_file)
print(mca)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# ----- set ident -------
print(colnames(mca@meta.data)[grep("_snn_res.", colnames(mca@meta.data))])
#
ident_col <- colnames(mca@meta.data)[grep(paste0("_snn_res.", res), colnames(mca@meta.data))]
ident_col <- ident_col[1]
Idents(object = mca) <- ident_col 
print(table(Idents(mca)))

# Filter out small clusters

ident_to_keep <- names(table(Idents(mca)))[which(table(Idents(mca)) > mincells)]
print(ident_to_keep)
cells_to_keep <- rownames(mca@meta.data)[which(mca@meta.data[[ident_col]] %in% ident_to_keep)]
print(paste("Cells to keep:", length(cells_to_keep)))
mca <- mca[, cells_to_keep]
print(mca)

# ----------------------------------------------------
#  Visualizing the selected cluster granularity
# ----------------------------------------------------

g_cluster <- Seurat::DimPlot(mca, 
		     reduction = "umap", 
		     pt.size = 0.01) +
			     ggplot2::theme(legend.position = "bottom") + 
			     ggtitle(paste0(jobname, "_", res))

# --- title
title <- ggdraw() + 
  draw_label(
    paste0(jobname, ' res: ', res),
    fontface = 'bold',
    hjust = 0.5
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "_", res), "_umap.jpeg"), 
     	    width = 980, height = 900, pointsize = 74, quality = 100)
    plot_grid(title, plot_grid(g_cluster, ncol = 1), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "_", res), "_umap.pdf"))
    plot_grid(title, plot_grid(g_cluster, ncol = 1), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

# ----------------------------------------------------
# Identifying the marker genes
# ----------------------------------------------------
# Mb*1024^2
l <-3000*1024^2
options(future.globals.maxSize=l)

# --------------------
# ROC method
# --------------------

mca.markers.roc <- FindAllMarkers(mca,
				  assay = `assay`,
				  only.pos = TRUE,
				  min.pct = min_cell_prop,
				  test.use = "roc",
				  logfc.threshold = min_logfc)

print(head(mca.markers.roc, 2))

write.table(mca.markers.roc,
	    paste0(outputfolder, '/', jobname, "_res_", 
		   gsub("\\.", "_", res), "_marker_genes_roc.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

cns <- colnames(mca.markers.roc)
cns[grep("pct|cluster", cns)] <- paste0(cns[grep("pct|cluster", cns)], "_auc")
colnames(mca.markers.roc) <- cns

# --------------------
# Wilcox method
# --------------------

mca.markers.wilcox <- FindAllMarkers(mca,
				     only.pos = TRUE, 
				     logfc.threshold = min_logfc,
				     assay = `assay`,
				     min.pct = min_cell_prop,
				     test.use = "wilcox")

print(head(mca.markers.wilcox, 2))

write.table(mca.markers.wilcox,
	    paste0(outputfolder, '/', jobname, "_res_", 
		   gsub("\\.", "_", res), "_marker_genes_wilcox.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

cns <- colnames(mca.markers.wilcox)
cns[grep("pct|cluster", cns)] <- paste0(cns[grep("pct|cluster", cns)], "_wilcox")
colnames(mca.markers.wilcox) <- cns

# ---------------------------
# Negative binomial method 
# --------------------------

mca.markers.bn <- FindAllMarkers(mca,
				 only.pos = TRUE, 
				 logfc.threshold = min_logfc,
				 assay = `assay`,
				 min.pct = min_cell_prop,
				 test.use = "negbinom")

print(head(mca.markers.bn, 2))

write.table(mca.markers.bn,
	    paste0(outputfolder, '/', jobname, "_res_", 
		   gsub("\\.", "_", res), "_marker_genes_negbin.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

cns <- colnames(mca.markers.bn)
cns[grep("pct|cluster|p_val|avg_logFC", cns)] <- paste0(cns[grep("pct|cluster|p_val|avg_logFC", cns)], "_negbin")
colnames(mca.markers.bn) <- cns

mca.markers <- merge(mca.markers.roc, mca.markers.wilcox, by = "gene", all = TRUE)
mca.markers <- merge(mca.markers, mca.markers.bn, by = "gene", all = TRUE)

# Check if stable gene ID was provided with or without version
g2test <- mca.markers[['gene']][1]

if(grepl("\\.\\d+$", g2test)) {
	mca.markers[['gene_version']] <- mca.markers[['gene']]
	mca.markers[['gene']] <- gsub("\\.\\d+$", "", mca.markers[['gene_version']])
} else {
	mca.markers[['gene_version']] <- mca.markers[['gene']]
}

print(head(mca.markers, 2))

# ------------------------------------------------------------------------------
# Marker expression visualization (umap)
# ------------------------------------------------------------------------------

if(specie == "hg") {

	refbiomart <- data.table::fread("./data/ensg_biomart_20191206.tsv")
	mca.markers <- merge(mca.markers, refbiomart, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE) 

} else if(specie == "mm"){

	refbiomart <- data.table::fread("./data/biomart_mm_20191206.tsv")
	mca.markers <- merge(mca.markers, refbiomart, by.x = "gene", by.y = "Gene stable ID", all.x = TRUE) 

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}

print(head(mca.markers))

write.table(mca.markers, paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_marker_genes.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

auc <- dplyr::filter(mca.markers, !is.na(cluster_auc))
auc <- split(auc, auc[['cluster_auc']]) 
auc <- lapply(auc, function(w) head(dplyr::arrange(w, -myAUC), 200))
topmarkers <- Reduce(rbind, lapply(auc, function(top) {
					   top <- dplyr::mutate(top, d = pct.1_auc - pct.2_auc)
					   head(dplyr::arrange(top, -d), 100)
	})
)

print(table(topmarkers$cluster_auc))
print(dim(topmarkers))

write.table(topmarkers, paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_top100_marker_genes.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

mkrs <- split(topmarkers, topmarkers$cluster_auc) 
mkrs <- Filter(Negate(is.null), mkrs)

if(any(grepl("gene_version", colnames(mkrs[[1]])))) {
	   mkrs <- lapply(mkrs, function(m) head(m$gene_version, 50))
} else {
	   mkrs <- lapply(mkrs, function(m) head(m$gene, 50))
	
}

parallel::mclapply(names(mkrs), function(cls) {

	g <- mkrs[[cls]]
	gnames <- filter(topmarkers, gene%in%g)[['Gene name']]

	ggs <- FeaturePlot(mca, 
			   features = unique(g), 
			   combine=FALSE)

	ggs <- lapply(seq(ggs), function(i) {
			      igene <- names(ggs[[i]]$data)[4]
			      ggs[[i]] + 
				      ggtitle(paste(unique(dplyr::filter(topmarkers,
							    gene_version==igene)[['Gene name']]),
				      collapse = '///')) 
			   })

	title <- ggdraw() + 
		draw_label(
			paste0(jobname, " C", cls),
			fontface = 'bold',
			#x = 0,
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
  	)

	jpeg(paste0(outputfolder, '/', jobname, "_C", cls, "_res_", gsub("\\.", "", res), "_auc_marker_genes.jpeg"), 
	    width = 1980, height = 1980, pointsize = 74, quality = 100)

	plot(
	    plot_grid(title,
		      plot_grid(plotlist = ggs, ncol = 5),
		      ncol = 1, rel_heights = c(0.01, 1))
	)

	dev.off()
}, mc.cores = ncores)

# ------------------------------------------------------------------------------
# Marker expression visualization (heatmap) [AUC]
# ------------------------------------------------------------------------------

if(length(colnames(mca)) > 10000) {
	ncell <- 10000
} else {
	ncell <- length(colnames(mca))
}

cells <- sample(colnames(mca), ncell, replace=FALSE)

print(head(cells))

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_auc_marker_genes_heatmap.pdf"), 
     width = 18, height = 15)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_auc_marker_genes_heatmap.jpeg"), 
 	    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()


#
# Wilcoxon
#

wilcox <- dplyr::filter(mca.markers, !is.na(cluster_wilcox))
wilcox <- split(wilcox, wilcox[['cluster_wilcox']]) 
wilcox <- lapply(wilcox, function(w) head(dplyr::arrange(w, p_val_adj), 200))
topmarkers <- Reduce(rbind, lapply(wilcox, function(top) {
					   top <- dplyr::mutate(top, d = pct.1_wilcox - pct.2_wilcox)
					   head(dplyr::arrange(top, -d), 100)
	})
)

print(table(topmarkers$cluster_wilcox))

write.table(topmarkers, paste0(outputfolder, '/', jobname, "_res_", 
gsub("\\.", "", res),"_wilcox_top100_marker_genes.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

mkrs <- split(topmarkers, topmarkers$cluster_auc) 
mkrs <- Filter(Negate(is.null), mkrs)

if(any(grepl("gene_version", colnames(mkrs[[1]])))) {

	   mkrs <- lapply(mkrs, function(m) head(m$gene_version, 50))

} else {

	   mkrs <- lapply(mkrs, function(m) head(m$gene, 50))
	
}


parallel::mclapply(names(mkrs), function(cls) {

	g <- mkrs[[cls]]
	gnames <- filter(topmarkers, gene%in%g)[['Gene name']]

	ggs <- FeaturePlot(mca, 
			   features = unique(g), 
			   combine=FALSE)

	ggs <- lapply(seq(ggs), function(i) {
			      igene <- names(ggs[[i]]$data)[4]
			      ggs[[i]] + 
				      ggtitle(paste(unique(dplyr::filter(topmarkers,
							    gene_version==igene)[['Gene name']]),
				      collapse = '///')) 
			   })

	title <- ggdraw() + 
		draw_label(
			paste0(jobname, " C", cls),
			fontface = 'bold',
			#x = 0,
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
  	)

	jpeg(paste0(outputfolder, '/', jobname, "_C", cls, "_res_", gsub("\\.", "", res), "_wilcox_marker_genes.jpeg"), 
	    width = 1980, height = 1980, pointsize = 74, quality = 100)

	plot(
	    plot_grid(title,
		      plot_grid(plotlist = ggs, ncol = 5),
		      ncol = 1, rel_heights = c(0.01, 1))
	)

	dev.off()
}, mc.cores = ncores)

# ------------------------------------------------------------------------------
# Marker expression visualization (heatmap)
# ------------------------------------------------------------------------------

if(length(colnames(mca)) > 10000) {
	ncell <- 10000
} else {
	ncell <- length(colnames(mca))
}

cells <- sample(colnames(mca), ncell, replace=FALSE)

print(head(cells))

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_wilcox_marker_genes_heatmap.pdf"), 
     width = 18, height = 15)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_wilcox_marker_genes_heatmap.jpeg"), 
 	    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

#
# Neg binomial
#

wilcox <- dplyr::filter(mca.markers, !is.na(cluster_negbin))
wilcox <- split(wilcox, wilcox[['cluster_negbin']]) 
wilcox <- lapply(wilcox, function(w) dplyr::filter(w, p_val_negbin < 0.05))
topmarkers <- Reduce(rbind, lapply(wilcox, function(top) {
					   top <- dplyr::mutate(top, d = pct.1_wilcox - pct.2_wilcox)
					   head(dplyr::arrange(top, -d), 100)
	})
)

print(table(topmarkers$cluster_negbin))

write.table(topmarkers, paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_negbin_top100_marker_genes.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)

mkrs <- split(topmarkers, topmarkers$cluster_negbin)
mkrs <- Filter(Negate(is.null), mkrs)

if(any(grepl("gene_version", colnames(mkrs[[1]])))) {

	   mkrs <- lapply(mkrs, function(m) head(m$gene_version, 50))

} else {

	   mkrs <- lapply(mkrs, function(m) head(m$gene, 50))
	
}


parallel::mclapply(names(mkrs), function(cls) {

	g <- mkrs[[cls]]
	gnames <- filter(topmarkers, gene%in%g)[['Gene name']]

	ggs <- FeaturePlot(mca, 
			   features = unique(g), 
			   combine=FALSE)

	ggs <- lapply(seq(ggs), function(i) {
			      igene <- names(ggs[[i]]$data)[4]
			      ggs[[i]] + 
				      ggtitle(paste(unique(dplyr::filter(topmarkers,
							    gene_version==igene)[['Gene name']]),
				      collapse = '///')) 
			   })

	title <- ggdraw() + 
		draw_label(
			paste0(jobname, " C", cls),
			fontface = 'bold',
			hjust = 0.5
		) +
	theme(
		plot.margin = margin(0, 0, 0, 0)
  	)

	jpeg(paste0(outputfolder, '/', jobname, "_C", cls, "_res_", gsub("\\.", "", res), "_negbin_marker_genes.jpeg"), 
	    width = 1980, height = 1980, pointsize = 74, quality = 100)

	plot(
	    plot_grid(title,
		      plot_grid(plotlist = ggs, ncol = 5),
		      ncol = 1, rel_heights = c(0.01, 1))
	)

	dev.off()
}, mc.cores = ncores)


# ------------------------------------------------------------------------------
# Marker expression visualization (heatmap)
# ------------------------------------------------------------------------------

if(length(colnames(mca)) > 10000) {
	ncell <- 10000
} else {
	ncell <- length(colnames(mca))
}

cells <- sample(colnames(mca), ncell, replace=FALSE)

print(head(cells))

pdf(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res), "_negbin_marker_genes_heatmap.pdf"), 
     width = 18, height = 15)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_res_", gsub("\\.", "", res),"_negbin_marker_genes_heatmap.jpeg"), 
 	    width = 1980, height = 980, pointsize = 74, quality = 100)
plot(
DoHeatmap(
  mca,
  features = unlist(mkrs),
  cells = cells,
  group.by = "ident",
  group.bar = TRUE,
  group.colors = NULL,
  disp.min = -2.5,
  disp.max = NULL,
  slot = "scale.data",
  assay = "SCT",
  label = TRUE,
  size = 5.5,
  hjust = 0,
  angle = 45,
  raster = TRUE,
  draw.lines = TRUE,
  lines.width = NULL,
  group.bar.height = 0.02,
  combine = TRUE
) + scale_fill_gradientn(colors = c("#92c5de", "#f7f7f7", "#b2182b"))
)
dev.off()
