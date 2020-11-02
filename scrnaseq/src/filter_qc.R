#
# Read and filter a cell_data_set object (Monocle3)
#

"Read and filter a cell_data_set object (Monocle3) 

Usage: filter_qc.R --scobject=<file> --specie=<value> --infolder=<folder> --outfolder=<folder> --minfeat=<value> --maxfeat=<value> --mincellfeat=<value> --maxcellfeat=<value> --mincountcell=<value> --maxcountcell=<value> --mingenecell=<value> --maxgenecell=<value> --pctmt=<value> --pctrb=<value> --mt=<logical> --rb=<logical> [--sys.rt=<logical>]

Options:
  -h --help              Show this screen.
  --version              00.99.01
  --scobject=<file>      *.rds file cell_data_set (Monocle3).
  --specie=<value>       Either 'hg' or 'mm'.
  --infolder=<file>      Path to the single_cell_data .rds file.
  --outfolder=<file>     Path to results folder.
  --minfeat=<value>      Minimum feature count-depth.
  --maxfeat=<value>      Maximum feature count-depth.
  --mincellfeat=<value>	 Minimum number of feature cell-depth.
  --maxcellfeat=<value>  Maximum number of feature cell-depth. 
  --mincountcell=<value> Minimum cell count-depth.
  --maxcountcell=<value> Maximum cell count-depth.
  --mingenecell=<value>  Minimum cell gene-depth.
  --maxgenecell=<value>  Maximum cell gene-depth.
  --pctmt=<value>        Maximum mitochondrial-count/total-count/cell rate allowed.
  --pctrb=<value>        Maximum ribosomal-count/total-count/cell rate allowed.
  --mt=<logical>         TRUE for filtering mitochondrial genes out.
  --rb=<logical>         TRUE for filtering ribosomal genes out.
  --sys.rt=<logical>     Optional. Hardcode to remove cells with diff. RT reaction efficiency.

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "monocle3")

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

# ---------------
# --- Parameters
# ---------------

rdsfile <- arguments$scobject
sp <- arguments$specie
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
minfeat <- as.numeric(arguments$minfeat)
maxfeat <- as.numeric(arguments$maxfeat)
mincellfeat <- as.numeric(arguments$mincellfeat)
maxcellfeat <- as.numeric(arguments$maxcellfeat)
mincountcell <- as.numeric(arguments$mincountcell)
maxcountcell <- as.numeric(arguments$maxcountcell)
mingenecell <- as.numeric(arguments$mingenecell)
maxgenecell <- as.numeric(arguments$maxgenecell)
pctmt <- as.numeric(arguments$pctmt)
pctrb <- as.numeric(arguments$pctrb)
mt <- as.logical(arguments$mt)
rb <- as.logical(arguments$rb)
jobname <- gsub("_RAW_qc.rds|_FILTERED_qc.rds", "", rdsfile)

if(is.null(arguments$sys.rt)) {
	sys.rt <- FALSE
} else {
	sys.rt <- as.logical(arguments$sys.rt)
}

# Read single_cell_object
cds <- readRDS(paste0(inputfolder, '/', rdsfile))

if(class(cds) != "cell_data_set") {
    stop(paste("Sorry", print(class(cds)), "objects are not supported!",
	       "Try cell_data_set (Monocle3) instead!"))
}


features <- data.frame(fData(cds)) %>%
	dplyr::select(id) %>%
	unlist %>%
	as.character

cells <- data.frame(pData(cds)) %>%
	dplyr::select(barcode) %>%
	unlist %>%
	as.character

if(mt) {

	if(sp == "mm") {

		# Mitochondrial genes
		mt_genes <- fData(cds)[grepl("^mt-", fData(cds)$gene_short_name), ]$id
		message("Mitochondrial genes detected:")
		print(fData(cds)[grepl("^mt-", fData(cds)$gene_short_name), 
		      ]$gene_short_name)
		if(length(mt_genes) == 0) {
			message("No mitochondrial genes detected to filter out")
		} else {
			features <- features[!features%in%mt_genes]
		}
	
	} else if (sp == "hg") {

		# Mitochondrial genes
		mt_genes <- fData(cds)[grepl("^MT-", fData(cds)$gene_short_name), ]$id
		message("Mitochondrial genes detected:")
		print(fData(cds)[grepl("^MT-", fData(cds)$gene_short_name), 
		      ]$gene_short_name)
		if(length(mt_genes) == 0) {
			message("No mitochondrial genes detected to filter out")
		} else {
			features <- features[!features%in%mt_genes]
		}

	} else {
	
		stop("Please provide a valid specie parameter")
	
	}
}

if(rb) {

	if(sp == "mm") {

	# Ribosomal genes
	rb_genes <- fData(cds)[grepl("^Rps|^Rpl", fData(cds)$gene_short_name), ]$id
	message("Ribosomal genes detected:")
	print(fData(cds)[grepl("^Rps|^Rpl", fData(cds)$gene_short_name), 
	      ]$gene_short_name)
	if(length(rb_genes) == 0) {
		message("No ribosomal genes detected to filter out")
	} else {
		features <- features[!features%in%rb_genes]
	}

} else if (sp == "hg") {

	# Ribosomal genes
	rb_genes <- fData(cds)[grepl("^RPS|^RPL", fData(cds)$gene_short_name), ]$id
	message("Ribosomal genes detected:")
	print(fData(cds)[grepl("^RPS|^RPL", fData(cds)$gene_short_name), 
	      ]$gene_short_name)
	if(length(rb_genes) == 0) {
		message("No ribosomal genes detected to filter out")
	} else {
			features <- features[!features%in%rb_genes]
		}

	} else {
	
		stop("Please provide a valid specie parameter")

	}
}

cds <- cds[features, cells]

features <- data.frame(fData(cds)) %>%
	dplyr::filter(as.numeric(feat_count_depth) > minfeat & as.numeric(feat_count_depth) < maxfeat & as.numeric(feat_cell_depth) > mincellfeat & as.numeric(feat_cell_depth) < maxcellfeat) %>%
	dplyr::select(id) %>%
	unlist %>%
	as.character

cells <- data.frame(pData(cds)) %>%
	dplyr::filter(as.numeric(cell_depth_counts) > mincountcell & as.numeric(cell_depth_counts) < maxcountcell & as.numeric(cell_depth_genes_detected) > mingenecell & as.numeric(cell_depth_genes_detected) < maxgenecell & as.numeric(pct_mt_counts) < pctmt & as.numeric(pct_rb_counts) < pctrb) %>%
	dplyr::select(barcode) %>%
	unlist %>%
	as.character

if(sys.rt) {

	cells <- cells[cells %in% keep_cells]

}

sub_cds <- cds[features, cells]

report <- data.frame('cds' = c('RAW', 'FILTER'),
		     'n_features' = c(nrow(cds), nrow(sub_cds)),
		     'n_cells' = c(ncol(cds), ncol(sub_cds)))

g1 <- ggplot(report, aes(x = cds, y = n_features)) +
	geom_bar(stat = "identity", alpha = 0.3) +
	geom_text(aes(label = n_features), vjust = -1) +
	theme_classic() +
	ggtitle("Features")

g2 <- ggplot(report, aes(x = cds, y = n_cells)) +
	geom_bar(stat = "identity", color = "#a6bddb", fill = "#a6bddb") +
	geom_text(aes(label = n_cells), vjust = -1) +
	theme_classic() +
	ggtitle("Cells")

# --- Write filter report plots
title <- ggdraw() + 
  draw_label(
    paste0(jobname, " filter report (", Sys.time(), ")"),
    fontface = 'bold',
    #x = 0,
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0)
  )


filter_file <- paste0(outputfolder, '/', jobname, '_', 'filter_report.pdf')

pdf(filter_file, height = 5.85, width = 8.27)
    plot_grid(title, plot_grid(g1, g2, ncol = 2), 
	      ncol = 1, rel_heights = c(0.06, 1))
dev.off()

message(paste("The filter report ", filter_file, "has been created."))

# --- Write the filtered cell_data_set object
cds_file <- paste0(outputfolder, '/', jobname, '_filtered.rds')
saveRDS(sub_cds, file = cds_file)

# -------------------------------------------
# --- Files for duplets ---------------------
# --- Write matrix for duplets identification
# -------------------------------------------

m <- exprs(sub_cds)#@assays@data$counts
Matrix::writeMM(m, gsub("rds", "mtx", cds_file))
# --- Write genes for duplets identification pipeline
features <- data.frame(Reduce(cbind, sub_cds@rowRanges@elementMetadata@listData))
write.table(features, paste0(outputfolder, '/', jobname, '_', 'filtered_features_metadata.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# --- Write barcodes for duplets identification pipeline
barcodes <- colData(sub_cds)
write.table(barcodes, paste0(outputfolder, '/', jobname, '_', 'filtered_barcodes_metadata.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message(paste("The filtered cell_data_set object (Monocle3)", 
      cds_file, "has been created."))
