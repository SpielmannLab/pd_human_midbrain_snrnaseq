#
# Filter cell_data_set objects (Monocle3) given the run_scrublet.py output
#

"Filter cell_data_set objects (Monocle3) given the run_scrublet.py output

Usage: filter_duplets.R --scobject=<file> --dupthr=<value> --infolder=<folder> --outfolder=<folder> --dataset=<value>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --scobject=<file>    *.rds file SingleCellObject.
  --dupthr=<value>     Duplet score threshold to retain singlets (score <= dupthr).
  --infolder=<file>    Path to the single_cell_data .rds file.
  --outfolder=<file>   Path to results folder.
  --dataset=<value>    Either RAW or FILTERED.

"-> doc

library(docopt)
arguments <- docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools', 'ggplot2',
          'cowplot', "dplyr", "ggwordcloud", "monocle3")

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
dupthr <- as.numeric(arguments$dupthr)
inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
ds <- arguments$dataset

if(ds == "RAW") {

	jobname <- gsub("\\.rds", "", rdsfile)
	
} else if(ds == "FILTERED") {

	jobname <- gsub("\\.rds|_filtered", "", rdsfile)
	
} else {

	stop(paste("Sorry, datasets", ds, "are not supported. Try",
		 "FILTERED or RAW."))
}

# Read single_cell_object
cds <- readRDS(paste0(inputfolder, '/', rdsfile))
if(class(cds) != "cell_data_set") {
    stop(paste("Sorry", print(class(cds)), "objects are not supported!",
	       "Try cell_data_set (Monocle3) instead!"))
}
print("cds dimensions genes x cells:")
print(dim(cds))

# Read scrublet output
# UMAP
umap_file <- paste0(inputfolder, '/', jobname, '_umap_scrublet.csv')
umap <- data.table::fread(umap_file)
colnames(umap) <- c("UMAP_1", "UMAP_2")
print("UMAP sapce from scrublet:")
print(head(umap))
print(dim(umap))

# Duplet observed-score
dup_obs_file <- paste0(inputfolder, '/', jobname, '_duplets_score.csv')
dup_obs <- data.table::fread(dup_obs_file)

g_dup <- ggplot(dup_obs, aes(V1)) +
	geom_histogram() +
	geom_vline(xintercept = dupthr) + 
	ggtitle("observed duplets") +
	theme_classic()

umap <- cbind(umap, dup_obs)
colnames(umap)[3] <- "obs_dup"

# Duplet simulated-score
dup_obs_file <- paste0(inputfolder, '/', jobname, '_sim_duplets_score.csv')
dup_obs <- data.table::fread(dup_obs_file)

g_sim_dup <- ggplot(dup_obs, aes(V1)) +
	geom_histogram() +
	geom_vline(xintercept = dupthr) + 
	ggtitle("simulated duplets") +
	theme_classic()

umap %>%
	dplyr::mutate(dup = ifelse(obs_dup >= dupthr, "duplet", "singlet")) -> umap 

g_umap <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2)) +
	geom_point(size = 0.8, aes(color = obs_dup), alpha = 0.4) +
	scale_colour_gradient(low = "grey", high = "#67001f") +
	theme_classic()

g_umap_dup <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2)) +
	geom_point(aes(color = dup)) +
	ggtitle(paste(table(umap$dup), collapse = "/")) +
	theme_classic()

# --- Write QC plots
title <- ggdraw() + 
  draw_label(
    paste("Histogram Scrublet scores", jobname),
    fontface = 'bold',
    #x = 0,
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0)
  )

qc_file <- paste0(outputfolder, '/', jobname, '_duplet_score.pdf')

pdf(qc_file, width = 11.69, height = 8.27)
plot(
    plot_grid(title, plot_grid(g_dup, g_sim_dup, g_umap, g_umap_dup,
			       ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
    )
dev.off()

sub_cds <- cds[, which(umap$dup == "singlet")]
colData(sub_cds) <- cbind(colData(sub_cds), dplyr::filter(umap, dup == "singlet"))

# --- Write cell_data_set object with QC metrics

cds_file <- paste0(outputfolder, '/', jobname, '_filtered.rds')
saveRDS(sub_cds, file = cds_file)

message(paste("The cell_data_set object (Monocle3)", 
	      cds_file, "has been created."))
