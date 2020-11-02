#
# Estimate cell type composition difference given a colData variable (groupvar), a cell-ontology and considering the inter-sample variability (Seurat3 objects only)
#

"  Estimate cell type composition difference given a colData variable (groupvar), a cell-ontology and considering the inter-sample variability (Seurat3 objects only)

Usage: diff_cellcomp.R --jobname=<value> [--mcafile=<file>] [--samplevar=<value>] --groupvar=<value> --resolution=<value> --infolder=<folder> --outfolder=<folder>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --mcafile=<file>     rds file with mca (Seurat3) object.
  --groupvar=<value>   Variable to gather cells eg. Condition: WT, KO.
  --samplevar=<value>  Sample variable.
  --resolution=<value> Clustering resolution (numeric) or a string if the name of the cell-ontology variable.
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
	  'Rtsne', 'MASS', 'monocle3', "betareg", "lmtest")

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
groupvar <- arguments$groupvar
resolution <- arguments$resolution

# --- Functions

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# --- Read cell_data_set object with merged samples
if(!is.null(arguments$mcafile)) {

	mca_file <- paste0(inputfolder, '/', arguments$mcafile)
	jobname <- gsub("\\.rds", "", arguments$mcafile)
	
} else {

	mca_file <- paste0(outputfolder, '/', jobname, '_seurat3_merged_clustered_resrange.rds')
	
}

if(!is.null(arguments$auxlabel)) {

	auxlabel <- arguments$auxlabel
	
} else {

	auxlabel <- NULL

}

# --- Read cell_data_set object with merged samples

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

cond <- DimPlot(mca, 
	  reduction = "umap", 
	  pt.size = 0.01, 
	  split.by = groupvar) + 
		ggplot2::theme(legend.position = "bottom") + 
		ggtitle(groupvar)

cond2 <- DimPlot(mca, 
	  reduction = "umap", 
	  pt.size = 0.01, 
	  group.by = groupvar) +
		ggplot2::theme(legend.position = "bottom") + 
		ggtitle(groupvar)

# --- title
title <- ggdraw() + 
  draw_label(
    paste0(jobname, ' ', groupvar),
    fontface = 'bold',
    hjust = 0.5
  ) +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

jpeg(paste0(outputfolder, '/', jobname, "_", "pickedup_clusters_umap.jpeg"), 
		    width = 1380, height = 900, pointsize = 74, quality = 100)
    plot_grid(title, plot_grid(cond, cond2, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

#
meta <- mca@meta.data
meta[['barcode']] <- rownames(meta)
#
print(head(meta, 2))
samples_levels <- unique(meta[[samplevar]]) 
print(samples_levels)

umap_dat <- data.frame(mca@reductions$umap@cell.embeddings)
dim(umap_dat)
umap_dat[['barcode']] <- rownames(umap_dat)
#
umap_dat <- merge(umap_dat, meta, by = "barcode")
rownames(umap_dat) <- umap_dat[['barcode']]
print(head(umap_dat[, 1:10], 2))
#
groupvar_levels <- unique(umap_dat[[groupvar]]) 
print(groupvar_levels)
#
samples_levels <- unique(umap_dat[[samplevar]]) 
print(samples_levels)

dens_plots <- Reduce(rbind,
		     lapply(groupvar_levels, function(l) {
			print(l)
			sub_umap_mat <- filter(umap_dat, !!rlang::sym(groupvar) == l)
			print(nrow(sub_umap_mat))
			l_cells <- sub_umap_mat[['barcode']]
			x <- umap_dat[l_cells, 'UMAP_1']
			y <- umap_dat[l_cells, 'UMAP_2']
			d <- get_density(y, x, n = 100)
			sub_umap_mat[['dens_2d_grvar']] <- d
			sub_umap_mat
		})
	)

g <- ggplot(dens_plots, aes(x = UMAP_1, y = UMAP_2)) +
	geom_point(aes(color = dens_2d_grvar), 
		   size = 0.01) +
	theme_classic() +
	scale_color_gradient2(low = "grey", 
			      high = "#005a32") +
	facet_wrap(as.formula(paste0("~", groupvar)))


pdf(paste0(outputfolder, '/', jobname, "_diff_density_", groupvar, ".pdf"),
    width = 8, height = 4)
	plot(g)
dev.off()

jpeg(paste0(outputfolder, '/', jobname, "_diff_density_", groupvar, ".jpeg"), 
	    width = 980, height = 440, pointsize = 74, quality = 100)
	plot(g)
dev.off()

# --------------------------------------------------
# Differential cell composition in umap space
# --------------------------------------------------

dens_plots[['2d_dens_all']] <- get_density(dens_plots[['UMAP_1']], dens_plots[['UMAP_2']], n = 100)
dens_plots[['diff_dens']] <-  dens_plots[['dens_2d_grvar']] - dens_plots[['2d_dens_all']]
dens_plots[['diff_dens_log']] <-  log2(dens_plots[['dens_2d_grvar']] / dens_plots[['2d_dens_all']])

head(dens_plots,2)

g <- ggplot(dens_plots, aes(x = UMAP_1, y = UMAP_2)) +
	geom_point(aes(color = diff_dens_log), 
		   size = 1,
		   alpha = 0.7) +
	theme_classic() +
	scale_color_gradient2(low = "#7ca2c7", 
			      mid = "#cccccc",
			      high = "#d24644",
			      midpoint = 0,
			      limits = c(-2, 2)) + 
	facet_wrap(as.formula(paste0("~", groupvar)))

pdf(paste0(outputfolder, '/', jobname, "_diff_density.pdf"),
    width = 25, height = 12)
	plot(g)
dev.off()


jpeg(paste0(outputfolder, '/', jobname, "_diff_density.jpeg"), 
	    width = 1500, height = 700, pointsize = 74, quality = 100)
	plot(g)
dev.off()

# ---------------------------------------------------------------------
# cell proportion differences among the levels of the grouping variable
# ---------------------------------------------------------------------

celltypes <- split(meta, meta[[groupvar]])

celltypes <- Reduce(rbind,
       lapply(names(celltypes), function(nsc) { 
			  celltypes[[nsc]] %>%
			  group_by(!!rlang::sym(cluster_var)) %>%
			  summarise(p = n()/nrow(celltypes[[nsc]]),
				    condition = nsc)
			  })
       )

g <- ggplot(celltypes, aes(x = !!rlang::sym(cluster_var), y = p, fill = !!rlang::sym(groupvar))) +
	geom_bar(stat = "identity", position = "dodge") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90))


pdf(paste0(outputfolder, '/', jobname, '_cell_types_propportion', cluster_var, '_', groupvar ,'barplot.pdf'))
	plot(g)
dev.off()

g <- ggplot(celltypes, aes(x = !!rlang::sym(cluster_var), y = p, fill = !!rlang::sym(groupvar))) +
	geom_bar(stat = "identity", position = "stack") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90)) 

pdf(paste0(outputfolder, '/', jobname, '_cell_types_propportion', cluster_var, '_', groupvar ,'barplot.pdf'))
	plot(g)
dev.off()

# ---------------------------------------------------------------------
# Cell type proportion per sample
# ---------------------------------------------------------------------

sample_cell <- split(meta, meta[[samplevar]])

sample_cell <- Reduce(rbind,
       lapply(names(sample_cell), function(nsc) { 
			  sample_cell[[nsc]] %>%
			  group_by(!!rlang::sym(cluster_var)) %>%
			  summarise(p = n()/nrow(sample_cell[[nsc]]),
				    sample = nsc)
			  })
       )

g <- ggplot(sample_cell, aes(x = !!rlang::sym(samplevar), y = p, fill = !!rlang::sym(cluster_var))) +
	geom_bar(stat = "identity") +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90))

pdf(paste0(outputfolder, '/', jobname, '_cell_per_cluster_per_sample_cell_clusters_barplot.pdf'))
	plot(g)
dev.off()

# ---------------------------------------------------------------------
# Cell proportions per class jitter 
# ---------------------------------------------------------------------

print(head(sample_cell, 2))

meta %>%
	group_by(!!rlang::sym(samplevar)) %>%
	summarize(condition = unique(!!rlang::sym(groupvar))) -> sample_condition

print(head(sample_condition, 2))

sample_cell <- merge(sample_cell, sample_condition,  by = samplevar)

print(head(sample_cell, 2))

g <- ggplot(sample_cell, aes(x = !!rlang::sym(groupvar), 
			     y = p, 
			     color = !!rlang::sym(groupvar), 
			     fill = !!rlang::sym(groupvar))) +
	geom_boxplot(alpha = 0.5, size = .1, 
		     outlier.shape = NA) +
	geom_jitter(alpha = 1, size = .9) +
	theme_classic() +
	facet_wrap(as.formula(paste0('~', !!rlang::sym(groupvar))), nrow = 1, scale = "free")

pdf(paste0(outputfolder, '/', jobname, '_cell_per_cluster_per_sample_cell_clusters_jitter_scalefree.pdf'), 
    height = 2, width = 10)
	plot(g)
dev.off()

# -----------------------------------------------------------
# t-test to compare the distributions
# -----------------------------------------------------------

ttest <- list()

for(i in names(celln)) {  
	print(paste('C', i))
	ttest[[i]] <- t.test(formula=p~condition, data=filter(sample_cell, cell_ontology == i))
}

tt_res <- do.call(rbind,
		  lapply(names(ttest), function(tt_name) {
			tt <- ttest[[tt_name]]
			data.frame('mean_control' = tt$estimate[1],
				   'mean_PD' = tt$estimate[2],
				   'statistic' = tt$statistic,
				   'parameter' = tt$parameter,
				   'p.value' = tt$p.value,
				   'CIL' = tt$conf.int[1],
				   'CIR' = tt$conf.int[2],
				   'null.value' = tt$null.value,
				   'cell_type' = tt_name)
			})
		  )

print(dim(tt_res))
print(head(tt_res, 2))

write.table(tt_res, paste0(outputfolder, '/', jobname, '_ttest_cell_type_percentage_differences.tsv'),
	    quote = FALSE, row.names = FALSE, sep = "\t")

# -----------------------------------------------------------
# Modelling proportions based on clinical features
# -----------------------------------------------------------

sample_cell <- split(meta, meta[[samplevar]])

sample_cell <- Reduce(rbind,
       lapply(names(sample_cell), function(nsc) { 
			  sample_cell[[nsc]] %>%
			  group_by(!!rlang::sym(cluster_var)) %>%
			  summarise(p = n()/nrow(sample_cell[[nsc]]),
				    sample = nsc)
			  })
       )

print(head(sample_cell, 2))

meta_sample_curated <- data.table::fread("./data/sample_metadata_curated.tsv")
meta_sample_curated <- merge(meta_sample_curated, sample_cell, by = samplevar)

meta %>%
	group_by(!!rlang::sym(samplevar))) %>%
	summarise('ncells' = n()) -> ncells

print(head(ncells))

meta_sample_curated <- merge(meta_sample_curated, ncells, by = samplevar)

print(head(meta_sample_curated, 2))
print(dim(meta_sample_curated))

model2df <- function(modelres) {
	m <- data.frame(coef(summary(modelres))$mean)
	m[["covar"]] <- rownames(m)
	m[["model_precision_Estimate"]] <- coef(summary(modelres))$precision[1]
	m[["model_precision_Std_Error"]] <- coef(summary(modelres))$precision[2]
	m[["model_precision_z_value"]] <- coef(summary(modelres))$precision[3]
	m[["model_precision_pr_higherthan_z"]] <- coef(summary(modelres))$precision[4]
	m
}

print(head(meta_sample_curated,2))

cellprop <- model2df(betareg(`p` ~ cell_ontology * condition + as.numeric(ncells) * as.numeric(pmi) * sex * as.numeric(age_at_death) , data = meta_sample_curated))
cellprop

write.table(cellprop, paste0(outputfolder, '/', jobname, '_beta_regression_all_cells_results.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE)

cellprop <- data.table::fread(paste0(outputfolder, '/', jobname, '_beta_regression_all_cells_results.tsv'))
print(head(cellprop))

betar <- lapply(unique(meta[[clustervar]]), function(cell_type) {

 	print(cell_type)

	sub_meta_sample <- filter(meta_sample_curated, cell_ontology == cell_type)
	
	print(dim(sub_meta_sample))
	print(head(sub_meta_sample))

	cellprop <- model2df(betareg(`p` ~ condition, data = sub_meta_sample))
	cellprop[["celltype"]] <- cell_type
	head(cellprop)

	cellprop1 <- model2df(betareg(`p` ~ condition + ncells, data = sub_meta_sample))
	cellprop1[["celltype"]] <- cell_type
	head(cellprop1)

	if(cell_type == "CADPS2+ neurons") {
		
		cellprop2 <- model2df(betareg(`p` ~ condition + ncells + pmi + age_at_death, data = sub_meta_sample))
		cellprop2[["celltype"]] <- cell_type
		head(cellprop2)

	} else {

		cellprop2 <- model2df(betareg(`p` ~ condition + ncells + pmi + age_at_death + sex, data = sub_meta_sample))
		cellprop2[["celltype"]] <- cell_type
		head(cellprop2)

	}
	
	list('covs' = cellprop,
	     'covs1' = cellprop1,
	     'covs2' = cellprop2)
})

# Naive model p ~ condition
naive <- do.call(rbind, lapply(betar, function(b) b[['covs']]))
print(dim(naive))
write.table(naive, paste0(outputfolder, '/', jobname, '_beta_regression_naive_results.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE)

# Complete model p ~ condition + ncells + pmi + age_at_death + sex
complete <- do.call(rbind, lapply(betar, function(b) b[['covs2']]))
print(dim(complete))
write.table(complete, paste0(outputfolder, '/', jobname, '_beta_regression_complete_results.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE)


# Visualization of the resulting coefficients
print(head(complete, 2))

# Split by cell type
gg <- ggplot(filter(complete, covar != "(Intercept)"), aes(x = covar, y = Estimate, colour = celltype)) +
	geom_point(size = 1.6) +
	geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error), width = 0.0, lwd = 0.8) +
	geom_hline(yintercept = 0, color='grey', size=0.05) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	scale_color_manual(values = mycolors) +
	scale_x_discrete(limits = c("conditionIPD", "pmi", "age_at_death", "sexM", "ncells")) +
	facet_wrap(~celltype, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_model_residuals_by_celltype.pdf"),
    width = 4.606*2, height = 1.3*2)
	plot(gg)
dev.off()

# Split by cofactor
gg <- ggplot(filter(complete, covar != "(Intercept)"), aes(x = celltype, y = Estimate, colour = celltype)) +
	geom_hline(yintercept = 0, color='grey', size=0.1) +
	geom_point() +
	geom_errorbar(aes(ymin = Estimate - Std..Error, ymax = Estimate + Std..Error)) +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
	scale_color_manual(values = mycolors) +
	facet_wrap(~covar, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_model_residuals_by_covar.pdf"),
    width = 4.706*2, height = 1.719*2)
	plot(gg)
dev.off()
