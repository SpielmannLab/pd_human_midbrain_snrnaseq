#
# Estimate expression difference given a colData variable and a clustering resolution (Seurat3 objects only)
#

" Estimate expression difference given a colData variable and a clustering resolution (Seurat3 objects only)

Usage: diff_exp_seurat2monocle.R --jobname=<value> --groupvar=<value> --samplevar=<value> --resolution=<value> --specie=<value> --infolder=<folder> --outfolder=<folder>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --jobname=<file>     Descriptive name for your experiment.
  --groupvar=<value>   Variable to do comparison.
  --samplevar=<value>  Sample variable. Ref for pseudo-bulk.
  --resolution=<value> Clustering resolution.
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
	  'Rtsne', 'MASS', 'htmlwidgets')

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

# Install last version of both (Beta version keep changing)

#devtools::install_github('cole-trapnell-lab/leidenbase')
#devtools::install_github('cole-trapnell-lab/monocle3', force=TRUE)


# -------------------------------
suppressMessages(library(monocle3))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(Seurat))
#suppressMessages(library(batchelor))
suppressMessages(library(sctransform))
suppressMessages(library(future))
suppressMessages(library(MASS))
suppressMessages(library(zinbwave))
suppressMessages(library(DESeq2))
suppressMessages(library(ggfortify))
suppressMessages(library(topconfects))
suppressMessages(library(plotly))
suppressMessages(library(edgeR))

#suppressMessages(library(dplyr))
#suppressMessages(library(ggwordcloud))
# ------------------------------


outputfolder <- "/project/zvi/cprada/hbrain_project/results/PD_CO_midbrain/PD_CO_hg_brain_seurat3_integrated_npcs_25_nvgs_4000_full_cca_seurat3_integrated_norm_seurat3_merged_clustered_resrange"
jobname <- "PD_CO_hg_brain"
ncores <- 5
groupvar <- "condition"
res <- 0.1
samplevar <- "sample"
specie <- "hg"

jobname <- arguments$jobname
groupvar <- arguments$groupvar
samplevar <- arguments$samplevar
res <- arguments$resolution
specie <- arguments$specie
inputfolder <- arguments$infolder
outputfolder <- arguments$outfolder

# --- Read cell_data_set object with merged samples
inputfolder <- "/project/zvi/cprada/hbrain_project/results/PD_CO_midbrain"
mca_file <- paste0(inputfolder, '/PD_CO_hg_brain_seurat3_integrated_npcs_25_nvgs_4000_full_cca_seurat3_integrated_norm_seurat3_merged_clustered_resrange.rds')
inputfolder <- "/project/zvi/cprada/hbrain_project/results/PD_CO_midbrain/PD_CO_hg_brain_seurat3_integrated_npcs_25_nvgs_4000_full_cca_seurat3_integrated_norm_seurat3_merged_clustered_resrange"

mca <- readRDS(mca_file)

message(paste("The cell_data_set object (Seurat3)", 
	      mca_file, "has been read."))

# Setting default cell identity based on clustering resolution


ident_col <- "cell_ontology"
print(table(mca@meta.data[[ident_col]]))
Idents(mca) <- ident_col

# Checking cell metadata variable
print(table(mca@meta.data[[groupvar]]))

# Checking sample metadata
print(table(mca@meta.data[[samplevar]]))

# Gene annotation
if(specie == "hg") {

	refbiomart <- data.table::fread("/project/zvi/cprada/hbrain_project/data/ensg_biomart_20191206.tsv")

} else if(specie == "mm"){

	refbiomart <- data.table::fread("/project/zvi/cprada/limb_project/data/biomart_mm_20191206.tsv")

} else {

	stop("Sorry this specie isnt supported, try mm or hg please")

}


# Counts
raw_eset <- mca@assays$RNA@counts 
print(raw_eset[1:5, 1:10])
print(dim(raw_eset))

row_df <- mca@assays$RNA@meta.features
row_df[["gene"]] <- rownames(row_df)
dim(row_df)
row_df <- merge(row_df, refbiomart, by.x ="gene", by.y = "Gene stable ID", all.x = TRUE) 
row_df <- filter(row_df, !duplicated(`gene`))
dim(row_df)
colnames(row_df)[7] <- "gene_short_name"
rownames(row_df) <- row_df[["gene"]]
head(row_df)

lapply(unique(mca@meta.data[[ident_col]]), function(cell) {

	       cells <- rownames(mca@meta.data)[which(mca@meta.data[[ident_col]] == cell)]
	       sub_mca <- mca[, cells]
	       seurat2monocleDE(sub_mca, paste0(gsub("\\+| ", "_", cell), '_', "DEG"))

	  })

seurat2monocleDE <- function(mca, jobname) {
	
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
	det_thr <- 0.2

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
				 model_formula_str = "~condition",
				 expression_family = "quasipoisson",
				 cores = 10)

	fit_coefs <- coefficient_table(cgene_fits)

	colnames(fit_coefs)
	table(fit_coefs$term)
	unique(fit_coefs$term)

	coef_res <- fit_coefs %>% 
	       filter(term == "conditionIPD") %>% 
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

# Sample metadata
meta <- mca@meta.data
print(head(meta,2))
meta %>%
	mutate('sample' = gsub("_filtered_FILTERED_qc_seurat3", "", `sample`)) %>%
	group_by(!!rlang::sym(samplevar)) %>%
	summarise('sex' = gsub("_", "", unique(sex.y)),
		  'age' = as.numeric(gsub("_", "", unique(age_at_death))),
		  'condition' = unique(condition)) %>%
	data.frame -> sample_meta

#
# read meta-data from Anne (setting curated meta-data)
#
metadata_file <- "/project/zvi/cprada/hbrain_project/data/sample_metadata_curated.tsv"
sample_meta <- data.table::fread(metadata_file)
head(sample_meta)
#
#

sample_meta <- DataFrame(sample_meta)
sample_meta[[groupvar]] <- as.factor(sample_meta[[groupvar]])
rownames(sample_meta) <- sample_meta[[samplevar]] 
print(head(sample_meta,2))

# Cell metadata
meta <- mca@meta.data
meta[['barcode']] <- rownames(meta)
print(head(meta,2))

# Minimum number of UMI counts to consider in the DE analysis
# (UMI counts per sample) 
mingenecount <- 200

# Top X% perturbed genes
metathr <- 0.05


# Sample pseudo-bulk PCA
# Sum counts per sample
ss <- unique(meta[[samplevar]])

sample_pbulk <- Reduce(cbind,
			lapply(ss, function(s) {
				print(s)
				cells <- rownames(meta)[which(meta[[samplevar]] == s)]
				print(head(cells))
				bulk <- data.frame(rowSums(raw_eset[, cells]))
				#sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
				sname <- s
				colnames(bulk) <- sname
				print(head(bulk, 2))
				bulk
			})
		)
	
print("Genes x samples matrix")
print(head(sample_pbulk))
print(dim(sample_pbulk))

length(which(rowMeans(sample_pbulk) > 50))

cell_exp <- sample_pbulk[which(rowMeans(sample_pbulk) > 50), ]
head(cell_exp)

mm <- sum(cell_exp)

norm_cell_exp <- apply(cell_exp, 2, function(x) log2((x/(mm/1000000))+1))
head(norm_cell_exp)


center_eset <- apply(norm_cell_exp, 1, function(x) (x-mean(x))/sd(x))


vs <- colnames(sample_meta)
pdf(paste0(outputfolder, "/", jobname, "_TPM_norm_metadata_PCA.pdf"), width = 4, height = 3)
for(v in vs) {
	plot(
		autoplot(prcomp(center_eset[sample_meta[['sample']], ]), 
			 data=data.frame(sample_meta), 
			 colour = v,
			 label = TRUE, 
			 label.size = 3) +
	     	theme_classic()
	)
}
dev.off()


# Check sample expression distributions
dist_data <- reshape2::melt(sample_pbulk)

gg_raw_dist <- ggplot(dist_data, aes(`value`)) +
	geom_density(aes(group = `variable`, fill = `variable`), alpha = .4) +
	theme_classic()

gg_raw_log_dist <- ggplot(dist_data, aes(log10(`value` + 1))) +
	geom_density(aes(group = `variable`, fill = `variable`), alpha = .4) +
	theme_classic()

dds <- DESeqDataSetFromMatrix(countData = sample_pbulk[, sample_meta[["sample"]]],
			      colData = sample_meta,
			      design = ~ condition)# + age_at_death)

print(dds)

print(quantile(rowSums(counts(dds))))
keep <- rowSums(counts(dds)) > 50
dds <- dds[keep, ]

#
# Normalization for PCA visualization ---------
#
vsd <- vst(dds, blind=FALSE)
dim(assay(vsd))
head(assay(vsd), 3)
vsd_mat <- assay(vsd)
colnames(vsd_mat) <- colnames(sample_pbulk)
vs <- colnames(sample_meta)
pdf(paste0(outputfolder, "/", jobname, "_DESeq2_vsd_norm_metadata_PCA_sample_pbulk.pdf"))
for(v in vs) {
	plot(
		autoplot(prcomp(t(vsd_mat)[sample_meta[['sample']], ]), 
			 data=data.frame(sample_meta), 
			 colour = v,
			 label = TRUE, 
			 label.size = 3) +
	     	theme_classic()
	)
}
dev.off()

# raw norm per seq.size
dds <- estimateSizeFactors(dds)
counts_norm <- counts(dds, normalized=TRUE)
pdf(paste0(outputfolder, "/", jobname, "_DESeq2_factorsize_norm_metadata_PCA_sample_pbulk.pdf"))
for(v in vs) {
	plot(
		autoplot(prcomp(t(counts_norm)[sample_meta[['sample']], ]), 
			 data=data.frame(sample_meta), 
			 colour = v,
			 label = TRUE, 
			 label.size = 3) +
	     	theme_classic()
	)
}
dev.off()
# --------------------------------


# Check sample distributions of this normalized expression
dist_data <- reshape2::melt(vsd_mat)

gg_vsd_dist <- ggplot(dist_data, aes(`value`)) +
	geom_density(aes(group = `Var2`, fill = `Var2`), alpha = .4) +
	theme_classic()

pdf(paste0(outputfolder, "/", jobname, "_DESeq2_expression_sample_pbulk_distribution_sample_pbulk.pdf"))
	plot(gg_raw_dist)
	plot(gg_raw_log_dist)
	plot(gg_vsd_dist)
dev.off()
	# -----------------------------------------





# Main function for pseudo-bulk differential expression
pbulk_diffexp <- function(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart) { # count matrix plus metadata
	# Sum counts per sample
	ss <- unique(meta[[samplevar]])

	sample_pbulk <- Reduce(cbind,
			lapply(ss, function(s) {
				print(s)
				cells <- rownames(meta)[which(meta[[samplevar]] == s)]
				print(head(cells))
				bulk <- data.frame(rowSums(raw_eset[, cells]))
				#sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
				sname <- s
				colnames(bulk) <- sname
				print(head(bulk, 2))
				bulk
			})
		)
	
	print("Genes x samples matrix")
	print(head(sample_pbulk))
	print(dim(sample_pbulk))

	# Check sample expression distributions
	dist_data <- reshape2::melt(sample_pbulk)

	gg_raw_dist <- ggplot(dist_data, aes(`value`)) +
		geom_density(aes(group = `variable`, fill = `variable`), alpha = .4) +
		theme_classic()

	gg_raw_log_dist <- ggplot(dist_data, aes(log10(`value` + 1))) +
		geom_density(aes(group = `variable`, fill = `variable`), alpha = .4) +
		theme_classic()

	dds <- DESeqDataSetFromMatrix(countData = sample_pbulk[, rownames(sample_meta)],
				      colData = sample_meta,
				      design = ~ condition)# + age_at_death)

	print(dds)

	print(quantile(rowSums(counts(dds))))
	keep <- rowSums(counts(dds)) > mingenecount
	dds <- dds[keep, ]

	#
	# Normalization for PCA visualization ---------
	#
	vsd <- vst(dds, blind=FALSE)
	dim(assay(vsd))
	head(assay(vsd), 3)
	vsd_mat <- assay(vsd)
	colnames(vsd_mat) <- colnames(sample_pbulk)
	vs <- colnames(sample_meta)
	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_vsd_norm_metadata_PCA.pdf"))
	for(v in vs) {
		plot(
			autoplot(prcomp(t(vsd_mat)[sample_meta[['sample']], ]), 
				 data=data.frame(sample_meta), 
				 colour = v,
				 label = TRUE, 
				 label.size = 3) +
		     	theme_classic()
		)
	}
	dev.off()

	# raw norm per seq.size
	dds <- estimateSizeFactors(dds)
	counts_norm <- counts(dds, normalized=TRUE)
	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_factorsize_norm_metadata_PCA.pdf"))
	for(v in vs) {
		plot(
			autoplot(prcomp(t(counts_norm)[sample_meta[['sample']], ]), 
				 data=data.frame(sample_meta), 
				 colour = v,
				 label = TRUE, 
				 label.size = 3) +
		     	theme_classic()
		)
	}
	dev.off()
	# --------------------------------



	# Check sample distributions of this normalized expression
	dist_data <- reshape2::melt(vsd_mat)

	gg_vsd_dist <- ggplot(dist_data, aes(`value`)) +
		geom_density(aes(group = `Var2`, fill = `Var2`), alpha = .4) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_expression_sample_pbulk_distribution.pdf"))
		plot(gg_raw_dist)
		plot(gg_raw_log_dist)
		plot(gg_vsd_dist)
	dev.off()
	# -----------------------------------------


	dds <- DESeq(object = dds, 
		     test = "LRT",
		     sfType="poscounts",
		     minmu=1e-6,
		     useT=TRUE,
		     quiet = TRUE,
		     parallel = FALSE,
		     reduced = ~1,
		     BPPARAM = bpparam())

	res <- results(dds, contrast=c("condition", "IPD", "Control"))

	res[['gene']] <- rownames(res)

	res <- data.frame(res)

	# Check if stable gene ID was provided with or without version
	g2test <- res[['gene']][1]

	if(grepl("\\.\\d+$", g2test)) {
		res[['gene_version']] <- res[['gene']]
		res[['gene']] <- gsub("\\.\\d+$", "", res[['gene_version']])
	} else {
		res[['gene_version']] <- res[['gene']]
	}

	res <- merge(res, refbiomart, by.x = "gene", by.y = "Gene stable ID") 
	
	res[['index']] <- seq(nrow(res))

	confects <- normal_confects(as.numeric(res$log2FoldChange),
				    se = as.numeric(res$lfcSE), 
				    fdr = 0.05, 
				    full = TRUE)

	print(head(confects$table, 3))

	res <- merge(res, 
		     dplyr::select(confects$table, c(index, `rank`)),
		     by = 'index', all = TRUE)


	irank <- quantile(res[['rank']], metathr)
	print(irank)

	res %>%
		dplyr::mutate(perturbation = ifelse(`rank` <= irank, 
						    ifelse(log2FoldChange < 0, "Down", "Up"),
						    "Unperturbed")) -> res 

	write.table(res, paste0(outputfolder, "/", jobname, "_deg_DESeq2_pbulk_results.tsv"),
		    sep = "\t", quote = FALSE, row.names = FALSE)

	# Volcano plot
	gg <- ggplot(dplyr::filter(res, !is.na(padj)), 
		     aes(x = log2FoldChange, y = -log10(padj), 
			 text = `Gene name`,
			 color = perturbation)) +
		geom_point(alpha = 0.5) +
		geom_errorbarh(aes(xmin = log2FoldChange - lfcSE, 
				   xmax = log2FoldChange + lfcSE), 
			       alpha = .2) +
		scale_color_manual(values = c("blue", "grey", "red")) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_samplewise_pbulk_volcano.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_DESeq2_volcanoplot.html"))
	# -

	# MA plot
	gg <- ggplot(dplyr::filter(res, !is.na(padj)), 
		     aes(x = baseMean, y = log2FoldChange, 
			 text = `Gene name`,
			 color = perturbation)) +
		geom_point(alpha = 0.5) +
		geom_errorbar(aes(ymin = log2FoldChange - lfcSE, 
				  ymax = log2FoldChange + lfcSE), 
			       alpha = .2) +
		scale_color_manual(values = c("blue", "grey", "red")) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_samplewise_pbulk_MAplot.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_DESeq2_MAplot.html"))
	# -

	vis_eset <- reshape2::melt(vsd_mat)
	vis_eset <- merge(vis_eset, sample_meta, 
			  by.x = "Var2", by.y = "sample")

	vis_eset <- merge(vis_eset, 
			  select(res, c(gene_version, `Gene name`)), 
			  by.x = "Var1", by.y = "gene_version")

	top50 <- head(arrange(filter(res, rank <= irank), 
			      rank), 50)[['gene_version']]
	print(head(top50))

	dplyr::filter(data.frame(vis_eset), 
		      Var1 %in% top50) %>%
		ggplot(aes(x = condition, y = value,
			   fill = condition, color = condition)) +
		geom_violin(alpha = .2) +
		geom_jitter(alpha = .5) +
		theme_classic() +
		facet_wrap(~`Gene.name`, nrow = 10) -> ggex

	pdf(paste0(outputfolder, "/", jobname, "_DESeq2_samplewise_pbulk_top_50gene_violin.pdf"),
		width = 20, height = 20)
		plot(ggex)
	dev.off()
}

# all cells 
pbulk_diffexp(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart)

# For each cluster independently
cls <- unique(meta[[ident_col]])

lapply(cls, function(cl) {
	       print(cl)

	       sub_meta <- dplyr::filter(meta, !!rlang::sym(ident_col) == cl)
	       print(dim(sub_meta))

	       cl_cells <- sub_meta[['barcode']]
	       rownames(sub_meta) <- cl_cells
	       print(head(sub_meta, 2))
	       print(dim(sub_meta))
	       sub_raw_eset <- raw_eset[, cl_cells]
	       print(dim(sub_raw_eset))

	       if(ncol(sub_raw_eset) < mingenecount) {
		       NULL
	       } else {

		       sub_jobname <- paste0(jobname, "_cluster_", cl)
		       
		       pbulk_diffexp(sub_raw_eset, sub_meta, sample_meta, sub_jobname, mingenecount, refbiomart)
	       }

})

# ----------

# differential expression using edgeR

pbulk_edger_diffexp <- function(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart) {
	# Sum counts per sample
	ss <- unique(meta[[samplevar]])
	print(ss)

	print(dim(raw_eset))


	sample_pbulk <- Reduce(cbind,
			lapply(ss, function(s) {
				cells <- rownames(meta)[which(meta[[samplevar]] == s)]
				print(head(cells))
				if(length(cells) == 1) {
					bulk <- data.frame(raw_eset[, cells])
				} else {
					bulk <- data.frame(rowSums(raw_eset[, cells]))
				}

				sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
				colnames(bulk) <- sname
				print(head(bulk, 2))
				bulk
			})
		)
	
	print("Genes x samples matrix")
	print(head(sample_pbulk,2))
	print(head(sample_meta))

	print(rownames(sample_meta))
	print(colnames(sample_pbulk))

 	y <- DGEList(counts = sample_pbulk[, rownames(sample_meta)],
				      samples = sample_meta)

	# Checking sample library sizes
	discarded <- scater::isOutlier(y$samples$lib.size, log=TRUE, type="lower")

	if(any(discarded)) {
		print("Notice that these samples were removed from the DE analysis:")
		print(colnames(y)[discarded])
		y <- y[, !discarded]
	} else {
		print("Sample library sizes are more/less homogeneus")
		print(y$samples$lib.size)
	}
	print("Discarded?:")
	summary(discarded)

	keep <- filterByExpr(y, group=sample_meta[[groupvar]],
			     min.total.count = mingenecount)
	table(keep)

	y <- y[keep, ]
	
	y <- calcNormFactors(y)
	y$samples

	# The actual design matrix /should it be customed for every dataset?/
	design <- model.matrix(~ as.numeric(age_at_death) + as.numeric(pmi) + as.factor(sex) + factor(condition),
			       y$samples)

	print(`design`)
	y <- estimateDisp(y, design)
	print(summary(y$trended.dispersion))

	pdf(paste0(outputfolder, "/", jobname, "_biological_component_geneVariance_edgeR.pdf"))
		plotBCV(y)
	dev.off()
	
	fit <- glmQLFit(y, design, robust=TRUE)
	print(summary(fit$var.prior))

	pdf(paste0(outputfolder, "/", jobname, "_likelihood_dispersion_edgeR.pdf"))
		plotQLDisp(fit)
	dev.off()

	res <- glmQLFTest(fit, coef=ncol(design))
	summary(decideTests(res))
	topTags(res)



	results <- glmQLFTest(fit)


	coef_df <- data.frame(results$coefficient)
	coef_df[["gene"]] <- rownames(coef_df)
	print(head(coef_df))
	
	head(arrange(coef_df, -abs(as.numeric.age_at_death.)))

	coef_thr <- 0.1
	head(filter(coef_df, abs(as.numeric.age_at_death.) < 0.1, abs(as.numeric.pmi.) < 0.1, abs(factor.condition.IPD) < 0.1,  abs(as.factor.sex.M) > 1) %>% arrange(-abs(as.factor.sex.M)))
	dim(filter(coef_df, abs(as.numeric.age_at_death.) < 0.1, abs(as.numeric.pmi.) < 0.1, abs(factor.condition.IPD) < 0.1,  abs(as.factor.sex.M) > 1) %>% arrange(-abs(as.factor.sex.M)))

	coef_thr <- 0.1
	head(filter(coef_df, abs(as.numeric.age_at_death.) < 0.1, abs(as.numeric.pmi.) < 0.1, abs(factor.condition.IPD) > 0.1,  abs(as.factor.sex.M) < 0.1) %>% arrange(-abs(factor.condition.IPD)))


	de_res <- data.frame(topTags(res, nrow(res)))
	de_res[['gene']] <- rownames(de_res)

	# Check if stable gene ID was provided with or without version
	g2test <- de_res[['gene']][1]

	if(grepl("\\.\\d+$", g2test)) {
		de_res[['gene_version']] <- de_res[['gene']]
		de_res[['gene']] <- gsub("\\.\\d+$", "", de_res[['gene_version']])
	} else {
		de_res[['gene_version']] <- de_res[['gene']]
	}

	de_res <- merge(de_res, refbiomart, 
			by.x = "gene", by.y = "Gene stable ID") 

	print(head(de_res, 3))

	de_res %>%
		mutate(pertb = abs(logFC)*(-log10(PValue))) %>%
		arrange(-pertb) -> de_res

	de_res[['rank']] <- seq(nrow(de_res))
	print(head(de_res))
	
	metathr <- 0.05
	irank <- de_res[round(nrow(de_res)*metathr), ][['rank']]

	de_res %>%
		dplyr::mutate(perturbation = ifelse(`rank` <= irank, 
						    ifelse(logFC < 0, "Down", "Up"),
						    "Unperturbed")) -> de_res 

	write.table(de_res, paste0(outputfolder, "/", jobname, "_deg_edgeR_pbulk_results.tsv"),
		    sep = "\t", quote = FALSE, row.names = FALSE)

	gg <- ggplot(de_res, 
		     aes(x = logFC, y = -log10(PValue), 
			 text = `Gene name`,
			 color = perturbation)) +
		geom_point(alpha = 0.5) +
		scale_color_manual(values = c("blue", "grey", "red")) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_volcano.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_edgeR_volcanoplot.html"))

	vis_eset <- data.frame(res$fitted.values)
	vis_eset[['gene']] <- rownames(vis_eset)
	vis_eset <- reshape2::melt(vis_eset)
	vis_eset <- merge(vis_eset, sample_meta, 
			  by.x = "variable", by.y = "sample")

	vis_eset <- merge(vis_eset, 
			  select(de_res, c(gene_version, `Gene name`)), 
			  by.x = "gene", by.y = "gene_version")

	top50 <- head(arrange(de_res, rank), 50)[['gene_version']]
	print(head(top50))

	dplyr::filter(data.frame(vis_eset), 
		      gene %in% top50) %>%
		ggplot(aes(x = condition, y = log10(value),
			   fill = condition, color = condition)) +
		geom_violin(alpha = .2) +
		geom_jitter(alpha = .5) +
		theme_classic() +
		facet_wrap(~`Gene.name`, nrow = 10) -> ggex

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_top_50gene_violin.pdf"),
		width = 20, height = 20)
		plot(ggex)
	dev.off()
}

# all cells 
pbulk_edger_diffexp(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart)

# For each cluster independently
cls <- unique(meta[[ident_col]])
cls <- cls[c(9,11)]
mingenecount<-50
lapply(cls, function(cl) {
	       print(cl)

	       sub_meta <- dplyr::filter(meta, !!rlang::sym(ident_col) == cl)
	       print(dim(sub_meta))

	       cl_cells <- sub_meta[['barcode']]
	       rownames(sub_meta) <- cl_cells
	       print(head(sub_meta, 2))
	       print(dim(sub_meta))
	       sub_raw_eset <- raw_eset[, cl_cells]
	       print(dim(sub_raw_eset))
	       head(sub_raw_eset)

	       if(ncol(sub_raw_eset) < mingenecount) {
		       NULL
	       } else {
		       sub_jobname <- paste0(jobname, "_cluster_", cl)

		       pbulk_edger_diffexp(sub_raw_eset, sub_meta, sample_meta[unique(sub_meta[["sample"]]),], sub_jobname, mingenecount, refbiomart)
	       }

	      
})


# Comparison edgeR AND DESeq2
deseq2 <- data.table::fread(paste0(outputfolder, "/", jobname, "_deg_DESeq2_pbulk_results.tsv"))
edger <- data.table::fread(paste0(outputfolder, "/", jobname, "_deg_edgeR_pbulk_results.tsv"))

r <- merge(deseq2, edger, by = 'gene_version')
print(head(r))

r %>%
	mutate(perturbation = ifelse(`perturbation.x` != "Unperturbed" & `perturbation.y` != "Unperturbed",
				     "Perturbed", "Unperturbed")) -> r 

gg <- ggplot(r, aes(x = log2FoldChange, y = logFC, 
		    text = `Gene name.x`,
		    color = `perturbation`)) +
	geom_point(alpha = .4) +
#	geom_density_2d(alpha = .5) +
	theme_classic()

pdf(paste0(outputfolder, "/", jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.pdf"))
	plot(gg)
dev.off()

htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
           paste0(normalizePath(outputfolder), 
           "/", jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.html"))

# for all clusters
lapply(cls, function(cl) {
		print(cl)
		sub_jobname <- paste0(jobname, "_cluster_", cl)
		deseq2 <- data.table::fread(paste0(outputfolder, "/", sub_jobname, "_deg_DESeq2_pbulk_results.tsv"))
		edger <- data.table::fread(paste0(outputfolder, "/", sub_jobname, "_deg_edgeR_pbulk_results.tsv"))

		r <- merge(deseq2, edger, by = 'gene_version')
		print(head(r))

		r %>%
			mutate(perturbation = ifelse(`perturbation.x` != "Unperturbed" & `perturbation.y` != "Unperturbed",
					     "Perturbed", "Unperturbed")) -> r 

		gg <- ggplot(r, aes(x = log2FoldChange, y = logFC, 
			    text = `Gene name.x`,
			    color = `perturbation`)) +
			geom_point(alpha = .4) +
		#	geom_density_2d(alpha = .5) +
			theme_classic()

		pdf(paste0(outputfolder, "/", sub_jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.pdf"))
			plot(gg)
		dev.off()

		htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
		           paste0(normalizePath(outputfolder), 
		           "/", sub_jobname, "_corr_edgeR_DESeq2_logFC_samplewise_pbulk.html"))
})


## Testing results relaibility 
# DXL1 gene
test_ds <- raw_eset["ENSG00000144355.15", ]
# CHI3L1
test_ds <- raw_eset["ENSG00000133048.13", ]
# DXL1 gene
test_ds <- raw_eset["ENSG00000144355", ]
# CHI3L1
test_ds <- raw_eset["ENSG00000133048", ]


co <- filter(meta, condition == "Control")[['barcode']]
pd <- filter(meta, condition == "IPD")[['barcode']]

mean(test_ds[co])
mean(test_ds[pd])

sum(test_ds[co])
sum(test_ds[pd])

length(test_ds[co])
length(test_ds[pd])



# Differential expresssion analysis modeling the zero-inflation
# zinbwave + DESeq2
#print(quantile(rowSums(raw_eset > 0)))
#keep <- rowSums(raw_eset > 0) >= 200
#table(keep)

#filter_raw_eset <- raw_eset[keep, ]

#se <- SummarizedExperiment(filter_raw_eset,
#			   colData = meta[colnames(filter_raw_eset), ])

#m <- zinbwave(se, X="~condition", BPPARAM=MulticoreParam(2))

#print(m)
# this looks non-scalable for this application 36K cells x 16K genes
# 2 cores overnight did not even create tje zinbawe object
#f It might be usefull for cluster-specific comparisons



# Clinical associations

##########################
# Aging
##########################

mca <- readRDS(paste0(inputfolder, '/PD_CO_hg_brain_seurat3_integrated_npcs_25_nvgs_4000_full_cca_seurat3_integrated_norm_seurat3_merged_clustered_resrange.rds'))  
ident_col <- "cell_ontology"
print(table(mca@meta.data[[ident_col]]))
Idents(mca) <- ident_col
meta <- mca@meta.data
samplevar <- "sample"
sample_meta <- read.delim("/home/prada/hbrain_project/data/sample_metadata_curated.tsv")
rownames(sample_meta) <- sample_meta[["sample"]]
mingenecount <- 50


## Counts
raw_eset <- mca@assays$RNA@counts 
print(raw_eset[1:5, 1:10])
print(dim(raw_eset))


# controls & IPD
sample_pattern <- "CO"
pbulk_edger_aging(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart, sample_pattern)

# For each cluster independently
cls <- unique(meta[[ident_col]])
cls 
mingenecount <- 50
cls[c(11:12)]

lapply(cls[c(11:12)], function(cl) {
	       
	       print(cl)

	       sub_meta <- dplyr::filter(meta, !!rlang::sym(ident_col) == cl)
	       print(dim(sub_meta))

	       cl_cells <- sub_meta[['barcode']]
	       rownames(sub_meta) <- cl_cells
	       print(head(sub_meta, 2))
	       print(dim(sub_meta))
	       sub_raw_eset <- raw_eset[, cl_cells]
	       print(dim(sub_raw_eset))
	       head(sub_raw_eset)

	       if(ncol(sub_raw_eset) < mingenecount) {

		       NULL

	       } else {

		       sub_jobname <- paste0(jobname, "aging_cluster_", cl)
		       pbulk_edger_aging(sub_raw_eset, sub_meta, sample_meta[unique(sub_meta[["sample"]]),], sub_jobname, mingenecount, refbiomart, sample_pattern)

	       }

})

# PD
sample_pattern <- "PD"
pbulk_edger_aging(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart, sample_pattern)


pbulk_edger_aging <- function(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart, sample_pattern) {
	# Sum counts per sample
	ss <- unique(meta[[samplevar]])

	ss <- ss[grep(sample_pattern, ss)]
   	print(ss)

	sample_meta <- sample_meta[grep(sample_pattern, rownames(sample_meta)), ] 
	print(ss)

	print(dim(raw_eset))

	sample_pbulk <- Reduce(cbind,
			lapply(ss, function(s) {
				cells <- rownames(meta)[which(meta[[samplevar]] == s)]
				print(head(cells))
				if(length(cells) == 1) {
					bulk <- data.frame(raw_eset[, cells])
				} else {
					bulk <- data.frame(rowSums(raw_eset[, cells]))
				}

				sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
				colnames(bulk) <- sname
				print(head(bulk, 2))
				bulk
			})
		)
	
	print("Genes x samples matrix")
	print(head(sample_pbulk, 2))
	print(head(sample_meta))

	print(rownames(sample_meta))
	print(colnames(sample_pbulk))

 	y <- DGEList(counts = sample_pbulk[, rownames(sample_meta)],
				      samples = sample_meta)

	# Checking sample library sizes
	discarded <- scater::isOutlier(y$samples$lib.size, log=TRUE, type="lower")

	if(any(discarded)) {
		print("Notice that these samples were removed from the DE analysis:")
		print(colnames(y)[discarded])
		y <- y[, !discarded]
	} else {
		print("Sample library sizes are more/less homogeneus")
		print(y$samples$lib.size)
	}

	print("Discarded?:")
	summary(discarded)

	keep <- filterByExpr(y, group=sample_meta[[groupvar]],
			     min.total.count = mingenecount)
	table(keep)

	y <- y[keep, ]
	
	y <- calcNormFactors(y)
	y$samples

	if(sample_pattern == "CO") {

		if(grepl("DaN", jobname)) {

			print(jobname)

			design <- model.matrix(~as.numeric(age_at_death),
				       y$samples)
			
		} else {

			# The actual design matrix /should it be customed for every dataset?/
			design <- model.matrix(~as.numeric(pmi) + as.factor(sex) + as.numeric(age_at_death),
					       y$samples)

		}

		jobname <- paste0(jobname, "_CO_individuals")

	} else if (sample_pattern == "PD"){

		design <- model.matrix(~as.factor(sex) + as.numeric(dis_dur) + as.numeric(age_at_death),
				       y$samples)

		jobname <- paste0(jobname, "_PD_patients")

	} else if (sample_pattern == "CO|PD") {

		design <- model.matrix(~as.factor(condition) + as.numeric(pmi) + as.factor(sex) + as.numeric(age_at_death),
				       y$samples)

		jobname <- paste0(jobname, "_PD_and_CO_individuals_aging")

	}

	print(`design`)
	y <- estimateDisp(y, design)
	print(summary(y$trended.dispersion))

	pdf(paste0(outputfolder, "/", jobname, "_biological_component_geneVariance_edgeR.pdf"))
		plotBCV(y)
	dev.off()
	
	fit <- glmQLFit(y, design, robust=TRUE)

	print(summary(fit$var.prior))

	pdf(paste0(outputfolder, "/", jobname, "_likelihood_dispersion_edgeR.pdf"))
		plotQLDisp(fit)
	dev.off()

	res <- glmQLFTest(fit, coef=ncol(design))
	summary(decideTests(res))
	topTags(res)


	coef_df <- data.frame(res$coefficient)
	coef_df[["gene"]] <- rownames(coef_df)
	print(head(coef_df))
	
	head(arrange(coef_df, -abs(as.numeric.age_at_death.)))

	de_res <- data.frame(topTags(res, nrow(res)))
	de_res[['gene']] <- rownames(de_res)

	# Check if stable gene ID was provided with or without version
	g2test <- de_res[['gene']][1]

	if(grepl("\\.\\d+$", g2test)) {
		de_res[['gene_version']] <- de_res[['gene']]
		de_res[['gene']] <- gsub("\\.\\d+$", "", de_res[['gene_version']])
	} else {
		de_res[['gene_version']] <- de_res[['gene']]
	}

	de_res <- merge(de_res, refbiomart, 
			by.x = "gene", by.y = "Gene stable ID") 

	print(head(de_res, 3))

	de_res %>%
		mutate(pertb = abs(logFC)*(-log10(PValue))) %>%
		arrange(-pertb) -> de_res

	de_res[['rank']] <- seq(nrow(de_res))
	print(head(de_res))
	
	metathr <- 0.05
	irank <- de_res[round(nrow(de_res)*metathr), ][['rank']]

	de_res %>%
		dplyr::mutate(perturbation = ifelse(`rank` <= irank, 
						    ifelse(logFC < 0, "Down", "Up"),
						    "Unperturbed")) -> de_res 

	write.table(de_res, paste0(outputfolder, "/", jobname, "_deg_edgeR_pbulk_results.tsv"),
		    sep = "\t", quote = FALSE, row.names = FALSE)

	gg <- ggplot(de_res, 
		     aes(x = logFC, y = -log10(PValue), 
			 text = `Gene name`,
			 color = perturbation)) +
		geom_point(alpha = 0.5) +
		scale_color_manual(values = c("blue", "grey", "red")) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_volcano.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_edgeR_volcanoplot.html"))

	vis_eset <- data.frame(res$fitted.values)
	vis_eset[['gene']] <- rownames(vis_eset)
	vis_eset <- reshape2::melt(vis_eset)
	vis_eset <- merge(vis_eset, sample_meta, 
			  by.x = "variable", by.y = "sample")

	vis_eset <- merge(vis_eset, 
			  select(de_res, c(gene_version, `Gene name`)), 
			  by.x = "gene", by.y = "gene_version")

	top50 <- head(arrange(de_res, rank), 50)[['gene_version']]
	print(head(top50))

	dplyr::filter(data.frame(vis_eset), 
		      gene %in% top50) %>%
		ggplot(aes(x = age_at_death, y = log10(value),
			   fill = age_at_death, color = age_at_death)) +
		geom_smooth(alpha = .2) +
		geom_point(alpha = .5) +
		theme_classic() +
		facet_wrap(~`Gene.name`, nrow = 10) -> ggex

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_top_50gene_violin.pdf"),
		width = 20, height = 20)
		plot(ggex)
	dev.off()
}

## monocle3 to validate and trajectory edgeR pseudobulk to identify aging signature

# Counts
raw_eset <- mca@assays$RNA@counts 
print(raw_eset[1:5, 1:10])
print(dim(raw_eset))

row_df <- mca@assays$RNA@meta.features
row_df[["gene"]] <- rownames(row_df)
dim(row_df)
row_df <- merge(row_df, refbiomart, by.x ="gene", by.y = "Gene stable ID", all.x = TRUE) 
row_df <- filter(row_df, !duplicated(`gene`))
dim(row_df)
colnames(row_df)[7] <- "gene_short_name"
rownames(row_df) <- row_df[["gene"]]
head(row_df)

## Monocle diff. exp. framework
cds <- new_cell_data_set(expression_data = mca@assays$RNA@counts,
			 cell_metadata = mca@meta.data,
			 gene_metadata = row_df[rownames(mca@assays$RNA@counts), ])
dim(cds)

cell_types <- unique(pData(cds)[['cell_ontology']])
cell_type <- "Astrocytes"
cell_type <- "DaNs"

cell_types
sample_pattern <- "PD"
sample_pattern
options(future.globals.maxSize=61943040000000)

for(cell_type in cell_types[c(1:9)]) {

	       sub_cells <- cds[, pData(cds)$cell_ontology == cell_type]

		if(sample_pattern == "CO") {

			sub_cells <- sub_cells[, pData(sub_cells)$condition == "Control"]

			reg_feat <- read.delim(paste0(inputfolder, "/PD_CO_hg_brainaging_cluster_", cell_type, "_CO_individuals_deg_edgeR_pbulk_results.tsv"))

		} else if(sample_pattern == "PD") {

			sub_cells <- sub_cells[, pData(sub_cells)$condition == "IPD"]
			
			reg_feat <- read.delim(paste0(inputfolder, "/PD_CO_hg_brainaging_cluster_", cell_type, "_PD_patients_deg_edgeR_pbulk_results.tsv"))

		}

	       print(dim(sub_cells))

#	       reg_feat <- read.delim(paste0(inputfolder, "/PD_CO_hg_brain_aging_PDpatientsaging_cluster_", cell_type, "_CO_individuals_deg_edgeR_pbulk_results.tsv"))
#	       reg_feat <- read.delim(paste0(inputfolder, "/PD_CO_hg_brainaging_cluster_", cell_type, "_PD_and_CO_individuals_aging_deg_edgeR_pbulk_results.tsv"))

	       print(head(reg_feat))

	       print(dim(filter(reg_feat, PValue < 0.05)))

	       aging_feat <- as.character(filter(reg_feat, PValue < 0.05)[['gene']])
	       aging_feat <- aging_feat[!aging_feat%in%c("ENSG00000205696", "ENSG00000233359")]
	       #fit_coefs <- coefficient_table(filter(cgene_fits, !gene%in%c("ENSG00000205696", "ENSG00000233359")))

	       head(aging_feat)

	       sub_cells <- sub_cells[aging_feat, ]

	       print(dim(sub_cells))
	       
	       print(quantile(rowSums(exprs(sub_cells))))

		sub_cells <- sub_cells[rowSums(exprs(sub_cells)) > mingenecount, ]

	       print(dim(sub_cells))
	       print(colnames(colData(sub_cells)))

		sub_cells <- preprocess_cds(sub_cells,
					    method = "PCA",
					    num_dim = 25,
					    norm_method = "log", 
					    use_genes = NULL,
					    scaling = TRUE,  
					    verbose = TRUE, 
					    cores = 10)

		print(cell_type)

		if(sample_pattern == "PD") {

			cgene_fits <- fit_models(sub_cells, 
						 model_formula_str = "~pmi+age_at_death+dis_dur+sex.y",
						 expression_family = "negbinomial",
						 cores = 5)

		} else if(sample_pattern == "CO") {

			cgene_fits <- fit_models(sub_cells, 
						 model_formula_str = "~pmi+age_at_death+sex.y",
						 expression_family = "negbinomial",
						 cores = 5)

		}

	       fit_coefs <- coefficient_table(cgene_fits)
	       colnames(fit_coefs)
	       table(fit_coefs$term)
	       head(fit_coefs)

	       coef_res <- fit_coefs %>% 
		       filter(term == "age_at_death") %>% 
		       filter (q_value < 0.05) %>%
		       select(gene, gene_short_name, term, q_value, estimate) %>% 
		       data.frame

	       print(head(coef_res))
	       print(dim(coef_res))

	       coef_res <- mutate(coef_res, `r` = -log10(q_value)*abs(estimate)) %>%
		       arrange(-r)


	       write.table(arrange(coef_res, q_value), paste0(outputfolder, "/monocle3_reg_aging_signature__", cell_type, "_", sample_pattern, ".tsv"),
			   sep = "\t", quote = FALSE, row.names = FALSE)

	      
	       pdf(paste0(outputfolder, "/monocle3_reg_aging_sign_", cell_type, "_", sample_pattern, ".pdf"), width = 20, height = 20)

	       plot(
	       		plot_genes_violin(sub_cells[head(arrange(coef_res, -estimate),25)[["gene"]], ], 
				 group_cells_by="age_at_death", ncol=5) 

			)

	       dev.off()

	       if(FALSE) {


		       sub_cells <- sub_cells[coef_res[["gene"]], ]
		       sub_cells <- align_cds(sub_cells, alignment_group = "sample")
		       sub_cells <- reduce_dimension(sub_cells, reduction_method = "Aligned")
		       sub_cells <- reduce_dimension(sub_cells, reduction_method = "UMAP", preprocess_method = "Aligned")

		       sub_cells <- cluster_cells(sub_cells, preprocess_method = "Aligned")
		       sub_cells <- learn_graph(sub_cells,
						use_partition = FALSE)

		       sub_cells <- order_cells(sub_cells, root_pr_nodes=get_earliest_principal_node(sub_cells))
	
		       pdf(paste0(outputfolder, '/', jobname, "_", cell_type, "_CO_aging_signature_trajectory.pdf"))
		       plot_cells(sub_cells,
				  color_cells_by = "sample",
				  label_groups_by_cluster=FALSE,
				  label_leaves=TRUE,
				  label_branch_points=TRUE)
		       plot_cells(sub_cells,
				  color_cells_by = "batch",
				  label_groups_by_cluster=FALSE,
				  label_leaves=FALSE,
				  label_branch_points=FALSE)
		       plot_cells(sub_cells,
				  color_cells_by = "age_at_death",
				  label_groups_by_cluster=FALSE,
				  label_leaves=FALSE,
				  label_branch_points=FALSE)
		       plot_cells(sub_cells,
				  color_cells_by = "pseudotime",
				  label_groups_by_cluster=FALSE,
				  label_leaves=FALSE,
				  label_branch_points=FALSE)
		       dev.off()

	 	       pdf(paste0(outputfolder, '/', jobname, "_", cell_type, "_CO_aging_signature_trajectory_genes_pseudotime.pdf"), height = 35)
		       plot_genes_in_pseudotime(sub_cells[rowData(sub_cells)$gene %in% head(arrange(coef_res, q_value),20)[["gene"]], ],
	                         color_cells_by="dis_dur",
	                         min_expr=0.5)
		       dev.off()
		}
}

#
# Visualizing aging signature
#

# PD

files <- list.files(path = outputfolder, pattern = "monocle3_reg_aging_signature__")
files <- files[-grep("COPD", files)]
files

# CO
co_files <- files[grep("CO.tsv", files)]
names(co_files) <- gsub("monocle3_reg_aging_signature__|\\.tsv", "", co_files)
co_files

aging <- lapply(co_files, function(f) read.delim(paste0(outputfolder, '/', f)))

aging_tab <- do.call(rbind,
		     lapply(names(aging), function(na) {
				    ag <- aging[[na]]
				    ag[["celltype"]] <- na
				    ag
				 })
		     )

aging_tab <- mutate(aging_tab, 
		    `sign` = ifelse(estimate < 0, "neg", "pos"))

head(aging_tab)
dim(aging_tab)

aging_tab <- filter(aging_tab, q_value < 0.05 & abs(estimate) > 0.05)

gg <- ggplot(aging_tab, aes(x = `sign`)) +
	geom_bar(aes(fill = `sign`)) +
	theme_classic() +
	scale_fill_manual(values = c("darkgreen", "darkblue")) +
	facet_wrap(~celltype, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_aging_signature_CO.pdf"))
	plot(gg)
dev.off()

gg <- ggplot(aging_tab, aes(x = `sign`, y = estimate)) +
	geom_violin(aes(fill = `sign`)) +
	theme_classic() +
	scale_fill_manual(values = c("darkgreen", "darkblue")) +
	facet_wrap(~celltype, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_aging_signature_estimate_CO.pdf"))
	plot(gg)
dev.off()

dim(aging_tab)
head(aging_tab)
write.table(aging_tab, paste0(outputfolder, '/', jobname, '_aging_signature_tab_CO.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE)

aging_mkrs <- split(aging_tab, aging_tab[['celltype']])
aging_mkrs_co <- lapply(aging_mkrs, function(a) as.character(a[['gene_short_name']]))
aging_mkrs_co

# PD
pd_files <- files[grep("_PD.tsv", files)]
names(pd_files) <- gsub("monocle3_reg_aging_signature__|\\.tsv", "", pd_files)
pd_files

aging <- lapply(pd_files, function(f) read.delim(paste0(outputfolder, '/', f)))

aging_tab <- do.call(rbind,
		     lapply(names(aging), function(na) {
				    ag <- aging[[na]]
				    ag[["celltype"]] <- na
				    ag
				 })
		     )

aging_tab <- mutate(aging_tab, 
		    `sign` = ifelse(estimate < 0, "neg", "pos"))

head(aging_tab)
dim(aging_tab)

aging_tab <- filter(aging_tab, q_value < 0.05 & abs(estimate) > 0.05)

gg <- ggplot(aging_tab, aes(x = `sign`)) +
	geom_bar(aes(fill = `sign`)) +
	theme_classic() +
	scale_fill_manual(values = c("darkgreen", "darkblue")) +
	facet_wrap(~celltype, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_aging_signature_PD.pdf"))
	plot(gg)
dev.off()

gg <- ggplot(aging_tab, aes(x = `sign`, y = estimate)) +
	geom_violin(aes(fill = `sign`)) +
	theme_classic() +
	scale_fill_manual(values = c("darkgreen", "darkblue")) +
	facet_wrap(~celltype, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_aging_signature_estimate_PD.pdf"))
	plot(gg)
dev.off()

dim(aging_tab)
head(aging_tab)
write.table(aging_tab, paste0(outputfolder, '/', jobname, '_aging_signature_tab_PD.tsv'),
	    sep = "\t", quote = FALSE, row.names = FALSE)

aging_mkrs <- split(aging_tab, aging_tab[['celltype']])
aging_mkrs_pd <- lapply(aging_mkrs, function(a) as.character(a[['gene_short_name']]))
aging_mkrs_pd


aging_mkrs <- append(aging_mkrs_co, aging_mkrs_pd)


calc_simil <- function(mods, type) { # take a list of string vectors and compare them against them. It returns the -log10(Fisher.pval)
  gene.universe <- as.vector(unlist(mods))
  result <- lapply(names(mods)[-length(names(mods))], function(mod) {
    i <- which(names(mods) == mod) + 1
    similarity <- c()
    while (i <= length(mods)) {
      if(type == "fisher") {
        s <- -log10(fisher_exact_test(mods[[mod]], mods[[i]], gene.universe))
      } else if (type == "intersect") {
        s <- length(intersect(mods[[mod]], mods[[i]]))
      } else if (type == "jaccard") {
        s <- set_similarity(mods[[mod]], mods[[i]], method = "Jaccard")
      }
      similarity <- append(similarity,s)
      i <- i+1
    }
    result = cbind("DEG1" = rep(mod, length(similarity) ), 
                   "DEG2" = names(mods)[(which(names(mods)==mod)+1):length(names(mods))], 
                   "Weight" = as.numeric(similarity))
    return(result)
  })
  return(result)
}

fisher_exact_test <- function(m1, m2, universe) { # Function: Calculate the Fisher Exact Test (Diogenes)
  isect <- length(intersect(m1, m2))
  m1uniq <- length(m1[!m1 %in% isect])
  m2uniq <- length(m2[!m2 %in% isect])
  tot <- length(universe) - length(unique(c(m1,m2)))
  mm <- c(isect, m1uniq, m2uniq, tot)
  to_fisher <- matrix(mm, nrow = 2)
  if (all(to_fisher > 0)) { 
    k <- fisher.test(to_fisher, alternative = 'greater')$p.value
    return(k)
  } else { # If any cell os the matrix is 0
    return(1)
  }
}

write_simil_matrix <- function(simil.3col) {
  simil.3col[which(simil.3col[,3]=="Inf"),3] <- max(as.numeric(simil.3col[-which(simil.3col[,3]=='Inf'),3]))
  vars <- unique(c(simil.3col[,1], simil.3col[,2]))
  nvars <- length(vars)
  simil.mat <- data.frame(matrix(0, nrow=nvars, ncol=nvars))
  colnames(simil.mat) <- vars
  row.names(simil.mat) <- vars
  for(s in 1:nrow(simil.3col)) {
    simil.mat[simil.3col[s,1], simil.3col[s,2]] <- simil.3col[s,3]
    simil.mat[simil.3col[s,2], simil.3col[s,1]] <- simil.3col[s,3]
  }
  return(simil.mat)
}

smilmat <- calc_simil(aging_mkrs, "intersect")
smilmat <- Reduce(rbind, smilmat)
smilmat <- write_simil_matrix(smilmat)
smilmat <- smilmat[names(aging_mkrs_co), names(aging_mkrs_pd)]

smilmat_f <- calc_simil(aging_mkrs, "fisher")
smilmat_f <- Reduce(rbind, smilmat_f)
smilmat_f <- write_simil_matrix(smilmat_f)
smilmat_f <- smilmat_f[names(aging_mkrs_co), names(aging_mkrs_pd)]


colfunc <- colorRampPalette(c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f"))(100)
pdf(paste0(outputfolder, '/', jobname, "_aging_signature_intersection_CO_vs_PD.pdf"),
    width = 5, height = 5)
  gplots::heatmap.2(data.matrix(smilmat_f), cellnote=ifelse(data.matrix(smilmat)==0, NA, data.matrix(smilmat)), Colv = FALSE, Rowv = FALSE, scale = "none", col = colfunc, trace = "none", notecol="black")
dev.off()


aging_mkrs <- split(filter(aging_tab, `sign` == 'neg'), filter(aging_tab, `sign` == 'neg')[['celltype']])
aging_mkrs <- lapply(aging_mkrs, function(a) as.character(a[['gene_short_name']]))
aging_mkrs

smilmat <- calc_simil(aging_mkrs, "intersect")
smilmat <- Reduce(rbind, smilmat)
smilmat <- write_simil_matrix(smilmat)
smilmat

colfunc <- colorRampPalette(c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f"))
pdf(paste0(outputfolder, '/', jobname, "_aging_signature_intersection_down.pdf"),
    width = 5, height = 5)
  gplots::heatmap.2(data.matrix(smilmat), cellnote=ifelse(data.matrix(smilmat)==0, NA, data.matrix(smilmat)), scale = "none", col = colfunc(100), trace = "none", notecol="black")
dev.off()

aging_mkrs <- split(filter(aging_tab, `sign` == 'pos'), filter(aging_tab, `sign` == 'pos')[['celltype']])
aging_mkrs <- lapply(aging_mkrs, function(a) as.character(a[['gene_short_name']]))
aging_mkrs

smilmat <- calc_simil(aging_mkrs, "intersect")
smilmat <- Reduce(rbind, smilmat)
smilmat <- write_simil_matrix(smilmat)
smilmat

colfunc <- colorRampPalette(c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f"))
pdf(paste0(outputfolder, '/', jobname, "_aging_signature_intersection_pos.pdf"),
    width = 5, height = 5)
  gplots::heatmap.2(data.matrix(smilmat), cellnote=ifelse(data.matrix(smilmat)==0, NA, data.matrix(smilmat)), scale = "none", col = colfunc(100), trace = "none", notecol="black")
dev.off()


# -----------------------------------------
# Disease duration
# -----------------------------------------


# PD
pbulk_edger_disdur <- function(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart, sample_pattern) {
	# Sum counts per sample
	ss <- unique(meta[[samplevar]])
	ss <- ss[grep(sample_pattern, ss)]
   	print(ss)

	sample_meta <- sample_meta[grep(sample_pattern, rownames(sample_meta)), ] 

	print(dim(raw_eset))

	sample_pbulk <- Reduce(cbind,
			lapply(ss, function(s) {
				cells <- rownames(meta)[which(meta[[samplevar]] == s)]
				print(head(cells))
				if(length(cells) == 1) {
					bulk <- data.frame(raw_eset[, cells])
				} else {
					bulk <- data.frame(rowSums(raw_eset[, cells]))
				}

				sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
				colnames(bulk) <- sname
				print(head(bulk, 2))
				bulk
			})
		)

	print(dim(sample_pbulk))
	sample_pbulk <- sample_pbulk[which(rowSums(sample_pbulk) > mingenecount), ]
	print(dim(sample_pbulk))

 	y <- DGEList(counts = sample_pbulk[, rownames(sample_meta)],
		      samples = sample_meta)

	# Checking sample library sizes
	discarded <- scater::isOutlier(y$samples$lib.size, log=TRUE, type="lower")

	if(any(discarded)) {
		print("Notice that these samples were removed from the DE analysis:")
		print(colnames(y)[discarded])
		y <- y[, !discarded]
	} else {
		print("Sample library sizes are more/less homogeneus")
		print(y$samples$lib.size)
	}

	print("Discarded?:")
	print(summary(discarded))

	y <- calcNormFactors(y)
	print(y$samples[,1:2])

	if (sample_pattern == "PD") {

		if(grepl("CADPS2|GABA", jobname)) {

			print(jobname)

			design <- model.matrix(~as.numeric(age_at_death) + as.numeric(dis_dur),
					       y$samples)
			print(`design`)
			y <- estimateDisp(y)

			} else {

				design <- model.matrix(~as.factor(sex) + as.numeric(age_at_death) + as.numeric(dis_dur),
						       y$samples)
				print(`design`)
				y <- estimateDisp(y, `design`)

			}

		jobname <- paste0(jobname, "_PD_patients_disease_duration")
	}

	print(summary(y$trended.dispersion))

	pdf(paste0(outputfolder, "/", jobname, "_biological_component_geneVariance_edgeR.pdf"))
		plotBCV(y)
	dev.off()

	print("Fitting:")
	
	fit <- glmQLFit(y, design, robust=TRUE)

	coef_fit <- data.frame(coef(fit))
	coef_fit[['gene']] <- rownames(coef_fit)
	coef_fit <- merge(coef_fit, refbiomart,
			by.x = "gene", by.y = "Gene stable ID")
	print(head(coef_fit))
	write.table(coef_fit, paste0(outputfolder, "/", jobname, "_deg_edgeR_pbulk_coeficients.tsv"),
		    sep = "\t", quote = FALSE, row.names = FALSE)

	print(summary(fit$var.prior))

	pdf(paste0(outputfolder, "/", jobname, "_likelihood_dispersion_edgeR.pdf"))
		plotQLDisp(fit)
	dev.off()

	res <- glmQLFTest(fit, coef=ncol(design))
	print(summary(decideTests(res)))
	print(topTags(res))

	
	de_res <- data.frame(topTags(res, nrow(res)))
	de_res[['gene']] <- rownames(de_res)

	# Check if stable gene ID was provided with or without version
	g2test <- de_res[['gene']][1]

	if(grepl("\\.\\d+$", g2test)) {
		de_res[['gene_version']] <- de_res[['gene']]
		de_res[['gene']] <- gsub("\\.\\d+$", "", de_res[['gene_version']])
	} else {
		de_res[['gene_version']] <- de_res[['gene']]
	}

	de_res <- merge(de_res, refbiomart, 
			by.x = "gene", by.y = "Gene stable ID") 

	print(head(de_res, 3))

	de_res %>%
		mutate(pertb = abs(logFC)*(-log10(PValue))) %>%
		arrange(-pertb) -> de_res

	de_res[['rank']] <- seq(nrow(de_res))
	print(head(de_res))
	
	metathr <- 0.05
	irank <- de_res[round(nrow(de_res)*metathr), ][['rank']]

	de_res %>%
		dplyr::mutate(perturbation = ifelse(`rank` <= irank, 
						    ifelse(logFC < 0, "Down", "Up"),
						    "Unperturbed")) -> de_res 

	write.table(de_res, paste0(outputfolder, "/", jobname, "_deg_edgeR_pbulk_results.tsv"),
		    sep = "\t", quote = FALSE, row.names = FALSE)

	gg <- ggplot(de_res, 
		     aes(x = logFC, y = -log10(PValue), 
			 text = `Gene name`,
			 color = perturbation)) +
		geom_point(alpha = 0.5) +
		scale_color_manual(values = c("blue", "grey", "red")) +
		theme_classic()

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_volcano.pdf"))
		plot(gg)
	dev.off()

	htmlwidgets::saveWidget(as_widget(ggplotly(gg)), 
            paste0(normalizePath(outputfolder), 
	           "/", jobname, "_edgeR_volcanoplot.html"))

	vis_eset <- data.frame(res$fitted.values)
	vis_eset[['gene']] <- rownames(vis_eset)
	vis_eset <- reshape2::melt(vis_eset)
	vis_eset <- merge(vis_eset, sample_meta, 
			  by.x = "variable", by.y = "sample")

	vis_eset <- merge(vis_eset, 
			  select(de_res, c(gene_version, `Gene name`)), 
			  by.x = "gene", by.y = "gene_version")

	top50 <- head(arrange(de_res, rank), 50)[['gene_version']]
	print(head(top50))

	dplyr::filter(data.frame(vis_eset), 
		      gene %in% top50) %>%
		ggplot(aes(x = dis_dur, y = log10(value),
			   fill = age_at_death, color = sex)) +
		geom_smooth(alpha = .2) +
		geom_point(alpha = .5) +
		theme_classic() +
		facet_wrap(~`Gene.name`, nrow = 10) -> ggex

	pdf(paste0(outputfolder, "/", jobname, "_edgeR_samplewise_pbulk_top_50gene_violin.pdf"),
		width = 20, height = 20)
		plot(ggex)
	dev.off()
}
## edgeR pseudobulk to identify aging signature
## monocle3 to validate and trajectory
cls <- unique(meta[[ident_col]])
cls 
mingenecount <- 100
sample_pattern <- "PD"
pbulk_edger_disdur(raw_eset, meta, sample_meta, jobname, mingenecount, refbiomart, sample_pattern)

# For each cluster independently

lapply(cls[c(1:11)], function(cl) {
	       print(cl)

	       sub_meta <- dplyr::filter(meta, !!rlang::sym(ident_col) == cl)

	       cl_cells <- sub_meta[['barcode']]
	       rownames(sub_meta) <- cl_cells
	       sub_raw_eset <- raw_eset[, cl_cells]

	       if(ncol(sub_raw_eset) < mingenecount) {
		       NULL
	       } else {
		       sub_jobname <- paste0(jobname, "_", cl)

		       pbulk_edger_disdur(sub_raw_eset, sub_meta, sample_meta[unique(sub_meta[["sample"]]),], sub_jobname, mingenecount, refbiomart, sample_pattern)
	       }
})

# Double filtering with monocle3
# Counts
raw_eset <- mca@assays$RNA@counts 
print(raw_eset[1:5, 1:10])
print(dim(raw_eset))

row_df <- mca@assays$RNA@meta.features
row_df[["gene"]] <- rownames(row_df)
dim(row_df)
row_df <- merge(row_df, refbiomart, by.x ="gene", by.y = "Gene stable ID", all.x = TRUE) 
row_df <- filter(row_df, !duplicated(`gene`))
dim(row_df)
colnames(row_df)[7] <- "gene_short_name"
rownames(row_df) <- row_df[["gene"]]
head(row_df)


## Monocle diff. exp. framework
cds <- new_cell_data_set(expression_data = mca@assays$RNA@counts,
			 cell_metadata = mca@meta.data,
			 gene_metadata = row_df[rownames(mca@assays$RNA@counts), ])
dim(cds)

cell_types <- unique(pData(cds)[['cell_ontology']])
cell_type <- "Astrocytes"
cell_type <- "GABA"

cell_types
options(future.globals.maxSize=61943040000000)

for(cell_type in cell_types[c(10:11)]) {

	       sub_cells <- cds[, pData(cds)$cell_ontology == cell_type]
	       sub_cells <- sub_cells[, pData(sub_cells)$condition == "IPD"]
	       print(dim(sub_cells))

	       reg_feat <- read.delim(paste0(inputfolder, "/PD_CO_hg_brain_PD_patients_disease_duration_", cell_type, "_PD_patients_disease_duration_deg_edgeR_pbulk_results.tsv"))
	       print(head(reg_feat))

	       print(dim(filter(reg_feat, PValue < 0.05)))

	       aging_feat <- as.character(filter(reg_feat, PValue < 0.05)[['gene']])
	       head(aging_feat)

	       sub_cells <- sub_cells[aging_feat, ]

	       print(dim(sub_cells))
	       
	       print(quantile(rowSums(exprs(sub_cells))))

		sub_cells <- sub_cells[rowSums(exprs(sub_cells)) > mingenecount, ]

	       print(dim(sub_cells))
	       print(colnames(colData(sub_cells)))

		sub_cells <- preprocess_cds(sub_cells,
					    method = "PCA",
					    num_dim = 25,
					    norm_method = "log", 
					    use_genes = NULL,
					    scaling = TRUE,  
					    verbose = TRUE, 
					    cores = 10)

		print(cell_type)


	       cgene_fits <- fit_models(sub_cells, 
			 model_formula_str = "~pmi+age_at_death+sex.y+dis_dur",
			 expression_family = "negbinomial",
			 cores = 5)

#	       fit_coefs <- coefficient_table(filter(cgene_fits, !gene%in%c("ENSG00000205696", "ENSG00000233359")))
	       
	       fit_coefs <- coefficient_table(cgene_fits)
	       colnames(fit_coefs)
	       table(fit_coefs$term)

	       coef_res <- fit_coefs %>% 
		       filter(term == "dis_dur") %>% 
		       filter (q_value < 0.05) %>%
		       select(gene, gene_short_name, term, q_value, estimate) %>% 
		       data.frame

	       if(nrow(coef_res) == 0) {
		       next
	       } else {

	       print(head(coef_res))
	       print(dim(coef_res))

#	       coef_res <- merge(coef_res, refbiomart, by.x = "gene", by.y = "Gene stable ID")

	       write.table(coef_res, paste0(outputfolder, "/monocle3_reg_disease_duration_", cell_type, ".tsv"),
			   sep = "\t", quote = FALSE, row.names = FALSE)

	       pdf(paste0(outputfolder, "/monocle3_reg_disease_duration_IPD_", cell_type, ".pdf"), width = 20, height = 20)

	       plot(
	       		plot_genes_violin(sub_cells[head(arrange(coef_res, q_value),49)[["gene"]], ], 
				 group_cells_by="dis_dur", ncol=7) 

			)

	       dev.off()

	       sub_cells <- sub_cells[coef_res[["gene"]], ]
#	       sub_cells <- align_cds(sub_cells, alignment_group = "sample")
#	       sub_cells <- reduce_dimension(sub_cells, reduction_method = "Aligned")
	       sub_cells <- reduce_dimension(sub_cells, 
					     reduction_method = "UMAP", 
					     preprocess_method = "PCA")

	       sub_cells <- cluster_cells(sub_cells, preprocess_method = "PCA")
	       sub_cells <- learn_graph(sub_cells,
					use_partition = FALSE)

	       sub_cells <- order_cells(sub_cells, root_pr_nodes=get_earliest_principal_node(sub_cells))

	       pdf(paste0(outputfolder, '/', jobname, "_", cell_type, "_IPD_disease_duration_trajectory.pdf"))
	       plot(
	       plot_cells(sub_cells,
			  color_cells_by = "sample",
			  label_groups_by_cluster=FALSE,
			  label_leaves=TRUE,
			  label_branch_points=TRUE)
	       )
	       plot(
	       plot_cells(sub_cells,
			  color_cells_by = "batch",
			  label_groups_by_cluster=FALSE,
			  label_leaves=FALSE,
			  label_branch_points=FALSE)
	       )
	       plot(
	       plot_cells(sub_cells,
			  color_cells_by = "age_at_death",
			  label_groups_by_cluster=FALSE,
			  label_leaves=FALSE,
			  label_branch_points=FALSE)
	       )
	       plot(
	       plot_cells(sub_cells,
			  color_cells_by = "dis_dur",
			  label_groups_by_cluster=FALSE,
			  label_leaves=FALSE,
			  label_branch_points=FALSE)
	       )
	       plot(
	       plot_cells(sub_cells,
			  color_cells_by = "pseudotime",
			  label_groups_by_cluster=FALSE,
			  label_leaves=FALSE,
			  label_branch_points=FALSE)
	       )
	       dev.off()

 	       pdf(paste0(outputfolder, '/', jobname, "_", cell_type, "_PD_disease_duration_trajectory_genes_pseudotime.pdf"), height = 35)
	       plot(
	       plot_genes_in_pseudotime(sub_cells[rowData(sub_cells)$gene %in% head(arrange(coef_res, q_value),20)[["gene"]], ],
                        color_cells_by="dis_dur",
                       min_expr=0.5)
	       )
	       dev.off()
	       
	       }

}

get_earliest_principal_node <- function(cds){
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex))))]
  root_pr_nodes

}

#
# visualize disease duration
#


files <- list.files(path = outputfolder, pattern = "monocle3_reg_disease_duration_")
files <- files[-grep("IPD|GABA|CAD", files)]
files

names(files) <- gsub("monocle3_reg_disease_duration_|\\.tsv", "", files)
files

aging <- lapply(files, function(f) read.delim(paste0(outputfolder, '/', f)))

aging_tab <- do.call(rbind,
		     lapply(names(aging), function(na) {
				    ag <- aging[[na]]
				    ag[["celltype"]] <- na
				    ag
				 })
		     )

aging_tab <- mutate(aging_tab, 
		    `sign` = ifelse(estimate < 0, "neg", "pos"))

head(aging_tab)
dim(aging_tab)

aging_tab <- filter(aging_tab, q_value < 0.05 & abs(estimate) > 0.05)
dim(aging_tab)

gg <- ggplot(aging_tab, aes(x = `sign`)) +
	geom_bar(aes(fill = `sign`)) +
	theme_classic() +
	scale_fill_manual(values = c("#bf812d", "#35978f")) +
	facet_wrap(~celltype, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_disease_duration_signature_PD.pdf"))
	plot(gg)
dev.off()

gg <- ggplot(aging_tab, aes(x = `sign`, y = estimate)) +
	geom_violin(aes(fill = `sign`)) +
	theme_classic() +
	scale_fill_manual(values = c("#bf812d", "#35978f")) +
	facet_wrap(~celltype, nrow = 1)

pdf(paste0(outputfolder, '/', jobname, "_disease_duration_signature_PD_estimate.pdf"))
	plot(gg)
dev.off()

write.table(aging_tab, paste0(outputfolder, "/monocle3_reg_disease_duration_cell_types_midbrain.tsv"),
	    sep = "\t", quote = FALSE, row.names = FALSE)



aging_mkrs <- split(aging_tab, aging_tab[["celltype"]])
aging_mkrs <- lapply(aging_mkrs, function(x) as.character(x[["gene_short_name"]]))
aging_mkrs


smilmat <- calc_simil(aging_mkrs, "intersect")
smilmat <- Reduce(rbind, smilmat)
smilmat <- write_simil_matrix(smilmat)

smilmat_f <- calc_simil(aging_mkrs, "fisher")
smilmat_f <- Reduce(rbind, smilmat_f)
smilmat_f <- write_simil_matrix(smilmat_f)

colfunc <- colorRampPalette(c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f"))(100)
pdf(paste0(outputfolder, '/', jobname, "_disease_duration_signature_intersection.pdf"),
    width = 5, height = 5)
  gplots::heatmap.2(data.matrix(smilmat_f), cellnote=ifelse(data.matrix(smilmat)==0, NA, data.matrix(smilmat)), symm=TRUE, scale = "none", col = colfunc, trace = "none", notecol="black")
dev.off()



# ----------------------------------------------------------------
# Cluster cell types by peusobulk correlation distances
# ----------------------------------------------------------------


samplevar <- "cell_ontology"
ss <- unique(meta[[samplevar]])
print(ss)


print(dim(raw_eset))

sample_pbulk <- Reduce(cbind,
		lapply(ss, function(s) {
			cells <- rownames(meta)[which(meta[[samplevar]] == s)]
			print(head(cells))
			if(length(cells) == 1) {
				bulk <- data.frame(raw_eset[, cells])
			} else {
				bulk <- data.frame(rowSums(raw_eset[, cells]))
			}
			sname <- gsub("_filtered_FILTERED_qc_seurat3", "", s)
			colnames(bulk) <- sname
			print(head(bulk, 2))
			bulk
		})
	)

dim(sample_pbulk)

head(sample_pbulk)

cell_exp <- sample_pbulk[mca@assays$SCT@var.features, ]
head(cell_exp)


mm <- sum(cell_exp)

norm_cell_exp <- apply(cell_exp, 2, function(x) log2((x/(mm/1000000))+1))
head(norm_cell_exp)


corr <- cor(norm_cell_exp)
dim(corr)

head(corr)
quantile(corr)
palette = colorRampPalette(c("grey90", "darkred"))(20)
pdf("~/PD_midbrain_tree.pdf")
gplots::heatmap.2(x = corr, col = palette, symm = TRUE)
dev.off()
