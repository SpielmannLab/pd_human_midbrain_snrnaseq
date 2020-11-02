#
# Quality Control of cell_data_set objects (Monocle3)
#

"Quality Control of cell_data_set objects (Monocle3) 

Usage: qc_sc.R --scobject=<file> --specie=<value> --infolder=<folder> --outfolder=<folder> --dataset=<value>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --scobject=<file>    *.rds file SingleCellObject.
  --specie=<value>     Either 'hg' or 'mm'.
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
          'cowplot', "dplyr", "ggwordcloud", "monocle3", "dplyr")

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

# --------------
# Main function
# --------------

plot_qc <- function(cds, mt_genes, rb_genes, ds, jobname) {

    # ---------------------------
    # Feature centered histograms
    # ---------------------------

    logcount <- log10(rowSums(exprs(cds)))
    features <- data.frame('feat_count_depth' = rowSums(exprs(cds)),
                           'feat_log10_count_depth' = ifelse(is.infinite(logcount),
						       0.001, logcount),
			   'feat_cell_depth' = rowSums(exprs(cds) != 0),
			   'gene_short_name' = fData(cds)[rownames(exprs(cds)), ]$gene_short_name
			   )
    features[['id']] <- rownames(features)

    librarysize <- sum(features$feat_count_depth)
    features$pctlib <- features$feat_count_depth/librarysize
    
    libsizecell <- colSums(exprs(cds))

    top20feat <- head(arrange(features, -feat_count_depth),20)$id
    top20exp <- exprs(cds)[top20feat, ]
    pctlibsizecell <- lapply(seq(ncol(top20exp)), function(i) top20exp[, i]/libsizecell[i])
    names(pctlibsizecell) <- colnames(top20exp)
    top20 <- do.call(rbind, 
	      lapply(names(pctlibsizecell), function(cell) {
		  data.frame('pctlibcellsize'=pctlibsizecell[[cell]],
			     'id'=names(pctlibsizecell[[cell]]),
			     'cell'=cell)
		  })
	     )

    print(head(top20))
    top20 <- merge(top20, select(features, c(`id`, gene_short_name)), by='id')
    print("Done")
    
    gtop20 <- ggplot(top20, aes(x = gene_short_name, y=pctlibcellsize)) +
	    geom_jitter(aes(group = gene_short_name), alpha=.2, color = "grey", shape=".") +
	    geom_violin(alpha=.3, fill = "grey") +
	    scale_x_discrete(limits=rev(head(arrange(features, -feat_count_depth),20)$gene_short_name)) +
	    coord_flip() +
	    theme_classic() +
	    theme(axis.text=element_text(size=5)) +
	    ggtitle(paste(ds, "Top 20 deepest genes"),
		    subtitle = paste0("They account for ~",
				     ceiling(sum(features[top20feat,]$pctlib)*100),
				     "% of total-counts"))

     g1 <- ggplot(`features`, aes(`feat_count_depth`)) +
	geom_histogram(color = "grey", fill = "grey", 
		       bins = ceiling(max(features$feat_count_depth)/100)) +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(features$feat_count_depth))), 
		 x = mean(features$feat_count_depth), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/100)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	ggplot2::annotate("text",
		 label = "> Top 1% deepest", 
		 x = ceiling(quantile(features$feat_count_depth, .99)), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/100)),
			      plot = FALSE)$counts)/2,
		 hjust = 0,
		 color = "grey"
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Feature count depth"))

    g1_1 <- ggplot(features, aes(feat_count_depth)) +
	  geom_histogram(color = "grey", fill = "grey", 
			 bins = (ceiling(max(features$feat_count_depth)/100))) +
	  theme_classic() +
	  xlim(c(-0.5, ceiling(quantile(features$feat_count_depth, 0.5)))) +
	  ggtitle(paste(ds, "Feature count depth"),
		  subtitle = "Quantile 50")

    g2_wc <- ggplot(head(arrange(features, -feat_count_depth), nrow(features)*0.01), 
		    aes(label = gene_short_name, 
				size = feat_count_depth, 
				color = feat_count_depth)) +
		geom_text_wordcloud(rm_outside = TRUE) +
		scale_size_area(max_size = 2.5) +
		scale_color_gradient(low = "lightgrey", high = "black") +
		theme_minimal() +
		ggtitle("Top 1% count-deepest features")


    g2 <- ggplot(features, aes(feat_log10_count_depth)) +
	geom_histogram(color = "grey", fill = "grey", bins=100) +
	geom_vline(xintercept = round(quantile(features$feat_log10_count_depth, 0.5), digits = 2),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       round(mean(features$feat_log10_count_depth), digits=2)), 
		 x = mean(features$feat_log10_count_depth), 
		 y = max(hist(features$feat_log10_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_log10_count_depth),
					   by = max(features$feat_log10_count_depth)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "log10(feature count depth)"))

    g1 <- ggplot(`features`, aes(`feat_count_depth`)) +
	geom_histogram(color = "grey", fill = "grey", 
		       bins = ceiling(max(features$feat_count_depth)/200)) +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(features$feat_count_depth, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(features$feat_count_depth))), 
		 x = mean(features$feat_count_depth), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/200)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	ggplot2::annotate("text",
		 label = "> Top 1% deepest", 
		 x = ceiling(quantile(features$feat_count_depth, .99)), 
		 y = max(hist(features$feat_count_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_count_depth),
					   by = max(features$feat_count_depth)/ceiling(max(features$feat_count_depth)/100)),
			      plot = FALSE)$counts)/2,
		 hjust = 0,
		 color = "grey"
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Feature count depth"))


    g2_1 <- ggplot(features, aes(feat_log10_count_depth)) +
	  geom_histogram(color = "grey", fill = "grey", binwidth = 0.01) +
          xlim(c(-.1, quantile(features$feat_log10_count_depth, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(ds, "log10(feature count depth)"),
		  subtitle = "Quantile 50")

    g2_2 <- ggplot(`features`, aes(`feat_cell_depth`)) +
	geom_histogram(color = "grey", fill = "grey", 
		       bins = ceiling(max(features$feat_cell_depth)/100)) +
	geom_vline(xintercept = ceiling(quantile(features$feat_cell_depth, 0.5)),
		   linetype = "dashed", color = "black") +
	geom_vline(xintercept = ceiling(quantile(features$feat_cell_depth, .99)),
		   linetype = "dotted", color = "grey") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(features$feat_cell_depth))), 
		 x = mean(features$feat_cell_depth), 
		 y = max(hist(features$feat_cell_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_cell_depth),
					   by = max(features$feat_cell_depth)/ceiling(max(features$feat_cell_depth)/100)),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	ggplot2::annotate("text",
		 label = "> Top 1% deepest", 
		 x = ceiling(quantile(features$feat_cell_depth, .99)), 
		 y = max(hist(features$feat_cell_depth, 
			      breaks = seq(from = 0, 
					   to = max(features$feat_cell_depth),
					   by = max(features$feat_cell_depth)/ceiling(max(features$feat_cell_depth)/100)),
			      plot = FALSE)$counts)/2,
		 hjust = 0,
		 color = "grey"
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Feature cell depth"))

    g2_2_wc <- ggplot(head(arrange(features, -feat_cell_depth), nrow(features)*0.01), 
		    aes(label = gene_short_name, 
				size = feat_cell_depth, 
				color = feat_cell_depth)) +
		geom_text_wordcloud(rm_outside = TRUE) +
		scale_size_area(max_size = 2.5) +
		scale_color_gradient(low = "lightgrey", high = "black") +
		theme_minimal() +
		ggtitle("Top 1% cell-deepest features")



    # ------------------
    # Cell based density
    # ------------------

    logcount <- log10(colSums(exprs(cds)))
    cell <- data.frame('cell_depth_counts' = colSums(exprs(cds)),
                       'cell_log10_depth_counts' = ifelse(is.infinite(logcount),
							  0.1, logcount),
		       'cell_depth_genes_detected' = colSums(exprs(cds) != 0)
		       )

    cell[['barcode']] <- rownames(cell)
    print(head(cell))

    # ---
    # Mitochondrial gene information
    # ---
    
    if(length(mt_genes) == 0) {
	    mtc <- setNames(rep(0, ncol(cds)), colnames(cds))
    } else {
	    mtc <- colSums(exprs(cds)[mt_genes, ])
    }
    mt <- data.frame(mt_counts = mtc)
    mt$pct_mt_counts <- mt$mt_counts/cell$cell_depth_counts
    mt$barcode <- rownames(mt)
    print(head(mt))
    cell <- merge(cell, mt, by = "barcode")

    print(head(cell))

    # ---
    # Ribosomal
    # ---
    
    if(length(rb_genes) == 0) {
	    rbc <- setNames(rep(0, ncol(cds)), colnames(cds))
    } else {
	    rbc <- colSums(exprs(cds)[rb_genes, ])
    }
    rb <- data.frame(rb_counts = rbc)
    rb$pct_rb_counts <- rb$rb_counts/cell$cell_depth_counts
    rb$barcode <- rownames(rb)

    print("Rb genes")
    print(head(cell))
    cell <- merge(cell, rb, by = "barcode")
    print(head(cell))

    cell %>%
	dplyr::arrange(-cell_depth_counts) %>%
	mutate(`rank` = seq(nrow(cell))) -> cell
    print(head(cell,2))
   
    rownames(cell) <- cell$barcode

    g3 <- ggplot(cell, aes(cell_depth_counts)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 1000) +
	geom_vline(xintercept = ceiling(quantile(cell$cell_depth_counts, 0.5)),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(cell$cell_depth_counts))), 
		 x = mean(cell$cell_depth_counts), 
		 y = max(hist(cell$cell_depth_counts, 
			      breaks = seq(from = 0, 
					   to = max(cell$cell_depth_counts),
					   by = max(cell$cell_depth_counts)/1000),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Count depth"))

    g3_1 <- ggplot(cell, aes(cell_depth_counts)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 500) +
	  xlim(c(-.5, ceiling(quantile(cell$cell_depth_counts, 0.5)))) +
	  theme_classic() +
	  ggtitle(paste(ds, "Count depth - Quantile 50"))

    g4 <- ggplot(cell, aes(cell_log10_depth_counts)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
	geom_vline(xintercept = quantile(cell$cell_log10_depth_counts, 0.5),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       round(mean(cell$cell_log10_depth_counts), digits=2)), 
		 x = mean(cell$cell_log10_depth_counts), 
		 y = max(hist(cell$cell_log10_depth_counts, 
			      breaks = seq(from = 0, 
					   to = max(cell$cell_log10_depth_counts),
					   by = max(cell$cell_log10_depth_counts)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "log10(count depth)"))

    g4_1 <- ggplot(cell, aes(cell_log10_depth_counts)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 50) +
	  xlim(c(-.1, quantile(cell$cell_log10_depth_counts, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(ds, "log10(count depth) - Quantile 50"))

    g5 <- ggplot(cell, aes(cell_depth_genes_detected)) +
	geom_histogram(color = "#a6bddb", fill = "#a6bddb", bins = 100) +
	geom_vline(xintercept = quantile(cell$cell_depth_genes_detected, 0.5),
		   linetype = "dashed", color = "black") +
	ggplot2::annotate("text",
		 label = paste("mean ~", 
			       ceiling(mean(cell$cell_depth_genes_detected))), 
		 x = mean(cell$cell_depth_genes_detected), 
		 y = max(hist(cell$cell_depth_genes_detected, 
			      breaks = seq(from = 0, 
					   to = max(cell$cell_depth_genes_detected),
					   by = max(cell$cell_depth_genes_detected)/100),
			      plot = FALSE)$counts),
		 hjust = 0
		 ) +
	theme_classic() +
	ggtitle(paste(ds, "Gene depth"))

    g5_1 <- ggplot(cell, aes(cell_depth_genes_detected)) +
	  geom_histogram(color = "#a6bddb", fill = "#a6bddb", binwidth = 1) +
	  xlim(c(-.1, quantile(cell$cell_depth_genes_detected, 0.5))) +
	  theme_classic() +
	  ggtitle(paste(ds, "Gene depth - Quantile 50"))

    g6 <- ggplot(cell, aes(x = `rank`, y = cell_depth_counts)) +
	geom_point(aes(color = pct_mt_counts), 
		   stat = "identity", alpha = .2, size = 0.5) +
	scale_color_gradient(low = "black",  high = "red") +
	ylim(c(0, max(cell$cell_depth_counts))) +
	theme_classic() +
	ggtitle(paste(ds, "Depth vs. rank"))

    g6_1 <- ggplot(cell, aes(x = `rank`, y = cell_depth_counts)) +
	  geom_point(aes(color = pct_rb_counts), 
		     stat = "identity", alpha = .2, size = 0.5) +
	  scale_color_gradient(low = "pink",  high = "black") +
	  ylim(c(0, max(cell$cell_depth_counts))) +
	  theme_classic() +
	  ggtitle(paste(ds, "Depth vs. rank"))

    g7 <- ggplot(cell, aes(x = cell_depth_counts, y = cell_depth_genes_detected)) +
	geom_point(aes(color = pct_mt_counts), stat = "identity", 
		   alpha = 0.2, size = 0.5) +
	scale_color_gradient(low = "black",  high = "red") +
	theme_classic() +
	ggtitle(paste(ds, "Genes-detected vs. depth"),
		subtitle = "% Mitochondrial-counts")

    g8 <- ggplot(cell, aes(x = cell_depth_counts, y = cell_depth_genes_detected)) +
	geom_point(aes(color = pct_rb_counts, fill = pct_rb_counts), 
		   stat = "identity", alpha = 0.2, size = 0.5) +
	scale_color_gradient(high = "black",  low = "pink") +
	scale_fill_gradient(high = "black",  low = "pink") +
	theme_classic() +
	ggtitle(paste(ds, "Genes-detected vs. depth"),
	       subtitle = "% Ribosomal-counts")

    g9 <- ggplot(cell, aes(x = cell_depth_counts, y = pct_mt_counts)) +
	geom_point(color = "red", alpha = .2, size = 0.5) +
	geom_density_2d(color = "black", alpha = 0.5) +
	theme_classic() +
	ggtitle(paste(ds, "% Mitochondrial-counts vs. depth"))


    if(max(cell$pct_mt_counts) == 0) {
	    g10 <-  ggplot(cell, aes(pct_mt_counts)) +
		        geom_density(fill = "red", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_mt_counts),
				   linetype = "dashed", color = "black") + 
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Mitochondrial-counts density"))
    } else {
	    g10 <- ggplot(cell, aes(pct_mt_counts)) +
			geom_density(fill = "red", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_mt_counts),
				   linetype = "dashed", color = "black") +
			ggplot2::annotate("text",
					  label = paste("mean ~", 
							round(mean(cell$pct_mt_counts), digits=2)), 
					  	 x = mean(cell$pct_mt_counts), 
						 y = max(hist(cell$pct_mt_counts, 
							      breaks = seq(from = 0, 
									   to = max(cell$pct_mt_counts),
									   by = max(cell$pct_mt_counts)/1000),
							      plot = FALSE)$density),
					  hjust = 1,
					  vjust = 0
					  ) +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Mitochondrial-counts density"))
    }

    g11 <- ggplot(cell, aes(x = cell_depth_counts, y = pct_rb_counts)) +
	geom_point(alpha = 0.4, color = "pink", size = 0.5) +
	geom_density_2d(color = "black", alpha = 0.6) +
	theme_classic() +
	ggtitle(paste(ds, "% Ribosomal-counts vs. depth"))


    if(max(cell$pct_rb_counts) == 0) {

	    g12 <-  ggplot(cell, aes(pct_rb_counts)) +
		        geom_density(fill = "pink", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_rb_counts),
				   linetype = "dashed", color = "black") + 
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Ribosomal-counts density"))

    } else {

	    g12 <- ggplot(cell, aes(pct_rb_counts)) +
			geom_density(fill = "pink", color = "white", alpha = 0.8) +
			geom_vline(xintercept = median(cell$pct_rb_counts),
				   linetype = "dashed", color = "black") +
			ggplot2::annotate("text",
					  label = paste("mean ~", 
							round(mean(cell$pct_rb_counts), digits=2)), 
					  	 x = mean(cell$pct_rb_counts), 
						 y = max(hist(cell$pct_rb_counts, 
							      breaks = seq(from = 0, 
									   to = max(cell$pct_rb_counts),
									   by = max(cell$pct_rb_counts)/1000),
							      plot = FALSE)$density),
					  hjust = 1,
					  vjust = 0
					  ) +
			theme_classic() +
			coord_flip() +
			ggtitle(paste(ds, "Ribosomal-counts density"))

    }

    # ------------------------------
    # Save QC metrics in cds object
    # ------------------------------
    if(ds == "FILTERED") {
	    oldcols <- seq(ncol(fData(cds)))[-c(1:2)]
	    colnames(fData(cds))[oldcols] <- paste0(colnames(fData(cds))[oldcols], '_RAW')
	    oldcols <- seq(ncol(pData(cds)))[-c(1:2)]
	    colnames(pData(cds))[oldcols] <- paste0(colnames(pData(cds))[oldcols], '_RAW')
    }
    pData(cds)[["barcode"]] <- rownames(pData(cds))
    print(head(pData(cds)))

    cell <- merge(pData(cds), cell, by = "barcode")
    rownames(cell) <- cell$barcode
    pData(cds) <- cell[rownames(pData(cds)), ]
    #colnames(fData(cds))[1] <- "id"
    features <- merge(fData(cds), 
		      dplyr::select(features, -gene_short_name), 
		      by = "id")
    rownames(features) <- features$id
    fData(cds) <- features[rownames(cds), ]

    print(head(fData(cds)))

    return(list("ps1" = list(g1, g2_wc, g2, gtop20, g2_2, g2_2_wc, g3, g3_1, g4, g4_1),
		"ps2" = list(g5, g5_1, g6, g7, g9, g10, g6_1, g8, g11, g12),
		"cds" = cds))
}

# ---------------
# --- Parameters
# ---------------

if(ds == "RAW") {
	jobname <- gsub("\\.rds", "", rdsfile)
} else if(ds == "FILTERED") {
	jobname <- gsub("\\.rds", "", rdsfile)
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

if(sp == "mm") {

	
	if(colnames(fData(cds))[1] == "gene_id") {

		colnames(fData(cds))[1] <- "id"
	}

	print(head(fData(cds)))

	# Mitochondrial genes
	mt_genes <- fData(cds)[grepl("^mt-", fData(cds)$gene_short_name), ]$id
	message("Mitochondrial genes detected:")
	print(fData(cds)[grepl("^mt-", fData(cds)$gene_short_name), ]$gene_short_name)
	if(length(mt_genes) == 0) {
		message("No mitochondrial genes detected") 
	}
	
	# Ribosomal genes
	rb_genes <- fData(cds)[grepl("^Rps|^Rpl", fData(cds)$gene_short_name), ]$id
	message("Ribosomal genes detected:")
	print(fData(cds)[grepl("^Rps|^Rpl", fData(cds)$gene_short_name), 
	      ]$gene_short_name)
	if(length(rb_genes) == 0) {
		message("Ribosomal genes detected")
	}

} else if (sp == "hg") {

	# Mitochondrial genes
	mt_genes <- fData(cds)[grepl("^MT-", fData(cds)$gene_short_name), ]$id
	message("Mitochondrial genes detected:")
	print(fData(cds)[grepl("^MT-", fData(cds)$gene_short_name), 
	      ]$gene_short_name)
	if(length(mt_genes) == 0) {
		message("No mitochondrial genes detected") 
	}
	
	# Ribosomal genes
	rb_genes <- fData(cds)[grepl("^RPS|^RPL", fData(cds)$gene_short_name), ]$id
	message("Ribosomal genes detected:")
	print(fData(cds)[grepl("^RPS|^RPL", fData(cds)$gene_short_name), 
	      ]$gene_short_name)
	if(length(rb_genes) == 0) {
		message("Ribosomal genes detected")
	}

} else {

	stop("Please provide a valid specie parameter")

}

# --- Run
res <- plot_qc(cds, mt_genes, rb_genes, ds, jobname)

# --- Write QC plots
title <- ggdraw() + 
  draw_label(
    jobname,
    fontface = 'bold',
    #x = 0,
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 0)
  )

qc_file <- paste0(outputfolder, '/', jobname, '_', ds, '_mapping_qc.pdf')

pdf(qc_file, height = 11.69, width = 8.27)
    plot_grid(title, plot_grid(plotlist = res$ps1, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
    plot_grid(title, plot_grid(plotlist = res$ps2, ncol = 2), 
	      ncol = 1, rel_heights = c(0.02, 1))
dev.off()

message(paste("The QC report ", qc_file, "has been created."))

# --- Write cell_data_set object with QC metrics

cds_file <- paste0(outputfolder, '/', jobname, '_', ds, '_qc.rds')
saveRDS(res$cds, file = cds_file)

message(paste("The cell_data_set object (Monocle3)", 
	      cds_file, "has been created."))
