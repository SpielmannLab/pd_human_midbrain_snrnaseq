#
# Read 10X - CellRanger output and write and SingleCellExperiment object
#

"CellRanger output to cell_data_set (Monocle3)

Usage: read_10x.R --infolder=<folder> --outfolder=<folder> --jobname=<value> --genome=<value>

Options:
  -h --help            Show this screen.
  --version            00.99.01
  --infolder=<file>    Path to ../outs parent folder.
  --outfolder=<file>   Path to results folder.
  --jobname=<value>    Descriptive name.
  --genome=<value>     Genome string eg. mm10.

"-> doc

library(docopt)
arguments<-docopt(doc, quoted_args=TRUE)
print(arguments)

# --- Dependencies

pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
	  'limma', 'S4Vectors', 'SingleCellExperiment',
	  'SummarizedExperiment', 'batchelor', 'devtools')

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

# --- Parameters

inputfolder <- arguments$infolder 
outputfolder <- arguments$outfolder
jobname <- arguments$jobname
`genome` <- arguments$genome
umi_cutoff <- 0

# --- Run

cds <- load_cellranger_data(inputfolder)

cds_file <- paste0(outputfolder, '/', jobname, '.rds')

saveRDS(cds, file = cds_file)

message(paste("The cell_data_set object (Monocle3)", 
	      cds_file, "has been created."))
