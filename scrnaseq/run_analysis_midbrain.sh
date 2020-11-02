#!/bin/bash
# ---------------------
# Reading raw data (CellRanger outs folder)
# ---------------------

jobname="PD_CO_hg_brain" 
infolder="./<CellRanger mapping results>"
outfolder="./results"

sp="hg"

IFS=$'\n'

cr2cds=(
	# IPD samples
	"Rscript ./src/read_10x.R --infolder=${infolder}/IPD1 --outfolder=${outfolder} --jobname=IPD1 --genome=${sp}"
	"Rscript ./src/read_10x.R --infolder=${infolder}/IPD2 --outfolder=${outfolder} --jobname=IPD2 --genome=${sp}"
	"Rscript ./src/read_10x.R --infolder=${infolder}/IPD3 --outfolder=${outfolder} --jobname=IPD3 --genome=${sp}"		
	"Rscript ./src/read_10x.R --infolder=${infolder}/IPD4 --outfolder=${outfolder} --jobname=IPD4 --genome=${sp}"
	"Rscript ./src/read_10x.R --infolder=${infolder}/IPD5 --outfolder=${outfolder} --jobname=IPD5 --genome=${sp}"
	# Control
	"Rscript ./src/read_10x.R --infolder=${infolder}/C1 --outfolder=${outfolder} --jobname=C1 --genome=${sp}"
	"Rscript ./src/read_10x.R --infolder=${infolder}/C2 --outfolder=${outfolder} --jobname=C2 --genome=${sp}"
	"Rscript ./src/read_10x.R --infolder=${infolder}/C3 --outfolder=${outfolder} --jobname=C3 --genome=${sp}"		
	"Rscript ./src/read_10x.R --infolder=${infolder}/C4 --outfolder=${outfolder} --jobname=C4 --genome=${sp}"
	"Rscript ./src/read_10x.R --infolder=${infolder}/C5 --outfolder=${outfolder} --jobname=C5 --genome=${sp}"
	"Rscript ./src/read_10x.R --infolder=${infolder}/C6 --outfolder=${outfolder} --jobname=C6 --genome=${sp}"
	)

parallel {} ::: ${cr2cds[*]}

# ---------------------
# QC
# ---------------------

infolder="./results"

dataset="RAW"

qcreport=(
	# IPD samples
	"Rscript ./src/qc_sc.R --scobject=IPD1.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD2.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD3.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD4.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD5.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	# Control
	"Rscript ./src/qc_sc.R --scobject=C1.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C2.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C3.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C4.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C5.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C6.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	)

parallel {} ::: ${qcreport[*]}

# ---------------------
# FILTERING
# ---------------------

minfeat=1
maxfeat="Inf"
mincellfeat=3
maxcellfeat="Inf"
mincountcell=1500
maxcountcell="Inf"
mingenecell=1000
maxgenecell="Inf"
pctmt=1
pctrb=1
mt="TRUE"
rb="TRUE"

filter=(
	# IPD samples
	"Rscript ./src/filter_qc.R --scobject=PD1_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=PD2_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=PD3_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=PD4_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=PD5_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	# Controls
	"Rscript ./src/filter_qc.R --scobject=C1_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=C2_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=C3_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=C4_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=C5_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	"Rscript ./src/filter_qc.R --scobject=C6_RAW_qc.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --minfeat=${minfeat} --maxfeat=${maxfeat} --mincellfeat=${mincellfeat} --maxcellfeat=${maxcellfeat} --mincountcell=${mincountcell} --maxcountcell=${maxcountcell} --mingenecell=${mingenecell} --maxgenecell=${maxgenecell} --pctmt=${pctmt} --pctrb=${pctrb} --mt=${mt} --rb=${rb}"
	)

parallel {} ::: ${filter[*]}

# -----------------------
# Duplets score Scrublet
# -----------------------

npcs=25
expected_db_rate=0.06
db_thr=0.25

duplets=(
	# IPD samples
	"python3 ./src/run_scrublet.py IPD1 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py IPD2 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py IPD3 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py IPD4 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py IPD5 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	# Controls
	"python3 ./src/run_scrublet.py C1 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py C2 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py C3 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py C4 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py C5 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	"python3 ./src/run_scrublet.py C6 ${infolder} ${outfolder} ${npcs} ${expected_db_rate} ${db_thr}"
	)

parallel {} ::: ${duplets[*]}

# ---------------------
# Filter duplets given a duplet-score Scrublet threshold
# ---------------------

dataset="FILTERED"
dupthr=0.15

filterdup=(
	# IPD samples
	"Rscript ./src/filter_duplets.R --scobject=IPD1_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=IPD2_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=IPD3_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=IPD4_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=IPD5_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	# Controls
	"Rscript ./src/filter_duplets.R --scobject=C1_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=C2_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=C3_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=C4_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=C5_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/filter_duplets.R --scobject=C6_filtered.rds --dupthr=${dupthr} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	)

parallel {} ::: ${filterdup[*]}

# ---------------------
# QC after filter
# ---------------------

dataset="FILTERED"

qcreport=(
	# IPD samples
	"Rscript ./src/qc_sc.R --scobject=IPD1_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD2_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD3_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD4_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=IPD5_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	# Controls
	"Rscript ./src/qc_sc.R --scobject=C1_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C2_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C3_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"	
	"Rscript ./src/qc_sc.R --scobject=C4_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C5_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"
	"Rscript ./src/qc_sc.R --scobject=C6_filtered.rds --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --dataset=${dataset}"				
	)

parallel {} ::: ${qcreport[*]}


# ---------------------
# Merge samples
# ---------------------

nhvg=4000
npcs=25
ncor=5
write="TRUE"
explore="TRUE"
testpcs="TRUE"
norm_method="SCT"
sample_sufix="_filtered_FILTERED_qc.rds"

# seurat3 step1
Rscript ./src/merge_sc_seurat3_step1.R --jobname=${jobname} --npcs=${npcs} --nhvg=${nhvg} --nneigh=30 --method=${norm_method} --ncores=${ncor} --testpcs=${testpcs} --specie=${sp} --write=${write} --infolder=${infolder} --outfolder=${outfolder} --explore=${explore} IPD1${sample_sufix} IPD2${sample_sufix} IPD3${sample_sufix} IPD4${sample_sufix} IPD5${sample_sufix} C1${sample_sufix} C2${sample_sufix} C3${sample_sufix} C4${sample_sufix} C5${sample_sufix} C6${sample_sufix}

samplemeta="./data/sample_metadata_curated.tsv"
sample_col="sample"
method="full_cca"
nhvg=4000
ref=NULL
norm_method="sct"
sample_sufix="_filtered_FILTERED_qc_"${npcs}"PCs_"${nhvg}"g_seurat3.rds"

# seurat step1.5
Rscript ./src/merge_sc_seurat3_step1_5.R --jobname=${jobname} --npcs=${npcs} --nhvg=${nhvg} --metafile=${samplemeta} --samplevar=${sample_col} --method=${method} --normethod=${norm_method} --ncores=${ncor} --testpcs=TRUE --infolder=${infolder} --outfolder=${outfolder} IPD1${sample_sufix} IPD2${sample_sufix} IPD3${sample_sufix} IPD4${sample_sufix} IPD5${sample_sufix} C1${sample_sufix} C2${sample_sufix} C3${sample_sufix} C4${sample_sufix} C5${sample_sufix} C6${sample_sufix}

# seurat3 step2
mcafile=${jobname}"_seurat3_integrated_npcs_"${npcs}"_nvgs_"${nhvg}"_"${method}"_"${norm_method}".rds"
Rscript ./src/merge_sc_seurat3_step2.R --jobname=${jobname} --mcafile=${mcafile} --npcs=${npcs} --nhvg=${nhvg} --colgroup=${sample_col} --ncores=${ncor} --testpcs=${testpcs} --infolder=${infolder} --outfolder=${outfolder}


# --------------------------------------------------
# Clustering
# --------------------------------------------------

mcafile=${jobname}"_seurat3_integrated_npcs_"${npcs}"_nvgs_"${nhvg}"_"${method}"_"${norm_method}"_seurat3_integrated_norm.rds"
assay="integrated"
Rscript ./src/cluster_seurat3.R --jobname=${jobname} --mcafile=${mcafile} --npcs=${npcs} --assay=${assay} --testpcs=FALSE --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}

# --------------------------------------------------
# Marker gene identification
# --------------------------------------------------

mcafile=${jobname}"_seurat3_integrated_npcs_"${npcs}"_nvgs_"${nhvg}"_"${method}"_"${norm_method}"_seurat3_integrated_norm.rds"
assay="RNA"
ncor=15
res=0.01
mincells=100
min_cell_prop=0.1
min_logfc=0

Rscript ./src/get_mkr_genes_seurat3.R --jobname=${jobname} --mcafile=${mcafile} --assay=${assay} --mincells=${mincells} --min_cell_prop=${min_cell_prop} --min_logfc=${min_logfc} --resolution=${res} --specie=${sp} --infolder=${infolder} --outfolder=${outfolder} --ncores=${ncor}

# --------------------------------------------------
# Differential cellular composition (two-level variables)
# --------------------------------------------------

mcafile=${jobname}"_seurat3_integrated_npcs_"${npcs}"_nvgs_"${nhvg}"_"${method}"_"${norm_method}"_seurat3_integrated_norm.rds"
groupvar="condition"
samplevar="sample"
res="cell_ontology"

Rscript ./src/diff_cellcomp.R --jobname=${jobname} --mcafile=${mcafile} --samplevar=${samplevar} --groupvar=${groupvar} --resolution=${res} --infolder=${infolder} --outfolder=${outfolder}

# --------------------------------------------------
# Differential gene expression (two-level variables)
# --------------------------------------------------

level2test="IPD"
distribution="quasipoisson"
ncor=15

Rscript ./src/diff_exp_seurat2monocle.R --jobname=${jobname} --mcafile=${mcafile} --groupvar=${groupvar} --level2test=${level2test} --samplevar=${samplevar} --resolution=${res} --distribution=${distribution} --det_thr=${min_cell_prop} --ncores=${ncor} --specie=${sp} --infolder=${infolder} --outfolder=${outfolder}

# --------------------------------------------------
# Sub-clustering (split dataset into the comprising cell-types
# --------------------------------------------------



# --------------------------------------------------
# Trajectory inference
# --------------------------------------------------



# 
