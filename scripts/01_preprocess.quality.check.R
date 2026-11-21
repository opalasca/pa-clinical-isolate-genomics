# Assembly quality check and sample metadata --------------------------------
# Combines the sample sheet, MLST types and QUAST assembly statistics for the
# P. aeruginosa isolates (two sequencing duplicates per isolate).
# Assemblies with a total length > 8 Mb are flagged (see Methods).

# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages-------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

# Send on-screen plots to a null device when run non-interactively (avoids Rplots.pdf)
if (!interactive()) pdf(NULL)

# Read input files ---------

sample.info <- read_excel(paste(datadir,"WGS\ samplesheet.xlsx",sep=""))
sample.info <- as.data.frame(sample.info)

samples.isolates <- sample.info
samples.isolates$isolate <- sub("[^0-9].*$", "", samples.isolates$isolate.replicate)

samples.isolates = samples.isolates[samples.isolates$isolate.replicate!="blank",]

samples.isolates$sample.id <- as.integer(sub(".*x", "", samples.isolates$sample))
samples.isolates$sample.id <- paste0("S", samples.isolates$sample.id)


# Read MLST type and save an xlsx with it ---------
mlst <- read.csv(paste0(datadir,"MLST_results_p_aeruginosa/MLST.output.all.samples.txt"), header = FALSE, sep=" ")
mlst <- mlst[c(1,3)]
colnames(mlst) <- c("sample.id","MLST_type")
mlst$sample.id <- sub("_.*$", "", mlst$sample.id)
samples.isolates.plus.mlst <- merge(samples.isolates, mlst, by="sample.id", all.x=TRUE)
write.xlsx(samples.isolates.plus.mlst, paste0(resultdir,"metadata_with_MLST_types_p_aeruginosa.xlsx"))


# Read quast report ---------
quast_report <- read.csv(paste0(datadir,"quast_transposed_report.tsv"), header = TRUE, sep="\t")
quast_report$sample <- sub("_.*$", "", quast_report$Assembly)

quast_report$Total.length.....0.bp. <- as.numeric(gsub(",", "", quast_report$Total.length.....0.bp.))
quast <- quast_report[c("sample","N50","Total.length","Total.length.....10000.bp.")]
colnames(quast)[4]="Total.length.10kbp"


# Merge QC with sample info and flag oversized assemblies ---------
samples.isolates.QC <- merge(samples.isolates.plus.mlst, quast, by.x="sample", by.y="sample", all.x=TRUE)
samples.isolates.QC$isolate <- as.numeric(samples.isolates.QC$isolate)
samples.isolates.QC <- samples.isolates.QC[order(samples.isolates.QC$isolate), ]
samples.isolates.QC$assembly.over.8Mb <- samples.isolates.QC$Total.length > 8e6

write.xlsx(samples.isolates.QC, file=paste0(resultdir,"assembly_QC_p_aeruginosa.xlsx"), rowNames=FALSE)


# Write list of P. aeruginosa samples ---------
samples.paeruginosa <- samples.isolates[samples.isolates$Species =="Pseudomonas aeruginosa", ]
samples.paeruginosa$id <- paste(samples.paeruginosa$sample, samples.paeruginosa$sample.id, sep="_")
samples.paeruginosa <- samples.paeruginosa[c("id")]
write.csv(samples.paeruginosa, file=paste0(resultdir,"samples_paeruginosa.txt"), row.names=FALSE, quote=FALSE)
