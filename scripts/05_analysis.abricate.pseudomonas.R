# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/05_virulence_genes/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures, recursive = TRUE), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages-------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx", "pheatmap","tibble","ComplexHeatmap","circlize")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

# Send on-screen plots to a null device when run non-interactively (avoids Rplots.pdf)
if (!interactive()) pdf(NULL)

# Read and process metadata ---------

sample.info <- read_excel(paste(datadir,"WGS\ samplesheet.xlsx",sep=""))
sample.info <- as.data.frame(sample.info)

samples.isolates <- sample.info
samples.isolates$isolate <- sub("[^0-9].*$", "", samples.isolates$isolate.replicate)

samples.isolates = samples.isolates[samples.isolates$isolate.replicate!="blank",] 
samples.isolates$sample.id <- as.integer(sub(".*x", "", samples.isolates$sample))
samples.isolates$sample.id <- paste0("S", samples.isolates$sample.id)

# Read MLST type 
mlst <- read.csv(paste0(datadir,"MLST_results_p_aeruginosa/MLST.output.all.samples.txt"), header = FALSE, sep=" ")
mlst <- mlst[c(1,3)]
colnames(mlst) <- c("sample.id","MLST_type")
mlst$sample.id <- sub("_.*$", "", mlst$sample.id)
samples.isolates <- merge(samples.isolates, mlst, by="sample.id", all.x=TRUE)

# Remove samples with too long genome size (therefore suspected contamination)
samples_to_remove = c("S68", "S69", "S39","S54")

# Define order based on phylogenetic tree (tip_order)
isolate_keep_order <-  c("32", "25", "22", "4", "3", "18", "5", "27", "11", "33", "19", "9", "26", "PAO1")

isolate_order <- isolate_keep_order

# Read and process abricate result file -------------------
# Read all abricate files, including those for reference genomes
abricatedir <- paste0(datadir,"ABRicate/")
files <- list.files(abricatedir, full.names = TRUE)

abricate <- lapply(files, function(f) {
  read.delim(f, header = TRUE, stringsAsFactors = FALSE)
})
abricate <- do.call(rbind, abricate)

abricate$sampleID <- sub("_.*$", "", abricate$X.FILE)
# Filter for coverage >= 80%
abricate = abricate[abricate$X.COVERAGE >=80, ]

abricate = merge(abricate, samples.isolates, by.x="sampleID",by.y="sample.id", all.x=TRUE)

# For reference genomes, set isolate to sampleID and Species to P aeruginosa (for consistency with our isolates)
ref.genomes =c("PAO1","PA7","PA14","LESB5","DK2")
ref.genomes.remove = c("PA7","PA14","LESB5","DK2")

abricate[abricate$sampleID %in% ref.genomes, ]$isolate <- abricate[abricate$sampleID %in% ref.genomes, ]$sampleID
abricate[abricate$sampleID %in% ref.genomes, ]$Species <- "Pseudomonas aeruginosa"

# Add column indicating whether gene is present in both duplicates of an isolate
abricate <- abricate %>%
  group_by(isolate, GENE, PRODUCT) %>%
  mutate(Present.in.duplicate = ifelse(n_distinct(sampleID) > 1, "yes", "no")) %>%
  ungroup()

abricate <- abricate[!abricate$isolate %in% ref.genomes.remove, ]  

# Remove genes not present in both duplicates for selected samples
abricate <- abricate %>%
  filter(!((Present.in.duplicate == "no") & (isolate!="PAO1")))

# Merge duplicates of the same isolate (take union of genes across duplicates, by keeping only first occurance when a gene is predicted in both duplicates)
abricate <- abricate %>%
  group_by(isolate, GENE, PRODUCT) %>%
  dplyr::slice(1) %>%        # keep one row per gene per isolate
  ungroup()

# Filter for Pseudomonas
abricate = abricate[which(abricate$Species == "Pseudomonas aeruginosa"),  ]

# Simplify data by keeping only the strictly necessary columns 
abricate <- abricate[c("isolate","GENE","PRODUCT","MLST_type")]

abricate[abricate$isolate=="PAO1",]$MLST_type <- "549"
#abricate[abricate$isolate=="PA14",]$MLST_type <- "253"
#abricate[abricate$isolate=="LESB5",]$MLST_type <- "146"
#abricate[abricate$isolate=="DK2",]$MLST_type <- "386"
#abricate[abricate$isolate=="PA7",]$MLST_type <- "1195"

abricate$MLST_type = paste0("ST",abricate$MLST_type)

# Colour for gene presence (blue used for presence/absence in all manuscript figures; green is used for mutations)
presence.col <- "#95AAD3"

# Extract functional annotation from PRODUCT column ----------
abricate$func.annot.full <- str_extract(abricate$PRODUCT, "\\[[^\\]]+\\]")
abricate$func.annot.full <- str_remove_all(abricate$func.annot.full, "[\\[\\]]")

abricate$func.annot <- str_replace(abricate$func.annot.full, "\\s*\\(VF[0-9]+\\)", "")
abricate$func.annot <- str_replace(abricate$func.annot, "\\s*\\(AI[0-9]+\\)", "")

abricate_annotation <- unique(abricate[c("GENE","func.annot.full","func.annot")])

# Heatmap of all virulence genes, grouped by VFDB functional category (supplementary overview) -----

# Build presence/absence matrix ---
gene_mat <- abricate %>%
  mutate(present = 1) %>%
  distinct(GENE, isolate, present) %>%
  pivot_wider(
    names_from = isolate,
    values_from = present,
    values_fill = 0
  ) %>%
  column_to_rownames("GENE") %>%
  as.matrix()

# Reorder isolates ---

isolate_keep <- isolate_keep_order

gene_mat <- gene_mat[, intersect(isolate_keep, colnames(gene_mat)), drop = FALSE]

# Functional annotation per gene ---
gene_func <- abricate %>%
  distinct(GENE, func.annot) %>%
  arrange(GENE)

gene_func <- gene_func[match(rownames(gene_mat), gene_func$GENE), ]
stopifnot(all(gene_func$GENE == rownames(gene_mat)))

# Identify unique functional groups ---
func_groups <- unique(gene_func$func.annot)

# Make sure rownames match
stopifnot(all(rownames(gene_mat) == gene_func$GENE))

# Order rows by functional annotation alphabetically
gene_order <- gene_func %>%
  arrange(func.annot, GENE) %>%  # first by functional annotation, then by gene name
  pull(GENE)

gene_mat_ordered <- gene_mat[gene_order, ]

# Make a factor for row titles in the same order
row_titles_ordered <- factor(gene_func$func.annot[match(rownames(gene_mat_ordered), gene_func$GENE)],
                             levels = sort(unique(gene_func$func.annot)))

ht_abricate <- Heatmap(
  gene_mat_ordered,
  name = "Presence",
  col = c("0" = "white", "1" = presence.col),
  cluster_rows = FALSE,            # no clustering
  cluster_columns = FALSE,
  row_split = row_titles_ordered,  # split by function
  cluster_row_slices = FALSE,      # don't cluster within group
  show_row_names = FALSE,
  row_title = levels(row_titles_ordered),
  row_title_gp = gpar(fontsize = 12),
  row_title_rot = 0,           # horizontal
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 13),
  row_names_gp = gpar(fontsize = 8),
  rect_gp = gpar(col = "grey60", lwd = 0.3)
)

pdf(paste0(figures, "heatmap_virulence_all_genes.pdf"), width = 10, height = 17)
draw(ht_abricate)
dev.off()

png(filename = paste0(figures, "heatmap_virulence_all_genes.png"),width = 3000,height = 4800,res = 300)
draw(ht_abricate)
dev.off()

# Heatmap of virulence genes with variable presence (genes found in all isolates removed), grouped by VFDB category -----
# Basis for Fig. 4A; the gene subset shown in the published figure was selected manually

# Build presence/absence matrix ---
gene_mat <- abricate %>%
  mutate(present = 1) %>%
  distinct(GENE, isolate, present) %>%
  pivot_wider(
    names_from = isolate,
    values_from = present,
    values_fill = 0
  ) %>%
  column_to_rownames("GENE") %>%
  as.matrix()

# Reorder isolates ---

isolate_keep <- isolate_keep_order

gene_mat <- gene_mat[, intersect(isolate_keep, colnames(gene_mat)), drop = FALSE]

# Remove genes found in all isolates
gene_mat <- gene_mat[rowSums(gene_mat == 1) < ncol(gene_mat), ]

# Functional annotation per gene ---
gene_func <- abricate %>%
  distinct(GENE, func.annot) %>%
  arrange(GENE)

gene_func <- gene_func[match(rownames(gene_mat), gene_func$GENE), ]
stopifnot(all(gene_func$GENE == rownames(gene_mat)))

# Identify unique functional groups ---
func_groups <- unique(gene_func$func.annot)

# Make sure rownames match
stopifnot(all(rownames(gene_mat) == gene_func$GENE))

# Order rows by functional annotation alphabetically
gene_order <- gene_func %>%
  arrange(func.annot, GENE) %>%  # first by functional annotation, then by gene name
  pull(GENE)

gene_mat_ordered <- gene_mat[gene_order, ]

# Make a factor for row titles in the same order
row_titles_ordered <- factor(gene_func$func.annot[match(rownames(gene_mat_ordered), gene_func$GENE)],
                             levels = sort(unique(gene_func$func.annot)))

ht_abricate <- Heatmap(
  gene_mat_ordered,
  name = "Presence",
  col = c("0" = "white", "1" = presence.col),
  cluster_rows = FALSE,            # no clustering
  cluster_columns = FALSE,
  row_split = row_titles_ordered,  # split by function
  cluster_row_slices = FALSE,      # don't cluster within group
  show_row_names = TRUE,
  row_title = levels(row_titles_ordered),
  row_title_gp = gpar(fontsize = 11),
  row_title_rot = 0,           # horizontal
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 13),
  row_names_gp = gpar(fontsize = 9),
  rect_gp = gpar(col = "grey60", lwd = 0.3)
)

pdf(paste0(figures, "heatmap_virulence_variable_genes.pdf"), width = 8.5, height = 9)
draw(ht_abricate)
dev.off()

# Save as xlsx ------

abricate.in.all.isolates <- abricate %>%
  group_by(GENE, PRODUCT, func.annot) %>%
  filter(n_distinct(isolate) == 14 ) %>%
  ungroup()

abricate.not.in.all.isolates <- abricate %>%
  group_by(GENE, PRODUCT, func.annot) %>%
  filter(n_distinct(isolate) < 14 ) %>%
  ungroup()

abricate.wide <- abricate.not.in.all.isolates %>%
  mutate(present = 1) %>%
  select(GENE, PRODUCT, func.annot,isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = func.annot) %>%
  arrange(func.annot, GENE)

abricate.wide.found.in.all <- abricate.in.all.isolates %>%
  mutate(present = 1) %>%
  select(GENE, PRODUCT, func.annot,isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = func.annot) %>%
  arrange(func.annot, GENE)

wb <- createWorkbook()
addWorksheet(wb, "Virulence genes")
writeData(wb, "Virulence genes", abricate.wide)

addWorksheet(wb, "Virulence genes found in all")
writeData(wb, "Virulence genes found in all", abricate.wide.found.in.all)

saveWorkbook(wb, "results/abricate_virulence_results.xlsx", overwrite = TRUE)
