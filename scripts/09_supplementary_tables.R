# supplementary_tables.R
# Generates supplementary tables for the manuscript.
# Outputs:
#   results/Supplementary_Resistance.xlsx  (tabs S1-S6)
#   results/Supplementary_Virulence.xlsx   (tabs S7-S8)

# Setup ----
dir <- "./"  # run all scripts from the repository root
datadir  <- paste0(dir, "data/")
resultdir <- paste0(dir, "results/")
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

list.of.packages <- c("dplyr", "tidyr", "readxl", "stringr", "readr", "openxlsx", "tibble")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, library, character.only = TRUE)

# Send on-screen plots to a null device when run non-interactively (avoids Rplots.pdf)
if (!interactive()) pdf(NULL)

# Metadata ----
sample.info <- read_excel(paste0(datadir, "WGS samplesheet.xlsx"))
sample.info <- as.data.frame(sample.info)
samples.isolates <- sample.info
samples.isolates$isolate <- sub("[^0-9].*$", "", samples.isolates$isolate.replicate)
samples.isolates <- samples.isolates[samples.isolates$isolate.replicate != "blank", ]
samples.isolates$sample.id <- as.integer(sub(".*x", "", samples.isolates$sample))
samples.isolates$sample.id <- paste0("S", samples.isolates$sample.id)

mlst <- read.csv(paste0(datadir, "MLST_results_p_aeruginosa/MLST.output.all.samples.txt"),
                 header = FALSE, sep = " ")
mlst <- mlst[c(1, 3)]
colnames(mlst) <- c("sample.id", "mlst_type")
mlst$sample.id <- sub("_.*$", "", mlst$sample.id)
samples.isolates <- merge(samples.isolates, mlst, by = "sample.id", all.x = TRUE)

samples_to_remove  <- c("S68", "S69", "S39", "S54")
isolate_keep_order <- c("32", "25", "22", "4", "3", "18", "5", "27", "11", "33", "19", "9", "26", "PAO1")
ref.genomes        <- c("PAO1", "PA7", "PA14", "LESB5", "DK2")
ref.genomes.remove <- c("PA7", "PA14", "LESB5", "DK2")

# Gene name harmonisation map (ResFinder/AMRFinder -> RGI names) ----
resfinder_to_rgi <- c(
  "aph(3')-VIa" = "aPH(3')-VIa",
  "aph(3'')-Ib" = "aPH(3'')-Ib",
  "aph(6)-Id"   = "aPH(6)-Id",
  "mph(E)"      = "mphE",
  "mph(F)"      = "mphF",
  "msr(E)"      = "msrE",
  "blaIMP-1b"   = "blaIMP-1",
  "blaTEM-1A"   = "blaTEM-1",
  "blaTEM-1B"   = "blaTEM-1",
  "ARR-3"       = "arr-3",
  "aadA1"       = "aadA",
  "aadA6"       = "aadA"
)

harmonise_gene_name <- function(x) {
  ifelse(x %in% names(resfinder_to_rgi), resfinder_to_rgi[x], x)
}

# Clean RGI gene names (add bla prefix, lower-case first letter) ----
clean_rgi_gene_name <- function(name) {
  bla.pattern <- "^(OXA|TEM|CTX-M|NDM|IMP|PDC|CMY|SHV|KPC|SRT|PER)"
  name <- ifelse(grepl(bla.pattern, name), paste0("bla", name), name)
  name <- gsub("Pseudomonas aeruginosa ", "", name)
  name <- paste0(tolower(substr(name, 1, 1)), substr(name, 2, nchar(name)))
  dplyr::case_when(
    grepl("gyrA", name) ~ "gyrA",
    grepl("EF-Tu", name) ~ "EF-Tu",
    TRUE ~ name
  )
}

# Helper: write a worksheet with a description row, blank row, then data ----
write_sheet <- function(wb, sheet_name, description, data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name,
            x = data.frame(Description = description),
            startRow = 1, colNames = FALSE)
  writeData(wb, sheet_name, x = data, startRow = 3, colNames = TRUE)
  addStyle(wb, sheet_name,
           createStyle(fontColour = "#666666", fontSize = 10, wrapText = TRUE),
           rows = 1, cols = 1)
  addStyle(wb, sheet_name,
           createStyle(textDecoration = "bold", border = "Bottom"),
           rows = 3, cols = seq_len(ncol(data)))
  setColWidths(wb, sheet_name, cols = 1, widths = 80)
}

# ============================================================
# RGI ----
# ============================================================
rgidir <- paste0(datadir, "RGI/")
files  <- list.files(rgidir, pattern = "fa(sta)?\\.txt$", full.names = TRUE)

rgi <- do.call(rbind, lapply(files, function(f) {
  df    <- read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  fname <- basename(f)
  df$sampleID <- if (grepl("_(S[0-9]+)_contigs", fname)) {
    sub(".*_(S[0-9]+)_contigs.*$", "\\1", fname)
  } else {
    sub(".*_([A-Za-z0-9]+)_genome.*$", "\\1", fname)
  }
  df
}))

rgi <- merge(rgi, samples.isolates, by.x = "sampleID", by.y = "sample.id", all.x = TRUE)
rgi[rgi$sampleID %in% ref.genomes, ]$isolate  <- rgi[rgi$sampleID %in% ref.genomes, ]$sampleID
rgi[rgi$sampleID %in% ref.genomes, ]$Species  <- "Pseudomonas aeruginosa"
rgi <- rgi[!rgi$sampleID %in% ref.genomes.remove, ]
rgi[rgi$isolate == "PAO1", ]$mlst_type <- "549"

rgi <- rgi %>%
  group_by(isolate, Best_Hit_ARO, Model_type, SNPs_in_Best_Hit_ARO) %>%
  mutate(confirmed_in_duplicate = ifelse(n_distinct(sampleID) > 1, "yes", "no")) %>%
  ungroup() %>%
  filter(!(confirmed_in_duplicate == "no" & sampleID %in% samples_to_remove)) %>%
  group_by(isolate, Best_Hit_ARO, Model_type, SNPs_in_Best_Hit_ARO) %>%
  dplyr::slice(1) %>%
  ungroup()

rgi <- rgi[rgi$isolate %in% isolate_keep_order, ]

rgi <- rgi %>%
  mutate(gene_name = clean_rgi_gene_name(Best_Hit_ARO)) %>%
  filter(!gene_name %in% "escherichia coli soxS with mutation conferring antibiotic resistance") %>%
  rename(
    aro_name             = Best_Hit_ARO,
    model_type           = Model_type,
    snps                 = SNPs_in_Best_Hit_ARO,
    drug_class           = Drug.Class,
    resistance_mechanism = Resistance.Mechanism,
    amr_gene_family      = AMR.Gene.Family,
    sequence_identity    = Best_Identities,
    sequence_coverage    = Percentage.Length.of.Reference.Sequence,
    antibiotic           = Antibiotic,
    rgi_cutoff           = Cut_Off
  ) %>%
  select(isolate, mlst_type, gene_name, aro_name, model_type, snps,
         drug_class, resistance_mechanism, amr_gene_family,
         antibiotic, rgi_cutoff, sequence_identity, sequence_coverage,
         confirmed_in_duplicate)

rgi.acquired  <- rgi %>% filter(model_type == "protein homolog model")
rgi.mutations <- rgi %>% filter(model_type %in% c("protein variant model", "protein overexpression model"))

rgi.mutations.wide <- rgi.mutations %>%
  mutate(present = 1) %>%
  select(gene_name, snps, model_type, drug_class, resistance_mechanism, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(intersect(isolate_keep_order, names(.))), .after = resistance_mechanism) %>%
  arrange(drug_class, gene_name)

# ============================================================
# ResFinder ----
# ============================================================
resfinder <- read.csv(paste0(datadir, "Resfinder_plasmidfinder.tsv.tabular"),
                      header = TRUE, sep = "\t")
resfinder$sample <- sub("_.*$", "", resfinder$sampleID)
resfinder <- merge(resfinder, samples.isolates, by.x = "sample", by.y = "sample.id", all.x = TRUE)
resfinder <- resfinder[resfinder$Species == "Pseudomonas aeruginosa", ]
resfinder <- resfinder[c("sample", "isolate", "data", "data_type",
                         "CGE_predicted_phenotype", "X.identity", "X.overlap", "mlst_type")]
colnames(resfinder)[3:7] <- c("gene_name", "data_type", "predicted_phenotype",
                               "sequence_identity", "sequence_coverage")

resfinder.ref <- read.csv(paste0(datadir, "Staramr_reference_genomes.tabular"),
                          header = TRUE, sep = "\t")
resfinder.ref$sample  <- sub("_.*$", "", resfinder.ref$Isolate.ID)
resfinder.ref$isolate <- resfinder.ref$sample
resfinder.ref <- resfinder.ref[!resfinder.ref$isolate %in% ref.genomes.remove, ]
resfinder.ref <- resfinder.ref[c("sample", "isolate", "Data", "Data.Type",
                                  "CGE.Predicted.Phenotype", "X.Identity", "X.Overlap")]
resfinder.ref$mlst_type <- NA
resfinder.ref[resfinder.ref$isolate == "PAO1", ]$mlst_type <- "549"
colnames(resfinder.ref)[3:7] <- c("gene_name", "data_type", "predicted_phenotype",
                                   "sequence_identity", "sequence_coverage")
resfinder <- bind_rows(resfinder, resfinder.ref)

# Fix aac(6')-Ib entry
resfinder[which(resfinder$sample == "S19" & resfinder$gene_name == "aac(6')-Ib"),
          "predicted_phenotype"] <-
  "Tobramycin, Dibekacin, Amikacin, Sisomicin, Netilmicin, Fluoroquinolone, Ciprofloxacin"
resfinder[which(resfinder$sample == "S19" & resfinder$gene_name == "aac(6')-Ib"),
          "gene_name"] <- "aac(6')-Ib-cr"

resfinder <- resfinder %>%
  group_by(isolate, gene_name, data_type, predicted_phenotype) %>%
  mutate(confirmed_in_duplicate = ifelse(n_distinct(sample) > 1, "yes", "no")) %>%
  ungroup() %>%
  filter(!(confirmed_in_duplicate == "no" & sample %in% samples_to_remove)) %>%
  group_by(isolate, gene_name, predicted_phenotype) %>%
  dplyr::slice(1) %>%
  ungroup()

resfinder.resistance <- resfinder %>%
  filter(data_type == "Resistance") %>%
  mutate(gene_name_harmonised = harmonise_gene_name(gene_name))

# ============================================================
# AMRFinder ----
# ============================================================
amrfinderdir <- paste0(datadir, "AMRFinderPlus/")
files <- list.files(amrfinderdir, full.names = TRUE)

amrfinder <- do.call(rbind, lapply(files, function(f) {
  df <- read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  df$sampleID <- sub("_contigs.*$", "", basename(f))
  if (!"Hierarchy.node" %in% colnames(df)) df$Hierarchy.node <- NA
  df
}))

amrfinder <- amrfinder[amrfinder$X..Coverage.of.reference.sequence >= 80, ]
amrfinder <- merge(amrfinder, samples.isolates, by.x = "sampleID", by.y = "sample.id", all.x = TRUE)
amrfinder[amrfinder$sampleID == "PAO1", "Species"]    <- "Pseudomonas aeruginosa"
amrfinder[amrfinder$sampleID == "PAO1", "isolate"]    <- "PAO1"
amrfinder[amrfinder$sampleID == "PAO1", "mlst_type"]  <- "549"
amrfinder <- amrfinder[amrfinder$Species == "Pseudomonas aeruginosa", ]

amrfinder <- amrfinder %>%
  group_by(isolate, Gene.symbol, Element.type, Element.subtype, Method, Class, Subclass) %>%
  mutate(confirmed_in_duplicate = ifelse(n_distinct(sampleID) > 1, "yes", "no")) %>%
  ungroup() %>%
  filter(!(confirmed_in_duplicate == "no" & sampleID %in% samples_to_remove)) %>%
  group_by(isolate, Gene.symbol, Element.type, Element.subtype, Class, Subclass) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  rename(
    gene_name        = Gene.symbol,
    element_type     = Element.type,
    element_subtype  = Element.subtype,
    detection_method = Method,
    sequence_name    = Sequence.name,
    drug_class       = Class,
    drug_subclass    = Subclass
  ) %>%
  mutate(gene_name_harmonised = harmonise_gene_name(gene_name)) %>%
  select(isolate, mlst_type, gene_name, gene_name_harmonised, element_type, element_subtype,
         detection_method, sequence_name, drug_class, drug_subclass, confirmed_in_duplicate)

amrfinder.acquired  <- amrfinder %>% filter(detection_method %in% c("EXACTX", "ALLELEX", "BLASTX"))
amrfinder.mutations <- amrfinder %>% filter(detection_method == "POINTX")

amrfinder.mutations.wide <- amrfinder.mutations %>%
  mutate(present = 1) %>%
  select(gene_name, sequence_name, drug_class, drug_subclass, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(intersect(isolate_keep_order, names(.))), .after = drug_subclass) %>%
  arrange(drug_class, gene_name)

# ============================================================
# Keep acquired genes confirmed in both duplicates by at least one tool ----
# ============================================================
# A gene is considered present in an isolate only if at least one tool (RGI,
# ResFinder or AMRFinderPlus) detected it in both sequencing duplicates (see Methods).
# PAO1 is a single reference genome and is kept as is.
confirmed.pairs <- bind_rows(
  rgi.acquired         %>% transmute(isolate, gene_name_harmonised = gene_name, confirmed_in_duplicate),
  resfinder.resistance %>% transmute(isolate, gene_name_harmonised, confirmed_in_duplicate),
  amrfinder.acquired   %>% transmute(isolate, gene_name_harmonised, confirmed_in_duplicate)
) %>%
  group_by(isolate, gene_name_harmonised) %>%
  filter(isolate == "PAO1" | any(confirmed_in_duplicate == "yes")) %>%
  ungroup() %>%
  distinct(isolate, gene_name_harmonised)

resfinder.resistance.all <- resfinder.resistance  # all calls, kept for the long-format table (Tab S5)

rgi.acquired <- rgi.acquired %>%
  semi_join(confirmed.pairs, by = c("isolate", "gene_name" = "gene_name_harmonised"))
resfinder.resistance <- resfinder.resistance %>%
  semi_join(confirmed.pairs, by = c("isolate", "gene_name_harmonised"))
amrfinder.acquired <- amrfinder.acquired %>%
  semi_join(confirmed.pairs, by = c("isolate", "gene_name_harmonised"))

# Per-tool presence/absence tables (Tabs S2-S4) ----
resfinder.wide <- resfinder.resistance %>%
  mutate(present = 1) %>%
  select(gene_name, predicted_phenotype, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(intersect(isolate_keep_order, names(.))), .after = predicted_phenotype) %>%
  arrange(predicted_phenotype, gene_name)

rgi.acquired.wide <- rgi.acquired %>%
  mutate(present = 1) %>%
  select(gene_name, drug_class, resistance_mechanism, amr_gene_family, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(intersect(isolate_keep_order, names(.))), .after = amr_gene_family) %>%
  arrange(drug_class, gene_name)

amrfinder.acquired.wide <- amrfinder.acquired %>%
  mutate(present = 1) %>%
  select(gene_name, element_type, element_subtype, drug_class, drug_subclass, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(intersect(isolate_keep_order, names(.))), .after = drug_subclass) %>%
  arrange(element_type, drug_class, gene_name)

# ============================================================
# Combined wide (Tab S1) ----
# ============================================================
rgi.pairs <- rgi.acquired %>%
  select(isolate, gene_name) %>% distinct() %>%
  mutate(detected_by_rgi = "yes")

resfinder.pairs <- resfinder.resistance %>%
  select(isolate, gene_name = gene_name_harmonised) %>% distinct() %>%
  mutate(detected_by_resfinder = "yes")

amrfinder.pairs <- amrfinder.acquired %>%
  select(isolate, gene_name = gene_name_harmonised) %>% distinct() %>%
  mutate(detected_by_amrfinder = "yes")

combined.long <- rgi.pairs %>%
  full_join(resfinder.pairs, by = c("isolate", "gene_name")) %>%
  full_join(amrfinder.pairs, by = c("isolate", "gene_name")) %>%
  mutate(
    detected_by_rgi       = replace_na(detected_by_rgi, "no"),
    detected_by_resfinder = replace_na(detected_by_resfinder, "no"),
    detected_by_amrfinder = replace_na(detected_by_amrfinder, "no")
  ) %>%
  filter(isolate %in% isolate_keep_order)

gene.tool.summary <- combined.long %>%
  group_by(gene_name) %>%
  summarise(
    detected_by_rgi       = ifelse(any(detected_by_rgi == "yes"), "yes", "no"),
    detected_by_resfinder = ifelse(any(detected_by_resfinder == "yes"), "yes", "no"),
    detected_by_amrfinder = ifelse(any(detected_by_amrfinder == "yes"), "yes", "no"),
    .groups = "drop"
  )

combined.wide <- combined.long %>%
  mutate(present = 1) %>%
  select(gene_name, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0)

combined.wide.supp <- gene.tool.summary %>%
  left_join(combined.wide, by = "gene_name") %>%
  relocate(all_of(intersect(isolate_keep_order, names(.))), .after = detected_by_amrfinder) %>%
  arrange(gene_name)

# ============================================================
# Combined long format (Tab S6) ----
# ============================================================
rgi.long.out <- rgi %>%
  mutate(tool = "RGI", predicted_phenotype = antibiotic) %>%
  select(tool, isolate, mlst_type, gene_name, model_type, drug_class,
         resistance_mechanism, predicted_phenotype,
         sequence_identity, sequence_coverage, confirmed_in_duplicate)

resfinder.long.out <- resfinder.resistance.all %>%
  mutate(tool = "ResFinder", model_type = NA,
         resistance_mechanism = NA, drug_class = NA) %>%
  select(tool, isolate, mlst_type, gene_name, model_type, drug_class,
         resistance_mechanism, predicted_phenotype,
         sequence_identity, sequence_coverage, confirmed_in_duplicate)

amrfinder.long.out <- amrfinder %>%
  mutate(tool = "AMRFinderPlus", model_type = element_subtype,
         resistance_mechanism = NA, predicted_phenotype = NA,
         sequence_identity = NA, sequence_coverage = NA) %>%
  rename(drug_class_amr = drug_class) %>%
  mutate(drug_class = drug_class_amr) %>%
  select(tool, isolate, mlst_type, gene_name, model_type, drug_class,
         resistance_mechanism, predicted_phenotype,
         sequence_identity, sequence_coverage, confirmed_in_duplicate)

combined.long.out <- bind_rows(rgi.long.out, resfinder.long.out, amrfinder.long.out) %>%
  filter(isolate %in% isolate_keep_order) %>%
  arrange(tool, isolate, gene_name)

# ============================================================
# ABRicate / VFDB (Virulence) ----
# ============================================================
abricatedir <- paste0(datadir, "ABRicate/")
files <- list.files(abricatedir, full.names = TRUE)

abricate <- do.call(rbind, lapply(files, function(f) {
  read.delim(f, header = TRUE, stringsAsFactors = FALSE)
}))

abricate$sampleID <- sub("_.*$", "", abricate$X.FILE)
abricate <- abricate[abricate$X.COVERAGE >= 80, ]
abricate <- merge(abricate, samples.isolates, by.x = "sampleID", by.y = "sample.id", all.x = TRUE)
abricate[abricate$sampleID %in% ref.genomes, "isolate"] <- abricate[abricate$sampleID %in% ref.genomes, "sampleID"]
abricate[abricate$sampleID %in% ref.genomes, "Species"] <- "Pseudomonas aeruginosa"

abricate <- abricate %>%
  group_by(isolate, GENE, PRODUCT) %>%
  mutate(confirmed_in_duplicate = ifelse(n_distinct(sampleID) > 1, "yes", "no")) %>%
  ungroup()

abricate <- abricate[!abricate$isolate %in% ref.genomes.remove, ]
abricate <- abricate %>%
  filter(!((confirmed_in_duplicate == "no") & (isolate != "PAO1"))) %>%
  group_by(isolate, GENE, PRODUCT) %>%
  dplyr::slice(1) %>%
  ungroup()

abricate <- abricate[abricate$Species == "Pseudomonas aeruginosa", ]
abricate[abricate$isolate == "PAO1", "mlst_type"] <- "549"

abricate <- abricate %>%
  mutate(
    func.annot.full      = str_remove_all(str_extract(PRODUCT, "\\[[^\\]]+\\]"), "[\\[\\]]"),
    functional_category  = str_replace(func.annot.full, "\\s*\\(VF[0-9]+\\)", ""),
    functional_category  = str_replace(functional_category, "\\s*\\(AI[0-9]+\\)", "")
  ) %>%
  rename(gene_name = GENE, product = PRODUCT) %>%
  select(isolate, mlst_type, gene_name, product, functional_category, confirmed_in_duplicate)

abricate.wide <- abricate %>%
  mutate(present = 1) %>%
  select(gene_name, product, functional_category, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(intersect(isolate_keep_order, names(.))), .after = functional_category) %>%
  arrange(functional_category, gene_name)

abricate.long <- abricate %>%
  filter(isolate %in% isolate_keep_order) %>%
  arrange(functional_category, gene_name, isolate)

# ============================================================
# Write Supplementary_Resistance_Acquired.xlsx ----
# ============================================================
wb_acq <- createWorkbook()

write_sheet(wb_acq,
  "S1 - All tools combined",
  paste0("Combined presence/absence of acquired resistance genes detected by ResFinder, RGI, and/or AMRFinderPlus. ",
         "'detected_by_*' columns indicate whether a gene was detected by each tool across the dataset. ",
         "Isolate columns show presence (1) or absence (0). ",
         "Only genes confirmed in both sequencing duplicates by at least one tool are included. ",
         "Gene names are harmonised across tools where possible. ",
         "Isolates are ordered according to the phylogeny."),
  combined.wide.supp)

write_sheet(wb_acq,
  "S2 - ResFinder",
  paste0("Presence/absence of acquired resistance genes detected by ResFinder (staramr v0.12.1). ",
         "Minimum nucleotide identity 90%, minimum coverage 60%. ",
         "Rows are genes; isolate columns show presence (1) or absence (0). ",
         "Isolates are ordered according to the phylogeny."),
  resfinder.wide)

write_sheet(wb_acq,
  "S3 - RGI acquired",
  paste0("Presence/absence of acquired resistance genes detected by RGI (v6.0.5, CARD v4.0.1), protein homolog model. ",
         "Only 'strict' or 'perfect' hits are included. ",
         "Rows are genes; isolate columns show presence (1) or absence (0). ",
         "Isolates are ordered according to the phylogeny."),
  rgi.acquired.wide)

write_sheet(wb_acq,
  "S4 - AMRFinder acquired",
  paste0("Presence/absence of acquired resistance genes detected by AMRFinderPlus (v3.12.8), ",
         "exact, allele, and BLAST methods (EXACTX/ALLELEX/BLASTX). ",
         "Minimum reference sequence coverage 80%. ",
         "Rows are genes; isolate columns show presence (1) or absence (0). ",
         "Isolates are ordered according to the phylogeny."),
  amrfinder.acquired.wide)

write_sheet(wb_acq,
  "S5 - All tools detail",
  paste0("Long-format table of all acquired resistance gene calls across ResFinder, RGI, and AMRFinderPlus. ",
         "Each row represents one gene detected in one isolate by one tool. ",
         "The 'confirmed_in_duplicate' column indicates whether the gene was detected in both sequencing duplicates."),
  combined.long.out)

saveWorkbook(wb_acq, paste0(resultdir, "Supplementary_Resistance_Acquired.xlsx"), overwrite = TRUE)
message("Saved Supplementary_Resistance_Acquired.xlsx")

# ============================================================
# Write Supplementary_Virulence.xlsx ----
# ============================================================
wb_vir <- createWorkbook()

write_sheet(wb_vir,
  "S9 - Virulence genes",
  paste0("Presence/absence of virulence-associated genes detected by ABRicate (v1.0.1) with the VFDB database. ",
         "Minimum nucleotide identity 90%, minimum coverage 80%. ",
         "Genes are grouped by functional category. ",
         "Rows are genes; isolate columns show presence (1) or absence (0). ",
         "Only genes confirmed in both sequencing duplicates are included (all isolates except PAO1). ",
         "Isolates are ordered according to the phylogeny."),
  abricate.wide)

write_sheet(wb_vir,
  "S10 - Virulence detail",
  paste0("Long-format table of all virulence-associated genes detected by ABRicate (v1.0.1) with the VFDB database. ",
         "Each row represents one gene detected in one isolate, ",
         "with associated product description and functional category."),
  abricate.long)

saveWorkbook(wb_vir, paste0(resultdir, "Supplementary_Virulence.xlsx"), overwrite = TRUE)
message("Saved Supplementary_Virulence.xlsx")

# ============================================================
# Mutation supplementary tables: VCF / GFF pipeline ----
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("vcfR", "rtracklayer")) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}
library(vcfR)
library(rtracklayer)

# GFF: gene ID → name lookup ----
gff.data  <- rtracklayer::import(paste0(datadir, "bacdist/ref.gff"))
gene.info <- as_tibble(gff.data)[c("ID", "gene")]
colnames(gene.info)[2] <- "gene.name"

# VCF: load, extract GT and INFO ----
vcf  <- read.vcfR(paste0(datadir, "bacdist/merged_Pa_raw.vcf.gz"), verbose = FALSE)
gt   <- extract.gt(vcf, element = "GT")
fix  <- getFIX(vcf) %>% as_tibble() %>% mutate(ID = rownames(gt))
info <- getINFO(vcf) %>% as_tibble()

# Parse SnpEff ANN field ----
ann <- bind_cols(fix, info) %>%
  mutate(ANN = str_extract(value, "ANN=[^;]+")) %>%
  filter(!is.na(ANN)) %>%
  mutate(ANN      = str_remove(ANN, "^ANN="),
         ANN_list = str_split(ANN, ",")) %>%
  unnest(ANN_list) %>%
  separate(ANN_list,
           into = c("alt","effect","impact","gene","gene_id","feature",
                    "feature_id","biotype","rank","hgvs_c","hgvs_p",
                    "cdna_pos","cds_pos","aa_pos","distance","errors"),
           sep = "\\|", fill = "right") %>%
  dplyr::select(ID, CHROM, POS, REF, ALT, gene, effect, impact, hgvs_p)

ann.mod.high <- ann[ann$impact %in% c("MODERATE", "HIGH"), ]

# GT to long format ----
gt.long <- gt %>%
  as_tibble(rownames = "ID") %>%
  pivot_longer(-ID, names_to = "sample.id", values_to = "GT") %>%
  mutate(mutated = !is.na(GT) & GT != "0/0")

# Sample metadata: VCF sample IDs are formatted as sampleXX_SYY ----
samples.vcf <- samples.isolates %>%
  mutate(
    sample.id.full = paste(sample, sample.id, sep = "_"),
    mlst_type = if_else(!is.na(mlst_type) & mlst_type == "773!", "773", mlst_type)
  ) %>%
  filter(Species == "Pseudomonas aeruginosa") %>%
  select(isolate, mlst_type, sample.id.full)

# HGVS 3-letter → 1-letter amino acid conversion ----
aa.map.vcf <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",
                His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",
                Ser="S",Thr="T",Trp="W",Tyr="Y",Val="V")

# Build combined mutation table ----
mutations.base <- gt.long %>%
  filter(mutated) %>%
  select(-mutated) %>%
  inner_join(ann.mod.high, by = "ID") %>%
  inner_join(samples.vcf, by = c("sample.id" = "sample.id.full")) %>%
  left_join(gene.info, by = c("gene" = "ID")) %>%
  mutate(gene.name = coalesce(gene.name, gene)) %>%
  group_by(isolate, CHROM, POS, REF, ALT, effect, GT, ID, hgvs_p, gene, gene.name) %>%
  mutate(in.both.replicates = ifelse(n_distinct(sample.id) == 2, "yes", "no")) %>%
  ungroup() %>%
  select(-sample.id) %>%
  distinct()

mutations.base$mutation <- sub("^p\\.", "", mutations.base$hgvs_p)
for (aa in names(aa.map.vcf)) {
  mutations.base$mutation <- gsub(aa, aa.map.vcf[aa], mutations.base$mutation, fixed = TRUE)
}

# Helper: collapse co-occurring mutations (same isolate set within each gene) ----
collapse_cooccurring <- function(muts) {
  muts <- muts %>%
    mutate(
      mut.key = paste(gene.name, POS, REF, ALT, sep = "_"),
      isolate = factor(isolate, levels = isolate_keep_order)
    ) %>%
    group_by(gene, mut.key, impact) %>%
    mutate(isolate.set = paste(sort(unique(as.character(isolate))), collapse = "_")) %>%
    ungroup()

  collapsed <- muts %>%
    distinct(gene, gene.name, mut.key, isolate.set, impact, mutation) %>%
    group_by(gene, gene.name, isolate.set, impact) %>%
    summarise(n.muts = n(), .groups = "drop") %>%
    group_by(gene.name) %>%
    arrange(desc(n.muts), isolate.set, .by_group = TRUE) %>%
    mutate(mut.id = row_number()) %>%
    ungroup() %>%
    mutate(mut.label.old = ifelse(
      n.muts > 1,
      paste0(gene.name, ".mut", mut.id, ".c", n.muts),
      paste0(gene.name, ".mut", mut.id)
    ))

  muts %>%
    left_join(
      collapsed %>% select(gene, gene.name, isolate.set, mut.label.old, n.muts),
      by = c("gene", "gene.name", "isolate.set")
    ) %>%
    group_by(mut.label.old) %>%
    mutate(
      mut.label = if (n.muts[1] == 1) {
        paste0(gene.name, ".", mutation)
      } else {
        aa.pos <- as.numeric(gsub("[A-Z*]", "", mutation))
        paste0(gene.name, ".", min(aa.pos, na.rm = TRUE), "-", max(aa.pos, na.rm = TRUE))
      }
    ) %>%
    ungroup()
}

# Helper: build wide table — one row per mut.label, isolates as 0/1 columns ----
mut_to_wide <- function(muts.collapsed, class_col = NULL) {
  annot <- muts.collapsed %>%
    group_by(gene.name, mut.label, mut.label.old) %>%
    summarise(
      mutations = paste(sort(unique(na.omit(mutation))), collapse = ", "),
      effect    = paste(sort(unique(na.omit(effect))),   collapse = ", "),
      impact    = ifelse(any(impact == "HIGH", na.rm = TRUE), "HIGH", "MODERATE"),
      .groups   = "drop"
    )

  if (!is.null(class_col)) {
    cls <- muts.collapsed %>%
      select(mut.label, all_of(class_col)) %>%
      distinct() %>%
      group_by(mut.label) %>%
      dplyr::slice(1) %>%
      ungroup()
    annot <- annot %>% left_join(cls, by = "mut.label")
  }

  presence <- muts.collapsed %>%
    mutate(present = 1) %>%
    select(mut.label, isolate, present) %>%
    distinct() %>%
    pivot_wider(names_from = isolate, values_from = present, values_fill = 0)

  annot %>%
    left_join(presence, by = "mut.label") %>%
    relocate(all_of(intersect(isolate_keep_order[isolate_keep_order != "PAO1"], names(.))),
             .after = impact)
}

# ============================================================
# Chromosomal resistance mutations ----
# ============================================================

efflux.genes.v  <- c("mexA","mexB","oprM","mexC","mexD","oprJ",
                      "mexE","mexF","oprN","mexX","PA2019","mexY","PA2018")
efflux.reg.v    <- c("mexR","nalC","nalD","nfxB","mexT","mexS","PA2491","mexZ","PA2020")
porin.v         <- c("oprD")
ampc.v          <- c("ampC","ampD","ampE","ampR","dacB","PA3047")
fq.v            <- c("gyrA","gyrB","parC","parE","mexR","nalC","nalD","mexT","nfxB","mexZ")
lps.v           <- c("pmrA","pmrB","phoP","phoQ","arnA","arnB","arnC","arnD","arnE","arnF","pagL")
pbp.v           <- c("pbpA","pbpB","pbpC","ftsI")
resistome.all.v <- unique(c(efflux.genes.v, efflux.reg.v, porin.v, ampc.v, fq.v, lps.v, pbp.v))

mutations.res <- mutations.base %>%
  filter(gene.name %in% resistome.all.v) %>%
  mutate(
    gene.name = case_when(
      gene.name == "PA3047" ~ "dacB",
      gene.name == "PA2491" ~ "mexS",
      gene.name == "PA2019" ~ "mexX",
      gene.name == "PA2018" ~ "mexY",
      gene.name == "PA2020" ~ "mexZ",
      TRUE ~ gene.name
    ),
    res.class = case_when(
      gene.name %in% c("mexA","mexB","oprM","mexC","mexD","oprJ",
                        "mexE","mexF","oprN","mexX","mexY")                       ~ "efflux",
      gene.name %in% c("mexR","nalC","nalD","nfxB","mexT","mexS","mexZ")         ~ "efflux.regulators",
      gene.name == "oprD"                                                         ~ "porins",
      gene.name %in% c("ampC","ampD","ampE","ampR","dacB")                       ~ "ampC",
      gene.name %in% c("gyrA","gyrB","parC","parE")                              ~ "fluoroquinolone",
      gene.name %in% c("pmrA","pmrB","phoP","phoQ",
                        "arnA","arnB","arnC","arnD","arnE","arnF","pagL")         ~ "lps",
      gene.name %in% c("pbpA","pbpB","pbpC","ftsI")                              ~ "pbp",
      TRUE ~ NA_character_
    )
  )

mut.collapsed.res <- collapse_cooccurring(mutations.res)
mut.wide.res <- mut_to_wide(mut.collapsed.res, class_col = "res.class") %>%
  arrange(res.class, gene.name, mut.label)

# ============================================================
# QS gene mutations ----
# ============================================================

qs.genes.v <- c(
  "lasR","lasI","lasA","lasB",
  "rhlR","rhlI","rhlA","rhlB","rhlC",
  paste0("phz", c("A1","B1","C1","D1","E1","F1","G1",
                   "A2","B2","C2","D2","E2","F2","G2")),
  "phzM","phzS", "mvfR",
  paste0("pqs", c("A","B","C","D","E","H","L","R"))
)

mutations.qs <- mutations.base %>%
  filter(gene.name %in% qs.genes.v) %>%
  mutate(
    qs.class = case_when(
      gene.name %in% c("lasR","lasI")                                               ~ "las",
      gene.name %in% c("lasA","lasB")                                               ~ "lasA-B",
      gene.name %in% c("rhlR","rhlI")                                               ~ "rhl",
      gene.name %in% c("rhlA","rhlB","rhlC")                                        ~ "rhlA-C",
      gene.name %in% c(paste0("pqs", c("A","B","C","D","E","H","L","R")), "mvfR")   ~ "pqs",
      gene.name %in% c(paste0("phz", c("A1","B1","C1","D1","E1","F1","G1",
                                         "A2","B2","C2","D2","E2","F2","G2")),
                        "phzM","phzS")                                               ~ "phz",
      TRUE ~ NA_character_
    )
  )

mut.collapsed.qs <- collapse_cooccurring(mutations.qs)
mut.wide.qs <- mut_to_wide(mut.collapsed.qs, class_col = "qs.class") %>%
  arrange(qs.class, gene.name, mut.label)

# ============================================================
# Write Supplementary_Resistance_Mutations.xlsx ----
# ============================================================
wb_rmut <- createWorkbook()

write_sheet(wb_rmut,
  "S6 - RGI mutations",
  paste0("Resistance-associated mutations detected by RGI (v6.0.5, CARD v4.0.1), ",
         "protein variant and protein overexpression models. ",
         "Only 'strict' or 'perfect' hits are included. ",
         "Rows are gene/mutation combinations; isolate columns show presence (1) or absence (0). ",
         "Isolates are ordered according to the phylogeny."),
  rgi.mutations.wide)

write_sheet(wb_rmut,
  "S7 - AMRFinder mutations",
  paste0("Resistance-associated point mutations detected by AMRFinderPlus (v3.12.8), POINTX method. ",
         "Minimum reference sequence coverage 80%. ",
         "Rows are gene/mutation combinations; isolate columns show presence (1) or absence (0). ",
         "Isolates are ordered according to the phylogeny."),
  amrfinder.mutations.wide)

write_sheet(wb_rmut,
  "S8 - Chromosomal mutations",
  paste0(
    "Chromosomal mutations in resistance-associated genes detected by SNP calling relative to PAO1 ",
    "(SnpEff-annotated VCF; MODERATE and HIGH impact variants only). ",
    "Mutations co-occurring in exactly the same set of isolates are collapsed to a single row; ",
    "individual amino acid changes are listed in the 'mutations' column. ",
    "Gene categories (res.class): efflux, efflux.regulators, porins, ampC, fluoroquinolone, lps, pbp. ",
    "Note: genes shared between categories (e.g. mexR in efflux.regulators and fluoroquinolone) ",
    "are assigned to efflux.regulators. ",
    "Isolate columns show presence (1) or absence (0). ",
    "Isolates are ordered according to the phylogeny. ",
    "PAO1 is excluded (reference genome used for SNP calling)."
  ),
  mut.wide.res)

saveWorkbook(wb_rmut, paste0(resultdir, "Supplementary_Resistance_Mutations.xlsx"), overwrite = TRUE)
message("Saved Supplementary_Resistance_Mutations.xlsx")

# ============================================================
# Write Supplementary_QS_Mutations.xlsx ----
# ============================================================
wb_qs <- createWorkbook()

write_sheet(wb_qs,
  "S11 - QS gene mutations",
  paste0(
    "Chromosomal mutations in quorum sensing (QS) genes detected by SNP calling relative to PAO1 ",
    "(SnpEff-annotated VCF; MODERATE and HIGH impact variants only). ",
    "Mutations co-occurring in exactly the same set of isolates are collapsed to a single row; ",
    "individual amino acid changes are listed in the 'mutations' column. ",
    "Gene subsystems (qs.class): las (lasR/lasI), lasA-B, rhl (rhlR/rhlI), rhlA-C, pqs, phz. ",
    "Isolate columns show presence (1) or absence (0). ",
    "Isolates are ordered according to the phylogeny. ",
    "PAO1 is excluded (reference genome used for SNP calling)."
  ),
  mut.wide.qs)

saveWorkbook(wb_qs, paste0(resultdir, "Supplementary_QS_Mutations.xlsx"), overwrite = TRUE)
message("Saved Supplementary_QS_Mutations.xlsx")
