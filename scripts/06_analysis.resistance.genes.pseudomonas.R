# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/06_resistance_genes/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures, recursive = TRUE), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages-------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx", "pheatmap","tibble","ComplexHeatmap","circlize","ggpubr")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

# Send on-screen plots to a null device when run non-interactively (avoids Rplots.pdf)
if (!interactive()) pdf(NULL)

theme_set(theme_minimal(base_family = "Calibri"))


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



##########################################################
# Read and process RGI results ---------------
rgidir <- paste0(datadir,"RGI/")
files <- list.files(rgidir,pattern = "fa(sta)?\\.txt$",full.names = TRUE)

rgi <- lapply(files, function(f) {
  df <- read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  fname <- basename(f)
  sampleID <- if (grepl("_(S[0-9]+)_contigs", fname)) {
    sub(".*_(S[0-9]+)_contigs.*$", "\\1", fname)
  } else {
    sub(".*_([A-Za-z0-9]+)_genome.*$", "\\1", fname)
  }
  df$sampleID <- sampleID
  df
})

rgi <- do.call(rbind, rgi)

rgi = merge(rgi, samples.isolates, by.x="sampleID",by.y="sample.id", all.x=TRUE)

# For reference genomes, set isolate to sampleID and Species to P aeruginosa (for consistency with our isolates)
ref.genomes =c("PAO1","PA7","PA14","LESB5","DK2")
ref.genomes.remove = c("PA7","PA14","LESB5","DK2")
rgi[rgi$sampleID %in% ref.genomes, ]$isolate <- rgi[rgi$sampleID %in% ref.genomes, ]$sampleID
rgi[rgi$sampleID %in% ref.genomes, ]$Species <- "Pseudomonas aeruginosa"
rgi <- rgi[!rgi$sampleID %in% ref.genomes.remove, ]  # remove reference genomes that are not PAO1, to avoid confusion in downstream analyses - can be added back in if wanted to compare with other reference genomes than PAO1

rgi <- rgi[c("sampleID","Contig","Start","Stop","Cut_Off","Best_Hit_ARO","Best_Identities", "Model_type", "SNPs_in_Best_Hit_ARO", "Drug.Class", 
             "Resistance.Mechanism", "AMR.Gene.Family", "Percentage.Length.of.Reference.Sequence","Antibiotic","isolate","MLST_type")]
rgi[rgi$isolate=="PAO1",]$MLST_type <- "549"

# Add column indicating whether gene is present in both duplicates of an isolate
rgi <- rgi %>%
  group_by(isolate, Best_Hit_ARO, Model_type,SNPs_in_Best_Hit_ARO) %>%
  mutate(Present.in.duplicate = ifelse(n_distinct(sampleID) > 1, "yes", "no")) %>%
  ungroup()

# Remove genes not present in both duplicates for selected samples
rgi <- rgi %>%
  filter(!(Present.in.duplicate == "no" & sampleID %in% samples_to_remove))

# Merge duplicates of the same isolate (take union of genes across duplicates, by keeping only first occurance when a gene is predicted in both duplicates)
rgi <- rgi %>%
  group_by(isolate, Best_Hit_ARO, Model_type,SNPs_in_Best_Hit_ARO) %>%
  dplyr::slice(1) %>%        # keep one row per gene per isolate
  ungroup()

# Wide summary with separate columns for yes and no
rgi.presence.in.duplicate.summary <- rgi %>%
  dplyr::count(isolate, sampleID, Present.in.duplicate) %>%
  tidyr::pivot_wider(
    names_from = Present.in.duplicate,
    values_from = n,
    values_fill = 0
  )


rgi$gene.name <- rgi$Best_Hit_ARO
# Add 'bla' prefix to common beta-lactamases
bla.pattern <- "^(OXA|TEM|CTX-M|NDM|IMP|PDC|CMY|SHV|KPC|SRT|PER)"
rgi <- rgi %>%
  mutate(gene.name = ifelse(grepl(bla.pattern, Best_Hit_ARO), paste0("bla", Best_Hit_ARO), gene.name))

# Capitalize first letter
rgi$gene.name = gsub("Pseudomonas aeruginosa ","",rgi$gene.name)
rgi$gene.name <- paste0(tolower(substr(rgi$gene.name, 1, 1)), substr(rgi$gene.name, 2, nchar(rgi$gene.name)))

# Simplify descriptive entries
rgi <- rgi %>%
  mutate(gene.name = case_when(
    grepl("gyrA", Best_Hit_ARO) ~ "gyrA",
    grepl("EF-Tu", Best_Hit_ARO) ~ "EF-Tu",
    TRUE ~ gene.name))

rgi <- rgi[c("isolate","MLST_type","gene.name","Best_Hit_ARO","Model_type", "SNPs_in_Best_Hit_ARO", "Drug.Class", 
             "Resistance.Mechanism" ,"Antibiotic","AMR.Gene.Family","Cut_Off","Best_Identities", "Percentage.Length.of.Reference.Sequence","Present.in.duplicate")]

rgi.full = rgi

rgi.in.all.isolates <- rgi.full %>%
  group_by(Best_Hit_ARO, Model_type, SNPs_in_Best_Hit_ARO) %>%
  filter(n_distinct(isolate) == 14 ) %>%
  ungroup()

rgi <- rgi %>%
  group_by(Best_Hit_ARO, Model_type, SNPs_in_Best_Hit_ARO) %>%
  filter(n_distinct(isolate) < 14) %>%
  ungroup()

rgi <- rgi[!rgi$gene.name %in% c("escherichia coli soxS with mutation conferring antibiotic resistance"), ]  # remove empty calls
rgi <- rgi[which(rgi$isolate %in% isolate_keep_order), ]

rgi.presence <- rgi[rgi$Model_type == "protein homolog model",]

rgi.variant <- rgi[rgi$Model_type == "protein variant model",]
rgi.overexpression <- rgi[rgi$Model_type == "protein overexpression model",]

rgi.duplicated = rgi[rgi$Present.in.duplicate =="yes" | rgi$isolate=="PAO1",]

rgi.wide <- rgi %>%
  mutate(present = 1) %>%
  select(gene.name,Model_type,SNPs_in_Best_Hit_ARO, Drug.Class, Resistance.Mechanism, Antibiotic, AMR.Gene.Family, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = AMR.Gene.Family) %>%
  arrange(Model_type, Drug.Class, gene.name)

rgi.wide.found.in.all <- rgi.in.all.isolates %>%
  mutate(present = 1) %>%
  select(gene.name,Model_type,SNPs_in_Best_Hit_ARO, Drug.Class, Resistance.Mechanism, Antibiotic, AMR.Gene.Family, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = AMR.Gene.Family) %>%
  arrange(Model_type, Drug.Class, gene.name)

#write.xlsx(rgi.wide, "results/rgi_results.xlsx",SheetName = "RGI not in all isolates", append = TRUE)
# Create workbook and write two sheets
wb <- createWorkbook()
addWorksheet(wb, "RGI aquired genes and mutations")
writeData(wb, "RGI aquired genes and mutations", rgi.wide)

addWorksheet(wb, "RGI genes in all isolates")
writeData(wb, "RGI genes in all isolates", rgi.wide.found.in.all)

addWorksheet(wb, "RGI long format")
writeData(wb, "RGI long format", rgi)

saveWorkbook(wb, "results/rgi_results.xlsx", overwrite = TRUE)


##########################################################
# Read and process resfinder results  -------------------

# Read and process resfinder result file on isolates
resfinder <- read.csv(paste0(datadir,"Resfinder_plasmidfinder.tsv.tabular"), header = TRUE, sep="\t")
# Change sample id
resfinder$sample <- sub("_.*$", "", resfinder$sampleID)
# Merge with metadata
resfinder <- merge(resfinder, samples.isolates,by.x="sample",by.y="sample.id", all.x=TRUE)
# Select for Pseudomonas 
resfinder <- resfinder[which(resfinder$Species=="Pseudomonas aeruginosa"), ]
resfinder <- resfinder[c("sample","isolate","data","data_type","CGE_predicted_phenotype", "X.identity","X.overlap")]
colnames(resfinder) <- c("sample","isolate","gene.name","data_type","CGE_predicted_phenotype","Identity.RF","Overlap.RF")

# Read and process resfinder file on reference genomes (keep only PAO1) 
resfinder.ref <- read.csv(paste0(datadir,"Staramr_reference_genomes.tabular"), header = TRUE, sep="\t")
# Change sample id
resfinder.ref$sample <- sub("_.*$", "", resfinder.ref$Isolate.ID)
resfinder.ref$isolate <- resfinder.ref$sample
# remove reference genomes that are not PAO1
resfinder.ref <- resfinder.ref[!resfinder.ref$isolate %in% ref.genomes.remove, ]  
# Select specific columns 
resfinder.ref <- resfinder.ref[c("sample","isolate","Data","Data.Type","CGE.Predicted.Phenotype", "X.Identity","X.Overlap")]
colnames(resfinder.ref) <- c("sample","isolate","gene.name","data_type","CGE_predicted_phenotype","Identity.RF","Overlap.RF")

# Combine isolate data with PAO1 data
resfinder <- bind_rows(resfinder, resfinder.ref)

resfinder[resfinder$sample=="S19" & resfinder$gene.name=="aac(6')-Ib",]$CGE_predicted_phenotype = "Tobramycin, Dibekacin, Amikacin, Sisomicin, Netilmicin, Fluoroquinolone, Ciprofloxacin"
resfinder[resfinder$sample=="S19" & resfinder$gene.name=="aac(6')-Ib",]$gene.name = "aac(6')-Ib-cr"


# Add column indicating whether gene is present in both duplicates of an isolate
resfinder <- resfinder %>%
  group_by(isolate, gene.name, data_type, CGE_predicted_phenotype) %>%
  mutate(Present.in.duplicate = ifelse(n_distinct(sample) > 1, "yes", "no")) %>%
  ungroup()

#resfinder.resistance[resfinder.resistance$sample=="S19" & resfinder.resistance$gene.name=="aac(6')-Ib",]$gene.name = "aac(6')-Ib-cr"

# Wide summary with separate columns for yes and no
resfinder.presence.in.duplicate.summary <- resfinder %>%
  dplyr::count(isolate, sample, Present.in.duplicate) %>%
  tidyr::pivot_wider(
    names_from = Present.in.duplicate,
    values_from = n,
    values_fill = 0
  )

# Remove genes not present in both duplicates for selected samples
resfinder <- resfinder %>%
  filter(!(Present.in.duplicate == "no" & sample %in% samples_to_remove))

# Take union of genes across duplicates, by keeping only first occurance when a gene is found in both duplicates
resfinder <- resfinder %>%
  group_by(isolate, gene.name, CGE_predicted_phenotype) %>%
  dplyr::slice(1) %>%        # keep one row per gene per isolate
  ungroup()

# Select for resistance 
resfinder.resistance <- resfinder[resfinder$data_type == "Resistance", ]
# Select for resistance 
resfinder.plasmid <- resfinder[resfinder$data_type == "Plasmid", ]

resfinder_to_rgi <- c(
  #"aac(6')-Ib"   = "aAC(6')-Ib10",
  "aph(3')-VIa"  = "aPH(3')-VIa",
  "aph(3'')-Ib"  = "aPH(3'')-Ib",
  "aph(6)-Id"    = "aPH(6)-Id",
  "mph(E)"       = "mphE",
  "mph(F)"       = "mphF",
  "msr(E)"       = "msrE",
  "blaIMP-1b"    = "blaIMP-1",
  "blaTEM-1A"    = "blaTEM-1",
  "blaTEM-1B"    = "blaTEM-1",
  "ARR-3"        = "arr-3",
  "aadA1"        = "aadA",
  "aadA6"        = "aadA"
)

# Add a column in ResFinder that matches RGI names
resfinder.resistance$gene.name.rgi <- ifelse(
  resfinder.resistance$gene.name %in% names(resfinder_to_rgi),
  resfinder_to_rgi[resfinder.resistance$gene.name],
  resfinder.resistance$gene.name
)


resfinder.resistance.full = resfinder.resistance

resfinder.resistance <- resfinder.resistance %>%
  group_by(gene.name) %>%
  filter(n_distinct(isolate) < 14) %>%
  ungroup()

resfinder.resistance.duplicated = resfinder.resistance[resfinder.resistance$Present.in.duplicate =="yes" | resfinder.resistance$isolate=="PAO1",]

resfinder.in.all.isolates <- resfinder.resistance.full %>%
  group_by(gene.name) %>%
  filter(n_distinct(isolate) == 14 ) %>%
  ungroup()

resfinder.wide <- resfinder.resistance.duplicated %>%
  mutate(present = 1) %>%
  select(gene.name,CGE_predicted_phenotype,isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = CGE_predicted_phenotype) %>%
  arrange(CGE_predicted_phenotype, gene.name)

resfinder.wide.found.in.all <- resfinder.in.all.isolates %>%
  mutate(present = 1) %>%
  select(gene.name,CGE_predicted_phenotype,isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = CGE_predicted_phenotype) %>%
  arrange(CGE_predicted_phenotype, gene.name)

#write.xlsx(resfinder.wide, "results/resfinder_results.xlsx")

# Create workbook and write two sheets
wb <- createWorkbook()
addWorksheet(wb, "Resfinder aquired genes")
writeData(wb, "Resfinder aquired genes", resfinder.wide)

addWorksheet(wb, "Resfinder genes in all isolates")
writeData(wb, "Resfinder genes in all isolates", resfinder.wide.found.in.all)

addWorksheet(wb, "Resfinder long format")
writeData(wb, "Resfinder long format", resfinder)

saveWorkbook(wb, "results/resfinder_results.xlsx", overwrite = TRUE)



##########################################################
# Read and process AMRFInderPlus result file on isolates  -------------------

amrfinderdir <- paste0(datadir,"AMRFinderPlus/")
files <- list.files(amrfinderdir, full.names = TRUE)

amrfinder <- lapply(files, function(f) {
  df <- read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  sampleID <- sub("_contigs.*$", "", basename(f))
  df$sampleID <- sampleID
  # Ensure Hierarchy.node column exists
  if(!"Hierarchy.node" %in% colnames(df)) {
    df$Hierarchy.node <- NA
  }
  df
})

amrfinder <- do.call(rbind, amrfinder)

# Filter for coverage >= 80%
amrfinder = amrfinder[amrfinder$X..Coverage.of.reference.sequence >=80, ]

# Merge with metadata
amrfinder = merge(amrfinder, samples.isolates, by.x="sampleID",by.y="sample.id", all.x=TRUE)

amrfinder[amrfinder$sampleID == "PAO1", ]$Species <- "Pseudomonas aeruginosa"
amrfinder[amrfinder$sampleID == "PAO1", ]$isolate <- "PAO1"
amrfinder[amrfinder$sampleID == "PAO1",]$MLST_type <- "549"

# Select for Pseudomonas - AMRfinderPLUS results for now are already only for pseudomonas
amrfinder <- amrfinder[which(amrfinder$Species=="Pseudomonas aeruginosa"), ]

# Add column indicating whether gene is present in both duplicates of an isolate
amrfinder <- amrfinder %>%
  group_by(isolate, Gene.symbol, Element.type, Element.subtype, Method, Class, Subclass) %>%
  mutate(Present.in.duplicate = ifelse(n_distinct(sampleID) > 1, "yes", "no")) %>%
  ungroup()

# Remove genes not present in both duplicates for selected samples
amrfinder <- amrfinder %>%
  filter(!(Present.in.duplicate == "no" & sampleID %in% samples_to_remove))

# Merge duplicates of the same isolate (take union of genes across duplicates, by keeping only first occurance when a gene is predicted in both duplicates)
amrfinder <- amrfinder %>%
  group_by(isolate, Gene.symbol, Element.type, Element.subtype, Class, Subclass) %>%
  dplyr::slice(1) %>%        # keep one row per gene per isolate
  ungroup()

amrfinder_to_rgi <- c(
  "aph(3')-VIa"  = "aPH(3')-VIa",
  "aph(3'')-Ib"  = "aPH(3'')-Ib",
  "aph(6)-Id"    = "aPH(6)-Id",
  "mph(E)"       = "mphE",
  "mph(F)"       = "mphF",
  "msr(E)"       = "msrE"
)

# Add a column in amrfinder that matches RGI names
amrfinder$gene.name.rgi <- ifelse(
  amrfinder$Gene.symbol %in% names(resfinder_to_rgi),
  amrfinder_to_rgi[amrfinder$Gene.symbol],
  amrfinder$Gene.symbol
)


# Select only the necessary columns, for simplicity 
amrfinder <- amrfinder[c("isolate","MLST_type", "Gene.symbol","gene.name.rgi", "Element.type","Element.subtype","Method",
              "Sequence.name","Hierarchy.node","Class","Subclass","Present.in.duplicate")]

amrfinder.full = amrfinder

#Remove genes found in all
amrfinder <- amrfinder %>%
  group_by( Gene.symbol, Element.type, Element.subtype, Class, Subclass) %>%
  filter(n_distinct(isolate) < 14) %>%
  ungroup()

amrfinder.in.all.isolates <- amrfinder.full %>%
  group_by( Gene.symbol, Element.type, Element.subtype, Class, Subclass) %>%
  filter(n_distinct(isolate) == 14 ) %>%
  ungroup()

amrfinder.presence <- amrfinder %>% filter(Method %in% c("EXACTX","ALLELEX","BLASTX"))
amrfinder.mutations <- amrfinder %>% filter(Method %in% c("POINTX"))

amrfinder.duplicated = amrfinder[amrfinder$Present.in.duplicate =="yes" | amrfinder$isolate=="PAO1",]


# Reshape mutation table to wide (1/0 per isolate) and export as xlsx
amrfinder.wide <- amrfinder.duplicated %>%
  mutate(present = 1) %>%
  select(Gene.symbol,Sequence.name, Element.subtype, Class, Subclass, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = Subclass) %>%
  arrange(Element.subtype, Class, Subclass, Gene.symbol)

amrfinder.wide.found.in.all <- amrfinder.in.all.isolates %>%
  mutate(present = 1) %>%
  select(Gene.symbol,Sequence.name, Element.subtype, Class,Subclass, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = Subclass) %>%
  arrange(Element.subtype, Class, Subclass, Gene.symbol)

# Create workbook and write two sheets
wb <- createWorkbook()
addWorksheet(wb, "AMRFinder genes and mutations")
writeData(wb, "AMRFinder genes and mutations", amrfinder.wide)
addWorksheet(wb, "AMRFinder genes in all isolates")
writeData(wb, "AMRFinder genes in all isolates", amrfinder.wide.found.in.all)
addWorksheet(wb, "AMRFinder long format")
writeData(wb, "AMRFinder long format", amrfinder)

saveWorkbook(wb, "results/amrfinder_results.xlsx", overwrite = TRUE)


##########################################################
# Merge/compare results between the 3 tools  -----------------------------------------

rgi.presence.less <- rgi.presence[c("isolate","gene.name","Present.in.duplicate")] 
resfinder.resistance.less <- resfinder.resistance[c("isolate","gene.name","gene.name.rgi","Present.in.duplicate" )]
amrfinder.presence.less <- amrfinder.presence[c("isolate","Gene.symbol","gene.name.rgi","Present.in.duplicate")]
colnames(rgi.presence.less)[3] = "duplicate.rgi"
colnames(resfinder.resistance.less)[4] = "duplicate.resfinder"
colnames(resfinder.resistance.less)[2] = "gene.resfinder"
colnames(amrfinder.presence.less)[2] = "gene.amrfinder"
colnames(amrfinder.presence.less)[4] = "duplicate.amrfinder"

rgi.resfinder <- merge(rgi.presence.less, resfinder.resistance.less, by.x=c("isolate","gene.name"), by.y=c("isolate","gene.name.rgi"), all=TRUE)
rgi.resfinder.amrfinder <- merge(rgi.resfinder, amrfinder.presence.less,  by.x=c("isolate","gene.name"), by.y=c("isolate","gene.name.rgi"), all=TRUE)

# Keep rows with at least one duplicate == "yes", retain all PAO1 rows
combined.filtered.duplicate <- rgi.resfinder.amrfinder %>%
  filter(
    isolate == "PAO1" |
      duplicate.rgi == "yes" |
      duplicate.resfinder == "yes" |
      duplicate.amrfinder == "yes"
  )

# Rows with no duplicate in any database (excluding PAO1)
combined.no.duplicate <- rgi.resfinder.amrfinder %>%
  filter(
    isolate != "PAO1" &
      (is.na(duplicate.rgi) | duplicate.rgi != "yes") &
      (is.na(duplicate.resfinder) | duplicate.resfinder != "yes") &
      (is.na(duplicate.amrfinder) | duplicate.amrfinder != "yes")
  )

gene.names.rgi = unique(na.omit(combined.filtered.duplicate$gene.name))
gene.names.resfinder = unique(na.omit(combined.filtered.duplicate$gene.resfinder))
gene.names.amrfinder = unique(na.omit(combined.filtered.duplicate$gene.amrfinder))

all.gene.names = unique(c(gene.names.rgi,gene.names.resfinder,gene.names.amrfinder))
all.gene.names = sort(all.gene.names)

#rgi.not.in.resfinder = setdiff(rgi.presence$gene.name, resfinder.resistance$gene.name.rgi)
#resfinder.not.in.rgi = setdiff(resfinder.resistance$gene.name.rgi, rgi.presence$gene.name)

# Detect genes to remove from initial resfinder results, due to not being present in duplicate
resfinder.genes.to.remove <- merge(resfinder.resistance,combined.no.duplicate[!is.na(combined.no.duplicate$gene.resfinder),], 
                          by.x= c("isolate","gene.name"),by.y=c("isolate","gene.resfinder"))
resfinder.genes.to.remove <- resfinder.genes.to.remove[c(1,2)]
resfinder.genes.to.remove <- resfinder.genes.to.remove[order(resfinder.genes.to.remove$gene.name, resfinder.genes.to.remove$gene.name),]

write.xlsx(resfinder.genes.to.remove, "results/resfinder_gene_isolate_pairs_to_remove.xlsx")


# Remove isolate–gene.name pairs listed for removal 
resfinder.resistance.filtered <- resfinder.resistance %>%
  anti_join(resfinder.genes.to.remove, by = c("isolate", "gene.name"))

a=merge(combined.filtered.duplicate,resfinder.resistance.filtered[c(2,9,3,5)],by.x=c("isolate","gene.name"),by.y=c("isolate","gene.name.rgi"),all.x=TRUE)

# From combined filtered list, keep the genes not found in resfinder, in order to add to resfinder heatmap

combined.filtered.duplicate.rgi.amrfinder <- combined.filtered.duplicate[is.na(combined.filtered.duplicate$gene.resfinder),]
pairs.unique.amrfinder <- combined.filtered.duplicate.rgi.amrfinder[is.na(combined.filtered.duplicate.rgi.amrfinder$gene.name),]
pairs.unique.rgi <- combined.filtered.duplicate.rgi.amrfinder[is.na(combined.filtered.duplicate.rgi.amrfinder$gene.amrfinder),]

# Keep rows in rgi that are present in combined.filtered.duplicate.rgi.amrfinder
rgi.data.extra <- rgi %>%
  semi_join(combined.filtered.duplicate.rgi.amrfinder,
            by = c("isolate", "gene.name"))

iso.order <- isolate_keep_order
plot.data <- rgi.data.extra
plot.data <- plot.data %>% mutate(isolate = factor(isolate, levels = iso.order))
plot.data$gene.name <- factor(plot.data$gene.name,
                              levels = rev(sort(unique(plot.data$gene.name))))

# Plot heatmap ------
ggplot(plot.data, aes(x = isolate, y = gene.name)) +
  geom_tile(fill = "#95AAD3",color = "white", linewidth=0.4) +
  labs(x = "Isolate", y = "Gene") +
  theme_minimal() +
  scale_y_discrete(position = "right")+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

ggsave(paste0(figures,"resistance_genes_rgi_minus_resfinder.pdf"), width = 4, height = 6)

rgi.extra.resfinder.wide <- rgi.data.extra %>%
  mutate(present = 1) %>%
  select(gene.name,Drug.Class, Resistance.Mechanism, Antibiotic, AMR.Gene.Family, isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order), .after = AMR.Gene.Family) %>%
  arrange(Drug.Class, gene.name)

write.xlsx(rgi.extra.resfinder.wide, "results/extra_resistance_genes_not_in_resfinder.xlsx")

# Get full table of presence/absence for combined rgi/resfinder ------

full.table.presence.res = rbind(resfinder.wide[c(1,3:16)], rgi.extra.resfinder.wide[c(1,6:19)])

# Heatmap of acquired resistance genes, grouped by class (basis for Fig. 3A) -------
# Genes confirmed in both sequencing duplicates (ResFinder genes plus genes detected only
# by RGI/AMRFinderPlus); genes present in all isolates are not shown.
# The gene subset shown in the published Fig. 3A was selected manually from this heatmap.
presence.col <- "#95AAD3"  # blue used for gene presence in all manuscript figures

gene.class <- dplyr::case_when(
  grepl("^bla(IMP|NDM|VIM|KPC|GES|SPM)", full.table.presence.res$gene.name) ~ "Carbapenemases",
  grepl("^bla", full.table.presence.res$gene.name)                         ~ "Other beta-lactamases",
  grepl("^(aac|aad|aph|ant|rmt)", full.table.presence.res$gene.name, ignore.case = TRUE) ~ "Aminoglycosides",
  grepl("^(qnr|crpP)", full.table.presence.res$gene.name)                  ~ "Fluoroquinolones",
  TRUE                                                                      ~ "Other"
)
gene.class <- factor(gene.class, levels = c("Carbapenemases", "Other beta-lactamases",
                                            "Aminoglycosides", "Fluoroquinolones", "Other"))

res.mat <- as.matrix(full.table.presence.res[, isolate_keep_order])
rownames(res.mat) <- full.table.presence.res$gene.name
res.mat <- res.mat[order(gene.class, rownames(res.mat)), ]
gene.class <- sort(gene.class)

# Isolates as rows (PAO1 on top), genes as columns, as in Fig. 3A
ht.res <- Heatmap(
  t(res.mat)[rev(isolate_keep_order), ],
  name = "Presence",
  col = c("0" = "white", "1" = presence.col),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_split = gene.class,
  column_title_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 9),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "grey60", lwd = 0.3),
  heatmap_legend_param = list(at = c(0, 1), labels = c("Not present", "Present"))
)

pdf(paste0(figures, "heatmap_resistance_genes_by_class.pdf"), width = 16, height = 5)
draw(ht.res)
dev.off()
# Sum presence (1) per isolate column
# Column sums excluding gene.name
total.genes.per.isolate.res <- data.frame(
  isolate = names(colSums(full.table.presence.res[,-1], na.rm = TRUE)),
  total.genes.resistance = colSums(full.table.presence.res[,-1], na.rm = TRUE)
)

virulence <- read.xlsx("results/abricate_virulence_results.xlsx")
full.table.presence.vir = virulence[c(1,4:17)]

total.genes.per.isolate.vir <- data.frame(
  isolate = names(colSums(full.table.presence.vir[,-1], na.rm = TRUE)),
  total.genes.virulence = colSums(full.table.presence.vir[,-1], na.rm = TRUE)
)

total.genes.res.vir = merge(total.genes.per.isolate.res,total.genes.per.isolate.vir,by="isolate")
total.genes.res.vir 

total.genes.res.vir <- total.genes.res.vir |>
  dplyr::mutate(isolate = factor(isolate, levels = isolate_keep_order)) |>
  dplyr::arrange(isolate)


# Pivot to long format for faceting
plot.df <- total.genes.res.vir |>
  pivot_longer(cols = c(total.genes.resistance, total.genes.virulence),
               names_to = "type",
               values_to = "count")

plot.df <- plot.df |>
  dplyr::mutate(type = dplyr::recode(
    type,
    total.genes.resistance = "Resistance genes",
    total.genes.virulence  = "Virulence genes"
  ))

# Keep custom isolate order if defined
#plot.df$isolate <- factor(plot.df$isolate, levels = isolate.keep.order)

# Faceted barplot
p<-ggplot(plot.df, aes(x = isolate, y = count)) +
  geom_bar(stat = "identity", width=0.8, fill = "#95AAD3") +
  facet_wrap(~ type, ncol = 1, scales = "free_y") +
  labs(x = "Isolate", y = "Number of genes") +
  theme_minimal()+
  theme(
    axis.text = element_text(size = 12),       # x/y axis tick labels
    axis.title = element_text(size = 14),      # x/y axis titles
    strip.text = element_text(size = 13)       # facet titles
  )


table.df <- as.data.frame(t(total.genes.res.vir))
colnames(table.df)=table.df[1,]
table.df = table.df[-1,]

rownames(table.df) = c("Resistance genes", "Virulence genes")

table.plot <- ggtexttable(
  table.df,
  theme = ttheme("light")
)

ggarrange(p, table.plot, ncol = 1, heights = c(2,1))
ggsave(paste0(figures,"barplot.table.resistance.virulence.genes.per.isolate.pdf"), width = 6.5, height = 5)

