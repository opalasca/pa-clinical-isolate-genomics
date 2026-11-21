# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/08_quorum_sensing/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures, recursive = TRUE), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages-------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx", "pheatmap","tibble",
                      "ComplexHeatmap","circlize", "vcfR","purrr","rtracklayer","patchwork")

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

# Prepare metadata for intersecting with mutation data -------
samples = samples.isolates
samples$sample.id.full = paste(samples$sample, samples$sample.id,sep="_")
samples = samples[samples$Species =="Pseudomonas aeruginosa",]
samples = samples[c("isolate","MLST_type", "sample.id.full" )]
samples[samples$MLST_type=="773!", "MLST_type"] <- "773"

# Remove samples with too long genome size (therefore suspected contamination)
samples_to_remove = c("S68", "S69", "S39","S54")

# Define order based on phylogenetic tree (tip_order)
isolate_keep_order <-  c("32", "25", "22", "4", "3", "18", "5", "27", "11", "33", "19", "9", "26", "PAO1")

# Read xlsx -------------------

pheno <- read.xlsx(paste0(datadir,"Phenotypic_data.xlsx"))
sheets <- getSheetNames(paste0(datadir, "Phenotypic_data.xlsx"))
pheno_list <- lapply(sheets, function(sheet) {
  read.xlsx(paste0(datadir, "Phenotypic_data.xlsx"), sheet = sheet)
})
names(pheno_list) <- sheets

las = as.data.frame(pheno_list[["P. a Quorum sensing(Las)"]])
rhl = as.data.frame(pheno_list[["P. a Quorum sensing(Rhl)"]])
pyocyanin = as.data.frame(pheno_list[["P. a Pyocyanin production"]])
rhamnolipid = as.data.frame(pheno_list[["P. a Rhamnolipid production"]])

# Create pheno data, adding phenotypic data as the different sheets are processed -----------
# Rerun for resetting

pheno.data = samples.isolates %>% filter(Species %in% c("Pseudomonas aeruginosa")) 
pheno.data = pheno.data[c("isolate","MLST_type") ]
pheno.data[pheno.data$MLST_type=="773!", "MLST_type"] <- 773
pheno.data <- pheno.data %>% add_row(isolate = "PAO1",MLST_type = "549")
pheno.data <- unique(pheno.data)

# LasR -----------
data = las; name = "las" 

isolates <- c("PAO1", "QSI", "3", "4", "5", "9", "11","18", "19", "22", "25", "26", "27", "32", "33")
isolate_vec <- rep(isolates, each = 3)
replicate_vec <- rep(1:3, times = length(isolates))
n_time <- nrow(data)

data.long <- data %>%
  pivot_longer(-1, values_to = "value") %>%
  mutate(
    isolate   = rep(isolate_vec, times = nrow(data)),
    replicate = rep(replicate_vec, times = nrow(data))
  ) %>%
  transmute(time = .[[1]], isolate, replicate, value)

data.median <- data.long %>%
  group_by(time, isolate) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

data.median.log <- data.long %>%
  group_by(time, isolate) %>%
  summarise(value = median(log2(value), na.rm = TRUE), .groups = "drop")

ggplot(data.median, aes(x = time, y = value, color = isolate)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(y = "QS reporter value", x = "Time (h)")
ggsave(paste0(figures, name, "_QS_reporter_over_time.pdf"), width = 7, height = 5)

ggplot(data.median.log, aes(x = time, y = value, color = isolate)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(y = "QS reporter value", x = "Time (h)")
ggsave(paste0(figures, name, "_QS_reporter_over_time_log2.pdf"), width = 7, height = 5)

selected.time.point <- data.median %>%
  filter(time == 9) %>%
  group_by(isolate) 

selected.time.point$isolate <- factor(selected.time.point$isolate, levels = isolates)
selected.time.point$value <- as.numeric(selected.time.point$value)

ggplot(selected.time.point, aes(isolate, value)) +
  geom_col() +
  theme_minimal() +
  labs(y = paste(name,"reporter value, hour 9"), x = "isolate")

ggsave(paste0(figures, name, "_reporter_hour9.pdf"), width = 7, height = 3)

pheno.data <- pheno.data %>%
  left_join(
    selected.time.point %>% dplyr::select(isolate, value),
    by = "isolate"
  ) %>%
  dplyr::rename(!!name := value)

  
# Rhl -----------
data = rhl; name = "rhl"

isolates <- c("PAO1", "QSI", "3", "4", "5", "9", "11","18", "19", "22", "25", "26", "27", "32", "33")
isolate_vec <- rep(isolates, each = 3)
replicate_vec <- rep(1:3, times = length(isolates))
n_time <- nrow(data)

data.long <- data %>%
  pivot_longer(-1, values_to = "value") %>%
  mutate(
    isolate   = rep(isolate_vec, times = nrow(data)),
    replicate = rep(replicate_vec, times = nrow(data))
  ) %>%
  transmute(time = .[[1]], isolate, replicate, value)

data.median <- data.long %>%
  group_by(time, isolate) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

data.median.log <- data.long %>%
  group_by(time, isolate) %>%
  summarise(value = median(log2(value), na.rm = TRUE), .groups = "drop")

ggplot(data.median, aes(x = time, y = value, color = isolate)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(y = "QS reporter value", x = "Time (h)")
ggsave(paste0(figures, name, "_QS_reporter_over_time.pdf"), width = 7, height = 5)

ggplot(data.median.log, aes(x = time, y = value, color = isolate)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(y = "QS reporter value", x = "Time (h)")
ggsave(paste0(figures, name, "_QS_reporter_over_time_log2.pdf"), width = 7, height = 5)

selected.time.point <- data.median %>%
  filter(time == 9) %>%
  group_by(isolate) 

selected.time.point$isolate <- factor(selected.time.point$isolate, levels = isolates)
selected.time.point$value <- as.numeric(selected.time.point$value)

ggplot(selected.time.point, aes(isolate, value)) +
  geom_col() +
  theme_minimal() +
  labs(y = paste(name,"reporter value, hour 9"), x = "isolate")

ggsave(paste0(figures, name, "_reporter_hour9.pdf"), width = 7, height = 3)

pheno.data <- pheno.data %>%
  left_join(
    selected.time.point %>% dplyr::select(isolate, value),
    by = "isolate"
  ) %>%
  dplyr::rename(!!name := value)

# Pyocyanin -----------
data = pyocyanin; name = "pyocyanin"
data = data[-1,-1]

isolates <- c("PAO1", "QSI", "3", "4", "5", "9", "11","18", "19", "22", "25", "26", "27", "32", "33")
isolate_vec <- rep(isolates, each = 3)
replicate_vec <- rep(1:3, times = length(isolates))

data.long <- data %>%
  pivot_longer(everything(), values_to = "value") %>%
  mutate(
    isolate   = isolate_vec,
    replicate = replicate_vec
  ) %>%
  dplyr::select(isolate, replicate, value)

data.long$value = as.numeric(data.long$value)

data.median <- data.long %>%
  group_by(isolate) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

data.median$isolate <- factor(data.median$isolate, levels = isolates)

ggplot(data.median, aes(isolate, value)) +
  geom_col() +
  theme_minimal() +
  labs(y = paste(name,"production"), x = "isolate")

ggsave(paste0(figures, name, "_production.pdf"), width = 7, height = 3)

pheno.data <- pheno.data %>%
  left_join(
    data.median %>% dplyr::select(isolate, value),
    by = "isolate"
  ) %>%
  dplyr::rename(!!name := value)

# Rhamnolipid -----------

data = rhamnolipid; name = "rhamnolipid"
data = data[-c(1,47)]

isolates <- c("PAO1", "QSI", "3", "4", "5", "9", "11","18", "19", "22", "25", "26", "27", "32", "33")
isolate_vec <- rep(isolates, each = 3)
replicate_vec <- rep(1:3, times = length(isolates))

data.long <- data %>%
  pivot_longer(everything(), values_to = "value") %>%
  mutate(
    isolate   = isolate_vec,
    replicate = replicate_vec
  ) %>%
  dplyr::select(isolate, replicate, value)

data.long$value = as.numeric(data.long$value)

data.median <- data.long %>%
  group_by(isolate) %>%
  summarise(value = median(value, na.rm = TRUE), .groups = "drop")

data.median$isolate <- factor(data.median$isolate, levels = isolates)

ggplot(data.median, aes(isolate, value)) +
  geom_col() +
  theme_minimal() +
  labs(y = paste(name,"production"), x = "isolate")

ggsave(paste0(figures, name, "_production.pdf"), width = 7, height = 3)

pheno.data <- pheno.data %>%
  left_join(
    data.median %>% dplyr::select(isolate, value),
    by = "isolate"
  ) %>%
  dplyr::rename(!!name := value)

# parse GFF, to be able to match gene ID from vcf to gene name ------
gff <- import(paste0(datadir,"bacdist/ref.gff"))
gene.info <- as_tibble(gff) 
gene.info <- gene.info[c("ID","gene")]
colnames(gene.info)[2] = "gene.name"

# Read vcf file ---------

#vcf <- read.vcfR(paste0(datadir,"bacdist/Pa_final.vcf"), verbose = FALSE)
vcf <- read.vcfR(paste0(datadir,"bacdist/merged_Pa_raw.vcf.gz"), verbose = FALSE)

fix <- getFIX(vcf) %>% as_tibble()
info <- getINFO(vcf) %>% as_tibble()
gt <- extract.gt(vcf, element = "GT") 
fix$ID = rownames(gt)
info.keep = info
info.keep$ID=rownames(gt)

ann <- bind_cols(fix, info) %>%
  mutate(ANN = str_extract(value, "ANN=[^;]+")) %>%
  filter(!is.na(ANN)) %>%
  mutate(ANN = str_remove(ANN, "^ANN="),           # remove ANN= prefix
         ANN_list = str_split(ANN, ",")) %>%      # split multiple effects
  unnest(ANN_list) %>%                             # create one row per effect
  separate(ANN_list,
           into = c("alt","effect","impact","gene","gene_id","feature",
                    "feature_id","biotype","rank","hgvs_c","hgvs_p",
                    "cdna_pos","cds_pos","aa_pos","distance","errors"),
           sep = "\\|", fill = "right") %>%
  dplyr::select(ID,CHROM, POS, REF, ALT, gene, effect, impact, hgvs_p)

ann.medium.high <- ann[ann$impact %in% c("MODERATE","HIGH"),]

gt.long <- extract.gt(vcf, element = "GT") %>%
  as_tibble(rownames = "ID") %>%
  pivot_longer(-ID, names_to = "sample.id", values_to = "GT") %>%
  mutate(mutated = !is.na(GT) & GT != "0/0")

gt.long.mutated <- gt.long[gt.long$mutated==TRUE,] 

mutations <- gt.long.mutated %>% inner_join(ann.medium.high, by = "ID")

# Remove "mutated" column
mutations <- mutations[, !names(mutations) %in% "mutated"]

# Add isolate info
mutations <- merge(mutations, samples, by.x="sample.id", by.y = "sample.id.full")

# Add gene name
mutations <- merge(mutations, gene.info, by.x="gene", by.y="ID")
mutations[is.na(mutations$gene.name),]$gene.name <- mutations[is.na(mutations$gene.name),]$gene

# flag mutations present in both replicates with identical GT
mutations <- mutations |>
  dplyr::group_by(isolate, CHROM, POS, REF, ALT, effect, GT, ID, hgvs_p, gene, gene.name) |>
  dplyr::mutate(in.both.replicates = ifelse(dplyr::n_distinct(sample.id) == 2, "yes", "no")) |>
  dplyr::ungroup()

# Remove sample id column
mutations <- mutations[, !names(mutations) %in% "sample.id"]

# Convert HGVS protein notation from 3-letter to 1-letter and remove "p."
aa.map <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
  Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
  Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
  Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V"
)

mutations$mutation <- sub("^p\\.", "", mutations$hgvs_p)

for(i in names(aa.map)){
  mutations$mutation <- gsub(i, aa.map[i], mutations$mutation, fixed=TRUE)
}

# Remove duplicate rows, coming from duplicates of the same isolate
mutations <- unique(mutations)

dim(mutations[which(mutations$in.both.replicates=="no"),])

# Identify genes of interest -------

# QS-related genes provided by expert, grouped as a single reference vector
qs.genes <- c(
  "lasR","lasI","lasA","lasB",
  "rhlR","rhlI","rhlA","rhlB","rhlC",
  paste0("phz", c("A1","B1","C1","D1","E1","F1","G1",
                  "A2","B2","C2","D2","E2","F2","G2")),
  "phzM","phzS",
  "mvfR", 
  paste0("pqs", c("A","B","C","D","E","H","L","R")))
  
additional.qs.genes <- c("toxA","hcnA","hcnB","hcnC", "lasA" ,"lasB")

# Subsets expected to impact individual phenotypes if mutated

# Las signal (3-oxo-C12-HSL output / las reporter)
las.genes <- c(
  "lasR","lasI","lasA","lasB",
  "rhlR","rhlI",
  "pqsR","pqsH","mvfR"
)

# Rhl signal (C4-HSL output / rhl reporter)
rhl.genes <- c(
  "rhlR","rhlI","rhlA","rhlB","rhlC",
  "lasR","lasI",
  "pqsR","mvfR"
)

# Pyocyanin / phenazines
pyocyanin.genes <- c(
  paste0("phz", c("A1","B1","C1","D1","E1","F1","G1",
                  "A2","B2","C2","D2","E2","F2","G2")),
  "phzM","phzS",
  "pqsA","pqsB","pqsC","pqsD","pqsE","pqsH","pqsR",
  "lasR","rhlR","mvfR"
)

# Rhamnolipids
rhamnolipid.genes <- c(
  "rhlA","rhlB","rhlC",
  "rhlR","rhlI",
  "lasR",
  "pqsR","mvfR"
)

# Additional QS-controlled virulence factors (optional, for completeness)
#toxA.genes <- c("toxA","lasR","pqsR")
#hcn.genes  <- c("hcnA","hcnB","hcnC","lasR","rhlR","pqsR")

# Select main genes to include in the mutation heatmap -------
#key.genes <- main.genes
key.genes <- qs.genes
#key.genes <- unique(c(key.genes, additional.qs.genes))

mutations.selected <- mutations[mutations$gene.name %in% key.genes, ]
 
# Organise data for plotting heatmap -------------

# Scale phenotypes 0-1, reorder, and rename las/rhl
pheno.long <- pheno.data %>%
  pivot_longer(cols = las:rhamnolipid, names_to = "phenotype", values_to = "value") %>%
  group_by(phenotype) %>%
  mutate(value_scaled = (value - min(value)) / (max(value) - min(value))) %>%  # min-max scaling
  ungroup() %>%
  mutate(
    phenotype = recode(phenotype, las = "las (9h)", rhl = "rhl (9h)"),
    phenotype = factor(phenotype, levels = c("las (9h)", "rhl (9h)", "pyocyanin", "rhamnolipid")))

mut.plot = mutations.selected
mut.plot <- mutations.selected %>%
  mutate(mut.key = paste(gene.name, POS, REF, ALT, sep = "_")) #%>%

mut.plot <- bind_rows(mut.plot, c(gene = NA,ID = NA,GT = NA,CHROM = NA,POS = NA,REF = NA,ALT = NA,effect = NA,impact = NA,
         hgvs_p = NA,isolate = "PAO1",MLST_type = NA,gene.name = NA,in.both.replicates = "no",mut.key = NA)) 

# Order everything by the order from the phylogenetic tree 
iso.order <- isolate_keep_order
pheno.long <- pheno.long %>% mutate(isolate = factor(isolate, levels = iso.order))  
mut.plot <- mut.plot %>% mutate(isolate = factor(isolate, levels = iso.order))

# Collapsing coocuring mutations -----

mut.signature <- mut.plot %>%
  group_by(gene, mut.key,impact ,mutation) %>%
  summarize(isolate.set = paste(sort(unique(as.character(isolate))), collapse = "_"),.groups = "drop")

cooccur <- mut.signature %>%
  distinct(gene, mut.key, mutation, impact,isolate.set) %>%
  group_by(gene, isolate.set,impact) %>%
  summarize(
    muts   = paste(mut.key, collapse = ";"),
    n.muts = n(),
    .groups = "drop"
  )

mut.plot <- mut.plot %>%
  group_by(gene, mut.key, impact) %>%
  mutate(isolate.set = paste(sort(unique(as.character(isolate))), collapse = "_")) %>%
  ungroup()

collapsed <- mut.plot %>%
  distinct(gene, gene.name, mut.key, isolate.set, impact, mutation) %>%
  group_by(gene, gene.name, isolate.set, impact) %>%
  summarize(
    n.muts = n(),
    .groups = "drop"
  )

collapsed <- collapsed %>%
  group_by(gene.name) %>%
  arrange(desc(n.muts), isolate.set, .by_group = TRUE) %>%
  mutate(mut.id = row_number()) %>%
  ungroup()

collapsed <- collapsed %>%
  mutate(
    mut.label = ifelse(
      n.muts > 1,
      paste0(gene.name, ".mut", mut.id, ".c", n.muts),
      paste0(gene.name, ".mut", mut.id)
    )
  )

# Replace mut.label only for PAO1 rows
collapsed <- collapsed %>%
  mutate(mut.label = ifelse(isolate.set == "PAO1", "z.dummy", mut.label))

mut.plot.collapsed <- mut.plot %>%
  left_join(
    collapsed %>% dplyr::select(gene, gene.name, isolate.set, mut.label,n.muts),
    by = c("gene","gene.name","isolate.set")
  )

mut.plot.collapsed$mut.label.old <- mut.plot.collapsed$mut.label

mut.plot.collapsed <- mut.plot.collapsed %>%
  group_by(mut.label.old) %>%
  mutate(
    mut.label = if(n.muts[1] == 1) {
      paste0(gene.name, ".", mutation)
    } else {
      # extract numeric position from mutation string
      aa.pos <- as.numeric(gsub("[A-Z]", "", mutation))
      paste0(gene.name, ".", min(aa.pos), "-", max(aa.pos))
    }
  ) %>%
  ungroup()

# create a combined label factor for plotting, ordered by gene
# you can keep mut.label within each gene as-is, or sort alphabetically
mut.plot.collapsed <- mut.plot.collapsed %>%
  mutate(gene.name = factor(gene.name, levels = rev(key.genes))) %>%
  arrange(gene.name) %>%   # sort by gene order
  mutate(mut.label = factor(mut.label, levels = unique(mut.label)))

# create a gene factor based on qs.genes
mut.plot.collapsed$gene.name <- factor(mut.plot.collapsed$gene.name,
                                       levels = rev(key.genes))

# Plot heatmaps long format ------------
p.pheno <- ggplot(pheno.long,aes(x = isolate, y = phenotype,fill = value_scaled)) +
  geom_tile(color = "white",linewidth=0.4) +
  scale_fill_gradient(low = "#fff5f0", high = "darkblue") +
  scale_y_discrete(position = "right")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10 ),
        axis.title.x = element_blank(),
        axis.text.y =element_text(size=13),
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
)

p.mut <- ggplot(mut.plot.collapsed,aes(x = isolate, y = mut.label, fill = impact)) +
  geom_tile(color = "white",linewidth=0.4) +
  scale_fill_manual(values = c(
    HIGH = "#335A30",
    MODERATE ="#94C47D"
  )) +
  scale_y_discrete(position = "right")+
  labs(x = "Isolate", y = "Mutation") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9))
        #panel.background = element_rect(fill = "white", color = NA),
        #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)))

p.pheno / p.mut +
  plot_layout(heights = c(0.8, 7))

ggsave(paste0(figures,"phenotype_mutation_heatmap_collapsed.pdf"), width = 6, height = 12)
ggsave(paste0(figures,"phenotype_mutation_heatmap_collapsed.png"), width = 6, height = 12, dpi = 600)

# Plot heatmaps long format without pqs genes------------

pqs.genes = c("pqsA","pqsB","pqsC","pqsD","pqsE","pqsH","pqsL")
mut.plot.collapsed.without.pqs <- mut.plot.collapsed[!mut.plot.collapsed$gene.name %in% pqs.genes,] 

p.pheno <- ggplot(pheno.long,aes(x = isolate, y = phenotype,fill = value_scaled)) +
  geom_tile(color = "white",linewidth=0.4) +
  scale_fill_gradient(low = "#fff5f0", high = "darkblue") +
  scale_y_discrete(position = "right")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 10 ),
        axis.title.x = element_blank(),
        axis.text.y =element_text(size=11),
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  )

p.mut <- ggplot(mut.plot.collapsed.without.pqs,aes(x = isolate, y = mut.label, fill = impact)) +
  geom_tile(color = "white",linewidth=0.4) +
  scale_fill_manual(values = c(
    HIGH = "#335A30",
    MODERATE ="#94C47D"
  )) +
  scale_y_discrete(position = "right")+
  labs(x = "Isolate", y = "Mutation") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 9))
#panel.background = element_rect(fill = "white", color = NA),
#theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)))

p.pheno / p.mut +
  plot_layout(heights = c(0.8, 5.5))

ggsave(paste0(figures,"phenotype_mutation_heatmap_collapsed_without_pqs.pdf"), width = 5.5, height = 9)
ggsave(paste0(figures,"phenotype_mutation_heatmap_collapsed_without_pqs.png"), width = 5, height = 9, dpi = 600)

# Cluster mutation rows based on isolate patterns -------
mut.mat <- mut.plot.collapsed %>%
  mutate(val = 1) %>%
  dplyr::select(mut.label, isolate, val) %>%
  pivot_wider(id_cols = mut.label, names_from = isolate, values_from = val, values_fn = max, values_fill = 0) %>%
  column_to_rownames("mut.label") %>%
  as.matrix()

mut.order <- rownames(mut.mat)[hclust(dist(mut.mat))$order]

mut.plot.clust <- mut.plot.collapsed %>%
  mutate(mut.label = factor(mut.label, levels = mut.order))

p.pheno <- ggplot(pheno.long,aes(x = isolate, y = phenotype,fill = value_scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#fff5f0", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank())

p.mut <- ggplot(mut.plot.clust ,aes(x = isolate, y = mut.label, fill = impact)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(
    HIGH = "#335A30",
    MODERATE = "#94C47D"
  )) +
  labs(x = "Isolate", y = "Mutation") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p.pheno / p.mut +
  plot_layout(heights = c(1, 7))

ggsave(paste0(figures,"phenotype_mutation_heatmap_collapsed_clustered_rows.pdf"), width = 6, height = 12)

# Plot pheno subsets -----------

plot.pheno.mut.heatmap <- function(gene.subset, phenotype.selected, name){
  
  mut.plot.subset <- mut.plot.collapsed %>%
    filter(gene.name %in% c(gene.subset) | isolate == "PAO1")
  
  #pheno.long.subset <- pheno.long %>%
  #  filter(phenotype %in% c(phenotype.selected))
  
  # isolate order based on phenotype value
  #isolate.order <- pheno.long.subset %>% arrange(value) %>%                  
  #  pull(isolate) %>% as.character() %>% unique()
  
  #pheno.long.subset <- pheno.long.subset %>%
  #  mutate(isolate = factor(isolate, levels = isolate.order))
  
  # isolate order based on selected phenotype only
  isolate.order <- pheno.long %>%
    filter(phenotype == phenotype.selected) %>% 
    arrange(value) %>%
    pull(isolate) %>%
    as.character() %>%
    unique()
  
  # apply order to full phenotype table
  pheno.long.subset <- pheno.long %>%
    mutate(isolate = factor(isolate, levels = isolate.order))
  
  
  mut.plot.subset <- mut.plot.subset %>%
    mutate(isolate = factor(isolate, levels = isolate.order))
  
  heatmap.height = length(unique(mut.plot.subset$mut.label)) * 0.2 + 1  # adjust height based on number of mutations
  

  # cluster mutation rows based on isolate patterns -------
  mut.mat <- mut.plot.subset %>%
    mutate(val = 1) %>%
    dplyr::select(mut.label, isolate, val) %>%
    pivot_wider(id_cols = mut.label, names_from = isolate, values_from = val, values_fn = max, values_fill = 0) %>%
    column_to_rownames("mut.label") %>%
    as.matrix()
  
  mut.order <- rownames(mut.mat)[hclust(dist(mut.mat))$order]
  
  mut.plot.subset <- mut.plot.subset %>%
    mutate(mut.label = factor(mut.label, levels = mut.order))
  # ----------
  
  p.pheno <- ggplot(pheno.long.subset,aes(x = isolate, y = phenotype,fill = value_scaled)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "#fff5f0", high = "darkblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank())
  
  p.mut <- ggplot(mut.plot.subset,aes(x = isolate, y = mut.label, fill = impact)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c(
      HIGH = "#335A30",
      MODERATE = "#94C47D"
    )) +
    labs(x = "Isolate", y = "Mutation") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  p.pheno + p.mut +
    plot_layout(heights = c(1, heatmap.height-1))
  
  ggsave(paste0(figures,"phenotype_mutation_heatmap_collapsed_",name, ".pdf"), width = 5, height = heatmap.height)
}

plot.pheno.mut.heatmap(las.genes, "las (9h)", "las")
plot.pheno.mut.heatmap(rhl.genes, "rhl (9h)", "rhl")
plot.pheno.mut.heatmap(pyocyanin.genes, "pyocyanin", "pyocyanin")
plot.pheno.mut.heatmap(rhamnolipid.genes, "rhamnolipid", "rhamnolipid")

# Correlation analysis -------

mut.mat <- mut.plot.collapsed %>%
  dplyr::select(isolate, mut.label) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from  = mut.label,
    values_from = present,
    values_fill = 0
  )

phenotype.selected <- "rhl (9h)"
phenotype.selected <- "las (9h)"
phenotype.selected <- "pyocyanin"
phenotype.selected <- "rhamnolipid"

pheno.vec <- pheno.long %>%
  filter(phenotype == phenotype.selected) %>%
  dplyr::select(isolate, value)

dat <- mut.mat %>%
  left_join(pheno.vec, by = "isolate")

assoc.results <- dat %>%
  pivot_longer(-c(isolate, value),
               names_to = "mut.label",
               values_to = "present") %>%
  group_by(mut.label) %>%
  summarize(
    n.mut   = sum(present),
    n.wt    = sum(present == 0),
    effect  = mean(value[present == 1], na.rm = TRUE) -
      mean(value[present == 0], na.rm = TRUE),
    p.value = ifelse(n.mut >= 2 & n.wt >= 2,
                     wilcox.test(value ~ present)$p.value,
                     NA_real_),
    .groups = "drop"
  ) %>%
  arrange(p.value)

hist(assoc.results$p.value, breaks = 20, main = paste("P-value distribution for", phenotype.selected),xlab = "P-value")
hist(assoc.results$effect, breaks = 20, main = paste("Effect size distribution for", phenotype.selected),xlab = "Effect size")

############################
# Plot heatmaps wide format ------------
p.pheno <- ggplot(pheno.long,aes(x = phenotype, y = isolate,fill = value_scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#fff5f0", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_blank())

p.mut <- ggplot(mut.plot.collapsed,aes(x = mut.label, y = isolate, fill = impact)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(
    HIGH = "#335A30",
    MODERATE = "#94C47D"
  )) +
  labs(x = "Isolate", y = "Mutation") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p.pheno + p.mut + plot_layout(widths = c(1, 5))

# Combine side by side, collect legends on the right
p.pheno + p.mut +
  plot_layout(widths = c(1, 7), guides = "collect") & 
  theme(legend.position = "right")

ggsave(paste0(figures,"phenotype_mutation_heatmap_wide.pdf"), width = 14, height = 6)

# Plot clustered heatmap  ------------

# cluster mutations (rows) and isolates (columns)
mat <- xtabs(~ mut.label + isolate, mut.plot.collapsed)
iso.order <- colnames(mat)[hclust(dist(t(mat)))$order]
mut.order <- rownames(mat)[hclust(dist(mat))$order]

mut.plot.clust <- mut.plot.collapsed

mut.plot.clust$isolate <- factor(mut.plot$isolate, levels = iso.order)
pheno.long.clust <- pheno.long
pheno.long.clust$isolate <- factor(pheno.long$isolate, levels = iso.order)

#mut.plot.clust$mut.label <- factor(mut.plot$mut.label, levels = mut.order)

p.mut <- ggplot(mut.plot.clust, aes(x = isolate, y = mut.label, fill = impact)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(HIGH = "#335A30", MODERATE = "#94C47D")) +
  labs(x = "Isolate", y = "Mutation") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p.pheno <- ggplot(pheno.long,aes(x = isolate, y = phenotype,fill = value_scaled)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "#fff5f0", high = "darkblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p.pheno / p.mut +
  plot_layout(heights = c(1, 5))

# Save mutations as xlsx ------
# Reshape mutation table to wide (1/0 per isolate) and export as xlsx

mut.wide <- mut.plot.collapsed %>%
  mutate(present = 1) %>%
  select(gene, gene.name, mutation, mut.label, mut.label.old, ID, effect, impact,  isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order[isolate_keep_order != "PAO1"]), .after = impact)

write.xlsx(mut.wide, paste0(resultdir,"mutations_QS.xlsx"))
# Get all mutations, for plotting general stats ----------

mutations.full <- gt.long.mutated %>% inner_join(ann, by = "ID")

# Remove "mutated" column
mutations.full <- mutations.full[, !names(mutations.full) %in% "mutated"]

# Add isolate info
mutations.full <- merge(mutations.full, samples, by.x="sample.id", by.y = "sample.id.full")

# Add gene name
mutations.full <- merge(mutations.full, gene.info, by.x="gene", by.y="ID",all.x=TRUE)
mutations.full[is.na(mutations.full$gene.name),]$gene.name <- mutations.full[is.na(mutations.full$gene.name),]$gene

# Remove sample id columns, such as to then easily remove duplicate rows, coming from duplicates of the same isolate
mutations.full <- mutations.full[, !names(mutations.full) %in% "sample.id"]
mutations.full <- unique(mutations.full)

# Count mutations per isolate and impact --------
mut.counts <- aggregate(
  list(counts = mutations.full$gene),
  by = list(isolate = mutations.full$isolate,
            impact  = mutations.full$impact),
  FUN = length
)

# set impact order: bottom → top in stacked bars
mut.counts$impact <- factor(
  mut.counts$impact,
  levels = rev(c("MODIFIER","LOW","MODERATE","HIGH"))
)

# set isolate order explicitly
mut.counts$isolate <- factor(
  mut.counts$isolate,
  levels = isolate_keep_order[-c(length(isolate_keep_order))]  # remove PAO1 from barplot
)

# Stacked barplot all mutations  -------

barplot.all <- ggplot(mut.counts, aes(x = isolate, y = counts, fill = impact)) +
    geom_col(color = NA) +  # remove any bar borders
    scale_fill_manual(values = c(
      MODIFIER = "grey80",
      LOW      = "grey65",
      #MODERATE = "#FDB863",
      #HIGH     = "#B2182B"
      MODERATE = "#94C47D",
      HIGH     = "#335A30"
    )) +
    ggtitle("All mutations") +
    theme_classic(base_size = 14) +  # clean minimal theme
    theme(
      plot.title = element_text(size = rel(0.8)),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.grid = element_blank(),       # remove gridlines
      axis.line = element_line(color = "black"),  # only axes in #335A30
      axis.ticks = element_line(color = "black"),
      axis.title.x     = element_blank()    # remove x-axis title
    )  
barplot.all  
ggsave(paste0(figures,"mutations_per_isolate_all.pdf"), width = 7, height = 4)

# Plot barplot with HIGH + MODERATE impact only ------
high.mod.mut <- subset(mut.counts, impact %in% c("HIGH","MODERATE"))

barplot.high.mod <- ggplot(high.mod.mut, aes(x = isolate, y = counts, fill = impact)) +
  geom_col(color = NA) +  # remove any bar borders
  scale_fill_manual(values = c(
    #MODERATE = "#FDB863",
    #HIGH     = "#B2182B"
    MODERATE = "#94C47D",
    HIGH     = "#335A30"
  )) +
  ggtitle("Moderate and high impact mutations") +
  theme_classic(base_size = 14) +  # clean minimal theme
  theme(
    plot.title = element_text(size = rel(0.8)),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),       # remove gridlines
    axis.line = element_line(color = "black"),  # only axes in #335A30
    axis.ticks = element_line(color = "black"),
    axis.title.x     = element_blank()    # remove x-axis title
  )  
barplot.high.mod
ggsave(paste0(figures,"mutations_per_isolate_high.mod.pdf"), width = 7, height = 4)

# Plot barplot with HIGH impact only ------

high.mut <- subset(mut.counts, impact == "HIGH")

# optional: ensure isolates keep the same order
high.mut$isolate <- factor(high.mut$isolate, levels = isolate_keep_order[-c(length(isolate_keep_order))])

# barplot for HIGH impact mutations
barplot.high <- ggplot(high.mut, aes(x = isolate, y = counts)) +
  #geom_col(fill = "#B2182B") +  # same HIGH color
  geom_col(fill = "#335A30") +  # same HIGH color
  labs(x = "Isolate", y = "counts") +
  ggtitle("High impact mutations") +
  theme_classic(base_size = 14) +  # clean minimal theme
  theme(
    plot.title = element_text(size = rel(0.8)),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.title.x     = element_blank()    # remove x-axis title
  )
barplot.high
ggsave(paste0(figures,"mutations_per_isolate_high_impact.pdf"), width = 7, height = 3)

barplot.all / barplot.high +
  plot_layout(heights = c(4, 1),guides = "collect")

ggsave(paste0(figures,"mutations_per_isolate_all_plus_high.pdf"), width = 7, height = 3)

barplot.all + barplot.high.mod + barplot.high +
  plot_layout(heights = c(4,2,1),guides = "collect")

ggsave(paste0(figures,"mutations_per_isolate_all_plus_mod_plus_high.pdf"), width = 7, height = 6)

# Count mutations in QS selected genes ---------

mutations.selected <- mutations[mutations$gene.name %in% qs.genes, ]

mut.counts.qs <- aggregate(
  list(counts = mutations.selected$gene),
  by = list(isolate = mutations.selected$isolate,
            impact  = mutations.selected$impact),
  FUN = length
)

# set impact order: bottom → top in stacked bars
mut.counts.qs$impact <- factor(
  mut.counts.qs$impact,
  levels = rev(c("MODIFIER","LOW","MODERATE","HIGH"))
)

# set isolate order explicitly
mut.counts.qs$isolate <- factor(
  mut.counts.qs$isolate,
  levels = isolate_keep_order[-c(length(isolate_keep_order))]  # remove PAO1 from barplot
)

qs.mut <- mut.counts.qs

# Plot barplot with QS related mutations ------
barplot.high.mod.QS <- ggplot(qs.mut, aes(x = isolate, y = counts, fill = impact)) +
  geom_col(color = NA) +  # remove any bar borders
  scale_fill_manual(values = c(
    #MODERATE = "#FDB863",
    #HIGH     = "#B2182B"
    MODERATE = "#94C47D",
    HIGH     = "#335A30"
  )) +
  theme_classic(base_size = 14) +  # clean minimal theme
  ggtitle("Mutations in QS-related genes") +
  scale_y_continuous(breaks = c(0,max(qs.mut$counts))) +
  
  theme(
    plot.title = element_text(size = rel(0.8)),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),       # remove gridlines
    axis.line = element_line(color = "black"),
    axis.title.x     = element_blank(),    # remove x-axis title
    axis.ticks = element_line(color = "black")
  )  
barplot.high.mod.QS
ggsave(paste0(figures,"mutations_per_isolate_QS.pdf"), width = 7, height = 2)

barplot.all + barplot.high.mod + barplot.high + barplot.high.mod.QS +
  plot_layout(heights = c(1.5,1.5,1,1),guides = "collect")

ggsave(paste0(figures,"mutations_per_isolate_all_plus_QS.pdf"), width = 8, height = 7)

  

# Count mutations in QS selected genes, separated into phz and pqs ---------

mutations.phz <- mutations[mutations$gene.name %in% paste0("phz", c("A1","B1","C1","D1","E1","F1","G1","A2","B2","C2","D2","E2","F2","G2","M","S")), ]
mutations.pqs <- mutations[mutations$gene.name %in% c("pqsA","pqsB","pqsC","pqsD","pqsE","pqsH","pqsR"), ]
mutations.rhl <- mutations[mutations$gene.name %in% c("rhlR","rhlI","rhlA","rhlB","rhlC"), ]
mutations.las <- mutations[mutations$gene.name %in% c("lasR","lasI"), ]

mutations.subset = mutations.phz; title.subset = "phz"
mutations.subset = mutations.pqs; title.subset = "pqs"
mutations.subset = mutations.rhl; title.subset = "rhl"
mutations.subset = mutations.las; title.subset = "las"

mut.counts.qs.subset <- aggregate(
  list(counts = mutations.subset$gene),
  by = list(isolate = mutations.subset$isolate,
            impact  = mutations.subset$impact),
  FUN = length
)

# set impact order: bottom → top in stacked bars
mut.counts.qs.subset$impact <- factor(
  mut.counts.qs.subset$impact,
  levels = rev(c("MODIFIER","LOW","MODERATE","HIGH"))
)

# set isolate order explicitly
mut.counts.qs.subset$isolate <- factor(
  mut.counts.qs.subset$isolate,
  levels = isolate_keep_order[-c(length(isolate_keep_order))]  # remove PAO1 from barplot
)

qs.mut.subset <- mut.counts.qs.subset

# Plot barplot with QS related mutations separated into phz and pqs ------
barplot.high.mod.QS.subset <- ggplot(qs.mut.subset, aes(x = isolate, y = counts, fill = impact)) +
  geom_col(color = NA) +  # remove any bar borders
  scale_fill_manual(values = c(
    MODERATE = "#94C47D",
    HIGH     = "#335A30"
  )) +
  ggtitle(paste("Mutations in",title.subset, "genes")) +
  theme_classic(base_size = 14) +  # clean minimal theme
  ggtitle("Mutations in QS-related genes") +
  theme(
    plot.title = element_text(size = rel(0.8)),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),       # remove gridlines
    axis.line = element_line(color = "black"),  # only axes in #335A30
    axis.ticks = element_line(color = "black")
  )  
barplot.high.mod.QS.subset
ggsave(paste0(figures,"mutations_per_isolate_QS_",title.subset,".pdf"), width = 7, height = 2)

barplot.all + barplot.high.mod + barplot.high + barplot.high.mod.QS +
  plot_layout(heights = c(1.5,1.5,1,1),guides = "collect")

ggsave(paste0(figures,"mutations_per_isolate_all_plus_QS.pdf"), width = 8, height = 7)

# Plot barplot with QS mutations by categories v2---------
mutations.qs <- mutations %>%
  mutate(qs.class = case_when(
    gene.name %in% paste0("phz", c("A1","B1","C1","D1","E1","F1","G1",
                                   "A2","B2","C2","D2","E2","F2","G2","M","S")) ~ "phz",
    gene.name %in% c("pqsA","pqsB","pqsC","pqsD","pqsE","pqsH","pqsR","mvfR")   ~ "pqs",
    gene.name %in% c("rhlR","rhlI") ~ "rhl",
    gene.name %in% c("rhlA","rhlB","rhlC") ~ "rhlA-C",
    gene.name %in% c("lasR","lasI")  ~ "las",
    gene.name %in% c("lasA","lasB") ~ "lasA-B",
    TRUE  ~ NA_character_
  )) %>%
  filter(!is.na(qs.class))

# count mutations
qs.mut.counts <- mutations.qs %>%
  group_by(isolate, impact, qs.class) %>%
  summarise(counts = n(), .groups = "drop")

qs.mut.counts <- qs.mut.counts %>%
  mutate(
    impact = factor(impact, levels = rev(c("MODIFIER","LOW","MODERATE","HIGH"))),
    isolate = factor(isolate, levels = isolate_keep_order[-length(isolate_keep_order)]),
    qs.class = factor(qs.class, levels = c("las","lasA-B","rhl","rhlA-C","pqs","phz"))
  )

# faceted barplot
p.qs <- ggplot(qs.mut.counts, aes(x = isolate, y = counts, fill = impact)) +
  geom_col(color = NA) +
  scale_fill_manual(values = c(
    MODERATE = "#94C47D",
    HIGH     = "#335A30"
  )) +
  scale_y_continuous(breaks = c(0,10)) +
  facet_wrap(~ qs.class, nrow = 6) +
  theme_classic(base_size = 14) +
  ggtitle("Mutations in subsets of QS-related genes") +
  theme(
    plot.title = element_text(size = rel(0.8)),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p.qs

ggsave(paste0(figures,"mutations_per_isolate_QS_subsets.pdf"), width = 7, height = 8)

barplot.all + barplot.high.mod + barplot.high.mod.QS + p.qs +
  plot_layout(heights = c(1,1,1,10),guides = "collect")

ggsave(paste0(figures,"barplot_mutations_all_plus_QS_subsets.pdf"), width = 6, height = 9)

# Fig. 6: all, high/moderate, resistome and QS variant counts per isolate -------
# Uses the resistome panels saved by script 07
res.panels <- readRDS(paste0(resultdir, "fig6_resistome_panels.rds"))

res.panels$barplot.all + res.panels$barplot.high.mod + res.panels$barplot.all.res + res.panels$p.res +
  barplot.high.mod.QS + p.qs +
  plot_layout(heights = c(1,1,1,12,1,10),guides = "collect")

ggsave(paste0(figures,"barplot_mutations_all_plus_resistome_plus_QS_subsets.pdf"), width = 7, height = 16)
ggsave(paste0(figures,"barplot_mutations_all_plus_resistome_plus_QS_subsets.png"), width = 7, height = 16)

