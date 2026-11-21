# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/07_resistance_mutations/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures, recursive = TRUE), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages-------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx", "pheatmap","tibble",
                      "ComplexHeatmap","circlize", "vcfR","purrr","rtracklayer","patchwork","openxlsx","gsubfn")

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

# Read disc diffusion results and combine with resistance predicted phenotypes --------
disc.diff <- read_excel(paste(datadir,"disc_diffusion_results.xlsx",sep=""))
#Select for P aeruginoa rows
disc.diff = disc.diff[disc.diff$Species=="P. aeruginosa", ]
# Remove columns with all NAs
disc.diff <- disc.diff %>%
  dplyr::select(where(~ !all(is.na(.))))
# Remove species column
disc.diff <- disc.diff[,-c(2) ]
# Set rownames to Isolate
disc.diff = as.data.frame(disc.diff)
colnames(disc.diff)[1] <- "Isolate"
rownames(disc.diff) <- disc.diff$Isolate
# Rename Piperacilin / Tazobactam column to match the name in resfinder predicted phenotype
colnames(disc.diff)[4] = c("Piperacillin+Tazobactam")
colnames(disc.diff)[7] = c("Ceftolozane+Tazobactam")

# Remove Isolate ID column
disc.diff <- disc.diff[,-c(1) ]
# Reorder antibiotics
abx_order <- rev(c("Ciprofloxacin","Ceftazidime","Piperacillin+Tazobactam","Aztreonam","Meropenem","Ceftolozane+Tazobactam"))
disc.diff <- disc.diff[,abx_order, drop = FALSE]
disc.diff <- disc.diff[isolate_keep_order, drop = FALSE,]

# parse GFF, to be able to match gene ID from vcf to gene name ------
gff <- import(paste0(datadir,"bacdist/ref.gff"))
gene.info <- as_tibble(gff) 
gene.info <- gene.info[c("ID","gene")]
colnames(gene.info)[2] = "gene.name"

# Read vcf file ---------

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

# Count number of mutations per sample.id
mutations.summary <- mutations %>%
  dplyr::count(sample.id, name = "n.mutations")

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

#mutations.including.not.in.both.replicates <- mutations
#mutations <- mutations[mutations$in.both.replicates=="yes",]

#dim(mutations.including.not.in.both.replicates[which(mutations.including.not.in.both.replicates$in.both.replicates=="no"),])

# Identify genes of interest -------

# Clinically relevant efflux pump structural genes in P. aeruginosa
efflux.genes <- c(
  "mexA","mexB","oprM",        # MexAB-OprM (major MDR pump)
  "mexC","mexD","oprJ",        # MexCD-OprJ
  "mexE","mexF","oprN",        # MexEF-OprN
  "mexX","PA2019",
  "mexY","PA2018",
  # MexXY-OprM (aminoglycosides + FQ)
  "oprM"                       # shared outer membrane component
)

# Regulators of major efflux systems
efflux.regulators <- c(
  "mexR","nalC","nalD",        
  "nfxB",                      
  "mexT", 
  "mexS", "PA2491",                    
  "mexZ","PA2020"                       
)

# Porins affecting beta-lactam and carbapenem uptake
porin.genes <- c(
  "oprD"                       # primary carbapenem porin
)

# AmpC and its regulatory pathway
ampc.genes <- c(
  "ampC",
  "ampD","ampE",
  "ampR",
  "dacB", "PA3047" #(dacB)
)

# Acquired beta-lactamase families to screen
#betalactamase.genes <- c(
#  "blaIMP","blaVIM","blaNDM",   # metallo-beta-lactamases
#  "blaGES","blaKPC",            # serine carbapenemases
#  "blaOXA","blaPER","blaVEB"   # ESBLs relevant to P. aeruginosa
#)

# Fluoroquinolone resistance determinants
fluoroquinolone.genes <- c(
  "gyrA","gyrB",
  "parC","parE",
  "mexR","nalC","nalD","mexT","nfxB","mexZ"
)

# LPS modification and peptide resistance
lps.genes <- c(
  "pmrA","pmrB",
  "phoP","phoQ",
  "arnA","arnB","arnC","arnD","arnE","arnF",
  "pagL"
)

# Penicillin-binding proteins and cell wall enzymes # PBP3 (ceftazidime)
pbp.genes <- c(
  "pbpA","pbpB","pbpC",
  "ftsI"                       
)

# Master resistome gene list for P. aeruginosa (clinically focused)
resistome.genes <- unique(c(
  efflux.genes,
  efflux.regulators,
  porin.genes,
  ampc.genes,
  fluoroquinolone.genes,
  lps.genes,
  pbp.genes
))

# Select main genes to include in the mutation heatmap -------
key.genes <- resistome.genes
mutations.selected <- mutations[mutations$gene.name %in% key.genes, ]
mutations.selected[mutations.selected$gene.name == "PA3047",]$gene.name <- "dacB"
mutations.selected[mutations.selected$gene.name == "PA2491",]$gene.name <- "mexS"
mutations.selected[mutations.selected$gene.name == "PA2019",]$gene.name <- "mexX"
mutations.selected[mutations.selected$gene.name == "PA2018",]$gene.name <- "mexY"
mutations.selected[mutations.selected$gene.name == "PA2020",]$gene.name <- "mexZ"

# Organise data for plotting heatmap -------------

# Prepare disc diffusion data
disc.long <- disc.diff %>%
  rownames_to_column("isolate") %>%
  pivot_longer(
    cols = -isolate,
    names_to = "antibiotic",
    values_to = "phenotype"
  ) %>%
  mutate(
    isolate = factor(isolate, levels = isolate_keep_order),
    antibiotic = factor(
      antibiotic,
      levels = c(
        "Ciprofloxacin","Ceftazidime", "Piperacillin+Tazobactam", "Aztreonam", "Meropenem","Ceftolozane+Tazobactam" )),
    phenotype = factor(phenotype, levels = c("S", "I", "R"))
  )

mut.plot = mutations.selected
mut.plot <- mutations.selected %>%
  mutate(mut.key = paste(gene.name, POS, REF, ALT, sep = "_")) #%>%

#Remove PAO1 from disc diff
disc.long <- disc.long %>% filter(isolate != "PAO1")

# Order everything by the order from the phylogenetic tree 
iso.order <- isolate_keep_order
disc.long <- disc.long %>% mutate(isolate = factor(isolate, levels = iso.order))  
mut.plot <- mut.plot %>% mutate(isolate = factor(isolate, levels = iso.order))

# Collapse co-occuring mutations -----------
# create a string of isolates per mutation, *within each gene*
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

# create a gene factor based on qs.genes
mut.plot.collapsed$gene.name <- factor(mut.plot.collapsed$gene.name,
                                       levels = rev(key.genes))

# create a combined label factor for plotting, ordered by gene
# you can keep mut.label within each gene as-is, or sort alphabetically
mut.plot.collapsed <- mut.plot.collapsed %>%
  arrange(gene.name) %>%   # sort by gene order
  mutate(mut.label = factor(mut.label, levels = unique(mut.label)))

# Plot heatmaps long format ------------
p.abx <- ggplot(disc.long, aes(x = isolate, y = antibiotic, fill = phenotype)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = c(
      S = "darkblue",
      I = "#fdae61",  # orange
      R = "#d73027"   # red
    ),
    drop = FALSE
  ) +
  labs(x = NULL, y = NULL, fill = "Phenotype") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

p.mut <- ggplot(mut.plot.collapsed,aes(x = isolate, y = mut.label, fill = impact)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c(
    HIGH = "#335A30",
    MODERATE = "#94C47D"
  )) +
  labs(x = "Isolate", y = "Mutation") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

p.abx / p.mut +
  plot_layout(heights = c(1, 10))

ggsave(paste0(figures,"mutations_resistance_disc_diff.pdf"), width = 8, height = 15)

# Plot subsets of resistance genes -----------

plot.mut.heatmap <- function(gene.subset, name){
  
  mut.plot.subset <- mut.plot.collapsed %>%
    filter(gene.name %in% c(gene.subset))

  heatmap.height = length(unique(mut.plot.subset$mut.label)) * 0.2 + 1  # adjust height based on number of mutations
  
  # cluster mutation rows based on isolate patterns -------
  mut.mat <- mut.plot.subset %>%
    mutate(val = 1) %>%
    dplyr::select(mut.label, isolate, val) %>%
    pivot_wider(id_cols = mut.label, names_from = isolate, values_from = val, values_fn = max, values_fill = 0) %>%
    column_to_rownames("mut.label") %>%
    as.matrix()
  
  # ----------
  
  p.abx <- ggplot(disc.long, aes(x = isolate, y = antibiotic, fill = phenotype)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c(
        S = "darkblue",
        I = "#fdae61",  # orange
        R = "#d73027"   # red
      ),
      drop = FALSE
    ) +
    labs(x = NULL, y = NULL, fill = "Phenotype") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  p.mut <- ggplot(mut.plot.subset,aes(x = isolate, y = mut.label, fill = impact)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c(
      HIGH = "#335A30",
      MODERATE = "#94C47D"
    )) +
    labs(x = "Isolate", y = "Mutation") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  
  p.abx + p.mut +
    plot_layout(heights = c(1, heatmap.height-1))
  
  ggsave(paste0(figures,"mutations_resistance_",name, ".pdf"), width = 6, height = heatmap.height)
}

plot.mut.heatmap(c(efflux.genes,efflux.regulators), "efflux.genes.and.regulators")
plot.mut.heatmap(c(porin.genes,ampc.genes), "porins.and.ampC")
plot.mut.heatmap(fluoroquinolone.genes, "fluoroquinolone.genes")
plot.mut.heatmap(c(lps.genes,pbp.genes), "lps.pbp")
plot.mut.heatmap(c(efflux.regulators), "efflux.regulators")

############################

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
  scale_y_continuous(breaks = c(0,90000)) +
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
  scale_y_continuous(breaks = c(0,max(high.mod.mut$counts))) +
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

# Plot barplot with resistance related mutations ------

mutations.selected <- mutations[mutations$gene.name %in% resistome.genes, ]

mut.counts.res <- aggregate(
  list(counts = mutations.selected$gene),
  by = list(isolate = mutations.selected$isolate,
            impact  = mutations.selected$impact),
  FUN = length
)

# set impact order: bottom → top in stacked bars
mut.counts.res$impact <- factor(
  mut.counts.res$impact,
  levels = rev(c("MODIFIER","LOW","MODERATE","HIGH"))
)

# set isolate order explicitly
mut.counts.res$isolate <- factor(
  mut.counts.res$isolate,
  levels = isolate_keep_order[-c(length(isolate_keep_order))]  # remove PAO1 from barplot
)

res.mut <- mut.counts.res

barplot.all.res <- ggplot(res.mut, aes(x = isolate, y = counts, fill = impact)) +
  geom_col(color = NA) +  # remove any bar borders
  scale_fill_manual(values = c(
    #MODERATE = "#FDB863",
    #HIGH     = "#B2182B"
    MODERATE = "#94C47D",
    HIGH     = "#335A30"
  )) +
  theme_classic(base_size = 14) +  # clean minimal theme
  ggtitle("Mutations resistome genes") +
  scale_y_continuous(breaks = c(0,max(mut.counts.res$counts))) +
  theme(
    plot.title = element_text(size = rel(0.8)),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),       # remove gridlines
    axis.line = element_line(color = "black"),  # only axes in #335A30
    axis.title.x     = element_blank(),    # remove x-axis title
    axis.ticks = element_line(color = "black")
  )  
barplot.all.res

# Plot barplot with resistome mutations by categories ---------
mutations.resistome <- mutations %>%
  mutate(res.class = case_when(
    gene.name %in% efflux.genes ~ "efflux",
    gene.name %in% efflux.regulators   ~ "efflux.regulators",
    gene.name %in% porin.genes ~ "porins",
    gene.name %in% ampc.genes ~ "ampC",
    gene.name %in% fluoroquinolone.genes ~ "fluoroquinolone",
    gene.name %in% lps.genes ~ "lps",
    gene.name %in% pbp.genes ~ "pbp",
    TRUE  ~ NA_character_
  )) %>%
  filter(!is.na(res.class))

# count mutations
res.mut.counts <- mutations.resistome %>%
  group_by(isolate, impact, res.class) %>%
  summarise(counts = n(), .groups = "drop")

res.mut.counts <- res.mut.counts %>%
  mutate(
    impact = factor(impact, levels = rev(c("MODIFIER","LOW","MODERATE","HIGH"))),
    isolate = factor(isolate, levels = isolate_keep_order[-length(isolate_keep_order)]),
    res.class = factor(res.class, levels = c("efflux","efflux.regulators","porins","ampC","fluoroquinolone","lps","pbp"))
  )

# faceted barplot
p.res <- ggplot(res.mut.counts, aes(x = isolate, y = counts, fill = impact)) +
  geom_col(color = NA) +
  scale_fill_manual(values = c(
    MODERATE = "#94C47D",
    HIGH     = "#335A30"
  )) +
  scale_y_continuous(breaks = c(0,20)) +
  facet_wrap(~ res.class, nrow = 7) +
  theme_classic(base_size = 14) +
  ggtitle("Mutations in subsets of resistome genes") +
  theme(
    plot.title = element_text(size = rel(0.8)),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

p.res

ggsave(paste0(figures,"barplot_mutations_resistome_subsets.pdf"), width = 7, height = 8)

barplot.all + barplot.high.mod + barplot.all.res + p.res +
  plot_layout(heights = c(1,1,1,10),guides = "collect")

ggsave(paste0(figures,"barplot_mutations_all_plus_resistome_subsets.pdf"), width = 6, height = 9)

# Save resistome panels; the combined Fig. 6 (with QS panels) is drawn at the end of script 08
saveRDS(list(barplot.all = barplot.all, barplot.high.mod = barplot.high.mod,
             barplot.all.res = barplot.all.res, p.res = p.res),
        paste0(resultdir, "fig6_resistome_panels.rds"))

# Save mutations as xlsx ------
# Reshape mutation table to wide (1/0 per isolate) and export as xlsx

mut.wide <- mut.plot.collapsed %>%
  mutate(present = 1) %>%
  select(gene, gene.name, mutation, mut.label, mut.label.old, ID, effect, impact,  isolate, present) %>%
  distinct() %>%
  pivot_wider(names_from = isolate, values_from = present, values_fill = 0) %>%
  relocate(all_of(isolate_keep_order[isolate_keep_order != "PAO1"]), .after = impact)

write.xlsx(mut.wide, paste0(resultdir,"mutations_resistome.xlsx"))

snps.card <- read.csv(paste0(datadir,"selected.snps.CARD.txt"),sep="\t", header=FALSE)
snps.card <- read.table(paste0(datadir,"selected.snps.CARD.overexpression.txt"),header = FALSE,stringsAsFactors = FALSE)

# Mutations also reported in CARD (hatched boxes in Fig. 3B)
mut.in.card <- merge(mut.wide, snps.card, by.x="mutation", by.y="V9")

