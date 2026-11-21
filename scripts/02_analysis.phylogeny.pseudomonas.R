# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/02_snp_phylogeny/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures, recursive = TRUE), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages-------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx", "pheatmap","tibble","ComplexHeatmap","circlize","ape","ggtree")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

lapply(list.of.packages, library, character.only=TRUE)

# Send on-screen plots to a null device when run non-interactively (avoids Rplots.pdf)
if (!interactive()) pdf(NULL)


save_pheatmap_pdf <- function(x, filename, width, height) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

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


# Prepare metadata -----
metadata=samples.isolates[samples.isolates$Species == "Pseudomonas aeruginosa", c("sample.id", "isolate","MLST_type")]
colnames(metadata)[1] <- "label"
metadata$MLST_type <- gsub("!", "", metadata$MLST_type)

metadata <- metadata %>%
  mutate(display_label = paste0(" ", isolate, " | ", label, " | ST", MLST_type))


# Read SNP-based RAxML tree including the PAO1 reference ------
# (BacDist: Snippy SNPs, recombinant regions identified by ClonalFrameML removed, RAxML GTRCAT)
tree <- read.tree(paste0(datadir,"bacdist/RAxML_bestTree.Pa_incl_reference"))
tree$tip.label[tree$tip.label == "AE004091"] <- "PAO1_reference"
tree$tip.label <- sub(".*_(S[0-9]+)$", "\\1", tree$tip.label)


# Plot tree with both sequencing duplicates of each isolate and PAO1 ------
metadata.with.PAO1 <- metadata %>%
  add_row(
    label = "PAO1_reference",
    isolate = "PAO1",
    MLST_type = "549",
    display_label = " PAO1_reference | ST549"
  )

ggtree(tree) %<+% metadata.with.PAO1 + 
  geom_tiplab(aes(label = display_label), size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.5))) 

ggsave(paste0(figures,"Pseudomonas_phylogeny_with_PAO1.pdf"), width = 8, height = 6)


# Plot tree with one representative assembly per isolate and PAO1 (Fig. 1B) -------
# The first sequencing duplicate of each isolate (S3-S33) is kept
metadata.without.duplicates <- metadata %>%
  mutate(S_num = as.numeric(sub("^S", "", label))) %>%
  filter(S_num <= 33) %>%
  dplyr::select(-S_num)

metadata.without.duplicates <- metadata.without.duplicates %>%
  mutate(display_label = paste0(" ", isolate, " | ST", MLST_type))

metadata.without.duplicates.with.PAO1 <- metadata.without.duplicates %>%
  add_row(label = "PAO1_reference",isolate = "PAO1", MLST_type = "549", display_label = " PAO1 | ST549")

tips_to_keep <- metadata.without.duplicates.with.PAO1$label
tree_pruned <- keep.tip(tree, tips_to_keep)

tree_pruned$edge.length <- tree_pruned$edge.length * 0.2  # scale branch lengths for plotting

# Horizontal layout, tips at the bottom, so the tree can sit on top of the heatmaps
p<-ggtree(tree_pruned, layout = "rectangular") %<+% metadata.without.duplicates.with.PAO1 +
  geom_tiplab(aes(label = display_label), size = 4.2, angle = 45, hjust = 1, vjust=2, fontface = "bold") +
  scale_x_reverse(expand = expansion(mult = c(0.2, 2))) +
  scale_y_reverse()+
  coord_flip(clip = "off") +
  theme(plot.margin = margin(t = 10, r = 10, b = 40, l = 10))  # expand bottom margin
p
ggsave(paste0(figures, "Pseudomonas_phylogeny_without_duplicates_with_PAO1_compact_flipped_bold_longer_branches.pdf"),width = 10, height = 5)


# Extract tip order used to sort isolates in the resistance/virulence/disc diffusion heatmaps -------
# (hard-coded as isolate_keep_order in scripts 05-09)
tip_order <- p$data %>%
  filter(isTip) %>%           # only the tips
  arrange(y) %>%              # sort by y-coordinate (top to bottom)
  pull(isolate)             

tip_order <- rev(tip_order)  # reverse the order to match the heatmap orientation
print(tip_order)
