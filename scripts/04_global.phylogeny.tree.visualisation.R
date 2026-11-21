# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/04_global_phylogeny/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures, recursive = TRUE), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages -------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx", "ape","tibble","ggtree","ggtreeExtra","readr","viridis","ggrepel",
                      "sf","rnaturalearth","countrycode","rnaturalearthdata","RColorBrewer","ggnewscale","Polychrome", "phangorn","colorspace","pals")

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

samples.PA = samples.isolates %>% filter(grepl("Pseudomonas aeruginosa", `Species`))
samples.PA$sample.name = paste(samples.PA$sample,samples.PA$sample.id,sep="_")

samples.PA <- samples.PA[c(8,1,2,6,7)]

# Run script 03 to obtain the meta.mlst data frame (NCBI metadata + MLST)
source(paste0(dir, "scripts/03_global.phylogeny.mlst.all.genomes.R"))
figures=paste0(dir,"figures/04_global_phylogeny/")  # reset after sourcing script 03

# Build rooted fan tree with custom labels, outgroup marking, top ST/country rings --------
# Select tree ------

tree <- read.tree(paste0(datadir, "global_phylogeny/PA_global_tree.treefile")); name="final"

# Plot tree ------

# Collapse near-zero branches into polytomies
tol <- 1e-6
tree$edge.length[tree$edge.length < tol] <- 0
tree <- di2multi(tree, tol = tol)

outgroup.root <- c("GCA_000007565.2", "GCA_000021045.1", "GCA_019704535.1")

tree <- root(tree, outgroup = outgroup.root, resolve.root = FALSE)

tip.df <- tibble(label = tree$tip.label)

meta.public <- meta.mlst %>%
  transmute(
    label = `Assembly Accession`,
    country = ifelse(is.na(country) | country == "" | country == "unknown", "unknown", country),
    Collection.date = `Assembly BioSample Collection date`,
    ST = as.character(ST),
    isolate = label,
    dataset = "public"
  )

samples.PA[samples.PA$MLST_type == "773!", "MLST_type"] <- "773"
my.isolates.ST = unique(samples.PA[c(4,5)])

# Build metadata for own isolates and attach the correct ST from my.isolates.ST
meta.own <- tip.df %>%
  filter(grepl("^THJ", label)) %>%
  transmute(
    label,
    country = str_extract(label, "(?<=\\|)[^|]+(?=\\|[^|]+$)"),
    Collection.date = NA_character_,
    isolate = str_extract(label, "(?<=\\|287\\|)\\d+(?=\\|)"),
    dataset = "my_isolates"
  ) %>%
  left_join(
    my.isolates.ST %>% transmute(isolate = as.character(isolate), ST = as.character(MLST_type)),
    by = "isolate"
  )

tip.meta <- tip.df %>%
  left_join(meta.public, by = "label") %>%
  rows_update(meta.own, by = "label") %>%
  mutate(
    dataset = ifelse(is.na(dataset), "public", dataset),
    country.plot = ifelse(is.na(country) | country == "", "unknown", country),
    ST.plot = ifelse(is.na(ST) | ST == "" | ST == "-", "unknown", ST),
    isolate = ifelse(is.na(isolate), label, isolate)
  )

ref.map <- c(
  "GCA_000006765.1" = "PAO1",
  "GCA_000014625.1" = "PA14",
  "GCA_000026645.1" = "LESB58",
  "GCA_000017205.1" = "PA7"
  #"GCA_003025345.2" = "ASM302534v2"
  )

outgroup.vec <- c(
  "GCA_000019225.1", "GCA_000021045.1", "GCA_019704535.1", "GCA_900215245.1",  "GCA_000007565.2"
  #"GCA_003025345.2",#"GCA_000017205.1", 
)

outgroup.vec <- c("GCA_000007565.2", "GCA_000021045.1", "GCA_019704535.1","GCA_003025345.2")

tip.meta <- tip.meta %>%
  mutate(
    isolate = case_when(
      label %in% names(ref.map) ~ unname(ref.map[label]),
      TRUE ~ isolate
    ),
    outgroup = label %in% outgroup.vec,
    ref.label = case_when(
      label %in% names(ref.map) ~ unname(ref.map[label]),
      TRUE ~ NA_character_
    ),
    tip.class = case_when(
      outgroup ~ "outgroup",
      label %in% names(ref.map) ~ "reference",
      TRUE ~ "ingroup"
    )
  )

top.country <- tip.meta %>% count(country.plot, sort = TRUE) %>% slice_head(n = 10) %>% pull(country.plot)
top.st <- tip.meta %>%   filter(ST.plot != "unknown") %>% count(ST.plot, sort = TRUE) %>% slice_head(n = 30) %>% pull(ST.plot)

tip.meta <- tip.meta %>%
  mutate(
    country.plot = ifelse(country.plot %in% top.country, country.plot, "Other"),
    ST.plot = ifelse(ST.plot %in% top.st, ST.plot, "Other")
  )

country.levels <- tip.meta %>% count(country.plot, sort = TRUE) %>% pull(country.plot)
country.levels <- c(setdiff(country.levels, "Other"), "Other")

st.levels <- tip.meta %>% count(ST.plot, sort = TRUE) %>% pull(ST.plot)
st.levels <- c(setdiff(st.levels, "Other"), "Other")

tip.meta <- tip.meta %>%
  mutate(
    country.plot = factor(country.plot, levels = country.levels),
    ST.plot = factor(ST.plot, levels = st.levels)
  )

country.levels <- levels(tip.meta$country.plot)
country.cols <- setNames(grDevices::hcl.colors(length(country.levels), "Dynamic"), country.levels)
#country.cols <- setNames(createPalette(length(country.levels), seedcolors = c("#5DA5DA", "#F15854")), country.levels)
#country_vivid <- all_cols[hcl_mat[, "C"] > 20]  # reuse if already computed above
#country.cols <- setNames(country_vivid[seq_along(country.levels)], country.levels)
all_cols <- polychrome(36)
hcl_mat  <- coords(as(hex2RGB(all_cols), "polarLUV"))
vivid    <- all_cols[hcl_mat[, "C"] > 20]
country_base <- vivid[round(seq(1, length(vivid), length.out = length(country.levels)))]
country.cols <- setNames(desaturate(country_base, amount = 0.2), country.levels)
country.cols <- setNames(qualitative_hcl(length(country.levels), h = c(10, 260), c = 65, l = 68),country.levels)
country.cols <- setNames(grDevices::hcl.colors(length(country.levels), "Dark 3"), country.levels)

country.cols["Other"] <- "grey80"

st.levels <- levels(tip.meta$ST.plot)
#st.cols <- setNames(grDevices::hcl.colors(length(st.levels), "Dark 3"), st.levels)
#st.cols <- setNames(createPalette(length(st.levels), seedcolors = c("#60BD68", "#B276B2")), st.levels)
#st.cols <- setNames(qualitative_hcl(length(st.levels), h = c(10, 320), c = 60, l = 72),st.levels)
#st.cols <- setNames(polychrome(length(st.levels)), st.levels)

all_cols <- polychrome(36)
hcl_mat  <- coords(as(hex2RGB(all_cols), "polarLUV"))
vivid    <- all_cols[hcl_mat[, "C"] > 20]   # drop low-chroma (grey) entries

base_cols <- vivid[seq_along(st.levels)]
st.cols <- setNames(desaturate(base_cols, amount = 0.2), st.levels)  # 0.1–0.3 range

st.cols["Other"] <- "grey80"

#p.base.bl <- ggtree(tree, layout = "fan", open.angle = 20, linewidth = 0.1)

p.base <- ggtree(tree, layout = "fan", branch.length = "none", open.angle = 20, linewidth = 0.1)

tip.plot <- p.base$data %>%
  filter(isTip) %>%
  left_join(tip.meta, by = "label")

ring.country <- tip.meta %>% transmute(label, country.ring = country.plot) %>% distinct()
ring.st <- tip.meta %>% transmute(label, ST.ring = ST.plot) %>% distinct()

overlap.labels <- c("3", "4", "25", "32")

p <- p.base +
  geom_tippoint(
    data = tip.plot %>% filter(dataset == "public" & !outgroup),
    color = "grey80", size = 0.1, alpha = 0.8
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(!is.na(ref.label) & dataset == "public"),
    color = "blue", size = 1.5, alpha = 0.9
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(dataset == "my_isolates"),
    color = "red", size = 1.5, alpha = 0.95
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(outgroup),
    shape = 17, color = "black", size = 1.5, alpha = 1
  ) +
  geom_tiplab2(
    data = tip.plot %>% filter(dataset == "my_isolates", !isolate %in% overlap.labels),
    aes(label = isolate),
    color = "red", size = 4, fontface = "plain",
    offset = 1.2, align = FALSE
  ) +
  geom_tiplab2(
    data = tip.plot %>% filter(dataset == "my_isolates", isolate == "3"),
    aes(label = isolate),
    color = "red", size = 4,
    offset = 1.2, align = FALSE
  ) +
  geom_tiplab2(
    data = tip.plot %>% filter(dataset == "my_isolates", isolate == "4"),
    aes(label = isolate),
    color = "red", size = 4,
    offset = 4.8, align = FALSE
  ) +
  geom_tiplab2(
    data = tip.plot %>% filter(dataset == "my_isolates", isolate == "25"),
    aes(label = isolate),
    color = "red", size = 4,
    offset = 1.2, align = FALSE
  ) +
  geom_tiplab2(
    data = tip.plot %>% filter(dataset == "my_isolates", isolate == "32"),
    aes(label = isolate),
    color = "red", size = 4,
    offset = 5.8, align = FALSE
  ) +
  geom_tiplab2(
    data = tip.plot %>% filter(!is.na(ref.label)),
    aes(label = ref.label),
    color = "blue4", size = 2, fontface = "bold",
    offset = 1.2, align = FALSE
  ) +
  geom_fruit(
    data = ring.st, geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = ST.ring),
    inherit.aes = FALSE, offset = 0.25, pwidth = 0.1, color = NA, width = 7
  ) +
  scale_fill_manual(values = st.cols, name = "ST", drop = FALSE) +
  ggnewscale::new_scale_fill() +
  geom_fruit(
    data = ring.country, geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = country.ring),
    inherit.aes = FALSE, offset = 0.08, pwidth = 0.1, color = NA, width = 3
  ) +
  scale_fill_manual(values = country.cols, name = "Country", drop = FALSE) +
  theme_void() +
  theme(
    legend.text = element_text( size = 11),
    legend.title = element_text(size=11,face="bold")
  )

p

ggsave(filename = paste0(figures, "PA_global_tree_",name,".pdf"),plot = p, width = 9, height = 9)

# Save same tree without any tip labels ----------
p <- p.base +
  geom_tippoint(
    data = tip.plot %>% filter(dataset == "public" & !outgroup),
    color = "grey80", size = 0.1, alpha = 0.8
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(!is.na(ref.label) & dataset == "public"),
    color = "blue", size = 1.5, alpha = 0.9
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(dataset == "my_isolates"),
    color = "red", size = 1.5, alpha = 0.95
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(outgroup),
    shape = 17, color = "black", size = 2.5, alpha = 1
  ) +
  geom_fruit(
    data = ring.st, geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = ST.ring),
    inherit.aes = FALSE, offset = 0.25, pwidth = 0.1, color = NA, width = 7
  ) +
  scale_fill_manual(values = st.cols, name = "ST", drop = FALSE) +
  ggnewscale::new_scale_fill() +
  geom_fruit(
    data = ring.country, geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = country.ring),
    inherit.aes = FALSE, offset = 0.08, pwidth = 0.1, color = NA, width = 3
  ) +
  scale_fill_manual(values = country.cols, name = "Country", drop = FALSE) +
  theme_void() +
  theme(
    legend.text = element_text( size = 11),
    legend.title = element_text(size=11,face="bold")
  )

p

ggsave(filename = paste0(figures, "PA_global_tree_without_labels_",name,".pdf"),plot = p, width = 9, height = 9)

# Rectangular phylogram with transformed branch lengths --------
tree.plot <- tree
tree.plot$edge.length <- sqrt(tree.plot$edge.length)

p.rect <- ggtree(tree.plot, layout = "rectangular", linewidth = 0.1)

tip.plot <- p.rect$data %>%
  filter(isTip) %>%
  left_join(tip.meta, by = "label")

p.rect <- p.rect +
  geom_tippoint(
    data = tip.plot %>% filter(dataset == "public" & !outgroup),
    color = "grey80", size = 0.5, alpha = 0.8
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(ST.plot == "549" & dataset == "public"),
    color = "blue", size = 1.6, alpha = 0.9
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(dataset == "my_isolates"),
    color = "red", size = 2.2, alpha = 0.95
  ) +
  geom_tippoint(
    data = tip.plot %>% filter(outgroup),
    shape = 17, color = "black", size = 2.5, alpha = 1
  ) +
  geom_tiplab(
    data = tip.plot %>% filter(dataset == "my_isolates"),
    aes(label = isolate),
    color = "red", size = 3, align = TRUE, linesize = 0.2
  ) +
  geom_tiplab(
    data = tip.plot %>% filter(!is.na(ref.label)),
    aes(label = ref.label),
    color = "blue4", size = 3, fontface = "bold", offset = 0.002
  ) +
  theme_tree2()

p.rect
#ggsave(filename = paste0(figures, "PA_tree_rectangular.",name,".pdf"),plot = p.rect, width = 8, height = 10)

# Rectangular tree with bootstrap support values on nodes ----

tree <- tree.plot

# Parse support from node labels
parse_support <- function(tree) {
  n_tips <- Ntip(tree)
  tibble(
    node = (n_tips + 1):(n_tips + Nnode(tree)),
    label = tree$node.label
  ) %>%
    mutate(
      label = ifelse(is.na(label) | label == "", "0/0", label),
      alrt  = as.numeric(gsub("/.*", "", label)),
      uboot = as.numeric(gsub(".*/", "", label))
    )
}

support_df <- parse_support(tree)

# Get node depths from root
n_tips <- Ntip(tree)

internal_nodes <- tibble(
  node = (n_tips + 1):(n_tips + Nnode(tree))
) %>%
  left_join(support_df, by = "node") %>%
  mutate(
    depth = node.depth(tree)[(n_tips + 1):(n_tips + Nnode(tree))],
    clade_size = sapply(node, function(n) {
      length(Descendants(tree, n, type = "tips")[[1]])
    })
  ) %>%
  arrange(depth)

# Define "main" as nodes subtending at least X% of ingroup tips
n_ingroup <- sum(!tree$tip.label %in% outgroup.vec)
threshold <- 0.03  # at least 3% of tips = ~45 tips in your 1500 genome tree

main_splits <- internal_nodes %>%
  filter(clade_size >= threshold * n_ingroup) %>%
  arrange(depth) %>%
  mutate(
    support_label = paste0(round(alrt), "/", round(uboot)),
    well_supported = alrt >= 80 & uboot >= 95
  )

# Join support data to tree plot data first
tree_plot_data <- p.rect$data %>%
  left_join(main_splits %>% 
              select(node, alrt, uboot, well_supported, 
                     clade_size, support_label), 
            by = "node")

p.rect.supp <- p.rect +
  geom_text(
    data = tree_plot_data %>%
      filter(node %in% main_splits$node,
             clade_size > 0.05 * n_ingroup),
    aes(x = x, y = y,
        label = round(uboot),
        color = alrt >= 80),
    size = 2, hjust = -0.1
  ) +
  scale_color_manual(values = c("TRUE"="black", "FALSE"="red"))
p.rect.supp
ggsave(filename = paste0(figures, "PA_tree_rectangular_with_support.",name,".pdf"),plot = p.rect.supp, width = 8, height = 10)

