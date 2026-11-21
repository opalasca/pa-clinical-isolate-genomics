# Set folders -----------
dir <- "./"  # run all scripts from the repository root

scripts=paste0(dir,"scripts/")
datadir=paste0(dir,"data/")
figures=paste0(dir,"figures/03_global_st_distribution/")
resultdir=paste0(dir,"results/")
ifelse(!dir.exists(figures), dir.create(figures, recursive = TRUE), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)

# Load packages -------------

list.of.packages <- c("ggplot2","dplyr","tidyr","readxl","stringr","readr","openxlsx", "ape","tibble","ggtree","ggtreeExtra","readr","viridis","ggrepel",
                      "sf","rnaturalearth","countrycode","rnaturalearthdata")

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
#write.table(samples.PA, paste0(resultdir,"metadata.PA.isolates.tsv"),sep="\t", quote=FALSE, row.names=FALSE)
#samples.PA$niche = "wound"

# Read full metadata ----

# Read a tab-separated metadata file and keep original column names
meta.full <- read.delim(paste0(datadir,"global_phylogeny/PA_all_ncbi_genomes_metadata.tsv"), sep = "\t", header = TRUE, check.names = FALSE, quote = "", fill = TRUE, na.strings = c("", "missing"))

table(meta.full$`Assembly Level`)

meta.complete <- meta.full[meta.full$`Assembly Level` == "Complete Genome", ]

# Count assemblies per year
table(format(as.Date(meta.complete$`Assembly Release Date`), "%Y"))

# Filter by quality
meta.filtered <- meta.full[meta.full$`Assembly Stats Total Sequence Length` > 5.8e6 & meta.full$`Assembly Stats Total Sequence Length` < 7.2e6, ]
meta.filtered <- meta.filtered[meta.filtered$`Assembly Stats Number of Contigs` <= 500 & meta.filtered$`Assembly Stats Contig N50` >= 20000,]
# Keep only Genbank assemblies
meta.filtered <- meta.filtered[grepl("^GCA_", meta.filtered$`Assembly Accession`), ]
# Save filtered accessions for downstream use on computerome
write.table(meta.filtered[,1], paste0(resultdir,"PA_accessions.all.quality.filtered.tsv"),sep="\t", quote=FALSE, row.names=FALSE)

# Clean country names from Assembly BioSample Geographic location -------
meta.filtered <- meta.filtered %>%
  mutate(
    #country.raw = `Assembly BioSample Geographic location`,
    country = `Assembly BioSample Geographic location`
  ) %>%
  mutate(
    country = trimws(country),
    country = na_if(country, "NA"),
    country = na_if(country, "unknown"),
    country = na_if(country, "Unknown"),
    country = na_if(country, "missing"),
    country = na_if(country, "Missing"),
    country = na_if(country, "")
  ) %>%
  mutate(
    country = sub(":.*", "", country),
    country = trimws(country)
  ) %>%
  mutate(
    country = case_when(
      country == "USA" ~ "United States",
      country == "Viet Nam" ~ "Vietnam",
      TRUE ~ country
    )
  ) %>%
  mutate(
    country = tidyr::replace_na(country, "unknown")
  )

meta.filtered[tolower(meta.filtered$country) %in% c("unknown","missing","not applicable","not determined","not provided","not collected","not available"), ]$country = "unknown"


# Read MLST results for all genomes ----

# Read a tab-separated metadata file and keep original column names
mlst.all <- read.delim(paste0(datadir,"global_phylogeny/mlst_ncbi_genomes.tsv"), sep = "\t", header = TRUE, check.names = FALSE, quote = "", fill = TRUE, na.strings = c("", "missing"))
mlst.PA.isolates <- read.delim(paste0(datadir,"global_phylogeny/mlst_our_isolates.tsv"), sep = "\t", header = TRUE, check.names = FALSE, quote = "", fill = TRUE, na.strings = c("", "missing"))

meta.mlst.all = merge(meta.filtered, mlst.all, by.x="Assembly Accession", by.y="genome", all.x=TRUE)
meta.mlst = meta.mlst.all %>% filter(species == "paeruginosa") 

ST.counts = meta.mlst %>%
  filter(!ST=="-") %>%
  group_by(ST) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

meta.mlst.ukraine <- meta.mlst %>% filter(country == "Ukraine" | country == "ukraine")

ST.counts.ukraine = meta.mlst.ukraine %>%
  filter(!ST=="-") %>%
  group_by(ST) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

meta.mlst.ourSTs = meta.mlst %>% filter(ST %in% mlst.PA.isolates$ST & !ST=="-" )

ST.counts.ourSTs = meta.mlst.ourSTs %>%
  filter(!ST=="-") %>%
  group_by(ST) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

meta.mlst.ourSTs.complete <- meta.mlst.ourSTs %>%
  filter(`Assembly Level` == "Complete Genome") 
  
meta.mlst.complete <- meta.mlst %>% filter(`Assembly Level` == "Complete Genome")

ST.counts.complete = meta.mlst.complete %>%
  filter(!ST=="-") %>%
  group_by(ST) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

# Select only assemblies with known country and clean up country names for mapping
meta.mlst.loc <- meta.mlst[!tolower(meta.mlst$country) %in% c("unknown","missing","not applicable","not determined","not provided","not collected","not available"), ]
meta.mlst.loc <- meta.mlst.loc %>% filter(!is.na(country))

meta.mlst.loc <- meta.mlst.loc %>%mutate(ST = ifelse(ST != "-", paste0("ST", ST), ST))

# Bubble map of total counts by continent for selected STs
# Plot world map with our STs (one colored bubble per ST around each continent center, sized by assembly count) ---------

our.STs <- paste0("ST",c("773","207","1047","253","235","369","1129","480","270","244"))

# Create regions with Asia split based on your data distribution
meta.cont <- meta.mlst.loc %>%
  mutate(
    continent = countrycode(country, "country.name", "continent"),
    region = case_when(
      continent == "Europe" ~ "Europe",
      continent == "Africa" ~ "Africa",
      continent == "Oceania" ~ "Oceania",
      
      # Americas split
      continent == "Americas" & country %in% c("United States","Canada","Mexico","Greenland") ~ "North America",
      continent == "Americas" ~ "South America",
      
      # Asia split (data-driven)
      continent == "Asia" & country %in% c(
        "China","Japan","South Korea","North Korea","Taiwan","Hong Kong",
        "Indonesia","Singapore","Thailand","Vietnam","Malaysia","Philippines","Myanmar"
      ) ~ "East/Southeast Asia",
      
      continent == "Asia" ~ "South/West Asia",
      
      TRUE ~ NA_character_
    )
  )


cont.centers <- data.frame(
  region = c("North America","South America","Europe","Africa",
             "East/Southeast Asia","South/West Asia","Oceania"),
  lon = c(-100, -60, 15, 20, 120, 60, 135),
  lat = c(45, -15, 55, 5, 30, 25, -25)
)


plot.df <- meta.cont %>%
  filter(ST %in% our.STs, !is.na(region)) %>%
  count(region, ST, name = "n") %>%
  mutate(ST = factor(ST, levels = our.STs)) %>%   
  left_join(cont.centers, by = "region") %>%
  group_by(region) %>%
  arrange(ST, .by_group = TRUE) %>%
  mutate(
    k = n(),
    i = row_number(),
    angle = 2 * pi * (i - 1) / k,
    lon.plot = lon + 15 * cos(angle),
    lat.plot = lat + 9 * sin(angle)
  ) %>%
  ungroup()

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_point(data = plot.df, aes(lon.plot, lat.plot, size = n, fill = ST), shape = 21, color = "black", alpha = 0.85) +
  geom_text(data = cont.centers, aes(lon, lat, label = region), fontface = "bold", size = 3) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  theme_void() +
  labs(size = "Assemblies", fill = "ST")

# Create circular legend for STs
legend.df <- data.frame(
  ST = factor(our.STs, levels = our.STs)
) %>%
  mutate(
    i = row_number(),
    k = n(),
    angle = 2 * pi * (i - 1) / k,
    lon = -130 + 35 * cos(angle),   # position legend (adjust)
    lat = -30 + 25  * sin(angle)
  )

p <- ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  
  # your data bubbles
  geom_point(data = plot.df,
             aes(lon.plot, lat.plot, size = n, fill = ST),
             shape = 21, color = "black", alpha = 0.85) +
  
  # circular legend
  geom_point(data = legend.df,
             aes(lon, lat, fill = ST),
             shape = 21, size = 4, color = "black") +
  
  geom_text(data = legend.df,
            aes(lon, lat, label = ST),
            size = 3.5, vjust = -1) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  scale_size_continuous(
    name = "Assemblies",
    trans = "sqrt",
    breaks = c(10, 50, 100, 250, 500, 1000),
    range = c(1, 6)
  ) +
  theme_void() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = "none")   # remove default legend

p
ggsave(p, filename = paste0(figures,"world.map.ourSTs.pdf"), width = 10, height = 6)

# Summary of ST by country, region and total counts -------

ST.by.country <- meta.cont %>%
  mutate(ST = as.character(ST)) %>%
  filter(!is.na(country), !is.na(region), !is.na(ST), ST != "-") %>%
  count(ST, country, region, name = "n") %>%   
  arrange(desc(n), country)

ST.by.region <- meta.cont %>%
  mutate(ST = as.character(ST)) %>%
  filter(!is.na(region), !is.na(ST), ST != "-") %>%
  count(ST, region, name = "n") %>%
  arrange(desc(n), region)

ST.counts <- meta.cont %>%
  filter(!is.na(ST), ST != "-") %>%
  count(ST, name = "n") %>%
  arrange(desc(n))

region.order <- meta.cont %>%
  filter(!is.na(region), ST != "-") %>%
  count(region, name = "total") %>%
  arrange(desc(total)) %>%
  pull(region)

ST.region.summary <- meta.cont %>%
  mutate(ST = as.character(ST)) %>%
  filter(!is.na(region), !is.na(ST), ST != "-") %>%
  count(ST, region, name = "n") %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(-ST))) %>%
  relocate(total, .after = ST) %>%
  select(ST, total, all_of(region.order)) %>%   # 👈 key line
  arrange(desc(total))

write.xlsx(ST.by.country, paste0(resultdir,"ST_by_country.xlsx"), rowNames = FALSE)
write.xlsx(ST.region.summary, paste0(resultdir,"ST_region_summary.xlsx"), rowNames = FALSE)

write.table(ST.by.country, paste0(resultdir,"ST_by_country.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(ST.region.summary, paste0(resultdir,"ST_region_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)



# Subsample assemblies for global phylogeny proportional to ST representation ---------
set.seed(1)

our.STs <- c("773","207","1047","253","235","369","1129","480","270","244")

meta.for.subsampling <- meta.mlst.loc %>% filter(!ST %in% c("_","-")) 

meta.sub <- meta.for.subsampling |>
  mutate(ST = as.character(ST)) |>
  filter(!is.na(ST), ST != "")

st.counts <- meta.sub |>
  count(ST, name = "count") |>
  mutate(quota = case_when(
    count < 5 ~ 0L,
    TRUE ~ pmax(1, floor(0.02 * count))
  ),
  quota = pmin(quota, count))

meta.sub.main <- meta.sub |>
  left_join(st.counts |> select(ST, count, quota), by = "ST") |>
  filter(count >= 5) |>
  group_by(ST) |>
  group_modify(~ slice_sample(.x, n = .x$quota[1])) |>
  ungroup() |>
  select(-count, -quota)

n.main <- nrow(meta.sub.main)

meta.sub.rare <- meta.sub |>
  left_join(st.counts |> select(ST, count), by = "ST") |>
  filter(count < 5) |>
  select(-count)

n.to.add <- max(0, 1500 - n.main)

meta.sub.rare.sampled <- if(n.to.add > 0){
  slice_sample(meta.sub.rare, n = min(n.to.add, nrow(meta.sub.rare)))
} else {
  meta.sub.rare[0, ]
}

meta.subsampled <- bind_rows(meta.sub.main, meta.sub.rare.sampled)

ST.counts.subsampled = meta.subsampled %>%
  filter(!ST=="-") %>%
  group_by(ST) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

a=merge(ST.counts, ST.counts.subsampled, by="ST", all.x=TRUE)

ref.strains =c("GCA_000006765.1", "GCA_000014625.1", "GCA_000026645.1")

outgroup.strains = c("GCA_000017205.1", "GCA_000019225.1", "GCA_000021045.1", "GCA_003025345.2", "GCA_019704535.1", "GCA_900215245.1")

sample.list = c(ref.strains, outgroup.strains, meta.subsampled$`Assembly Accession`)

write.table(sample.list, paste0(resultdir,"PA_accessions.subsampled.1500.plus.ref.tsv"),sep="\t", quote=FALSE, row.names=FALSE,col.names = FALSE)

