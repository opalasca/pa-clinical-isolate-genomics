# pa-clinical-isolate-genomics

R scripts and processed data for the genomic and phenotypic analysis of *Pseudomonas aeruginosa* clinical isolates described in:

> Palasca O, Schwartz CH, Kristensen R, Koshova I, Murashevych B, Thomsen TR, Hertz FB, Bjarnsholt T, Jakobsen TH.
> **Molecular epidemiology and virulence profiles of MDR and XDR *Pseudomonas aeruginosa* from combat-related wounds in Ukraine.**
> *New Microbes and New Infections* (2026). https://doi.org/10.1016/j.nmni.2026.101880

Raw sequencing reads and assemblies are available at NCBI under BioProject **PRJNA1473102**.

## Contents

```
run_all.R                           runs all scripts in order
data/                               processed inputs (P. aeruginosa isolates only)
scripts/                            R scripts, numbered in the order they should be run
pa-clinical-isolate-genomics.Rproj  RStudio project (sets the working directory)
```

`figures/`, `results/` and `logs/` are created when the scripts run. Each script writes its figures to its own subfolder (e.g. `figures/05_virulence_genes/`); the combined variant barplot (Fig. 6) is written by script 08 to `figures/08_quorum_sensing/`.

## Scripts

To run the full analysis, open `pa-clinical-isolate-genomics.Rproj` in RStudio and source `run_all.R`, or run from the repository root:

```
Rscript run_all.R
```

`run_all.R` runs the scripts in numerical order and writes one log per script to `logs/` (about 3 minutes in total).
Individual scripts can also be run on their own from the repository root, as long as the scripts they depend on have been run first.
Missing packages are installed automatically via `BiocManager`.
Scripts 04 and 08 need script 03 and 07 outputs, respectively; script 06 needs the output of script 05.

Notes:
- ggtree 3.x (the version available for R 4.4) does not work with ggplot2 ≥ 4.0. Install ggplot2 3.5.2 before running:
  `remotes::install_version("ggplot2", version = "3.5.2", repos = "https://cloud.r-project.org")`
  (`run_all.R` checks this and stops with this message if needed.)

| Script | Analysis | Figure / table |
|---|---|---|
| `01_preprocess.quality.check.R` | Sample metadata, MLST types and QUAST assembly statistics; flags assemblies > 8 Mb | Methods |
| `02_analysis.phylogeny.pseudomonas.R` | SNP-based phylogeny (BacDist: Snippy, ClonalFrameML, RAxML) of the isolates and PAO1: one tree with both sequencing duplicates, and one with a single assembly per isolate; defines the isolate order used in all heatmaps | Fig. 1B, 3, 5 |
| `03_global.phylogeny.mlst.all.genomes.R` | Quality filtering of public NCBI genomes, ST/country summaries, world map of the study STs, ST-proportional subsampling of 1,500 genomes | Fig. 2A–B |
| `04_global.phylogeny.tree.visualisation.R` | Global IQ-TREE phylogeny (40 marker genes) annotated with ST and country (sources script 03) | Fig. 2C |
| `05_analysis.abricate.pseudomonas.R` | Virulence genes (ABRicate/VFDB), duplicate-confirmed presence/absence; heatmaps of all genes and of genes with variable presence, grouped by VFDB category | Fig. 4A (basis), Suppl. |
| `06_analysis.resistance.genes.pseudomonas.R` | Acquired resistance genes: RGI, ResFinder (StarAMR) and AMRFinderPlus, keeping genes detected in both sequencing duplicates by at least one tool; heatmap of these genes grouped by class; resistance and virulence gene counts per isolate (needs output of script 05) | Fig. 3A (basis), 4B |
| `07_analysis.mutations.resistance.R` | Disc-diffusion susceptibility; moderate/high-impact variants in resistance-associated genes relative to PAO1; variant counts per isolate | Fig. 1B, 3B, 6 |
| `08_QS.phenotype.and.mutations.R` | Las/Rhl reporter assays, pyocyanin and rhamnolipid production; variants in QS-associated genes; combined variant barplots (uses panels saved by script 07) | Fig. 5, 6 |
| `09_supplementary_tables.R` | Supplementary tables: acquired resistance genes (same duplicate rule as script 06), virulence genes, resistance and QS mutations | Suppl. tables |

The final figures were assembled and edited outside R; for some panels (Fig. 3A, 4A) the subset of genes shown was selected manually from the heatmaps produced here.
Statistical comparisons of phenotypes with PAO1 (one-way ANOVA with Dunnett's test) were not done in these scripts.

## Data

All files are restricted to the *P. aeruginosa* isolates in the study (13 isolates, each sequenced in duplicate; sample IDs `S3`–`S69`) and reference genomes.

| Path | Content |
|---|---|
| `WGS samplesheet.xlsx` | Sequencing sample ↔ isolate mapping: `sample`, `plate.well`, `isolate.replicate` (`3a`/`3b` = duplicates of isolate 3), `Species` |
| `quast_transposed_report.tsv` | QUAST assembly statistics (Shovill assemblies) |
| `MLST_results_p_aeruginosa/` | MLST (CGE MLST 2.0) per sample |
| `RGI/` | RGI (CARD) results per assembly, incl. reference genomes |
| `Resfinder_plasmidfinder.tsv.tabular`, `Staramr_reference_genomes.tabular` | StarAMR (ResFinder/PlasmidFinder) results for isolates and reference genomes |
| `AMRFinderPlus/` | AMRFinderPlus results per assembly |
| `ABRicate/` | ABRicate (VFDB) results per assembly and reference genome |
| `selected.snps.CARD*.txt` | Resistance-associated variants reported in CARD |
| `disc_diffusion_results.xlsx` | Disc-diffusion susceptibility categories (EUCAST) |
| `Phenotypic_data.xlsx` | Las and Rhl reporter assays, pyocyanin and rhamnolipid measurements |
| `bacdist/` | BacDist outputs: SNP-based RAxML tree incl. PAO1 (recombinant regions identified by ClonalFrameML removed), merged raw Snippy VCF (gzipped), PAO1 annotation (GFF) |
| `global_phylogeny/PA_all_ncbi_genomes_metadata.tsv` | NCBI metadata for all *P. aeruginosa* assemblies (downloaded 05/03/2026) |
| `global_phylogeny/mlst_ncbi_genomes.tsv` | MLST (`mlst` command-line tool, run on Computerome) for all quality-filtered public genomes, computed **before** subsampling so that genomes could be subsampled proportionally to ST |
| `global_phylogeny/mlst_our_isolates.tsv` | MLST (`mlst` command-line tool, run on Computerome) for the study isolates, for consistency with the public genomes |
| `global_phylogeny/PA_global_tree.treefile` | Final IQ-TREE maximum-likelihood tree (Fig. 2C; IQ-TREE seed 1111) |

## Upstream processing

Read processing, assembly (Shovill v1.4.2) and SNP calling (BacDist) were run on a high-performance computing cluster.
RGI was run via the CARD web interface; StarAMR, AMRFinderPlus and ABRicate via Galaxy.
For the global phylogeny, sequence types were first assigned to all downloaded public genomes and to the study isolates with the `mlst` command-line tool on Computerome (`data/global_phylogeny/mlst_*.tsv`). Script 03 uses these STs to subsample 1,500 genomes proportionally to ST abundance, and the subsampled genomes were then used to build the tree (Prodigal, fetchMGs, MAFFT, trimAl, IQ-TREE) with [core-genome-phylogeny](https://github.com/opalasca/core-genome-phylogeny).
Tool versions and parameters are given in the Methods section of the paper.

## Software

Tested with R 4.4.3. Main packages: ggplot2 3.5.2, dplyr 1.2.0, tidyr 1.3.2, readxl 1.4.5, openxlsx 4.2.8, ComplexHeatmap 2.22.0, circlize 0.4.16, ggtree 3.14.0, ape 5.8.1, vcfR 1.15.0, rtracklayer 1.66.0, patchwork 1.3.0.
