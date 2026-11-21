# Run all analysis scripts in order ------------------------------------------
# Usage: open pa-clinical-isolate-genomics.Rproj in RStudio and source this file,
#        or from the repository root: Rscript run_all.R
#
# The output of each script is written to logs/.
# Figures are written to figures/ and tables to results/.

scripts <- c(
  "01_preprocess.quality.check.R",
  "02_analysis.phylogeny.pseudomonas.R",
  "03_global.phylogeny.mlst.all.genomes.R",
  "04_global.phylogeny.tree.visualisation.R",   # sources script 03
  "05_analysis.abricate.pseudomonas.R",
  "06_analysis.resistance.genes.pseudomonas.R",    # needs results of script 05
  "07_analysis.mutations.resistance.R",
  "08_QS.phenotype.and.mutations.R",             # needs results of script 07
  "09_supplementary_tables.R"
)

if (!file.exists("pa-clinical-isolate-genomics.Rproj")) {
  stop("Run this script from the repository root (or open the .Rproj file first).")
}

# ggtree 3.x and patchwork < 1.3.1 do not work with ggplot2 >= 4.0
if (requireNamespace("ggplot2", quietly = TRUE) && requireNamespace("ggtree", quietly = TRUE) &&
    packageVersion("ggplot2") >= "4.0.0" && packageVersion("ggtree") < "4.0.0") {
  stop("ggplot2 ", packageVersion("ggplot2"), " is not compatible with ggtree ", packageVersion("ggtree"), ".\n",
       "Install ggplot2 3.5.2 and restart R:\n",
       '  remotes::install_version("ggplot2", version = "3.5.2", repos = "https://cloud.r-project.org")')
}

dir.create("logs", showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

for (s in scripts) {
  log.file <- file.path("logs", sub("\\.R$", ".log", s))
  message(sprintf("Running %-48s", s), appendLF = FALSE)
  t0 <- Sys.time()
  status <- system2(rscript, file.path("scripts", s), stdout = log.file, stderr = log.file)
  elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (status != 0) {
    message(sprintf("FAILED after %ds", elapsed))
    stop("Script ", s, " failed. See ", log.file, " for details.")
  }
  message(sprintf("done (%ds)", elapsed))
}

message("\nAll scripts finished. Figures are in figures/, tables in results/.")
