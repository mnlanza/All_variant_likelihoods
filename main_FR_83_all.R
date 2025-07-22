suppressPackageStartupMessages({
  library(Biostrings)
  library(ggplot2)
  library(optparse)
})

# ---- CLI options ----
option_list <- list(
  make_option(c("-f", "--for_fasta"), type="character", help="Forward FASTA file path"),
  make_option(c("--rev_fasta"), type="character", help="Reverse complement FASTA file path"),
  make_option(c("-r", "--for_raw_dir"), type="character", help="Forward logits folder path"),
  make_option(c("--rev_raw_dir"), type="character", help="Reverse complement logits folder path"),
  make_option(c("-s", "--save"), action="store_true", default=FALSE, help="Save output PDFs"),
  make_option(c("-l", "--highlight"), type="character", default="247,248,249",
              help="Comma-separated positions to highlight"),
  make_option(c("-o", "--output_dir"), type="character", default="output",
              help="Output directory for plots")
)
opt <- parse_args(OptionParser(option_list = option_list))

# sets raw directory for reverse strand to same as raw data for forward strand if no input
if (is.null(opt$rev_raw_dir)) {
  opt$rev_raw_dir <- opt$for_raw_dir
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output_dir)) {
  dir.create(opt$output_dir, recursive = TRUE)
}

# Prevent creation of Rplots.pdf
if (!interactive()) {
  pdf(NULL)
}

# ---- source helper scripts ----
script_dir <- file.path("scripts")
source(file.path(script_dir, "handle_EVO2_output.r"))  # defines save_as_pdf()
source(file.path(script_dir, "evo2_analysis_functions.R"))  # defines initialize_df()
source(file.path(script_dir, "83_all_fx.R"))

for_variants <- load_variants(opt$for_fasta)
rev_variants <- load_variants(opt$rev_fasta)

for_seq_names <- names(for_variants)
rev_seq_names <- names(rev_variants)

for_gene_dfs <- get_dfs(for_seq_names, opt$for_fasta, opt$for_raw_dir)
rev_gene_dfs <- get_dfs(rev_seq_names, opt$rev_fasta, opt$rev_raw_dir)

for_tot_lls <- get_tot_lls(for_gene_dfs)
rev_tot_lls <- get_tot_lls(rev_gene_dfs)

for_tot_ls_minusC <- 2^(for_tot_lls - max(for_tot_lls))
rev_tot_ls_minusC <- 2^(rev_tot_lls - max(rev_tot_lls))

for_rel_ls <- for_tot_ls_minusC / sum(for_tot_ls_minusC)
rev_rel_ls <- rev_tot_ls_minusC / sum(rev_tot_ls_minusC)

if (opt$save) {
  # Save all plots using the updated plot_FR_all function
  plot_FR_all(for_rel_ls, rev_rel_ls, for_tot_lls, rev_tot_lls, out_dir = opt$output_dir)
  cat(sprintf("✅ Plots saved to %s/\n", opt$output_dir))
} else {
  cat("ℹ️  Use --save to write plots to PDF files.\n")
}

# Close any open graphics devices
while (dev.cur() > 1) {
  dev.off()
}



