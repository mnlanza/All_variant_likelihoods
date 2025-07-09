# === evo2_variant_main.R ===
library(Biostrings)
library(ggplot2)
library(optparse)

parse_args   <- optparse::parse_args

# ---- CLI options ----
option_list <- list(
  make_option(c("-f", "--fasta"),      type="character", help="FASTA file"),
  make_option(c("-r", "--raw_dir"),    type="character", help="Logits folder"),
  make_option(c("-s", "--save"),       action="store_true", default=FALSE,
              help="Save output PDFs"),
  make_option(c("-l", "--highlight"),  type="character",
              default="247,248,249",
              help="Comma-separated positions to highlight")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- source helper scripts ----
source("handle_EVO2_output_copy.r")           # defines save_as_pdf()
source("evo2_analysis_functions_copy.R")      # defines initialize_df()
source("gyrae_83_all.R")

# ---- load data & build dfs ----
codon_table_11 <- getGeneticCode("11")

all_variants <- load_variants(opt$fasta)
seq_names    <- names(all_variants)

gene_dfs <- get_dfs(seq_names, opt$fasta, opt$raw_dir)    # << pass args
tot_lls  <- get_tot_lls(gene_dfs)

rel_ls <- 2^(tot_lls - max(tot_lls))
rel_ls <- rel_ls / sum(rel_ls)

# ---- plot ----
print(plot_norm_ll(rel_ls))
print(plot_aa_ls(rel_ls, codon_table_11, seq_names))

if (opt$save) {
  plot_all(rel_ls, codon_table_11, seq_names)             # << pass args
  cat("✅ Plots saved to My_figures/\n")
} else {
  cat("ℹ️  Use --save to write plots to PDF files.\n")
}