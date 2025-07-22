library(Biostrings)  # readDNAStringSet()
library(ggplot2)     # plotting

# ------------------------------------------------------------------
# 1. Data loaders ---------------------------------------------------
# ------------------------------------------------------------------

# Read all variant sequences from FASTA
load_variants <- function(fasta_path) {
  readDNAStringSet(fasta_path)
}

# Create one data-frame per variant by calling initialize_df()
# Outputs a list of dataframes, one for each variant
get_dfs <- function(seq_names, fasta_path, raw_dir) {
  all_variants <- load_variants(fasta_path)
  dfs <- lapply(seq_names, function(name) {
    npy_fn <- file.path(raw_dir, sprintf("input_%s_logits.npy", name))
    initialize_df(npy_fn, name, all_variants)
  })
  names(dfs) <- seq_names
  dfs
}

# Sum log-likelihood column for every df and return a named numeric vector
# Input a list of dataframes, one for each variant (output of get_dfs())
get_tot_lls <- function(gene_dfs) {
  sapply(gene_dfs, function(df) sum(df$log_likelihood))
}

# ------------------------------------------------------------------
# 2. Plot helpers ----------------------------------------------------
# ------------------------------------------------------------------

# Bar plot of relative likelihood per variant
# takes relative likelihoods as input (output of 2^get_tot_lls() normalized)
plot_norm_ll <- function(rel_ls) {
  ggplot(data.frame(gene = names(rel_ls), rel_l = rel_ls),
         aes(x = gene, y = rel_l)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(x = "Gene variant", y = "Relative likelihood",
         title = "Normalized likelihood per gene") +
    theme(axis.text.y = element_text(size = 8))
}

# Bar plot of relative likelihood summed by amino-acid
# takes relative likelihoods as input (output of 2^get_tot_lls() normalized)
# hard coded codon table 11 because only working with bacteria
plot_aa_ls <- function(rel_ls) {
  codon_table_11 <- getGeneticCode("11")
  aa_rel_ls <- setNames(numeric(length(unique(codon_table_11))),
                        unique(codon_table_11))
  for (name in names(rel_ls)) {
    codon <- substr(name, 4, 6)
    aa <- codon_table_11[codon]
    aa_rel_ls[aa] <- aa_rel_ls[aa] + rel_ls[name]
  }
  ggplot(data.frame(amino_acids = names(aa_rel_ls), val = aa_rel_ls),
         aes(x = val, y = amino_acids)) +
    geom_col(fill = "steelblue") +
    theme_minimal(base_size = 14) +
    labs(x = "Relative AA likelihood", y = "Amino acids",
         title = "Normalized likelihood per amino-acid")
}

plot_FR_comp_ls <- function(for_rel_ls, rev_rel_ls) {
  plot_df <- data.frame(
    forward = for_rel_ls, 
    reverse = rev_rel_ls,
    codon = substr(names(for_rel_ls), 4, 6)  # Extract codon from name (e.g., "83_AAA" -> "AAA")
  )
  ggplot(plot_df, aes(x = forward, y = reverse)) + 
    geom_point(color = "steelblue", alpha = 0.6) +
    geom_text(aes(label = codon)) +
    labs(x = "Relative Likelihood (Forward)", 
         y = "Relative Likelihood (RC)",
         title = "Relative Likelihood: Forward vs Reverse Complement Variants (Position 83)") +
    theme_minimal()
}

plot_FR_lls <- function(for_tot_lls, rev_tot_lls) {
  # Load required package
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install the ggrepel package: install.packages('ggrepel')")
  }
  
  # Create data frame with codon names as a separate column
  plot_df <- data.frame(
    forward = for_tot_lls,  # Use total log likelihoods directly
    reverse = rev_tot_lls,  # Use total log likelihoods directly
    codon = substr(names(for_tot_lls), 4, 6)  # Extract codon from name (e.g., "83_AAA" -> "AAA")
  )
  ggplot(plot_df, aes(x = forward, y = reverse)) + 
    geom_point(color = "steelblue", alpha = 0.6) +
    ggrepel::geom_text_repel(
      aes(label = codon),
      size = 3,
      box.padding = 0.5,     # Padding around the labels
      point.padding = 0.2,   # Padding between the label and its point
      min.segment.length = 0, # Allow shorter connecting segments
      max.overlaps = Inf,    # Don't hide any labels
      seed = 42,             # For reproducible label placement
      force = 2,             # Increase repulsion force between labels
      segment.color = "grey50", # Color of the connecting lines
      segment.alpha = 0.5     # Transparency of connecting lines
    ) +
    labs(x = "Log Likelihood (Forward)", 
         y = "Log Likelihood (RC)",
         title = "Log Likelihood: Forward vs Reverse Complement Variants (Position 83)") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 11),
      axis.title = element_text(size = 10)
    )
}

# Save both plots via save_as_pdf() (defined in evo2_analysis_functions)
plot_all <- function(rel_ls, out_dir = "output") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  save_as_pdf(plot_norm_ll(rel_ls), "Norm_codon_ls.pdf", out_dir)
  save_as_pdf(plot_aa_ls(rel_ls), "Norm_aa_ls.pdf", out_dir)
}

plot_FR_all <- function(for_rel_ls, rev_rel_ls, for_tot_lls, rev_tot_lls, out_dir = "output") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  # Save normalized likelihood plots
  plot_all(for_rel_ls, out_dir)
  save_as_pdf(plot_norm_ll(rev_rel_ls), "Rev_norm_codon_ls.pdf", out_dir)
  save_as_pdf(plot_aa_ls(rev_rel_ls), "Rev_norm_aa_ls.pdf", out_dir)
  save_as_pdf(plot_FR_comp_ls(for_rel_ls, rev_rel_ls), "FR_comp_ls.pdf", out_dir)
  save_as_pdf(plot_FR_lls(for_tot_lls, rev_tot_lls), "FR_tot_lls.pdf", out_dir)

  # Display plots
  print(plot_norm_ll(for_rel_ls))
  print(plot_aa_ls(for_rel_ls))
  print(plot_norm_ll(rev_rel_ls))
  print(plot_aa_ls(rev_rel_ls))
  print(plot_FR_comp_ls(for_rel_ls, rev_rel_ls))
  print(plot_FR_lls(for_tot_lls, rev_tot_lls))
}


