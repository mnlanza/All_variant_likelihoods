# initializing
library(Biostrings)
library(ggplot2)

# Suppress warnings for externally defined functions
utils::globalVariables(c("initialize_df", "save_as_pdf"))

# Source external functions first
source("~/Desktop/Relman_lab/EVO2_Exercise/handle_EVO2_output.r")  # Contains save_as_pdf
source("~/Desktop/Relman_lab/EVO2_Exercise/ML_EVO2_Exercise.R")    # Contains initialize_df

codon_table_11 <- getGeneticCode('11')
setwd("~/Desktop/Relman_lab/All_variant_likelihoods")

all_variants <- readDNAStringSet("codon83_variants.fasta") 
seq_names <- names(all_variants)

# don't reuse has specific file name
get_dfs <- function(){
  dfs <- list()
  for (name in seq_names){
    npy_fn <- paste0("raw_data/input_", name, "_logits.npy")
    dfs[[name]] <- initialize_df(npy_fn, name)
  }
  return(dfs)
}

# returns a vector with total ll for each gene
get_tot_lls <- function(all_gene_dfs){
  # initialize vec of length all_genes_dfs to store total ll 
  total_ll_vec <- numeric(length(all_gene_dfs))
  names(total_ll_vec) <- names(all_gene_dfs)
  for (i in 1:length(all_gene_dfs)) { 
    total_ll_vec[i] <- sum(all_gene_dfs[[i]]$log_likelihood)
  }
  return(total_ll_vec)
}

# storing all gene dfs
all_gene_dfs <- get_dfs()

# vector with names and total_ll for each gene
tot_lls <- get_tot_lls(all_gene_dfs)

# subtracting lowest ll from all lls
lowered_ll <- tot_lls - max(tot_lls)

# turn ll in to likelihood
gene_ls_vec <- 2^lowered_ll

# relative likelihoods
rel_ls <- gene_ls_vec / sum(gene_ls_vec)


plot_norm_ll <- function(){

  # turn into a dataframe for plotting
  df_norm <- data.frame(
    gene   = names(rel_ls),
    rel_l = rel_ls
  )
  
  # ggplot bar‐plot
  ggplot(df_norm, aes(x = gene, y = rel_l)) +
      geom_col(fill = "steelblue") +
      coord_flip() + 
      theme_minimal(base_size = 16) +
      labs(x = "Gene variant", y = "Relative likelihood", 
           title = "Normalized (relative) likelihood per gene") +
      theme(axis.text.y = element_text(size = 8))
    
}

# plotting normalized log-likelihoods
print(plot_norm_ll())

###########################
# sorting by amino acid ###
###########################
rel_ls["83_ATG"]
which.max(rel_ls)

# getting vector with amino acids
amino_acids <- unique(as.vector(codon_table_11))
# making a new vector with amino acids as names and 0's for relative likelihoods
aa_rel_ls <- setNames(numeric(length(amino_acids)), amino_acids)
# iterating to add relative likelihood values to aa_counts
for (name in seq_names){
  codon <- substr(name, 4, 6)
  aa <- codon_table_11[codon]
  aa_rel_ls[aa] <- aa_rel_ls[aa] + rel_ls[name]
}

plot_aa_ls <- function(){
  df_norm_aa <- data.frame(
    amino_acids = names(aa_rel_ls),
    aa_rel_ls = aa_rel_ls
  )
  ggplot(df_norm_aa, aes(x = aa_rel_ls, y = amino_acids)) + 
      geom_col(fill = "steelblue") +     
      theme_minimal(base_size = 14) +
      labs(x = "Relative amino acid likelihood", y = "Amino acids", 
             title = "Normalized (relative) amino acids per gene") 
}

# plotting amino acid likelihoods
print(plot_aa_ls())

plot_all <- function(){
  save_as_pdf(plot_norm_ll(), "Norm_codon_ls.pdf", out_dir = "My_figures", width = 10, height = 6.7)
  save_as_pdf(plot_aa_ls(), "Norm_aa_ls.pdf", out_dir = "My_figures")
}
#plot_all()